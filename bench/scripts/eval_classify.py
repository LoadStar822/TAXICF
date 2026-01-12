#!/usr/bin/env python3
"""
Rank-aware evaluation for TAXICF classify output using NCBI taxonomy (ete3).

This script is copied/adapted from TAXICF/bench/scripts/eval_classify.py.

Assumptions:
- Prediction TSV: columns <read_id> <taxid:score ...>; first taxid is used. "unclassified" means no call.
- Truth TSV: header contains "#anonymous_contig_id" or "read_id" and "tax_id".
- Taxonomy: provided via ete3 NCBITaxa. Uses an existing local sqlite DB if present.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    from ete3 import NCBITaxa
except ImportError as e:  # pragma: no cover - runtime check
    print("[error] ete3 is required. Install with `pip install ete3`.", file=sys.stderr)
    raise


def _split_whitespace_row(row: List[str]) -> List[str]:
    # Sometimes the "TSV" is actually space-delimited (e.g., CAMI *_gsa_mapping_new.tsv).
    if len(row) == 1 and row[0] and ("\t" not in row[0]) and (" " in row[0]):
        return row[0].split()
    return row


def load_truth(path: Path) -> Dict[str, str]:
    """Load truth mapping.

    Supported formats:
    - CAMI gsa_mapping.tsv: has header with "tax_id" and "#anonymous_contig_id" (or "read_id").
    - Simple mapping TSV: header "read_id\ttax_id".
    - SIMULATED mapping: no header; typically 3 columns <read_id> <ignored> <taxid> (taxid constant per file).
    """
    truth: Dict[str, str] = {}
    with path.open() as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return truth
        header = _split_whitespace_row(header)

        # Headered format (Cami / generic).
        if "tax_id" in header:
            try:
                id_idx = header.index("#anonymous_contig_id")
            except ValueError:
                try:
                    id_idx = header.index("read_id")
                except ValueError:
                    raise SystemExit(
                        f"truth file missing id column (#anonymous_contig_id/read_id): {path}"
                    )
            tax_idx = header.index("tax_id")
            for row in reader:
                row = _split_whitespace_row(row)
                if len(row) <= tax_idx or len(row) <= id_idx:
                    continue
                truth[row[id_idx]] = row[tax_idx]
            return truth

        # Headerless format (SIMULATED-like). Treat the first line as data.
        first = header
        if len(first) == 2:
            truth[first[0]] = first[1]
        elif len(first) >= 3:
            truth[first[0]] = first[2]

        for row in reader:
            row = _split_whitespace_row(row)
            if len(row) == 2:
                truth[row[0]] = row[1]
            elif len(row) >= 3:
                truth[row[0]] = row[2]
    return truth


def load_pred(path: Path) -> Dict[str, str]:
    preds: Dict[str, str] = {}
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rid = parts[0]
            second = parts[1]
            if second == "unclassified":
                preds[rid] = "unclassified"
            else:
                taxid = second.split(":", 1)[0]
                preds[rid] = taxid
    return preds


def build_ncbi(dbfile: Optional[Path], taxdump: Optional[Path]) -> NCBITaxa:
    """Instantiate NCBITaxa, building the sqlite DB from a local taxdump (no download)."""
    ncbi = NCBITaxa(dbfile=str(dbfile) if dbfile else None)

    need_build = dbfile and not dbfile.exists()
    need_build = need_build or (dbfile is None and not Path(ncbi.dbfile).exists())

    if need_build:
        if not taxdump or not taxdump.is_dir():
            raise SystemExit("taxdump directory required to build ete3 database")
        with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tf:
            with tarfile.open(tf.name, "w:gz") as tar:
                for fname in taxdump.iterdir():
                    if fname.is_file():
                        tar.add(fname, arcname=fname.name)
            tar_path = tf.name
        ncbi.update_taxonomy_database(taxdump_file=tar_path)
    return ncbi


def ancestors_at_ranks(
    ncbi: NCBITaxa, taxids: Iterable[int], ranks: List[str]
) -> Dict[int, Dict[str, Optional[int]]]:
    """Return mapping taxid -> {rank: ancestor_taxid_at_rank or None}."""
    cache: Dict[int, Dict[str, Optional[int]]] = {}
    for tid in taxids:
        try:
            lineage = ncbi.get_lineage(tid)
        except Exception:
            cache[tid] = {r: None for r in ranks}
            continue
        rank_map = ncbi.get_rank(lineage)
        rank_lookup = {}
        for anc in lineage:
            r = rank_map.get(anc)
            if r in ranks and r not in rank_lookup:
                rank_lookup[r] = anc
        cache[tid] = {r: rank_lookup.get(r) for r in ranks}
    return cache


def build_lineage_cache(ncbi: NCBITaxa, taxids: Iterable[int]) -> Dict[int, set[int]]:
    """Return mapping taxid -> set(lineage taxids)."""
    cache: Dict[int, set[int]] = {}
    for tid in taxids:
        try:
            lineage = set(ncbi.get_lineage(tid))
        except Exception:
            lineage = set()
        cache[tid] = lineage
    return cache


def is_descendant(pred_tid: int, true_tid: int, lineage_cache: Dict[int, set[int]]) -> bool:
    lineage = lineage_cache.get(pred_tid)
    if lineage:
        return true_tid in lineage
    return pred_tid == true_tid


def descendant_counts(
    truth: Dict[str, str],
    preds: Dict[str, str],
    lineage_cache: Dict[int, set[int]],
) -> Tuple[int, int, int]:
    tp = fp = fn = 0
    for rid, true_tid_str in truth.items():
        try:
            true_tid = int(true_tid_str)
        except ValueError:
            continue
        pred_tid_str = preds.get(rid)
        if pred_tid_str is None or pred_tid_str == "unclassified":
            fn += 1
            continue
        try:
            pred_tid = int(pred_tid_str)
        except ValueError:
            fp += 1
            fn += 1
            continue
        if is_descendant(pred_tid, true_tid, lineage_cache):
            tp += 1
        else:
            fp += 1
            fn += 1
    return tp, fp, fn


def prf(tp: int, fp: int, fn: int) -> Dict[str, float]:
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"precision": prec, "recall": rec, "f1": f1}


def compute_desc_metrics(
    truth: Dict[str, str],
    preds: Dict[str, str],
    lineage_cache: Dict[int, set[int]],
) -> Dict[str, float]:
    tp, fp, fn = descendant_counts(truth, preds, lineage_cache)
    base = prf(tp, fp, fn)
    accuracy = tp / len(truth) if truth else 0.0
    return {
        "precision": base["precision"],
        "recall": base["recall"],
        "f1": base["f1"],
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "accuracy": accuracy,
    }


def score_from_abundance(data: dict) -> Tuple[str, float]:
    if "desc_presence_f1" in data:
        return "abundance-presence-f1", float(data["desc_presence_f1"])
    if "presence_f1" in data:
        return "abundance-presence-f1", float(data["presence_f1"])
    if "l1_distance_pct" in data:
        return "abundance-l1", -float(data["l1_distance_pct"])
    return "none", float("-inf")


def score_from_per_read(data: dict) -> Tuple[str, float]:
    if "desc_f1" in data:
        return "per-read-f1", float(data["desc_f1"])
    if "f1" in data:
        return "per-read-f1", float(data["f1"])
    return "none", float("-inf")


def main() -> None:
    ap = argparse.ArgumentParser(description="Rank-aware evaluation of TAXICF classification output.")
    ap.add_argument("--pred", required=True, type=Path, help="TAXICF classify TSV output")
    ap.add_argument("--truth", required=True, type=Path, help="Ground-truth mapping TSV")
    ap.add_argument("--out", required=True, type=Path, help="metrics.json output path")
    ap.add_argument(
        "--taxdump",
        type=Path,
        default=Path("/mnt/sda/tianqinzhong/project/SimDataset/taxdump"),
        help="Directory containing nodes.dmp & names.dmp",
    )
    ap.add_argument("--ete-db", type=Path, default=None, help="ete3 sqlite DB path (will be created if missing)")
    ap.add_argument(
        "--ranks",
        nargs="*",
        default=["species", "genus", "family", "order", "class", "phylum", "superkingdom"],
        help="Ranks to evaluate",
    )
    ap.add_argument(
        "--primary-rank",
        choices=["auto", "species", "genus"],
        default="auto",
        help="Top-level rank reported as precision/recall/F1 (auto prefers species, falls back to first rank)",
    )
    ap.add_argument(
        "--collapse-to-rank",
        choices=["none", "species", "genus"],
        default="none",
        help="If set, collapse both truth and prediction to this rank before scoring",
    )
    args = ap.parse_args()

    truth = load_truth(args.truth)
    preds = load_pred(args.pred)

    ncbi = build_ncbi(args.ete_db, args.taxdump)

    all_taxids = {int(t) for t in truth.values()} | {
        int(t) for t in preds.values() if t not in ("unclassified", None)
    }

    rank_cache = ancestors_at_ranks(ncbi, all_taxids, args.ranks)
    lineage_cache = build_lineage_cache(ncbi, all_taxids)

    # per-rank counters
    counters = {r: {"tp": 0, "fp": 0, "fn": 0, "support": 0} for r in args.ranks}
    unclassified = 0
    missing = 0

    for rid, true_tid_str in truth.items():
        true_tid = int(true_tid_str)
        pred_tid_str = preds.get(rid)
        if pred_tid_str is None:
            missing += 1
            pred_tid = None
        elif pred_tid_str == "unclassified":
            pred_tid = None
            unclassified += 1
        else:
            pred_tid = int(pred_tid_str)

        ranks_to_use = [args.collapse_to_rank] if args.collapse_to_rank != "none" else args.ranks

        for rank in ranks_to_use:
            true_at = rank_cache.get(true_tid, {}).get(rank)
            if true_at is None:
                continue  # skip if truth missing at this rank
            counters[rank]["support"] += 1

            if pred_tid is None:
                counters[rank]["fn"] += 1
                continue

            pred_at = rank_cache.get(pred_tid, {}).get(rank)
            if pred_at is None:
                counters[rank]["fp"] += 1
                counters[rank]["fn"] += 1
            elif pred_at == true_at:
                counters[rank]["tp"] += 1
            else:
                counters[rank]["fp"] += 1
                counters[rank]["fn"] += 1

    per_rank = {}
    for rank, c in counters.items():
        support = c["support"]
        per_rank[rank] = {
            **c,
            **prf(c["tp"], c["fp"], c["fn"]),
            "accuracy": c["tp"] / support if support else 0.0,
        }

    if args.primary_rank == "auto":
        primary_rank = "species" if per_rank.get("species", {}).get("support", 0) else args.ranks[0]
    else:
        primary_rank = args.primary_rank
    top = per_rank.get(
        primary_rank,
        {"tp": 0, "fp": 0, "fn": 0, "support": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0},
    )
    desc = compute_desc_metrics(truth, preds, lineage_cache)

    exact_accuracy = top["tp"] / top["support"] if top["support"] else 0.0
    metrics = {
        "metric_version": "descendant-aware-v1",
        "reads_total": len(truth),
        "reads_pred": len(preds),
        "truth_entries": len(truth),
        "unclassified": unclassified,
        "missing_pred": missing,
        "per_rank": per_rank,
        "collapse_to_rank": None if args.collapse_to_rank == "none" else args.collapse_to_rank,
        "precision": desc["precision"],
        "recall": desc["recall"],
        "f1": desc["f1"],
        "tp": desc["tp"],
        "fp": desc["fp"],
        "fn": desc["fn"],
        "accuracy": desc["accuracy"],
        "desc_precision": desc["precision"],
        "desc_recall": desc["recall"],
        "desc_f1": desc["f1"],
        "desc_tp": desc["tp"],
        "desc_fp": desc["fp"],
        "desc_fn": desc["fn"],
        "desc_accuracy": desc["accuracy"],
        "exact_precision": top["precision"],
        "exact_recall": top["recall"],
        "exact_f1": top["f1"],
        "exact_tp": top["tp"],
        "exact_fp": top["fp"],
        "exact_fn": top["fn"],
        "exact_accuracy": exact_accuracy,
        "primary_rank": primary_rank,
        "desc": {
            "precision": desc["precision"],
            "recall": desc["recall"],
            "f1": desc["f1"],
            "tp": desc["tp"],
            "fp": desc["fp"],
            "fn": desc["fn"],
            "accuracy": desc["accuracy"],
        },
        "exact": {
            "rank": primary_rank,
            "precision": top["precision"],
            "recall": top["recall"],
            "f1": top["f1"],
            "tp": top["tp"],
            "fp": top["fp"],
            "fn": top["fn"],
            "support": top.get("support", 0),
            "accuracy": exact_accuracy,
        },
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        json.dump(metrics, f, indent=2)

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
