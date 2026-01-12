#!/usr/bin/env python3
"""
Summarize TAXICF pipeline runs into a single TSV.

This scans `bench/runs/pipeline/**/run.json` and, if present, loads the
corresponding `metrics.json`.

Note: output is intended to be local-only (bench/results is gitignored).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional


def _safe_float(x: Any) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None


def _safe_int(x: Any) -> Optional[int]:
    try:
        return int(x)
    except Exception:
        return None


def _load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _db_stem(db_path: str) -> str:
    p = Path(db_path)
    return p.stem


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize TAXICF pipeline runs under bench/runs/pipeline.")
    ap.add_argument("--runs-root", type=Path, default=Path("bench/runs/pipeline"), help="Runs root to scan.")
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("bench/results/pipeline_summary.tsv"),
        help="Output TSV path (default is gitignored).",
    )
    args = ap.parse_args()

    runs_root: Path = args.runs_root
    out_path: Path = args.out

    rows: List[Dict[str, Any]] = []
    for run_json in sorted(runs_root.rglob("run.json")):
        try:
            meta = _load_json(run_json)
        except Exception:
            continue

        run_dir = run_json.parent
        metrics_path = Path(meta.get("metrics") or (run_dir / "metrics.json"))
        metrics: Dict[str, Any] = {}
        if metrics_path.is_file():
            try:
                metrics = _load_json(metrics_path)
            except Exception:
                metrics = {}

        desc = metrics.get("desc") or {}
        exact = metrics.get("exact") or {}

        row: Dict[str, Any] = {
            "run_dir": str(run_dir.resolve()),
            "timestamp": meta.get("timestamp"),
            "db": meta.get("db"),
            "db_stem": _db_stem(str(meta.get("db") or "")),
            "input": meta.get("input"),
            "truth": meta.get("truth"),
            "threads": _safe_int(meta.get("threads")),
            "elapsed_sec": _safe_float(meta.get("elapsed_sec")),
            "returncode": _safe_int(meta.get("returncode")),
            "eval_returncode": _safe_int(meta.get("eval_returncode")),
            "reads_total": _safe_int(metrics.get("reads_total")),
            "reads_pred": _safe_int(metrics.get("reads_pred")),
            "unclassified": _safe_int(metrics.get("unclassified")),
            "desc_precision": _safe_float(desc.get("precision")),
            "desc_recall": _safe_float(desc.get("recall")),
            "desc_f1": _safe_float(desc.get("f1")),
            "exact_rank": exact.get("rank"),
            "exact_precision": _safe_float(exact.get("precision")),
            "exact_recall": _safe_float(exact.get("recall")),
            "exact_f1": _safe_float(exact.get("f1")),
        }
        rows.append(row)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "timestamp",
        "run_dir",
        "db_stem",
        "input",
        "truth",
        "threads",
        "elapsed_sec",
        "returncode",
        "eval_returncode",
        "reads_total",
        "reads_pred",
        "unclassified",
        "desc_precision",
        "desc_recall",
        "desc_f1",
        "exact_rank",
        "exact_precision",
        "exact_recall",
        "exact_f1",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})

    print(f"[OK] wrote {len(rows)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

