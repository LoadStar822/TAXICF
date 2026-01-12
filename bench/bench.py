#!/usr/bin/env python3
"""
TAXICF benchmark driver (inspired by TAXICF/bench).

Goals:
- Provide a unified entry to run `TAXICF build` / `TAXICF classify`
- Record command/env/logs under `bench/runs/`
- Keep large artifacts out of git (see .gitignore)

This script itself only uses Python stdlib. Optional metric evaluation may
invoke `bench/scripts/eval_classify.py` (requires `ete3`).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
BENCH_ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = BENCH_ROOT / "scripts"

DEFAULT_BIN = REPO_ROOT / "build" / "TAXICF"
DEFAULT_RUN_ROOT = BENCH_ROOT / "runs"
DEFAULT_DB_ROOT = BENCH_ROOT / "DB"
DEFAULT_THREADS = 32
DEFAULT_EVAL_SCRIPT = SCRIPTS_DIR / "eval_classify.py"

# Dataset roots/patterns (local machine; do not commit data).
CAMI_MARINE_LONG_FASTA = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMARINE_LONG/fasta/"
    "marmgCAMI2_sample_{sid}_contigs_anonymous_gsa.fasta"
)
CAMI_MARINE_LONG_TRUTH = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMARINE_LONG/mapping/"
    "marmgCAMI2_sample_{sid}_contigs_gsa_mapping.tsv"
)
CAMI_MARINE_SHORT_FASTA = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMARINE_SHORT/output_dir/"
    "marmgCAMI2_sample_{sid}_contigs_simulation_short_read_2018.08.15_09.49.32_sample_{sid}_contigs_anonymous_gsa.fasta"
)
CAMI_MARINE_SHORT_TRUTH = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMARINE_SHORT/output_dir/"
    "marmgCAMI2_sample_{sid}_contigs_simulation_short_read_2018.08.15_09.49.32_sample_{sid}_contigs_gsa_mapping.tsv"
)

CAMI_MOUSE_LONG_FASTA = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMOUSE_LONG/output_dir/"
    "2018.02.13_14.02.01_sample_{sid}_contigs_2018.02.13_14.02.01_sample_{sid}_contigs_anonymous_gsa.fasta"
)
CAMI_MOUSE_LONG_TRUTH = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMOUSE_LONG/output_dir/"
    "2018.02.13_14.02.01_sample_{sid}_contigs_2018.02.13_14.02.01_sample_{sid}_contigs_gsa_mapping.tsv"
)
CAMI_MOUSE_SHORT_FASTA = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMOUSE_SHORT/output_dir/"
    "2017.12.29_11.37.26_sample_{sid}_contigs_2017.12.29_11.37.26_sample_{sid}_contigs_anonymous_gsa.fasta"
)
CAMI_MOUSE_SHORT_TRUTH = Path(
    "/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMOUSE_SHORT/output_dir/"
    "2017.12.29_11.37.26_sample_{sid}_contigs_2017.12.29_11.37.26_sample_{sid}_contigs_gsa_mapping.tsv"
)

SIMULATED_ROOT = Path("/mnt/sda/tianqinzhong/project/taxicfBench/dataset/SIMULATED")


def _now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def _md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_short_rev(repo_root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=repo_root, text=True
        ).strip()
    except Exception:
        return "unknown"


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _run_cmd(
    cmd: List[str],
    *,
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> int:
    with open(stdout_path, "w", encoding="utf-8") as out, open(
        stderr_path, "w", encoding="utf-8"
    ) as err:
        p = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=err)
        return int(p.returncode)


def _normalize_db_path(db: str, db_root: Path) -> Path:
    p = Path(db)
    if p.is_absolute() or p.parent != Path("."):
        # user explicitly provided a path (absolute or contains path separators)
        if p.suffix == ".icf":
            return p
        return p.with_suffix(".icf")

    # treat as a DB name under DB_ROOT
    name = db
    if name.endswith(".icf"):
        name = name[: -len(".icf")]
    return (db_root / name).with_suffix(".icf")


@dataclass(frozen=True)
class BenchConfig:
    taxicf_bin: Path
    threads: int
    run_root: Path
    db_root: Path

    @staticmethod
    def from_env() -> "BenchConfig":
        return BenchConfig(
            taxicf_bin=Path(os.environ.get("TAXICF_BIN", str(DEFAULT_BIN))),
            threads=int(os.environ.get("THREADS", str(DEFAULT_THREADS))),
            run_root=Path(os.environ.get("RUN_ROOT", str(DEFAULT_RUN_ROOT))),
            db_root=Path(os.environ.get("DB_ROOT", str(DEFAULT_DB_ROOT))),
        )

    def ensure_bin(self) -> None:
        if not self.taxicf_bin.exists():
            raise FileNotFoundError(
                "TAXICF binary not found. Build it first or set TAXICF_BIN.\n"
                f"  default: {DEFAULT_BIN}\n"
                f"  current: {self.taxicf_bin}"
            )
        if not os.access(self.taxicf_bin, os.X_OK):
            raise PermissionError(f"TAXICF binary is not executable: {self.taxicf_bin}")


def _prepare_run_dir(cfg: BenchConfig, kind: str, name: str, tag: str) -> Path:
    run_dir = cfg.run_root / kind / name / f"{_timestamp()}_{tag}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def _write_env(cfg: BenchConfig, run_dir: Path) -> None:
    data = {
        "timestamp": _now_iso(),
        "git": _git_short_rev(REPO_ROOT),
        "taxicf_bin": str(cfg.taxicf_bin),
        "taxicf_md5": _md5sum(cfg.taxicf_bin) if cfg.taxicf_bin.exists() else "missing",
        "threads": cfg.threads,
        "hostname": os.uname().nodename if hasattr(os, "uname") else "unknown",
        "platform": sys.platform,
        "python": sys.version.replace("\n", " "),
    }
    _write_text(run_dir / "env.json", json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def _write_cmd(run_dir: Path, cmd: Iterable[str]) -> None:
    _write_text(run_dir / "command.txt", " ".join(cmd) + "\n")


def _split_fields(line: str) -> List[str]:
    """Split a line into fields (tab-delimited or whitespace-delimited variants)."""
    s = line.rstrip("\n")
    if "\t" in s:
        return s.split("\t")
    return s.split()


def _truth_has_header(path: Path) -> bool:
    """Best-effort check whether the truth file is headered (contains tax_id column)."""
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip():
                    continue
                return "tax_id" in _split_fields(line)
    except Exception:
        return False
    return False


def _merge_truth_files(truth_files: Sequence[Path], out_path: Path) -> Path:
    """Merge multiple truth TSVs into one file for eval.

    Supported:
    - CAMI gsa_mapping.tsv style (headered): keep header once, concatenate data rows.
    - SIMULATED style (no header): concatenate whole files.
    """
    if not truth_files:
        raise ValueError("truth_files is empty")

    headered = _truth_has_header(truth_files[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("wb") as out:
        for i, src in enumerate(truth_files):
            with src.open("rb") as f:
                first = f.readline()
                if not first:
                    continue

                if headered:
                    if i == 0:
                        header = "\t".join(_split_fields(first.decode("utf-8", errors="replace"))).encode("utf-8")
                        out.write(header + b"\n")
                    # Always skip the first line for headered inputs.
                else:
                    # Headerless: keep the first data line.
                    out.write(first)

                shutil.copyfileobj(f, out, length=1024 * 1024)

    return out_path


def _run_eval_script(
    *,
    pred_path: Path,
    truth_path: Path,
    out_path: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> int:
    script = Path(os.environ.get("EVAL_SCRIPT", str(DEFAULT_EVAL_SCRIPT))).resolve()
    if not script.is_file():
        raise FileNotFoundError(f"Eval script not found: {script}")

    cmd = [
        sys.executable,
        str(script),
        "--pred",
        str(pred_path),
        "--truth",
        str(truth_path),
        "--out",
        str(out_path),
    ]
    with open(stdout_path, "w", encoding="utf-8") as out, open(
        stderr_path, "a", encoding="utf-8"
    ) as err:
        p = subprocess.run(cmd, stdout=out, stderr=err)
        return int(p.returncode)


def _summarize_metrics(metrics_path: Path) -> str:
    try:
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        desc = metrics.get("desc") or {}
        exact = metrics.get("exact") or {}
        s = []
        if desc:
            s.append(
                "desc: "
                f"P={float(desc.get('precision', 0.0)):.4f} "
                f"R={float(desc.get('recall', 0.0)):.4f} "
                f"F1={float(desc.get('f1', 0.0)):.4f} "
                f"(tp={int(desc.get('tp', 0))} fp={int(desc.get('fp', 0))} fn={int(desc.get('fn', 0))})"
            )
        if exact:
            rank = str(exact.get("rank") or metrics.get("primary_rank") or "primary")
            s.append(
                f"exact@{rank}: "
                f"P={float(exact.get('precision', 0.0)):.4f} "
                f"R={float(exact.get('recall', 0.0)):.4f} "
                f"F1={float(exact.get('f1', 0.0)):.4f} "
                f"(tp={int(exact.get('tp', 0))} fp={int(exact.get('fp', 0))} fn={int(exact.get('fn', 0))}"
                f" support={int(exact.get('support', 0))})"
            )
        return "\n".join(s).strip()
    except Exception:
        return ""


@dataclass(frozen=True)
class DatasetJob:
    key: str
    reads: Path
    truth: Path


@dataclass(frozen=True)
class DatasetBatchJob:
    key: str
    reads: List[Path]
    truths: List[Path]


def _jobs_cami(pattern_reads: Path, pattern_truth: Path, sample_ids: List[int], key_prefix: str) -> List[DatasetJob]:
    jobs: List[DatasetJob] = []
    for sid in sample_ids:
        reads = Path(str(pattern_reads).format(sid=sid))
        truth = Path(str(pattern_truth).format(sid=sid))
        jobs.append(DatasetJob(key=f"{key_prefix}_sample{sid}", reads=reads, truth=truth))
    return jobs


def _jobs_simulated(limit: Optional[int]) -> List[DatasetJob]:
    if not SIMULATED_ROOT.is_dir():
        return []
    fastas = sorted(SIMULATED_ROOT.glob("*.fasta"))
    if limit is not None:
        fastas = fastas[: max(0, int(limit))]
    jobs: List[DatasetJob] = []
    for fasta in fastas:
        truth = fasta.with_suffix(".tsv")
        if not truth.is_file():
            continue
        jobs.append(DatasetJob(key=f"SIMULATED_{fasta.stem}", reads=fasta, truth=truth))
    return jobs


def _batch_cami(pattern_reads: Path, pattern_truth: Path, sample_ids: List[int], key: str) -> DatasetBatchJob:
    reads = [Path(str(pattern_reads).format(sid=sid)) for sid in sample_ids]
    truths = [Path(str(pattern_truth).format(sid=sid)) for sid in sample_ids]
    return DatasetBatchJob(key=key, reads=reads, truths=truths)


def _batch_simulated(limit: Optional[int]) -> DatasetBatchJob:
    jobs = _jobs_simulated(limit)
    reads = [j.reads for j in jobs]
    truths = [j.truth for j in jobs]
    return DatasetBatchJob(key=f"SIMULATED_n{len(reads)}", reads=reads, truths=truths)


def _safe_name(s: str, max_len: int = 140) -> str:
    # Keep run_dir paths short and portable.
    out = []
    for ch in s:
        if ch.isalnum() or ch in ("-", "_", "."):
            out.append(ch)
        else:
            out.append("_")
    name = "".join(out)
    while "__" in name:
        name = name.replace("__", "_")
    return name[:max_len].strip("_") or "run"


def cmd_build(args: argparse.Namespace) -> int:
    cfg = BenchConfig.from_env()
    cfg.ensure_bin()
    cfg.db_root.mkdir(parents=True, exist_ok=True)
    cfg.run_root.mkdir(parents=True, exist_ok=True)

    db_out = _normalize_db_path(args.db, cfg.db_root).resolve()
    name = args.name or db_out.stem
    run_dir = _prepare_run_dir(cfg, "build", name, args.tag)

    input_path = Path(args.input).resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"--input not found: {input_path}")

    cmd = [
        str(cfg.taxicf_bin),
        "build",
        "-i",
        str(input_path),
        "-o",
        str(db_out.with_suffix("")),  # TAXICF build will append .icf
        "-t",
        str(args.threads if args.threads is not None else cfg.threads),
    ]
    if args.extra:
        cmd.extend(args.extra)

    _write_env(cfg, run_dir)
    _write_cmd(run_dir, cmd)

    start = time.perf_counter()
    rc = _run_cmd(cmd, cwd=run_dir, stdout_path=run_dir / "stdout.log", stderr_path=run_dir / "stderr.log")
    elapsed = time.perf_counter() - start

    metrics = {
        "tool": "TAXICF",
        "mode": "build",
        "timestamp": _now_iso(),
        "elapsed_sec": elapsed,
        "returncode": rc,
        "input": str(input_path),
        "db_out": str(db_out),
        "threads": int(args.threads if args.threads is not None else cfg.threads),
    }
    # Keep runtime info separate from TAXICF-style evaluation metrics (metrics.json).
    _write_text(run_dir / "run.json", json.dumps(metrics, indent=2, ensure_ascii=False) + "\n")

    if rc != 0:
        return rc

    if not db_out.is_file():
        raise RuntimeError(
            "Build finished but expected DB file not found:\n"
            f"  expected: {db_out}\n"
            f"  run_dir:  {run_dir}"
        )

    print(f"[OK] build -> {db_out}")
    print(f"run_dir: {run_dir}")
    return 0


def cmd_classify(args: argparse.Namespace) -> int:
    cfg = BenchConfig.from_env()
    cfg.ensure_bin()
    cfg.run_root.mkdir(parents=True, exist_ok=True)

    db_path = _normalize_db_path(args.db, cfg.db_root).resolve()
    if not db_path.is_file():
        raise FileNotFoundError(
            f"Database not found: {db_path}\n"
            "Tip: run `bench/bench.py build ...` first or set DB_ROOT."
        )

    name = args.name or Path(args.input).stem
    run_dir = _prepare_run_dir(cfg, "classify", name, args.tag)

    input_path = Path(args.input).resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"--input not found: {input_path}")

    out_tsv = (run_dir / "pred.tsv").resolve()
    cmd = [
        str(cfg.taxicf_bin),
        "classify",
        "-i",
        str(input_path),
        "-d",
        str(db_path),
        "-o",
        str(out_tsv),
        "-t",
        str(args.threads if args.threads is not None else cfg.threads),
    ]
    if args.extra:
        cmd.extend(args.extra)

    _write_env(cfg, run_dir)
    _write_cmd(run_dir, cmd)

    start = time.perf_counter()
    rc = _run_cmd(cmd, cwd=run_dir, stdout_path=run_dir / "stdout.log", stderr_path=run_dir / "stderr.log")
    elapsed = time.perf_counter() - start

    metrics = {
        "tool": "TAXICF",
        "mode": "classify",
        "timestamp": _now_iso(),
        "elapsed_sec": elapsed,
        "returncode": rc,
        "input": str(input_path),
        "db": str(db_path),
        "output": str(out_tsv),
        "threads": int(args.threads if args.threads is not None else cfg.threads),
    }
    # Keep runtime info separate from TAXICF-style evaluation metrics (metrics.json).
    _write_text(run_dir / "run.json", json.dumps(metrics, indent=2, ensure_ascii=False) + "\n")

    if rc != 0:
        return rc

    if not out_tsv.is_file():
        raise RuntimeError(
            "Classify finished but expected output TSV not found:\n"
            f"  expected: {out_tsv}\n"
            f"  run_dir:  {run_dir}"
        )

    print(f"[OK] classify -> {out_tsv}")
    print(f"run_dir: {run_dir}")
    return 0


def cmd_pipeline(args: argparse.Namespace) -> int:
    """classify + eval in one run_dir (TAXICF-style)."""
    cfg = BenchConfig.from_env()
    cfg.ensure_bin()
    cfg.run_root.mkdir(parents=True, exist_ok=True)

    db_path = _normalize_db_path(args.db, cfg.db_root).resolve()
    if not db_path.is_file():
        raise FileNotFoundError(f"Database not found: {db_path}")

    input_specs = args.input if isinstance(args.input, list) else [args.input]
    truth_specs = args.truth if isinstance(args.truth, list) else [args.truth]

    input_paths = [Path(p).resolve() for p in input_specs]
    if not input_paths:
        raise SystemExit("--input is empty")
    for p in input_paths:
        if not p.is_file():
            raise FileNotFoundError(f"--input not found: {p}")

    truth_paths = [Path(p).resolve() for p in truth_specs]
    if not truth_paths:
        raise SystemExit("--truth is empty")
    for p in truth_paths:
        if not p.is_file():
            raise FileNotFoundError(f"--truth not found: {p}")

    if len(truth_paths) not in (1, len(input_paths)):
        raise SystemExit(
            "pipeline truth must be either:\n"
            "  - one merged truth file, OR\n"
            "  - the same number of truth files as inputs\n"
            f"got inputs={len(input_paths)} truths={len(truth_paths)}"
        )

    default_name = (
        _safe_name(f"{input_paths[0].stem}__{db_path.stem}")
        if len(input_paths) == 1
        else _safe_name(f"{len(input_paths)}inputs__{db_path.stem}")
    )
    name = args.name or default_name
    run_dir = _prepare_run_dir(cfg, "pipeline", name, args.tag)

    out_tsv = (run_dir / "pred.tsv").resolve()

    truth_path = truth_paths[0]
    if len(truth_paths) > 1:
        truth_path = _merge_truth_files(truth_paths, (run_dir / "truth.tsv").resolve())
        _write_text(run_dir / "truth_inputs.txt", "\n".join(str(p) for p in truth_paths) + "\n")

    cmd = [str(cfg.taxicf_bin), "classify"]
    for p in input_paths:
        cmd.extend(["-i", str(p)])
    cmd.extend(
        [
            "-d",
            str(db_path),
            "-o",
            str(out_tsv),
            "-t",
            str(args.threads if args.threads is not None else cfg.threads),
        ]
    )
    if args.extra:
        cmd.extend(args.extra)

    _write_env(cfg, run_dir)
    _write_cmd(run_dir, cmd)
    _write_text(run_dir / "inputs.txt", "\n".join(str(p) for p in input_paths) + "\n")

    start = time.perf_counter()
    rc = _run_cmd(cmd, cwd=run_dir, stdout_path=run_dir / "stdout.log", stderr_path=run_dir / "stderr.log")
    classify_elapsed = time.perf_counter() - start

    run_meta = {
        "tool": "TAXICF",
        "mode": "pipeline",
        "timestamp": _now_iso(),
        "returncode": rc,
        "input": str(input_paths[0]),
        "inputs": [str(p) for p in input_paths],
        "db": str(db_path),
        "output": str(out_tsv),
        "truth": str(truth_path),
        "truth_inputs": [str(p) for p in truth_paths],
        "threads": int(args.threads if args.threads is not None else cfg.threads),
        "elapsed_sec": classify_elapsed,
    }

    if rc != 0:
        _write_text(run_dir / "run.json", json.dumps(run_meta, indent=2, ensure_ascii=False) + "\n")
        return rc

    if not out_tsv.is_file():
        raise RuntimeError(
            "Classify finished but expected output TSV not found:\n"
            f"  expected: {out_tsv}\n"
            f"  run_dir:  {run_dir}"
        )

    metrics_path = (run_dir / "metrics.json").resolve()
    eval_stdout = run_dir / "metrics.txt"
    eval_rc = 0
    try:
        eval_rc = _run_eval_script(
            pred_path=out_tsv,
            truth_path=truth_path,
            out_path=metrics_path,
            stdout_path=eval_stdout,
            stderr_path=run_dir / "stderr.log",
        )
    except Exception as e:
        # Non-fatal for pipeline orchestration; preserve classify results.
        eval_rc = 127
        with open(run_dir / "stderr.log", "a", encoding="utf-8") as err:
            err.write(f"[pipeline] eval failed: {e}\n")

    run_meta["eval_returncode"] = eval_rc
    if metrics_path.is_file():
        run_meta["metrics"] = str(metrics_path)
    _write_text(run_dir / "run.json", json.dumps(run_meta, indent=2, ensure_ascii=False) + "\n")

    print(f"[OK] pipeline -> {out_tsv}")
    print(f"run_dir: {run_dir}")
    if metrics_path.is_file():
        summary = _summarize_metrics(metrics_path)
        if summary:
            print(summary)
    return 0


def cmd_run_all(args: argparse.Namespace) -> int:
    cfg = BenchConfig.from_env()
    cfg.ensure_bin()
    cfg.run_root.mkdir(parents=True, exist_ok=True)

    # DB list (names under DB_ROOT or absolute paths).
    db_specs = args.dbs or [
        "completeONE_v2_tag8",
        "completeONE_v2_tag16",
        "complete_v2_tag8",
        "complete_v2_tag16",
    ]
    db_paths = []
    for db in db_specs:
        p = _normalize_db_path(db, cfg.db_root).resolve()
        if not p.is_file():
            raise FileNotFoundError(f"Database not found: {p}")
        db_paths.append(p)

    datasets = set(args.datasets or [])
    if not datasets:
        datasets = {
            "cami-marine-long",
            "cami-marine-short",
            "cami-mouse-long",
            "cami-mouse-short",
            "simulated",
        }

    marine_samples = args.marine_samples or list(range(10))
    mouse_samples = args.mouse_samples or list(range(64))

    batch_mode = args.batch_mode
    batch_jobs: List[DatasetBatchJob] = []
    jobs: List[DatasetJob] = []

    if batch_mode == "group":
        if "cami-marine-long" in datasets:
            batch_jobs.append(
                _batch_cami(
                    CAMI_MARINE_LONG_FASTA,
                    CAMI_MARINE_LONG_TRUTH,
                    marine_samples,
                    f"CAMIMARINE_LONG_n{len(marine_samples)}",
                )
            )
        if "cami-marine-short" in datasets:
            batch_jobs.append(
                _batch_cami(
                    CAMI_MARINE_SHORT_FASTA,
                    CAMI_MARINE_SHORT_TRUTH,
                    marine_samples,
                    f"CAMIMARINE_SHORT_n{len(marine_samples)}",
                )
            )
        if "cami-mouse-long" in datasets:
            batch_jobs.append(
                _batch_cami(
                    CAMI_MOUSE_LONG_FASTA,
                    CAMI_MOUSE_LONG_TRUTH,
                    mouse_samples,
                    f"CAMIMOUSE_LONG_n{len(mouse_samples)}",
                )
            )
        if "cami-mouse-short" in datasets:
            batch_jobs.append(
                _batch_cami(
                    CAMI_MOUSE_SHORT_FASTA,
                    CAMI_MOUSE_SHORT_TRUTH,
                    mouse_samples,
                    f"CAMIMOUSE_SHORT_n{len(mouse_samples)}",
                )
            )
        if "simulated" in datasets:
            batch_jobs.append(_batch_simulated(args.simulated_limit))
    else:
        if "cami-marine-long" in datasets:
            jobs.extend(_jobs_cami(CAMI_MARINE_LONG_FASTA, CAMI_MARINE_LONG_TRUTH, marine_samples, "CAMIMARINE_LONG"))
        if "cami-marine-short" in datasets:
            jobs.extend(_jobs_cami(CAMI_MARINE_SHORT_FASTA, CAMI_MARINE_SHORT_TRUTH, marine_samples, "CAMIMARINE_SHORT"))
        if "cami-mouse-long" in datasets:
            jobs.extend(_jobs_cami(CAMI_MOUSE_LONG_FASTA, CAMI_MOUSE_LONG_TRUTH, mouse_samples, "CAMIMOUSE_LONG"))
        if "cami-mouse-short" in datasets:
            jobs.extend(_jobs_cami(CAMI_MOUSE_SHORT_FASTA, CAMI_MOUSE_SHORT_TRUTH, mouse_samples, "CAMIMOUSE_SHORT"))
        if "simulated" in datasets:
            jobs.extend(_jobs_simulated(args.simulated_limit))

    total = len(db_paths) * (len(batch_jobs) if batch_mode == "group" else len(jobs))
    if total == 0:
        print("[skip] no jobs to run (datasets or paths missing).")
        return 0

    idx = 0
    for db_path in db_paths:
        if batch_mode == "group":
            for job in batch_jobs:
                idx += 1
                missing = False
                for p in job.reads:
                    if not p.is_file():
                        print(f"[skip] missing reads: {p}")
                        missing = True
                for p in job.truths:
                    if not p.is_file():
                        print(f"[skip] missing truth: {p}")
                        missing = True
                if missing:
                    continue

                name = _safe_name(f"{job.key}__{db_path.stem}")
                print(f"[{idx}/{total}] {name} (batch inputs={len(job.reads)})")

                ns = argparse.Namespace(
                    input=[str(p) for p in job.reads],
                    truth=[str(p) for p in job.truths],
                    db=str(db_path),
                    name=name,
                    tag=args.tag,
                    threads=args.threads,
                    extra=args.extra,
                )
                rc = cmd_pipeline(ns)
                if rc != 0:
                    print(f"[warn] job failed (rc={rc}): {name}")
        else:
            for job in jobs:
                idx += 1
                if not job.reads.is_file():
                    print(f"[skip] missing reads: {job.reads}")
                    continue
                if not job.truth.is_file():
                    print(f"[skip] missing truth: {job.truth}")
                    continue

                name = _safe_name(f"{job.key}__{db_path.stem}")
                print(f"[{idx}/{total}] {name}")

                ns = argparse.Namespace(
                    input=[str(job.reads)],
                    truth=[str(job.truth)],
                    db=str(db_path),
                    name=name,
                    tag=args.tag,
                    threads=args.threads,
                    extra=args.extra,
                )
                rc = cmd_pipeline(ns)
                if rc != 0:
                    print(f"[warn] job failed (rc={rc}): {name}")

    return 0


def cmd_eval(args: argparse.Namespace) -> int:
    script = Path(os.environ.get("EVAL_SCRIPT", str(DEFAULT_EVAL_SCRIPT))).resolve()
    if not script.is_file():
        raise FileNotFoundError(f"Eval script not found: {script}")

    pred_path = Path(args.pred).resolve()
    truth_path = Path(args.truth).resolve()
    out_path = Path(args.out).resolve()

    if not pred_path.is_file():
        raise FileNotFoundError(f"--pred not found: {pred_path}")
    if not truth_path.is_file():
        raise FileNotFoundError(f"--truth not found: {truth_path}")

    cmd = [
        sys.executable,
        str(script),
        "--pred",
        str(pred_path),
        "--truth",
        str(truth_path),
        "--out",
        str(out_path),
    ]
    if args.extra:
        cmd.extend(args.extra)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)

    if not out_path.is_file():
        raise RuntimeError(f"Eval finished but metrics file not found: {out_path}")

    # Print a concise summary for quick inspection.
    print(f"[OK] eval -> {out_path}")
    summary = _summarize_metrics(out_path)
    if summary:
        print(summary)

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bench/bench.py", description="Benchmark driver for TAXICF.")

    parser.add_argument(
        "--version",
        action="version",
        version="TAXICF bench (repo-local)",
    )

    sub = parser.add_subparsers(dest="cmd", required=True)

    p_build = sub.add_parser("build", help="Run `TAXICF build` and record logs/metrics.")
    p_build.add_argument("-i", "--input", required=True, help="Input target TSV (path -> taxid).")
    p_build.add_argument("--db", required=True, help="DB name (under DB_ROOT) or a path.")
    p_build.add_argument("--name", default=None, help="Run name (default: db stem).")
    p_build.add_argument("--tag", default="default", help="Tag appended to run directory.")
    p_build.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help=f"Threads (default: env THREADS or {DEFAULT_THREADS}).",
    )
    p_build.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Extra args passed to TAXICF after '--extra'. Example: --extra --load-factor 0.9",
    )
    p_build.set_defaults(func=cmd_build)

    p_classify = sub.add_parser(
        "classify", help="Run `TAXICF classify` and record logs/metrics."
    )
    p_classify.add_argument("-i", "--input", required=True, help="Input reads (FASTA/FASTQ, optionally compressed).")
    p_classify.add_argument("--db", required=True, help="DB name (under DB_ROOT) or a path.")
    p_classify.add_argument("--name", default=None, help="Run name (default: input stem).")
    p_classify.add_argument("--tag", default="default", help="Tag appended to run directory.")
    p_classify.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help=f"Threads (default: env THREADS or {DEFAULT_THREADS}).",
    )
    p_classify.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Extra args passed to TAXICF after '--extra'. Example: --extra --post-thres 0.6",
    )
    p_classify.set_defaults(func=cmd_classify)

    p_eval = sub.add_parser("eval", help="Evaluate a classify TSV against truth (TAXICF-style metrics).")
    p_eval.add_argument("--pred", required=True, help="Prediction TSV (TAXICF classify output).")
    p_eval.add_argument("--truth", required=True, help="Truth mapping TSV (e.g., CAMI gsa_mapping.tsv).")
    p_eval.add_argument("--out", required=True, help="Output metrics JSON path.")
    p_eval.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Extra args passed to eval script after '--extra'. Example: --extra --primary-rank genus",
    )
    p_eval.set_defaults(func=cmd_eval)

    p_pipe = sub.add_parser("pipeline", help="Run `TAXICF classify` + eval and record logs/metrics.")
    p_pipe.add_argument(
        "-i",
        "--input",
        action="append",
        required=True,
        help="Input reads (FASTA/FASTQ). Repeat -i for multiple files (DB loaded once).",
    )
    p_pipe.add_argument(
        "--truth",
        action="append",
        required=True,
        help="Truth mapping TSV. Either one merged file, or repeat --truth to match each -i.",
    )
    p_pipe.add_argument("--db", required=True, help="DB name (under DB_ROOT) or a path.")
    p_pipe.add_argument("--name", default=None, help="Run name (default: <input-stem>__<db-stem>).")
    p_pipe.add_argument("--tag", default="default", help="Tag appended to run directory.")
    p_pipe.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help=f"Threads (default: env THREADS or {DEFAULT_THREADS}).",
    )
    p_pipe.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Extra args passed to TAXICF classify after '--extra'. Example: --extra --post-thres 0.6",
    )
    p_pipe.set_defaults(func=cmd_pipeline)

    p_all = sub.add_parser("run-all", help="Run pipeline for all datasets (completeONE/complete + CAMI/MOUSE/SIMULATED).")
    p_all.add_argument(
        "--dbs",
        nargs="*",
        default=[],
        help="DB names (under DB_ROOT) or paths. Default: completeONE_v2_tag8/tag16 + complete_v2_tag8/tag16.",
    )
    p_all.add_argument(
        "--datasets",
        nargs="*",
        default=[],
        choices=["cami-marine-long", "cami-marine-short", "cami-mouse-long", "cami-mouse-short", "simulated"],
        help="Which dataset groups to run (default: all).",
    )
    p_all.add_argument("--marine-samples", nargs="*", type=int, default=[], help="CAMI marine sample IDs (default: 0..9).")
    p_all.add_argument("--mouse-samples", nargs="*", type=int, default=[], help="CAMI mouse sample IDs (default: 0..63).")
    p_all.add_argument("--simulated-limit", type=int, default=None, help="Limit number of SIMULATED fasta files (default: all).")
    p_all.add_argument(
        "--batch-mode",
        choices=["group", "file"],
        default="group",
        help="Batch mode: group=one classify per dataset group (fast, DB loads once); file=per file (slow).",
    )
    p_all.add_argument("--tag", default="all", help="Tag appended to run directory.")
    p_all.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help=f"Threads (default: env THREADS or {DEFAULT_THREADS}).",
    )
    p_all.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Extra args passed to TAXICF classify after '--extra' (applies to ALL jobs).",
    )
    p_all.set_defaults(func=cmd_run_all)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
