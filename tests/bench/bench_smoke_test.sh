#!/usr/bin/env bash
set -euo pipefail

python3 bench/bench.py --help >/dev/null
python3 bench/bench.py build --help >/dev/null
python3 bench/bench.py classify --help >/dev/null
python3 bench/bench.py eval --help >/dev/null
python3 bench/bench.py pipeline --help >/dev/null
python3 bench/bench.py run-all --help >/dev/null

# Default threads should be 32 unless overridden.
python3 - <<'PY'
import os
import runpy

os.environ.pop("THREADS", None)
mod = runpy.run_path("bench/bench.py")
cfg = mod["BenchConfig"].from_env()
assert cfg.threads == 32, f"expected default THREADS=32, got {cfg.threads}"
PY

THREADS=12 python3 - <<'PY'
import runpy

mod = runpy.run_path("bench/bench.py")
cfg = mod["BenchConfig"].from_env()
assert cfg.threads == 12, f"expected THREADS=12 override, got {cfg.threads}"
PY

# Basic eval sanity check (uses local ete3 taxonomy DB if available).
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

cat >"$tmpdir/truth.tsv" <<'TSV'
#anonymous_contig_id	genome_id	tax_id	contig_id	number_reads	start_position	end_position
readA	genomeA	9606	contigA	1	0	0
readB	genomeB	562	contigB	1	0	0
TSV

cat >"$tmpdir/pred.tsv" <<'TSV'
readA	9606:1.0
readB	unclassified
TSV

python3 bench/bench.py eval --pred "$tmpdir/pred.tsv" --truth "$tmpdir/truth.tsv" --out "$tmpdir/metrics.json" >/dev/null

METRICS_JSON="$tmpdir/metrics.json" python3 - <<'PY'
import json
import os
from pathlib import Path

metrics = json.loads(Path(os.environ["METRICS_JSON"]).read_text())
assert metrics["desc_tp"] == 1, metrics
assert metrics["desc_fp"] == 0, metrics
assert metrics["desc_fn"] == 1, metrics

# Precision/recall should be present (paper-style reporting).
assert "desc" in metrics and isinstance(metrics["desc"], dict), metrics
assert metrics["desc"]["tp"] == metrics["desc_tp"], metrics
assert metrics["desc"]["fp"] == metrics["desc_fp"], metrics
assert metrics["desc"]["fn"] == metrics["desc_fn"], metrics
assert abs(metrics["desc"]["precision"] - 1.0) < 1e-9, metrics
assert abs(metrics["desc"]["recall"] - 0.5) < 1e-9, metrics

assert "exact" in metrics and isinstance(metrics["exact"], dict), metrics
assert metrics["exact"]["rank"] in ("species", "genus"), metrics
assert "precision" in metrics["exact"] and "recall" in metrics["exact"] and "f1" in metrics["exact"], metrics
PY

# Headerless truth format (SIMULATED-style): <read_id> <ignored> <taxid>
cat >"$tmpdir/truth_no_header.tsv" <<'TSV'
readA	0	9606
readB	0	562
TSV
python3 bench/bench.py eval --pred "$tmpdir/pred.tsv" --truth "$tmpdir/truth_no_header.tsv" --out "$tmpdir/metrics_no_header.json" >/dev/null

METRICS_JSON="$tmpdir/metrics_no_header.json" python3 - <<'PY'
import json
import os
from pathlib import Path

metrics = json.loads(Path(os.environ["METRICS_JSON"]).read_text())
assert metrics["desc_tp"] == 1, metrics
assert metrics["desc_fp"] == 0, metrics
assert metrics["desc_fn"] == 1, metrics
PY

echo "bench smoke test: OK"
