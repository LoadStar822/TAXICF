# TAXICF

TAXICF is a lightweight metagenomic sequence classification tool.

Current version: **0.1.0**

Core idea: build an **ICF (Interleaved Cuckoo Filter)** database from per-taxid syncmer hash sets derived from reference sequences; during classification, compute syncmer hashes for each read and count hits in the ICF to produce candidate taxids.

---

## Features

- `build`: build an ICF database from an input manifest (outputs `*.icf`)
- `classify`: classify FASTA/FASTQ reads (compressed input supported) (outputs `*.tsv`)

---

## Build

### Dependencies

- CMake ≥ 3.10
- A C++20-capable compiler (GCC / Clang)
- zlib, bzip2 (for compressed inputs; typically required by SeqAn3)
- OpenMP (optional, to accelerate parts of the pipeline)

Ubuntu example:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git libbz2-dev zlib1g-dev
```

### Compile

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Built executable: `build/TAXICF`

---

## Usage

### 1) Build database (`build`)

Prepare a two-column TSV manifest (one reference file path and its taxid per line):

```
/path/to/ref1.fna.gz    562
/path/to/ref2.fasta     1280
```

Run:

```bash
./build/TAXICF build -i target.tsv -o TAXICFDB
```

Output database file: `TAXICFDB.icf`

#### DB size optimization (recommended)

By default, `build` uses the **v2 multi-block ICF** format: taxids are grouped into multiple ICF blocks by the `binSize` they need (still ICF-only). This avoids the large database blow-up caused by a global `maxHashes` forcing an oversized `binSize` for everyone.

You can control the grouping granularity with `--block-factor` (smaller = finer groups and smaller DB, but more blocks; start with `1.5`):

```bash
./build/TAXICF build -i target.tsv -o TAXICFDB --block-factor 1.5
```

To force the legacy single-block v1 format (usually larger, but useful for comparison), set it to `<= 0`:

```bash
./build/TAXICF build -i target.tsv -o TAXICFDB --block-factor 0
```

### 2) Classify reads (`classify`)

Single-end:

```bash
./build/TAXICF classify -i reads.fq -d TAXICFDB -o TAXICFClassify.tsv
```

Paired-end (paired files; must pass an even number of `-p` arguments):

```bash
./build/TAXICF classify -p R1.fq.gz -p R2.fq.gz -d TAXICFDB -o TAXICFClassify.tsv
```

Notes:
- `-d/--database` accepts `TAXICFDB` or `TAXICFDB.icf` (the program will append `.icf` automatically)
- Output is TSV; by default it will have a `.tsv` suffix
  - `classify` supports both v1 and v2 databases

---

## Temporary files

During `build`, TAXICF creates a temporary directory in the **current working directory** (e.g., `__taxicf_tmp_<output>`) to store intermediate artifacts, and removes it automatically on successful completion.

It is recommended to run builds under the benchmark run directory (this repo’s `bench/bench.py` is organized that way).

---

## Tests (optional)

```bash
./build/icf_filter_tests
./build/syncmer_tests
```

---

## License

MIT
