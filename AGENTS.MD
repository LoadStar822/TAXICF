# TAXICF Agent Guide

## 项目定位
- **软件名**：TAXICF
- **版本**：0.1.0
- **核心数据结构**：仅使用 **ICF（Interleaved Cuckoo Filter）**。

## 重要约束（必须遵守）
- **禁止提交/上传**任何私有材料或“AI 相关”的内部内容。
  - `docs/` 目录被视为私有资料：不要提交（已在 `.gitignore` 中忽略）。
- **不要提交大体积产物**：
  - 构建数据库（`*.icf`）、分类输出（`*.tsv`）、临时目录与 benchmark 运行结果都必须留在本地（已在 `.gitignore` 中忽略）。
- **除非用户明确要求**，不要新增第三方依赖；优先删除未使用依赖而不是继续堆叠。

## 代码结构速览
- `TAXICF.cpp`：CLI 入口（子命令：`build` / `classify`）
- `src/build/`：数据库构建逻辑（ICF）
- `src/classify/`：分类逻辑
- `src/utils/`：syncmer、EM/VEM 等工具模块
- `bench/bench.py`：benchmark 入口脚本（用于落盘命令、日志、metrics；不属于核心运行路径）
- `tests/`：轻量测试（ICF / syncmer / bench 冒烟）

## 构建与测试（本地可复现）
### 编译
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### 运行测试
```bash
./build/icf_filter_tests
./build/syncmer_tests
bash tests/bench/bench_smoke_test.sh
```

## Benchmark 约定（参考 TAXICF/bench 的思路）
- 统一入口：`bench/bench.py`
  - build：记录 `TAXICF build` 的运行日志与指标
  - classify：记录 `TAXICF classify` 的运行日志与指标
- 产物目录（默认）：
  - 数据库：`bench/DB/`
  - 运行日志：`bench/runs/`
- 可通过环境变量覆盖：
  - `TAXICF_BIN`（默认：`build/TAXICF`）
  - `DB_ROOT`（默认：`bench/DB`）
  - `RUN_ROOT`（默认：`bench/runs`）
  - `THREADS`（默认：32）

## 数据集路径（本机只读）
> 这些路径用于本机 bench/build/classify 的复现实验；不要把任何数据集/数据库/结果提交到 GitHub。

### Build 输入 target.tsv
- RefSeq complete（路径为**绝对路径**，可在任意工作目录使用）：
  - `/mnt/sda/CommonData/TaxonomyDataset/complete/target.tsv`
  - （可选分组版本）`/mnt/sda/CommonData/TaxonomyDataset/complete/target_group.tsv`
- RefSeq completeONE（`target.tsv` 内是**相对路径**，需要注意工作目录）：
  - `/mnt/sda/CommonData/TaxonomyDataset/completeONE/target.tsv`
  - 其中每行类似 `completeONE/.../GCF_*.fna.gz <taxid>`，相对路径以 `/mnt/sda/CommonData/TaxonomyDataset` 为基准。

### CAMI 海洋数据集（用于测试 classify）
- 长数据集（10 个样本，id=0..9）：
  - 根目录：`/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMARINE_LONG`
  - 输入 FASTA：`.../fasta/marmgCAMI2_sample_<id>_contigs_anonymous_gsa.fasta`
  - 金标准 mapping：`.../mapping/marmgCAMI2_sample_<id>_contigs_gsa_mapping.tsv`
- 短数据集（10 个样本，id=0..9）：
  - 根目录：`/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMARINE_SHORT`
  - 输入/金标准（已解压在）：`.../output_dir/`
  - 输入 FASTA：`.../output_dir/marmgCAMI2_sample_<id>_contigs_simulation_short_read_2018.08.15_09.49.32_sample_<id>_contigs_anonymous_gsa.fasta`
  - 金标准 mapping：`.../output_dir/marmgCAMI2_sample_<id>_contigs_simulation_short_read_2018.08.15_09.49.32_sample_<id>_contigs_gsa_mapping.tsv`

### CAMI 老鼠数据集（用于测试 classify）
- 长数据集（64 个样本，id=0..63）：
  - 根目录：`/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMOUSE_LONG`
  - 输入/金标准（已解压在）：`.../output_dir/`
  - 输入 FASTA：`.../output_dir/2018.02.13_14.02.01_sample_<id>_contigs_2018.02.13_14.02.01_sample_<id>_contigs_anonymous_gsa.fasta`
  - 金标准 mapping：`.../output_dir/2018.02.13_14.02.01_sample_<id>_contigs_2018.02.13_14.02.01_sample_<id>_contigs_gsa_mapping.tsv`
- 短数据集（64 个样本，id=0..63）：
  - 根目录：`/mnt/sda/tianqinzhong/project/taxicfBench/dataset/CAMIMOUSE_SHORT`
  - 输入/金标准（已解压在）：`.../output_dir/`
  - 输入 FASTA：`.../output_dir/2017.12.29_11.37.26_sample_<id>_contigs_2017.12.29_11.37.26_sample_<id>_contigs_anonymous_gsa.fasta`
  - 金标准 mapping（推荐用这个：TSV header 正常）：`.../output_dir/2017.12.29_11.37.26_sample_<id>_contigs_2017.12.29_11.37.26_sample_<id>_contigs_gsa_mapping.tsv`
  - 备注：同目录下存在 `*_gsa_mapping_new.tsv`，但其 header 是空格分隔（不是 tab），当前 `bench/scripts/eval_classify.py` 会解析失败。

### SIMULATED 模拟数据集（用于测试 classify）
- 根目录：`/mnt/sda/tianqinzhong/project/taxicfBench/dataset/SIMULATED`
- 特点：每个 `.fasta` 文件只包含**一个物种**的 reads，物种名可从文件名开头解析（例如 `Aeropyrum_pernix_... .fasta` 对应物种 `Aeropyrum_pernix`）。
- “金标准”需要特殊处理：
  - 目录下通常有同名 `.tsv`（无 header），其第 1 列为 read id，第 3 列为 taxid（整文件同一 taxid）。
  - 或者在 bench/eval 中直接把同一 `.fasta` 内的所有 reads 视为同一物种（taxid 由文件名映射得到）。

## 改动原则（面向长期维护）
- 优先做“删减/收敛”而不是扩张：只保留 ICF + build/classify 的清晰边界。
- 任何影响 CLI 行为/输出格式的修改，都要同步更新 `README.md`（并保持示例可跑）。
- 修改后至少确保：
  - `cmake --build build -j` 可编译通过；
  - `./build/icf_filter_tests` 与 `./build/syncmer_tests` 通过；
  - `tests/bench/bench_smoke_test.sh` 通过。
