# TAXICF

TAXICF 是一个轻量的宏基因组序列分类工具（metagenomic classifier）。

当前版本：**0.1.0**

核心思路：对参考序列按 taxid 生成 syncmer 哈希集合，并构建 **ICF（Interleaved Cuckoo Filter）** 数据库；分类时对每条 read 计算 syncmer 哈希并在 ICF 中统计命中，输出候选 taxid。

---

## 功能

- `build`：从输入清单构建 ICF 数据库（输出 `*.icf`）
- `classify`：对 FASTA/FASTQ（支持压缩）进行分类（输出 `*.tsv`）

---

## 构建

### 依赖

- CMake ≥ 3.10
- 支持 C++20 的编译器（GCC / Clang）
- zlib、bzip2（用于压缩输入，通常由 SeqAn3 依赖）
- OpenMP（可选，用于加速部分流程）

Ubuntu 示例：

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git libbz2-dev zlib1g-dev
```

### 编译

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

生成可执行文件：`build/TAXICF`

---

## 使用

### 1) 构建数据库（build）

准备一个两列的 TSV 清单（每行一个参考文件与其 taxid）：

```
/path/to/ref1.fna.gz    562
/path/to/ref2.fasta     1280
```

运行：

```bash
./build/TAXICF build -i target.tsv -o TAXICFDB
```

输出数据库文件：`TAXICFDB.icf`

#### DB 体积优化（推荐）

默认情况下，`build` 会使用 **v2 多块 ICF** 格式：把不同 taxid 按各自需要的 `binSize` 分组到多个 ICF block 里（仍然只使用 ICF），从根上避免“全局 maxHashes 绑死 binSize”导致的超大数据库。

你可以用 `--block-factor` 控制分组粒度（越小越细、越省空间，但 block 数会变多；建议从 `1.5` 开始）：

```bash
./build/TAXICF build -i target.tsv -o TAXICFDB --block-factor 1.5
```

如果你想强制输出旧的单块 v1 格式（通常会更大，但便于对比），把它设为 `<=0`：

```bash
./build/TAXICF build -i target.tsv -o TAXICFDB --block-factor 0
```

### 2) 序列分类（classify）

单端：

```bash
./build/TAXICF classify -i reads.fq -d TAXICFDB -o TAXICFClassify.tsv
```

双端（成对文件，必须偶数个参数）：

```bash
./build/TAXICF classify -p R1.fq.gz -p R2.fq.gz -d TAXICFDB -o TAXICFClassify.tsv
```

说明：
- `-d/--database` 可以写 `TAXICFDB` 或 `TAXICFDB.icf`（程序会自动补全 `.icf`）
- 输出为 TSV；默认会带 `.tsv` 后缀
  - `classify` 同时支持读取 v1 / v2 数据库

---

## 临时文件说明

`build` 阶段会在**当前工作目录**创建一个临时目录（形如 `__taxicf_tmp_<output>`）用于存放中间产物，并在成功结束后自动删除。
建议在 benchmark 脚本的 run 目录下执行 build（本仓库的 `bench/bench.py` 已按这个方式组织）。

---

## 测试（可选）

```bash
./build/icf_filter_tests
./build/syncmer_tests
```

---

## License

MIT
