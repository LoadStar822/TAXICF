/*
 * -----------------------------------------------------------------------------
 * Filename:      TAXICFBuild.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-30
 *
 * Last Modified: 2025-12-30
 *
 * Description:
 *  TAXICF build module (ICF-only).
 *
 * -----------------------------------------------------------------------------
 */
#include <TAXICFBuild.hpp>

#include <interleaved-cuckoo-filter.h>
#include <icf_db.hpp>
#include <robin_hood.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <dna4_traits.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <span>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <utils/Syncmer.hpp>

namespace TAXICFBuild {
namespace {
static inline bool file_is_compressed(const std::filesystem::path &filepath) {
  std::string extension = filepath.extension().string();
  return extension == ".gz" || extension == ".bgzf" || extension == ".bz2";
}

void print_build_time(long long milliseconds) {
  long long total_seconds = milliseconds / 1000;
  long long seconds = total_seconds % 60;
  long long total_minutes = total_seconds / 60;
  long long minutes = total_minutes % 60;
  long long hours = total_minutes / 60;

  if (hours > 0) {
    std::cout << hours << "h " << minutes << "min " << seconds << "s "
              << milliseconds % 1000 << "ms" << std::endl;
  } else if (minutes > 0) {
    std::cout << minutes << "min " << seconds << "s " << milliseconds % 1000
              << "ms" << std::endl;
  } else {
    std::cout << seconds << "s " << milliseconds % 1000 << "ms" << std::endl;
  }
}

void create_or_reset_directory(const std::string &dir,
                               const BuildConfig &config) {
  std::filesystem::path directoryPath = dir;
  if (std::filesystem::exists(directoryPath)) {
    if (!std::filesystem::is_directory(directoryPath)) {
      throw std::runtime_error("'" + dir +
                               "' exists but is not a directory.");
    }
    std::filesystem::remove_all(directoryPath);
    if (config.verbose) {
      std::cout << "Directory '" << dir << "' existed and was removed."
                << std::endl;
    }
  }
  std::filesystem::create_directory(directoryPath);
  if (config.verbose) {
    std::cout << "Directory '" << dir << "' created successfully." << std::endl;
  }
}

void parse_input_file(
    const std::string &filePath,
    robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    FileInfo &fileInfo) {
  std::ifstream inputFile(filePath);
  if (!inputFile.is_open()) {
    throw std::runtime_error("Failed to open input file: " + filePath);
  }
  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    std::string seqPath;
    std::string taxidStr;
    if (!(iss >> seqPath >> taxidStr)) {
      fileInfo.invalidNum++;
      continue;
    }
    hashCount[taxidStr] = 0;
    inputFiles[taxidStr].push_back(seqPath);
    fileInfo.fileNum++;
  }
}

void syncmer_count(
    BuildConfig &config,
    robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles,
    robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount,
    FileInfo &fileInfo,
    const std::string &tmpDir) {
  const std::array<size_t, 1> syncmer_positions{
      static_cast<size_t>(config.syncmer_position)};
  const std::span<const size_t> syncmer_pos_span(syncmer_positions);
  const uint64_t syncmer_seed = TAXICFBuild::adjust_seed(config.kmer_size);

  struct TaxidTask {
    std::string taxid;
    const std::vector<std::string> *files;
  };

  std::vector<TaxidTask> tasks;
  tasks.reserve(inputFiles.size());
  for (auto &[taxid, files] : inputFiles) {
    // Sorting within each taxid improves sequential access patterns on HDDs.
    std::sort(files.begin(), files.end());
    tasks.push_back(TaxidTask{taxid, &files});
  }

  std::vector<uint64_t> taxidHashCounts(tasks.size(), 0);
  std::mutex fileInfo_mutex;

#pragma omp parallel for schedule(dynamic)
  for (size_t idx = 0; idx < tasks.size(); ++idx) {
    const auto &task = tasks[idx];
    const std::string &taxid = task.taxid;
    const std::vector<std::string> &taxid_files = *task.files;

    FileInfo localFileInfo{};
    size_t max_hashes = config.max_hashes_per_taxid;
    uint64_t current_count = 0;

    robin_hood::unordered_flat_map<uint64_t, uint8_t> local_hash_counts;
    std::vector<uint64_t> hashes;
    hashes.reserve(4096);
    std::vector<uint64_t> filtered_hashes;
    filtered_hashes.reserve(4096);

    std::ofstream ofile;
    std::string output_filename = tmpDir + "/" + taxid + ".sync";
    ofile.open(output_filename, std::ios::binary | std::ios::app);
    if (!ofile.is_open()) {
      std::lock_guard<std::mutex> g(fileInfo_mutex);
      std::cerr << "Unable to open the syncmer file: " << output_filename
                << std::endl;
    }

    for (const auto &filename : taxid_files) {
      if (max_hashes > 0 && current_count >= max_hashes) {
        localFileInfo.skippedNum++;
        continue;
      }

      uint8_t cutoff = 0;
      if (config.adaptive_cutoff) {
        std::filesystem::path filepath(filename);
        size_t filesize = std::filesystem::file_size(filepath);
        bool is_compressed = file_is_compressed(filepath);
        size_t adjusted_filesize = filesize * 2 / (is_compressed ? 1 : 3);

        if (adjusted_filesize <= 314'572'800ULL)
          cutoff = 1;
        else if (adjusted_filesize <= 524'288'000ULL)
          cutoff = 3;
        else if (adjusted_filesize <= 1'073'741'824ULL)
          cutoff = 10;
        else if (adjusted_filesize <= 3'221'225'472ULL)
          cutoff = 20;
        else
          cutoff = 50;
      }

      local_hash_counts.clear();

      seqan3::sequence_file_input<raptor::dna4_traits,
                                 seqan3::fields<seqan3::field::id,
                                               seqan3::field::seq>>
          fin{filename};
      for (auto &record : fin) {
        auto &seq = record.sequence();
        if (seq.size() < config.min_length) {
          localFileInfo.skippedSeqNum++;
          continue;
        }
        localFileInfo.sequenceNum++;
        localFileInfo.bpLength += seq.size();

        hashes.clear();
        taxicf::syncmer::compute_hashes_append(
            seq, config.smer_size, config.kmer_size, syncmer_pos_span,
            syncmer_seed, true, hashes);
        for (uint64_t hash : hashes) {
          uint8_t &count = local_hash_counts[hash];
          if (count < 255) {
            ++count;
          }
        }
      }

      filtered_hashes.clear();
      filtered_hashes.reserve(local_hash_counts.size());
      for (const auto &[hash, count] : local_hash_counts) {
        if (count >= cutoff) {
          filtered_hashes.emplace_back(hash);
        }
      }
      local_hash_counts.clear();

      if (max_hashes > 0) {
        size_t remaining = static_cast<size_t>(max_hashes - current_count);
        if (filtered_hashes.size() > remaining) {
          filtered_hashes.resize(remaining);
        }
      }

      if (!filtered_hashes.empty()) {
        if (ofile.is_open()) {
          ofile.write(reinterpret_cast<const char *>(filtered_hashes.data()),
                      static_cast<std::streamsize>(filtered_hashes.size() *
                                                   sizeof(uint64_t)));
          if (!ofile.good()) {
            std::lock_guard<std::mutex> g(fileInfo_mutex);
            std::cerr << "Failed to write syncmer hashes: " << output_filename
                      << std::endl;
          }
        }
        current_count += filtered_hashes.size();
      }
    }

    taxidHashCounts[idx] = current_count;

    {
      std::lock_guard<std::mutex> lock(fileInfo_mutex);
      fileInfo.skippedNum += localFileInfo.skippedNum;
      fileInfo.skippedSeqNum += localFileInfo.skippedSeqNum;
      fileInfo.sequenceNum += localFileInfo.sequenceNum;
      fileInfo.bpLength += localFileInfo.bpLength;
    }
  }

  // Update hashCount map sequentially (avoid concurrent writes to the hash table).
  for (size_t idx = 0; idx < tasks.size(); ++idx) {
    hashCount[tasks[idx].taxid] = taxidHashCounts[idx];
  }
}

uint64_t get_max_value(
    const robin_hood::unordered_flat_map<std::string, uint64_t> &hashCount) {
  uint64_t maxValue = 0;
  for (const auto &kv : hashCount) {
    maxValue = std::max(maxValue, kv.second);
  }
  return maxValue;
}

std::vector<std::string> sorted_taxids(
    const robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
        &inputFiles) {
  std::vector<std::string> taxids;
  taxids.reserve(inputFiles.size());
  for (const auto &kv : inputFiles) {
    taxids.push_back(kv.first);
  }
  std::sort(taxids.begin(), taxids.end());
  return taxids;
}

void save_icf(const taxicf::InterleavedCuckooFilter &icf,
              const std::string &output_base,
              const std::vector<std::string> &indexToTaxid,
              const ICFConfig &config) {
  const std::string outPath = output_base + ".icf";
  std::ofstream os(outPath, std::ios::binary);
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open file: " + outPath);
  }

  cereal::BinaryOutputArchive archive(os);
  archive(icf);
  archive(indexToTaxid);
  archive(config);
}

size_t round_up_geometric(size_t value, double factor) {
  if (value <= 1) {
    return 1;
  }
  if (!(factor > 1.0)) {
    return value;
  }
  size_t cur = 1;
  while (cur < value) {
    const double next_d = std::ceil(static_cast<double>(cur) * factor);
    size_t next = static_cast<size_t>(next_d);
    if (next <= cur) {
      next = cur + 1;
    }
    cur = next;
  }
  return cur;
}
} // namespace

void run(BuildConfig config) {
  if (config.threads == 0) {
    unsigned int hardwareThreads = std::thread::hardware_concurrency();
    if (hardwareThreads == 0) {
      hardwareThreads = 1;
    }
    const auto maxThreads =
        static_cast<unsigned int>(std::numeric_limits<uint16_t>::max());
    if (hardwareThreads > maxThreads) {
      hardwareThreads = maxThreads;
    }
    config.threads = static_cast<uint16_t>(hardwareThreads);
  }

  if (config.smer_size == 0) {
    throw std::runtime_error("Syncmer s-mer size must be greater than 0.");
  }
  if (config.smer_size >= config.kmer_size) {
    throw std::runtime_error("Syncmer s-mer size must be smaller than k-mer.");
  }
  const uint16_t syncmer_span =
      static_cast<uint16_t>(config.kmer_size - config.smer_size + 1);
  if (config.syncmer_position >= syncmer_span) {
    throw std::runtime_error("Syncmer offset must satisfy 0 <= pos < k - s + 1.");
  }
  if (!(config.load_factor > 0.0 && config.load_factor < 1.0)) {
    throw std::runtime_error("--load-factor must be in (0, 1).");
  }

  if (config.verbose) {
    std::cout << config << std::endl;
  }

#ifdef _OPENMP
  omp_set_num_threads(config.threads);
#endif

  auto build_start = std::chrono::high_resolution_clock::now();

  std::cout << "Reading input files..." << std::endl;
  FileInfo fileInfo{};
  robin_hood::unordered_flat_map<std::string, uint64_t> hashCount;
  robin_hood::unordered_flat_map<std::string, std::vector<std::string>>
      inputFiles;
  parse_input_file(config.input_file, inputFiles, hashCount, fileInfo);

  const std::string tmpDir =
      "__taxicf_tmp_" + std::filesystem::path(config.output_file).filename().string();
  create_or_reset_directory(tmpDir, config);

  auto calculate_start = std::chrono::high_resolution_clock::now();
  std::cout << "Calculating syncmers..." << std::endl;
  syncmer_count(config, inputFiles, hashCount, fileInfo, tmpDir);
  auto calculate_end = std::chrono::high_resolution_clock::now();
  if (config.verbose) {
    std::cout << "Calculate time: ";
    print_build_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                         calculate_end - calculate_start)
                         .count());
  }

  std::vector<std::string> indexToTaxid = sorted_taxids(inputFiles);
  if (indexToTaxid.empty()) {
    throw std::runtime_error("No taxids found in input file.");
  }

  uint64_t maxHashes = get_max_value(hashCount);
  constexpr size_t kTagNum = 4;
  const double denom = config.load_factor * static_cast<double>(kTagNum);

  const bool use_multiblock = config.block_factor > 0.0;
  if (!use_multiblock) {
    size_t binNum = indexToTaxid.size();
    size_t binSize = 1;
    if (maxHashes > 0) {
      binSize = std::max<size_t>(
          1, static_cast<size_t>(
                 std::ceil(static_cast<double>(maxHashes) / denom)));
    }

    std::cout << "Building ICF... (bins=" << binNum << ", binSize=" << binSize
              << ", tagBits=" << static_cast<int>(config.tag_bits)
              << ", maxHashesPerBin=" << maxHashes << ")" << std::endl;

    taxicf::InterleavedCuckooFilter icf(binNum, binSize, config.tag_bits);

    auto insert_one_taxid = [&](size_t binIndex) {
      const std::string &taxid = indexToTaxid[binIndex];
      std::ifstream ifile(tmpDir + "/" + taxid + ".sync", std::ios::binary);
      if (!ifile.is_open()) {
        return; // taxid may have 0 hashes
      }

      try {
        std::vector<uint64_t> buffer(128 * 1024);
        while (ifile) {
          ifile.read(reinterpret_cast<char *>(buffer.data()),
                     static_cast<std::streamsize>(buffer.size() *
                                                  sizeof(uint64_t)));
          std::streamsize bytesRead = ifile.gcount();
          if (bytesRead <= 0) {
            break;
          }
          size_t count = static_cast<size_t>(bytesRead) / sizeof(uint64_t);
          for (size_t i = 0; i < count; ++i) {
            icf.insertTag(binIndex, buffer[i]);
          }
        }
      } catch (const std::exception &e) {
        throw std::runtime_error("ICF insert failed (taxid=" + taxid +
                                 ", bin=" + std::to_string(binIndex) +
                                 "): " + e.what());
      }
    };

    std::atomic<bool> failed{false};
    std::exception_ptr eptr = nullptr;

    if (config.tag_bits == 8) {
      const size_t pairs = (indexToTaxid.size() + 1) / 2;
#pragma omp parallel for schedule(dynamic)
      for (size_t pair = 0; pair < pairs; ++pair) {
        if (failed.load(std::memory_order_relaxed)) {
          continue;
        }
        const size_t bin0 = pair * 2;
        const size_t bin1 = bin0 + 1;
        try {
          insert_one_taxid(bin0);
          if (bin1 < indexToTaxid.size()) {
            insert_one_taxid(bin1);
          }
        } catch (...) {
          failed.store(true, std::memory_order_relaxed);
#pragma omp critical
          {
            if (!eptr) {
              eptr = std::current_exception();
            }
          }
        }
      }
    } else {
#pragma omp parallel for schedule(dynamic)
      for (size_t binIndex = 0; binIndex < indexToTaxid.size(); ++binIndex) {
        if (failed.load(std::memory_order_relaxed)) {
          continue;
        }
        try {
          insert_one_taxid(binIndex);
        } catch (...) {
          failed.store(true, std::memory_order_relaxed);
#pragma omp critical
          {
            if (!eptr) {
              eptr = std::current_exception();
            }
          }
        }
      }
    }

    if (eptr) {
      std::rethrow_exception(eptr);
    }

    ICFConfig icfConfig{};
    icfConfig.dbVersion = ICFConfig::CurrentDbVersion;
    icfConfig.binNum = binNum;
    icfConfig.binSize = binSize;
    icfConfig.kmerSize = config.kmer_size;
    icfConfig.smerSize = config.smer_size;
    icfConfig.syncmerPosition = config.syncmer_position;
    icfConfig.loadFactor = config.load_factor;
    icfConfig.seed64 = TAXICFBuild::adjust_seed(config.kmer_size);

    auto save_start = std::chrono::high_resolution_clock::now();
    std::cout << "Saving ICF... " << std::endl;
    save_icf(icf, config.output_file, indexToTaxid, icfConfig);
    auto save_end = std::chrono::high_resolution_clock::now();
    if (config.verbose) {
      std::cout << "Save time: ";
      print_build_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                           save_end - save_start)
                           .count());
    }
  } else {
    std::cout << "Building multi-block ICF... (bins=" << indexToTaxid.size()
              << ", tagBits=" << static_cast<int>(config.tag_bits)
              << ", maxHashesPerBin=" << maxHashes
              << ", blockFactor=" << config.block_factor << ")" << std::endl;

    robin_hood::unordered_flat_map<size_t, std::vector<std::string>> groups;
    groups.reserve(indexToTaxid.size());

    for (const auto &taxid : indexToTaxid) {
      uint64_t cnt = 0;
      auto it = hashCount.find(taxid);
      if (it != hashCount.end()) {
        cnt = it->second;
      }
      size_t needBuckets = 1;
      if (cnt > 0) {
        needBuckets = std::max<size_t>(
            1, static_cast<size_t>(
                   std::ceil(static_cast<double>(cnt) / denom)));
      }
      const size_t groupBinSize = round_up_geometric(
          needBuckets, config.block_factor);
      groups[groupBinSize].push_back(taxid);
    }

    std::vector<size_t> groupKeys;
    groupKeys.reserve(groups.size());
    for (const auto &kv : groups) {
      groupKeys.push_back(kv.first);
    }
    std::sort(groupKeys.begin(), groupKeys.end());

    IcfDbV2 db;
    db.config.kmerSize = config.kmer_size;
    db.config.smerSize = config.smer_size;
    db.config.syncmerPosition = config.syncmer_position;
    db.config.loadFactor = config.load_factor;
    db.config.seed64 = TAXICFBuild::adjust_seed(config.kmer_size);
    db.config.tagBits = config.tag_bits;

    db.blocks.reserve(groupKeys.size());
    for (size_t binSize : groupKeys) {
      const auto &taxids = groups[binSize];
      if (taxids.empty()) {
        continue;
      }
      IcfDbBlockV2 block;
      block.binSize = binSize;
      block.indexToTaxid = taxids;
      block.icf = taxicf::InterleavedCuckooFilter(taxids.size(), binSize,
                                                   config.tag_bits);
      db.blocks.push_back(std::move(block));
    }

    std::cout << "ICF blocks: " << db.blocks.size() << std::endl;
    if (config.verbose) {
      for (const auto &block : db.blocks) {
        std::cout << "  - binSize=" << block.binSize
                  << " bins=" << block.indexToTaxid.size() << std::endl;
      }
    }

    struct InsertTask {
      size_t blockIdx;
      size_t bin0;
      size_t bin1;
      bool hasBin1;
    };

    auto insert_one_taxid_in_block = [&](size_t blockIdx, size_t binIndex) {
      auto &block = db.blocks[blockIdx];
      const std::string &taxid = block.indexToTaxid[binIndex];
      std::ifstream ifile(tmpDir + "/" + taxid + ".sync", std::ios::binary);
      if (!ifile.is_open()) {
        return;
      }

      try {
        std::vector<uint64_t> buffer(128 * 1024);
        while (ifile) {
          ifile.read(reinterpret_cast<char *>(buffer.data()),
                     static_cast<std::streamsize>(buffer.size() *
                                                  sizeof(uint64_t)));
          std::streamsize bytesRead = ifile.gcount();
          if (bytesRead <= 0) {
            break;
          }
          size_t count = static_cast<size_t>(bytesRead) / sizeof(uint64_t);
          for (size_t i = 0; i < count; ++i) {
            block.icf.insertTag(binIndex, buffer[i]);
          }
        }
      } catch (const std::exception &e) {
        throw std::runtime_error("ICF insert failed (taxid=" + taxid +
                                 ", bin=" + std::to_string(binIndex) +
                                 ", binSize=" + std::to_string(block.binSize) +
                                 "): " + e.what());
      }
    };

    std::vector<InsertTask> tasks;
    tasks.reserve(indexToTaxid.size());
    for (size_t blockIdx = 0; blockIdx < db.blocks.size(); ++blockIdx) {
      const auto &taxids = db.blocks[blockIdx].indexToTaxid;
      if (config.tag_bits == 8) {
        const size_t pairs = (taxids.size() + 1) / 2;
        for (size_t pair = 0; pair < pairs; ++pair) {
          const size_t bin0 = pair * 2;
          const size_t bin1 = bin0 + 1;
          tasks.push_back(
              InsertTask{blockIdx, bin0, bin1, bin1 < taxids.size()});
        }
      } else {
        for (size_t bin = 0; bin < taxids.size(); ++bin) {
          tasks.push_back(InsertTask{blockIdx, bin, 0, false});
        }
      }
    }

    std::atomic<bool> failed{false};
    std::exception_ptr eptr = nullptr;

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < tasks.size(); ++i) {
      if (failed.load(std::memory_order_relaxed)) {
        continue;
      }
      const auto &t = tasks[i];
      try {
        insert_one_taxid_in_block(t.blockIdx, t.bin0);
        if (t.hasBin1) {
          insert_one_taxid_in_block(t.blockIdx, t.bin1);
        }
      } catch (...) {
        failed.store(true, std::memory_order_relaxed);
#pragma omp critical
        {
          if (!eptr) {
            eptr = std::current_exception();
          }
        }
      }
    }

    if (eptr) {
      std::rethrow_exception(eptr);
    }

    auto save_start = std::chrono::high_resolution_clock::now();
    std::cout << "Saving ICF... " << std::endl;
    TAXICFBuild::write_v2_db(config.output_file, db);
    auto save_end = std::chrono::high_resolution_clock::now();
    if (config.verbose) {
      std::cout << "Save time: ";
      print_build_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                           save_end - save_start)
                           .count());
    }
  }

  if (config.verbose) {
    std::cout << "Remove temporary files..." << std::endl;
  }
  std::filesystem::remove_all(tmpDir);

  auto build_end = std::chrono::high_resolution_clock::now();
  if (config.verbose) {
    std::cout << "Total build time: ";
    print_build_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                         build_end - build_start)
                         .count());
  }
}
} // namespace TAXICFBuild
