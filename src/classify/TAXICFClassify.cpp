/*
 * -----------------------------------------------------------------------------
 * Filename:      TAXICFClassify.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-08-10
 *
 * Last Modified: 2025-12-30
 *
 * Description:
 *  TAXICF classify module (ICF-only).
 * -----------------------------------------------------------------------------
 */
#include <TAXICFClassify.hpp>

#include <EM.hpp>
#include <VEM.hpp>
#include <buildConfig.hpp>
#include <dna4_traits.hpp>
#include <icf_db.hpp>
#include <interleaved-cuckoo-filter.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/contrib/stream/bgzf.hpp>

#include <utils/Syncmer.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace TAXICFClassify {
namespace {
constexpr const char *kUnclassified = "unclassified";

struct ICFDatabase {
  struct RuntimeConfig {
    uint8_t kmerSize{};
    uint16_t smerSize{};
    uint16_t syncmerPosition{0};
    double loadFactor{0.0};
    uint64_t seed64{0};
    uint8_t tagBits{8};
    uint8_t dbVersion{0};
  };

  RuntimeConfig config;

  // v1 (legacy single-block)
  taxicf::InterleavedCuckooFilter icf;
  std::vector<std::string> indexToTaxid; // bin -> taxid

  // v2 (multi-block)
  TAXICFBuild::IcfDbV2 db2;
  bool multiblock{false};
};

void print_classify_time(long long milliseconds) {
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

ICFDatabase load_icf_db(const std::string &dbFile) {
  std::filesystem::path archivePath{dbFile};
  if (!std::filesystem::exists(archivePath)) {
    archivePath = std::filesystem::path{dbFile + ".icf"};
  }
  if (!std::filesystem::exists(archivePath)) {
    throw std::runtime_error("无法找到 ICF 数据库文件: " + dbFile);
  }

  std::ifstream is(archivePath, std::ios::binary);
  if (!is.is_open()) {
    throw std::runtime_error("无法打开 ICF 数据库文件: " + archivePath.string());
  }

  ICFDatabase db;

  // Detect v2 DB by raw magic header.
  {
    uint8_t version = 0;
    if (TAXICFBuild::try_read_v2_header(is, version)) {
      if (version != TAXICFBuild::kIcfDbVersionV2) {
        throw std::runtime_error("Unsupported ICF v2 db version: " +
                                 std::to_string(version));
      }
      cereal::BinaryInputArchive archive(is);
      archive(db.db2);
      if (db.db2.config.dbVersion != TAXICFBuild::kIcfDbVersionV2) {
        throw std::runtime_error("Invalid v2 DB payload version: " +
                                 std::to_string(db.db2.config.dbVersion));
      }
      db.multiblock = true;

      db.config.kmerSize = db.db2.config.kmerSize;
      db.config.smerSize = db.db2.config.smerSize;
      db.config.syncmerPosition = db.db2.config.syncmerPosition;
      db.config.loadFactor = db.db2.config.loadFactor;
      db.config.seed64 = db.db2.config.seed64;
      db.config.tagBits = db.db2.config.tagBits;
      db.config.dbVersion = db.db2.config.dbVersion;
      return db;
    }
  }

  // Legacy v1 DB (no magic header).
  cereal::BinaryInputArchive archive(is);
  archive(db.icf);
  archive(db.indexToTaxid);
  TAXICFBuild::ICFConfig cfg_v1;
  archive(cfg_v1);

  if (cfg_v1.dbVersion != TAXICFBuild::ICFConfig::CurrentDbVersion) {
    throw std::runtime_error("ICF 数据库版本不兼容：db_version=" +
                             std::to_string(cfg_v1.dbVersion));
  }
  if (db.indexToTaxid.size() != cfg_v1.binNum) {
    throw std::runtime_error("ICF 数据库损坏：indexToTaxid.size != binNum");
  }

  db.config.kmerSize = cfg_v1.kmerSize;
  db.config.smerSize = cfg_v1.smerSize;
  db.config.syncmerPosition = cfg_v1.syncmerPosition;
  db.config.loadFactor = cfg_v1.loadFactor;
  db.config.seed64 = cfg_v1.seed64;
  db.config.dbVersion = cfg_v1.dbVersion;

  return db;
}

template <typename T>
static inline void update_minmax(T value, T &minVal, T &maxVal) {
  minVal = std::min(minVal, value);
  maxVal = std::max(maxVal, value);
}

classifyResult classify_one(const std::string &id,
                            const std::vector<uint64_t> &hashes,
                            const ClassifyConfig &config,
                            ICFDatabase &db) {
  classifyResult result;
  result.id = id;
  result.evaluated = static_cast<double>(hashes.size());

  if (hashes.empty()) {
    result.taxidCount = {{kUnclassified, 1}};
    return result;
  }

  const size_t hashNum = hashes.size();
  const size_t threshold = static_cast<size_t>(
      std::ceil(config.shotThreshold * static_cast<double>(hashNum)));

  // Keep only TopK hits to avoid keeping huge per-read vector capacities (tag8 can
  // produce many low-count hits due to higher FP rate).
  const size_t topK = config.preEmTopK;
  std::vector<std::pair<std::string, size_t>> hits;
  if (topK > 0) {
    hits.reserve(topK);
  } else {
    hits.reserve(128);
  }

  auto add_candidate = [&](std::string_view taxid, int count) {
    if (count <= 0) {
      return;
    }
    const size_t c = static_cast<size_t>(count);
    if (threshold > 0 && c < threshold) {
      return;
    }

    if (topK == 0) {
      hits.emplace_back(std::string(taxid), c);
      return;
    }

    if (hits.size() < topK) {
      hits.emplace_back(std::string(taxid), c);
      return;
    }

    auto min_it = std::min_element(
        hits.begin(), hits.end(),
        [](const auto &a, const auto &b) { return a.second < b.second; });
    if (min_it == hits.end() || c <= min_it->second) {
      return;
    }
    *min_it = {std::string(taxid), c};
  };

  if (!db.multiblock) {
    auto counts = db.icf.bulk_count(hashes);
    for (size_t bin = 0; bin < db.indexToTaxid.size(); ++bin) {
      add_candidate(db.indexToTaxid[bin], counts[bin]);
    }
  } else {
    for (auto &block : db.db2.blocks) {
      auto counts = block.icf.bulk_count(hashes);
      for (size_t bin = 0; bin < block.indexToTaxid.size(); ++bin) {
        add_candidate(block.indexToTaxid[bin], counts[bin]);
      }
    }
  }

  if (hits.empty()) {
    result.taxidCount = {{kUnclassified, 1}};
    return result;
  }

  std::sort(hits.begin(), hits.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });

  result.taxidCount = std::move(hits);
  return result;
}

void saveResult(const std::vector<classifyResult> &classifyResults,
                const ClassifyConfig &config) {
  std::string outputFile = config.outputFile;
  if (std::filesystem::path(outputFile).extension() != ".tsv") {
    outputFile += ".tsv";
  }
  std::ofstream os(outputFile, std::ios::out);
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open file: " + config.outputFile);
  }

  for (const auto &result : classifyResults) {
    os << result.id << '\t';

    bool handled = false;
    if (!result.taxidCount.empty() &&
        result.taxidCount.front().first == kUnclassified) {
      os << kUnclassified;
      handled = true;
    } else if (config.output_posterior && !result.posteriors.empty()) {
      auto oldFlags = os.flags();
      auto oldPrecision = os.precision();
      os.setf(std::ios::fixed, std::ios::floatfield);
      os << std::setprecision(4) << result.posteriors.front().first << ':'
         << result.posteriors.front().second;
      if (result.posteriors.size() > 1) {
        os << '\t' << "POST_TOP2=" << result.posteriors[1].first << ':'
           << result.posteriors[1].second;
      }
      os.flags(oldFlags);
      os.precision(oldPrecision);
      handled = true;
    }

    if (!handled) {
      for (const auto &[taxid, count] : result.taxidCount) {
        if (taxid == kUnclassified) {
          os << taxid;
          continue;
        }
        os << taxid << ':' << count << '\t';
      }
    }
    os << '\n';
  }
}
} // namespace

void postEmDecision(std::vector<classifyResult> &results,
                    const DecisionConfig &decisionConfig,
                    const std::unordered_map<std::string, double> &classWeights) {
  auto prune_by_global_pi = [&](classifyResult &res, double pi_min) {
    if (pi_min <= 0.0 || res.posteriors.empty() || classWeights.empty()) {
      return;
    }
    std::vector<std::pair<std::string, double>> kept;
    kept.reserve(res.posteriors.size());
    double sum = 0.0;
    for (const auto &kv : res.posteriors) {
      auto weight_it = classWeights.find(kv.first);
      double w = (weight_it != classWeights.end()) ? weight_it->second : 0.0;
      if (w >= pi_min) {
        kept.push_back(kv);
        sum += kv.second;
      }
    }
    if (kept.empty()) {
      res.posteriors.clear();
      return;
    }
    if (sum > 0.0) {
      for (auto &kv : kept) {
        kv.second /= sum;
      }
    }
    res.posteriors.swap(kept);
  };

  for (auto &result : results) {
    if (decisionConfig.min_class_weight > 0.0) {
      prune_by_global_pi(result, decisionConfig.min_class_weight);
    }

    if (result.posteriors.empty()) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(kUnclassified, 1);
      continue;
    }

    auto top1 = result.posteriors.front();
    auto top2 = (result.posteriors.size() > 1) ? result.posteriors[1]
                                               : std::make_pair(std::string{}, 0.0);

    bool pass = true;
    if (decisionConfig.posterior_threshold > 0.0 &&
        top1.second < decisionConfig.posterior_threshold) {
      pass = false;
    }
    if (pass && decisionConfig.margin_delta > 0.0 &&
        (top1.second - top2.second) < decisionConfig.margin_delta) {
      pass = false;
    }
    if (pass && decisionConfig.margin_ratio > 0.0 && top2.second > 0.0) {
      double ratio = top1.second / top2.second;
      if (ratio < decisionConfig.margin_ratio) {
        pass = false;
      }
    }

    if (!pass) {
      result.taxidCount.clear();
      result.taxidCount.emplace_back(kUnclassified, 1);
      continue;
    }

    result.taxidCount.clear();
    result.taxidCount.emplace_back(top1.first, 0);
  }
}

void run(ClassifyConfig config) {
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

  if (!config.em && !config.vem) {
    config.em = true;
  }
  if (!(config.post_ratio > 0.0)) {
    config.post_ratio = std::numeric_limits<double>::quiet_NaN();
  }

  if (config.verbose) {
    std::cout << config << std::endl;
  }

  seqan3::contrib::bgzf_thread_count = config.threads;
#ifdef _OPENMP
  omp_set_num_threads(static_cast<int>(config.threads));
#endif

  auto totalStart = std::chrono::high_resolution_clock::now();
  auto loadStart = std::chrono::high_resolution_clock::now();
  std::cout << "Loading ICF database..." << std::endl;
  ICFDatabase db = load_icf_db(config.dbFile);
  auto loadEnd = std::chrono::high_resolution_clock::now();
  if (config.verbose) {
    std::cout << "Load time: ";
    print_classify_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                            loadEnd - loadStart)
                            .count());
  }

  const std::array<size_t, 1> syncmer_positions{
      static_cast<size_t>(db.config.syncmerPosition)};
  const std::span<const size_t> syncmer_pos_span(syncmer_positions);

  FileInfo fileInfo{};
  fileInfo.minLen = std::numeric_limits<size_t>::max();
  fileInfo.maxLen = 0;

  std::vector<classifyResult> classifyResults;
  classifyResults.reserve(1024);

  auto classifyStart = std::chrono::high_resolution_clock::now();
  std::cout << "Classifying sequences by icf..." << std::endl;

  auto handle_record = [&](FileInfo &info,
                           std::vector<classifyResult> &out,
                           const std::string &id,
                           const std::vector<seqan3::dna4> &seq1,
                           const std::vector<seqan3::dna4> *seq2 = nullptr) {
    std::vector<uint64_t> hashes;
    size_t readLen = seq1.size() + (seq2 ? seq2->size() : 0);
    update_minmax(readLen, info.minLen, info.maxLen);
    info.bpLength += readLen;
    info.sequenceNum += 1;

    if (seq1.size() >= db.config.kmerSize) {
      hashes = taxicf::syncmer::compute_hashes(seq1, db.config.smerSize,
                                                db.config.kmerSize,
                                                syncmer_pos_span,
                                                db.config.seed64, true);
    }
    if (seq2 && seq2->size() >= db.config.kmerSize) {
      auto hashes2 = taxicf::syncmer::compute_hashes(*seq2, db.config.smerSize,
                                                      db.config.kmerSize,
                                                      syncmer_pos_span,
                                                      db.config.seed64, true);
      hashes.insert(hashes.end(), hashes2.begin(), hashes2.end());
    }
    if (hashes.size() > 2048) {
      std::sort(hashes.begin(), hashes.end());
      hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());
    }

    auto res = classify_one(id, hashes, config, db);
    if (!res.taxidCount.empty() && res.taxidCount.front().first == kUnclassified) {
      info.unclassifiedNum += 1;
    } else {
      info.classifiedNum += 1;
    }

    out.emplace_back(std::move(res));
  };

  if (!config.pairedFiles.empty()) {
    if (config.pairedFiles.size() % 2 != 0) {
      throw std::runtime_error("--paired 必须提供偶数个文件");
    }

    const size_t pairCount = config.pairedFiles.size() / 2;
#pragma omp parallel
    {
      FileInfo localInfo{};
      localInfo.minLen = std::numeric_limits<size_t>::max();
      localInfo.maxLen = 0;
      std::vector<classifyResult> localResults;
      localResults.reserve(1024);

#pragma omp for schedule(dynamic)
      for (size_t pairIdx = 0; pairIdx < pairCount; ++pairIdx) {
        const auto &p1 = config.pairedFiles[pairIdx * 2];
        const auto &p2 = config.pairedFiles[pairIdx * 2 + 1];
        seqan3::sequence_file_input<raptor::dna4_traits,
                                   seqan3::fields<seqan3::field::id,
                                                 seqan3::field::seq>>
            fin1{p1};
        seqan3::sequence_file_input<raptor::dna4_traits,
                                   seqan3::fields<seqan3::field::id,
                                                 seqan3::field::seq>>
            fin2{p2};

        auto it1 = fin1.begin();
        auto it2 = fin2.begin();
        for (; it1 != fin1.end() && it2 != fin2.end(); ++it1, ++it2) {
          const auto &rec1 = *it1;
          const auto &rec2 = *it2;
          handle_record(localInfo, localResults, std::string{rec1.id()},
                        rec1.sequence(), &rec2.sequence());
        }
      }

#pragma omp critical
      {
        if (localInfo.sequenceNum > 0) {
          fileInfo.sequenceNum += localInfo.sequenceNum;
          fileInfo.classifiedNum += localInfo.classifiedNum;
          fileInfo.unclassifiedNum += localInfo.unclassifiedNum;
          fileInfo.bpLength += localInfo.bpLength;
          fileInfo.minLen = std::min(fileInfo.minLen, localInfo.minLen);
          fileInfo.maxLen = std::max(fileInfo.maxLen, localInfo.maxLen);
        }
        classifyResults.reserve(classifyResults.size() + localResults.size());
        std::move(localResults.begin(), localResults.end(),
                  std::back_inserter(classifyResults));
      }
    }
  } else {
#pragma omp parallel
    {
      FileInfo localInfo{};
      localInfo.minLen = std::numeric_limits<size_t>::max();
      localInfo.maxLen = 0;
      std::vector<classifyResult> localResults;
      localResults.reserve(1024);

#pragma omp for schedule(dynamic)
      for (size_t fileIdx = 0; fileIdx < config.singleFiles.size(); ++fileIdx) {
        const auto &path = config.singleFiles[fileIdx];
        seqan3::sequence_file_input<raptor::dna4_traits,
                                   seqan3::fields<seqan3::field::id,
                                                 seqan3::field::seq>>
            fin{path};
        for (auto &record : fin) {
          handle_record(localInfo, localResults, std::string{record.id()},
                        record.sequence(), nullptr);
        }
      }

#pragma omp critical
      {
        if (localInfo.sequenceNum > 0) {
          fileInfo.sequenceNum += localInfo.sequenceNum;
          fileInfo.classifiedNum += localInfo.classifiedNum;
          fileInfo.unclassifiedNum += localInfo.unclassifiedNum;
          fileInfo.bpLength += localInfo.bpLength;
          fileInfo.minLen = std::min(fileInfo.minLen, localInfo.minLen);
          fileInfo.maxLen = std::max(fileInfo.maxLen, localInfo.maxLen);
        }
        classifyResults.reserve(classifyResults.size() + localResults.size());
        std::move(localResults.begin(), localResults.end(),
                  std::back_inserter(classifyResults));
      }
    }
  }

  auto classifyEnd = std::chrono::high_resolution_clock::now();
  if (config.verbose) {
    std::cout << "Classify time: ";
    print_classify_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                            classifyEnd - classifyStart)
                            .count());
  }

  std::unordered_map<std::string, double> classWeights;
  bool posteriorModelUsed = false;

  if (config.em) {
    auto EMstart = std::chrono::high_resolution_clock::now();
    std::cout << "Running EM algorithm..." << std::endl;
    EMOptions options;
    auto [posterior, weights] = EMAlgorithm(std::move(classifyResults), config.emIter,
                                            config.emThreshold, options);
    classifyResults = std::move(posterior);
    classWeights = std::move(weights);
    posteriorModelUsed = true;
    auto EMend = std::chrono::high_resolution_clock::now();
    if (config.verbose) {
      std::cout << "EM time: ";
      print_classify_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                              EMend - EMstart)
                              .count());
    }
  }

  if (config.vem) {
    auto VEMstart = std::chrono::high_resolution_clock::now();
    std::cout << "Running VEM algorithm..." << std::endl;
    VEMOptions options;
    auto [posterior, weights] = VEMAlgorithm(std::move(classifyResults), config.emIter,
                                             config.emThreshold, options);
    classifyResults = std::move(posterior);
    classWeights = std::move(weights);
    posteriorModelUsed = true;
    auto VEMend = std::chrono::high_resolution_clock::now();
    if (config.verbose) {
      std::cout << "VEM time: ";
      print_classify_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                              VEMend - VEMstart)
                              .count());
    }
  }

  if (posteriorModelUsed) {
    DecisionConfig decisionConfig;
    decisionConfig.posterior_threshold = config.post_thres;
    decisionConfig.margin_delta = config.post_margin;
    decisionConfig.margin_ratio = config.post_ratio;
    decisionConfig.min_class_weight = config.post_pi_min;

    postEmDecision(classifyResults, decisionConfig, classWeights);

    fileInfo.classifiedNum = 0;
    fileInfo.unclassifiedNum = 0;
    for (const auto &result : classifyResults) {
      if (!result.taxidCount.empty() &&
          result.taxidCount.front().first == kUnclassified) {
        ++fileInfo.unclassifiedNum;
      } else {
        ++fileInfo.classifiedNum;
      }
    }
  }

  auto saveStart = std::chrono::high_resolution_clock::now();
  std::cout << "Saving classification results..." << std::endl;
  saveResult(classifyResults, config);
  auto saveEnd = std::chrono::high_resolution_clock::now();
  if (config.verbose) {
    std::cout << "Save time: ";
    print_classify_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                            saveEnd - saveStart)
                            .count());
    std::cout << "Total sequences: " << fileInfo.sequenceNum << std::endl;
    std::cout << "Classified sequences: " << fileInfo.classifiedNum << std::endl;
    std::cout << "Unclassified sequences: " << fileInfo.unclassifiedNum
              << std::endl;
  }

  auto totalEnd = std::chrono::high_resolution_clock::now();
  if (config.verbose) {
    std::cout << "Total time: ";
    print_classify_time(std::chrono::duration_cast<std::chrono::milliseconds>(
                            totalEnd - totalStart)
                            .count());
  }
}
} // namespace TAXICFClassify
