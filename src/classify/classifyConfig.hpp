/*
 * -----------------------------------------------------------------------------
 * Filename:      classifyConfig.hpp
 *
 * Description:
 *  TAXICF classify configuration (ICF-only, simplified).
 * -----------------------------------------------------------------------------
 */
#pragma once

#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace TAXICFClassify {
struct ClassifyConfig {
  std::vector<std::string> singleFiles;
  std::vector<std::string> pairedFiles;
  std::string outputFile{"TAXICFClassify"};
  std::string dbFile;

  double shotThreshold = 0.62;
  size_t preEmTopK = 32;

  uint16_t threads{0};
  bool verbose{true};

  // Optional posterior refinement.
  bool em{false};
  bool vem{false};
  double emThreshold{0.001};
  size_t emIter{80};

  // Posterior decision thresholds (used when EM/VEM enabled).
  double post_thres{0.56};
  double post_margin{0.03};
  double post_ratio{1.30};
  double post_pi_min{1e-4};

  // Output
  bool output_posterior{true};
};

inline std::ostream &operator<<(std::ostream &os, const ClassifyConfig &config) {
  os << std::string(40, '=') << '\n';
  os << " Classify Configuration " << '\n';
  os << std::string(40, '=') << '\n';

  os << std::left << std::setw(20) << "Single files:" << '\n';
  for (const auto &file : config.singleFiles) {
    os << std::setw(20) << "" << file << '\n';
  }
  os << std::setw(20) << "Paired files:" << '\n';
  for (const auto &file : config.pairedFiles) {
    os << std::setw(20) << "" << file << '\n';
  }

  os << std::setw(20) << "Output file:" << config.outputFile << '\n'
     << std::setw(20) << "Database file:" << config.dbFile << '\n'
     << std::setw(20) << "Shot threshold:" << config.shotThreshold << '\n'
     << std::setw(20) << "Pre-EM topK:" << config.preEmTopK << '\n'
     << std::setw(20) << "EM:" << config.em << '\n'
     << std::setw(20) << "VEM:" << config.vem << '\n'
     << std::setw(20) << "Threads:" << config.threads << '\n'
     << std::setw(20) << "Verbose:" << config.verbose << '\n'
     << std::setw(20) << "Posterior thres:" << config.post_thres << '\n'
     << std::setw(20) << "Posterior margin:" << config.post_margin << '\n'
     << std::setw(20) << "Posterior ratio:"
     << (std::isnan(config.post_ratio) ? std::string("nan")
                                       : std::to_string(config.post_ratio))
     << '\n'
     << std::setw(20) << "Posterior pi min:" << config.post_pi_min << '\n'
     << std::setw(20) << "Output posterior:" << config.output_posterior << '\n';

  os << std::string(40, '=') << '\n';
  return os;
}

struct FileInfo {
  size_t sequenceNum{0};
  size_t classifiedNum{0};
  size_t unclassifiedNum{0};
  size_t minLen{0};
  size_t maxLen{0};
  size_t bpLength{0};

  std::unordered_set<std::string> uniqueTaxids;
  std::unordered_map<std::string, size_t> taxidTotalMatches;
};

struct classifyResult {
  std::string id;
  std::vector<std::pair<std::string, size_t>> taxidCount;
  std::vector<std::pair<std::string, double>> posteriors;
  double evaluated{0.0}; // 实际参与判别的 syncmer 数，用于归一化
};

struct DecisionConfig {
  double posterior_threshold{0.56};
  double margin_delta{0.03};
  double margin_ratio{1.30};
  double min_class_weight{1e-4};
};
} // namespace TAXICFClassify
