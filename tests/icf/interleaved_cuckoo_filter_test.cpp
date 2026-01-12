#include <interleaved-cuckoo-filter.h>

#include <cereal/archives/binary.hpp>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
template <typename T, typename U>
void expect_equal(const T &lhs, const U &rhs, const std::string &message) {
  if (!(lhs == rhs)) {
    throw std::runtime_error(message + "：实际值=" + std::to_string(lhs) +
                             "，期望值=" + std::to_string(rhs));
  }
}

void expect_true(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

uint64_t splitmix64(uint64_t &state) {
  uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

std::vector<uint64_t> make_values(size_t n, uint64_t seed) {
  std::vector<uint64_t> values;
  values.reserve(n);
  uint64_t state = seed;
  for (size_t i = 0; i < n; ++i) {
    uint64_t v = splitmix64(state);
    if (v == 0) {
      v = 1;
    }
    values.push_back(v);
  }
  return values;
}

struct TestCase {
  std::string name;
  std::function<void()> body;
};

class TestRunner {
public:
  void add(std::string name, std::function<void()> body) {
    tests_.push_back(TestCase{std::move(name), std::move(body)});
  }

  int run() const {
    auto format_duration = [](double ms) {
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(3) << ms;
      return oss.str();
    };

    struct Result {
      std::string name;
      double elapsed_ms;
      bool success;
      std::string message;
    };

    std::vector<Result> results;
    results.reserve(tests_.size());

    int failures = 0;
    double total_ms = 0.0;
    for (const auto &test : tests_) {
      auto start = std::chrono::steady_clock::now();
      bool success = true;
      std::string message;
      try {
        test.body();
      } catch (const std::exception &ex) {
        success = false;
        message = ex.what();
        ++failures;
      }
      auto end = std::chrono::steady_clock::now();
      double elapsed_ms =
          std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
              end - start)
              .count();
      total_ms += elapsed_ms;
      results.push_back(Result{test.name, elapsed_ms, success, message});
    }

    size_t nameWidth = std::string("用例").size();
    for (const auto &result : results) {
      nameWidth = std::max(nameWidth, result.name.size());
    }
    const size_t statusWidth = std::string("结果").size();
    const size_t timeWidth = std::string("耗时 (ms)").size();
    const size_t totalWidth =
        nameWidth + statusWidth + timeWidth + 8; // padding and separators
    std::string divider(totalWidth, '-');

    std::cout << "\n测试结果\n" << divider << '\n';
    std::cout << std::left << std::setw(nameWidth) << "用例"
              << " | " << std::setw(statusWidth) << "结果"
              << " | " << std::right << std::setw(timeWidth) << "耗时 (ms)"
              << '\n';
    std::cout << divider << '\n';

    for (const auto &result : results) {
      std::string timeString = format_duration(result.elapsed_ms);
      std::cout << std::left << std::setw(nameWidth) << result.name << " | "
                << std::setw(statusWidth) << (result.success ? "通过" : "失败")
                << " | " << std::right << std::setw(timeWidth) << timeString
                << '\n';
      std::cout << std::left;
    }
    std::cout << divider << '\n';

    if (failures > 0) {
      std::cout << "\n失败用例详情：" << '\n';
      for (const auto &result : results) {
        if (!result.success) {
          std::cout << " - " << result.name << "：" << result.message << '\n';
        }
      }
    }

    std::cout << "\n共 " << results.size() << " 个用例，其中 "
              << (results.size() - failures) << " 个通过，用时 "
              << format_duration(total_ms) << " ms" << std::endl;

    return failures == 0 ? 0 : 1;
  }

private:
  std::vector<TestCase> tests_;
};
} // namespace

int main() {
  TestRunner runner;

  runner.add("ICF bulk_count：空输入返回全 0", []() {
    constexpr size_t bins = 8;
    taxicf::InterleavedCuckooFilter icf(bins, /*bin_size=*/64, /*bitNum=*/16);

    std::vector<uint64_t> empty;
    auto counts = icf.bulk_count(empty);

    expect_equal(static_cast<size_t>(counts.size()), bins, "counts bins size mismatch");
    for (size_t i = 0; i < bins; ++i) {
      expect_equal(counts[i], 0, "empty query should yield zero hits");
    }
  });

  runner.add("ICF 插入：同一 value 可写入不同 bin", []() {
    constexpr size_t bins = 6;
    taxicf::InterleavedCuckooFilter icf(bins, /*bin_size=*/256, /*bitNum=*/16);

    auto values = make_values(/*n=*/128, /*seed=*/0x1234);
    for (auto v : values) {
      icf.insertTag(/*binIndex=*/1, v);
      icf.insertTag(/*binIndex=*/4, v);
    }

    auto counts = icf.bulk_count(values);
    expect_equal(static_cast<size_t>(counts.size()), bins, "counts bins size mismatch");
    expect_equal(counts[1], static_cast<int>(values.size()), "bin 1 hit count mismatch");
    expect_equal(counts[4], static_cast<int>(values.size()), "bin 4 hit count mismatch");
  });

  runner.add("ICF bulk_count：按 bin 查询无漏报", []() {
    constexpr size_t bins = 4;
    taxicf::InterleavedCuckooFilter icf(bins, /*bin_size=*/512, /*bitNum=*/16);

    std::vector<std::vector<uint64_t>> per_bin;
    per_bin.reserve(bins);
    for (size_t bin = 0; bin < bins; ++bin) {
      per_bin.push_back(make_values(/*n=*/256, /*seed=*/0xBEEF + bin * 17));
      for (auto v : per_bin.back()) {
        icf.insertTag(bin, v);
      }
    }

    for (size_t bin = 0; bin < bins; ++bin) {
      auto counts = icf.bulk_count(per_bin[bin]);
      expect_equal(static_cast<size_t>(counts.size()), bins, "counts bins size mismatch");
      expect_equal(counts[bin], static_cast<int>(per_bin[bin].size()), "hit count mismatch for bin query");
    }
  });

  runner.add("ICF 序列化往返：插入值不丢失", []() {
    constexpr size_t bins = 4;
    taxicf::InterleavedCuckooFilter icf(bins, /*bin_size=*/256, /*bitNum=*/16);

    auto values0 = make_values(/*n=*/128, /*seed=*/0xABCDEF01);
    auto values2 = make_values(/*n=*/128, /*seed=*/0xABCDEF02);
    for (auto v : values0) {
      icf.insertTag(/*binIndex=*/0, v);
    }
    for (auto v : values2) {
      icf.insertTag(/*binIndex=*/2, v);
    }

    std::stringstream buffer;
    {
      cereal::BinaryOutputArchive archive(buffer);
      archive(icf);
    }

    taxicf::InterleavedCuckooFilter restored;
    {
      cereal::BinaryInputArchive archive(buffer);
      archive(restored);
    }

    {
      auto counts = restored.bulk_count(values0);
      expect_equal(static_cast<size_t>(counts.size()), bins, "counts bins size mismatch");
      expect_equal(counts[0], static_cast<int>(values0.size()), "restored bin 0 hit count mismatch");
    }
    {
      auto counts = restored.bulk_count(values2);
      expect_equal(static_cast<size_t>(counts.size()), bins, "counts bins size mismatch");
      expect_equal(counts[2], static_cast<int>(values2.size()), "restored bin 2 hit count mismatch");
    }
  });

  runner.add("ICF bitNum=8：基本插入与查询", []() {
    constexpr size_t bins = 3;
    taxicf::InterleavedCuckooFilter icf(bins, /*bin_size=*/256, /*bitNum=*/8);

    auto values = make_values(/*n=*/256, /*seed=*/0x777);
    for (auto v : values) {
      icf.insertTag(/*binIndex=*/2, v);
    }

    auto counts = icf.bulk_count(values);
    expect_equal(static_cast<size_t>(counts.size()), bins, "counts bins size mismatch");
    expect_equal(counts[2], static_cast<int>(values.size()), "bin 2 hit count mismatch");
  });

  return runner.run();
}
