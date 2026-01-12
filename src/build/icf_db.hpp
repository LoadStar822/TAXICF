/*
 * -----------------------------------------------------------------------------
 * Filename:      icf_db.hpp
 *
 * Description:
 *  TAXICF ICF database IO helpers.
 *
 *  - v1 (legacy): cereal archive containing (icf, indexToTaxid, ICFConfig)
 *  - v2 (current): file starts with a raw magic header, followed by a cereal
 *    archive containing a multi-block database.
 * -----------------------------------------------------------------------------
 */
#pragma once

#include <buildConfig.hpp>
#include <interleaved-cuckoo-filter.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace TAXICFBuild {

inline constexpr std::array<char, 8> kIcfDbMagic{{'T', 'A', 'X', 'I', 'C', 'F', 'D', 'B'}};
inline constexpr uint8_t kIcfDbVersionV2 = 2;

struct IcfDbGlobalConfigV2 {
  uint8_t dbVersion{kIcfDbVersionV2};
  uint8_t kmerSize{};
  uint16_t smerSize{};
  uint16_t syncmerPosition{0};
  double loadFactor{0.0};
  uint64_t seed64{0};
  uint8_t tagBits{8};

  template <class Archive>
  void serialize(Archive &archive) {
    archive(dbVersion, kmerSize, smerSize, syncmerPosition, loadFactor, seed64,
            tagBits);
  }
};

struct IcfDbBlockV2 {
  size_t binSize{};
  taxicf::InterleavedCuckooFilter icf;
  std::vector<std::string> indexToTaxid;

  template <class Archive>
  void serialize(Archive &archive) {
    archive(binSize, icf, indexToTaxid);
  }
};

struct IcfDbV2 {
  IcfDbGlobalConfigV2 config;
  std::vector<IcfDbBlockV2> blocks;

  template <class Archive>
  void serialize(Archive &archive) {
    archive(config, blocks);
  }
};

inline bool try_read_v2_header(std::istream &is, uint8_t &version_out) {
  std::array<char, 8> magic{};
  is.read(magic.data(), static_cast<std::streamsize>(magic.size()));
  if (!is.good()) {
    return false;
  }
  if (magic != kIcfDbMagic) {
    is.clear();
    is.seekg(0);
    return false;
  }
  uint8_t version = 0;
  is.read(reinterpret_cast<char *>(&version), 1);
  if (!is.good()) {
    throw std::runtime_error("ICF DB header truncated (missing version byte).");
  }
  version_out = version;
  return true;
}

inline void write_v2_db(const std::string &output_base,
                        const IcfDbV2 &db) {
  const std::string outPath = output_base + ".icf";
  std::ofstream os(outPath, std::ios::binary);
  if (!os.is_open()) {
    throw std::runtime_error("Failed to open file: " + outPath);
  }

  os.write(kIcfDbMagic.data(), static_cast<std::streamsize>(kIcfDbMagic.size()));
  const uint8_t ver = kIcfDbVersionV2;
  os.write(reinterpret_cast<const char *>(&ver), 1);
  if (!os.good()) {
    throw std::runtime_error("Failed to write ICF DB header: " + outPath);
  }

  cereal::BinaryOutputArchive archive(os);
  archive(db);
}

inline IcfDbV2 read_v2_db(std::istream &is) {
  uint8_t version = 0;
  if (!try_read_v2_header(is, version)) {
    throw std::runtime_error("Not a v2 ICF DB (missing magic header).");
  }
  if (version != kIcfDbVersionV2) {
    throw std::runtime_error("Unsupported ICF DB version: " +
                             std::to_string(version));
  }
  cereal::BinaryInputArchive archive(is);
  IcfDbV2 db;
  archive(db);
  if (db.config.dbVersion != kIcfDbVersionV2) {
    throw std::runtime_error("Invalid v2 DB payload version: " +
                             std::to_string(db.config.dbVersion));
  }
  return db;
}

} // namespace TAXICFBuild
