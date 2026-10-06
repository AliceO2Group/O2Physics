// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FlatLutEntry.h
/// \brief Flat LUT data structures and buffer handling for the ALICE3 fast smearing backend.

#ifndef ALICE3_CORE_FLATLUTENTRY_H_
#define ALICE3_CORE_FLATLUTENTRY_H_

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <span>
#include <vector>

static constexpr int LUTCOVM_VERSION = 20210801;

namespace o2::fastsim
{

constexpr int kNumCovarianceTerms = 15;
constexpr int kNumEigenModes = 5;
typedef std::array<float, kNumCovarianceTerms> CovarianceArray;
typedef std::array<float, kNumEigenModes> EigenArray;
typedef std::array<double, kNumEigenModes> EigenArrayDouble;
typedef std::array<std::array<float, kNumEigenModes>, kNumEigenModes> EigenMatrix;

/**
 * @brief Flat LUT entry structure.
 */
struct LutEntry {
  float nch = 0.f;
  float eta = 0.f;
  float pt = 0.f;
  bool valid = false;
  float eff = 0.f;
  float eff2 = 0.f;
  float itof = 0.f;
  float otof = 0.f;
  CovarianceArray covm = {0.f};
  EigenArray eigval = {0.f};
  EigenMatrix eigvec = {{{0.f}}};
  EigenMatrix eiginv = {{{0.f}}};

  void print() const;
};

using lutEntry_t = LutEntry;

/**
 * @brief Binning map.
 */
struct Map {
  int nbins = 1;
  float min = 0.f;
  float max = 1.e6f;
  bool log = false;

  [[nodiscard]] float eval(int bin) const
  {
    float width = (max - min) / nbins;
    float val = min + (bin + 0.5f) * width;
    if (log) {
      return std::pow(10.f, val);
    }
    return val;
  }

  [[nodiscard]] float fracPositionWithinBin(float val) const;
  [[nodiscard]] int find(float val) const;
  void print() const;
};

using map_t = Map;

/**
 * @brief LUT header.
 */
struct LutHeader {
  int version = LUTCOVM_VERSION;
  int pdg = 0;
  float mass = 0.f;
  float field = 0.f;
  map_t nchmap;
  map_t radmap;
  map_t etamap;
  map_t ptmap;

  [[nodiscard]] bool checkVersion() const;
  void print() const;
};

using lutHeader_t = LutHeader;

/**
 * @brief Flat LUT data container - single contiguous buffer.
 * Memory layout: [header][entry_0][entry_1]...[entry_N]
 *
 * All entries stored sequentially in a single allocation.
 * Can be directly mapped from file or shared memory without copying.
 */
class FlatLutData
{
 public:
  FlatLutData() = default;
  FlatLutData(FlatLutData&&) = default;
  FlatLutData& operator=(FlatLutData&&) = default;

  /**
   * @brief Initialize from binning information
   * Pre-allocates contiguous memory for all entries
   */
  void initialize(const lutHeader_t& header);

  [[nodiscard]] size_t getEntryIndex(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const
  {
    // Linear index: nch varies slowest, pt varies fastest
    // idx = nch * (rad*eta*pt) + rad * (eta*pt) + eta * pt + pt
    return static_cast<size_t>(nch_bin) * (mRadBins * mEtaBins * mPtBins) + static_cast<size_t>(rad_bin) * (mEtaBins * mPtBins) + static_cast<size_t>(eta_bin) * mPtBins + static_cast<size_t>(pt_bin);
  }

  /**
   * @brief Get LUT entry by bin indices (view)
   */
  [[nodiscard]] const lutEntry_t* getEntryRef(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const;

  /**
   * @brief Get LUT entry by bin indices (owned)
   */
  lutEntry_t* getEntry(int nch_bin, int rad_bin, int eta_bin, int pt_bin);

  /**
   * @brief Get LUT header (view)
   */
  [[nodiscard]] const lutHeader_t& getHeaderRef() const;

  /**
   * @brief Get LUT header (owned)
   */
  lutHeader_t& getHeader();

  /**
   * @brief Get raw data buffer
   */
  uint8_t* data() { return mData.data(); }                              // owned
  [[nodiscard]] const uint8_t* data() const { return mDataRef.data(); } // view

  /**
   * @brief Total size in bytes
   */
  [[nodiscard]] size_t bytes() const { return mDataRef.size(); }

  /**
   * @brief Construct a new FlatLutData from external buffer as a copy
   */
  static FlatLutData AdoptFromBuffer(const uint8_t* buffer, size_t size);

  /**
   * @brief Construct a new FlatLutData from external buffer as a view
   */
  static FlatLutData ViewFromBuffer(const uint8_t* buffer, size_t size);

  /**
   * @brief Construct a new FlatLutData from external span as a view
   */
  static FlatLutData ViewFromBuffer(std::span<std::byte> const& span);

  /**
   * @brief Construct a new FlatLutData from a file
   */
  static FlatLutData loadFromFile(std::ifstream& file, const char* filename);

  /**
   * @brief Preview buffer header for version and other compatibility checks
   */
  static lutHeader_t previewHeader(const uint8_t* buffer, size_t size);

  /**
   * @brief Preview file-stored header for version and other compatibility checks
   */
  static lutHeader_t previewHeader(std::ifstream& file, const char* filename);

  /**
   * @brief Check if the LUT is loaded
   */
  [[nodiscard]] bool isLoaded() const;

  /**
   * @brief Reset LUT to empty
   */
  void reset();

 private:
  /**
   * @brief Linear index calculation for entry access
   */
  [[nodiscard]] size_t getEntryOffset(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const;

  /**
   * @brief Update dimensions from the current header
   */
  void cacheDimensions();
  /**
   * @brief Reset dimensions to 0
   */
  void resetDimensions();

  /**
   * @brief Adopt a buffer by copying
   */
  void adopt(const uint8_t* buffer, size_t size);
  /**
   * @brief Adopt a buffer as a view
   */
  void view(const uint8_t* buffer, size_t size);
  /**
   * @brief Update mDataRef from mData
   */
  void updateRef();
  /**
   * @brief Validate external buffer
   */
  static void validateBuffer(const uint8_t* buffer, size_t size);

  std::vector<uint8_t> mData;
  std::span<uint8_t const> mDataRef;

  // Cache dimensions for quick access
  int mNchBins = 0;
  int mRadBins = 0;
  int mEtaBins = 0;
  int mPtBins = 0;
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_FLATLUTENTRY_H_
