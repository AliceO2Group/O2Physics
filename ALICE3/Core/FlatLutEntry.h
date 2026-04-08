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

#ifndef ALICE3_CORE_FLATLUTENTRY_H_
#define ALICE3_CORE_FLATLUTENTRY_H_

#include <cmath>
#include <cstdint>
#include <vector>

#define LUTCOVM_VERSION 20260408

namespace o2::delphes
{

/**
 * @brief Flat LUT entry structure
 */
struct lutEntry_t {
  float nch = 0.f;
  float eta = 0.f;
  float pt = 0.f;
  bool valid = false;
  float eff = 0.f;
  float eff2 = 0.f;
  float itof = 0.f;
  float otof = 0.f;
  float covm[15] = {0.f};
  float eigval[5] = {0.f};
  float eigvec[5][5] = {{0.f}};
  float eiginv[5][5] = {{0.f}};

  void print() const;
};

/**
 * @brief Binning map
 */
struct map_t {
  int nbins = 1;
  float min = 0.f;
  float max = 1.e6f;
  bool log = false;

  float eval(int bin) const
  {
    float width = (max - min) / nbins;
    float val = min + (bin + 0.5f) * width;
    if (log)
      return std::pow(10.f, val);
    return val;
  }

  float fracPositionWithinBin(float val) const;
  int find(float val) const;
  void print() const;
};

/**
 * @brief LUT header
 */
struct lutHeader_t {
  int version = LUTCOVM_VERSION;
  int pdg = 0;
  float mass = 0.f;
  float field = 0.f;
  map_t nchmap;
  map_t radmap;
  map_t etamap;
  map_t ptmap;

  bool check_version() const;
  void print() const;
};

/**
 * @brief Flat LUT data container - single contiguous buffer
 * Memory layout: [header][entry_0][entry_1]...[entry_N]
 *
 * All entries stored sequentially in a single allocation.
 * Can be directly mapped from file or shared memory without copying.
 */
class FlatLutData
{
 public:
  FlatLutData() = default;

  /**
   * @brief Initialize from binning information
   * Pre-allocates contiguous memory for all entries
   */
  void initialize(const lutHeader_t& header);

  /**
   * @brief Get LUT entry by bin indices
   * O(1) access via linear index calculation
   */
  lutEntry_t* getEntry(int nch_bin, int rad_bin, int eta_bin, int pt_bin);
  const lutEntry_t* getEntry(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const;

  /**
   * @brief Get raw data buffer (for serialization/shared memory)
   */
  uint8_t* data() { return mData.data(); }
  const uint8_t* data() const { return mData.data(); }

  /**
   * @brief Total size in bytes
   */
  size_t bytes() const { return mData.size(); }

  /**
   * @brief Construct from external buffer (e.g., shared memory or file mapping)
   */
  static FlatLutData fromBuffer(const uint8_t* buffer, size_t size);

  /**
   * @brief Reference-based access without copying
   * Useful when data is already in shared memory
   */
  static FlatLutData fromExternalBuffer(uint8_t* buffer, size_t size);

  const lutHeader_t& header() const
  {
    return *reinterpret_cast<const lutHeader_t*>(mData.data());
  }

  lutHeader_t& header()
  {
    return *reinterpret_cast<lutHeader_t*>(mData.data());
  }

 private:
  /**
   * @brief Linear index calculation for entry access
   */
  size_t getEntryOffset(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const;

  std::vector<uint8_t> mData;

  // Cache dimensions for quick access
  int mNchBins = 0;
  int mRadBins = 0;
  int mEtaBins = 0;
  int mPtBins = 0;
};

} // namespace o2::delphes

#endif // ALICE3_CORE_FLATLUTENTRY_H_
