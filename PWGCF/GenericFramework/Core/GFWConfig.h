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

/// \file GFWConfig.h
/// \brief Utility classes for Generic Framework
/// \author Emil Gorm Nielsen, NBI, emil.gorm.nielsen@cern.ch

#ifndef PWGCF_GENERICFRAMEWORK_CORE_GFWCONFIG_H_
#define PWGCF_GENERICFRAMEWORK_CORE_GFWCONFIG_H_

#include "GFW.h"

#include "Framework/Logger.h"

#include <TMath.h>
#include <TObject.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace o2
{

namespace analysis::genericframework
{

template <typename T>
int CheckSameSize(const std::vector<T>& first)
{
  return first.size();
}

template <typename T, typename... Args>
int CheckSameSize(const std::vector<T>& first, const std::vector<Args>&... rest)
{
  size_t size = first.size();
  bool allSameSize = ((size == rest.size()) && ...);

  return allSameSize ? size : -1;
}

/// \class GFWBinningCuts
/// \brief Class which implements configurable acceptance cuts
///
class GFWBinningCuts
{
 public:
  GFWBinningCuts(int vtxzbins_ = 40, int etabins_ = 16, int phibins_ = 72, int nchbins_ = 300, float nchmin_ = 0, float nchmax_ = 3000,
                 float ptpoimin_ = 0.2, float ptpoimax_ = 10., float ptrefmin_ = 0.2, float ptrefmax_ = 3.,
                 std::vector<double> ptbinning_ = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                                                   0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                                                   1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                                   2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10},
                 std::vector<double> centbinning_ = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                                                     21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0,
                                                     40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0,
                                                     60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0,
                                                     80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0}) : mVtxZbins{vtxzbins_}, mEtabins{etabins_}, mPhibins{phibins_}, mNchbins{nchbins_}, mNchmin{nchmin_}, mNchmax{nchmax_}, mPTpoimin{ptpoimin_}, mPTpoimax{ptpoimax_}, mPTrefmin{ptrefmin_}, mPTrefmax{ptrefmax_}, mPTbinning{std::move(ptbinning_)}, mCentbinning{std::move(centbinning_)} {};

  auto Print() const
  {
    LOGF(info, "Vz: %d | Eta: %d | Phi: %d | Pt POI: %.2f, %.2f | Pt Ref: %.2f, %.2f | Nch: %d, %.1f, %.1f", mVtxZbins, mEtabins, mPhibins, mPTpoimin, mPTpoimax, mPTrefmin, mPTrefmax, mNchbins, mNchmin, mNchmax);
    return;
  }

  void SetVtxZBinning(int vtxbins)
  {
    mVtxZbins = vtxbins;
  }
  const auto& GetVtxZbins() const { return mVtxZbins; }

  void SetPtPOI(float ptpoimin, float ptpoimax)
  {
    mPTpoimin = ptpoimin;
    mPTpoimax = ptpoimax;
  }
  const auto& GetPtPOImin() const { return mPTpoimin; }
  const auto& GetPtPOImax() const { return mPTpoimax; }

  void SetPtBinning(std::vector<double> ptbinning_) { mPTbinning = std::move(ptbinning_); }
  const auto& GetPtBinning() const { return mPTbinning; }

  void SetEtaBinning(int etabins)
  {
    mEtabins = etabins;
  }
  const auto& GetEtaBins() const { return mEtabins; }

  void SetPhiBins(int phibins) { mPhibins = phibins; }
  const auto& GetPhiBins() const { return mPhibins; }

  void SetPtRef(float ptrefmin, float ptrefmax)
  {
    mPTrefmin = ptrefmin;
    mPTrefmax = ptrefmax;
  }
  const auto& GetPtRefMin() const { return mPTrefmin; }
  const auto& GetPtRefMax() const { return mPTrefmax; }

  const float& GetPtMin() const { return std::min(mPTrefmin, mPTpoimin); }
  const float& GetPtMax() const { return std::max(mPTrefmax, mPTpoimax); }

  void SetNchBinning(int nchbins, float nchmin, float nchmax)
  {
    mNchbins = nchbins;
    mNchmin = nchmin;
    mNchmax = nchmax;
  }
  const auto& GetNchBins() const { return mNchbins; }
  const auto& GetNchMin() const { return mNchmin; }
  const auto& GetNchMax() const { return mNchmax; }

  void SetCentBinning(std::vector<double> centbinning_) { mCentbinning = std::move(centbinning_); }
  const auto& GetCentBinning() const { return mCentbinning; }

 private:
  int mVtxZbins;
  int mEtabins;
  int mPhibins;
  int mNchbins;
  float mNchmin;
  float mNchmax;
  float mPTpoimin;
  float mPTpoimax;
  float mPTrefmin;
  float mPTrefmax;
  std::vector<double> mPTbinning;
  std::vector<double> mCentbinning;
  ClassDefNV(GFWBinningCuts, 1);
};

/// \class GFWRegions
/// \brief Class which implements configurable regions for the GFW
///

class GFWRegions
{
 public:
  GFWRegions(std::vector<std::string> names_ = {"refN", "refP", "refFull"}, std::vector<float> etaminvals_ = {-0.8, 0.4, -0.8},
             std::vector<float> etamaxvals_ = {-0.4, 0.8, 0.8},
             std::vector<int> pTDifs_ = {0, 0, 0},
             std::vector<int> bitmasks_ = {1, 1, 1}) : names{std::move(names_)}, etaminvals{std::move(etaminvals_)}, etamaxvals{std::move(etamaxvals_)}, pTDifs{std::move(pTDifs_)}, bitmasks{std::move(bitmasks_)} {};

  auto Print() const
  {
    for (size_t i = 0; i < names.size(); ++i) {
      LOGF(info, "{%s, %.1f, %.1f, %d, %d}", names[i].c_str(), etaminvals[i], etamaxvals[i], pTDifs[i], bitmasks[i]);
    }
    return;
  }

  auto GetSize() const { return CheckSameSize(names, etaminvals, etamaxvals, pTDifs, bitmasks); }

  void SetNames(std::vector<std::string> names_) { names = std::move(names_); }
  const auto& GetNames() const { return names; }

  void SetEtaMin(std::vector<float> etaminvals_) { etaminvals = std::move(etaminvals_); }
  const auto& GetEtaMin() const { return etaminvals; }

  void SetEtaMax(std::vector<float> etamaxvals_) { etamaxvals = std::move(etamaxvals_); }
  const auto& GetEtaMax() const { return etamaxvals; }

  void SetpTDifs(std::vector<int> pTDifs_) { pTDifs = std::move(pTDifs_); }
  const auto& GetpTDifs() const { return pTDifs; }

  void SetBitmasks(std::vector<int> bitmasks_) { bitmasks = std::move(bitmasks_); }
  const auto& GetBitmasks() const { return bitmasks; }

 private:
  std::vector<std::string> names;
  std::vector<float> etaminvals;
  std::vector<float> etamaxvals;
  std::vector<int> pTDifs;
  std::vector<int> bitmasks;
  ClassDefNV(GFWRegions, 1);
};

/// \class GFWCorrConfigs
/// \brief Class which implements configurable correlations for GFW
///
class GFWCorrConfigs
{
 public:
  GFWCorrConfigs(std::vector<std::string> corrs_ = {"refP {2} refN {-2}", "refP {3} refN {-3}", "refP {4} refN {-4}", "refFull {2 -2}",
                                                    "refFull {2 2 -2 -2}"},
                 std::vector<std::string> heads_ = {"ChGap22", "ChGap32", "ChGap42", "ChFull22", "ChFull24"},
                 std::vector<int> pTDifs_ = {0, 0, 0, 0, 0}, std::vector<int> pTCorrMasks_ = {15, 1, 1, 0, 0}) : corrs{std::move(corrs_)}, heads{std::move(heads_)}, pTDifs{std::move(pTDifs_)}, pTCorrMasks{std::move(pTCorrMasks_)} {};

  auto Print() const
  {
    for (size_t i = 0; i < corrs.size(); ++i) {
      LOGF(info, "{%s,%s,%d,%d}", heads[i].c_str(), corrs[i].c_str(), pTDifs[i], pTCorrMasks[i]);
    }
    return;
  }

  auto GetSize() const { return CheckSameSize(corrs, heads, pTDifs); }

  void SetCorrs(std::vector<std::string> corrs_) { corrs = std::move(corrs_); }
  const auto& GetCorrs() const { return corrs; }

  void SetHeads(std::vector<std::string> heads_) { heads = std::move(heads_); }
  const auto& GetHeads() const { return heads; }

  void SetpTDifs(std::vector<int> pTDifs_) { pTDifs = std::move(pTDifs_); }
  const auto& GetpTDifs() const { return pTDifs; }

  void SetpTCorrMasks(std::vector<int> pTCorrMasks_) { pTCorrMasks = std::move(pTCorrMasks_); }
  const auto& GetpTCorrMasks() const { return pTCorrMasks; }

 private:
  std::vector<std::string> corrs;
  std::vector<std::string> heads;
  std::vector<int> pTDifs;
  std::vector<int> pTCorrMasks;
  ClassDefNV(GFWCorrConfigs, 1);
};

} // namespace analysis::genericframework
} // namespace o2
#endif // PWGCF_GENERICFRAMEWORK_CORE_GFWCONFIG_H_
