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

/// \file utilsAnalysis.h
/// \brief Utilities for HF analyses

#ifndef PWGHF_UTILS_UTILSANALYSIS_H_
#define PWGHF_UTILS_UTILSANALYSIS_H_

#include <algorithm> // std::upper_bound
#include <iterator>  // std::distance
#include <string>    //std::string
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/PhysicsConstants.h"

namespace o2::analysis
{
/// Finds pT bin in an array.
/// \param bins  array of pT bins
/// \param value  pT
/// \return index of the pT bin
/// \note Accounts for the offset so that pt bin array can be used to also configure a histogram axis.
template <typename T1, typename T2>
int findBin(T1 const& binsPt, T2 value)
{
  if (value < binsPt->front()) {
    return -1;
  }
  if (value >= binsPt->back()) {
    return -1;
  }
  return std::distance(binsPt->begin(), std::upper_bound(binsPt->begin(), binsPt->end(), value)) - 1;
}

/// Single-track cut on DCAxy and DCAz
/// \param binsPt pT bins
/// \param cuts cut configuration
/// \param pt is the prong pT
/// \param dcaXY is the prong dcaXY
/// \param dcaZ is the prong dcaZ
/// \return true if track passes all cuts
template <typename T1, typename T2>
bool isSelectedTrackDca(T1 const& binsPt, T2 const& cuts, const float pt, const float dcaXY, const float dcaZ)
{
  auto binPt = findBin(binsPt, pt);
  if (binPt == -1) {
    return false;
  }
  if (std::abs(dcaXY) < cuts->get(binPt, "min_dcaxytoprimary")) {
    return false; // minimum DCAxy
  }
  if (std::abs(dcaXY) > cuts->get(binPt, "max_dcaxytoprimary")) {
    return false; // maximum DCAxy
  }
  if (std::abs(dcaZ) < cuts->get(binPt, "min_dcaztoprimary")) {
    return false; // minimum DCAz
  }
  if (std::abs(dcaZ) > cuts->get(binPt, "max_dcaztoprimary")) {
    return false; // maximum DCAz
  }
  return true;
}

/// Single-track cut on ITS track properties
/// \param track track that has to satisfy the selection criteria
/// \param itsNClustersFoundMin is the minimum number of ITS clusters
/// \param itsChi2PerClusterMax is the maximum value of chi2 fit over ITS clusters
/// \return true if track passes all cuts
template <typename T>
bool isSelectedTrackItsQuality(T const& track, const int itsNClustersFoundMin, const float itsChi2PerClusterMax)
{
  if (track.itsNCls() < itsNClustersFoundMin) {
    return false;
  }
  if (track.itsChi2NCl() > itsChi2PerClusterMax) {
    return false;
  }
  return true;
}

/// Single-track cut on TPC track properties
/// \param track track that has to satisfy the selection criteria
/// \param tpcNClustersFoundMin is the minimum number of TPC clusters
/// \param tpcNCrossedRowsMin is the minimum number of crossed TPC rows
/// \param tpcNCrossedRowsOverFindableClustersMin is the minimum of TPC CrossedRows/FindableClusters value
/// \param tpcChi2PerClusterMax is the maximum value of chi2 fit over TPC clusters
/// \return true if track passes all cuts
template <typename T>
bool isSelectedTrackTpcQuality(T const& track, const int tpcNClustersFoundMin, const int tpcNCrossedRowsMin, const float tpcNCrossedRowsOverFindableClustersMin, const float tpcChi2PerClusterMax)
{
  if (track.tpcNClsFound() < tpcNClustersFoundMin) {
    return false;
  }
  if (track.tpcNClsCrossedRows() < tpcNCrossedRowsMin) {
    return false;
  }
  if (track.tpcCrossedRowsOverFindableCls() < tpcNCrossedRowsOverFindableClustersMin) {
    return false;
  }
  if (track.tpcChi2NCl() > tpcChi2PerClusterMax) {
    return false;
  }
  return true;
}

struct HfTriggerCuts : o2::framework::ConfigurableGroup {
  std::string prefix = "hfMassCut"; // JSON group name

  static constexpr float defaultDeltaMassPars3Prong[1][2] = {{-0.0025f, 0.0001f}};
  static constexpr float defaultSigmaPars3Prong[1][2] = {{0.00796f, 0.00176f}};
  o2::framework::Configurable<float> nSigma{"nSigma", 2.5, "Number of Sigmas for pT differential mass cut"};
  o2::framework::Configurable<float> pTMaxDeltaMass{"pTMaxDeltaMass", 10., "Max pT to apply delta mass shift to PDG mass value"};
  o2::framework::Configurable<float> pTMaxMassCut{"pTMaxMassCut", 8., "Max pT to apply pT differential cut"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> deltaMassPars3Prong{"deltaMassPars3Prong", {defaultDeltaMassPars3Prong[0], 2, {"constant", "linear"}}, "delta mass parameters for HF 3 prong trigger mass cut"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> sigmaPars3Prong{"sigmaPars3Prong", {defaultSigmaPars3Prong[0], 2, {"constant", "linear"}}, "sigma parameters for HF 3 prong trigger mass cut"};

  /// Mass selection of D+ candidates to build B0 candidates
  /// \param pTrackSameChargeFirst is the first same-charge track momentum
  /// \param pTrackSameChargeFirst is the second same-charge track momentum
  /// \param pTrackSameChargeFirst is the opposite charge track momentum
  /// \param ptD is the pt of the D+ meson candidate
  /// \return true if D+ passes selection
  bool isSelectedDplusInMassRange(const float& invMassDplus, const float& ptD)
  {
    float peakMean = (ptD < pTMaxDeltaMass) ? ((o2::constants::physics::MassDPlus + deltaMassPars3Prong->get("constant")) + deltaMassPars3Prong->get("linear") * ptD) : o2::constants::physics::MassDPlus;
    float peakWidth = sigmaPars3Prong->get("constant") + sigmaPars3Prong->get("linear") * ptD;

    if (std::fabs(invMassDplus - peakMean) > nSigma * peakWidth && ptD < pTMaxMassCut) {
      return false;
    }
    return true;
  }
};
} // namespace o2::analysis

#endif // PWGHF_UTILS_UTILSANALYSIS_H_
