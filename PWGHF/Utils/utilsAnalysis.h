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

namespace o2::analysis
{
/// Finds pt bin in an array.
/// \param bins  array of pt bins
/// \param value  pt
/// \return index of the pt bin
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
/// \param binsPt pt bins
/// \param cuts cut configuration
/// \param pt is the prong pt
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

/// Configurable group to apply trigger specific cuts for HF analysis
struct HfTriggerCuts : o2::framework::ConfigurableGroup {
  std::string prefix = "hfTriggerCuts"; // JSON group name

  static constexpr float defaultDeltaMassPars3Prong[1][2] = {{-0.0025f, 0.0001f}};
  static constexpr float defaultSigmaPars3Prong[1][2] = {{0.00796f, 0.00176f}};
  static constexpr float defaultDeltaMassPars2Prong[1][2] = {{-0.0025f, 0.0001f}};
  static constexpr float defaultSigmaPars2Prong[1][2] = {{0.01424f, 0.00178f}};
  o2::framework::Configurable<float> nSigma3ProngMax{"nSigma3ProngMax", 2, "Maximum number of sigmas for pT-differential mass cut for 3-prong candidates"};
  o2::framework::Configurable<float> nSigma2ProngMax{"nSigma2ProngMax", 2, "Maximum number of sigmas for pT-differential mass cut for 2-prong candidates"};
  o2::framework::Configurable<float> ptDeltaMass3ProngMax{"ptDeltaMass3ProngMax", 10., "Max pT to apply delta mass shift to PDG mass value for 3-prong candidates"};
  o2::framework::Configurable<float> ptDeltaMass2ProngMax{"ptDeltaMass2ProngMax", 10., "Max pT to apply delta mass shift to PDG mass value for 2-prong candidates"};
  o2::framework::Configurable<float> ptMassCut3ProngMax{"ptMassCut3ProngMax", 8., "Max pT to apply pT-differential cut for 3-prong candidates"};
  o2::framework::Configurable<float> ptMassCut2ProngMax{"ptMassCut2ProngMax", 8., "Max pT to apply pT-differential cut for 2-prong candidates"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> deltaMassPars3Prong{"deltaMassPars3Prong", {defaultDeltaMassPars3Prong[0], 2, {"constant", "linear"}}, "delta mass parameters for HF 3-prong trigger mass cut"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> deltaMassPars2Prong{"deltaMassPars2Prong", {defaultDeltaMassPars2Prong[0], 2, {"constant", "linear"}}, "delta mass parameters for HF 2-prong trigger mass cut"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> sigmaPars3Prong{"sigmaPars3Prong", {defaultSigmaPars3Prong[0], 2, {"constant", "linear"}}, "sigma parameters for HF 3-prong trigger mass cut"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> sigmaPars2Prong{"sigmaPars2Prong", {defaultSigmaPars2Prong[0], 2, {"constant", "linear"}}, "sigma parameters for HF 2-prong trigger mass cut"};

  /// Mass selection of 2 or 3 prong canidates in triggered data analysis
  /// \param invMass is the invariant mass of the candidate
  /// \param pdgMass is the pdg Mass of the candidate particle
  /// \param pt is the pt of the candidate
  /// \return true if candidate passes selection
  template <bool is3Prong>
  bool isCandidateInMassRange(const float& invMass, const float& pdgMass, const float& pt)
  {
    float peakMean{0.};
    float peakWidth{0.};
    float ptMassCutMax{0.};
    if constexpr (is3Prong) {
      peakMean = (pt < ptDeltaMass3ProngMax) ? ((pdgMass + deltaMassPars3Prong->get("constant")) + deltaMassPars3Prong->get("linear") * pt) : pdgMass;
      peakWidth = sigmaPars3Prong->get("constant") + sigmaPars3Prong->get("linear") * pt;
      ptMassCutMax = ptMassCut3ProngMax;
    } else {
      float peakMean = (pt < ptDeltaMass2ProngMax) ? ((pdgMass + deltaMassPars2Prong->get("constant")) + deltaMassPars2Prong->get("linear") * pt) : pdgMass;
      float peakWidth = sigmaPars2Prong->get("constant") + sigmaPars2Prong->get("linear") * pt;
      ptMassCutMax = ptMassCut2ProngMax;
    }
    return (!(std::fabs(invMass - peakMean) > nSigma3ProngMax * peakWidth && pt < ptMassCutMax));
  }
};
} // namespace o2::analysis

#endif // PWGHF_UTILS_UTILSANALYSIS_H_
