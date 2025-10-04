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
///
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_UTILS_UTILSANALYSIS_H_
#define PWGHF_UTILS_UTILSANALYSIS_H_

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>

#include <algorithm> // std::upper_bound
#include <cstdlib>
#include <iterator> // std::distance
#include <string>   //std::string

namespace o2::analysis
{
enum BHadMothers { NotMatched = 0,
                   BPlus,
                   BZero,
                   Bs,
                   LambdaBZero };

/// Convert the B hadron mother PDG for non prompt candidates to a flag
/// \param flagBHad pdg of the b hadron mother
/// \return integer map to specific mothers' PDG codes
inline BHadMothers getBHadMotherFlag(const int flagBHad)
{
  if (std::abs(flagBHad) == o2::constants::physics::kBPlus) {
    return BHadMothers::BPlus;
  }
  if (std::abs(flagBHad) == o2::constants::physics::kB0) {
    return BHadMothers::BZero;
  }
  if (std::abs(flagBHad) == o2::constants::physics::kBS) {
    return BHadMothers::Bs;
  }
  if (std::abs(flagBHad) == o2::constants::physics::kLambdaB0) {
    return BHadMothers::LambdaBZero;
  }
  return BHadMothers::NotMatched;
}

/// Finds bin in an array that contains a value.
/// \param bins  array of bins
/// \param value  value to find
/// \return index of the pT bin
/// \note Accounts for the offset so that pT bin array can be used to also configure a histogram axis.
template <typename TArrayPt, typename TNumber>
int findBin(TArrayPt const& bins, TNumber const value)
{
  if (value < bins->front()) {
    return -1;
  }
  if (value >= bins->back()) {
    return -1;
  }
  return std::distance(bins->begin(), std::upper_bound(bins->begin(), bins->end(), value)) - 1;
}

/// Single-track cut on DCAxy and DCAz
/// \param binsPt pT bins
/// \param cuts cut configuration
/// \param pt is the prong pT
/// \param dcaXY is the prong dcaXY
/// \param dcaZ is the prong dcaZ
/// \return true if track passes all cuts
template <typename TArrayPt, typename TArrayCuts>
bool isSelectedTrackDca(TArrayPt const& binsPt,
                        TArrayCuts const& cuts,
                        const float pt,
                        const float dcaXY,
                        const float dcaZ)
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
template <typename TTrack>
bool isSelectedTrackItsQuality(TTrack const& track,
                               const int itsNClustersFoundMin,
                               const float itsChi2PerClusterMax)
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
template <typename TTrack>
bool isSelectedTrackTpcQuality(TTrack const& track,
                               const int tpcNClustersFoundMin,
                               const int tpcNCrossedRowsMin,
                               const float tpcNCrossedRowsOverFindableClustersMin,
                               const float tpcChi2PerClusterMax)
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

/// Mass selection of 2 or 3 prong candidates in triggered data analysis
/// \tparam nProngs switch between 2-prong and 3-prong selection
/// \param invMass is the invariant mass of the candidate
/// \param pdgMass is the pdg Mass of the candidate particle
/// \param pt is the pT of the candidate
/// \param cutConfig is the struct with the pt-dependent mass configurations
/// \return true if candidate passes selection
template <typename Config>
bool isCandidateInMassRange(const float invMass,
                            const double pdgMass,
                            const float pt,
                            Config const& cutConfig)
{
  const float peakMean = (pt < cutConfig.ptDeltaMassMax.value) ? ((pdgMass + cutConfig.deltaMassPars->get("constant")) + cutConfig.deltaMassPars->get("linear") * pt) : pdgMass;
  const float peakWidth = cutConfig.sigmaPars->get("constant") + cutConfig.sigmaPars->get("linear") * pt;

  return (!(std::abs(invMass - peakMean) > cutConfig.nSigmaMax.value * peakWidth && pt < cutConfig.ptMassCutMax.value));
}

/// Configurable group to apply trigger specific cuts for 2-prong HF analysis
struct HfTrigger2ProngCuts : o2::framework::ConfigurableGroup {
  std::string prefix = "hfTrigger2ProngCuts"; // JSON group name

  static constexpr float DefaultDeltaMassPars[1][2] = {{-0.0025f, 0.0001f}};
  static constexpr float DefaultSigmaPars[1][2] = {{0.01424f, 0.00178f}};
  o2::framework::Configurable<float> nSigmaMax{"nSigmaMax", 2.f, "Maximum number of sigmas for pT-differential mass cut for 2-prong candidates"};
  o2::framework::Configurable<float> ptDeltaMassMax{"ptDeltaMassMax", 10.f, "Max pT to apply delta mass shift to PDG mass value for 2-prong candidates"};
  o2::framework::Configurable<float> ptMassCutMax{"ptMassCutMax", 9999.f, "Max pT to apply pT-differential cut for 2-prong candidates"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> deltaMassPars{"deltaMassPars", {DefaultDeltaMassPars[0], 2, {"constant", "linear"}}, "delta mass parameters for HF 2-prong trigger mass cut"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> sigmaPars{"sigmaPars", {DefaultSigmaPars[0], 2, {"constant", "linear"}}, "sigma parameters for HF 2-prong trigger mass cut"};
};

/// Configurable group to apply trigger specific cuts for 3-prong HF analysis
struct HfTrigger3ProngCuts : o2::framework::ConfigurableGroup {
  std::string prefix = "hfTrigger3ProngCuts"; // JSON group name

  static constexpr float DefaultDeltaMassPars[1][2] = {{-0.0025f, 0.0001f}};
  static constexpr float DefaultSigmaPars[1][2] = {{0.00796f, 0.00176f}};
  o2::framework::Configurable<float> nSigmaMax{"nSigmaMax", 2, "Maximum number of sigmas for pT-differential mass cut for 3-prong candidates"};
  o2::framework::Configurable<float> ptDeltaMassMax{"ptDeltaMassMax", 10.f, "Max pT to apply delta mass shift to PDG mass value for 3-prong candidates"};
  o2::framework::Configurable<float> ptMassCutMax{"ptMassCutMax", 9999.f, "Max pT to apply pT-differential cut for 3-prong candidates"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> deltaMassPars{"deltaMassPars", {DefaultDeltaMassPars[0], 2, {"constant", "linear"}}, "delta mass parameters for HF 3-prong trigger mass cut"};
  o2::framework::Configurable<o2::framework::LabeledArray<float>> sigmaPars{"sigmaPars", {DefaultSigmaPars[0], 2, {"constant", "linear"}}, "sigma parameters for HF 3-prong trigger mass cut"};
};

} // namespace o2::analysis

#endif // PWGHF_UTILS_UTILSANALYSIS_H_
