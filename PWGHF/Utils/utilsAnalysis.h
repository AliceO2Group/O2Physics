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
/// \return true if track passes all cuts
template <typename Track>
bool isItsQualityCuts(Track const& track, int itsNClsFoundMin, float itsChi2PerClusterMax)
{
  if (track.itsNCls() < itsNClsFoundMin) {
    return false;
  }
  if (track.itsChi2NCl() > itsChi2PerClusterMax) {
    return false;
  }
  return true;
}

/// Single-track cut on TPC track properties
/// \param track track that has to satisfy the selection criteria
/// \return true if track passes all cuts
template <typename Track>
bool isTpcQualityCuts(Track const& track, int tpcNClsFoundMin, int tpcNCrossedRowsMin, float tpcNCrossedRowsOverFindableClustersMin, float tpcChi2PerClusterMax)
{
  if (track.tpcNClsFound() < tpcNClsFoundMin) {
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

} // namespace o2::analysis

#endif // PWGHF_UTILS_UTILSANALYSIS_H_
