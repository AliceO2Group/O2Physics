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

/// \file pidUtils.h
/// \brief Common constants and functions for PID.
///
/// \author Marek Mytkowski <marek.mytkowski@cern.ch>

#ifndef TOOLS_PIDML_PIDUTILS_H_
#define TOOLS_PIDML_PIDUTILS_H_

#include <cstdlib>

namespace pidml::pidutils
{
// magic number had been changed for Run3
constexpr double kRun2TRDMissingSignal = 0.0f;
constexpr double kRun3TRDMissingSignal = -999.0f;
constexpr double kTOFMissingSignal = -999.0f;
constexpr double kTOFMissingBeta = -999.0f;
constexpr double kEpsilon = 1e-6f;

bool almostEqual(double a, double b, double eps = kEpsilon)
{
  return std::abs(a - b) <= eps;
}

template <typename T>
bool trdMissing(const T& track, bool isRun3 = true)
{
  return almostEqual(track.trdSignal(), isRun3 ? kRun3TRDMissingSignal : kRun2TRDMissingSignal);
}

template <typename T>
bool tofMissing(const T& track)
{
  // Because of run3 data we use also TOF beta value to determine if signal is present
  return almostEqual(track.tofSignal(), kTOFMissingSignal, kEpsilon) || almostEqual(track.beta(), kTOFMissingBeta);
}

template <typename T>
bool inPLimit(const T& track, double pLimit)
{
  return track.p() >= pLimit;
}
} // namespace pidml::pidutils

#endif // TOOLS_PIDML_PIDUTILS_H_
