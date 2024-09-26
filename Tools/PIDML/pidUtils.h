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
constexpr double kTOFMissingSignal = -999.0f;
constexpr double kGlobalEtaCut = 0.8f;

template <typename T>
bool trdMissing(const T& track)
{
  return !track.hasTRD();
}

template <typename T>
bool tofMissing(const T& track)
{
  return !track.hasTOF();
}

template <typename T>
bool inPLimit(const T& track, double pLimit)
{
  return track.p() >= pLimit;
}
} // namespace pidml::pidutils

#endif // TOOLS_PIDML_PIDUTILS_H_
