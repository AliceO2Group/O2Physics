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

/// \file utilsPid.h
/// \brief PID utilities for HF analyses
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, CERN

#ifndef PWGHF_UTILS_UTILSPID_H_
#define PWGHF_UTILS_UTILSPID_H_

#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <Framework/Logger.h>

#include <cstdint>

namespace o2::aod::pid_tpc_tof_utils
{
/// @brief Species of HF-candidate daughter tracks
enum HfProngSpecies : uint8_t {
  Pion = 0,
  Kaon,
  Proton,
  Deuteron,
  Triton,
  Helium,
  NHfProngSpecies
};

/// @brief PID methods used for HF-candidate daughter tracks
enum PidMethod {
  NoPid = 0, // none
  TpcOrTof,  // TPC or TOF
  TpcAndTof, // TPC and TOF
  NPidMethods
};

/// Function to combine TPC and TOF nSigma
/// \tparam Tiny switch between full and tiny (binned) PID tables
/// \param nSigmaTpc is the (binned) nSigma separation in TPC (if Tiny = true)
/// \param nSigmaTof is the (binned) nSigma separation in TOF (if Tiny = true)
/// \return combined nSigma of TPC and TOF
template <bool Tiny, typename TNumber>
TNumber combineNSigma(TNumber nSigmaTpc, TNumber nSigmaTof)
{
  static constexpr float NSigmaToleranceDefault = .1f;
  static constexpr float NSigmaDefault = -999.f + NSigmaToleranceDefault; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

  if constexpr (Tiny) {
    nSigmaTpc *= aod::pidtpc_tiny::binning::bin_width;
    nSigmaTof *= aod::pidtof_tiny::binning::bin_width;
  }

  if ((nSigmaTpc > NSigmaDefault) && (nSigmaTof > NSigmaDefault)) { // TPC and TOF
    return std::sqrt(.5f * (nSigmaTpc * nSigmaTpc + nSigmaTof * nSigmaTof));
  }
  if (nSigmaTpc > NSigmaDefault) { // only TPC
    return std::abs(nSigmaTpc);
  }
  if (nSigmaTof > NSigmaDefault) { // only TOF
    return std::abs(nSigmaTof);
  }
  return nSigmaTof; // no TPC nor TOF
}

/// \brief Function to fill tables with HF prong PID information
/// \tparam specPid particle species
/// \tparam TTrack datatype of the prong track
/// \tparam TCursor datatype of the cursor of the prong PID table to fill
/// \param track prong track
/// \param rowPid cursor of the prong PID table to fill
template <HfProngSpecies SpecPid, typename TTrack, typename TCursor>
void fillProngPid(TTrack const& track, TCursor& rowPid)
{

  // get PID information for the daughter tracks
  // TODO: add here the code for a possible PID post-calibrations in MC
  float nSigTpc = -999.f;
  float nSigTof = -999.f;
  if constexpr (SpecPid == HfProngSpecies::Pion) {
    // pion PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaPi();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaPi();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Kaon) {
    // kaon PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaKa();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaKa();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Proton) {
    // proton PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaPr();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaPr();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Deuteron) {
    // deuteron PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaDe();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaDe();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Triton) {
    // triton PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaTr();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaTr();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Helium) {
    // triton PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaHe();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaHe();
    }
  } else {
    LOG(fatal) << "Unsupported PID. Supported species in HF framework: HfProngSpecies::Pion, HfProngSpecies::Kaon, HfProngSpecies::Proton, HfProngSpecies::Deuteron, HfProngSpecies::Triton, HfProngSpecies::Helium";
  }

  // fill candidate prong PID rows
  rowPid(nSigTpc, nSigTof);
}
} // namespace o2::aod::pid_tpc_tof_utils

#endif // PWGHF_UTILS_UTILSPID_H_
