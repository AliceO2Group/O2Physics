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

namespace o2::aod::pid_tpc_tof_utils
{
enum HfProngSpecies : uint8_t { Pion = 0,
                                Kaon,
                                Proton };

/// Function to combine TPC and TOF NSigma
/// \param tiny switch between full and tiny (binned) PID tables
/// \param tpcNSigma is the (binned) NSigma separation in TPC (if tiny = true)
/// \param tofNSigma is the (binned) NSigma separation in TOF (if tiny = true)
/// \return combined NSigma of TPC and TOF
template <bool tiny, typename T1>
T1 combineNSigma(T1 tpcNSigma, T1 tofNSigma)
{
  static constexpr float DefaultNSigmaTolerance = .1f;
  static constexpr float DefaultNSigma = -999.f + DefaultNSigmaTolerance; // -999.f is the default value set in TPCPIDResponse.h and PIDTOF.h

  if constexpr (tiny) {
    tpcNSigma *= aod::pidtpc_tiny::binning::bin_width;
    tofNSigma *= aod::pidtof_tiny::binning::bin_width;
  }

  if ((tpcNSigma > DefaultNSigma) && (tofNSigma > DefaultNSigma)) { // TPC and TOF
    return std::sqrt(.5f * (tpcNSigma * tpcNSigma + tofNSigma * tofNSigma));
  }
  if (tpcNSigma > DefaultNSigma) { // only TPC
    return std::abs(tpcNSigma);
  }
  if (tofNSigma > DefaultNSigma) { // only TOF
    return std::abs(tofNSigma);
  }
  return tofNSigma; // no TPC nor TOF
}

/// @brief Function to fill tables with HF prong PID information
/// @tparam TRK datatype of the prong track
/// @tparam ROW datatype of the prong PID table to fill
/// @tparam specPid particle species
/// @param track prong track
/// @param rowPid cursor of the prong PID table to fill
template <HfProngSpecies specPid, typename TRK, typename ROW>
void fillProngPid(TRK const& track, ROW& rowPid)
{

  // get PID information for the daughter tracks
  // TODO: add here the code for a possible PID post-calibrations in MC
  float nSigTpc = -999.f;
  float nSigTof = -999.f;
  if constexpr (specPid == HfProngSpecies::Pion) {
    // pion PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaPi();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaPi();
    }
  } else if constexpr (specPid == HfProngSpecies::Kaon) {
    // kaon PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaKa();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaKa();
    }
  } else if constexpr (specPid == HfProngSpecies::Proton) {
    // proton PID
    if (track.hasTPC()) {
      nSigTpc = track.tpcNSigmaPr();
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaPr();
    }
  } else {
    LOG(fatal) << "Unsupported PID. Supported species in HF framework: HfProngSpecies::Pion, HfProngSpecies::Kaon, HfProngSpecies::Proton";
  }

  // fill candidate prong PID rows
  rowPid(nSigTpc, nSigTof);
}
} // namespace o2::aod::pid_tpc_tof_utils

#endif // PWGHF_UTILS_UTILSPID_H_
