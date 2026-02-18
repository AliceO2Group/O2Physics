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

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>
#include <MathUtils/BetheBlochAleph.h>

#include <cstdint>

using namespace o2::constants::physics;
using namespace o2::framework;

namespace o2::aod::pid_tpc_tof_utils
{
/// @brief Species of HF-candidate daughter tracks
enum HfProngSpecies : uint8_t {
  Pion = 0,
  Kaon,
  Proton,
  Deuteron,
  Triton,
  Helium3,
  Alpha,
  NHfProngSpecies
};

/// @brief PID methods used for HF-candidate daughter tracks
enum PidMethod {
  NoPid = 0, // none
  TpcOrTof,  // TPC or TOF
  TpcAndTof, // TPC and TOF
  NPidMethods
};

/// Compute TPC nσ for light nuclei (De/Tr/He3/Al) using a Bethe–Bloch parameter configuration (BB-based PID).
///
/// \tparam TrackType   Track/ASoA row type providing TPC accessors.
/// \param track        Track to be tested.
/// \param lightnuclei  Species selector: 3=Deuteron, 4=Triton, 5=Helium3, 6=Alpha.
/// \param bbParams is Bethe–Bloch  parameters
/// \return             TPC nσ for the chosen nucleus hypothesis (or -999 if not applicable).
template <typename TrackType>
float getTPCNSigmaLightNucleiBetheBloch(const TrackType& track, HfProngSpecies lightnuclei,
                                        const Configurable<LabeledArray<float>>& bbParams)
{
  if (!track.hasTPC()) {
    return -999.f;
  }
  const int row = static_cast<int>(lightnuclei) - static_cast<int>(HfProngSpecies::Deuteron);

  if (row < 0 || row >= HfProngSpecies::NHfProngSpecies) {
    return -999.f;
  }

  // Columns: [0..4] BB params, [5] relative resolution (sigma/mean)
  const double bb0 = bbParams->get(row, 0u);
  const double bb1 = bbParams->get(row, 1u);
  const double bb2 = bbParams->get(row, 2u);
  const double bb3 = bbParams->get(row, 3u);
  const double bb4 = bbParams->get(row, 4u);
  const double relRes = bbParams->get(row, 5u);

  if (relRes <= 0.) {
    return -999.f;
  }

  double mass = 0.;
  switch (lightnuclei) {
    case HfProngSpecies::Deuteron:
      mass = MassDeuteron;
      break;
    case HfProngSpecies::Triton:
      mass = MassTriton;
      break;
    case HfProngSpecies::Helium3:
      mass = MassHelium3;
      break;
    case HfProngSpecies::Alpha:
      mass = MassAlpha;
      break;
    default:
      LOG(fatal) << "Unhandled HfProngSpecies " << static_cast<int>(lightnuclei);
  }

  const int charge = (lightnuclei == HfProngSpecies::Helium3 || lightnuclei == HfProngSpecies::Alpha) ? 2 : 1;
  const float rigidity = track.tpcInnerParam(); // p/|q|

  const double x = static_cast<double>(charge) * static_cast<double>(rigidity) / mass;
  const double expBethe = common::BetheBlochAleph(x, bb0, bb1, bb2, bb3, bb4);
  const double expSigma = expBethe * relRes;

  if (expSigma <= 0.) {
    return -999.f;
  }

  return static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
}

/// \brief Function to fill tables with HF prong PID information
/// \tparam specPid particle species
/// \tparam TTrack datatype of the prong track
/// \tparam TCursor datatype of the cursor of the prong PID table to fill
/// \param track prong track
/// \param bbParams is Bethe–Bloch  parameters (only for light nuclei)
/// \param rowPid cursor of the prong PID table to fill
template <HfProngSpecies SpecPid, typename TTrack, typename TCursor>
void fillProngPidLightNuclei(TTrack const& track, TCursor& rowPid,
                             const Configurable<LabeledArray<float>>& bbParams)
{

  // get PID information for the daughter tracks
  // TODO: add here the code for a possible PID post-calibrations in MC
  float nSigTpc = -999.f;
  float nSigTof = -999.f;
  if constexpr (SpecPid == HfProngSpecies::Deuteron) {
    // deuteron PID
    if (track.hasTPC()) {
      nSigTpc = getTPCNSigmaLightNucleiBetheBloch(track, SpecPid, bbParams);
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaDe();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Triton) {
    // triton PID
    if (track.hasTPC()) {
      nSigTpc = getTPCNSigmaLightNucleiBetheBloch(track, SpecPid, bbParams);
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaTr();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Helium3) {
    // helium3 PID
    if (track.hasTPC()) {
      nSigTpc = getTPCNSigmaLightNucleiBetheBloch(track, SpecPid, bbParams);
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaHe();
    }
  } else if constexpr (SpecPid == HfProngSpecies::Alpha) {
    // helium4 PID
    if (track.hasTPC()) {
      nSigTpc = getTPCNSigmaLightNucleiBetheBloch(track, SpecPid, bbParams);
    }
    if (track.hasTOF()) {
      nSigTof = track.tofNSigmaAl();
    }
  } else {
    LOG(fatal) << "Unsupported PID. Supported species in HF framework:  HfProngSpecies::Deuteron, HfProngSpecies::Triton, HfProngSpecies::Helium3, HfProngSpecies::Alpha";
  }

  // fill candidate prong PID rows
  rowPid(nSigTpc, nSigTof);
}

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
  } else {
    LOG(fatal) << "Unsupported PID. Supported species in HF framework: HfProngSpecies::Pion, HfProngSpecies::Kaon, HfProngSpecies::Proton";
  }

  // fill candidate prong PID rows
  rowPid(nSigTpc, nSigTof);
}
} // namespace o2::aod::pid_tpc_tof_utils

#endif // PWGHF_UTILS_UTILSPID_H_
