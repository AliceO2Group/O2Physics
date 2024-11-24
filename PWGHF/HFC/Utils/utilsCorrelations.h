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

/// \file utilsCorrelations.h
/// \brief Utilities for HFC analyses
/// \author Xu Wang <waxu@cern.ch>

#ifndef PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_
#define PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_

#include <cmath>
#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"

namespace o2::analysis::hf_correlations
{
enum Region {
  Default = 0,
  Toward,
  Away,
  Transverse
};

template <typename T>
Region getRegion(T const deltaPhi)
{
  if (std::abs(deltaPhi) < o2::constants::math::PIThird) {
    return Toward;
  } else if (deltaPhi > 2. * o2::constants::math::PIThird && deltaPhi < 4. * o2::constants::math::PIThird) {
    return Away;
  } else {
    return Transverse;
  }
}

// ========= Find Leading Particle ==============
template <typename TTracks, typename T1, typename T2, typename T3>
int findLeadingParticle(TTracks const& tracks, T1 const dcaXYTrackMax, T2 const dcaZTrackMax, T3 const etaTrackMax)
{
  auto leadingParticle = tracks.begin();
  for (auto const& track : tracks) {
    if (std::abs(track.dcaXY()) >= dcaXYTrackMax || std::abs(track.dcaZ()) >= dcaZTrackMax) {
      continue;
    }
    if (std::abs(track.eta()) > etaTrackMax) {
      continue;
    }
    if (track.pt() > leadingParticle.pt()) {
      leadingParticle = track;
    }
  }
  return leadingParticle.globalIndex();
}

// ======= Find Leading Particle for McGen ============
template <typename TMcParticles, typename T1, typename T2>
int findLeadingParticleMcGen(TMcParticles const& mcParticles, T1 const etaTrackMax, T2 const ptTrackMin)
{
  auto leadingParticle = mcParticles.begin();
  for (auto const& mcParticle : mcParticles) {
    if (std::abs(mcParticle.eta()) > etaTrackMax) {
      continue;
    }
    if (mcParticle.pt() < ptTrackMin) {
      continue;
    }
    if ((std::abs(mcParticle.pdgCode()) != kElectron) && (std::abs(mcParticle.pdgCode()) != kMuonMinus) && (std::abs(mcParticle.pdgCode()) != kPiPlus) && (std::abs(mcParticle.pdgCode()) != kKPlus) && (std::abs(mcParticle.pdgCode()) != kProton)) {
      continue;
    }
    if (mcParticle.pt() > leadingParticle.pt()) {
      leadingParticle = mcParticle;
    }
  }
  return leadingParticle.globalIndex();
}
} // namespace o2::analysis::hf_correlations
#endif // PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_
