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
/// \brief Xu Wang <waxu@cern.ch>

#ifndef PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_
#define PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_

#include <cmath>

namespace o2::analysis::hf_correlations
{
enum Region {
  Default = 0,
  Toward,
  Away,
  Transverse
};

template <typename T>
Region getRegion(T deltaPhi)
{
  if (std::abs(deltaPhi) < o2::constants::math::PI / 3.) {
    return Toward;
  } else if (deltaPhi > 2. * o2::constants::math::PI / 3. && deltaPhi < 4. * o2::constants::math::PI / 3.) {
    return Away;
  } else {
    return Transverse;
  }
}

// Find Leading Particle
template <typename TTracks>
int findLeadingParticle(TTracks const& tracks)
{
  auto leadingParticle = tracks.begin();
  for (auto const& track : tracks) {
    if (std::abs(track.dcaXY()) >= 1. || std::abs(track.dcaZ()) >= 1.) {
      continue;
    }
    if (track.pt() > leadingParticle.pt()) {
      leadingParticle = track;
    }
  }
  int leadingIndex = leadingParticle.globalIndex();
  return leadingIndex;
}

// ======= Find Leading Particle for McGen ============
template <typename TMcParticles>
int findLeadingParticleMcGen(TMcParticles const& mcParticles, const float etaTrackMax, const float ptTrackMin)
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
