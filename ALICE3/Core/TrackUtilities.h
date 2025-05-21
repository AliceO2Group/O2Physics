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
///
/// \file TrackUtilities.h
///
/// \brief Set of utilities for the ALICE3 track handling
///
/// \since  May 21, 2025
///

#ifndef ALICE3_CORE_TRACKUTILITIES_H_
#define ALICE3_CORE_TRACKUTILITIES_H_

#include "ReconstructionDataFormats/Track.h"

namespace o2::upgrade
{

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
template <typename McParticleType>
void convertMCParticleToO2Track(McParticleType& particle,
                                o2::track::TrackParCov& o2track,
                                const o2::framework::Service<o2::framework::O2DatabasePDG>& pdg)
{
  auto pdgInfo = pdg->GetParticle(particle.pdgCode());
  int charge = 0;
  if (pdgInfo != nullptr) {
    charge = pdgInfo->Charge() / 3;
  }
  std::array<float, 5> params;
  std::array<float, 15> covm = {0.};
  float s, c, x;
  o2::math_utils::sincos(particle.phi(), s, c);
  o2::math_utils::rotateZInv(particle.vx(), particle.vy(), x, params[0], s, c);
  params[1] = particle.vz();
  params[2] = 0.; // since alpha = phi
  auto theta = 2. * std::atan(std::exp(-particle.eta()));
  params[3] = 1. / std::tan(theta);
  params[4] = charge / particle.pt();

  // Initialize TrackParCov in-place
  new (&o2track)(o2::track::TrackParCov)(x, particle.phi(), params, covm);
}

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
template <typename McParticleType>
o2::track::TrackParCov convertMCParticleToO2Track(McParticleType& particle,
                                                  const o2::framework::Service<o2::framework::O2DatabasePDG>& pdg)
{
  TrackParCov o2track;
  convertMCParticleToO2Track(particle, o2track, pdg);
}

} // namespace o2::upgrade

#endif // ALICE3_CORE_TRACKUTILITIES_H_
