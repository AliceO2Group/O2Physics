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

#include <vector>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/AnalysisHelpers.h"
#include "TLorentzVector.h"

namespace o2::upgrade
{

/// Function to convert a TLorentzVector into a perfect Track
/// \param charge particle charge (integer)
/// \param particle the particle to convert (TLorentzVector)
/// \param productionVertex where the particle was produced
/// \param o2track the address of the resulting TrackParCov
void convertTLorentzVectorToO2Track(const int charge,
                                    const TLorentzVector particle,
                                    const std::vector<double> productionVertex,
                                    o2::track::TrackParCov& o2track)
{
  std::array<float, 5> params;
  std::array<float, 15> covm = {0.};
  float s, c, x;
  o2::math_utils::sincos(static_cast<float>(particle.Phi()), s, c);
  o2::math_utils::rotateZInv(static_cast<float>(productionVertex[0]), static_cast<float>(productionVertex[1]), x, params[0], s, c);
  params[1] = static_cast<float>(productionVertex[2]);
  params[2] = 0.; // since alpha = phi
  const auto theta = 2. * std::atan(std::exp(-particle.PseudoRapidity()));
  params[3] = 1. / std::tan(theta);
  params[4] = charge / particle.Pt();

  // Initialize TrackParCov in-place
  new (&o2track)(o2::track::TrackParCov)(x, particle.Phi(), params, covm);
}

/// Function to convert a TLorentzVector into a perfect Track
/// \param pdgCode particle pdg
/// \param particle the particle to convert (TLorentzVector)
/// \param productionVertex where the particle was produced
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
void convertTLorentzVectorToO2Track(int pdgCode,
                                    TLorentzVector particle,
                                    std::vector<double> productionVertex,
                                    o2::track::TrackParCov& o2track,
                                    const o2::framework::Service<o2::framework::O2DatabasePDG>& pdg)
{
  const auto pdgInfo = pdg->GetParticle(pdgCode);
  int charge = 0;
  if (pdgInfo != nullptr) {
    charge = pdgInfo->Charge() / 3;
  }
  convertTLorentzVectorToO2Track(charge, particle, productionVertex, o2track);
}

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename McParticleType>
void convertMCParticleToO2Track(McParticleType& particle,
                                o2::track::TrackParCov& o2track,
                                const o2::framework::Service<o2::framework::O2DatabasePDG>& pdg)
{
  static TLorentzVector tlv;
  tlv.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.e());
  convertTLorentzVectorToO2Track(particle.pdgCode(), tlv, {particle.vx(), particle.vy(), particle.vz()}, o2track, pdg);
}

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename McParticleType>
o2::track::TrackParCov convertMCParticleToO2Track(McParticleType& particle,
                                                  const o2::framework::Service<o2::framework::O2DatabasePDG>& pdg)
{
  o2::track::TrackParCov o2track;
  convertMCParticleToO2Track(particle, o2track, pdg);
  return o2track;
}

} // namespace o2::upgrade

#endif // ALICE3_CORE_TRACKUTILITIES_H_
