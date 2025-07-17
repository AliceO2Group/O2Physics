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
/// \file   TrackUtilities.h
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
/// \brief  Set of utilities for the ALICE3 track handling
/// \since  May 21, 2025
///

#ifndef ALICE3_CORE_TRACKUTILITIES_H_
#define ALICE3_CORE_TRACKUTILITIES_H_

#include "ReconstructionDataFormats/Track.h"

#include "TLorentzVector.h"

#include <vector>

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
                                    o2::track::TrackParCov& o2track);

/// Function to convert a TLorentzVector into a perfect Track
/// \param pdgCode particle pdg
/// \param particle the particle to convert (TLorentzVector)
/// \param productionVertex where the particle was produced
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename PdgService>
void convertTLorentzVectorToO2Track(int pdgCode,
                                    TLorentzVector particle,
                                    std::vector<double> productionVertex,
                                    o2::track::TrackParCov& o2track,
                                    const PdgService& pdg)
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
template <typename McParticleType, typename PdgService>
void convertMCParticleToO2Track(McParticleType& particle,
                                o2::track::TrackParCov& o2track,
                                const PdgService& pdg)
{
  static TLorentzVector tlv;
  tlv.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.e());
  convertTLorentzVectorToO2Track(particle.pdgCode(), tlv, {particle.vx(), particle.vy(), particle.vz()}, o2track, pdg);
}

/// Function to convert a McParticle into a perfect Track
/// \param particle the particle to convert (mcParticle)
/// \param o2track the address of the resulting TrackParCov
/// \param pdg the pdg service
template <typename McParticleType, typename PdgService>
o2::track::TrackParCov convertMCParticleToO2Track(McParticleType& particle,
                                                  const PdgService& pdg)
{
  o2::track::TrackParCov o2track;
  convertMCParticleToO2Track(particle, o2track, pdg);
  return o2track;
}

} // namespace o2::upgrade

#endif // ALICE3_CORE_TRACKUTILITIES_H_
