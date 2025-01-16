// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUtils.h
/// \brief Utilities for the FemtoDream framework
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMUTILS_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMUTILS_H_

#include <vector>
#include <string>
#include "CommonConstants/PhysicsConstants.h"
#include "PWGCF/DataModel/FemtoDerived.h"

namespace o2::analysis::femtoDream
{

/// function for getting the mass of a particle depending on the pdg code
/// \param pdgCode pdg code of the particle
/// \return mass of the particle
inline float getMass(int pdgCode)
{
  // use this function instead of TDatabasePDG to return masses defined in the PhysicsConstants.h header
  // this approach saves a lot of memory and important partilces like deuteron are missing in TDatabasePDG anyway
  float mass = 0;
  // add new particles if necessary here
  switch (std::abs(pdgCode)) {
    case kPiPlus:
      mass = o2::constants::physics::MassPiPlus;
      break;
    case kKPlus:
      mass = o2::constants::physics::MassKPlus;
      break;
    case kProton:
      mass = o2::constants::physics::MassProton;
      break;
    case kLambda0:
      mass = o2::constants::physics::MassLambda;
      break;
    case o2::constants::physics::Pdg::kPhi:
      mass = o2::constants::physics::MassPhi;
      break;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      mass = o2::constants::physics::MassLambdaCPlus;
      break;
    case o2::constants::physics::Pdg::kDeuteron:
      mass = o2::constants::physics::MassDeuteron;
      break;
    case o2::constants::physics::Pdg::kTriton:
      mass = o2::constants::physics::MassTriton;
      break;
    case o2::constants::physics::Pdg::kHelium3:
      mass = o2::constants::physics::MassHelium3;
      break;
    default:
      LOG(fatal) << "PDG code is not suppored";
  }
  return mass;
}

inline int checkDaughterType(o2::aod::femtodreamparticle::ParticleType partType, int motherPDG)
{
  int partOrigin = 0;
  if (partType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
    switch (abs(motherPDG)) {
      case kLambda0:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterLambda;
        break;
      case kSigmaPlus:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterSigmaplus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    } // switch

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kV0) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kV0Child) {
    switch (abs(motherPDG)) {
      case kLambda0:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterLambda;
        break;
      case kSigmaPlus:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterSigmaplus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    } // switch

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
    partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
  }
  return partOrigin;
};

template <typename T, typename R>
inline bool containsNameValuePair(const std::vector<T>& myVector, const std::string& name, R value)
{
  for (const auto& obj : myVector) {
    if (obj.name == name) {
      if (std::abs(static_cast<float>((obj.defaultValue.template get<R>() - value))) < 1e-2) {
        return true; // Found a match
      }
    }
  }
  return false; // No match found
}

template <typename T>
float itsSignal(T const& track)
{
  uint32_t clsizeflag = track.itsClusterSizes();
  auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
  auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
  auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
  auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
  auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
  auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
  auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;
  int numLayers = 7;
  int sumClusterSizes = clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6 + clSizeLayer0;
  float cosLamnda = 1. / std::cosh(track.eta());
  return (static_cast<float>(sumClusterSizes) / numLayers) * cosLamnda;
};

} // namespace o2::analysis::femtoDream
#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMUTILS_H_
