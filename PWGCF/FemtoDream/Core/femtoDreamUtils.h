// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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

#include "PWGCF/DataModel/FemtoDerived.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>

#include <TMCProcess.h>
#include <TPDGCode.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

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
    case kXiMinus:
      mass = o2::constants::physics::MassXiMinus;
      break;
    case o2::constants::physics::Pdg::kPhi:
      mass = o2::constants::physics::MassPhi;
      break;
    case o2::constants::physics::Pdg::kD0:
      mass = o2::constants::physics::MassD0;
      break;
    case o2::constants::physics::Pdg::kDStar:
      mass = o2::constants::physics::MassDStar;
      break;
    case o2::constants::physics::Pdg::kDPlus:
      mass = o2::constants::physics::MassDPlus;
      break;
    case o2::constants::physics::Pdg::kDS:
      mass = o2::constants::physics::MassDS;
      break;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      mass = o2::constants::physics::MassLambdaCPlus;
      break;
    case o2::constants::physics::Pdg::kXiCPlus:
      mass = o2::constants::physics::MassXiCPlus;
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
    // case o2::constants::physics::Pdg::kOmegaMinus:
    case kOmegaMinus:
      mass = o2::constants::physics::MassOmegaMinus;
      break;
    case o2::constants::physics::Pdg::kK0Star892:
      mass = o2::constants::physics::MassK0Star892;
      break;
    case 310: /// K0Short is not implemented in o2::physics::constants::Pdg
      mass = o2::constants::physics::MassK0Short;
      break;
    default:
      LOG(fatal) << "PDG code is not supported";
  }
  return mass;
}

/// Identify conventional charm/beauty mesons and baryons from their PDG quark digits.
inline bool isHeavyFlavorHadron(int pdgCode)
{
  constexpr int DigitBase = 10;
  constexpr int HadronLimit = 10000000;
  const auto code = std::abs(pdgCode);
  if (code >= HadronLimit || code % DigitBase == 0) {
    return false;
  }
  const auto quark1 = (code / RecoDecay::PdgDivisorBaryon) % DigitBase;
  const auto quark2 = (code / RecoDecay::PdgDivisorMeson) % DigitBase;
  const auto quark3 = (code / DigitBase) % DigitBase;
  // Quarks, diquarks and generator-specific entries are not hadron decays.
  if (quark1 > kTop || quark2 < kDown || quark2 > kTop || quark3 < kDown || quark3 > kTop) {
    return false;
  }
  return quark1 == kCharm || quark1 == kBottom || quark2 == kCharm || quark2 == kBottom || quark3 == kCharm || quark3 == kBottom;
}

template <typename T>
inline bool hasHeavyFlavorHadronAncestor(T const& particles, typename T::iterator const& particle)
{
  auto ancestor = particle;
  std::vector<int64_t> motherIds;
  while (ancestor.has_mothers()) {
    motherIds.clear();
    RecoDecay::getOrderedMotherIds(ancestor, motherIds);
    if (motherIds.empty()) {
      break;
    }
    ancestor = particles.rawIteratorAt(motherIds.front() - particles.offset());
    if (isHeavyFlavorHadron(ancestor.pdgCode())) {
      return true;
    }
  }
  return false;
}

/// Classify true kaons from the correct collision for the DCA templates.
template <typename T>
inline int getKaonOrigin(T const& particles, typename T::iterator const& particle)
{
  using namespace o2::aod::femtodreamMCparticle;
  constexpr int GenFromTransport = -1;
  if (particle.getProcess() == TMCProcess::kPHInhelastic && particle.getGenStatusCode() == GenFromTransport) {
    return kMaterial;
  }
  // HF decay products can be physical primaries. Check their ancestry first.
  if (hasHeavyFlavorHadronAncestor(particles, particle)) {
    return kSecondaryFromHF;
  }
  if (particle.isPhysicalPrimary()) {
    return kPrimary;
  }
  if (particle.getProcess() == TMCProcess::kPDecay && particle.getGenStatusCode() == GenFromTransport && particle.has_mothers()) {
    return kSecondaryFromLF;
  }
  return kElse;
}

inline int checkDaughterType(o2::aod::femtodreamparticle::ParticleType partType, int motherPDG, int daughterPDG = 0)
{
  int partOrigin = 0;
  const auto absMotherPDG = std::abs(motherPDG);
  const auto absDaughterPDG = std::abs(daughterPDG);
  const bool isTrackLike = partType == o2::aod::femtodreamparticle::ParticleType::kTrack ||
                           partType == o2::aod::femtodreamparticle::ParticleType::kV0Child ||
                           partType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child ||
                           partType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor;

  if (absDaughterPDG == kKPlus && isTrackLike) {
    switch (absMotherPDG) {
      case kOmegaMinus:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterOmegaMinus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    }
    return partOrigin;
  }

  if (partType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
    switch (absMotherPDG) {
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
    switch (absMotherPDG) {
      case kSigma0:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterSigma0;
        break;
      case kXiMinus:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterXiMinus;
        break;
      case o2::constants::physics::Pdg::kXi0:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterXi0;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    }

  } else if (partType == o2::aod::femtodreamparticle::ParticleType::kV0Child || partType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child || partType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
    switch (absMotherPDG) {
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
    switch (absMotherPDG) {
      case kOmegaMinus:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterOmegaMinus;
        break;
      case 3324:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterXistar0;
        break;
      case 3314:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondaryDaughterXistarMinus;
        break;
      default:
        partOrigin = aod::femtodreamMCparticle::ParticleOriginMCTruth::kSecondary;
    }
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
