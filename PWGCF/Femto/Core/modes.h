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

/// \file modes.h
/// \brief common modes
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_MODES_H_
#define PWGCF_FEMTO_CORE_MODES_H_

#include "PWGCF/Femto/Core/dataTypes.h"

#include <Rtypes.h>

#include <cstdint>
#include <type_traits>

namespace o2::analysis::femto
{
namespace modes
{

// check if flag is set
template <typename T>
constexpr bool isFlagSet(T value, T flag)
{
  using U = std::underlying_type_t<T>;
  return (static_cast<U>(value) & static_cast<U>(flag)) != 0;
}

// check if flag is equal
template <typename T>
constexpr bool isEqual(T lhs, T rhs)
{
  using U = std::underlying_type_t<T>;
  return static_cast<U>(lhs) == static_cast<U>(rhs);
}

enum class Mode : uint32_t {
  kAnalysis = BIT(0),
  kQa = BIT(1),
  kMc = BIT(2),
  kAnalysis_Qa = kAnalysis | kQa,
  kAnalysis_Mc = kAnalysis | kMc,
  kAnalysis_Qa_Mc = kAnalysis | kQa | kMc,
};

enum class System : uint32_t {
  kPP = BIT(0),
  kPbPb = BIT(1),
  kMC = BIT(2),
  kRun3 = BIT(3),
  kRun2 = BIT(4),
  kPP_Run3 = kPP | kRun3,
  kPP_Run3_MC = kPP | kRun3 | kMC,
  kPP_Run2 = kPP | kRun2,
  kPbPb_Run3 = kPbPb | kRun3,
  kPbPb_Run2 = kPbPb | kRun2,
};

enum class MomentumType : o2::aod::femtodatatypes::MomentumType {
  kPt = 0,    // transverse momentum
  kPAtPv = 1, // momentum at primary vertex
  kPTpc = 2,  // momentum at inner wall of tpc
};

enum class TransverseMassType : o2::aod::femtodatatypes::TransverseMassType {
  kAveragePdgMass = 0,
  kReducedPdgMass = 1,
  kMt4Vector = 2
};

enum class Particle : o2::aod::femtodatatypes::ParticleType {
  kTrack = 0,
  kTwoTrackResonance = 1,
  kV0 = 2,
  kKink = 3,
  kCascade = 4,
};

enum class McOrigin : o2::aod::femtodatatypes::McOriginType {
  kNoMcParticle = 0,       // no associated mc particle normally indicated a wrongly reconstruced partilce
  kFromWrongCollision = 1, // partilce originates from the wrong collision or a collision which was wrongly reconstructed (like a split vertex)
  kPhysicalPrimary = 2,    // primary particle
  kFromSecondaryDecay = 3, // particle from secondary decay
  kFromMaterial = 4,       // partilce orginates from material
  kMissidentified = 5,     // partilce was kMissidentified (also know as fake)
  kMcOriginLast = 6
  // kFromFakeRecoCollision,
  // kFromUnkown
};

constexpr const char* mcOriginToString(McOrigin origin)
{
  switch (origin) {
    case McOrigin::kNoMcParticle:
      return "NoMcParticle";
    case McOrigin::kFromWrongCollision:
      return "FromWrongCollision";
    case McOrigin::kPhysicalPrimary:
      return "PhysicalPrimary";
    case McOrigin::kFromSecondaryDecay:
      return "FromSecondaryDecay";
    case McOrigin::kFromMaterial:
      return "FromMaterial";
    case McOrigin::kMissidentified:
      return "Missidentified";
    default:
      return "UnknownMcOrigin";
  }
}

constexpr bool
  hasMass(Particle p)
{
  switch (p) {
    case Particle::kV0:
    case Particle::kTwoTrackResonance:
    case Particle::kKink:
    case Particle::kCascade:
      return true;
    default:
      return false;
  }
}

enum class Track : o2::aod::femtodatatypes::TrackType {
  kTrack,
  kV0Daughter,
  kCascadeBachelor,
  kResonanceDaughter,
  kKinkDaughter
};

enum class V0 : o2::aod::femtodatatypes::V0Type {
  kLambda,
  kAntiLambda,
  kK0short
};

enum class Kink : o2::aod::femtodatatypes::KinkType {
  kSigma,
  kSigmaPlus
};

enum class Cascade : o2::aod::femtodatatypes::CascadeType {
  kXi,
  kOmega
};

// enum of supported resonances
enum class TwoTrackResonance : o2::aod::femtodatatypes::TwoTrackResonanceType {
  kRho0,
  kPhi,
  kKstar0,
  kKstar0Bar
};

}; // namespace modes
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_MODES_H_
