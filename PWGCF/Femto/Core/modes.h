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

#include <cstdint>
#include <type_traits>

namespace o2::analysis::femto::modes
{

#define BIT(n) (1ULL << (n))

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
  kReco = BIT(0),
  kQa = BIT(1),
  kMc = BIT(2),
  kSe = BIT(3),
  kMe = BIT(4),
  kReco_Qa = kReco | kQa,
  kReco_Mc = kReco | kMc,
  kReco_Qa_Mc = kReco | kQa | kMc,
  kSe_Reco = kSe | kReco,
  kMe_Reco = kMe | kReco,
  kSe_Reco_Mc = kSe | kReco | kMc,
  kMe_Reco_Mc = kMe | kReco | kMc,
  kSe_Mc = kSe | kMc,
  kMe_Mc = kMe | kMc,
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
  kPbPb_Run3_MC = kPbPb | kRun3 | kMC,
  kPbPb_Run2 = kPbPb | kRun2,
};

enum class MomentumType : o2::analysis::femto::datatypes::MomentumType {
  kPt = 0,    // transverse momentum
  kPAtPv = 1, // momentum at primary vertex
  kPTpc = 2,  // momentum at inner wall of tpc
};

enum class TransverseMassType : o2::analysis::femto::datatypes::TransverseMassType {
  kAveragePdgMass = 0,
  kReducedPdgMass = 1,
  kMt4Vector = 2
};

enum class Particle : o2::analysis::femto::datatypes::ParticleType {
  mcParticle = 0,
  kTrack = 1,
  kTwoTrackResonance = 2,
  kV0 = 3,
  kKink = 4,
  kCascade = 5,
};

enum class McOrigin : o2::analysis::femto::datatypes::McOriginType {
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

enum class Track : o2::analysis::femto::datatypes::TrackType {
  kTrack,
  kV0Daughter,
  kCascadeBachelor,
  kResonanceDaughter,
  kKinkDaughter
};

enum class V0 : o2::analysis::femto::datatypes::V0Type {
  kLambda,
  kAntiLambda,
  kK0short
};

enum class Kink : o2::analysis::femto::datatypes::KinkType {
  kSigma,
  kSigmaPlus
};

enum class Cascade : o2::analysis::femto::datatypes::CascadeType {
  kXi,
  kOmega
};

// enum of supported resonances
enum class TwoTrackResonance : o2::analysis::femto::datatypes::TwoTrackResonanceType {
  kRho0,
  kPhi,
  kKstar0,
  kKstar0Bar
};

}; // namespace o2::analysis::femto::modes
#endif // PWGCF_FEMTO_CORE_MODES_H_
