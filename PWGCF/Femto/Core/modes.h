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
  kNoCentCal = BIT(5),
  kPP_Run3 = kPP | kRun3,
  kPP_Run2 = kPP | kRun2,
  kPP_NoCentCal_Run3 = kPP | kRun3 | kNoCentCal,
  kPbPb_Run3 = kPbPb | kRun3,
  kPbPb_Run2 = kPbPb | kRun2,
};

enum class Track : o2::aod::femtodatatypes::TrackType {
  kPrimaryTrack,
  kV0Daughter,
  kCascadeBachelor,
  kResonanceDaughter
};

enum class V0 : o2::aod::femtodatatypes::V0Type {
  kLambda,
  kAntiLambda,
  kK0short
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

enum class Pairs : o2::aod::femtodatatypes::PairType {
  kTrackTrack,
  kTrackV0,
  kTrackResonance,
  kTrackCascade
};

enum class TrackPairs : o2::aod::femtodatatypes::PairType {
  kTrackTrack,
  kTrackPosDaughter,
  kTrackNegDaughter,
  kTrackBachelor
};

}; // namespace modes
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_MODES_H_
