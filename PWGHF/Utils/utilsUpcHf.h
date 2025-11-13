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

/// \file utilsUpcHf.h
/// \brief Utility functions for Ultra-Peripheral Collision (UPC) analysis in Heavy Flavor physics
///
/// \author Minjung Kim <minjung.kim@cern.ch>, CERN

#ifndef PWGHF_UTILS_UTILSUPCHF_H_
#define PWGHF_UTILS_UTILSUPCHF_H_

#include "PWGUD/Core/SGCutParHolder.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UDHelpers.h"

namespace o2::analysis::hf_upc
{

/// \brief Use TrueGap enum from SGSelector for gap type classification
using o2::aod::sgselector::TrueGap;

/// \brief Default thresholds for gap determination
namespace defaults
{
constexpr float AmplitudeThresholdFV0A = 100.0f; ///< Amplitude threshold for FV0-A (a.u.)
constexpr float AmplitudeThresholdFT0A = 100.0f; ///< Amplitude threshold for FT0-A (a.u.)
constexpr float AmplitudeThresholdFT0C = 50.0f;  ///< Amplitude threshold for FT0-C (a.u.)
constexpr float MaxFITTime = 4.0f;               ///< Maximum FIT time (ns)
constexpr int NDtColl = 1000;                    ///< Time window for BC range (ns)
constexpr int MinNBCs = 7;                       ///< Minimum number of BCs to check
constexpr int MinNTracks = 0;                    ///< Minimum number of tracks
constexpr int MaxNTracks = 100;                  ///< Maximum number of tracks
} // namespace defaults

/// \brief Determine gap type using SGSelector with BC range checking
/// \tparam TCollision Collision type
/// \tparam TBCs BC table type
/// \param collision Collision object
/// \param bcs BC table
/// \param amplitudeThresholdFV0A Threshold for FV0-A (default: 100.0)
/// \param amplitudeThresholdFT0A Threshold for FT0-A (default: 100.0)
/// \param amplitudeThresholdFT0C Threshold for FT0-C (default: 50.0)
/// \return SelectionResult with gap type value and BC pointer
template <typename TCollision, typename TBCs>
inline auto determineGapType(TCollision const& collision,
                             TBCs const& bcs,
                             float amplitudeThresholdFV0A = defaults::AmplitudeThresholdFV0A,
                             float amplitudeThresholdFT0A = defaults::AmplitudeThresholdFT0A,
                             float amplitudeThresholdFT0C = defaults::AmplitudeThresholdFT0C)
{
  using BCType = std::decay_t<decltype(collision.template foundBC_as<TBCs>())>;

  // Configure SGSelector thresholds
  SGCutParHolder sgCuts;
  sgCuts.SetNDtcoll(defaults::NDtColl);
  sgCuts.SetMinNBCs(defaults::MinNBCs);
  sgCuts.SetNTracks(defaults::MinNTracks, defaults::MaxNTracks);
  sgCuts.SetMaxFITtime(defaults::MaxFITTime);
  sgCuts.SetFITAmpLimits({amplitudeThresholdFV0A, amplitudeThresholdFT0A, amplitudeThresholdFT0C});

  // Get BC and BC range
  if (!collision.has_foundBC()) {
    return SelectionResult<BCType>{TrueGap::NoGap, nullptr};
  }

  const auto bc = collision.template foundBC_as<TBCs>();
  const auto bcRange = udhelpers::compatibleBCs(collision, sgCuts.NDtcoll(), bcs, sgCuts.minNBCs());

  // Create SGSelector instance and determine gap type with BC range checking
  SGSelector sgSelector;
  const auto sgResult = sgSelector.IsSelected(sgCuts, collision, bcRange, bc);

  return sgResult;
}

/// \brief Check if the gap type is a single-sided gap (SingleGapA or SingleGapC)
/// \param gap TrueGap enum value
/// \return true if single-sided gap, false otherwise
constexpr bool isSingleSidedGap(int gap) noexcept
{
  return (gap == TrueGap::SingleGapA || gap == TrueGap::SingleGapC);
}

/// \brief Get gap type name as string
/// \param gap TrueGap enum value
/// \return String representation of gap type
constexpr const char* getGapTypeName(int gap) noexcept
{
  switch (gap) {
    case TrueGap::NoGap:
      return "NoGap";
    case TrueGap::SingleGapA:
      return "SingleGapA";
    case TrueGap::SingleGapC:
      return "SingleGapC";
    case TrueGap::DoubleGap:
      return "DoubleGap";
    case TrueGap::NoUpc:
      return "NoUpc";
    case TrueGap::TrkOutOfRange:
      return "TrkOutOfRange";
    case TrueGap::BadDoubleGap:
      return "BadDoubleGap";
    default:
      return "Unknown";
  }
}

/// \brief Convert gap type to integer for histogram filling
/// \param gap TrueGap enum value
/// \return Integer representation (-1, 0, 1, 2, 3, 4, 5)
/// \note This is a pass-through function for consistency with other utility functions
constexpr int gapTypeToInt(int gap) noexcept
{
  return gap;
}

} // namespace o2::analysis::hf_upc

#endif // PWGHF_UTILS_UTILSUPCHF_H_
