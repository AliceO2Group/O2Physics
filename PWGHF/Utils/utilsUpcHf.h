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

#include <cstdint>

namespace o2::analysis::hf_upc
{

/// \brief Gap type classification for UPC events
enum class GapType : uint8_t {
  GapA = 0,     ///< Gap on A-side (C-side active)
  GapC = 1,     ///< Gap on C-side (A-side active)
  DoubleGap = 2 ///< Double gap (both sides empty)
};

/// \brief Default thresholds for gap determination
namespace defaults
{
constexpr float AmplitudeThresholdFT0A = 100.0f; ///< Amplitude threshold for FT0-A (a.u.)
constexpr float AmplitudeThresholdFT0C = 50.0f;  ///< Amplitude threshold for FT0-C (a.u.)
constexpr float EnergyThresholdZDC = 1.0f;       ///< Energy threshold for ZDC (a.u.)
} // namespace defaults

/// \brief Determine gap type based on FIT and ZDC signals
/// \param ft0A FT0-A amplitude
/// \param ft0C FT0-C amplitude
/// \param zdcA ZDC-A (ZNA) common energy
/// \param zdcC ZDC-C (ZNC) common energy
/// \param amplitudeThresholdFT0A Threshold for FT0-A (default: 100.0)
/// \param amplitudeThresholdFT0C Threshold for FT0-C (default: 50.0)
/// \param energyThresholdZDC Threshold for ZDC (default: 1.0)
/// \return Gap type classification
inline GapType determineGapType(float ft0A, float ft0C, float zdcA, float zdcC,
                                float amplitudeThresholdFT0A = defaults::AmplitudeThresholdFT0A,
                                float amplitudeThresholdFT0C = defaults::AmplitudeThresholdFT0C,
                                float energyThresholdZDC = defaults::EnergyThresholdZDC)
{
  // Gap on A-side: FT0-A empty, FT0-C active, ZNA empty, ZNC active
  if (ft0A < amplitudeThresholdFT0A && ft0C > amplitudeThresholdFT0C &&
      zdcA < energyThresholdZDC && zdcC > energyThresholdZDC) {
    return GapType::GapA;
  }

  // Gap on C-side: FT0-A active, FT0-C empty, ZNA active, ZNC empty
  if (ft0A > amplitudeThresholdFT0A && ft0C < amplitudeThresholdFT0C &&
      zdcA > energyThresholdZDC && zdcC < energyThresholdZDC) {
    return GapType::GapC;
  }

  // Default: Double gap (or no clear gap)
  return GapType::DoubleGap;
}

/// \brief Check if the gap type is a single-sided gap (GapA or GapC)
/// \param gap Gap type
/// \return true if single-sided gap, false otherwise
inline bool isSingleSidedGap(GapType gap)
{
  return (gap == GapType::GapA || gap == GapType::GapC);
}

/// \brief Get gap type name as string
/// \param gap Gap type
/// \return String representation of gap type
inline const char* getGapTypeName(GapType gap)
{
  switch (gap) {
    case GapType::GapA:
      return "GapA";
    case GapType::GapC:
      return "GapC";
    case GapType::DoubleGap:
      return "DoubleGap";
    default:
      return "Unknown";
  }
}

/// \brief Convert gap type to integer for histogram filling
/// \param gap Gap type
/// \return Integer representation (0=GapA, 1=GapC, 2=DoubleGap)
inline int gapTypeToInt(GapType gap)
{
  return static_cast<int>(gap);
}

} // namespace o2::analysis::hf_upc

#endif // PWGHF_UTILS_UTILSUPCHF_H_
