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

/// \file FemtoDreamTrackBitmpask.h
/// \brief Definition of the FemtoDreamTrackBitmpask
/// \author Luca Barioglio, TUM, luca.barioglio@cern.ch

#ifndef ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMTRACKBITMASK_H_
#define ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMTRACKBITMASK_H_

#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;

/// \class FemtoDreamBitmpask.h
/// \brief Bitmaks to be used to chose cuts
class FemtoDreamTrackBitmask
{
 public:
  FemtoDreamTrackBitmask()
  {
    resetBitmask();
  }

  // Full set of bitmasks

  // sign
  static constexpr uint32_t kSignMinusMask = 1;
  static constexpr uint32_t kSignPlusMask = 1 << 1;

  // pT
  static constexpr uint32_t kMinPtLowMask = 1 << 2;
  static constexpr uint32_t kMinPtMidMask = 1 << 3;
  static constexpr uint32_t kMinPtMinHighMask = 1 << 4;

  // eta
  static constexpr uint32_t kMaxEtaHighMask = 1 << 5;
  static constexpr uint32_t kMaxEtaMidMask = 1 << 6;
  static constexpr uint32_t kMaxEtaLowMask = 1 << 7;

  // TPC cls
  static constexpr uint32_t kMinTPCclsLowMask = 1 << 8;
  static constexpr uint32_t kMinTPCclsMidMask = 1 << 9;
  static constexpr uint32_t kMinTPCclsHighMask = 1 << 10;

  // fraction of findable TPC cls
  static constexpr uint32_t kMinFracTPCclsLowMask = 1 << 11;
  static constexpr uint32_t kMinFracTPCclsMidMask = 1 << 12;
  static constexpr uint32_t kMinFracTPCclshighMask = 1 << 13;

  // shared TPC clusters
  static constexpr uint32_t kMaxSharedTPCclsHighMask = 1 << 14;
  static constexpr uint32_t kMaxSharedTPCclsLowMask = 1 << 15;

  // DCA xy
  static constexpr uint32_t kMaxDCAxyHighMask = 1 << 16;
  static constexpr uint32_t kMaxkMaxDCAxyLowMask = 1 << 17;

  // DCA z
  static constexpr uint32_t kMaxDCAzHighMask = 1 << 18;
  static constexpr uint32_t kMaxkMaxDCAzLowMask = 1 << 19;

  // n sigma PID
  static constexpr uint32_t kMaxNsigmaPidHighMask = 1 << 20;
  static constexpr uint32_t kMaxNsigmaPidMidMask = 1 << 21;
  static constexpr uint32_t kMaxNsigmaPidLowMask = 1 << 22;

  // Fill bitmasks for reset

  static constexpr uint32_t kSignFullMask = kSignMinusMask + kSignPlusMask;
  static constexpr uint32_t kMinPtFullMask = kMinPtLowMask + kMinPtMidMask + kMinPtMinHighMask;
  static constexpr uint32_t kMaxEtaFullMask = kMaxEtaHighMask + kMaxEtaMidMask + kMaxEtaLowMask;
  static constexpr uint32_t kMinTPCclsFullMask = kMinTPCclsLowMask + kMinTPCclsMidMask + kMinTPCclsHighMask;
  static constexpr uint32_t kMaxSharedTPCclsFullMask = kMaxSharedTPCclsHighMask + kMaxSharedTPCclsLowMask;
  static constexpr uint32_t kMaxDCAxyFullMask = kMaxDCAxyHighMask + kMaxkMaxDCAxyLowMask;
  static constexpr uint32_t kMaxDCAzFullMask = kMaxDCAzHighMask + kMaxkMaxDCAzLowMask;
  static constexpr uint32_t kMaxNsigmaPidFullMask = kMaxNsigmaPidHighMask + kMaxNsigmaPidMidMask + kMaxNsigmaPidLowMask;

  template <typename T>
  bool acceptTrack(const T& track) const
  {
    return (mTotalBitmask & track.cut()) == mTotalBitmask;
  }

  uint32_t getBitmask() const
  {
    return mTotalBitmask;
  }

  void resetBitmask()
  {
    mTotalBitmask = kMinPtLowMask | kMaxEtaHighMask | kMinTPCclsLowMask | kMinFracTPCclsLowMask | kMaxSharedTPCclsHighMask | kMaxDCAxyHighMask | kMaxDCAzHighMask;
  }
  /// Returns true if the track passes the selections encoded in the bitmask

 private:
  uint32_t mTotalBitmask;
};

#endif /* ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMTRACKBITMASK_H_ */
