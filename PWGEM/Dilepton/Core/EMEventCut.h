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

//
// Class for em event selection
//

#ifndef PWGEM_DILEPTON_CORE_EMEVENTCUT_H_
#define PWGEM_DILEPTON_CORE_EMEVENTCUT_H_

#include "TNamed.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

using namespace std;

class EMEventCut : public TNamed
{
 public:
  EMEventCut() = default;
  EMEventCut(const char* name, const char* title) : TNamed(name, title) {}
  ~EMEventCut() {}

  enum class EMEventCuts : int {
    kSel8 = 0,
    kFT0AND,
    kZvtx,
    kNoTFB,     // no time frame border
    kNoITSROFB, // no ITS read out frame border
    kNoSameBunchPileup,
    kIsVertexITSTPC,
    kIsGoodZvtxFT0vsPV,
    kOccupancy,
    kNoCollInTimeRangeStandard,
    kNoCollInTimeRangeNarrow,
    kNCuts
  };

  template <typename T>
  bool IsSelected(T const& collision) const
  {
    if (mRequireSel8 && !IsSelected(collision, EMEventCuts::kSel8)) {
      return false;
    }
    if (mRequireFT0AND && !IsSelected(collision, EMEventCuts::kFT0AND)) {
      return false;
    }
    if (!IsSelected(collision, EMEventCuts::kZvtx)) {
      return false;
    }
    if (mRequireNoTFB && !IsSelected(collision, EMEventCuts::kNoTFB)) {
      return false;
    }
    if (mRequireNoITSROFB && !IsSelected(collision, EMEventCuts::kNoITSROFB)) {
      return false;
    }
    if (mRequireNoSameBunchPileup && !IsSelected(collision, EMEventCuts::kNoSameBunchPileup)) {
      return false;
    }
    if (mRequireVertexITSTPC && !IsSelected(collision, EMEventCuts::kIsVertexITSTPC)) {
      return false;
    }
    if (mRequireGoodZvtxFT0vsPV && !IsSelected(collision, EMEventCuts::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (!IsSelected(collision, EMEventCuts::kOccupancy)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeStandard && !IsSelected(collision, EMEventCuts::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeNarrow && !IsSelected(collision, EMEventCuts::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool IsSelected(T const& collision, const EMEventCuts& cut) const
  {
    switch (cut) {
      case EMEventCuts::kSel8:
        return collision.sel8();

      case EMEventCuts::kFT0AND:
        return collision.selection_bit(o2::aod::evsel::kIsTriggerTVX);

      case EMEventCuts::kZvtx:
        return mMinZvtx < collision.posZ() && collision.posZ() < mMaxZvtx;

      case EMEventCuts::kNoTFB:
        return collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);

      case EMEventCuts::kNoITSROFB:
        return collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);

      case EMEventCuts::kNoSameBunchPileup:
        return collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);

      case EMEventCuts::kIsVertexITSTPC:
        return collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC);

      case EMEventCuts::kIsGoodZvtxFT0vsPV:
        return collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);

      case EMEventCuts::kOccupancy:
        return mMinOccupancy <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < mMaxOccupancy;

      case EMEventCuts::kNoCollInTimeRangeStandard:
        return collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard);

      case EMEventCuts::kNoCollInTimeRangeNarrow:
        return collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow);

      default:
        return true;
    }
  }

  // Setters
  void SetRequireSel8(bool flag);
  void SetRequireFT0AND(bool flag);
  void SetZvtxRange(float min, float max);
  void SetOccupancyRange(int min, int max);
  void SetRequireNoTFB(bool flag);
  void SetRequireNoITSROFB(bool flag);
  void SetRequireNoSameBunchPileup(bool flag);
  void SetRequireVertexITSTPC(bool flag);
  void SetRequireGoodZvtxFT0vsPV(bool flag);
  void SetRequireNoCollInTimeRangeStandard(bool flag);
  void SetRequireNoCollInTimeRangeNarrow(bool flag);

 private:
  bool mRequireSel8{true};
  bool mRequireFT0AND{true};
  float mMinZvtx{-10.f}, mMaxZvtx{+10.f};
  int mMinOccupancy{static_cast<int>(-1e+9)}, mMaxOccupancy{static_cast<int>(+1e+9)};
  bool mRequireNoTFB{true};
  bool mRequireNoITSROFB{true};
  bool mRequireNoSameBunchPileup{false};
  bool mRequireVertexITSTPC{false};
  bool mRequireGoodZvtxFT0vsPV{false};
  bool mRequireNoCollInTimeRangeStandard{false};
  bool mRequireNoCollInTimeRangeNarrow{false};

  ClassDef(EMEventCut, 1);
};

#endif // PWGEM_DILEPTON_CORE_EMEVENTCUT_H_
