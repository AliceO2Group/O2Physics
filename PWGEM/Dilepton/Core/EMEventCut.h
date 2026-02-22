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

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

#include "TNamed.h"

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
    kIsVertexTOFmatched,
    kIsGoodZvtxFT0vsPV,
    kNoCollInTimeRangeStandard,
    kNoCollInTimeRangeStrict,
    kNoCollInITSROFStandard,
    kNoCollInITSROFStrict,
    kNoHighMultCollInPrevRof,
    kIsGoodITSLayer3,
    kIsGoodITSLayer0123,
    kIsGoodITSLayersAll,
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
    if (mRequireVertexTOFmatched && !IsSelected(collision, EMEventCuts::kIsVertexTOFmatched)) {
      return false;
    }
    if (mRequireGoodZvtxFT0vsPV && !IsSelected(collision, EMEventCuts::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeStandard && !IsSelected(collision, EMEventCuts::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeStrict && !IsSelected(collision, EMEventCuts::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (mRequireNoCollInITSROFStandard && !IsSelected(collision, EMEventCuts::kNoCollInITSROFStandard)) {
      return false;
    }
    if (mRequireNoCollInITSROFStrict && !IsSelected(collision, EMEventCuts::kNoCollInITSROFStrict)) {
      return false;
    }
    if (mRequireNoHighMultCollInPrevRof && !IsSelected(collision, EMEventCuts::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if (mRequireGoodITSLayer3 && !IsSelected(collision, EMEventCuts::kIsGoodITSLayer3)) {
      return false;
    }
    if (mRequireGoodITSLayer0123 && !IsSelected(collision, EMEventCuts::kIsGoodITSLayer0123)) {
      return false;
    }
    if (mRequireGoodITSLayersAll && !IsSelected(collision, EMEventCuts::kIsGoodITSLayersAll)) {
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
        return collision.selection_bit(o2::aod::emevsel::kIsTriggerTVX);

      case EMEventCuts::kZvtx:
        return mMinZvtx < collision.posZ() && collision.posZ() < mMaxZvtx;

      case EMEventCuts::kNoTFB:
        return collision.selection_bit(o2::aod::emevsel::kNoTimeFrameBorder);

      case EMEventCuts::kNoITSROFB:
        return collision.selection_bit(o2::aod::emevsel::kNoITSROFrameBorder);

      case EMEventCuts::kNoSameBunchPileup:
        return collision.selection_bit(o2::aod::emevsel::kNoSameBunchPileup);

      case EMEventCuts::kIsVertexITSTPC:
        return collision.selection_bit(o2::aod::emevsel::kIsVertexITSTPC);

      case EMEventCuts::kIsVertexTOFmatched:
        return collision.selection_bit(o2::aod::emevsel::kIsVertexTOFmatched);

      case EMEventCuts::kIsGoodZvtxFT0vsPV:
        return collision.selection_bit(o2::aod::emevsel::kIsGoodZvtxFT0vsPV);

      case EMEventCuts::kNoCollInTimeRangeStandard:
        return collision.selection_bit(o2::aod::emevsel::kNoCollInTimeRangeStandard);

      case EMEventCuts::kNoCollInTimeRangeStrict:
        return collision.selection_bit(o2::aod::emevsel::kNoCollInTimeRangeStrict);

      case EMEventCuts::kNoCollInITSROFStandard:
        return collision.selection_bit(o2::aod::emevsel::kNoCollInRofStandard);

      case EMEventCuts::kNoCollInITSROFStrict:
        return collision.selection_bit(o2::aod::emevsel::kNoCollInRofStrict);

      case EMEventCuts::kNoHighMultCollInPrevRof:
        return collision.selection_bit(o2::aod::emevsel::kNoHighMultCollInPrevRof);

      case EMEventCuts::kIsGoodITSLayer3:
        return collision.selection_bit(o2::aod::emevsel::kIsGoodITSLayer3);

      case EMEventCuts::kIsGoodITSLayer0123:
        return collision.selection_bit(o2::aod::emevsel::kIsGoodITSLayer0123);

      case EMEventCuts::kIsGoodITSLayersAll:
        return collision.selection_bit(o2::aod::emevsel::kIsGoodITSLayersAll);

      default:
        return true;
    }
  }

  // Setters
  void SetRequireSel8(bool flag);
  void SetRequireFT0AND(bool flag);
  void SetZvtxRange(float min, float max);
  void SetRequireNoTFB(bool flag);
  void SetRequireNoITSROFB(bool flag);
  void SetRequireNoSameBunchPileup(bool flag);
  void SetRequireVertexITSTPC(bool flag);
  void SetRequireVertexTOFmatched(bool flag);
  void SetRequireGoodZvtxFT0vsPV(bool flag);
  void SetRequireNoCollInTimeRangeStandard(bool flag);
  void SetRequireNoCollInTimeRangeStrict(bool flag);
  void SetRequireNoCollInITSROFStandard(bool flag);
  void SetRequireNoCollInITSROFStrict(bool flag);
  void SetRequireNoHighMultCollInPrevRof(bool flag);
  void SetRequireGoodITSLayer3(bool flag);
  void SetRequireGoodITSLayer0123(bool flag);
  void SetRequireGoodITSLayersAll(bool flag);

 private:
  bool mRequireSel8{false};
  bool mRequireFT0AND{true};
  float mMinZvtx{-10.f}, mMaxZvtx{+10.f};
  bool mRequireNoTFB{false};
  bool mRequireNoITSROFB{false};
  bool mRequireNoSameBunchPileup{false};
  bool mRequireVertexITSTPC{false};
  bool mRequireVertexTOFmatched{false};
  bool mRequireGoodZvtxFT0vsPV{false};
  bool mRequireNoCollInTimeRangeStandard{false};
  bool mRequireNoCollInTimeRangeStrict{false};
  bool mRequireNoCollInITSROFStandard{false};
  bool mRequireNoCollInITSROFStrict{false};
  bool mRequireNoHighMultCollInPrevRof{false};
  bool mRequireGoodITSLayer3{false};
  bool mRequireGoodITSLayer0123{false};
  bool mRequireGoodITSLayersAll{false};

  ClassDef(EMEventCut, 1);
};

#endif // PWGEM_DILEPTON_CORE_EMEVENTCUT_H_
