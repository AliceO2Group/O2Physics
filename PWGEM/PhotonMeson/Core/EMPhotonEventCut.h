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
// Class for em photon event selection
//

#ifndef PWGEM_PHOTONMESON_CORE_EMPHOTONEVENTCUT_H_
#define PWGEM_PHOTONMESON_CORE_EMPHOTONEVENTCUT_H_

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

#include "TNamed.h"

using namespace std;

class EMPhotonEventCut : public TNamed
{
 public:
  EMPhotonEventCut() = default;
  EMPhotonEventCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class EMPhotonEventCuts : int {
    kSel8 = 0,
    kFT0AND,
    kZvtx,
    kNoTFB,
    kNoITSROFB,
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
    kEMCReadoutInMB,
    kEMCHardwareTriggered,
    kNCuts
  };

  template <typename T>
  bool IsSelected(T const& collision) const
  {
    if (mRequireSel8 && !IsSelected(collision, EMPhotonEventCuts::kSel8)) {
      return false;
    }
    if (mRequireFT0AND && !IsSelected(collision, EMPhotonEventCuts::kFT0AND)) {
      return false;
    }
    if (!IsSelected(collision, EMPhotonEventCuts::kZvtx)) {
      return false;
    }
    if (mRequireNoTFB && !IsSelected(collision, EMPhotonEventCuts::kNoTFB)) {
      return false;
    }
    if (mRequireNoITSROFB && !IsSelected(collision, EMPhotonEventCuts::kNoITSROFB)) {
      return false;
    }
    if (mRequireNoSameBunchPileup && !IsSelected(collision, EMPhotonEventCuts::kNoSameBunchPileup)) {
      return false;
    }
    if (mRequireVertexITSTPC && !IsSelected(collision, EMPhotonEventCuts::kIsVertexITSTPC)) {
      return false;
    }
    if (mRequireVertexTOFmatched && !IsSelected(collision, EMPhotonEventCuts::kIsVertexTOFmatched)) {
      return false;
    }
    if (mRequireGoodZvtxFT0vsPV && !IsSelected(collision, EMPhotonEventCuts::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeStandard && !IsSelected(collision, EMPhotonEventCuts::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (mRequireNoCollInTimeRangeStrict && !IsSelected(collision, EMPhotonEventCuts::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (mRequireNoCollInITSROFStandard && !IsSelected(collision, EMPhotonEventCuts::kNoCollInITSROFStandard)) {
      return false;
    }
    if (mRequireNoCollInITSROFStrict && !IsSelected(collision, EMPhotonEventCuts::kNoCollInITSROFStrict)) {
      return false;
    }
    if (mRequireNoHighMultCollInPrevRof && !IsSelected(collision, EMPhotonEventCuts::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if (mRequireGoodITSLayer3 && !IsSelected(collision, EMPhotonEventCuts::kIsGoodITSLayer3)) {
      return false;
    }
    if (mRequireGoodITSLayer0123 && !IsSelected(collision, EMPhotonEventCuts::kIsGoodITSLayer0123)) {
      return false;
    }
    if (mRequireGoodITSLayersAll && !IsSelected(collision, EMPhotonEventCuts::kIsGoodITSLayersAll)) {
      return false;
    }
    if (mRequireEMCReadoutInMB && !IsSelected(collision, EMPhotonEventCuts::kEMCReadoutInMB)) {
      return false;
    }
    if (mRequireEMCHardwareTriggered && !IsSelected(collision, EMPhotonEventCuts::kEMCHardwareTriggered)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool IsSelected(T const& collision, const EMPhotonEventCuts& cut) const
  {
    switch (cut) {
      case EMPhotonEventCuts::kSel8:
        return collision.sel8();

      case EMPhotonEventCuts::kFT0AND:
        return collision.selection_bit(o2::aod::evsel::kIsTriggerTVX);

      case EMPhotonEventCuts::kZvtx:
        return mMinZvtx < collision.posZ() && collision.posZ() < mMaxZvtx;

      case EMPhotonEventCuts::kNoTFB:
        return collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);

      case EMPhotonEventCuts::kNoITSROFB:
        return collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);

      case EMPhotonEventCuts::kNoSameBunchPileup:
        return collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup);

      case EMPhotonEventCuts::kIsVertexITSTPC:
        return collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC);

      case EMPhotonEventCuts::kIsVertexTOFmatched:
        return collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched);

      case EMPhotonEventCuts::kIsGoodZvtxFT0vsPV:
        return collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV);

      case EMPhotonEventCuts::kNoCollInTimeRangeStandard:
        return collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard);

      case EMPhotonEventCuts::kNoCollInTimeRangeStrict:
        return collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict);

      case EMPhotonEventCuts::kNoCollInITSROFStandard:
        return collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard);

      case EMPhotonEventCuts::kNoCollInITSROFStrict:
        return collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict);

      case EMPhotonEventCuts::kNoHighMultCollInPrevRof:
        return collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof);

      case EMPhotonEventCuts::kIsGoodITSLayer3:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer3);

      case EMPhotonEventCuts::kIsGoodITSLayer0123:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123);

      case EMPhotonEventCuts::kIsGoodITSLayersAll:
        return collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll);

      case EMPhotonEventCuts::kEMCReadoutInMB:
        return (collision.alias_bit(kTVXinEMC));

      case EMPhotonEventCuts::kEMCHardwareTriggered:
        return (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7));

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
  void SetRequireEMCReadoutInMB(bool flag);
  void SetRequireEMCHardwareTriggered(bool flag);

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
  bool mRequireEMCReadoutInMB{false};
  bool mRequireEMCHardwareTriggered{false};

  ClassDef(EMPhotonEventCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_EMPHOTONEVENTCUT_H_
