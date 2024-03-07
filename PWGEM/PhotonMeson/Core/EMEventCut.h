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

#ifndef PWGEM_PHOTONMESON_CORE_EMEVENTCUT_H_
#define PWGEM_PHOTONMESON_CORE_EMEVENTCUT_H_

#include "TNamed.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

using namespace std;

class EMEventCut : public TNamed
{
 public:
  EMEventCut() = default;
  EMEventCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class EMEventCuts : int {
    kFT0AND = 0, // i.e. sel8
    kZvtx,
    kNoTFB,     // no time frame border
    kNoITSROFB, // no ITS read out frame border
    kNCuts
  };

  static const char* mCutNames[static_cast<int>(EMEventCuts::kNCuts)];

  template <typename T>
  bool IsSelected(T const& collision) const
  {
    if (!IsSelected(collision, EMEventCuts::kFT0AND)) {
      return false;
    }
    if (!IsSelected(collision, EMEventCuts::kZvtx)) {
      return false;
    }
    if (!IsSelected(collision, EMEventCuts::kNoTFB)) {
      return false;
    }
    if (!IsSelected(collision, EMEventCuts::kNoITSROFB)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool IsSelected(T const& collision, const EMEventCuts& cut) const
  {
    switch (cut) {
      case EMEventCuts::kFT0AND:
        return collision.sel8();
        // return collision.selection_bit(o2::aod::evsel::kIsTriggerTVX); // alternative way.

      case EMEventCuts::kZvtx:
        return mMinZvtx < collision.posZ() && collision.posZ() < mMaxZvtx;

      case EMEventCuts::kNoTFB:
        return collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);

      case EMEventCuts::kNoITSROFB:
        return collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);

      default:
        return false;
    }
  }

  // Setters
  void SetRequireFT0AND(bool flag);
  void SetZvtxRange(float min, float max);
  void SetRequireNoTFB(bool flag);
  void SetRequireNoITSROFB(bool flag);

  /// @brief Print the track selection
  void print() const;

 private:
  bool mRequireFT0AND{true};
  float mMinZvtx{-10.f}, mMaxZvtx{+10.f};
  bool mRequireNoTFB{true};
  bool mRequireNoITSROFB{true};

  ClassDef(EMEventCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_EMEVENTCUT_H_
