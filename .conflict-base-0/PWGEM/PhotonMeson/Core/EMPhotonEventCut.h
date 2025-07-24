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

#include "PWGEM/Dilepton/Core/EMEventCut.h"

using namespace std;

class EMPhotonEventCut : public EMEventCut
{
 public:
  EMPhotonEventCut() = default;
  EMPhotonEventCut(const char* name, const char* title) : EMEventCut(name, title) {}

  enum class EMPhotonEventCuts : int {
    kEMCReadoutInMB = 0,
    kEMCHardwareTriggered,
    kNCuts
  };

  template <typename T>
  bool IsSelected(T const& collision) const
  {
    if (!EMEventCut::IsSelected(collision)) {
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
      case EMPhotonEventCuts::kEMCReadoutInMB:
        return (collision.alias_bit(kTVXinEMC));

      case EMPhotonEventCuts::kEMCHardwareTriggered:
        return (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7));

      default:
        return true;
    }
  }

  // Setters
  void SetRequireEMCReadoutInMB(bool flag);
  void SetRequireEMCHardwareTriggered(bool flag);

 private:
  bool mRequireEMCReadoutInMB{false};
  bool mRequireEMCHardwareTriggered{false};

  ClassDef(EMPhotonEventCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_EMPHOTONEVENTCUT_H_
