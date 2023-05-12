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
// Class for v0 photon selection
//

#ifndef PWGEM_PHOTONMESON_CORE_PHOSPHOTONCUT_H_
#define PWGEM_PHOTONMESON_CORE_PHOSPHOTONCUT_H_

#include <set>
#include <vector>
#include <utility>
#include <string>
#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "Rtypes.h"
#include "TNamed.h"

class PHOSPhotonCut : public TNamed
{
 public:
  PHOSPhotonCut() = default;
  PHOSPhotonCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class PHOSPhotonCuts : int {
    kEnergy = 0,
    kDispersion,
    kCPV,
    kNCuts
  };

  static const char* mCutNames[static_cast<int>(PHOSPhotonCuts::kNCuts)];

  // Temporary function to check if track passes selection criteria. To be replaced by framework filters.
  template <typename T, typename Cluster>
  bool IsSelected(Cluster const& cluster) const
  {
    // auto track = cluster.template MatchedTrack_as<T>(); //please implement a column to point matched track index (DECLARE_SOA_ARRAY_INDEX_COLUMN) in SkimPHOSClusters table.
    if (!IsSelectedCluster(cluster, PHOSPhotonCuts::kEnergy)) {
      return false;
    }
    if (!IsSelectedCluster(cluster, PHOSPhotonCuts::kDispersion)) {
      return false;
    }
    if (!IsSelectedCluster(cluster, PHOSPhotonCuts::kCPV)) {
      return false;
    }
    return true;
  }

  // Temporary function to check if track passes a given selection criteria. To be replaced by framework filters.
  template <typename Cluster>
  bool IsSelectedCluster(Cluster const& cls, const PHOSPhotonCuts& cut) const
  {
    switch (cut) {
      case PHOSPhotonCuts::kEnergy:
        return cls.e() >= mMinEnergy && cls.e() <= mMaxEnergy;

      case PHOSPhotonCuts::kDispersion:
        return true;

      case PHOSPhotonCuts::kCPV:
        return true;

      default:
        return false;
    }
  }

  // Setters
  void SetEnergyRange(float min = 0.1f, float max = 1e+10f);

  /// @brief Print the PHOS selection
  void print() const;

 private:
  float mMinEnergy{0.1f}, mMaxEnergy{1e+10f};

  ClassDef(PHOSPhotonCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_PHOSPHOTONCUT_H_
