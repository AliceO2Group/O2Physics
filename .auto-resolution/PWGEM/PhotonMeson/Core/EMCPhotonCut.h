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
// Class for emcal photon selection
//

#ifndef PWGEM_PHOTONMESON_CORE_EMCPHOTONCUT_H_
#define PWGEM_PHOTONMESON_CORE_EMCPHOTONCUT_H_

#include <set>
#include <vector>
#include <utility>
#include <string>
#include <optional>
#include "Framework/Logger.h"
#include "Framework/DataTypes.h"
#include "Rtypes.h"
#include "TNamed.h"

class EMCPhotonCut : public TNamed
{
 public:
  EMCPhotonCut() = default;
  EMCPhotonCut(const char* name, const char* title) : TNamed(name, title) {}

  enum class EMCPhotonCuts : int {
    // cluster cut
    kDefinition = 0,
    kEnergy,
    kNCell,
    kM02,
    kTiming,
    kTM,
    kExotic,
    kNCuts
  };

  static const char* mCutNames[static_cast<int>(EMCPhotonCuts::kNCuts)];

  // Temporary function to check if cluster passes selection criteria. To be replaced by framework filters.
  template <typename T, typename Cluster>
  bool IsSelected(Cluster const& cluster) const
  {
    if (!IsSelectedEMCal(EMCPhotonCuts::kDefinition, cluster)) {
      return false;
    }
    if (!IsSelectedEMCal(EMCPhotonCuts::kEnergy, cluster)) {
      return false;
    }
    if (!IsSelectedEMCal(EMCPhotonCuts::kNCell, cluster)) {
      return false;
    }
    if (!IsSelectedEMCal(EMCPhotonCuts::kM02, cluster)) {
      return false;
    }
    if (!IsSelectedEMCal(EMCPhotonCuts::kTiming, cluster)) {
      return false;
    }
    if (mUseTM && (!IsSelectedEMCal(EMCPhotonCuts::kTM, cluster))) {
      return false;
    }
    if (!IsSelectedEMCal(EMCPhotonCuts::kExotic, cluster)) {
      return false;
    }
    return true;
  }

  // Returns true if a cluster survives the cuts!
  template <typename Cluster>
  bool IsSelectedEMCal(const EMCPhotonCuts& cut, Cluster const& cluster) const
  {
    switch (cut) {
      case EMCPhotonCuts::kDefinition:
        return cluster.definition() == mDefinition;

      case EMCPhotonCuts::kEnergy:
        return cluster.e() > mMinE;

      case EMCPhotonCuts::kNCell:
        return cluster.nCells() >= mMinNCell;

      case EMCPhotonCuts::kM02:
        return (cluster.nCells() == 1 || (mMinM02 <= cluster.m02() && cluster.m02() <= mMaxM02));

      case EMCPhotonCuts::kTiming:
        return mMinTime <= cluster.time() && cluster.time() <= mMaxTime;

      case EMCPhotonCuts::kTM: {
        auto trackseta = cluster.tracketa(); // std:vector<float>
        auto tracksphi = cluster.trackphi(); // std:vector<float>
        auto trackspt = cluster.trackpt();   // std:vector<float>
        auto tracksp = cluster.trackp();     // std:vector<float>
        int ntrack = tracksp.size();
        for (int itr = 0; itr < ntrack; itr++) {
          float dEta = fabs(trackseta[itr] - cluster.eta());
          float dPhi = fabs(tracksphi[itr] - cluster.phi());
          bool result = (dEta > mTrackMatchingEta(trackspt[itr])) || (dPhi > mTrackMatchingPhi(trackspt[itr])) || (cluster.e() / tracksp[itr] >= mMinEoverP);
          if (!result) {
            return false;
          }
        }
        return true; // when we don't have any tracks the cluster should always survive the TM cut!
      }

      case EMCPhotonCuts::kExotic:
        return mUseExoticCut ? !cluster.isExotic() : true;

      default:
        return false;
    }
  }

  // Setters
  void SetClusterizer(std::string clusterDefinitionString = "kV3Default");
  void SetMinE(float min = 0.7f);
  void SetMinNCell(int min = 1);
  void SetM02Range(float min = 0.1f, float max = 0.7f);
  void SetTimeRange(float min = -20.f, float max = 25.f);
  void SetTrackMatchingEta(std::function<float(float)> funcTM);
  void SetTrackMatchingPhi(std::function<float(float)> funcTM);
  void SetMinEoverP(float min = 0.7f);
  void SetUseExoticCut(bool flag = true);
  void SetUseTM(bool flag = true);

  /// @brief Print the cluster selection
  void print() const;

 private:
  // EMCal cluster cuts
  int mDefinition{10};      ///< clusterizer definition
  float mMinE{0.7f};        ///< minimum energy
  int mMinNCell{1};         ///< minimum number of cells per cluster
  float mMinM02{0.1f};      ///< minimum M02 for a cluster
  float mMaxM02{0.7f};      ///< maximum M02 for a cluster
  float mMinTime{-20.f};    ///< minimum cluster timing
  float mMaxTime{25.f};     ///< maximum cluster timing
  float mMinEoverP{1.75f};  ///< minimum cluster energy over track momentum ratio needed for the pair to be considered matched
  bool mUseExoticCut{true}; ///< flag to decide if the exotic cluster cut is to be checked or not
  bool mUseTM{true};        ///< flag to decide if track matching cut is to be checek or not

  std::function<float(float)> mTrackMatchingEta{}; ///< function to get check if a pre matched track and cluster pair is considered an actual match for eta
  std::function<float(float)> mTrackMatchingPhi{}; ///< function to get check if a pre matched track and cluster pair is considered an actual match for phi

  ClassDef(EMCPhotonCut, 1);
};

#endif // PWGEM_PHOTONMESON_CORE_EMCPHOTONCUT_H_
