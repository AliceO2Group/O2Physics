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

/// \file EMCPhotonCut.h
/// \brief Header of class for emcal photon selection.
/// \author M. Hemmer, marvin.hemmer@cern.ch; N. Strangmann, nicolas.strangmann@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_EMCPHOTONCUT_H_
#define PWGEM_PHOTONMESON_CORE_EMCPHOTONCUT_H_

#include <TNamed.h>

#include <Rtypes.h>

#include <cmath>
#include <string>

struct TrackMatchingParams {
  float a{0.01f};
  float b{4.07f};
  float c{-2.5f};
};

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
    kSecondaryTM,
    kExotic,
    kNCuts
  };

  static const char* mCutNames[static_cast<int>(EMCPhotonCuts::kNCuts)];

  /// \brief check if given cluster survives all cuts
  /// \param cluster cluster to check
  /// \return true if cluster survives all cuts else false
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
    if (mUseSecondaryTM && (!IsSelectedEMCal(EMCPhotonCuts::kSecondaryTM, cluster))) {
      return false;
    }
    if (!IsSelectedEMCal(EMCPhotonCuts::kExotic, cluster)) {
      return false;
    }
    return true;
  }

  /// \brief check if given cluster survives a given cut
  /// \param cut enum of the cluster cut to check
  /// \param cluster cluster to check
  /// \return true if cluster survives cut else false
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
        auto dEtas = cluster.deltaEta();   // std:vector<float>
        auto dPhis = cluster.deltaPhi();   // std:vector<float>
        auto trackspt = cluster.trackpt(); // std:vector<float>
        auto tracksp = cluster.trackp();   // std:vector<float>
        int ntrack = tracksp.size();
        for (int itr = 0; itr < ntrack; itr++) {
          float dEta = std::fabs(dEtas[itr]);
          float dPhi = std::fabs(dPhis[itr]);
          bool result = (dEta > GetTrackMatchingEta(trackspt[itr])) || (dPhi > GetTrackMatchingPhi(trackspt[itr])) || (cluster.e() / tracksp[itr] >= mMinEoverP);
          if (!result) {
            return false;
          }
        }
        return true; // when we don't have any tracks the cluster should always survive the TM cut!
      }
      case EMCPhotonCuts::kSecondaryTM: {
        auto dEtas = cluster.deltaEtaSec();   // std:vector<float>
        auto dPhis = cluster.deltaPhiSec();   // std:vector<float>
        auto trackspt = cluster.trackptSec(); // std:vector<float>
        auto tracksp = cluster.trackpSec();   // std:vector<float>
        int ntrack = tracksp.size();
        for (int itr = 0; itr < ntrack; itr++) {
          float dEta = std::fabs(dEtas[itr]);
          float dPhi = std::fabs(dPhis[itr]);
          bool result = (dEta > GetSecTrackMatchingEta(trackspt[itr])) || (dPhi > GetSecTrackMatchingPhi(trackspt[itr]));
          if (!result) {
            return false;
          }
        }
        return true; // when we don't have any secondary tracks the cluster should always survive the secondary TM cut!
      }

      case EMCPhotonCuts::kExotic:
        return mUseExoticCut ? !cluster.isExotic() : true;

      default:
        return false;
    }
  }

  // Setters
  /// \brief Set clusterizer
  /// \param clusterDefinitionString name of the clusterizer
  void SetClusterizer(std::string clusterDefinitionString = "kV3Default");

  /// \brief Set minimum cluster energy
  /// \param min minimum cluster energy
  void SetMinE(float min = 0.7f);

  /// \brief Set minimum number of cells per cluster
  /// \param min minimum number of cells per cluster
  void SetMinNCell(int min = 1);

  /// \brief Set cluster M02 range to select
  /// \param min minimum allowed cluster M02
  /// \param max maximum allowed cluster M02
  void SetM02Range(float min = 0.1f, float max = 0.7f);

  /// \brief Set cluster time range to select
  /// \param min minimum allowed cluster time
  /// \param max maximum allowed cluster time
  void SetTimeRange(float min = -20.f, float max = 25.f);

  /// \brief Set minimum cluster E over track momentum for track matching
  /// \param min minimum allowed E/p
  void SetMinEoverP(float min = 0.7f);

  /// \brief Set flag to reject exotic cluster
  /// \param flag flag to reject exotic cluster
  void SetUseExoticCut(bool flag = true);

  /// \brief Set flag to use track matching
  /// \param flag flag to use track matching
  void SetUseTM(bool flag = true);

  /// \brief Set flag to use secondary track matching
  /// \param flag flag to use secondary track matching
  void SetUseSecondaryTM(bool flag = false);

  /// \brief Set parameters for track matching delta eta = a + (pT + b)^c
  /// \param a a in a + (pT + b)^c
  /// \param b b in a + (pT + b)^c
  /// \param c c in a + (pT + b)^c
  void SetTrackMatchingEtaParams(float a, float b, float c)
  {
    mTrackMatchingEtaParams = {a, b, c};
  }

  /// \brief Set parameters for track matching delta phi = a + (pT + b)^c
  /// \param a a in a + (pT + b)^c
  /// \param b b in a + (pT + b)^c
  /// \param c c in a + (pT + b)^c
  void SetTrackMatchingPhiParams(float a, float b, float c)
  {
    mTrackMatchingPhiParams = {a, b, c};
  }

  /// \brief Set parameters for secondary track matching delta eta = a + (pT + b)^c
  /// \param a a in a + (pT + b)^c
  /// \param b b in a + (pT + b)^c
  /// \param c c in a + (pT + b)^c
  void SetSecTrackMatchingEtaParams(float a, float b, float c)
  {
    mSecTrackMatchingEtaParams = {a, b, c};
  }

  /// \brief Set parameters for secondary track matching delta phi = a + (pT + b)^c
  /// \param a a in a + (pT + b)^c
  /// \param b b in a + (pT + b)^c
  /// \param c c in a + (pT + b)^c
  void SetSecTrackMatchingPhiParams(float a, float b, float c)
  {
    mSecTrackMatchingPhiParams = {a, b, c};
  }

  /// \brief calculate delta eta for track matching at given track pT
  /// \param pT track pT
  float GetTrackMatchingEta(float pT) const
  {
    return mTrackMatchingEtaParams.a + std::pow(pT + mTrackMatchingEtaParams.b, mTrackMatchingEtaParams.c);
  }

  /// \brief calculate delta phi for track matching at given track pT
  /// \param pT track pT
  float GetTrackMatchingPhi(float pT) const
  {
    return mTrackMatchingPhiParams.a + std::pow(pT + mTrackMatchingPhiParams.b, mTrackMatchingPhiParams.c);
  }

  /// \brief calculate delta eta for secondary track matching at given track pT
  /// \param pT track pT
  float GetSecTrackMatchingEta(float pT) const
  {
    return mSecTrackMatchingEtaParams.a + std::pow(pT + mSecTrackMatchingEtaParams.b, mSecTrackMatchingEtaParams.c);
  }

  /// \brief calculate delta phi for secondary track matching at given track pT
  /// \param pT track pT
  float GetSecTrackMatchingPhi(float pT) const
  {
    return mSecTrackMatchingPhiParams.a + std::pow(pT + mSecTrackMatchingPhiParams.b, mSecTrackMatchingPhiParams.c);
  }

  /// \brief Print the cluster selection
  void print() const;

 private:
  // EMCal cluster cuts
  int mDefinition{10};         ///< clusterizer definition
  float mMinE{0.7f};           ///< minimum energy
  int mMinNCell{1};            ///< minimum number of cells per cluster
  float mMinM02{0.1f};         ///< minimum M02 for a cluster
  float mMaxM02{0.7f};         ///< maximum M02 for a cluster
  float mMinTime{-20.f};       ///< minimum cluster timing
  float mMaxTime{25.f};        ///< maximum cluster timing
  float mMinEoverP{1.75f};     ///< minimum cluster energy over track momentum ratio needed for the pair to be considered matched
  bool mUseExoticCut{true};    ///< flag to decide if the exotic cluster cut is to be checked or not
  bool mUseTM{true};           ///< flag to decide if track matching cut is to be checek or not
  bool mUseSecondaryTM{false}; ///< flag to decide if seconary track matching cut is to be checek or not

  TrackMatchingParams mTrackMatchingEtaParams = {-1.f, 0.f, 0.f};
  TrackMatchingParams mTrackMatchingPhiParams = {-1.f, 0.f, 0.f};
  TrackMatchingParams mSecTrackMatchingEtaParams = {-1.f, 0.f, 0.f};
  TrackMatchingParams mSecTrackMatchingPhiParams = {-1.f, 0.f, 0.f};

  ClassDef(EMCPhotonCut, 2);
};

#endif // PWGEM_PHOTONMESON_CORE_EMCPHOTONCUT_H_
