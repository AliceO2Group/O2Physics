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

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include <Framework/ASoA.h>

#include <TNamed.h>

#include <Rtypes.h>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <string>
#include <vector>

template <typename T>
concept IsTrackIterator = o2::soa::is_iterator<T> && requires(T t) {
  // Check that the *elements* of the container have the required methods:
  { t.deltaEta() } -> std::same_as<float>;
  { t.deltaPhi() } -> std::same_as<float>;
  { t.trackPt() } -> std::same_as<float>;
  { t.trackP() } -> std::same_as<float>;
};

template <typename T>
concept IsTrackContainer = o2::soa::is_table<T> && requires(T t) {
  // Check that the *elements* of the container have the required methods:
  { t.begin().deltaEta() } -> std::same_as<float>;
  { t.begin().deltaPhi() } -> std::same_as<float>;
  { t.begin().trackPt() } -> std::same_as<float>;
  { t.begin().trackP() } -> std::same_as<float>;
};

template <typename Cluster>
concept HasTrackMatching = requires(Cluster cluster) {
  // requires that the following are valid calls:
  { cluster.deltaEta() } -> std::convertible_to<std::vector<float>>;
  { cluster.deltaPhi() } -> std::convertible_to<std::vector<float>>;
  { cluster.trackpt() } -> std::convertible_to<std::vector<float>>;
  { cluster.trackp() } -> std::convertible_to<std::vector<float>>;
};

template <typename Cluster>
concept HasSecondaryMatching = requires(Cluster cluster) {
  // requires that the following are valid calls:
  { cluster.deltaEtaSec() } -> std::convertible_to<std::vector<float>>;
  { cluster.deltaPhiSec() } -> std::convertible_to<std::vector<float>>;
  { cluster.trackptSec() } -> std::convertible_to<std::vector<float>>;
  { cluster.trackpSec() } -> std::convertible_to<std::vector<float>>;
};

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

  constexpr auto getClusterId(o2::soa::is_iterator auto const& t) const
  {
    if constexpr (requires { t.emEmcClusterId(); }) {
      return t.emEmcClusterId();
    } else if constexpr (requires { t.minClusterId(); }) {
      return t.minClusterId();
    } else {
      return -1;
    }
  }

  /// \brief performs check if track is matched with given cluster
  /// \param cluster cluster to be checked
  /// \param emcmatchedtrack matched track iterator
  /// \param emcmatchedtrackEnd matched track end iterator
  /// \param GetEtaCut lambda to get the eta cut value
  /// \param GetPhiCut lambda to get the phi cut value
  /// \param applyEoverP bool to check if E/p should be checked (for secondaries we do not check this!)
  bool checkTrackMatching(o2::soa::is_iterator auto const& cluster, IsTrackIterator auto& emcmatchedtrack, o2::soa::RowViewSentinel const emcmatchedtrackEnd,
                          bool applyEoverP, auto GetEtaCut, auto GetPhiCut) const
  {
    // advance to cluster
    while (emcmatchedtrack != emcmatchedtrackEnd && getClusterId(emcmatchedtrack) < cluster.globalIndex()) {
      ++emcmatchedtrack;
    }
    // all matched tracks have been checked
    if (emcmatchedtrack == emcmatchedtrackEnd) {
      return true;
    }
    // if all remaining tracks are beyond this cluster, it survives
    if (getClusterId(emcmatchedtrack) > cluster.globalIndex()) {
      return true;
    }
    // iterate over tracks belonging to this cluster
    while (emcmatchedtrack != emcmatchedtrackEnd && getClusterId(emcmatchedtrack) == cluster.globalIndex()) {
      auto dEta = std::fabs(emcmatchedtrack.deltaEta());
      auto dPhi = std::fabs(emcmatchedtrack.deltaPhi());
      auto trackpt = emcmatchedtrack.trackPt();
      auto trackp = emcmatchedtrack.trackP();
      bool fail = (dEta > GetEtaCut(trackpt)) ||
                  (dPhi > GetPhiCut(trackpt)) ||
                  (applyEoverP && cluster.e() / trackp >= mMinEoverP);
      if (!fail) {
        return false; // cluster got a track matche to it
      }
      ++emcmatchedtrack;
    }
    return true; // all tracks checked, cluster survives
  }

  /// \brief check if given clusters survives all cuts
  /// \param flags EMBitFlags where results will be stored
  /// \param cluster cluster table to check
  /// \param matchedTracks matched primary tracks table
  /// \param matchedSecondaries matched secondary tracks table
  void AreSelectedRunning(EMBitFlags& flags, o2::soa::is_table auto const& clusters, IsTrackContainer auto const& emcmatchedtracks, IsTrackContainer auto const& secondaries) const
  {
    auto emcmatchedtrackIter = emcmatchedtracks.begin();
    auto emcmatchedtrackEnd = emcmatchedtracks.end();
    auto secondaryIter = secondaries.begin();
    auto secondaryEnd = secondaries.end();
    size_t iCluster = 0;
    for (const auto& cluster : clusters) {
      if (!IsSelectedRunning(cluster, emcmatchedtrackIter, emcmatchedtrackEnd, secondaryIter, secondaryEnd)) {
        flags.set(iCluster);
      }
      ++iCluster;
    }
  }

  /// \brief check if given cluster survives all cuts
  /// \param cluster cluster to check
  /// \param emcmatchedtrackIter current iterator of matched primary tracks
  /// \param emcmatchedtrackEnd end iterator of matched primary tracks
  /// \param secondaryIter current iterator of matched secondary tracks
  /// \param secondaryEnd end iterator of matched secondary tracks
  /// \return true if cluster survives all cuts else false
  bool IsSelectedRunning(o2::soa::is_iterator auto const& cluster, IsTrackIterator auto& emcmatchedtrackIter, o2::soa::RowViewSentinel const emcmatchedtrackEnd, IsTrackIterator auto& secondaryIter, o2::soa::RowViewSentinel const secondaryEnd) const
  {
    if (!IsSelectedEMCalRunning(EMCPhotonCuts::kDefinition, cluster)) {
      return false;
    }
    if (!IsSelectedEMCalRunning(EMCPhotonCuts::kEnergy, cluster)) {
      return false;
    }
    if (!IsSelectedEMCalRunning(EMCPhotonCuts::kNCell, cluster)) {
      return false;
    }
    if (!IsSelectedEMCalRunning(EMCPhotonCuts::kM02, cluster)) {
      return false;
    }
    if (!IsSelectedEMCalRunning(EMCPhotonCuts::kTiming, cluster)) {
      return false;
    }
    if (mUseTM && (!IsSelectedEMCalRunning(EMCPhotonCuts::kTM, cluster, emcmatchedtrackIter, emcmatchedtrackEnd))) {
      return false;
    }
    if (mUseSecondaryTM && (!IsSelectedEMCalRunning(EMCPhotonCuts::kSecondaryTM, cluster, secondaryIter, secondaryEnd))) {
      return false;
    }
    if (!IsSelectedEMCalRunning(EMCPhotonCuts::kExotic, cluster)) {
      return false;
    }
    return true;
  }

  /// \brief check if given cluster survives a given cut
  /// \param cut enum of the cluster cut to check
  /// \param cluster cluster to check
  /// \return true if cluster survives cut else false
  bool IsSelectedEMCalRunning(const EMCPhotonCuts& cut, o2::soa::is_iterator auto const& cluster) const
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

      case EMCPhotonCuts::kTM:
        return false;

      case EMCPhotonCuts::kSecondaryTM:
        return false;

      case EMCPhotonCuts::kExotic:
        return mUseExoticCut ? !cluster.isExotic() : true;

      default:
        return true;
    }
  }

  /// \brief check if given cluster survives a given cut
  /// \param cut enum of the cluster cut to check
  /// \param cluster cluster to check
  /// \param matchedTrackIter current iterator of matched primary or secondary tracks
  /// \param matchedTrackEnd end iterator of matched primary or secondary tracks
  /// \return true if cluster survives cut else false
  bool IsSelectedEMCalRunning(const EMCPhotonCuts& cut, o2::soa::is_iterator auto const& cluster, IsTrackIterator auto& matchedTrackIter, o2::soa::RowViewSentinel const matchedTrackEnd) const
  {
    switch (cut) {
      case EMCPhotonCuts::kTM:
        return checkTrackMatching(cluster, matchedTrackIter, matchedTrackEnd, true, [this](float pt) { return GetTrackMatchingEta(pt); }, [this](float pt) { return GetTrackMatchingPhi(pt); });

      case EMCPhotonCuts::kSecondaryTM:
        return checkTrackMatching(cluster, matchedTrackIter, matchedTrackEnd, false, [this](float pt) { return GetSecTrackMatchingEta(pt); }, [this](float pt) { return GetSecTrackMatchingPhi(pt); });

      default:
        return true;
    }
  }

  /// \brief check if given cluster survives all cuts
  /// \param cluster cluster to check
  /// \param matchedTracks subtable of the matched primary tracks (optional)
  /// \param matchedSecondaries subtable of the matched secondary tracks (optional)
  /// \return true if cluster survives all cuts else false
  template <o2::soa::is_iterator Cluster, typename TMatchedTracks = std::nullptr_t, typename TMatchedSecondaries = std::nullptr_t>
  bool IsSelected(Cluster const& cluster, TMatchedTracks const& emcmatchedtracks = nullptr, TMatchedSecondaries const& secondaries = nullptr) const
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
    if (mUseTM && (!IsSelectedEMCal(EMCPhotonCuts::kTM, cluster, emcmatchedtracks))) {
      return false;
    }
    if (mUseSecondaryTM && (!IsSelectedEMCal(EMCPhotonCuts::kSecondaryTM, cluster, secondaries))) {
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
  /// \param matchedTracks subtable of the matched primary tracks (optional)
  /// \param matchedSecondaries subtable of the matched secondary tracks (optional)
  /// \return true if cluster survives cut else false
  template <o2::soa::is_iterator Cluster, typename TMatchedTracks = std::nullptr_t>
  bool IsSelectedEMCal(const EMCPhotonCuts& cut, Cluster const& cluster, TMatchedTracks const& emcmatchedtracks = nullptr) const
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
        if constexpr (IsTrackContainer<TMatchedTracks>) {
          for (const auto& emcmatchedtrack : emcmatchedtracks) {
            auto dEta = std::fabs(emcmatchedtrack.deltaEta());
            auto dPhi = std::fabs(emcmatchedtrack.deltaPhi());
            auto trackpt = emcmatchedtrack.trackPt();
            auto trackp = emcmatchedtrack.trackP();
            bool result = (dEta > GetTrackMatchingEta(trackpt)) || (dPhi > GetTrackMatchingPhi(trackpt)) || (cluster.e() / trackp >= mMinEoverP);
            if (!result) {
              return false;
            }
          }
        } else if constexpr (HasTrackMatching<Cluster>) {
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
        } else {
          return true;
        }
        return true; // when we don't have any tracks the cluster should always survive the TM cut!
      }
      case EMCPhotonCuts::kSecondaryTM: {
        if constexpr (IsTrackContainer<TMatchedTracks>) {
          for (const auto& emcmatchedtrack : emcmatchedtracks) {
            auto dEta = std::fabs(emcmatchedtrack.deltaEta());
            auto dPhi = std::fabs(emcmatchedtrack.deltaPhi());
            auto trackpt = emcmatchedtrack.trackPt();
            auto trackp = emcmatchedtrack.trackP();
            bool result = (dEta > GetSecTrackMatchingEta(trackpt)) || (dPhi > GetSecTrackMatchingPhi(trackpt)) || (cluster.e() / trackp >= mMinEoverP);
            if (!result) {
              return false;
            }
          }
        } else if constexpr (HasSecondaryMatching<Cluster>) {
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
        } else {
          return true;
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

  ClassDef(EMCPhotonCut, 3);
};

#endif // PWGEM_PHOTONMESON_CORE_EMCPHOTONCUT_H_
