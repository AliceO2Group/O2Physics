// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file tripletCleaner.h
/// \brief triplet cleaner class
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TRIPLETCLEANER_H_
#define PWGCF_FEMTO_CORE_TRIPLETCLEANER_H_

#include "PWGCF/Femto/Core/pairCleaner.h"
namespace o2::analysis::femto::tripletcleaner
{

constexpr char PrefixTripletCleanerTrackTrackTrackSe[] = "TrackTrackTrackCleaner/SE/";
constexpr char PrefixTripletCleanerTrackTrackTrackMe[] = "TrackTrackTrackCleaner/ME/";
constexpr char PrefixTripletCleanerTrackTrackV0Se[] = "TrackTrackV0Cleaner/SE/";
constexpr char PrefixTripletCleanerTrackTrackV0Me[] = "TrackTrackV0Cleaner/ME/";
constexpr char PrefixTripletCleanerTrackTrackCascadeSe[] = "TrackTrackCascadeCleaner/SE/";
constexpr char PrefixTripletCleanerTrackTrackCascadeMe[] = "TrackTrackCascadeCleaner/ME/";

template <auto& prefix>
class TrackTrackTrackTripletCleaner : public paircleaner::BasePairCleaner<prefix>
{
 public:
  TrackTrackTrackTripletCleaner() = default;
  ~TrackTrackTrackTripletCleaner() override = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& track3, T4 const& /*trackTable*/, T5 const& tripletHistManager) const
  {
    this->fillAll();
    bool isClean = this->isCleanParticlePair(track1, track2) &&
                   this->isCleanParticlePair(track2, track3) &&
                   this->isCleanParticlePair(track1, track3);
    if (!isClean) {
      this->fillBlocked(tripletHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& track3, T4 const& trackTable, T5 const& mcParticles, T6 const& partonicMothers, T7 const& tripletHistManager) const
  {
    if (!this->isCleanTriplet(track1, track2, track3, trackTable, tripletHistManager)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, track2, mcParticles, partonicMothers) &&
             this->pairHasCommonAncestor(track2, track3, mcParticles, partonicMothers) &&
             this->pairHasCommonAncestor(track1, track3, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, track2, mcParticles, partonicMothers) &&
             this->pairHasNonCommonAncestor(track2, track3, mcParticles, partonicMothers) &&
             this->pairHasNonCommonAncestor(track1, track3, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class TrackTrackV0TripletCleaner : public paircleaner::BasePairCleaner<prefix>
{
 public:
  TrackTrackV0TripletCleaner() = default;
  ~TrackTrackV0TripletCleaner() override = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& v0, T4 const& trackTable, T5 const& tripletHistManager) const
  {
    this->fillAll();
    auto posDaughter = trackTable.rawIteratorAt(v0.posDauId() - trackTable.offset());
    auto negDaughter = trackTable.rawIteratorAt(v0.negDauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(track1, track2) &&
                   this->isCleanParticlePair(track1, posDaughter) &&
                   this->isCleanParticlePair(track1, negDaughter) &&
                   this->isCleanParticlePair(track2, posDaughter) &&
                   this->isCleanParticlePair(track2, negDaughter);
    if (!isClean) {
      this->fillBlocked(tripletHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& v0, T4 const& trackTable, T5 const& mcParticles, T6 const& partonicMothers, T7 const& tripletHistManager) const
  {
    if (!this->isCleanTriplet(track1, track2, v0, trackTable, tripletHistManager)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, track2, mcParticles, partonicMothers) &&
             this->pairHasCommonAncestor(track1, v0, mcParticles, partonicMothers) &&
             this->pairHasCommonAncestor(track2, v0, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, track2, mcParticles, partonicMothers) &&
             this->pairHasNonCommonAncestor(track1, v0, mcParticles, partonicMothers) &&
             this->pairHasNonCommonAncestor(track2, v0, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefix>
class TrackTrackCascadeTripletCleaner : public paircleaner::BasePairCleaner<prefix>
{
 public:
  TrackTrackCascadeTripletCleaner() = default;
  ~TrackTrackCascadeTripletCleaner() override = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& cascade, T4 const& trackTable, T5 const& tripletHistManager) const
  {
    this->fillAll();
    auto bachelor = trackTable.rawIteratorAt(cascade.bachelorId() - trackTable.offset());
    auto posDaughter = trackTable.rawIteratorAt(cascade.posDauId() - trackTable.offset());
    auto negDaughter = trackTable.rawIteratorAt(cascade.negDauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(track1, track2) &&
                   this->isCleanParticlePair(track1, posDaughter) &&
                   this->isCleanParticlePair(track1, negDaughter) &&
                   this->isCleanParticlePair(track1, bachelor) &&
                   this->isCleanParticlePair(track2, posDaughter) &&
                   this->isCleanParticlePair(track2, negDaughter) &&
                   this->isCleanParticlePair(track2, bachelor);
    if (!isClean) {
      this->fillBlocked(tripletHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& cascade, T4 const& trackTable, T5 const& mcParticles, T6 const& partonicMothers, T7 const& tripletHistManager) const
  {
    if (!this->isCleanTriplet(track1, track2, cascade, trackTable, tripletHistManager)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, track2, mcParticles, partonicMothers) &&
             this->pairHasCommonAncestor(track1, cascade, mcParticles, partonicMothers) &&
             this->pairHasCommonAncestor(track2, cascade, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, track2, mcParticles, partonicMothers) &&
             this->pairHasNonCommonAncestor(track1, cascade, mcParticles, partonicMothers) &&
             this->pairHasNonCommonAncestor(track2, cascade, mcParticles, partonicMothers);
    }
    return true;
  }
};
} // namespace o2::analysis::femto::tripletcleaner

#endif // PWGCF_FEMTO_CORE_TRIPLETCLEANER_H_
