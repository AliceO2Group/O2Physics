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

#include "PWGCF/Femto/Core/pairBuilder.h"

namespace o2::analysis::femto
{
namespace tripletcleaner
{

class TrackTrackTrackTripletCleaner : public paircleaner::BasePairCleaner
{
 public:
  TrackTrackTrackTripletCleaner() = default;
  ~TrackTrackTrackTripletCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& track3, T4 const& /*trackTable*/) const
  {
    return this->isCleanTrackPair(track1, track2) &&
           this->isCleanTrackPair(track2, track3) &&
           this->isCleanTrackPair(track1, track3);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& track3, T4 const& trackTable, T5 const& partonicMothers) const
  {
    if (!this->isCleanTriplet(track1, track2, track3, trackTable)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, track2, partonicMothers) &&
             this->pairHasCommonAncestor(track2, track3, partonicMothers) &&
             this->pairHasCommonAncestor(track1, track3, partonicMothers);
    }
    if (mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, track2, partonicMothers) &&
             this->pairHasNonCommonAncestor(track2, track3, partonicMothers) &&
             this->pairHasNonCommonAncestor(track1, track3, partonicMothers);
    }
    return true;
  }
};

class TrackTrackV0TripletCleaner : public paircleaner::BasePairCleaner
{
 public:
  TrackTrackV0TripletCleaner() = default;
  ~TrackTrackV0TripletCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& v0, T4 const& trackTable) const
  {
    auto posDaughter = trackTable.rawIteratorAt(v0.posDauId() - trackTable.offset());
    auto negDaughter = trackTable.rawIteratorAt(v0.negDauId() - trackTable.offset());
    return this->isCleanTrackPair(track1, track2) &&
           this->isCleanTrackPair(track1, posDaughter) &&
           this->isCleanTrackPair(track1, negDaughter) &&
           this->isCleanTrackPair(track2, posDaughter) &&
           this->isCleanTrackPair(track2, negDaughter);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool isCleanTriplet(T1 const& track1, T2 const& track2, T3 const& v0, T4 const& trackTable, T5 const& partonicMothers) const
  {
    if (!this->isCleanTriplet(track1, track2, v0, trackTable)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, track2, partonicMothers) &&
             this->pairHasCommonAncestor(track1, v0, partonicMothers) &&
             this->pairHasCommonAncestor(track2, v0, partonicMothers);
    }
    if (mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, track2, partonicMothers) &&
             this->pairHasNonCommonAncestor(track1, v0, partonicMothers) &&
             this->pairHasNonCommonAncestor(track2, v0, partonicMothers);
    }
    return true;
  }
};

} // namespace tripletcleaner
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_TRIPLETCLEANER_H_
