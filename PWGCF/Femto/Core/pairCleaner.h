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

/// \file pairCleaner.h
/// \brief pair cleaner class
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PAIRCLEANER_H_
#define PWGCF_FEMTO_CORE_PAIRCLEANER_H_

#include "PWGCF/Femto/Core/modes.h"

#include "fairlogger/Logger.h"

namespace o2::analysis::femto
{
namespace paircleaner
{

class BasePairCleaner
{
 public:
  BasePairCleaner() = default;
  virtual ~BasePairCleaner() = default;

  template <modes::Mode mode, typename T>
  void init(T const& PairCuts)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kMc)) {
      mMixPairsWithCommonAncestor = PairCuts.mixOnlyCommonAncestor.value;
      mMixPairsWithNonCommonAncestor = PairCuts.mixOnlyNonCommonAncestor.value;
      if (mMixPairsWithCommonAncestor && mMixPairsWithNonCommonAncestor) {
        LOG(fatal) << "Both mixing with common and non-common ancestor is activated. Breaking...";
      }
    }
  }

 protected:
  template <typename T1, typename T2>
  bool isCleanTrackPair(T1 const& track1, T2 const& track2) const
  {
    return track1.globalIndex() != track2.globalIndex();
  };

  template <typename T1, typename T2, typename T3>
  bool pairHasCommonAncestor(T1 const& particle1, T2 const& particle2, T3 const& /*partonicMothers*/) const
  {
    // if one of the two particles has no associated partonic mother, we cannot know if they have a common anchestor, so we break out with false
    if (!particle1.has_fMcPartMoth() || !particle2.has_fMcPartMoth()) {
      return false;
    }
    // get mc particles
    auto partonicMother1 = particle1.template fMcPartMoth_as<T3>();
    auto partonicMother2 = particle2.template fMcPartMoth_as<T3>();
    // get partonic mothers
    return partonicMother1.globalIndex() == partonicMother2.globalIndex();
  };

  template <typename T1, typename T2, typename T3>
  bool pairHasNonCommonAncestor(T1 const& particle1, T2 const& particle2, T3 const& /*partonicMothers*/) const
  {
    // if one of the two particles has no associated partonic mother, we cannot now if they have a non-common anchestor, so we break out with false
    if (!particle1.has_fMcPartMoth() || !particle2.has_fMcPartMoth()) {
      return false;
    }
    // get mc particles
    auto partonicMother1 = particle1.template fMcPartMoth_as<T3>();
    auto partonicMother2 = particle2.template fMcPartMoth_as<T3>();
    // get partonic mothers
    return partonicMother1.globalIndex() != partonicMother2.globalIndex();
  };

  bool mMixPairsWithCommonAncestor = false;
  bool mMixPairsWithNonCommonAncestor = false;
};

class TrackTrackPairCleaner : public BasePairCleaner
{
 public:
  TrackTrackPairCleaner() = default;

  template <typename T1, typename T2, typename T3>
  bool isCleanPair(T1 const& track1, T2 const& track2, T3 const& /*trackTable*/) const
  {
    return this->isCleanTrackPair(track1, track2);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track1, T2 const& track2, T3 const& trackTable, T4 const& partonicMothers) const
  {
    if (!this->isCleanPair(track1, track2, trackTable)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, track2, partonicMothers);
    }
    if (mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, track2, partonicMothers);
    }
    return true;
  }
};

class V0V0PairCleaner : public BasePairCleaner
{
 public:
  V0V0PairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& v01, const T2& v02, const T3& trackTable) const
  {
    auto posDaughter1 = trackTable.rawIteratorAt(v01.posDauId() - trackTable.offset());
    auto negDaughter1 = trackTable.rawIteratorAt(v01.negDauId() - trackTable.offset());
    auto posDaughter2 = trackTable.rawIteratorAt(v02.posDauId() - trackTable.offset());
    auto negDaughter2 = trackTable.rawIteratorAt(v02.negDauId() - trackTable.offset());
    return this->isCleanTrackPair(posDaughter1, posDaughter2) && this->isCleanTrackPair(negDaughter1, negDaughter2);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& v01, T2 const& v02, T3 const& trackTable, T4 const& partonicMothers) const
  {
    if (!this->isCleanPair(v01, v02, trackTable)) {
      return false;
    }
    // pair is clean
    // no check if we require common or non-common ancestry
    if (mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(v01, v02, partonicMothers);
    }
    if (mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(v01, v02, partonicMothers);
    }
    return true;
  }
};

class TrackV0PairCleaner : public BasePairCleaner // also works for particles decaying into a positive and negative daughter, like resonances
{
 public:
  TrackV0PairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& track, const T2& v0, const T3& trackTable) const
  {
    auto posDaughter = trackTable.rawIteratorAt(v0.posDauId() - trackTable.offset());
    auto negDaughter = trackTable.rawIteratorAt(v0.negDauId() - trackTable.offset());
    return (this->isCleanTrackPair(posDaughter, track) && this->isCleanTrackPair(negDaughter, track));
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track1, T2 const& v0, T3 const& trackTable, T4 const& partonicMothers) const
  {
    if (!this->isCleanPair(track1, v0, trackTable)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, v0, partonicMothers);
    }
    if (mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, v0, partonicMothers);
    }
    return true;
  }
};

class TrackKinkPairCleaner : public BasePairCleaner
{
 public:
  TrackKinkPairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& track, const T2& kink, const T3& trackTable) const
  {
    auto chaDaughter = trackTable.rawIteratorAt(kink.chaDauId() - trackTable.offset());
    return this->isCleanTrackPair(chaDaughter, track);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track1, T2 const& kink, T3 const& trackTable, T4 const& partonicMothers) const
  {
    if (!this->isCleanPair(track1, kink, trackTable)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, kink, partonicMothers);
    }
    if (mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, kink, partonicMothers);
    }
    return true;
  }
};

class TrackCascadePairCleaner : public BasePairCleaner
{
 public:
  TrackCascadePairCleaner() = default;
  template <typename T1, typename T2, typename T3>
  bool isCleanPair(const T1& track, const T2& cascade, const T3& trackTable) const
  {
    auto bachelor = trackTable.rawIteratorAt(cascade.bachelorId() - trackTable.offset());
    auto posDaughter = trackTable.rawIteratorAt(cascade.posDauId() - trackTable.offset());
    auto negDaughter = trackTable.rawIteratorAt(cascade.negDauId() - trackTable.offset());
    return (this->isCleanTrackPair(bachelor, track) && this->isCleanTrackPair(posDaughter, track) && this->isCleanTrackPair(negDaughter, track));
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& track1, T2 const& cascade, T3 const& trackTable, T4 const& partonicMothers) const
  {
    if (!this->isCleanPair(track1, cascade, trackTable)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(track1, cascade, partonicMothers);
    }
    if (mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(track1, cascade, partonicMothers);
    }
    return true;
  }
};
} // namespace paircleaner
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_PAIRCLEANER_H_
