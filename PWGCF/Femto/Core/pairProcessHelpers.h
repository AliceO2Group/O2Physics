// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file pairProcessHelpers.h
/// \brief process functions used in pair tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PAIRPROCESSHELPERS_H_
#define PWGCF_FEMTO_CORE_PAIRPROCESSHELPERS_H_

#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoAHelpers.h"

#include <random>

namespace o2::analysis::femto
{
namespace pairprocesshelpers
{

// process same event for identical tracks
template <typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5>
void processSameEvent(const T1& SliceParticle,
                      T2& ParticleHistManager,
                      T3& PairHistManager,
                      T4& CprManager,
                      T5& rng,
                      bool randomize)
{
  // Fill single particle histograms
  for (auto const& part : SliceParticle) {
    ParticleHistManager.fill(part);
  }

  std::uniform_real_distribution<float> dist(0.f, 1.f);

  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle))) {

    // Close pair rejection
    if (CprManager.isActivated()) {
      CprManager.setPair(p1, p2);
      if (CprManager.isClosePair()) {
        continue;
      }
    }
    CprManager.fill();

    // Randomize pair order if enabled
    float threshold = 0.5f;
    bool swapPair = randomize ? (dist(rng) > threshold) : false;
    if (swapPair) {
      PairHistManager.setPair(p2, p1);
    } else {
      PairHistManager.setPair(p1, p2);
    }
    PairHistManager.fill();
  }
}

// process same event for non-identical tracks
template <typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7>
void processSameEvent(T1& SliceParticle1,
                      T2& SliceParticle2,
                      T3& ParticleHistManager1,
                      T4& ParticleHistManager2,
                      T5& PairHistManager,
                      T6& CprManager,
                      T7& PcManager)
{
  // Fill single particle histograms
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.fill(part);
  }

  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.fill(part);
  }

  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2))) {
    // pair cleaning
    if (!PcManager.isCleanPair(p1, p2)) {
      continue;
    }
    // Close pair rejection
    if (CprManager.isActivated()) {
      CprManager.setPair(p1, p2);
      if (CprManager.isClosePair()) {
        continue;
      }
    }
    CprManager.fill();
    PairHistManager.setPair(p1, p2);
    PairHistManager.fill();
  }
}

// process same event for tracks and particles decaying into tracks
template <typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8>
void processSameEvent(T1& SliceParticle1,
                      T2& SliceParticle2,
                      T3& TrackTable,
                      T4& ParticleHistManager1,
                      T5& ParticleHistManager2,
                      T6& PairHistManager,
                      T7& CprManager,
                      T8& PcManager)
{
  // Fill single particle histograms
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.fill(part);
  }

  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.fill(part, TrackTable);
  }

  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2))) {
    // pair cleaning
    if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
      continue;
    }
    // Close pair rejection
    if (CprManager.isActivated()) {
      CprManager.setPair(p1, p2, TrackTable);
      if (CprManager.isClosePair()) {
        continue;
      }
    }
    CprManager.fill();
    PairHistManager.setPair(p1, p2);
    PairHistManager.fill();
  }
}

// process mixed event identical tracks
template <typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8>
void processMixedEvent(T1& Collisions,
                       T2& Partition,
                       T3& cache,
                       T4& policy,
                       T5& depth,
                       T6& PairHistManager,
                       T7& CprManager,
                       T8& PcManager)
{
  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {
    if (!(std::fabs(collision1.magField() - collision2.magField()) < o2::constants::math::Epsilon)) {
      continue;
    }
    CprManager.setMagField(collision1.magField());
    auto sliceParticle1 = Partition->sliceByCached(o2::aod::femtobase::stored::collisionId, collision1.globalIndex(), cache);
    auto sliceParticle2 = Partition->sliceByCached(o2::aod::femtobase::stored::collisionId, collision2.globalIndex(), cache);
    if (sliceParticle1.size() == 0 || sliceParticle2.size() == 0) {
      continue;
    }
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(sliceParticle1, sliceParticle2))) {
      // pair cleaning
      if (!PcManager.isCleanPair(p1, p2)) {
        continue;
      }
      // Close pair rejection
      if (CprManager.isActivated()) {
        CprManager.setPair(p1, p2);
        if (CprManager.isClosePair()) {
          continue;
        }
      }
      CprManager.fill();
      PairHistManager.setPair(p1, p2);
      PairHistManager.fill();
    }
  }
}

// process mixed event different tracks
template <typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9>
void processMixedEvent(T1& Collisions,
                       T2& Partition1,
                       T3& Partition2,
                       T4& cache,
                       T5& policy,
                       T6& depth,
                       T7& PairHistManager,
                       T8& CprManager,
                       T9& PcManager)
{
  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {
    if (!(std::fabs(collision1.magField() - collision2.magField()) < o2::constants::math::Epsilon)) {
      continue;
    }
    CprManager.setMagField(collision1.magField());
    auto sliceParticle1 = Partition1->sliceByCached(o2::aod::femtobase::stored::collisionId, collision1.globalIndex(), cache);
    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::collisionId, collision2.globalIndex(), cache);
    if (sliceParticle1.size() == 0 || sliceParticle2.size() == 0) {
      continue;
    }
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(sliceParticle1, sliceParticle2))) {
      // pair cleaning
      if (!PcManager.isCleanPair(p1, p2)) {
        continue;
      }
      // Close pair rejection
      if (CprManager.isActivated()) {
        CprManager.setPair(p1, p2);
        if (CprManager.isClosePair()) {
          continue;
        }
      }
      CprManager.fill();
      PairHistManager.setPair(p1, p2);
      PairHistManager.fill();
    }
  }
}

// process mixed event for track and particles decaying into tracks
template <typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9,
          typename T10>
void processMixedEvent(T1& Collisions,
                       T2& Partition1,
                       T3& Partition2,
                       T4& TrackTable,
                       T5& cache,
                       T6& policy,
                       T7& depth,
                       T8& PairHistManager,
                       T9& CprManager,
                       T10& PcManager)
{
  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {
    if (!(std::fabs(collision1.magField() - collision2.magField()) < o2::constants::math::Epsilon)) {
      continue;
    }
    CprManager.setMagField(collision1.magField());
    auto sliceParticle1 = Partition1->sliceByCached(o2::aod::femtobase::stored::collisionId, collision1.globalIndex(), cache);
    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::collisionId, collision2.globalIndex(), cache);
    if (sliceParticle1.size() == 0 || sliceParticle2.size() == 0) {
      continue;
    }
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(sliceParticle1, sliceParticle2))) {
      // pair cleaning
      if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
        continue;
      }
      // Close pair rejection
      if (CprManager.isActivated()) {
        CprManager.setPair(p1, p2, TrackTable);
        if (CprManager.isClosePair()) {
          continue;
        }
      }
      CprManager.fill();
      PairHistManager.setPair(p1, p2);
      PairHistManager.fill();
    }
  }
}
} // namespace pairprocesshelpers
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_PAIRPROCESSHELPERS_H_
