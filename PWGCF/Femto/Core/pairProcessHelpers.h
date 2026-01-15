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

#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoAHelpers.h"

namespace o2::analysis::femto
{
namespace pairprocesshelpers
{

enum PairOrder : uint8_t {
  kOrder12,
  kOrder21
};

// process same event for identical particles
template <modes::Mode mode,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7>
void processSameEvent(T1 const& SliceParticle,
                      T2 const& TrackTable,
                      T3 const& Collision,
                      T4& ParticleHistManager,
                      T5& PairHistManager,
                      T6& CprManager,
                      T7& PcManager,
                      PairOrder pairOrder)
{
  for (auto const& part : SliceParticle) {
    ParticleHistManager.template fill<mode>(part, TrackTable);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle))) {
    // check if pair is clean
    if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
      continue;
    }
    // check if pair is close
    CprManager.setPair(p1, p2, TrackTable);
    if (CprManager.isClosePair()) {
      continue;
    }
    // Randomize pair order if enabled
    switch (pairOrder) {
      case kOrder12:
        PairHistManager.setPair(p1, p2, Collision);
        break;
      case kOrder21:
        PairHistManager.setPair(p2, p1, Collision);
        break;
      default:
        PairHistManager.setPair(p1, p2, Collision);
    }
    // fill deta-dphi histograms with kstar cutoff
    CprManager.fill(PairHistManager.getKstar());
    // if pair cuts are configured check them before filling
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
    }
  }
}

// process same event for identical particles with mc information
template <modes::Mode mode,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9,
          typename T10,
          typename T11>
void processSameEvent(T1 const& SliceParticle,
                      T2 const& TrackTable,
                      T3 const& mcParticles,
                      T4 const& mcMothers,
                      T5 const& mcPartonicMothers,
                      T6 const& Collision,
                      T7 const& mcCollisions,
                      T8& ParticleHistManager,
                      T9& PairHistManager,
                      T10& CprManager,
                      T11& PcManager,
                      PairOrder pairOrder)
{
  for (auto const& part : SliceParticle) {
    ParticleHistManager.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle))) {
    // check if pair is clean
    if (!PcManager.isCleanPair(p1, p2, TrackTable, mcPartonicMothers)) {
      continue;
    }
    // check if pair is close
    CprManager.setPair(p1, p2, TrackTable);
    if (CprManager.isClosePair()) {
      continue;
    }
    // Randomize pair order if enabled
    switch (pairOrder) {
      case kOrder12:
        PairHistManager.setPairMc(p1, p2, mcParticles, Collision, mcCollisions);
        break;
      case kOrder21:
        PairHistManager.setPairMc(p2, p1, mcParticles, Collision, mcCollisions);
        break;
      default:
        PairHistManager.setPairMc(p1, p2, mcParticles, Collision, mcCollisions);
    }
    // fill deta-dphi histograms with kstar cutoff
    CprManager.fill(PairHistManager.getKstar());
    // if pair cuts are configured check them before filling
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
    }
  }
}

// process same event for non-identical particles
template <modes::Mode mode,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9>
void processSameEvent(T1 const& SliceParticle1,
                      T2 const& SliceParticle2,
                      T3 const& TrackTable,
                      T4 const& Collision,
                      T5& ParticleHistManager1,
                      T6& ParticleHistManager2,
                      T7& PairHistManager,
                      T8& CprManager,
                      T9& PcManager)
{
  // Fill single particle histograms
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable);
  }
  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.template fill<mode>(part, TrackTable);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2))) {
    // pair cleaning
    if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
      continue;
    }
    // Close pair rejection
    CprManager.setPair(p1, p2, TrackTable);
    if (CprManager.isClosePair()) {
      continue;
    }
    PairHistManager.setPair(p1, p2, Collision);
    CprManager.fill(PairHistManager.getKstar());
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
    }
  }
}

// process same event for non-identical particles with mc information
template <modes::Mode mode,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9,
          typename T10,
          typename T11,
          typename T12,
          typename T13>
void processSameEvent(T1 const& SliceParticle1,
                      T2 const& SliceParticle2,
                      T3 const& TrackTable,
                      T4 const& mcParticles,
                      T5 const& mcMothers,
                      T6 const& mcPartonicMothers,
                      T7 const& Collision,
                      T8 const& mcCollisions,
                      T9& ParticleHistManager1,
                      T10& ParticleHistManager2,
                      T11& PairHistManager,
                      T12& CprManager,
                      T13& PcManager)
{
  // Fill single particle histograms
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2))) {
    // pair cleaning
    if (!PcManager.isCleanPair(p1, p2, TrackTable, mcPartonicMothers)) {
      continue;
    }
    // Close pair rejection
    CprManager.setPair(p1, p2, TrackTable);
    if (CprManager.isClosePair()) {
      continue;
    }
    PairHistManager.setPairMc(p1, p2, mcParticles, Collision, mcCollisions);
    CprManager.fill(PairHistManager.getKstar());
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
    }
  }
}

// process mixed event
template <modes::Mode mode,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9,
          typename T10>
void processMixedEvent(T1 const& Collisions,
                       T2& Partition1,
                       T3& Partition2,
                       T4 const& TrackTable,
                       T5& cache,
                       T6 const& policy,
                       T7 const& depth,
                       T8& PairHistManager,
                       T9& CprManager,
                       T10& PcManager)
{
  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {
    if (collision1.magField() != collision2.magField()) {
      continue;
    }
    CprManager.setMagField(collision1.magField());
    auto sliceParticle1 = Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache);
    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache);
    if (sliceParticle1.size() == 0 || sliceParticle2.size() == 0) {
      continue;
    }
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(sliceParticle1, sliceParticle2))) {
      // pair cleaning
      if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
        continue;
      }
      // Close pair rejection
      CprManager.setPair(p1, p2, TrackTable);
      if (CprManager.isClosePair()) {
        continue;
      }
      PairHistManager.setPair(p1, p2, collision1, collision2);
      CprManager.fill(PairHistManager.getKstar());
      if (PairHistManager.checkPairCuts()) {
        PairHistManager.template fill<mode>();
      }
    }
  }
}

// process mixed event  with mc information
template <modes::Mode mode,
          typename T1,
          typename T2,
          typename T3,
          typename T4,
          typename T5,
          typename T6,
          typename T7,
          typename T8,
          typename T9,
          typename T10,
          typename T12,
          typename T11>
void processMixedEvent(T1 const& Collisions,
                       T2 const& mcCollisions,
                       T3& Partition1,
                       T4& Partition2,
                       T5 const& TrackTable,
                       T6 const& mcParticles,
                       T7& cache,
                       T8 const& policy,
                       T9 const& depth,
                       T10& PairHistManager,
                       T11& CprManager,
                       T12& PcManager)
{
  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {
    if (collision1.magField() != collision2.magField()) {
      continue;
    }
    CprManager.setMagField(collision1.magField());
    auto sliceParticle1 = Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache);
    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache);
    if (sliceParticle1.size() == 0 || sliceParticle2.size() == 0) {
      continue;
    }
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(sliceParticle1, sliceParticle2))) {
      // pair cleaning
      if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
        continue;
      }
      // Close pair rejection
      CprManager.setPair(p1, p2, TrackTable);
      if (CprManager.isClosePair()) {
        continue;
      }
      PairHistManager.setPairMc(p1, p2, mcParticles, collision1, collision2, mcCollisions);
      CprManager.fill(PairHistManager.getKstar());
      if (PairHistManager.checkPairCuts()) {
        PairHistManager.template fill<mode>();
      }
    }
  }
}

} // namespace pairprocesshelpers
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_PAIRPROCESSHELPERS_H_
