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

#include <Framework/ASoAHelpers.h>
#include <Framework/Logger.h>

#include <cstdint>
#include <optional>

namespace o2::analysis::femto::pairprocesshelpers
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
  PairHistManager.resetTrackedParticlesPerEvent();
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
        PairHistManager.setPair(p1, p2, TrackTable, Collision);
        break;
      case kOrder21:
        PairHistManager.setPair(p2, p1, TrackTable, Collision);
        break;
      default:
        PairHistManager.setPair(p1, p2, TrackTable, Collision);
    }
    // fill deta-dphi histograms with kstar cutoff
    CprManager.fill(PairHistManager.getKstar());
    // if pair cuts are configured check them before filling
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
      PairHistManager.trackParticlesPerEvent(p1, p2);
    }
  }
  PairHistManager.fillMixingQaSe();
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
          typename T11,
          typename T12>
void processSameEvent(T1 const& SliceParticle,
                      T2 const& TrackTable,
                      T3 const& mcParticles,
                      T4 const& mcMothers,
                      T5 const& mcPartonicMothers,
                      T6 const& Collision,
                      T7 const& mcCollisions,
                      T8& ParticleHistManager,
                      T9& PairHistManager,
                      T10& ParticleCleaner,
                      T11& CprManager,
                      T12& PcManager,
                      PairOrder pairOrder)
{
  PairHistManager.resetTrackedParticlesPerEvent();
  for (auto const& part : SliceParticle) {
    if (!ParticleCleaner.isClean(part, mcParticles, mcMothers, mcPartonicMothers)) {
      continue;
    }
    ParticleHistManager.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle))) {
    // check if particles are clean
    if (!ParticleCleaner.isClean(p1, mcParticles, mcMothers, mcPartonicMothers) ||
        !ParticleCleaner.isClean(p2, mcParticles, mcMothers, mcPartonicMothers)) {
      continue;
    }
    // check if pair is clean
    if (!PcManager.isCleanPair(p1, p2, TrackTable, mcParticles, mcPartonicMothers)) {
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
        PairHistManager.setPairMc(p1, p2, TrackTable, mcParticles, Collision, mcCollisions);
        break;
      case kOrder21:
        PairHistManager.setPairMc(p2, p1, TrackTable, mcParticles, Collision, mcCollisions);
        break;
      default:
        PairHistManager.setPairMc(p1, p2, TrackTable, mcParticles, Collision, mcCollisions);
    }
    // fill deta-dphi histograms with kstar cutoff
    CprManager.fill(PairHistManager.getKstar());
    // if pair cuts are configured check them before filling
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
      PairHistManager.trackParticlesPerEvent(p1, p2);
    }
  }
  PairHistManager.fillMixingQaSe();
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
  PairHistManager.resetTrackedParticlesPerEvent();
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
    PairHistManager.setPair(p1, p2, TrackTable, Collision);
    CprManager.fill(PairHistManager.getKstar());
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
      PairHistManager.trackParticlesPerEvent(p1, p2);
    }
  }
  PairHistManager.fillMixingQaSe();
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
          typename T13,
          typename T14,
          typename T15>
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
                      T12& ParticleCleaner1,
                      T13& ParticleCleaner2,
                      T14& CprManager,
                      T15& PcManager)
{
  PairHistManager.resetTrackedParticlesPerEvent();
  // Fill single particle histograms
  for (auto const& part : SliceParticle1) {
    if (!ParticleCleaner1.isClean(part, mcParticles, mcMothers, mcPartonicMothers)) {
      continue;
    }
    ParticleHistManager1.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle2) {
    if (!ParticleCleaner2.isClean(part, mcParticles, mcMothers, mcPartonicMothers)) {
      continue;
    }
    ParticleHistManager2.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2))) {
    // check if particles are clean
    if (!ParticleCleaner1.isClean(p1, mcParticles, mcMothers, mcPartonicMothers) ||
        !ParticleCleaner2.isClean(p2, mcParticles, mcMothers, mcPartonicMothers)) {
      continue;
    }
    // pair cleaning
    if (!PcManager.isCleanPair(p1, p2, TrackTable, mcParticles, mcPartonicMothers)) {
      continue;
    }
    // Close pair rejection
    CprManager.setPair(p1, p2, TrackTable);
    if (CprManager.isClosePair()) {
      continue;
    }
    PairHistManager.setPairMc(p1, p2, TrackTable, mcParticles, Collision, mcCollisions);
    CprManager.fill(PairHistManager.getKstar());
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
      PairHistManager.trackParticlesPerEvent(p1, p2);
    }
  }
  PairHistManager.fillMixingQaSe();
}
// process same event for identical particles, mc truth only (no track table, no reco collisions)
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
void processSameEvent(T1 const& SliceParticle,
                      T2 const& /*mcParticles*/,
                      T3 const& mcMothers,
                      T4 const& mcPartonicMothers,
                      T5 const& Collision,
                      T6& ParticleHistManager,
                      T7& PairHistManager,
                      T8& ParticleCleaner,
                      T9& CprManager,
                      T10& PcManager,
                      PairOrder pairOrder)
{
  PairHistManager.resetTrackedParticlesPerEvent();
  for (auto const& part : SliceParticle) {
    if (!ParticleCleaner.isClean(part, mcMothers, mcPartonicMothers)) {
      continue;
    }
    ParticleHistManager.fill(part, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle))) {
    if (!ParticleCleaner.isClean(p1, mcMothers, mcPartonicMothers) ||
        !ParticleCleaner.isClean(p2, mcMothers, mcPartonicMothers)) {
      continue;
    }
    if (!PcManager.isCleanPair(p1, p2, mcPartonicMothers)) {
      continue;
    }
    CprManager.setPair(p1, p2);
    if (CprManager.isClosePair()) {
      continue;
    }
    switch (pairOrder) {
      case kOrder12:
        PairHistManager.setPairMcTruth(p1, p2, Collision);
        break;
      case kOrder21:
        PairHistManager.setPairMcTruth(p2, p1, Collision);
        break;
      default:
        PairHistManager.setPairMcTruth(p1, p2, Collision);
    }
    CprManager.fill(PairHistManager.getKstar());
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
      PairHistManager.trackParticlesPerEvent(p1, p2);
    }
  }
  PairHistManager.fillMixingQaSe();
}

// process same event for non-identical particles, mc truth only
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
                      T3 const& /*mcParticles*/,
                      T4 const& mcMothers,
                      T5 const& mcPartonicMothers,
                      T6 const& Collision,
                      T7& ParticleHistManager1,
                      T8& ParticleHistManager2,
                      T9& PairHistManager,
                      T10& ParticleCleaner1,
                      T11& ParticleCleaner2,
                      T12& CprManager,
                      T13& PcManager)
{
  PairHistManager.resetTrackedParticlesPerEvent();
  for (auto const& part : SliceParticle1) {
    if (!ParticleCleaner1.isClean(part, mcMothers, mcPartonicMothers)) {
      continue;
    }
    ParticleHistManager1.fill(part, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle2) {
    if (!ParticleCleaner2.isClean(part, mcMothers, mcPartonicMothers)) {
      continue;
    }
    ParticleHistManager2.fill(part, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2))) {
    if (!ParticleCleaner1.isClean(p1, mcMothers, mcPartonicMothers) ||
        !ParticleCleaner2.isClean(p2, mcMothers, mcPartonicMothers)) {
      continue;
    }
    if (!PcManager.isCleanPair(p1, p2, mcPartonicMothers)) {
      continue;
    }
    CprManager.setPair(p1, p2);
    if (CprManager.isClosePair()) {
      continue;
    }
    PairHistManager.setPairMcTruth(p1, p2, Collision);
    CprManager.fill(PairHistManager.getKstar());
    if (PairHistManager.checkPairCuts()) {
      PairHistManager.template fill<mode>();
      PairHistManager.trackParticlesPerEvent(p1, p2);
    }
  }
  PairHistManager.fillMixingQaSe();
}

// mixed event in data
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
  int64_t lastCollisionIndex = -1;
  int windowSizeRaw = 0;
  int windowSizeEffective = 0;

  // collision1 is fixed across each mixing window, so its track slice is
  // materialized once per window and reused for every mixing partner, instead
  // of being re-sliced (a fresh arrow Slice + selection copy, the dominant cost
  // on the heaviest femto trains) on every (collision1, collision2) pair.
  std::optional<decltype(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, 0, cache))> sliceParticle1;

  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {

    // --- new window ---
    if (collision1.globalIndex() != lastCollisionIndex) {
      if (lastCollisionIndex != -1) {
        PairHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
      }
      windowSizeRaw = 0;
      windowSizeEffective = 0;
      lastCollisionIndex = collision1.globalIndex();
      sliceParticle1.emplace(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache));
    }

    ++windowSizeRaw;

    if (collision1.magField() != collision2.magField()) {
      LOG(warn) << "Tried mixing events with different magnetic field.";
      continue;
    }

    CprManager.setMagField(collision1.magField());

    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache);

    PairHistManager.resetTrackedParticlesPerEvent();

    if (sliceParticle1->size() == 0 || sliceParticle2.size() == 0) {
      PairHistManager.fillMixingQaMePerEvent();
      continue;
    }

    bool hasValidPair = false;
    PairHistManager.fillMixingQaMe(collision1, collision2);
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(*sliceParticle1, sliceParticle2))) {

      if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
        continue;
      }

      CprManager.setPair(p1, p2, TrackTable);
      if (CprManager.isClosePair()) {
        continue;
      }

      PairHistManager.setPair(p1, p2, TrackTable, collision1, collision2);
      CprManager.fill(PairHistManager.getKstar());

      if (PairHistManager.checkPairCuts()) {
        hasValidPair = true;
        PairHistManager.trackParticlesPerEvent(p1, p2);
        PairHistManager.template fill<mode>();
      }
    }

    if (hasValidPair) {
      ++windowSizeEffective;
    }

    PairHistManager.fillMixingQaMePerEvent();
  }

  // --- final window ---
  if (windowSizeRaw > 0) {
    PairHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
  }
}

// process mixed event in mc
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
          typename T13,
          typename T14,
          typename T15,
          typename T16>
void processMixedEvent(T1 const& Collisions,
                       T2 const& mcCollisions,
                       T3& Partition1,
                       T4& Partition2,
                       T5 const& TrackTable,
                       T6 const& mcParticles,
                       T7 const& mcMothers,
                       T8 const& mcPartonicMothers,
                       T9& cache,
                       T10 const& policy,
                       T11 const& depth,
                       T12& PairHistManager,
                       T13& ParticleCleaner1,
                       T14& ParticleCleaner2,
                       T15& CprManager,
                       T16& PcManager)
{
  int64_t lastCollisionIndex = -1;
  int windowSizeRaw = 0;
  int windowSizeEffective = 0;

  // collision1 is fixed across each mixing window, so its track slice is
  // materialized once per window and reused for every mixing partner, instead
  // of being re-sliced (a fresh arrow Slice + selection copy, the dominant cost
  // on the heaviest femto trains) on every (collision1, collision2) pair.
  std::optional<decltype(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, 0, cache))> sliceParticle1;

  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {
    if (collision1.globalIndex() != lastCollisionIndex) {
      if (lastCollisionIndex != -1) {
        PairHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
      }
      windowSizeRaw = 0;
      windowSizeEffective = 0;
      lastCollisionIndex = collision1.globalIndex();
      sliceParticle1.emplace(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache));
    }

    ++windowSizeRaw;

    if (collision1.magField() != collision2.magField()) {
      LOG(warn) << "Tried mixing events with different magnetic field.";
      continue;
    }

    CprManager.setMagField(collision1.magField());

    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache);

    PairHistManager.resetTrackedParticlesPerEvent();

    if (sliceParticle1->size() == 0 || sliceParticle2.size() == 0) {
      PairHistManager.fillMixingQaMePerEvent();
      continue;
    }

    bool hasValidPair = false;
    PairHistManager.fillMixingQaMe(collision1, collision2);

    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(*sliceParticle1, sliceParticle2))) {

      if (!ParticleCleaner1.isClean(p1, mcParticles, mcMothers, mcPartonicMothers) ||
          !ParticleCleaner2.isClean(p2, mcParticles, mcMothers, mcPartonicMothers)) {
        continue;
      }

      if (!PcManager.isCleanPair(p1, p2, TrackTable)) {
        continue;
      }

      CprManager.setPair(p1, p2, TrackTable);
      if (CprManager.isClosePair()) {
        continue;
      }

      PairHistManager.setPairMc(p1, p2, TrackTable, mcParticles, collision1, collision2, mcCollisions);

      CprManager.fill(PairHistManager.getKstar());

      if (PairHistManager.checkPairCuts()) {
        hasValidPair = true;
        PairHistManager.trackParticlesPerEvent(p1, p2);
        PairHistManager.template fill<mode>();
      }
    }

    if (hasValidPair) {
      ++windowSizeEffective;
    }

    PairHistManager.fillMixingQaMePerEvent();
  }

  if (windowSizeRaw > 0) {
    PairHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
  }
}
// process mixed event, mc truth only (no track table, collisions already mc truth so no separate mcCollisions)
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
          typename T13,
          typename T14>
void processMixedEvent(T1 const& Collisions,
                       T2& Partition1,
                       T3& Partition2,
                       T4 const& /*mcParticles*/,
                       T5 const& mcMothers,
                       T6 const& mcPartonicMothers,
                       T7& cache,
                       T8 const& policy,
                       T9 const& depth,
                       T10& PairHistManager,
                       T11& ParticleCleaner1,
                       T12& ParticleCleaner2,
                       T13& CprManager,
                       T14& PcManager)
{
  int64_t lastCollisionIndex = -1;
  int windowSizeRaw = 0;
  int windowSizeEffective = 0;

  std::optional<decltype(Partition1->sliceByCached(o2::aod::femtomcparticle::fMcColId, 0, cache))> sliceParticle1;

  for (auto const& [collision1, collision2] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions)) {

    if (collision1.globalIndex() != lastCollisionIndex) {
      if (lastCollisionIndex != -1) {
        PairHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
      }
      windowSizeRaw = 0;
      windowSizeEffective = 0;
      lastCollisionIndex = collision1.globalIndex();
      sliceParticle1.emplace(Partition1->sliceByCached(o2::aod::femtomcparticle::fMcColId, collision1.globalIndex(), cache));
    }

    ++windowSizeRaw;

    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtomcparticle::fMcColId, collision2.globalIndex(), cache);

    PairHistManager.resetTrackedParticlesPerEvent();

    if (sliceParticle1->size() == 0 || sliceParticle2.size() == 0) {
      PairHistManager.fillMixingQaMePerEvent();
      continue;
    }

    bool hasValidPair = false;
    PairHistManager.fillMixingQaMe(collision1, collision2);

    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(*sliceParticle1, sliceParticle2))) {

      if (!ParticleCleaner1.isClean(p1, mcMothers, mcPartonicMothers) ||
          !ParticleCleaner2.isClean(p2, mcMothers, mcPartonicMothers)) {
        continue;
      }

      if (!PcManager.isCleanPair(p1, p2, mcPartonicMothers)) {
        continue;
      }

      CprManager.setPair(p1, p2);
      if (CprManager.isClosePair()) {
        continue;
      }

      PairHistManager.setPairMcTruth(p1, p2, collision1, collision2);

      CprManager.fill(PairHistManager.getKstar());

      if (PairHistManager.checkPairCuts()) {
        hasValidPair = true;
        PairHistManager.trackParticlesPerEvent(p1, p2);
        PairHistManager.template fill<mode>();
      }
    }

    if (hasValidPair) {
      ++windowSizeEffective;
    }

    PairHistManager.fillMixingQaMePerEvent();
  }

  if (windowSizeRaw > 0) {
    PairHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
  }
}

} // namespace o2::analysis::femto::pairprocesshelpers

#endif // PWGCF_FEMTO_CORE_PAIRPROCESSHELPERS_H_
