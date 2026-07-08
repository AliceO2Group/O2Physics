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

/// \file tripletProcessHelpers.h
/// \brief process functions used in triplet tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TRIPLETPROCESSHELPERS_H_
#define PWGCF_FEMTO_CORE_TRIPLETPROCESSHELPERS_H_

#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoAHelpers.h>
#include <Framework/Logger.h>

#include <cstdint>
#include <optional>

namespace o2::analysis::femto::tripletprocesshelpers
{
enum TripletOrder : uint8_t {
  kOrder123, // no swap
  kOrder213, // swap 1&2: for the case that particle 1 & 2 are the same species, particle 3 is something else
  kOrder132, // swap 2&3
  kOrder321, // swap 1&2&3
};

// process same event for identical 3 particles
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
                      T5& TripletHistManager,
                      T6& CtrManager,
                      T7& TcManager,
                      TripletOrder tripletOrder)
{
  TripletHistManager.resetTrackedParticlesPerEvent();

  for (auto const& part : SliceParticle) {
    ParticleHistManager.template fill<mode>(part, TrackTable);
  }

  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle, SliceParticle))) {

    // check if triplet is clean
    if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable)) {
      continue;
    }

    // check if triplet is close
    CtrManager.setTriplet(p1, p2, p3, TrackTable);
    if (CtrManager.isCloseTriplet()) {
      continue;
    }

    // Randomize pair order if enabled
    switch (tripletOrder) {
      case kOrder123:
        TripletHistManager.setTriplet(p1, p2, p3, Collision);
        break;
      case kOrder213:
        TripletHistManager.setTriplet(p2, p1, p3, Collision);
        break;
      case kOrder132:
        TripletHistManager.setTriplet(p1, p3, p2, Collision);
        break;
      case kOrder321:
        TripletHistManager.setTriplet(p3, p2, p1, Collision);
        break;
      default:
        TripletHistManager.setTriplet(p1, p2, p3, Collision);
        break;
    }

    // fill deta-dphi histograms with q3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());

    // if triplet cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
      TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
    }
  }

  TripletHistManager.fillMixingQaSe();
}

// process same event for identical 2 particles and 1 other particle
template <modes::Mode mode,
          typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename T7, typename T8, typename T9>
void processSameEvent(T1 const& SliceParticle1, // 1&2 have same species
                      T2 const& SliceParticle3,
                      T3 const& TrackTable,
                      T4 const& Collision,
                      T5& ParticleHistManager1,
                      T6& ParticleHistManager3,
                      T7& TripletHistManager,
                      T8& CtrManager,
                      T9& TcManager,
                      TripletOrder tripletOrder)
{
  TripletHistManager.resetTrackedParticlesPerEvent();

  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable);
  }

  for (auto const& p3 : SliceParticle3) {
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle1, SliceParticle1))) {

      // check if triplet is clean
      if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable)) {
        continue;
      }

      // check if triplet is close
      CtrManager.setTriplet(p1, p2, p3, TrackTable);
      if (CtrManager.isCloseTriplet()) {
        continue;
      }

      // Randomize triplet order if enabled
      // only kOrder123 and kOrder213 are meaningful here since particle 1 & 2 are the same species
      switch (tripletOrder) {
        case kOrder213:
          TripletHistManager.setTriplet(p2, p1, p3, Collision);
          break;
        case kOrder123:
        default:
          TripletHistManager.setTriplet(p1, p2, p3, Collision);
          break;
      }

      // fill deta-dphi histograms with q3 cutoff
      CtrManager.fill(TripletHistManager.getQ3());

      // if triplet cuts are configured check them before filling
      if (TripletHistManager.checkTripletCuts()) {
        TripletHistManager.template fill<mode>();
        TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
      }
    }
  }

  TripletHistManager.fillMixingQaSe();
}

// process same event for 3 different particles
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
void processSameEvent(T1 const& SliceParticle1,
                      T2 const& SliceParticle2,
                      T3 const& SliceParticle3,
                      T4 const& TrackTable,
                      T5 const& Collision,
                      T6& ParticleHistManager1,
                      T7& ParticleHistManager2,
                      T8& ParticleHistManager3,
                      T9& TripletHistManager,
                      T10& CtrManager,
                      T11& TcManager)
{
  TripletHistManager.resetTrackedParticlesPerEvent();

  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable);
  }
  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.template fill<mode>(part, TrackTable);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable);
  }

  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2, SliceParticle3))) {

    // check if triplet is clean
    if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable)) {
      continue;
    }

    // check if triplet is close
    CtrManager.setTriplet(p1, p2, p3, TrackTable);
    if (CtrManager.isCloseTriplet()) {
      continue;
    }

    TripletHistManager.setTriplet(p1, p2, p3, Collision);

    // fill deta-dphi histograms with q3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());

    // if triplet cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
      TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
    }
  }

  TripletHistManager.fillMixingQaSe();
}

// process same event for 3 identical particles with mc information
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
                      T9& TripletHistManager,
                      T10& CtrManager,
                      T11& TcManager,
                      TripletOrder tripletOrder)
{
  TripletHistManager.resetTrackedParticlesPerEvent();

  for (auto const& part : SliceParticle) {
    ParticleHistManager.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }

  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle, SliceParticle))) {
    // check if triplet is clean
    if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable, mcParticles, mcPartonicMothers)) {
      continue;
    }
    // check if triplet is close
    CtrManager.setTriplet(p1, p2, p3, TrackTable);
    if (CtrManager.isCloseTriplet()) {
      continue;
    }
    // Randomize triplet order if enabled
    switch (tripletOrder) {
      case kOrder123:
        TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
        break;
      case kOrder213:
        TripletHistManager.setTripletMc(p2, p1, p3, mcParticles, Collision, mcCollisions);
        break;
      case kOrder132:
        TripletHistManager.setTripletMc(p1, p3, p2, mcParticles, Collision, mcCollisions);
        break;
      case kOrder321:
        TripletHistManager.setTripletMc(p3, p2, p1, mcParticles, Collision, mcCollisions);
        break;
      default:
        TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
        break;
    }
    // fill deta-dphi histograms with q3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());
    // if triplet cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
      TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
    }
  }

  TripletHistManager.fillMixingQaSe();
}

// process same event for 2 identical particles and one other with mc information
template <modes::Mode mode,
          typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
          typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13>
void processSameEvent(T1 const& SliceParticle1,
                      T2 const& SliceParticle3,
                      T3 const& TrackTable,
                      T4 const& mcParticles,
                      T5 const& mcMothers,
                      T6 const& mcPartonicMothers,
                      T7 const& Collision,
                      T8 const& mcCollisions,
                      T9& ParticleHistManager1,
                      T10& ParticleHistManager3,
                      T11& TripletHistManager,
                      T12& CtrManager,
                      T13& TcManager,
                      TripletOrder tripletOrder)
{
  TripletHistManager.resetTrackedParticlesPerEvent();

  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }

  for (auto const& p3 : SliceParticle3) {
    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle1, SliceParticle1))) {
      // check if triplet is clean
      if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable, mcParticles, mcPartonicMothers)) {
        continue;
      }
      // check if triplet is close
      CtrManager.setTriplet(p1, p2, p3, TrackTable);
      if (CtrManager.isCloseTriplet()) {
        continue;
      }
      // Randomize triplet order if enabled
      // only kOrder123 and kOrder213 are meaningful here since particle 1 & 2 are the same species
      switch (tripletOrder) {
        case kOrder213:
          TripletHistManager.setTripletMc(p2, p1, p3, mcParticles, Collision, mcCollisions);
          break;
        case kOrder123:
        default:
          TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
          break;
      }
      // fill deta-dphi histograms with q3 cutoff
      CtrManager.fill(TripletHistManager.getQ3());
      // if triplet cuts are configured check them before filling
      if (TripletHistManager.checkTripletCuts()) {
        TripletHistManager.template fill<mode>();
        TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
      }
    }
  }

  TripletHistManager.fillMixingQaSe();
}

// process same event for 3 different particles with mc information
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
                      T3 const& SliceParticle3,
                      T4 const& TrackTable,
                      T5 const& mcParticles,
                      T6 const& mcMothers,
                      T7 const& mcPartonicMothers,
                      T8 const& Collision,
                      T9 const& mcCollisions,
                      T10& ParticleHistManager1,
                      T11& ParticleHistManager2,
                      T12& ParticleHistManager3,
                      T13& TripletHistManager,
                      T14& CtrManager,
                      T15& TcManager)
{
  TripletHistManager.resetTrackedParticlesPerEvent();

  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }

  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(SliceParticle1, SliceParticle2, SliceParticle3))) {
    // check if triplet is clean
    if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable, mcParticles, mcPartonicMothers)) {
      continue;
    }
    // check if triplet is close
    CtrManager.setTriplet(p1, p2, p3, TrackTable);
    if (CtrManager.isCloseTriplet()) {
      continue;
    }
    TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
    // fill deta-dphi histograms with q3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());
    // if triplet cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
      TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
    }
  }

  TripletHistManager.fillMixingQaSe();
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
          typename T10,
          typename T11>
void processMixedEvent(T1 const& Collisions,
                       T2& Partition1,
                       T3& Partition2,
                       T4& Partition3,
                       T5 const& TrackTable,
                       T6& cache,
                       T7 const& policy,
                       T8 const& depth,
                       T9& TripletHistManager,
                       T10& CtrManager,
                       T11& TcManager)
{
  int64_t lastCollisionIndex1 = -1;
  int64_t lastCollisionIndex2 = -1;
  int windowSizeRaw = 0;
  int windowSizeEffective = 0;

  // collision1 stays fixed across the outer mixing window, and collision2 stays fixed across each inner sub-window (it only advances once collision3 wraps
  // back to its start Both slices are therefore materialized once per window/sub-window and reused, instead of being re-sliced, i.e. a fresh arrow Slice + selection copy
  std::optional<decltype(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, 0, cache))> sliceParticle1;
  std::optional<decltype(Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, 0, cache))> sliceParticle2;

  for (auto const& [collision1, collision2, collision3] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions, Collisions)) {

    // outer window
    if (collision1.globalIndex() != lastCollisionIndex1) {
      if (lastCollisionIndex1 != -1) {
        TripletHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
      }
      windowSizeRaw = 0;
      windowSizeEffective = 0;
      lastCollisionIndex1 = collision1.globalIndex();
      lastCollisionIndex2 = -1; // force sliceParticle2 to refresh below
      sliceParticle1.emplace(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache));
    }

    // inner sub-window
    if (collision2.globalIndex() != lastCollisionIndex2) {
      lastCollisionIndex2 = collision2.globalIndex();
      sliceParticle2.emplace(Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache));
    }

    ++windowSizeRaw;

    if (collision1.magField() != collision2.magField() ||
        collision2.magField() != collision3.magField() ||
        collision1.magField() != collision3.magField()) {
      LOG(warn) << "Tried mixing events with different magnetic field.";
      continue;
    }

    CtrManager.setMagField(collision1.magField());

    auto sliceParticle3 = Partition3->sliceByCached(o2::aod::femtobase::stored::fColId, collision3.globalIndex(), cache);

    TripletHistManager.resetTrackedParticlesPerEvent();

    if (sliceParticle1->size() == 0 || sliceParticle2->size() == 0 || sliceParticle3.size() == 0) {
      TripletHistManager.fillMixingQaMePerEvent();
      continue;
    }

    bool hasValidTriplet = false;
    TripletHistManager.fillMixingQaMe(collision1, collision2, collision3);

    for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(*sliceParticle1, *sliceParticle2, sliceParticle3))) {
      // pair cleaning
      if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable)) {
        continue;
      }
      // Close pair rejection
      CtrManager.setTriplet(p1, p2, p3, TrackTable);
      if (CtrManager.isCloseTriplet()) {
        continue;
      }

      TripletHistManager.setTriplet(p1, p2, p3, collision1, collision2, collision3);
      CtrManager.fill(TripletHistManager.getQ3());

      if (TripletHistManager.checkTripletCuts()) {
        hasValidTriplet = true;
        TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
        TripletHistManager.template fill<mode>();
      }
    }

    if (hasValidTriplet) {
      ++windowSizeEffective;
    }

    TripletHistManager.fillMixingQaMePerEvent();
  }

  // final window
  if (windowSizeRaw > 0) {
    TripletHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
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
          typename T13>
void processMixedEvent(T1 const& Collisions,
                       T2 const& mcCollisions,
                       T3& Partition1,
                       T4& Partition2,
                       T5& Partition3,
                       T6 const& TrackTable,
                       T7 const& mcParticles,
                       T8& cache,
                       T9 const& policy,
                       T10 const& depth,
                       T11& TripletHistManager,
                       T12& CtrManager,
                       T13& TcManager)
{
  int64_t lastCollisionIndex1 = -1;
  int64_t lastCollisionIndex2 = -1;
  int windowSizeRaw = 0;
  int windowSizeEffective = 0;

  // collision1 stays fixed across the outer mixing window, and collision2 stays fixed across each inner sub-window (it only advances once collision3 wraps
  // back to its start Both slices are therefore materialized once per window/sub-window and reused, instead of being re-sliced, i.e. a fresh arrow Slice + selection copy
  std::optional<decltype(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, 0, cache))> sliceParticle1;
  std::optional<decltype(Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, 0, cache))> sliceParticle2;

  for (auto const& [collision1, collision2, collision3] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions, Collisions)) {

    // outer window
    if (collision1.globalIndex() != lastCollisionIndex1) {
      if (lastCollisionIndex1 != -1) {
        TripletHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
      }
      windowSizeRaw = 0;
      windowSizeEffective = 0;
      lastCollisionIndex1 = collision1.globalIndex();
      lastCollisionIndex2 = -1; // force sliceParticle2 to refresh below
      sliceParticle1.emplace(Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache));
    }

    // inner sub-window
    if (collision2.globalIndex() != lastCollisionIndex2) {
      lastCollisionIndex2 = collision2.globalIndex();
      sliceParticle2.emplace(Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache));
    }

    ++windowSizeRaw;

    if (collision1.magField() != collision2.magField() ||
        collision2.magField() != collision3.magField() ||
        collision1.magField() != collision3.magField()) {
      LOG(warn) << "Tried mixing events with different magnetic field.";
      continue;
    }

    CtrManager.setMagField(collision1.magField());

    auto sliceParticle3 = Partition3->sliceByCached(o2::aod::femtobase::stored::fColId, collision3.globalIndex(), cache);

    TripletHistManager.resetTrackedParticlesPerEvent();

    if (sliceParticle1->size() == 0 || sliceParticle2->size() == 0 || sliceParticle3.size() == 0) {
      TripletHistManager.fillMixingQaMePerEvent();
      continue;
    }

    bool hasValidTriplet = false;
    TripletHistManager.fillMixingQaMe(collision1, collision2, collision3);

    for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(*sliceParticle1, *sliceParticle2, sliceParticle3))) {
      // pair cleaning
      if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable)) {
        continue;
      }
      // Close pair rejection
      CtrManager.setTriplet(p1, p2, p3, TrackTable);
      if (CtrManager.isCloseTriplet()) {
        continue;
      }

      TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, collision1, collision2, collision3, mcCollisions);
      CtrManager.fill(TripletHistManager.getQ3());

      if (TripletHistManager.checkTripletCuts()) {
        hasValidTriplet = true;
        TripletHistManager.trackParticlesPerEvent(p1, p2, p3);
        TripletHistManager.template fill<mode>();
      }
    }

    if (hasValidTriplet) {
      ++windowSizeEffective;
    }

    TripletHistManager.fillMixingQaMePerEvent();
  }

  // --- final window ---
  if (windowSizeRaw > 0) {
    TripletHistManager.fillMixingQaMePerMixingBin(windowSizeRaw, windowSizeEffective);
  }
}

} // namespace o2::analysis::femto::tripletprocesshelpers

#endif // PWGCF_FEMTO_CORE_TRIPLETPROCESSHELPERS_H_
