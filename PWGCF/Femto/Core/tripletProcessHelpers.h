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
/// \brief process functions used in pair tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TRIPLETPROCESSHELPERS_H_
#define PWGCF_FEMTO_CORE_TRIPLETPROCESSHELPERS_H_

#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/ASoAHelpers.h"

namespace o2::analysis::femto
{
namespace tripletprocesshelpers
{

enum TripletOrder : uint8_t {
  kOrder123, // no swap
  kOrder213, // first two swap 1&2 so we can use them for the case that particle 1 & 2 are the same species, particle 3 is something else
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

    // fill deta-dphi histograms with !3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());

    // if triplet cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
    }
  }
}

// process same event for identical 2 particles and 1 other particle
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
void processSameEvent(T1 const& SliceParticle1, // 1&2 have same species
                      T2 const& SliceParticle3,
                      T3 const& TrackTable,
                      T4 const& Collision,
                      T5& ParticleHistManager1,
                      T6& ParticleHistManager3,
                      T7& TripletHistManager,
                      T8& CtrManager,
                      T9& TcManager,
                      TripletOrder triplerOrder)
{
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable);
  }

  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle1, SliceParticle1, SliceParticle3))) {

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
    switch (triplerOrder) {
      case kOrder123:
        TripletHistManager.setTriplet(p1, p2, p3, Collision);
        break;
      case kOrder213:
        TripletHistManager.setTriplet(p1, p2, p3, Collision);
        break;
      default:
        TripletHistManager.setTriplet(p1, p2, p3, Collision);
        break;
    }

    // fill deta-dphi histograms with !3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());

    // if triplet cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
    }
  }
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
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable);
  }
  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.template fill<mode>(part, TrackTable);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable);
  }

  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle1, SliceParticle2, SliceParticle3))) {

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

    // fill deta-dphi histograms with !3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());

    // if triplet cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
    }
  }
}

// // process same event for 3 identical particles with mc information
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
                      int swapTriplet)
{
  for (auto const& part : SliceParticle) {
    ParticleHistManager.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle, SliceParticle, SliceParticle))) {
    // check if pair is clean
    if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable, mcPartonicMothers)) {
      continue;
    }
    // check if pair is close
    CtrManager.setTriplet(p1, p2, p3, TrackTable);
    if (CtrManager.isCloseTriplet()) {
      continue;
    }
    // Randomize pair order if enabled
    //
    switch (swapTriplet) {
      case 3:
        TripletHistManager.setTripletMc(p1, p3, p2, mcParticles, Collision, mcCollisions);
        break;
      case 2:
        TripletHistManager.setTripletMc(p2, p1, p3, mcParticles, Collision, mcCollisions);
        break;
      case 1:
        TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
        break;
      default:
        TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
        break;
    }
    // float threshold = 0.5f;
    // fill deta-dphi histograms with q3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());
    // if pair cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
    }
  }
}

// process same event for 2 identical particles and one other with mc information
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
                      int swapTriplet)
{
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle1, SliceParticle1, SliceParticle3))) {
    // check if pair is clean
    if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable, mcPartonicMothers)) {
      continue;
    }
    // check if pair is close
    CtrManager.setTriplet(p1, p2, p3, TrackTable);
    if (CtrManager.isCloseTriplet()) {
      continue;
    }
    // Randomize pair order if enabled
    //
    switch (swapTriplet) {
      case 2:
        TripletHistManager.setTripletMc(p2, p1, p3, mcParticles, Collision, mcCollisions);
        break;
      case 1:
        TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
        break;
      default:
        TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
        break;
    }
    // fill deta-dphi histograms with q3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());
    // if pair cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
    }
  }
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
  for (auto const& part : SliceParticle1) {
    ParticleHistManager1.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle2) {
    ParticleHistManager2.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& part : SliceParticle3) {
    ParticleHistManager3.template fill<mode>(part, TrackTable, mcParticles, mcMothers, mcPartonicMothers);
  }
  for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(SliceParticle1, SliceParticle2, SliceParticle3))) {
    // check if pair is clean
    if (!TcManager.isCleanTriplet(p1, p2, p3, TrackTable, mcPartonicMothers)) {
      continue;
    }
    // check if pair is close
    CtrManager.setTriplet(p1, p2, p3, TrackTable);
    if (CtrManager.isCloseTriplet()) {
      continue;
    }
    TripletHistManager.setTripletMc(p1, p2, p3, mcParticles, Collision, mcCollisions);
    // fill deta-dphi histograms with q3 cutoff
    CtrManager.fill(TripletHistManager.getQ3());
    // if pair cuts are configured check them before filling
    if (TripletHistManager.checkTripletCuts()) {
      TripletHistManager.template fill<mode>();
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
  for (auto const& [collision1, collision2, collision3] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions, Collisions)) {
    if (collision1.magField() != collision2.magField() ||
        collision2.magField() != collision3.magField() ||
        collision1.magField() != collision3.magField()) {
      continue;
    }
    CtrManager.setMagField(collision1.magField());
    auto sliceParticle1 = Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache);
    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache);
    auto sliceParticle3 = Partition3->sliceByCached(o2::aod::femtobase::stored::fColId, collision3.globalIndex(), cache);
    if (sliceParticle1.size() == 0 || sliceParticle2.size() == 0 || sliceParticle3.size() == 0) {
      continue;
    }
    for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(sliceParticle1, sliceParticle2, sliceParticle3))) {
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
        TripletHistManager.template fill<mode>();
      }
    }
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
  for (auto const& [collision1, collision2, collision3] : o2::soa::selfCombinations(policy, depth, -1, Collisions, Collisions, Collisions)) {
    if (collision1.magField() != collision2.magField() ||
        collision2.magField() != collision3.magField() ||
        collision1.magField() != collision3.magField()) {
      continue;
    }
    CtrManager.setMagField(collision1.magField());
    auto sliceParticle1 = Partition1->sliceByCached(o2::aod::femtobase::stored::fColId, collision1.globalIndex(), cache);
    auto sliceParticle2 = Partition2->sliceByCached(o2::aod::femtobase::stored::fColId, collision2.globalIndex(), cache);
    auto sliceParticle3 = Partition3->sliceByCached(o2::aod::femtobase::stored::fColId, collision3.globalIndex(), cache);
    if (sliceParticle1.size() == 0 || sliceParticle2.size() == 0 || sliceParticle3.size() == 0) {
      continue;
    }
    for (auto const& [p1, p2, p3] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(sliceParticle1, sliceParticle2, sliceParticle3))) {
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
        TripletHistManager.template fill<mode>();
      }
    }
  }
}

} // namespace tripletprocesshelpers
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_TRIPLETPROCESSHELPERS_H_
