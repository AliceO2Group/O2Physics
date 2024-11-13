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

/// \file HfDerivedData.h
/// \brief Class with helper functions for HF derived data format
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN

#ifndef PWGHF_CORE_HFDERIVEDDATA_H_
#define PWGHF_CORE_HFDERIVEDDATA_H_

#include "Common/Core/RecoDecay.h"

class HfDerivedData
{
 public:
  /// Default constructor
  HfDerivedData() = default;

  /// Default destructor
  ~HfDerivedData() = default;

  template <bool isMC, typename T, typename U, typename V, typename M>
  void fillCollTables(const T& collision, bool fillCollBase, U& collBaseTable, bool fillCollId, V& collIdTable, bool fillMcRCollId, M& matchedCollisions)
  {
    if (fillCollBase) {
      collBaseTable(
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        collision.numContrib(),
        collision.centFT0A(),
        collision.centFT0C(),
        collision.centFT0M(),
        collision.centFV0A(),
        collision.multZeqNTracksPV());
    }
    if (fillCollId) {
      collIdTable(collision.globalIndex());
    }

    if constexpr (isMC) {
      if (fillMcRCollId && collision.has_mcCollision()) {
        LOGF(debug, "Rec. collision %d: Filling derived-collision index %d for MC collision %d", collision.globalIndex(), collBaseTable.lastIndex(), collision.mcCollisionId());
        matchedCollisions[collision.mcCollisionId()].push_back(collBaseTable.lastIndex()); // [] inserts an empty element if it does not exist
      }
    }
  }

  template <typename T, typename U, typename V, typename M, typename N>
  void fillCollMcTables(const T& mcCollision, bool fillMcCollBase, U& mcCollBaseTable, bool fillMcCollId, V& mcCollIdTable, bool fillMcRCollId, M& matchedCollisions, N& rmcRCollIdTable)
  {
    if (fillMcCollBase) {
      mcCollBaseTable(
        mcCollision.posX(),
        mcCollision.posY(),
        mcCollision.posZ());
    }
    if (fillMcCollId) {
      mcCollIdTable(
        mcCollision.globalIndex());
    }
    if (fillMcRCollId) {
      rmcRCollIdTable(
        matchedCollisions[mcCollision.globalIndex()]);
    }
  }

  template <typename T, typename U, typename V>
  void fillCandidateTables(const T& candidate, bool fillCandidateBase, U& candidateBaseTable, V& collBaseTable, double invMass, double y)
  {
    if (fillCandidateBase) {
      candidateBaseTable(
        collBaseTable.lastIndex(), // lastIndex is not marked as const so collBaseTable cannot be passed as const either. Should be ok here though
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        invMass,
        y);
    }
  }

  template <typename T, typename U, typename V, typename M, typename N>
  void fillParticleTables(const T& particle, bool fillParticleBase, U& particleBaseTable, V& mcCollBaseTable, M mass, bool fillParticleId, N& particleIdTable)
  {
    if (fillParticleBase) {
      particleBaseTable(
        mcCollBaseTable.lastIndex(), // lastIndex is not marked as const so mcCollBaseTable cannot be passed as const either. Should be ok here though
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecayPtEtaPhi::y(particle.pt(), particle.eta(), mass),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
    if (fillParticleId) {
      particleIdTable(
        particle.mcCollisionId(),
        particle.globalIndex());
    }
  }

 private:
};

#endif // PWGHF_CORE_HFDERIVEDDATA_H_
