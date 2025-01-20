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

/// \file utilsDerivedData.h
/// \brief Utilities for derived-data creators
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#ifndef PWGHF_UTILS_UTILSDERIVEDDATA_H_
#define PWGHF_UTILS_UTILSDERIVEDDATA_H_

#include <vector>
#include <map>

#include "fairlogger/Logger.h"

#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"
#include "Framework/ASoA.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/DerivedTables.h"

// Macro to store nSigma for prong _id_ with PID hypothesis _hyp_ in an array
#define GET_N_SIGMA_PRONG(_array_, _candidate_, _id_, _hyp_) \
  _array_[0] = _candidate_.nSigTpc##_hyp_##_id_();           \
  _array_[1] = _candidate_.nSigTof##_hyp_##_id_();           \
  _array_[2] = _candidate_.tpcTofNSigma##_hyp_##_id_();

namespace o2::analysis::hf_derived
{
/// Reserve space in the filled table for all entries in the source table.
/// \param cursor  cursor of the filled table
/// \param enabled  switch for filling the table
/// \param size  size of the source table
template <typename T>
void reserveTable(T& cursor, const o2::framework::Configurable<bool>& enabled, const uint64_t size)
{
  if (enabled.value) {
    cursor.reserve(size);
  }
}

struct HfConfigurableDerivedData : o2::framework::ConfigurableGroup {
  // Candidates
  o2::framework::Configurable<bool> fillCandidateBase{"fillCandidateBase", true, "Fill candidate base properties"};
  // Collisions
  o2::framework::Configurable<bool> fillCollBase{"fillCollBase", true, "Fill collision base properties"};
  o2::framework::Configurable<bool> fillCollId{"fillCollId", true, "Fill original collision indices"};
  // MC collisions
  o2::framework::Configurable<bool> fillMcCollBase{"fillMcCollBase", true, "Fill MC collision base properties"};
  o2::framework::Configurable<bool> fillMcCollId{"fillMcCollId", true, "Fill original MC collision indices"};
  o2::framework::Configurable<bool> fillMcRCollId{"fillMcRCollId", true, "Fill indices of saved derived reconstructed collisions matched to saved derived MC collisions"};
  // MC particles
  o2::framework::Configurable<bool> fillParticleBase{"fillParticleBase", true, "Fill MC particle properties"};
  o2::framework::Configurable<bool> fillParticleId{"fillParticleId", true, "Fill original MC indices"};
};

template <
  typename HfBases,
  typename HfCollBases,
  typename HfCollIds,
  typename HfMcCollBases,
  typename HfMcCollIds,
  typename HfMcRCollIds,
  typename HfPBases,
  typename HfPIds>
struct HfProducesDerivedData : o2::framework::ProducesGroup {
  // Candidates
  o2::framework::Produces<HfBases> rowCandidateBase;
  // Collisions
  o2::framework::Produces<HfCollBases> rowCollBase;
  o2::framework::Produces<HfCollIds> rowCollId;
  // MC collisions
  o2::framework::Produces<HfMcCollBases> rowMcCollBase;
  o2::framework::Produces<HfMcCollIds> rowMcCollId;
  o2::framework::Produces<HfMcRCollIds> rowMcRCollId;
  // MC particles
  o2::framework::Produces<HfPBases> rowParticleBase;
  o2::framework::Produces<HfPIds> rowParticleId;

  HfConfigurableDerivedData const* conf;
  std::map<int, std::vector<int>> matchedCollisions; // indices of derived reconstructed collisions matched to the global indices of MC collisions
  std::map<int, bool> hasMcParticles;                // flags for MC collisions with HF particles

  void init(HfConfigurableDerivedData const& c)
  {
    conf = &c;
  }

  template <typename T>
  void reserveTablesCandidates(T size)
  {
    o2::analysis::hf_derived::reserveTable(rowCandidateBase, conf->fillCandidateBase, size);
  }

  template <typename T>
  void reserveTablesColl(T size)
  {
    o2::analysis::hf_derived::reserveTable(rowCollBase, conf->fillCollBase, size);
    o2::analysis::hf_derived::reserveTable(rowCollId, conf->fillCollId, size);
  }

  template <typename T>
  void reserveTablesMcColl(T size)
  {
    o2::analysis::hf_derived::reserveTable(rowMcCollBase, conf->fillMcCollBase, size);
    o2::analysis::hf_derived::reserveTable(rowMcCollId, conf->fillMcCollId, size);
    o2::analysis::hf_derived::reserveTable(rowMcRCollId, conf->fillMcRCollId, size);
  }

  template <typename T>
  void reserveTablesParticles(T size)
  {
    o2::analysis::hf_derived::reserveTable(rowParticleBase, conf->fillParticleBase, size);
    o2::analysis::hf_derived::reserveTable(rowParticleId, conf->fillParticleId, size);
  }

  template <typename T>
  void fillTablesCandidate(const T& candidate, double invMass, double y)
  {
    if (conf->fillCandidateBase.value) {
      rowCandidateBase(
        rowCollBase.lastIndex(),
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        invMass,
        y);
    }
  }

  template <bool isMC, typename T>
  // void fillTablesCollision(const T& collision, int isEventReject, int runNumber)
  void fillTablesCollision(const T& collision)
  {
    if (conf->fillCollBase.value) {
      rowCollBase(
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
    if (conf->fillCollId.value) {
      rowCollId(
        collision.globalIndex());
    }
    if constexpr (isMC) {
      if (conf->fillMcRCollId.value && collision.has_mcCollision()) {
        // Save rowCollBase.lastIndex() at key collision.mcCollisionId()
        LOGF(debug, "Rec. collision %d: Filling derived-collision index %d for MC collision %d", collision.globalIndex(), rowCollBase.lastIndex(), collision.mcCollisionId());
        matchedCollisions[collision.mcCollisionId()].push_back(rowCollBase.lastIndex()); // [] inserts an empty element if it does not exist
      }
    }
  }

  template <typename T>
  void fillTablesMcCollision(const T& mcCollision)
  {
    if (conf->fillMcCollBase.value) {
      rowMcCollBase(
        mcCollision.posX(),
        mcCollision.posY(),
        mcCollision.posZ(),
        mcCollision.centFT0M());
    }
    if (conf->fillMcCollId.value) {
      rowMcCollId(
        mcCollision.globalIndex());
    }
    if (conf->fillMcRCollId.value) {
      // Fill the table with the vector of indices of derived reconstructed collisions matched to mcCollision.globalIndex()
      rowMcRCollId(
        matchedCollisions[mcCollision.globalIndex()]);
    }
  }

  template <typename T, typename U>
  void fillTablesParticle(const T& particle, U mass)
  {
    if (conf->fillParticleBase.value) {
      rowParticleBase(
        rowMcCollBase.lastIndex(),
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecayPtEtaPhi::y(particle.pt(), particle.eta(), mass),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
    if (conf->fillParticleId.value) {
      rowParticleId(
        particle.mcCollisionId(),
        particle.globalIndex());
    }
  }

  template <typename CollisionType, typename ParticleType>
  void preProcessMcCollisions(CollisionType const& mcCollisions,
                              o2::framework::Preslice<ParticleType> const& mcParticlesPerMcCollision,
                              ParticleType const& mcParticles)
  {
    if (!conf->fillMcRCollId.value) {
      return;
    }
    hasMcParticles.clear();
    // Fill MC collision flags
    for (const auto& mcCollision : mcCollisions) {
      auto thisMcCollId = mcCollision.globalIndex();
      auto particlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, thisMcCollId);
      LOGF(debug, "MC collision %d has %d MC particles (preprocess)", thisMcCollId, particlesThisMcColl.size());
      hasMcParticles[thisMcCollId] = (particlesThisMcColl.size() > 0);
    }
  }

  template <typename CollisionType, typename ParticleType, typename TMass>
  void processMcParticles(CollisionType const& mcCollisions,
                          o2::framework::Preslice<ParticleType> const& mcParticlesPerMcCollision,
                          ParticleType const& mcParticles,
                          TMass const massParticle)
  {
    // Fill MC collision properties
    auto sizeTableMcColl = mcCollisions.size();
    reserveTablesMcColl(sizeTableMcColl);
    for (const auto& mcCollision : mcCollisions) {
      auto thisMcCollId = mcCollision.globalIndex();
      auto particlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, thisMcCollId);
      auto sizeTablePart = particlesThisMcColl.size();
      LOGF(debug, "MC collision %d has %d MC particles", thisMcCollId, sizeTablePart);
      // Skip MC collisions without HF particles (and without HF candidates in matched reconstructed collisions if saving indices of reconstructed collisions matched to MC collisions)
      LOGF(debug, "MC collision %d has %d saved derived rec. collisions", thisMcCollId, matchedCollisions[thisMcCollId].size());
      if (sizeTablePart == 0 && (!conf->fillMcRCollId.value || matchedCollisions[thisMcCollId].empty())) {
        LOGF(debug, "Skipping MC collision %d", thisMcCollId);
        continue;
      }
      LOGF(debug, "Filling MC collision %d at derived index %d", thisMcCollId, rowMcCollBase.lastIndex() + 1);
      fillTablesMcCollision(mcCollision);

      // Fill MC particle properties
      reserveTablesParticles(sizeTablePart);
      for (const auto& particle : particlesThisMcColl) {
        fillTablesParticle(particle, massParticle);
      }
    }
  }
};
} // namespace o2::analysis::hf_derived

#endif // PWGHF_UTILS_UTILSDERIVEDDATA_H_
