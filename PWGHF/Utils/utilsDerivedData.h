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

#include <fairlogger/Logger.h>

#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>

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

struct ConfigurableHfDerivedData : o2::framework::ConfigurableGroup {
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

struct ProducesHfDerivedData : o2::framework::ProducesGroup {
  // Candidates
  o2::framework::Produces<o2::aod::HfBases> rowCandidateBase;
  // Collisions
  o2::framework::Produces<o2::aod::HfCollBases> rowCollBase;
  o2::framework::Produces<o2::aod::HfCollIds> rowCollId;
  // MC collisions
  o2::framework::Produces<o2::aod::HfMcCollBases> rowMcCollBase;
  o2::framework::Produces<o2::aod::HfMcCollIds> rowMcCollId;
  o2::framework::Produces<o2::aod::HfMcRCollIds> rowMcRCollId;
  // MC particles
  o2::framework::Produces<o2::aod::HfPBases> rowParticleBase;
  o2::framework::Produces<o2::aod::HfPIds> rowParticleId;

  ConfigurableHfDerivedData const* conf;
  std::map<int, std::vector<int>> matchedCollisions; // indices of derived reconstructed collisions matched to the global indices of MC collisions

  void init(ConfigurableHfDerivedData const& c) {
    conf = &c;
  }

  template <typename T>
  void reserveTablesCandidates(T size) {
    o2::analysis::hf_derived::reserveTable(rowCandidateBase, conf->fillCandidateBase, size);
  }

  template <typename T>
  void reserveTablesColl(T size) {
    o2::analysis::hf_derived::reserveTable(rowCollBase, conf->fillCollBase, size);
    o2::analysis::hf_derived::reserveTable(rowCollId, conf->fillCollId, size);
  }

  template <typename T>
  void reserveTablesMcColl(T size) {
    o2::analysis::hf_derived::reserveTable(rowMcCollBase, conf->fillMcCollBase, size);
    o2::analysis::hf_derived::reserveTable(rowMcCollId, conf->fillMcCollId, size);
    o2::analysis::hf_derived::reserveTable(rowMcRCollId, conf->fillMcRCollId, size);
  }

  template <typename T>
  void reserveTablesParticles(T size) {
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
        mcCollision.posZ());
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
};
} // namespace o2::analysis::hf_derived

#endif // PWGHF_UTILS_UTILSDERIVEDDATA_H_
