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

/// \file mcBuilder.h
/// \brief monte carlo builder
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_MCBUILDER_H_
#define PWGCF_FEMTO_CORE_MCBUILDER_H_

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "fairlogger/Logger.h"

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{
namespace mcbuilder
{

struct ConfMc : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("MonteCarlo");
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "Passthrough all MC collisions and particles"};
};

struct McBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FMcCols> producedMcCollisions;
  o2::framework::Produces<o2::aod::FMcParticles> producedMcParticles;
  o2::framework::Produces<o2::aod::FMcMothers> producedMothers;
  o2::framework::Produces<o2::aod::FMcPartMoths> producedPartonicMothers;

  o2::framework::Produces<o2::aod::FColLabels> producedCollisionLabels;
  o2::framework::Produces<o2::aod::FTrackLabels> producedTrackLabels;
  o2::framework::Produces<o2::aod::FLambdaLabels> producedLambdaLabels;
  o2::framework::Produces<o2::aod::FK0shortLabels> producedK0shortLabels;
  o2::framework::Produces<o2::aod::FSigmaLabels> producedSigmaLabels;
  o2::framework::Produces<o2::aod::FSigmaPlusLabels> producedSigmaPlusLabels;
  o2::framework::Produces<o2::aod::FXiLabels> producedXiLabels;
  o2::framework::Produces<o2::aod::FOmegaLabels> producedOmegaLabels;
};

struct ConfMcTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("McTables");
  o2::framework::Configurable<int> produceMcCollisions{"produceMcCollisions", -1, "Produce MC collisions (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMcParticles{"produceMcParticles", -1, "Produce MC particles (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMcMothers{"produceMcMothers", -1, "Produce MC mother particles (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMcPartonicMothers{"produceMcPartonicMothers", -1, "Produce MC partonic mother particles (-1: auto; 0 off; 1 on)"};

  o2::framework::Configurable<int> producedCollisionLabels{"producedCollisionLabels", -1, "Produce MC partonic mother particles (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedTrackLabels{"producedTrackLabels", -1, "Produce track labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedLambdaLabels{"producedLambdaLabels", -1, "Produce lambda labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedK0shortLabels{"producedK0shortLabels", -1, "Produce k0short labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedSigmaLabels{"producedSigmaLabels", -1, "Produce k0short labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedSigmaPlusLabels{"producedSigmaPlusLabels", -1, "Produce k0short labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedXiLabels{"producedXiLabels", -1, "Produce xi labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedOmegaLabels{"producedOmegaLabels", -1, "Produce omega labels (-1: auto; 0 off; 1 on)"};
};

class McBuilder
{
 public:
  McBuilder() = default;
  ~McBuilder() = default;

  template <typename T1, typename T2, typename T3>
  void init(T1& config, T2& table, T3& initContext)
  {
    LOG(info) << "Initialize monte carlo builder...";

    mProduceMcCollisions = utils::enableTable("FMcCols_001", table.produceMcCollisions.value, initContext);
    mProduceMcParticles = utils::enableTable("FMcParticles_001", table.produceMcParticles.value, initContext);
    mProduceMcMothers = utils::enableTable("FMcMothers_001", table.produceMcMothers.value, initContext);
    mProduceMcPartonicMothers = utils::enableTable("FMcPartMoths_001", table.produceMcPartonicMothers.value, initContext);

    mProduceCollisionLabels = utils::enableTable("FColLabels", table.producedCollisionLabels.value, initContext);
    mProduceTrackLabels = utils::enableTable("FTrackLabels", table.producedTrackLabels.value, initContext);
    mProduceLambdaLabels = utils::enableTable("FLambdaLabels", table.producedLambdaLabels.value, initContext);
    mProduceK0shortLabels = utils::enableTable("FK0shortLabels", table.producedK0shortLabels.value, initContext);
    mProduceSigmaLabels = utils::enableTable("FSigmaLabels", table.producedSigmaLabels.value, initContext);
    mProduceSigmaPlusLabels = utils::enableTable("FSigmaPlusLabels", table.producedSigmaPlusLabels.value, initContext);
    mProduceXiLabels = utils::enableTable("FXiLabels", table.producedXiLabels.value, initContext);
    mProduceOmegaLabels = utils::enableTable("FOmegaLabels", table.producedOmegaLabels.value, initContext);

    if (mProduceMcCollisions || mProduceCollisionLabels ||
        mProduceMcParticles || mProduceMcMothers || mProduceMcPartonicMothers ||
        mProduceTrackLabels ||
        mProduceLambdaLabels || mProduceK0shortLabels ||
        mProduceSigmaLabels || mProduceSigmaPlusLabels ||
        mProduceXiLabels || mProduceOmegaLabels) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mPassThrough = config.passThrough.value;
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3>
  void fillMcCollisionWithLabel(T1& mcProducts, T2 const& col, T3 const& /*mcCols*/)
  {
    if (!mProduceCollisionLabels) {
      return;
    }
    // Case: This reconstructed collision has an MC collision
    if (col.has_mcCollision()) {
      const auto originalIndex = col.mcCollisionId();
      // check if MC collision already exists in femto table
      auto it = mCollisionMap.find(originalIndex);
      if (it == mCollisionMap.end()) {
        // Not yet created → create it
        auto mcCol = col.template mcCollision_as<T3>();
        this->fillMcCollision<system>(mcCol, mcProducts);
        it = mCollisionMap.find(originalIndex);
      }
      // Add label
      mcProducts.producedCollisionLabels(it->second);
    } else {
      // If no MC collision associated, fill empty label
      mcProducts.producedCollisionLabels(-1);
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcTrackWithLabel(T1 const& col, T2 const& mcCols, T3 const& track, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceTrackLabels) {
      mcProducts.producedTrackLabels(-1, -1, -1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, track, mcParticles, mcProducts, [](auto& prod, int64_t p, int64_t m, int64_t pm) { prod.producedTrackLabels(p, m, pm); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcLambdaWithLabel(T1 const& col, T2 const& mcCols, T3 const& lambda, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceLambdaLabels) {
      mcProducts.producedLambdaLabels(-1, -1, -1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, lambda, mcParticles, mcProducts, [](auto& prod, int64_t p, int64_t m, int64_t pm) { prod.producedLambdaLabels(p, m, pm); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcK0shortWithLabel(T1 const& col, T2 const& mcCols, T3 const& k0short, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceK0shortLabels) {
      mcProducts.producedK0shortLabels(-1, -1, -1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, k0short, mcParticles, mcProducts, [](auto& prod, int64_t p, int64_t m, int64_t pm) { prod.producedK0shortLabels(p, m, pm); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcSigmaWithLabel(T1 const& col, T2 const& mcCols, T3 const& sigmaDaughter, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceSigmaLabels) {
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, sigmaDaughter, mcParticles, mcProducts, [](auto& prod, int64_t p, int64_t m, int64_t pm) { prod.producedSigmaLabels(p, m, pm); }, true);
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcSigmaPlusWithLabel(T1 const& col, T2 const& mcCols, T3 const& sigmaPlusDaughter, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceSigmaPlusLabels) {
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, sigmaPlusDaughter, mcParticles, mcProducts, [](auto& prod, int64_t p, int64_t m, int64_t pm) { prod.producedSigmaPlusLabels(p, m, pm); }, true);
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcXiWithLabel(T1 const& col, T2 const& mcCols, T3 const& xi, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceXiLabels) {
      mcProducts.producedXiLabels(-1, -1, -1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, xi, mcParticles, mcProducts, [](auto& prod, int64_t p, int64_t m, int64_t pm) { prod.producedXiLabels(p, m, pm); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcOmegaWithLabel(T1 const& col, T2 const& mcCols, T3 const& omega, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceOmegaLabels) {
      mcProducts.producedOmegaLabels(-1, -1, -1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, omega, mcParticles, mcProducts, [](auto& prod, int64_t p, int64_t m, int64_t pm) { prod.producedOmegaLabels(p, m, pm); });
  }

  bool fillAnyTable() const { return mFillAnyTable; }

  template <typename T1, typename T2>
  void reset(T1 const& mcCollisions, T2 const& mcParticles)
  {
    mCollisionMap.clear();
    mCollisionMap.reserve(mcCollisions.size());
    mMcParticleMap.clear();
    mMcParticleMap.reserve(mcParticles.size());
    mMcMotherMap.clear();
    mMcMotherMap.reserve(mcParticles.size());
    mMcPartonicMotherMap.clear();
    mMcPartonicMotherMap.reserve(mcParticles.size());
  }

 private:
  template <modes::System system, typename T1, typename T2>
  void fillMcCollision(T1 const& mcCol, T2& mcProducts)
  {
    mcProducts.producedMcCollisions(
      mcCol.multMCNParticlesEta08(),
      mcCol.centFT0M());
    mCollisionMap.emplace(mcCol.globalIndex(), mcProducts.producedMcCollisions.lastIndex());
  }

  template <typename T1, typename T2, typename T3>
  modes::McOrigin getOrigin(T1 const& col, T2 const& /*mcCols*/, T3 const& mcParticle)
  {
    // constants
    const int producedByDecay = 4;
    // check if reconstructed collision has a generated collision
    if (col.has_mcCollision()) {
      // now check  collision ids, if they do not match, then the track belongs to another collision
      if (col.mcCollisionId() != mcParticle.mcCollisionId()) {
        return modes::McOrigin::kFromWrongCollision;
      }
      if (mcParticle.isPhysicalPrimary()) {
        return modes::McOrigin::kPhysicalPrimary;
      } else if (mcParticle.has_mothers() && mcParticle.getProcess() == producedByDecay) {
        return modes::McOrigin::kFromSecondaryDecay;
      } else {
        // not a primary and not from a decay, we label as material
        return modes::McOrigin::kFromMaterial;
      }
    }
    return modes::McOrigin::kFromWrongCollision;
  }

  template <modes::System system, typename T1, typename T2>
  int64_t getMcColId(T1 const& mcCol, T2& mcProducts)
  {
    auto gid = mcCol.globalIndex();
    // Find or create
    auto it = mCollisionMap.find(gid);
    if (it == mCollisionMap.end()) {
      this->fillMcCollision<system>(mcCol, mcProducts);
      it = mCollisionMap.find(gid);
      return it->second;
    }
    return it->second;
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fillMcLabelGeneric(T1 const& col,
                          T2 const& mcCols,
                          T3 const& particle,
                          T4 const& mcParticles,
                          T5& mcProducts,
                          T6 writeLabels,
                          bool startFromMotherParticle = false)
  {
    if (!particle.has_mcParticle()) {
      writeLabels(mcProducts, -1, -1, -1);
      return;
    }

    auto mcParticle = particle.template mcParticle_as<T4>();
    auto mcCol = mcParticle.template mcCollision_as<T2>();

    if (startFromMotherParticle) {
      // in case of e.g. sigmas we do not reconstruct the mother but the daughter, so here we want to start from the mother particle
      auto mcDaughterParticle = particle.template mcParticle_as<T4>();
      if (!mcDaughterParticle.has_mothers()) {
        writeLabels(mcProducts, -1, -1, -1);
        return;
      }
      auto mothersOfDaughter = mcDaughterParticle.template mothers_as<T4>();
      mcParticle = mothersOfDaughter.front();
    }

    int64_t mcParticleIndex = mcParticle.globalIndex();
    int64_t mcParticleRow = -1;

    // MC particle
    auto itP = mMcParticleMap.find(mcParticleIndex);
    if (itP != mMcParticleMap.end()) {
      mcParticleRow = itP->second;
    } else {
      auto origin = this->getOrigin(col, mcCols, mcParticle);
      int64_t mcColId = this->getMcColId<system>(mcCol, mcProducts);

      mcProducts.producedMcParticles(
        mcColId,
        static_cast<aod::femtodatatypes::McOriginType>(origin),
        mcParticle.pdgCode(),
        mcParticle.pt() * utils::signum(mcParticle.pdgCode()),
        mcParticle.eta(),
        mcParticle.phi());

      mcParticleRow = mcProducts.producedMcParticles.lastIndex();
      mMcParticleMap[mcParticleIndex] = mcParticleRow;
    }

    // MC mother
    int64_t mcMotherRow = -1;
    if (mcParticle.has_mothers()) {
      auto mothers = mcParticle.template mothers_as<T4>();
      auto mcMotherIndex = mothers.front().globalIndex();

      auto itM = mMcMotherMap.find(mcMotherIndex);
      if (itM != mMcMotherMap.end()) {
        mcMotherRow = itM->second;
      } else {
        mcProducts.producedMothers(mothers.front().pdgCode());
        mcMotherRow = mcProducts.producedMothers.lastIndex();
        mMcMotherMap[mcMotherIndex] = mcMotherRow;
      }
    }

    // Partonic mother
    int64_t mcPartonicMotherRow = -1;
    auto mcPartonicMotherIndex = this->findFirstPartonicMother(mcParticle, mcParticles);
    if (mcPartonicMotherIndex >= 0) {
      auto itPM = mMcPartonicMotherMap.find(mcPartonicMotherIndex);
      if (itPM != mMcPartonicMotherMap.end()) {
        mcPartonicMotherRow = itPM->second;
      } else {
        mcProducts.producedPartonicMothers(
          mcParticles.iteratorAt(mcPartonicMotherIndex).pdgCode());
        mcPartonicMotherRow = mcProducts.producedPartonicMothers.lastIndex();
        mMcPartonicMotherMap[mcPartonicMotherIndex] = mcPartonicMotherRow;
      }
    }

    writeLabels(mcProducts, mcParticleRow, mcMotherRow, mcPartonicMotherRow);
  }

  template <typename T1, typename T2>
  int findFirstPartonicMother(const T1& mcParticle, const T2& mcParticles)
  {
    if (!mcParticle.has_mothers()) {
      return -1;
    }

    auto motherIds = mcParticle.mothersIds();

    // adapted these checks from MCUtils in PWGEM
    const int defaultMotherSize = 2;
    std::vector<int> allMotherIds;
    if (motherIds.size() == defaultMotherSize && motherIds[1] > motherIds[0]) {
      // keep for now, might be needed later
      // && ((80 < std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) && std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) < 90) || (100 < std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) && std::abs(o2::mcgenstatus::getGenStatusCode(p.statusCode())) < 110))) {
      // Pythia: mother range
      for (int i = motherIds[0]; i <= motherIds[1]; i++)
        allMotherIds.push_back(i);
    } else {
      // Otherwise just use them as given
      for (const int& id : motherIds)
        allMotherIds.push_back(id);
    }

    // Loop over all mothers
    for (const int& i : allMotherIds) {

      if (i < 0 || i >= mcParticles.size())
        continue;

      const auto& mother = mcParticles.iteratorAt(i);
      int pdgAbs = std::abs(mother.pdgCode());

      // Is it a parton? (quark or gluon)
      if (pdgAbs <= PDG_t::kTop || pdgAbs == PDG_t::kGluon) {
        return i; // Found a parton → return index
      }

      // Recurse upward
      int found = this->findFirstPartonicMother(mother, mcParticles);
      if (found != -1)
        return found;
    }

    // No partonic ancestor found
    return -1;
  }

  bool mPassThrough = false;
  bool mFillAnyTable = false;
  bool mProduceMcCollisions = false;
  bool mProduceMcParticles = false;
  bool mProduceMcMothers = false;
  bool mProduceMcPartonicMothers = false;

  bool mProduceCollisionLabels = false;
  bool mProduceTrackLabels = false;
  bool mProduceLambdaLabels = false;
  bool mProduceK0shortLabels = false;
  bool mProduceSigmaLabels = false;
  bool mProduceSigmaPlusLabels = false;
  bool mProduceXiLabels = false;
  bool mProduceOmegaLabels = false;

  std::unordered_map<int64_t, int64_t> mCollisionMap;

  std::unordered_map<int64_t, int64_t> mMcParticleMap;
  std::unordered_map<int64_t, int64_t> mMcMotherMap;
  std::unordered_map<int64_t, int64_t> mMcPartonicMotherMap;
};

} // namespace mcbuilder
//
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_MCBUILDER_H_
