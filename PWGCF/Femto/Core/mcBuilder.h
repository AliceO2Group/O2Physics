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

#include "CommonConstants/MathConstants.h"
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
};

struct ConfMcTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("McTables");
  o2::framework::Configurable<int> produceMcCollisions{"produceMcCollisions", -1, "Produce MC collisions (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMcParticles{"produceMcParticles", -1, "Produce MC particles (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMcMothers{"produceMcMothers", -1, "Produce MC mother particles (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMcPartonicMothers{"produceMcPartonicMothers", -1, "Produce MC partonic mother particles (-1: auto; 0 off; 1 on)"};

  o2::framework::Configurable<int> producedCollisionLabels{"producedCollisionLabels", -1, "Produce MC partonic mother particles (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedTrackLabels{"producedTrackLabels", -1, "Produce MC partonic mother particles (-1: auto; 0 off; 1 on)"};
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

    if (mProduceMcCollisions || mProduceMcParticles || mProduceMcMothers || mProduceMcPartonicMothers || mProduceCollisionLabels || mProduceTrackLabels) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mPassThrough = config.passThrough.value;
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2>
  void fillMcCollision(T1 const& mcCol, T2& mcProducts)
  {
    mcProducts.producedMcCollisions(
      mcCol.multMCNParticlesEta08(),
      mcCol.centFT0M());
    mCollisionMap.emplace(mcCol.globalIndex(), mcProducts.producedMcCollisions.lastIndex());
  }

  template <modes::System system, typename T1, typename T2, typename T3>
  void fillMcCollisionWithLabel(T1& mcProducts, T2 const& col, T3 const& /*mcCols*/)
  {
    if (!mFillAnyTable) {
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
      // Case: No MC collision associated
      mcProducts.producedCollisionLabels(-1);
    }
  }

  template <typename T1, typename T2, typename T3>
  modes::McOrigin getOrigin(T1 const& col, T2 const& /*mcCols*/, T3 const& mcParticle)
  {
    // constants
    // const int producedDuringTransport = -1;
    // const int producedByDecay = 4;
    // const int producedByInelaticHadronicScattering = 4;
    // check if reconstructed collision has a generated collision
    if (col.has_mcCollision()) {
      // now check  collision ids,  if  they  not  match,  then the  track  does belong   to  the  wrong collision
      if (col.mcCollisionId() != mcParticle.mcCollisionId()) {
        return modes::McOrigin::kFromWrongCollision;
      }
      if (mcParticle.isPhysicalPrimary()) {
        return modes::McOrigin::kPhysicalPrimary;
        // } else if ((mcParticle.getProcess() == producedByDecay) && (mcParticle.getGenStatusCode() == producedDuringTransport) && mcParticle.has_mothers()) {
      } else if (mcParticle.has_mothers()) {
        return modes::McOrigin::kFromSecondaryDecay;
        // } else if ((mcParticle.getProcess() == producedByInelaticHadronicScattering) && (mcParticle.getGenStatusCode() == producedDuringTransport) && mcParticle.has_mothers()) {
      } else {
        return modes::McOrigin::kFromMaterial;
        // } else {
        //   return modes::McOrigin::kFromUnkown;
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

  template <modes::System system, modes::Track trackType, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcTrackWithLabel(T1 const& col, T2 const& mcCols, T3 const& track, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mFillAnyTable) {
      return;
    }

    if (!track.has_mcParticle()) {
      mcProducts.producedTrackLabels(-1, -1, -1);
      return;
    }

    auto mcParticle = track.template mcParticle_as<T4>();
    auto mcCol = mcParticle.template mcCollision_as<T2>();

    int64_t particleIndex = mcParticle.globalIndex();
    int64_t trackIndex = track.globalIndex();

    int64_t mcParticleRow = -1;

    // MC PARTICLE
    auto itP = mMcTrackToMcParticleMap.find(particleIndex);
    if (itP != mMcTrackToMcParticleMap.end()) {
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
      mMcTrackToMcParticleMap[particleIndex] = mcParticleRow;
    }

    // MOTHERS (fill only if exists)
    int64_t mothersRow = -1;
    auto itM = mMcTrackToMcMotherMap.find(trackIndex);

    if (itM != mMcTrackToMcMotherMap.end()) {
      mothersRow = itM->second;
    } else {

      auto mothers = mcParticle.template mothers_as<T4>();
      bool motherExists = mcParticle.has_mothers() && !mothers.empty();

      if (motherExists) {
        int motherPdg = mothers.front().pdgCode(); // PDG code is ALWAYS valid if the mother exists

        mcProducts.producedMothers(motherPdg);
        mothersRow = mcProducts.producedMothers.lastIndex();
        mMcTrackToMcMotherMap[trackIndex] = mothersRow;
      }
    }

    // PARTONIC MOTHER (fill only if exists)
    int64_t partonicRow = -1;
    auto itPM = mMcTrackToMcPartonicMap.find(trackIndex);

    if (itPM != mMcTrackToMcPartonicMap.end()) {
      partonicRow = itPM->second;
    } else {
      int partIdx = -1;
      if (mcParticle.has_mothers()) {
        partIdx = this->findFirstPartonicMother(mcParticle, mcParticles);
      }

      bool partonicExists = (partIdx >= 0);

      if (partonicExists) {
        int partonicPdg = mcParticles.iteratorAt(partIdx).pdgCode();

        mcProducts.producedPartonicMothers(partonicPdg);
        partonicRow = mcProducts.producedPartonicMothers.lastIndex();
        mMcTrackToMcPartonicMap[trackIndex] = partonicRow;
      }
    }

    // FINAL LABELS
    mcProducts.producedTrackLabels(mcParticleRow, mothersRow, partonicRow);
  }

  template <typename TParticle, typename TContainer>
  int findFirstPartonicMother(const TParticle& p, const TContainer& mcParticles)
  {
    if (!p.has_mothers()) {
      return -1;
    }

    auto motherIds = p.mothersIds();

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
      int apdg = std::abs(mother.pdgCode());

      // Is it a parton? (quark or gluon)
      if (apdg <= PDG_t::kTop || apdg == PDG_t::kGluon) {
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

  void reset()
  {
    mCollisionMap.clear();

    mMcTrackToMcParticleMap.clear();
    mMcTrackToMcMotherMap.clear();
    mMcTrackToMcPartonicMap.clear();

    mV0Map.clear();
  }

 private:
  bool mPassThrough = false;
  bool mFillAnyTable = false;
  bool mProduceMcCollisions = false;
  bool mProduceMcParticles = false;
  bool mProduceMcMothers = false;
  bool mProduceMcPartonicMothers = false;

  bool mProduceCollisionLabels = false;
  bool mProduceTrackLabels = false;

  std::unordered_map<int64_t, int64_t> mCollisionMap;

  std::unordered_map<int64_t, int64_t> mMcTrackToMcParticleMap;
  std::unordered_map<int64_t, int64_t> mMcTrackToMcMotherMap;
  std::unordered_map<int64_t, int64_t> mMcTrackToMcPartonicMap;

  std::unordered_map<int64_t, int64_t> mV0Map;
};

} // namespace mcbuilder
//
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_MCBUILDER_H_
