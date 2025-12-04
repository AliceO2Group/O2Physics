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
    const auto originalIndex = mcCol.globalIndex();
    auto it = mCollisionMap.find(originalIndex);

    // If already filled, stop here
    if (it != mCollisionMap.end()) {
      return;
    }

    // Fill MC collision table
    mcProducts.producedMcCollisions(
      1, 1);
    // mcCol.multMCNParticlesEta08(),
    // mcCol.centFT0M());

    // Record the new index
    const int newIndex = mcProducts.producedMcCollisions.lastIndex();
    mCollisionMap.emplace(originalIndex, newIndex);
  }

  template <modes::System system, typename T1, typename T2, typename T3>
  void fillMcCollisionWithLabel(T1& mcProducts, T2 const& col, T3 const& mcCols)
  {
    if (!mFillAnyTable) {
      return;
    }

    // Case: This reconstructed collision has an MC collision
    if (col.has_mcCollision()) {

      // auto genCol = col.template mcCollision_as<T3>();
      auto genCol = mcCols.rawIteratorAt(col.mcCollisionId() - mcCols.offset());
      const auto originalIndex = genCol.globalIndex();

      // Ensure MC collision exists
      auto it = mCollisionMap.find(originalIndex);
      if (it == mCollisionMap.end()) {
        // Not yet created → create it
        this->fillMcCollision<system>(genCol, mcProducts);
        it = mCollisionMap.find(originalIndex);
      }
      // Add label
      mcProducts.producedCollisionLabels(it->second);
      return;
    }
    // Case: No MC collision associated
    mcProducts.producedCollisionLabels(-1);
  }

  template <typename T1, typename T2, typename T3>
  modes::McOrigin getOrigin(T1 const& col, T2 const& /*mcCols*/, T3& mcParticle)
  {
    // check if reconstructed collision has a generated collision
    if (col.has_mcCollision()) {
      // now check  collision ids,  if  they  not  match,  then the  track  does belong   to  the  wrong collision
      if (col.mcCollisionId() != mcParticle.mcCollisionId()) {
        return modes::McOrigin::kFromWrongCollision;
      }
      if (mcParticle.isPhysicalPrimary()) {
        return modes::McOrigin::kPhysicalPrimary;
        //  getGenStatusCode  ==  -1  means  production  during  transport
        //  getProcess  ==  4  means production through  decay
        //  getProcess  ==  23  means production through  inelastic  hadronic  interaction
      } else if ((mcParticle.getProcess() == 4) && (mcParticle.getGenStatusCode() == -1) && mcParticle.has_mothers()) {
        return modes::McOrigin::kFromSecondaryDecay;
      } else if ((mcParticle.getProcess() == 23) && (mcParticle.getGenStatusCode() == -1) && mcParticle.has_mothers()) {
        return modes::McOrigin::kFromMaterial;
      } else {
        return modes::McOrigin::kFromUnkown;
      }
    }
    return modes::McOrigin::kFromFakeRecoCollision;
  }

  template <modes::System system, typename T1, typename T2, typename T3>
  int64_t getMcColId(T1 const& /*mcCols*/, T2 const& mcParticle, T3& mcProducts)
  {
    auto mcCol = mcParticle.mcCollision(); // <-- no template param
    const auto gid = mcCol.globalIndex();

    // Find or create
    auto it = mCollisionMap.find(gid);
    if (it == mCollisionMap.end()) {
      this->fillMcCollision<system>(mcCol, mcProducts);
      it = mCollisionMap.find(gid);
    }
    return it->second;
    // auto mcCol = mcParticle.template mcCollision_as<T1>();
    // this->fillMcCollision<system>(mcCol, mcProducts);
    // return mCollisionMap.at(mcCol.globalIndex());
  }
  template <modes::System system, modes::Track trackType, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcTrackWithLabel(T1 const& col, T2 const& mcCols, T3 const& track, T4 const& /*mcParticles*/, T5& mcProducts)
  {
    if (!mFillAnyTable) {
      return;
    }

    // No MC particle → label = -1
    if (!track.has_mcParticle()) {
      mcProducts.producedTrackLabels(-1);
      return;
    }

    auto mcParticle = track.template mcParticle_as<T4>();
    auto mcGlob = mcParticle.globalIndex();

    int mcIdx = -1;

    // Already processed?
    auto it = mMcTrackMap.find(mcGlob);
    if (it == mMcTrackMap.end()) {

      // ---- Fill new MC particle entry ----
      auto origin = this->getOrigin(col, mcCols, mcParticle);
      int64_t mcColId = this->getMcColId<system>(mcCols, mcParticle, mcProducts);

      mcProducts.producedMcParticles(
        mcColId,
        static_cast<aod::femtodatatypes::McOriginType>(origin),
        mcParticle.pdgCode(),
        mcParticle.pt(),
        mcParticle.eta(),
        mcParticle.phi());

      mcIdx = mcProducts.producedMcParticles.lastIndex();
      mMcTrackMap.emplace(mcGlob, mcIdx);

      // ---- Mothers ----
      int firstMotherPdg = 0;
      int partonicMotherPdg = 0;

      if (mcParticle.has_mothers()) {
        auto Mothers = mcParticle.template mothers_as<T4>();

        bool firstDone = false;
        for (auto const& m : Mothers) {
          if (!firstDone) {
            firstMotherPdg = m.pdgCode();
            firstDone = true;
          }
          if (isPartonic(m)) {
            partonicMotherPdg = m.pdgCode();
            break;
          }
        }
      }

      mcProducts.producedMothers(firstMotherPdg);
      mcProducts.producedPartonicMothers(partonicMotherPdg);

    } else {
      mcIdx = it->second;
    }

    // ---- Label track ----
    mcProducts.producedTrackLabels(mcIdx);
  }

  // template <modes::System system, modes::Track trackType, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  // void fillMcTrackWithLabel(T1 const& col, T2 const& mcCols, T3 const& track, T4& trackProducts, T5 const& /*mcParticles*/, T6& mcProducts)
  // {
  //   if (!mFillAnyTable) {
  //     return;
  //   }
  //   if (track.has_mcParticle()) {
  //     // check if track has a corresponding mc particle
  //     auto mcParticle = track.template mcParticle_as<T4>();
  //     // check if we filled the mc particle already
  //     if (mMcTrackMap.find(mcParticle.globalIndex()) == mMcTrackMap.end()) {
  //       auto origin = this->getOrigin(col, mcCols, mcParticle);
  //       int64_t mcColId = this->getMcColId<system>(mcCols, mcParticle, mcProducts);
  //       mcProducts.producedMcParticles(
  //         mcColId,
  //         static_cast<aod::femtodatatypes::McOriginType>(origin),
  //         mcParticle.pdgCode(),
  //         mcParticle.pt(),
  //         mcParticle.eta(),
  //         mcParticle.phi());
  //       mMcTrackMap.emplace(mcParticle.globalIndex(), mcProducts.producedMcParticles.lastIndex());
  //       mcProducts.producedTrackLabels(mcProducts.producedMcParticles.lastIndex());
  //       mMcTrackLabelMap.emplace(trackProducts.lastIndex(), mcProducts.producedMcParticles.lastIndex());
  //
  //       if (mcParticle.has_mothers()) {
  //         auto Mothers = mcParticle.template mothers_as<T4>();
  //         bool isFirstMother = true;
  //         bool foundPartonicMother = false;
  //         for (auto const& mother : Mothers) {
  //           if (isFirstMother) {
  //             if (!isPartonic(mother)) {
  //               mcProducts.producedMothers(mother.pdgCode());
  //             } else {
  //               mcProducts.producedMothers(0);
  //             }
  //             isFirstMother = false;
  //           }
  //
  //           if (isPartonic(mother)) {
  //             mcProducts.producedPartonicMothers(mother.pdgCode());
  //             foundPartonicMother = true;
  //             break;
  //           }
  //         }
  //         if (!foundPartonicMother) {
  //           mcProducts.producedPartonicMothers(0);
  //         }
  //       } else {
  //         mcProducts.producedMothers(mcColId, 0);
  //         mcProducts.producedPartonicMothers(mcColId, 0);
  //       }
  //     } else {
  //       // if already filled, recall the index
  //       mcProducts.producedTrackLabels(mMcTrackMap.at(track.globalIndex()));
  //       mMcTrackLabelMap.emplace(trackProducts.lastIndex(), mcProducts.producedMcParticles.lastIndex());
  //     }
  //   } else {
  //     // if no corresponding mc particle exists, we fill empty label
  //     mcProducts.producedTrackLabels(-1);
  //     mMcTrackLabelMap.emplace(trackProducts.lastIndex(), mcProducts.producedMcParticles.lastIndex());
  //   }
  // }

  template <typename T1>
  inline bool isPartonic(T1 const& mcParticle)
  {
    int pdgCode = std::abs(mcParticle.pdgCode());
    return (pdgCode >= 1 && pdgCode <= 6) || (pdgCode == 21);
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
  std::unordered_map<int64_t, int64_t> mMcTrackMap;
  std::unordered_map<int64_t, int64_t> mMcTrackLabelMap;
  std::unordered_map<int64_t, int64_t> mV0Map;
};

} // namespace mcbuilder
//
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_MCBUILDER_H_
