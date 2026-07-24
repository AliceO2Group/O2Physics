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

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>

#include <TPDGCode.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto::mcbuilder
{

constexpr int ProducedByDecay = 4;

struct ConfMc : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("MonteCarlo");
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "Passthrough all MC collisions and particles"};
  o2::framework::Configurable<bool> findLastPartonicMother{"findLastPartonicMother", true, "If true, the partonic mother will be the first parton directly after the initial collision. If false, the partonic mother will be the last parton before hadronization"};
  o2::framework::Configurable<float> etaAcceptanceMcOnly{"etaAcceptanceMcOnly", 0.8, "For MC ONLY processing. |eta| acceptance for estimating primary track multiplicity"};
};

struct McBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FMcCols> producedMcCollisions;
  o2::framework::Produces<o2::aod::FMcParticles> producedMcParticles;
  o2::framework::Produces<o2::aod::FMcMothers> producedMothers;
  o2::framework::Produces<o2::aod::FMcPartMoths> producedPartonicMothers;
  o2::framework::Produces<o2::aod::FMcMotherLabels> producedMcMotherLabels;

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
  o2::framework::Configurable<int> producedMcMotherLabels{"producedMcMotherLabels", -1, "Produce mother/partonic-mother labels (-1: auto; 0 off; 1 on)"};

  o2::framework::Configurable<int> producedCollisionLabels{"producedCollisionLabels", -1, "Produce MC collision labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedTrackLabels{"producedTrackLabels", -1, "Produce track labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedLambdaLabels{"producedLambdaLabels", -1, "Produce lambda labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedK0shortLabels{"producedK0shortLabels", -1, "Produce k0short labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedSigmaLabels{"producedSigmaLabels", -1, "Produce k0short labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedSigmaPlusLabels{"producedSigmaPlusLabels", -1, "Produce k0short labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedXiLabels{"producedXiLabels", -1, "Produce xi labels (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producedOmegaLabels{"producedOmegaLabels", -1, "Produce omega labels (-1: auto; 0 off; 1 on)"};
};

// filter/selections for mc collision and mc particles

struct ConfMcCollisionFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("McCollisionFilter");
  o2::framework::Configurable<float> vtxZMin{"vtxZMin", -10.f, "Minimum vertex Z position (cm)"};
  o2::framework::Configurable<float> vtxZMax{"vtxZMax", 10.f, "Maximum vertex Z position (cm)"};
  o2::framework::Configurable<float> multMin{"multMin", 0.f, "Minimum multiplicity"};
  o2::framework::Configurable<float> multMax{"multMax", 5000.f, "Maximum multiplicity"};
  o2::framework::Configurable<float> centMin{"centMin", 0.f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> centMax{"centMax", 100.f, "Maximum centrality (multiplicity percentile)"};
};

template <auto& Prefix>
struct ConfMcParticleSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string(Prefix);
  // kinematic cuts for filtering tracks
  o2::framework::Configurable<float> ptMin{"ptMin", 0.2f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", 2212, "Absolute value of PDG code. Set sign of charge to -1 for antiparticle."};
  o2::framework::Configurable<int> chargeSign{"chargeSign", 1, "Particle charge sign: +1 for positive, -1 for negative, 0 for both"};
};

constexpr const char PrefixMcParticleSelection1[] = "McParticleSelection1";
constexpr const char PrefixMcParticleSelection2[] = "McParticleSelection2";

using ConfMcParticleSelection1 = ConfMcParticleSelection<PrefixMcParticleSelection1>;
using ConfMcParticleSelection2 = ConfMcParticleSelection<PrefixMcParticleSelection2>;

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
    mProduceMcMotherLabels = utils::enableTable("FMcMotherLabels", table.producedMcMotherLabels.value, initContext);

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
        mProduceMcMotherLabels ||
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
    mEtaAcceptanceMcOnly = config.etaAcceptanceMcOnly.value;
    mFindLastPartonicMother = config.findLastPartonicMother.value;
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3>
  void fillMcCollisionWithLabel(T1& mcProducts, T2 const& col, T3 const& /*mcCols*/)
  {
    if (!mProduceCollisionLabels) {
      mcProducts.producedCollisionLabels(-1);
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
      }
      // Add label
      mcProducts.producedCollisionLabels(mCollisionMap.at(originalIndex)); // mc collsions has been added so we can now safely retrieve the index
    } else {
      // If no MC collision associated, fill empty label
      mcProducts.producedCollisionLabels(-1);
    }
  }

  template <modes::System system, typename T1, typename T2>
  void fillMcCollision(T1 const& mcCol, T2& mcProducts)
  {
    float centrality = -1;
    float multiplicity = -1;
    if constexpr (modes::isFlagSet(system, modes::System::kPP)) {
      centrality = mcCol.centFT0M();
      multiplicity = mcCol.multMCNParticlesEta08();
    }
    if constexpr (modes::isFlagSet(system, modes::System::kPbPb)) {
      centrality = mcCol.centFT0C();
      multiplicity = mcCol.multMCNParticlesEta08();
    }

    mcProducts.producedMcCollisions(
      mcCol.posZ(),
      multiplicity,
      centrality);
    mCollisionMap.emplace(mcCol.globalIndex(), mcProducts.producedMcCollisions.lastIndex());
  }

  // for mc only
  template <typename T1, typename T2, typename T3>
  void fillMcCollision(T1 const& mcCol, T2 const& mcParticles, T3& mcProducts)
  {
    float centrality = 0;   // no centrality estimator for mc only, so set to 0
    float multiplicity = 0; // no multiplicity estimator for mc only

    // define multiplicity ourselves by counting primary particles for |eta|,0.8
    // this is similar to how define it in data
    for (auto const& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && (std::fabs(mcParticle.eta()) < mEtaAcceptanceMcOnly)) {
        multiplicity += 1;
      }
    }

    mcProducts.producedMcCollisions(
      mcCol.posZ(),
      multiplicity,
      centrality);
    mCollisionMap.emplace(mcCol.globalIndex(), mcProducts.producedMcCollisions.lastIndex());
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4>
  void fillMcParticle(T1 const& mcParticle, T2 const& mcParticles, T3 const& mcCol, T4& mcProducts)
  {
    this->getOrCreateMcParticleRow<system>(mcParticle, mcParticles, mcCol, mcProducts);
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcTrackWithLabel(T1 const& col, T2 const& mcCols, T3 const& track, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceTrackLabels) {
      mcProducts.producedTrackLabels(-1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, track, mcParticles, mcProducts, [](auto& prod, int64_t p) { prod.producedTrackLabels(p); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcLambdaWithLabel(T1 const& col, T2 const& mcCols, T3 const& lambda, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceLambdaLabels) {
      mcProducts.producedLambdaLabels(-1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, lambda, mcParticles, mcProducts, [](auto& prod, int64_t p) { prod.producedLambdaLabels(p); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcK0shortWithLabel(T1 const& col, T2 const& mcCols, T3 const& k0short, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceK0shortLabels) {
      mcProducts.producedK0shortLabels(-1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, k0short, mcParticles, mcProducts, [](auto& prod, int64_t p) { prod.producedK0shortLabels(p); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcSigmaWithLabel(T1 const& col, T2 const& mcCols, T3 const& sigmaDaughter, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceSigmaLabels) {
      mcProducts.producedSigmaLabels(-1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, sigmaDaughter, mcParticles, mcProducts, [](auto& prod, int64_t p) { prod.producedSigmaLabels(p); }, true);
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcSigmaPlusWithLabel(T1 const& col, T2 const& mcCols, T3 const& sigmaPlusDaughter, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceSigmaPlusLabels) {
      mcProducts.producedSigmaPlusLabels(-1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, sigmaPlusDaughter, mcParticles, mcProducts, [](auto& prod, int64_t p) { prod.producedSigmaPlusLabels(p); }, true);
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcXiWithLabel(T1 const& col, T2 const& mcCols, T3 const& xi, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceXiLabels) {
      mcProducts.producedXiLabels(-1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, xi, mcParticles, mcProducts, [](auto& prod, int64_t p) { prod.producedXiLabels(p); });
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillMcOmegaWithLabel(T1 const& col, T2 const& mcCols, T3 const& omega, T4 const& mcParticles, T5& mcProducts)
  {
    if (!mProduceOmegaLabels) {
      mcProducts.producedOmegaLabels(-1);
      return;
    }
    fillMcLabelGeneric<system>(col, mcCols, omega, mcParticles, mcProducts, [](auto& prod, int64_t p) { prod.producedOmegaLabels(p); });
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

  // mc only, then there is only 1 mc collision
  template <typename T>
  void reset(T const& mcParticles)
  {
    mCollisionMap.clear();
    mMcParticleMap.clear();
    mMcParticleMap.reserve(mcParticles.size());
    mMcMotherMap.clear();
    mMcMotherMap.reserve(mcParticles.size());
    mMcPartonicMotherMap.clear();
    mMcPartonicMotherMap.reserve(mcParticles.size());
  }

 private:
  template <typename T1, typename T2, typename T3>
  modes::McOrigin getOrigin(T1 const& col, T2 const& /*mcCols*/, T3 const& mcParticle)
  {
    // whether a particle is misidentified or not can only be checked by qa/pair task later so it is not set here

    // check if reconstructed collision has a generated collision
    if (!col.has_mcCollision()) {
      return modes::McOrigin::kFromWrongCollision;
    }

    // now check collision ids, if they do not match, then the track belongs to another collision
    if (col.mcCollisionId() != mcParticle.mcCollisionId()) {
      return modes::McOrigin::kFromWrongCollision;
    }

    if (mcParticle.isPhysicalPrimary()) {
      return modes::McOrigin::kPhysicalPrimary;
    }

    if (mcParticle.has_mothers() && mcParticle.getProcess() == ProducedByDecay) {
      return modes::McOrigin::kFromSecondaryDecay;
    }

    // not a primary and not from a decay and not from a wrong collision, we label as material
    return modes::McOrigin::kFromMaterial;
  }

  template <typename T1>
  modes::McOrigin getOrigin(T1 const& mcParticle)
  {
    if (mcParticle.isPhysicalPrimary()) {
      return modes::McOrigin::kPhysicalPrimary;
    }
    if (mcParticle.has_mothers() && mcParticle.getProcess() == ProducedByDecay) {
      return modes::McOrigin::kFromSecondaryDecay;
    }
    return modes::McOrigin::kFromMaterial;
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

  /// Mc-only entry point: no reconstructed collision to match against, so origin
  /// is derived purely from the mc particle itself.
  template <modes::System system, typename T1, typename T2, typename T3, typename T4>
  int64_t getOrCreateMcParticleRow(T1 const& mcParticle, T2 const& mcParticles, T3 const& mcCol, T4& mcProducts)
  {
    auto origin = this->getOrigin(mcParticle);
    return this->buildMcParticleRow<system>(mcParticle, mcParticles, mcCol, origin, mcProducts);
  }

  /// Reco-matched entry point: origin is derived by comparing the reconstructed
  /// collision against the mc collision.
  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  int64_t getOrCreateMcParticleRow(T1 const& col, T2 const& mcCols, T3 const& mcParticle, T4 const& mcParticles, T5 const& mcCol, T6& mcProducts)
  {
    auto origin = this->getOrigin(col, mcCols, mcParticle);
    return this->buildMcParticleRow<system>(mcParticle, mcParticles, mcCol, origin, mcProducts);
  }

  /// Find-or-create the FMcParticles row for mcParticle. On first creation, also
  /// resolves the mother / partonic mother (with kinematics) and writes the single
  /// corresponding FMcMotherLabels row, so it always stays in lockstep, one row
  /// per FMcParticles row.
  template <modes::System system, typename T1, typename T2, typename T3, typename T4>
  int64_t buildMcParticleRow(T1 const& mcParticle, T2 const& mcParticles, T3 const& mcCol, modes::McOrigin origin, T4& mcProducts)
  {
    int64_t mcParticleIndex = mcParticle.globalIndex();

    auto itP = mMcParticleMap.find(mcParticleIndex);
    if (itP != mMcParticleMap.end()) {
      return itP->second;
    }

    int64_t mcColId = this->getMcColId<system>(mcCol, mcProducts);

    mcProducts.producedMcParticles(
      mcColId,
      static_cast<datatypes::McOriginType>(origin),
      mcParticle.pdgCode(),
      mcParticle.pt() * utils::signum(mcParticle.pdgCode()),
      mcParticle.eta(),
      mcParticle.phi());

    int64_t mcParticleRow = mcProducts.producedMcParticles.lastIndex();
    mMcParticleMap[mcParticleIndex] = mcParticleRow;

    // --- mother ---
    int64_t mcMotherRow = -1;
    if (mcParticle.has_mothers()) {
      auto mothers = mcParticle.template mothers_as<T2>();
      auto motherParticle = mothers.front();
      auto mcMotherIndex = motherParticle.globalIndex();

      auto itM = mMcMotherMap.find(mcMotherIndex);
      if (itM != mMcMotherMap.end()) {
        mcMotherRow = itM->second;
      } else {
        auto motherOrigin = this->getOrigin(motherParticle);
        mcProducts.producedMothers(
          static_cast<datatypes::McOriginType>(motherOrigin),
          motherParticle.pdgCode(),
          motherParticle.pt() * utils::signum(motherParticle.pdgCode()),
          motherParticle.eta(),
          motherParticle.phi());
        mcMotherRow = mcProducts.producedMothers.lastIndex();
        mMcMotherMap[mcMotherIndex] = mcMotherRow;
      }
    }

    // --- partonic mother ---
    int64_t mcPartonicMotherRow = -1;
    int64_t mcPartonicMotherIndex = mFindLastPartonicMother
                                      ? this->findLastPartonicMother(mcParticle, mcParticles)
                                      : this->findFirstPartonicMother(mcParticle, mcParticles);
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

    // exactly one FMcMotherLabels row per FMcParticles row, written here and only here
    if (mProduceMcMotherLabels) {
      mcProducts.producedMcMotherLabels(mcMotherRow, mcPartonicMotherRow);
    }

    return mcParticleRow;
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
      writeLabels(mcProducts, -1);
      return;
    }

    auto mcParticle = particle.template mcParticle_as<T4>();
    auto mcCol = mcParticle.template mcCollision_as<T2>();

    if (startFromMotherParticle) {
      // in case of e.g. sigmas we do not reconstruct the mother but the daughter, so here we want to start from the mother particle
      auto mcDaughterParticle = particle.template mcParticle_as<T4>();
      if (!mcDaughterParticle.has_mothers()) {
        writeLabels(mcProducts, -1);
        return;
      }
      auto mothersOfDaughter = mcDaughterParticle.template mothers_as<T4>();
      mcParticle = mothersOfDaughter.front();
    }

    int64_t mcParticleRow = this->getOrCreateMcParticleRow<system>(col, mcCols, mcParticle, mcParticles, mcCol, mcProducts);

    writeLabels(mcProducts, mcParticleRow);
  }

  template <typename T1, typename T2>
  int64_t findFirstPartonicMother(const T1& mcParticle, const T2& mcParticles)
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
      for (int i = motherIds[0]; i <= motherIds[1]; i++) {
        allMotherIds.push_back(i);
      }
    } else {
      // Otherwise just use them as given
      for (const int& id : motherIds) {
        allMotherIds.push_back(id);
      }
    }
    // Loop over all mothers
    for (const int& i : allMotherIds) {

      if (i < 0 || i >= mcParticles.size()) {
        continue;
      }
      const auto& mother = mcParticles.iteratorAt(i);
      int pdgAbs = std::abs(mother.pdgCode());
      // Is it a parton? (quark or gluon)
      if (pdgAbs <= PDG_t::kTop || pdgAbs == PDG_t::kGluon) {
        return i; // Found a parton → return index
      }
      // Recurse upward
      int64_t found = this->findFirstPartonicMother(mother, mcParticles);
      if (found != -1) {
        return found;
      }
    }
    // No partonic ancestor found
    return -1;
  }

  template <typename T1, typename T2>
  int64_t findLastPartonicMother(const T1& mcParticle, const T2& mcParticles)
  {
    int64_t lastPartonIndex = -1;
    int64_t currentIndex = mcParticle.globalIndex();
    while (currentIndex >= 0 && currentIndex < mcParticles.size()) {
      const auto& current = mcParticles.iteratorAt(currentIndex);
      if (!current.has_mothers()) {
        break;
      }
      auto motherIds = current.mothersIds();
      int nextIndex = -1;
      const int defaultMotherSize = 2;
      if (motherIds.size() == defaultMotherSize && motherIds[1] > motherIds[0]) {
        nextIndex = motherIds[0];
      } else {
        for (const int& id : motherIds) {
          if (id >= 0 && id < mcParticles.size()) {
            nextIndex = id;
            break;
          }
        }
      }

      if (nextIndex < 0 || nextIndex >= mcParticles.size()) {
        break;
      }
      const auto& mother = mcParticles.iteratorAt(nextIndex);
      int pdgAbs = std::abs(mother.pdgCode());
      int status = std::abs(o2::mcgenstatus::getGenStatusCode(mother.statusCode()));
      bool isParton = (pdgAbs <= PDG_t::kTop || pdgAbs == PDG_t::kGluon);
      const int isBeamParticleLowerLimit = 11;
      const int isBeamParticleUpperLimit = 19;
      bool isBeam = (status >= isBeamParticleLowerLimit && status <= isBeamParticleUpperLimit);
      if (isBeam) {
        return lastPartonIndex;
      }
      if (isParton) {
        lastPartonIndex = nextIndex;
      }
      currentIndex = nextIndex;
    }
    return -1;
  }

  bool mPassThrough = false;
  bool mFindLastPartonicMother = false;
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
  bool mProduceMcMotherLabels = false;

  float mEtaAcceptanceMcOnly = 0.8;

  std::unordered_map<int64_t, int64_t> mCollisionMap;

  std::unordered_map<int64_t, int64_t> mMcParticleMap;
  std::unordered_map<int64_t, int64_t> mMcMotherMap;
  std::unordered_map<int64_t, int64_t> mMcPartonicMotherMap;
};
} // namespace o2::analysis::femto::mcbuilder

#endif // PWGCF_FEMTO_CORE_MCBUILDER_H_
