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

/// \file cascadeBuilder.h
/// \brief cascade builder
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_CASCADEBUILDER_H_
#define PWGCF_FEMTO_CORE_CASCADEBUILDER_H_

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{
namespace cascadebuilder
{

struct ConfCascadeFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massXiMin{"massXiMin", 1.2f, "Minimum Xi mass"};
  o2::framework::Configurable<float> massXiMax{"massXiMax", 1.4f, "Maximum Xi mass"};
  o2::framework::Configurable<float> rejectMassXiMin{"rejectMassXiMin", 1.317f, "Reject Minimum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> rejectMassXiMax{"rejectMassXiMax", 1.325f, "Rejection Maximum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> massOmegaMin{"massOmegaMin", 1.5f, "Minimum Omega mass"};
  o2::framework::Configurable<float> massOmegaMax{"massOmegaMax", 1.9f, "Maximum Omega mass"};
  o2::framework::Configurable<float> rejectMassOmegaMin{"rejectMassOmegaMin", 1.668f, "Reject minimum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> rejectMassOmegaMax{"rejectMassOmegaMax", 1.676f, "Reject maximum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> massLambdaMin{"massLambdaMin", 1.0f, "Minimum Lambda mass"};
  o2::framework::Configurable<float> massLambdaMax{"massLambdaMax", 1.2f, "Maximum Lambda mass"};
};

#define CASCADE_DEFAULT_BITS                                                                                                                               \
  o2::framework::Configurable<std::vector<float>> cascadeCpaMin{"cascadeCpaMin", {0.95f}, "Minimum cosine of pointing angle"};                             \
  o2::framework::Configurable<std::vector<float>> cascadeTransRadMin{"cascadeTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> cascadeDcaDauMax{"cascadeDcaDauMax", {0.25f}, "Maximum DCA between the daughters at decay vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> lambdaCpaMin{"lambdaCpaMin", {0.78f}, "Minimum cosine of pointing angle"};                               \
  o2::framework::Configurable<std::vector<float>> lambdaTransRadMin{"lambdaTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                        \
  o2::framework::Configurable<std::vector<float>> lambdaDcaDauMax{"lambdaDcaDauMax", {0.5f}, "Maximum DCA between the daughters at decay vertex (cm)"};    \
  o2::framework::Configurable<std::vector<float>> lambdaDcaToPvMin{"lambdaDcaToPvMin", {0.3f}, "Minimum DCA between the lambda and primary vertex"};       \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.8f}, "Minimum DCA of the daughters from primary vertex (cm)"};           \
  o2::framework::Configurable<std::vector<float>> dauDcaMin{"dauDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};                \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};    \
  o2::framework::Configurable<std::vector<float>> posDauTpc{"posDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for positive daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> negDauTpc{"negDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for negative daughter tracks"};

struct ConfXiBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("XiBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcPion{"bachelorTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for bachelor tracks"};
};

struct ConfOmegaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcKaon{"bachelorTpcKaon", {5.f}, "Maximum |nsimga_Kaon| TPC for bachelor tracks"};
};

#undef CASCADE_DEFAULT_BITS

#define CASCADE_DEFAULT_SELECTION(defaultMassMin, defaultMassMax, defaultPdgCode)                              \
  o2::framework::Configurable<int> pdgCode{"pdgCode", defaultPdgCode, "Track PDG code"};                       \
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the charge of the Cascade "};                      \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                        \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                      \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                   \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                    \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};                                     \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};        \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Cascade"}; \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Cascade"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::CascadeMaskType> mask{"mask", 0, "Bitmask for cascade selection"};

struct ConfXiSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("XiSelection");
  CASCADE_DEFAULT_SELECTION(1.22, 1.42, 3312)
};

struct ConfOmegaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaSelection");
  CASCADE_DEFAULT_SELECTION(1.57, 1.77, 3334)
};

/// The different selections this task is capable of doing
enum CascadeSels {
  // selections for cascades
  kCascadeCpaMin,      ///< Min. CPA (cosine pointing angle)
  kCascadeDcaDaughMax, ///< Max. DCA of the daughers at decay vertex
  kCascadeTransRadMin, ///< max. transverse radius

  // selection for lambda daughter
  kLambdaCpaMin,      ///< Min. DCA of the lambda daughers at primary vertex
  kLambdaDcaDauMax,   ///< TPC PID for daughters (Pion/Proton)
  kLambdaTransRadMin, ///< Min. number of TPC clusters of daughter
  kLambdaDcaToPvMin,  ///< Min. DCA to primary vertex of daughter lambda

  // selection for bachelor/daugthers
  kDauAbsEtaMax, ///< Min. DCA of the daughers/bachelor at primary vertex
  kDauTpcClsMin, ///< Min. number of TPC clusters of daughters/bachelor
  kDauDcaMin,    ///< TPC Pion PID for negative daughter

  // PID selection for cascade bachelor
  kBachelorTpcPion, ///< TPC Pion PID for bachelor
  kBachelorTpcKaon, ///< TPC Kaon PID for bachelor
                    ///
  // PID selection for lambda daughers
  kPosDauTpc, ///< TPC PID for positive daughter
  kNegDauTpc, ///< TPC PID for negative daughter

  kCascadeSelsMax
};

const char cascadeSelsName[] = "Cascade Selection Object";
const std::unordered_map<CascadeSels, std::string> cascadeSelsToString = {
  {kCascadeCpaMin, "Cascade CPA Min"},
  {kCascadeDcaDaughMax, "Cascade DCA Daughters Max"},
  {kCascadeTransRadMin, "Cascade Transverse Radius Min"},

  {kLambdaCpaMin, "Lambda CPA Min"},
  {kLambdaDcaDauMax, "Lambda DCA Daughter Max"},
  {kLambdaTransRadMin, "Lambda Transverse Radius Min"},
  {kLambdaDcaToPvMin, "Lambda DCA to PV Min"},

  {kDauAbsEtaMax, "Daughter Abs Eta Max"},
  {kDauTpcClsMin, "Daughter TPC Clusters Min"},
  {kDauDcaMin, "Daughter DCA Min"},

  {kBachelorTpcPion, "Bachelor TPC Pion PID"},
  {kBachelorTpcKaon, "Bachelor TPC Kaon PID"},

  {kPosDauTpc, "Positive Daughter TPC PID"},
  {kNegDauTpc, "Negative Daughter TPC PID"},

  {kCascadeSelsMax, "Cascade Selections Max"}};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <modes::Cascade cascadeType>
class CascadeSelection : public BaseSelection<float, o2::aod::femtodatatypes::CascadeMaskType, kCascadeSelsMax>
{
 public:
  CascadeSelection() {}
  virtual ~CascadeSelection() = default;

  template <typename T1, typename T2>
  void configure(T1 const& config, T2 const& filter)
  {
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      mXiMassLowerLimit = filter.massXiMin.value;
      mXiMassUpperLimit = filter.massXiMax.value;
      mOmegaMassLowerLimit = filter.rejectMassOmegaMin.value;
      mOmegaMassUpperLimit = filter.rejectMassOmegaMax.value;
      this->addSelection(config.bachelorTpcPion.value, kBachelorTpcPion, limits::kAbsUpperLimit, true, true);
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      mOmegaMassLowerLimit = filter.massOmegaMin.value;
      mOmegaMassUpperLimit = filter.massOmegaMax.value;
      mXiMassLowerLimit = filter.rejectMassXiMin.value;
      mXiMassUpperLimit = filter.rejectMassXiMax.value;
      this->addSelection(config.bachelorTpcKaon.value, kBachelorTpcKaon, limits::kAbsUpperLimit, true, true);
    }

    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;
    mLambdaMassMin = filter.massLambdaMin.value;
    mLambdaMassMax = filter.massLambdaMax.value;

    this->addSelection(config.posDauTpc.value, kPosDauTpc, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.negDauTpc.value, kNegDauTpc, limits::kAbsUpperLimit, true, true);

    this->addSelection(config.cascadeCpaMin.value, kCascadeCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.cascadeTransRadMin.value, kCascadeTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.cascadeDcaDauMax.value, kCascadeDcaDaughMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.lambdaCpaMin.value, kLambdaCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.lambdaTransRadMin.value, kLambdaTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.lambdaDcaDauMax.value, kLambdaDcaDauMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.lambdaDcaToPvMin.value, kLambdaDcaToPvMin, limits::kLowerLimit, true, true);
    this->addSelection(config.dauAbsEtaMax.value, kDauAbsEtaMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauDcaMin.value, kDauDcaMin, limits::kAbsLowerLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, kDauTpcClsMin, limits::kLowerLimit, true, true);
  };

  template <typename T1, typename T2, typename T3>
  void applySelections(T1 const& cascade, T2 const& /*tracks*/, T3 const& col)
  {
    this->reset();
    // cascade selections
    this->evaluateObservable(kCascadeCpaMin, cascade.casccosPA(col.posX(), col.posY(), col.posZ()));
    this->evaluateObservable(kCascadeDcaDaughMax, cascade.dcacascdaughters());
    this->evaluateObservable(kCascadeTransRadMin, cascade.cascradius());

    // lambda selection
    this->evaluateObservable(kLambdaCpaMin, cascade.v0cosPA(col.posX(), col.posY(), col.posZ()));
    this->evaluateObservable(kLambdaDcaDauMax, cascade.dcaV0daughters());
    this->evaluateObservable(kLambdaTransRadMin, cascade.v0radius());
    this->evaluateObservable(kLambdaDcaToPvMin, cascade.dcav0topv(col.posX(), col.posY(), col.posZ()));

    auto bachelor = cascade.template bachelor_as<T2>();
    auto posDaughter = cascade.template posTrack_as<T2>();
    auto negDaughter = cascade.template negTrack_as<T2>();

    // daughter selections
    std::array<float, 3> etaDaughters = {std::fabs(bachelor.eta()), std::fabs(posDaughter.eta()), std::fabs(negDaughter.eta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaDaughters.begin(), etaDaughters.end()));

    std::array<float, 3> dcaDaughters = {std::hypot(bachelor.dcaXY(), bachelor.dcaZ()),
                                         std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()),
                                         std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ())};
    this->evaluateObservable(kDauDcaMin, *std::min_element(dcaDaughters.begin(), dcaDaughters.end()));

    std::array<float, 3> clustersDaughters = {1.f * bachelor.tpcNClsFound(), 1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauTpcClsMin, *std::min_element(clustersDaughters.begin(), clustersDaughters.end()));

    // bachelor pid selection
    // check both pion and kaon PID for xi and omega
    this->evaluateObservable(kBachelorTpcPion, bachelor.tpcNSigmaPi());
    this->evaluateObservable(kBachelorTpcKaon, bachelor.tpcNSigmaKa());

    // depending on the charge, we check lambda or antilambda hypothesis
    if (cascade.sign() < 0) {
      this->evaluateObservable(kPosDauTpc, posDaughter.tpcNSigmaPr());
      this->evaluateObservable(kNegDauTpc, negDaughter.tpcNSigmaPi());
    } else if (cascade.sign() > 0) {
      this->evaluateObservable(kPosDauTpc, posDaughter.tpcNSigmaPi());
      this->evaluateObservable(kNegDauTpc, negDaughter.tpcNSigmaPr());
    } else {
      LOG(warn) << "Encountered Cascade candidate with 0 charge";
    }

    this->assembleBitmask();
  };

  template <typename T>
  bool checkFilters(const T& cascade) const
  {
    return ((cascade.pt() > mPtMin && cascade.pt() < mPtMax) &&
            (cascade.eta() > mEtaMin && cascade.eta() < mEtaMax) &&
            (cascade.phi() > mPhiMin && cascade.phi() < mPhiMax) &&
            (cascade.mLambda() > mLambdaMassMin && cascade.mLambda() < mLambdaMassMax));
  }

  template <typename T>
  bool checkHypothesis(T const& cascadeCandidate)
  {
    // no need to check PID of the bachelor/daughters here, they are set as minimal cuts
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      return (mXiMassLowerLimit < cascadeCandidate.mXi() && mXiMassUpperLimit > cascadeCandidate.mXi()) &&           // inside xi mass window
             (cascadeCandidate.mOmega() < mOmegaMassLowerLimit || cascadeCandidate.mOmega() > mOmegaMassUpperLimit); // outside omega mass window
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      return (mOmegaMassLowerLimit < cascadeCandidate.mOmega() && mOmegaMassUpperLimit > cascadeCandidate.mOmega()) && // inside omega mass window
             (cascadeCandidate.mXi() < mXiMassLowerLimit || cascadeCandidate.mXi() > mXiMassUpperLimit);               // outside xi mass window
    }
    return false;
  }

 protected:
  float mXiMassLowerLimit = 0.f;
  float mXiMassUpperLimit = 999.f;

  float mOmegaMassLowerLimit = 0.f;
  float mOmegaMassUpperLimit = 999.f;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -0.9f;
  float mEtaMax = 0.9f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;
  float mLambdaMassMin = 1.f;
  float mLambdaMassMax = 1.2f;
};

struct CascadeBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FXis> producedXis;
  o2::framework::Produces<o2::aod::FXiMasks> producedXiMasks;
  o2::framework::Produces<o2::aod::FXiExtras> producedXiExtras;
  o2::framework::Produces<o2::aod::FOmegas> producedOmegas;
  o2::framework::Produces<o2::aod::FOmegaMasks> producedOmegaMasks;
  o2::framework::Produces<o2::aod::FOmegaExtras> producedOmegaExtras;
};

struct ConfCascadeTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeTables");
  o2::framework::Configurable<int> produceXis{"produceXis", -1, "Produce Xis (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceXiMasks{"produceXiMasks", -1, "Produce XiMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceXiExtras{"produceXiExtras", -1, "Produce XiExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOmegas{"produceOmegas", -1, "Produce Omegas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOmegaMasks{"produceOmegaMasks", -1, "Produce OmegaMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOmegaExtras{"produceOmegaExtras", -1, "Produce OmegaExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::Cascade cascadeType>
class CascadeBuilder
{
 public:
  CascadeBuilder() {}
  virtual ~CascadeBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(T1& config, T2& filter, T3& table, T4& initContext)
  {
    mCascadeSelection.configure(config, filter);
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      LOG(info) << "Initialize femto Xi builder...";
      mProduceXis = utils::enableTable("FXis_001", table.produceXis.value, initContext);
      mProduceXiMasks = utils::enableTable("FXiMasks_001", table.produceXiMasks.value, initContext);
      mProduceXiExtras = utils::enableTable("FXiExtras_001", table.produceXiExtras.value, initContext);
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      LOG(info) << "Initialize femto Omega builder...";
      mProduceOmegas = utils::enableTable("FOmegas_001", table.produceOmegas.value, initContext);
      mProduceOmegaMasks = utils::enableTable("FOmegaMasks_001", table.produceOmegaMasks.value, initContext);
      mProduceOmegaExtras = utils::enableTable("FOmegaExtras_001", table.produceOmegaExtras.value, initContext);
    }

    if (mProduceXis || mProduceXiExtras || mProduceXiMasks || mProduceOmegas || mProduceOmegaMasks || mProduceOmegaExtras) {
      mFillAnyTable = true;
      mCascadeSelection.printSelections(cascadeSelsName, cascadeSelsToString);
    } else {
      LOG(info) << "No tables configured";
    }
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillCascades(T1& collisionProducts, T2& trackProducts, T3& cascadeProducts, T4 const& fullCascades, T5 const& fullTracks, T6 const& col, T7& trackBuilder, T8& indexMap)
  {
    if (!mFillAnyTable) {
      return;
    }

    int64_t bachelorIndex = 0;
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& cascade : fullCascades) {
      if (!mCascadeSelection.checkFilters(cascade)) {
        continue;
      }
      mCascadeSelection.applySelections(cascade, fullTracks, col);
      if (mCascadeSelection.passesAllRequiredSelections() && mCascadeSelection.checkHypothesis(cascade)) {

        auto bachelor = cascade.template bachelor_as<T5>();
        auto posDaughter = cascade.template posTrack_as<T5>();
        auto negDaughter = cascade.template negTrack_as<T5>();

        bachelorIndex = trackBuilder.template getDaughterIndex<modes::Track::kCascadeBachelor>(bachelor, trackProducts, collisionProducts, indexMap);
        posDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(posDaughter, trackProducts, collisionProducts, indexMap);
        negDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(negDaughter, trackProducts, collisionProducts, indexMap);

        fillCascade(collisionProducts, cascadeProducts, cascade, col, bachelorIndex, posDaughterIndex, negDaughterIndex);
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillCascade(T1& collisionProducts, T2& cascadeProducts, T3 const& cascade, T4 const& col, int bachelorIndex, int posDaughterIndex, int negDaughterIndex)
  {
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      if (mProduceXis) {
        cascadeProducts.producedXis(collisionProducts.producedCollision.lastIndex(),
                                    cascade.sign() * cascade.pt(),
                                    cascade.eta(),
                                    cascade.phi(),
                                    cascade.mXi(),
                                    bachelorIndex,
                                    posDaughterIndex,
                                    negDaughterIndex);
      }
      if (mProduceXiMasks) {
        cascadeProducts.producedXiMasks(mCascadeSelection.getBitmask());
      }
      if (mProduceXiExtras) {
        cascadeProducts.producedXiExtras(
          cascade.mOmega(),
          cascade.casccosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcacascdaughters(),
          cascade.cascradius(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posY(), col.posY(), col.posZ()));
      }
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      if (mProduceOmegas) {
        cascadeProducts.producedOmegas(collisionProducts.producedCollision.lastIndex(),
                                       cascade.sign() * cascade.pt(),
                                       cascade.eta(),
                                       cascade.phi(),
                                       cascade.mOmega(),
                                       bachelorIndex,
                                       posDaughterIndex,
                                       negDaughterIndex);
      }
      if (mProduceOmegaMasks) {
        cascadeProducts.producedOmegaMasks(mCascadeSelection.getBitmask());
      }
      if (mProduceOmegaExtras) {
        cascadeProducts.producedOmegaExtras(
          cascade.mXi(),
          cascade.casccosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcacascdaughters(),
          cascade.cascradius(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posY(), col.posY(), col.posZ()));
      }
    }
  }

  bool fillAnyTable() { return mFillAnyTable; }

 private:
  CascadeSelection<cascadeType> mCascadeSelection;
  bool mFillAnyTable = false;
  bool mProduceXis = false;
  bool mProduceXiMasks = false;
  bool mProduceXiExtras = false;
  bool mProduceOmegas = false;
  bool mProduceOmegaMasks = false;
  bool mProduceOmegaExtras = false;
};

} // namespace cascadebuilder
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_CASCADEBUILDER_H_
