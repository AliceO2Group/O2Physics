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

/// \file v0Builder.h
/// \brief v0 builder
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_V0BUILDER_H_
#define PWGCF_FEMTO_CORE_V0BUILDER_H_

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
namespace v0builder
{

// filters applied in the producer task
struct ConfV0Filters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("V0Filters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMinLambda{"massMinLambda", 1.f, "Minimum mass for Lambda hypothesis"};
  o2::framework::Configurable<float> massMaxLambda{"massMaxLambda", 1.2f, "Maximum mass for Lambda hypothesis"};
  o2::framework::Configurable<float> massMinK0short{"massMinK0short", 0.45f, "Minimum mass for K0Short hypothesis"};
  o2::framework::Configurable<float> massMaxK0short{"massMaxK0short", 0.53f, "Maximum mass for K0Short hypothesis"};
  o2::framework::Configurable<float> rejectMassMinLambda{"rejectMassMinLambda", 1.11f, "Minimum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMaxLambda{"rejectMassMaxLambda", 1.12f, "Maximum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMinK0short{"rejectMassMinK0short", 0.48f, "Minimum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMaxK0short{"rejectMassMaxK0short", 0.5f, "Maximum mass to rejection K0short hypothesis for Lambda candidates"};
};

// selections bits for all v0s
#define V0_DEFAULT_BITS                                                                                                                                          \
  o2::framework::Configurable<std::vector<float>> dcaDauMax{"dcaDauMax", {1.5f}, "Maximum DCA between the daughters at decay vertex (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> cpaMin{"cpaMin", {0.99f}, "Minimum cosine of pointing angle"};                                                 \
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2f}, "Minimum transverse radius (cm)"};                                          \
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};                                         \
  o2::framework::Configurable<std::vector<float>> decayVtxMax{"decayVtxMax", {100.f}, "Maximum distance in x,y,z of the decay vertex from primary vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.8f}, "Maximum |eta| for daughter tracks"};                                     \
  o2::framework::Configurable<std::vector<float>> dauDcaMin{"dauDcaMin", {0.05f}, "Minimum DCA of the daughters from primary vertex (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};

// derived selection bits for lambda
struct ConfLambdaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("LambdaBits");
  V0_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDauTpcProton{"posDauTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcProton{"negDauTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC negative for daughter tracks"};
};

// derived selection bits for K0Short
struct ConfK0shortBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("K0shortBits");
  V0_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
};

#undef V0_DEFAULT_BITS

// base selection for analysis task for v0s
#define V0_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                                             \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", defaultPdgCode, "PDG code. Set sign to -1 for antiparticle"}; \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                   \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                                 \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                              \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                               \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};                                                \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                   \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Lambda"};             \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Lambda"};             \
  o2::framework::Configurable<o2::aod::femtodatatypes::V0MaskType> mask{"mask", 0, "Bitmask for v0 selection"};

// base selection for analysis task for lambdas
template <const char* Prefix>
struct ConfLambdaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_SELECTIONS(1.0, 1.2, 3122)
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Lambda (+1: Lambda; -1: Antilambda; 0: both)"};
};

// base selection for analysis task for k0short
template <const char* Prefix>
struct ConfK0shortSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_SELECTIONS(0.47, 0.51, 310)
  o2::framework::Configurable<int> sign{"sign", 0, "Dummy value. For compatability with Lambda selection"};
};

#undef V0_DEFAULT_SELECTIONS

constexpr const char PrefixLambdaSelection1[] = "LambdaSelection1";
constexpr const char PrefixLambdaSelection2[] = "LambdaSelection2";
using ConfLambdaSelection1 = ConfLambdaSelection<PrefixLambdaSelection1>;
using ConfLambdaSelection2 = ConfLambdaSelection<PrefixLambdaSelection2>;
constexpr const char PrefixK0shortSelection1[] = "K0shortSelection1";
constexpr const char PrefixK0shortSelection2[] = "K0shortSelection2";
using ConfK0shortSelection1 = ConfK0shortSelection<PrefixK0shortSelection1>;
using ConfK0shortSelection2 = ConfK0shortSelection<PrefixK0shortSelection2>;

/// The different selections for v0s
enum V0Sels {
  // selections for lambdas
  kCpaMin,      ///< Min. CPA (cosine pointing angle)
  kDcaDaughMax, ///< Max. DCA of the daughters at decay vertex
  kDecayVtxMax, ///< Max. distance of decay vertex in x,y,z
  kTransRadMin, ///< Min. transverse radius
  kTransRadMax, ///< max. transverse radius

  // selection for daughter
  kDauAbsEtaMax, ///< Max. absolute pseudo rapidity
  kDauDcaMin,    ///< Min. DCA of the daughters at primary vertex
  kDauTpcClsMin, ///< Min. number of TPC clusters of daughter

  // pid selection for daughters
  kPosDaughTpcPion,   ///< TPC Pion PID for positive daughter
  kPosDaughTpcProton, ///< TPC Proton PID for positive daughter
  kNegDaughTpcPion,   ///< TPC Pion PID for negative daughter
  kNegDaughTpcProton, ///< TPC Proton PID for negative daughter

  kV0SelsMax
};

constexpr char LambdaSelHistName[] = "hLambdaSelection";
constexpr char AntilambdaSelHistName[] = "hAntiLambdaSelection";
constexpr char K0shortSelHistName[] = "hK0shortSelection";
constexpr char V0SelsName[] = "V0 selection object";
const std::unordered_map<V0Sels, std::string> v0SelectionNames = {
  {kCpaMin, "Min. CPA (cosine pointing angle)"},
  {kDcaDaughMax, "Max. DCA of the daughters at decay vertex"},
  {kDecayVtxMax, "Max. distance of decay vertex in x,y,z"},
  {kTransRadMin, "Min. transverse radius"},
  {kTransRadMax, "Max. transverse radius"},

  {kDauAbsEtaMax, "Max. absolute pseudo rapidity of daughters"},
  {kDauDcaMin, "Min. DCA of the daughters at primary vertex"},
  {kDauTpcClsMin, "Min. number of TPC clusters of daughters"},

  {kPosDaughTpcPion, "TPC Pion PID for positive daughter"},
  {kPosDaughTpcProton, "TPC Proton PID for positive daughter"},
  {kNegDaughTpcPion, "TPC Pion PID for negative daughter"},
  {kNegDaughTpcProton, "TPC Proton PID for negative daughter"}};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <modes::V0 v0Type, const char* HistName>
class V0Selection : public BaseSelection<float, o2::aod::femtodatatypes::V0MaskType, kV0SelsMax>
{
 public:
  V0Selection() = default;
  ~V0Selection() = default;

  template <typename T1, typename T2>
  void configure(o2::framework::HistogramRegistry* registry, T1& config, T2& filter)
  {
    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;

    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda) || modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      mMassLambdaLowerLimit = filter.massMinLambda.value;
      mMassLambdaUpperLimit = filter.massMaxLambda.value;
      mMassK0shortLowerLimit = filter.rejectMassMinK0short.value;
      mMassK0shortUpperLimit = filter.rejectMassMaxK0short.value;

      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        this->addSelection(kPosDaughTpcProton, v0SelectionNames.at(kPosDaughTpcProton), config.posDauTpcProton.value, limits::kAbsUpperLimit, true, true, false);
        this->addSelection(kNegDaughTpcPion, v0SelectionNames.at(kNegDaughTpcPion), config.negDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
      }

      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        this->addSelection(kPosDaughTpcPion, v0SelectionNames.at(kPosDaughTpcPion), config.posDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
        this->addSelection(kNegDaughTpcProton, v0SelectionNames.at(kNegDaughTpcProton), config.negDauTpcProton.value, limits::kAbsUpperLimit, true, true, false);
      }
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      mMassK0shortLowerLimit = filter.massMinK0short.value;
      mMassK0shortUpperLimit = filter.massMaxK0short.value;
      mMassLambdaLowerLimit = filter.rejectMassMinLambda.value;
      mMassLambdaUpperLimit = filter.rejectMassMaxLambda.value;
      this->addSelection(kPosDaughTpcPion, v0SelectionNames.at(kPosDaughTpcPion), config.posDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
      this->addSelection(kNegDaughTpcPion, v0SelectionNames.at(kNegDaughTpcPion), config.negDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
    }

    this->addSelection(kDcaDaughMax, v0SelectionNames.at(kDcaDaughMax), config.dcaDauMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kCpaMin, v0SelectionNames.at(kCpaMin), config.cpaMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTransRadMin, v0SelectionNames.at(kTransRadMin), config.transRadMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTransRadMax, v0SelectionNames.at(kTransRadMax), config.transRadMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kDauAbsEtaMax, v0SelectionNames.at(kDauAbsEtaMax), config.dauAbsEtaMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kDauDcaMin, v0SelectionNames.at(kDauDcaMin), config.dauDcaMin.value, limits::kAbsLowerFunctionLimit, true, true, false);
    this->addSelection(kDauTpcClsMin, v0SelectionNames.at(kDauTpcClsMin), config.dauTpcClustersMin.value, limits::kLowerLimit, true, true, false);

    this->setupContainers<HistName>(registry);
  }

  template <typename T1, typename T2>
  void applySelections(T1 const& v0candidate, T2 const& /*tracks*/)
  {
    this->reset();
    // v0 selections
    this->evaluateObservable(kCpaMin, v0candidate.v0cosPA());
    this->evaluateObservable(kDcaDaughMax, v0candidate.dcaV0daughters());
    // for decay vertex, the x,y and z coordinate have to be below a certain threshold
    // compare the largest of the 3 to the limit set by the bit
    std::array<float, 3> decayCoordinates = {std::fabs(v0candidate.x()), std::fabs(v0candidate.y()), std::fabs(v0candidate.z())};
    this->evaluateObservable(kDecayVtxMax, *std::max_element(decayCoordinates.begin(), decayCoordinates.end()));
    this->evaluateObservable(kTransRadMin, v0candidate.v0radius());
    this->evaluateObservable(kTransRadMax, v0candidate.v0radius());

    // daughter selection
    // for daughter selections, both have to fit the same track quality selection, so we store only one bit for both
    // take largest/smallest from both daughters and evaluate the observable with this value
    auto posDaughter = v0candidate.template posTrack_as<T2>();
    auto negDaughter = v0candidate.template negTrack_as<T2>();

    std::array<float, 2> etaDaughters = {std::fabs(posDaughter.eta()), std::fabs(negDaughter.eta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaDaughters.begin(), etaDaughters.end()));

    std::array<float, 2> dcaDaughters = {std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()), std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ())};
    this->evaluateObservable(kDauDcaMin, *std::min_element(dcaDaughters.begin(), dcaDaughters.end()));

    std::array<float, 2> clustersDaughters = {1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauTpcClsMin, *std::min_element(clustersDaughters.begin(), clustersDaughters.end()));

    // daughter pid selections
    this->evaluateObservable(kPosDaughTpcPion, posDaughter.tpcNSigmaPi());
    this->evaluateObservable(kPosDaughTpcProton, posDaughter.tpcNSigmaPr());
    this->evaluateObservable(kNegDaughTpcPion, negDaughter.tpcNSigmaPi());
    this->evaluateObservable(kNegDaughTpcProton, negDaughter.tpcNSigmaPr());

    this->assembleBitmask<HistName>();
  }

  template <typename T>
  bool checkFilters(const T& v0) const
  {
    // check kinematics first
    const bool kinematicsOk =
      (v0.pt() > mPtMin && v0.pt() < mPtMax) &&
      (v0.eta() > mEtaMin && v0.eta() < mEtaMax) &&
      (v0.phi() > mPhiMin && v0.phi() < mPhiMax);
    if (!kinematicsOk) {
      return false;
    }
    // now check mass hypothesis
    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
      return (v0.mLambda() > mMassLambdaLowerLimit && v0.mLambda() < mMassLambdaUpperLimit) &&   // inside Λ
             (v0.mK0Short() < mMassK0shortLowerLimit || v0.mK0Short() > mMassK0shortUpperLimit); // outside K0s
    }

    if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      return (v0.mAntiLambda() > mMassLambdaLowerLimit && v0.mAntiLambda() < mMassLambdaUpperLimit) && // inside Λbar
             (v0.mK0Short() < mMassK0shortLowerLimit || v0.mK0Short() > mMassK0shortUpperLimit);       // outside K0s
    }

    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      return (v0.mK0Short() > mMassK0shortLowerLimit && v0.mK0Short() < mMassK0shortUpperLimit) &&   // inside K0s
             (v0.mLambda() < mMassLambdaLowerLimit || v0.mLambda() > mMassLambdaUpperLimit) &&       // outside Λ
             (v0.mAntiLambda() < mMassLambdaLowerLimit || v0.mAntiLambda() > mMassLambdaUpperLimit); // outside Λbar
    }
    return false;
  }

 protected:
  float mMassK0shortLowerLimit = 0.483f;
  float mMassK0shortUpperLimit = 0.503f;

  float mMassLambdaLowerLimit = 1.105f;
  float mMassLambdaUpperLimit = 1.125f;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -1.f;
  float mEtaMax = 1.f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;
};

struct V0BuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FLambdas> producedLambdas;
  o2::framework::Produces<o2::aod::FLambdaMasks> producedLambdaMasks;
  o2::framework::Produces<o2::aod::FLambdaExtras> producedLambdaExtras;
  o2::framework::Produces<o2::aod::FK0shorts> producedK0shorts;
  o2::framework::Produces<o2::aod::FK0shortMasks> producedK0shortMasks;
  o2::framework::Produces<o2::aod::FK0shortExtras> producedK0shortExtras;
};

struct ConfV0Tables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("V0Tables");
  o2::framework::Configurable<int> produceLambdas{"produceLambdas", -1, "Produce Lambdas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLambdaMasks{"produceLambdaMasks", -1, "Produce LambdaMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLambdaExtras{"produceLambdaExtras", -1, "Produce LambdaExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceK0shorts{"produceK0shorts", -1, "Produce K0shorts (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceK0shortMasks{"produceK0shortMasks", -1, "Produce K0shortMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceK0shortExtras{"produceK0shortExtras", -1, "Produce K0shortExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::V0 v0Type, const char* HistName>
class V0Builder
{
 public:
  V0Builder() = default;
  ~V0Builder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext)
  {
    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda) || modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        LOG(info) << "Initialize femto Lambda builder...";
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        LOG(info) << "Initialize femto AntiLambda builder...";
      }
      mProduceLambdas = utils::enableTable("FLambdas_001", table.produceLambdas.value, initContext);
      mProduceLambdaMasks = utils::enableTable("FLambdaMasks_001", table.produceLambdaMasks.value, initContext);
      mProduceLambdaExtras = utils::enableTable("FLambdaExtras_001", table.produceLambdaExtras.value, initContext);
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      LOG(info) << "Initialize femto K0short builder...";
      mProduceK0shorts = utils::enableTable("FK0shorts_001", table.produceK0shorts.value, initContext);
      mProduceK0shortMasks = utils::enableTable("FK0shortMasks_001", table.produceK0shortMasks.value, initContext);
      mProduceK0shortExtras = utils::enableTable("FK0shortExtras_001", table.produceK0shortExtras.value, initContext);
    }
    if (mProduceLambdas || mProduceLambdaMasks || mProduceLambdaExtras || mProduceK0shorts || mProduceK0shortMasks || mProduceK0shortExtras) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mV0Selection.configure(registry, config, filter);
    mV0Selection.printSelections(V0SelsName);
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void fillV0s(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts, T5& v0Products, T6 const& v0s, T7 const& tracks, T8 const& tracksWithItsPid, T9& trackBuilder)
  {
    if (!mFillAnyTable) {
      return;
    }
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& v0 : v0s) {
      if (!mV0Selection.checkFilters(v0)) {
        continue;
      }
      mV0Selection.applySelections(v0, tracks);
      if (!mV0Selection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillCollision<system>(collisionProducts, col);

      // cleaner, but without ITS pid: auto posDaughter = v0.template posTrack_as<T7>();
      auto posDaughter = tracksWithItsPid.iteratorAt(v0.posTrackId() - tracksWithItsPid.offset());
      // cleaner, but without ITS pid: auto negDaughter = v0.template negTrack_as<T7>();
      auto negDaughter = tracksWithItsPid.iteratorAt(v0.negTrackId() - tracksWithItsPid.offset());

      posDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(posDaughter, trackProducts, collisionProducts);
      negDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(negDaughter, trackProducts, collisionProducts);

      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        fillLambda(collisionProducts, v0Products, v0, 1.f, posDaughterIndex, negDaughterIndex);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        fillLambda(collisionProducts, v0Products, v0, -1.f, posDaughterIndex, negDaughterIndex);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
        fillK0short(collisionProducts, v0Products, v0, posDaughterIndex, negDaughterIndex);
      }
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13>
  void fillMcV0s(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts, T6& v0Products, T7 const& v0s, T8 const& tracks, T9 const& tracksWithItsPid, T10& trackBuilder, T11 const& mcParticles, T12& mcBuilder, T13& mcProducts)
  {

    if (!mFillAnyTable) {
      return;
    }
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& v0 : v0s) {
      if (!mV0Selection.checkFilters(v0)) {
        continue;
      }
      mV0Selection.applySelections(v0, tracks);
      if (!mV0Selection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);

      auto posDaughter = tracks.iteratorAt(v0.posTrackId() - tracksWithItsPid.offset());
      auto posDaughterWithItsPid = tracksWithItsPid.iteratorAt(v0.posTrackId() - tracksWithItsPid.offset());
      posDaughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kV0Daughter>(col, collisionProducts, mcCols, posDaughter, posDaughterWithItsPid, trackProducts, mcParticles, mcBuilder, mcProducts);

      auto negDaughter = tracks.iteratorAt(v0.negTrackId() - tracksWithItsPid.offset());
      auto negDaughterWithItsPid = tracksWithItsPid.iteratorAt(v0.negTrackId() - tracksWithItsPid.offset());
      negDaughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kV0Daughter>(col, collisionProducts, mcCols, negDaughter, negDaughterWithItsPid, trackProducts, mcParticles, mcBuilder, mcProducts);

      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        fillLambda(collisionProducts, v0Products, v0, 1.f, posDaughterIndex, negDaughterIndex);
        mcBuilder.template fillMcLambdaWithLabel<system>(col, mcCols, v0, mcParticles, mcProducts);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        fillLambda(collisionProducts, v0Products, v0, -1.f, posDaughterIndex, negDaughterIndex);
        mcBuilder.template fillMcLambdaWithLabel<system>(col, mcCols, v0, mcParticles, mcProducts);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
        fillK0short(collisionProducts, v0Products, v0, posDaughterIndex, negDaughterIndex);
        mcBuilder.template fillMcK0shortWithLabel<system>(col, mcCols, v0, mcParticles, mcProducts);
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillLambda(T1& collisionProducts, T2& v0Products, T3 const& v0, float sign, int64_t posDaughterIndex, int64_t negDaughterIndex)
  {
    float mass, massAnti;
    if (sign > 0.f) {
      mass = v0.mLambda();
      massAnti = v0.mAntiLambda();
    } else {
      mass = v0.mAntiLambda();
      massAnti = v0.mLambda();
    }
    if (mProduceLambdas) {
      v0Products.producedLambdas(collisionProducts.producedCollision.lastIndex(),
                                 sign * v0.pt(),
                                 v0.eta(),
                                 v0.phi(),
                                 mass,
                                 posDaughterIndex,
                                 negDaughterIndex);
    }
    if (mProduceLambdaMasks) {
      v0Products.producedLambdaMasks(mV0Selection.getBitmask());
    }
    if (mProduceLambdaExtras) {
      v0Products.producedLambdaExtras(
        massAnti,
        v0.mK0Short(),
        v0.v0cosPA(),
        v0.dcaV0daughters(),
        v0.v0radius(),
        v0.x(),
        v0.y(),
        v0.z());
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillK0short(T1& collisionProducts, T2& v0Products, T3 const& v0, int64_t posDaughterIndex, int64_t negDaughterIndex)
  {
    if (mProduceK0shorts) {
      v0Products.producedK0shorts(collisionProducts.producedCollision.lastIndex(),
                                  v0.pt(),
                                  v0.eta(),
                                  v0.phi(),
                                  v0.mK0Short(),
                                  posDaughterIndex,
                                  negDaughterIndex);
    }
    if (mProduceK0shortMasks) {
      v0Products.producedK0shortMasks(mV0Selection.getBitmask());
    }
    if (mProduceK0shortExtras) {
      v0Products.producedK0shortExtras(
        v0.mLambda(),
        v0.mAntiLambda(),
        v0.v0cosPA(),
        v0.dcaV0daughters(),
        v0.v0radius(),
        v0.x(),
        v0.y(),
        v0.z());
    }
  }

  bool fillAnyTable() const { return mFillAnyTable; }

 private:
  V0Selection<v0Type, HistName> mV0Selection;
  bool mFillAnyTable = false;
  bool mProduceLambdas = false;
  bool mProduceLambdaMasks = false;
  bool mProduceLambdaExtras = false;
  bool mProduceK0shorts = false;
  bool mProduceK0shortMasks = false;
  bool mProduceK0shortExtras = false;
};

struct ConfV0TablesDerivedToDerived : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("V0Tables");
  o2::framework::Configurable<int> limitLambda{"limitLambda", 1, "At least this many lambdas need to be in the collision"};
  o2::framework::Configurable<int> limitK0short{"limitK0short", 0, "At least this many k0short need to be in the collision"};
};

struct V0BuilderDerivedToDerivedProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::StoredFLambdas> producedLambdas;
  o2::framework::Produces<o2::aod::StoredFLambdaMasks> producedLambdaMasks;
  o2::framework::Produces<o2::aod::StoredFK0shorts> producedK0shorts;
  o2::framework::Produces<o2::aod::StoredFK0shortMasks> producedK0shortMasks;
};

class V0BuilderDerivedToDerived
{
 public:
  V0BuilderDerivedToDerived() = default;
  ~V0BuilderDerivedToDerived() = default;

  template <typename T>
  void init(T& config)
  {
    mLimitLambda = config.limitLambda.value;
    mLimitK0short = config.limitK0short.value;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewLambdas(T1& col, T2& /*lambdaTable*/, T3& partitionLambda, T4& cache)
  {
    auto lambdaSlice = partitionLambda->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (lambdaSlice.size() >= mLimitLambda) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewK0shorts(T1& col, T2& /*k0shortTable*/, T3& partitionK0short, T4& cache)
  {
    auto k0shortSlice = partitionK0short->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (k0shortSlice.size() >= mLimitK0short) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processLambdas(T1& col, T2& /*lambdaTable*/, T3& /*oldTrackTable*/, T4& partitionLambda, T5& trackBuilder, T6& cache, T7& newLambdaTable, T8& newTrackTable, T9& newCollisionTable)
  {
    auto lambdaSlice = partitionLambda->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& lambda : lambdaSlice) {

      auto posDaughter = lambda.template posDau_as<T3>();
      auto negDaughter = lambda.template negDau_as<T3>();

      int posDaughterIndex = trackBuilder.getDaughterIndex(posDaughter, newTrackTable, newCollisionTable);
      int negDaughterIndex = trackBuilder.getDaughterIndex(negDaughter, newTrackTable, newCollisionTable);

      newLambdaTable.producedLambdas(newCollisionTable.producedCollision.lastIndex(),
                                     lambda.signedPt(),
                                     lambda.eta(),
                                     lambda.phi(),
                                     lambda.mass(),
                                     posDaughterIndex,
                                     negDaughterIndex);
      newLambdaTable.producedLambdaMasks(lambda.mask());
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processK0shorts(T1& col, T2& /*k0shortTable*/, T3& /*oldTrackTable*/, T4& partitionK0short, T5& trackBuilder, T6& cache, T7& newK0shortTable, T8& newTrackTable, T9& newCollisionTable)
  {
    auto k0shortSlice = partitionK0short->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& k0short : k0shortSlice) {

      auto posDaughter = k0short.template posDau_as<T3>();
      auto negDaughter = k0short.template negDau_as<T3>();

      int posDaughterIndex = trackBuilder.getDaughterIndex(posDaughter, newTrackTable, newCollisionTable);
      int negDaughterIndex = trackBuilder.getDaughterIndex(negDaughter, newTrackTable, newCollisionTable);

      newK0shortTable.producedK0shorts(newCollisionTable.producedCollision.lastIndex(),
                                       k0short.pt(),
                                       k0short.eta(),
                                       k0short.phi(),
                                       k0short.mass(),
                                       posDaughterIndex,
                                       negDaughterIndex);
      newK0shortTable.producedK0shortMasks(k0short.mask());
    }
  }

 private:
  int mLimitLambda = 0;
  int mLimitK0short = 0;
};

} // namespace v0builder
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_V0BUILDER_H_
