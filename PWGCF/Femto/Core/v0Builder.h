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
#define V0_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                                 \
  o2::framework::Configurable<int> pdgCode{"pdgCode", defaultPdgCode, "V0 PDG code"};                         \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                       \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                     \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                  \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                   \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};                                    \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};       \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Lambda"}; \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Lambda"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::V0MaskType> mask{"mask", 0, "Bitmask for v0 selection"};

// base selection for analysis task for lambdas
template <const char* Prefix>
struct ConfLambdaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_SELECTIONS(1.0, 1.2, 3122)
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Lambda (+1 for Lambda and -1 for Antilambda"};
};

// base selection for analysis task for k0short
template <const char* Prefix>
struct ConfK0shortSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_SELECTIONS(0.47, 0.51, 310)
};

#undef V0_DEFAULT_SELECTIONS

constexpr const char PrefixLambdaSelection1[] = "LambdaSelection1";
using ConfLambdaSelection1 = ConfLambdaSelection<PrefixLambdaSelection1>;
constexpr const char PrefixK0shortSelection1[] = "K0shortSelection1";
using ConfK0shortSelection1 = ConfK0shortSelection<PrefixK0shortSelection1>;

/// The different selections for v0s
enum V0Seles {
  // selections for lambdas
  kCpaMin,      ///< Min. CPA (cosine pointing angle)
  kDcaDaughMax, ///< Max. DCA of the daughters at decay vertex
  kDecayVtxMax, ///< Max. distance of decay vertex in x,y,z
  kTransRadMin, ///< Min. transverse radius
  kTransRadMax, ///< max. transverse radius

  // selection for daughter
  kDauAbsEtaMax, ///< Max. absolute pseudo rapidity
  kDauDcaMin,    ///< Min. DCA of the positive daughters at primary vertex
  kDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter

  // pid selection for daughters
  kPosDaughTpcPion,   ///< TPC Pion PID for positive daughter
  kPosDaughTpcProton, ///< TPC Proton PID for positive daughter
  kNegDaughTpcPion,   ///< TPC Pion PID for negative daughter
  kNegDaughTpcProton, ///< TPC Proton PID for negative daughter

  kV0SelsMax
};

const char v0SelsName[] = "K0short selection object";
const std::unordered_map<V0Seles, std::string> v0SelsToStrings = {
  {kCpaMin, "Min. CPA (cosine pointing angle)"},
  {kDcaDaughMax, "Max. DCA of the daughters at decay vertex"},
  {kDecayVtxMax, "Max. distance of decay vertex in x,y,z"},
  {kTransRadMin, "Min. transverse radius"},
  {kTransRadMax, "Max. transverse radius"},

  {kDauAbsEtaMax, "Max. absolute pseudo rapidity"},
  {kDauDcaMin, "Min. DCA of the positive daughters at primary vertex"},
  {kDauTpcClsMin, "Min. number of TPC clusters of positive daughter"},

  {kPosDaughTpcPion, "TPC Pion PID for positive daughter"},
  {kPosDaughTpcProton, "TPC Proton PID for positive daughter"},
  {kNegDaughTpcPion, "TPC Pion PID for negative daughter"},
  {kNegDaughTpcProton, "TPC Proton PID for negative daughter"}};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <modes::V0 v0Type>
class V0Selection : public BaseSelection<float, o2::aod::femtodatatypes::V0MaskType, kV0SelsMax>
{
 public:
  V0Selection() {}
  virtual ~V0Selection() = default;

  template <typename T1, typename T2>
  void configure(T1& config, T2& filter)
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
        this->addSelection(config.posDauTpcProton.value, kPosDaughTpcProton, limits::kAbsUpperLimit, true, true);
        this->addSelection(config.negDauTpcPion.value, kNegDaughTpcPion, limits::kAbsUpperLimit, true, true);
      }

      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        this->addSelection(config.posDauTpcPion.value, kPosDaughTpcPion, limits::kAbsUpperLimit, true, true);
        this->addSelection(config.negDauTpcProton.value, kNegDaughTpcProton, limits::kAbsUpperLimit, true, true);
      }
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      mMassK0shortLowerLimit = filter.massMinK0short.value;
      mMassK0shortUpperLimit = filter.massMaxK0short.value;
      mMassLambdaLowerLimit = filter.rejectMassMinLambda.value;
      mMassLambdaUpperLimit = filter.rejectMassMaxLambda.value;
      this->addSelection(config.posDauTpcPion.value, kPosDaughTpcPion, limits::kAbsUpperLimit, true, true);
      this->addSelection(config.negDauTpcPion.value, kNegDaughTpcPion, limits::kAbsUpperLimit, true, true);
    }

    this->addSelection(config.dcaDauMax.value, kDcaDaughMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.cpaMin.value, kCpaMin, limits::kLowerLimit, true, true);
    this->addSelection(config.transRadMin.value, kTransRadMin, limits::kLowerLimit, true, true);
    this->addSelection(config.transRadMax.value, kTransRadMax, limits::kUpperLimit, true, true);
    this->addSelection(config.dauAbsEtaMax.value, kDauAbsEtaMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauDcaMin.value, kDauDcaMin, limits::kAbsLowerFunctionLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, kDauTpcClsMin, limits::kLowerLimit, true, true);
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

    this->assembleBitmask();
  }

  template <typename T>
  bool checkFilters(const T& v0) const
  {
    return ((v0.pt() > mPtMin && v0.pt() < mPtMax) &&
            (v0.eta() > mEtaMin && v0.eta() < mEtaMax) &&
            (v0.phi() > mPhiMin && v0.phi() < mPhiMax));
  }

  template <typename T>
  bool checkHypothesis(T const& v0candidate) const
  {
    // no need to check PID of the daughters here, they are set as minimal cuts
    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
      return (v0candidate.mLambda() > mMassLambdaLowerLimit && v0candidate.mLambda() < mMassLambdaUpperLimit) &&   // inside Lambda window
             (v0candidate.mK0Short() < mMassK0shortLowerLimit || v0candidate.mK0Short() > mMassK0shortUpperLimit); // outside K0short window
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      return                                                                                                        // check PID for daughters
        (v0candidate.mAntiLambda() > mMassLambdaLowerLimit && v0candidate.mAntiLambda() < mMassLambdaUpperLimit) && // inside AntiLambda window
        (v0candidate.mK0Short() < mMassK0shortLowerLimit || v0candidate.mK0Short() > mMassK0shortUpperLimit);       // outside K0short window
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      return (v0candidate.mK0Short() > mMassK0shortLowerLimit && v0candidate.mK0Short() < mMassK0shortUpperLimit) &&   // inside K0short window
             (v0candidate.mLambda() < mMassLambdaLowerLimit || v0candidate.mLambda() > mMassLambdaUpperLimit) &&       // outside Lambda window
             (v0candidate.mAntiLambda() < mMassLambdaLowerLimit || v0candidate.mAntiLambda() > mMassLambdaUpperLimit); // outside AntiLambda window
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

template <modes::V0 v0Type>
class V0Builder
{
 public:
  V0Builder() {}
  virtual ~V0Builder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(T1& config, T2& filter, T3& table, T4& initContext)
  {
    v0Selection.configure(config, filter);
    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda) || modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        LOG(info) << "Initialize femto Lambda builder...";
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        LOG(info) << "Initialize femto AntiLambda builder...";
      }
      produceLambdas = utils::enableTable("FLambdas_001", table.produceLambdas.value, initContext);
      produceLambdaMasks = utils::enableTable("FLambdaMasks_001", table.produceLambdaMasks.value, initContext);
      produceLambdaExtras = utils::enableTable("FLambdaExtras_001", table.produceLambdaExtras.value, initContext);
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      LOG(info) << "Initialize femto K0short builder...";
      produceK0shorts = utils::enableTable("FK0shorts_001", table.produceK0shorts.value, initContext);
      produceK0shortMasks = utils::enableTable("FK0shortMasks_001", table.produceK0shortMasks.value, initContext);
      produceK0shortExtras = utils::enableTable("FK0shortExtras_001", table.produceK0shortExtras.value, initContext);
    }
    if (produceLambdas || produceLambdaMasks || produceLambdaExtras || produceK0shorts || produceK0shortMasks || produceK0shortExtras) {
      mFillAnyTable = true;
      v0Selection.printSelections(v0SelsName, v0SelsToStrings);
    } else {
      LOG(info) << "No tables configured";
    }
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void fillV0s(T1& collisionProducts, T2& trackProducts, T3& v0products, T4 const& v0s, T5 const& tracks, T6& trackBuilder, T7& indexMap)
  {
    if (!mFillAnyTable) {
      return;
    }
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& v0 : v0s) {
      if (!v0Selection.checkFilters(v0)) {
        continue;
      }
      v0Selection.applySelections(v0, tracks);
      if (v0Selection.passesAllRequiredSelections() && v0Selection.checkHypothesis(v0)) {
        auto posDaughter = v0.template posTrack_as<T5>();
        auto negDaughter = v0.template negTrack_as<T5>();
        posDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(posDaughter, trackProducts, collisionProducts, indexMap);
        negDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(negDaughter, trackProducts, collisionProducts, indexMap);
        if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
          fillLambda(collisionProducts, v0products, v0, 1.f, posDaughterIndex, negDaughterIndex);
        }
        if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
          fillLambda(collisionProducts, v0products, v0, -1.f, posDaughterIndex, negDaughterIndex);
        }
        if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
          fillK0short(collisionProducts, v0products, v0, posDaughterIndex, negDaughterIndex);
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillLambda(T1& collisionProducts, T2& v0products, T3 const& v0, float sign, int posDaughterIndex, int negDaughterIndex)
  {
    float mass, massAnti;
    if (sign > 0.f) {
      mass = v0.mLambda();
      massAnti = v0.mAntiLambda();
    } else {
      mass = v0.mAntiLambda();
      massAnti = v0.mLambda();
    }
    if (produceLambdas) {
      v0products.producedLambdas(collisionProducts.producedCollision.lastIndex(),
                                 sign * v0.pt(),
                                 v0.eta(),
                                 v0.phi(),
                                 mass,
                                 posDaughterIndex,
                                 negDaughterIndex);
    }
    if (produceLambdaMasks) {
      v0products.producedLambdaMasks(v0Selection.getBitmask());
    }
    if (produceLambdaExtras) {
      v0products.producedLambdaExtras(
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
  void fillK0short(T1& collisionProducts, T2& v0products, T3 const& v0, int posDaughterIndex, int negDaughterIndex)
  {
    if (produceK0shorts) {
      v0products.producedK0shorts(collisionProducts.producedCollision.lastIndex(),
                                  v0.pt(),
                                  v0.eta(),
                                  v0.phi(),
                                  v0.mK0Short(),
                                  posDaughterIndex,
                                  negDaughterIndex);
    }
    if (produceK0shortMasks) {
      v0products.producedK0shortMasks(v0Selection.getBitmask());
    }
    if (produceK0shortExtras) {
      v0products.producedK0shortExtras(
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

  bool fillAnyTable() { return mFillAnyTable; }

 private:
  V0Selection<v0Type> v0Selection;
  bool mFillAnyTable = false;
  bool produceLambdas = false;
  bool produceLambdaMasks = false;
  bool produceLambdaExtras = false;
  bool produceK0shorts = false;
  bool produceK0shortMasks = false;
  bool produceK0shortExtras = false;
};
} // namespace v0builder
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_V0BUILDER_H_
