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

/// \file v0Selection.h
/// \brief Definition of v0 selections
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_V0SELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_V0SELECTION_H_

#include "PWGCF/FemtoUnited/Core/baseSelection.h"
#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/Core/modes.h"

#include "Framework/Configurable.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace v0selection
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
  o2::framework::Configurable<float> rejectMassMinLambda{"rejectMassMinLambda", 0.111f, "Minimum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMaxLambda{"rejectMassMaxLambda", 0.112f, "Maximum mass to rejection K0short hypothesis for Lambda candidates"};
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

#include <string>
#include <unordered_map>

const std::string lambdaSelsName = std::string("Lambda selection object");
const std::string antiLambdaSelsName = std::string("AntiLambda selection object");
const std::string k0shortSelsName = std::string("K0short selection object");
const std::unordered_map<V0Seles, std::string> V0SelesNames = {
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
template <o2::analysis::femtounited::modes::V0 v0Type>
class V0Selection : public BaseSelection<float, o2::aod::femtodatatypes::V0MaskType, kV0SelsMax>
{
 public:
  V0Selection() {}
  virtual ~V0Selection() = default;

  template <typename T1, typename T2>
  void configure(T1 const& config, T2 const& filter)
  {
    if constexpr (o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kLambda) ||
                  o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kAntiLambda)) {
      mMassLambdaLowerLimit = filter.massMinLambda.value;
      mMassLambdaUpperLimit = filter.massMaxLambda.value;
      mMassK0shortLowerLimit = filter.rejectMassMinK0short.value;
      mMassK0shortUpperLimit = filter.rejectMassMaxK0short.value;

      if constexpr (o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kLambda)) {
        this->addSelection(config.posDauTpcProton.value, kPosDaughTpcProton, limits::kAbsUpperLimit, true, true);
        this->addSelection(config.negDauTpcPion.value, kNegDaughTpcPion, limits::kAbsUpperLimit, true, true);
      }

      if constexpr (o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kAntiLambda)) {
        this->addSelection(config.posDauTpcPion.value, kPosDaughTpcPion, limits::kAbsUpperLimit, true, true);
        this->addSelection(config.negDauTpcProton.value, kNegDaughTpcProton, limits::kAbsUpperLimit, true, true);
      }
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kK0short)) {
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
  bool checkHypothesis(T const& v0candidate) const
  {
    // no need to check PID of the daughters here, they are set as minimal cuts
    if constexpr (o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kLambda)) {
      return (v0candidate.mLambda() > mMassLambdaLowerLimit && v0candidate.mLambda() < mMassLambdaUpperLimit) &&   // inside Lambda window
             (v0candidate.mK0Short() < mMassK0shortLowerLimit || v0candidate.mK0Short() > mMassK0shortUpperLimit); // outside K0short window
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kAntiLambda)) {
      return                                                                                                        // check PID for daughters
        (v0candidate.mAntiLambda() > mMassLambdaLowerLimit && v0candidate.mAntiLambda() < mMassLambdaUpperLimit) && // inside AntiLambda window
        (v0candidate.mK0Short() < mMassK0shortLowerLimit || v0candidate.mK0Short() > mMassK0shortUpperLimit);       // outside K0short window
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(v0Type, o2::analysis::femtounited::modes::V0::kK0short)) {
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
}; // namespace v0selection
} // namespace v0selection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_V0SELECTION_H_
