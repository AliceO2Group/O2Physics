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

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto::v0builder
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
  o2::framework::Configurable<bool> rejectHypothesisK0short{"rejectHypothesisK0short", true, "Rejection of K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMinK0short{"rejectMassMinK0short", 0.475f, "Minimum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMaxK0short{"rejectMassMaxK0short", 0.515f, "Maximum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> massMinK0short{"massMinK0short", 0.45f, "Minimum mass for K0Short hypothesis"};
  o2::framework::Configurable<float> massMaxK0short{"massMaxK0short", 0.53f, "Maximum mass for K0Short hypothesis"};
  o2::framework::Configurable<bool> rejectHypothesisLambda{"rejectHypothesisLambda", true, "Rejection of Lambda hypothesis for K0short candidates"};
  o2::framework::Configurable<float> rejectMassMinLambda{"rejectMassMinLambda", 1.11f, "Minimum mass to rejection K0short hypothesis for Lambda candidates"};
  o2::framework::Configurable<float> rejectMassMaxLambda{"rejectMassMaxLambda", 1.12f, "Maximum mass to rejection K0short hypothesis for Lambda candidates"};
};

// selections bits for all v0s
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define V0_DEFAULT_BITS                                                                                                                                          \
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, all V0s are passed through. Bits for all selections are stored."};               \
  o2::framework::Configurable<std::vector<float>> dcaDauMax{"dcaDauMax", {1.5f}, "Maximum DCA between the daughters at V0 decay vertex (cm)"};                   \
  o2::framework::Configurable<std::vector<float>> cpaMin{"cpaMin", {0.99f}, "Minimum cosine of pointing angle"};                                                 \
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2f}, "Minimum transverse radius (cm)"};                                          \
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};                                         \
  o2::framework::Configurable<std::vector<float>> decayVtxMax{"decayVtxMax", {100.f}, "Maximum distance in x,y,z of the decay vertex from primary vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.8f}, "Maximum |eta| for daughter tracks"};                                     \
  o2::framework::Configurable<std::vector<float>> dauAbsDcaxyMin{"dauAbsDcaxyMin", {0.05f}, "Minimum DCAxy of the daughters from primary vertex (cm)"};          \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};          \
  o2::framework::Configurable<bool> requireTof{"requireTof", false, "If true, TOF PID is a minimal selection. If false, TOF PID is an optional"};                \
  o2::framework::Configurable<bool> keepTracksWithoutTof{"keepTracksWithoutTof", true, "If true, the bit mask for the TOF selection will be ture for all limits if the daugher track has no TOF"};

// derived selection bits for lambda
struct ConfLambdaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("LambdaBits");
  V0_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDauTpcProton{"posDauTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcProton{"negDauTpcProton", {5.f}, "Maximum |nsimga_Proton| TPC negative for daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDauTofPion{"posDauTofPion", {}, "Maximum |nsigma_Pion| TOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDauTofProton{"posDauTofProton", {}, "Maximum |nsigma_Proton| TOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTofPion{"negDauTofPion", {}, "Maximum |nsigma_Pion| TOF for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTofProton{"negDauTofProton", {}, "Maximum |nsigma_Proton| TOF for negative daughter tracks"};
};

// derived selection bits for K0Short
struct ConfK0shortBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("K0shortBits");
  V0_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDauTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  o2::framework::Configurable<std::vector<float>> posDauTofPion{"posDauTofPion", {}, "Maximum |nsigma_Pion| TOF for positive daughter tracks"};
  o2::framework::Configurable<std::vector<float>> negDauTofPion{"negDauTofPion", {}, "Maximum |nsigma_Pion| TOF for negative daughter tracks"};
};

#undef V0_DEFAULT_BITS

// base selection for analysis task for v0s
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define V0_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                                               \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", (defaultPdgCode), "PDG code. Set sign to -1 for antiparticle"}; \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                     \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                                   \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                                \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                                 \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                  \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                     \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass"};                        \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass"};                        \
  o2::framework::Configurable<datatypes::V0MaskType> mask{"mask", 0, "Bitmask for v0 selection"};

// base selection for analysis task for lambdas
template <auto& Prefix>
struct ConfLambdaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  V0_DEFAULT_SELECTIONS(1.0, 1.2, 3122)
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Lambda (+1: Lambda; -1: Antilambda; 0: both)"};
};

// base selection for analysis task for k0short
template <auto& Prefix>
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
  kDauAbsEtaMax,   ///< Max. absolute pseudo rapidity
  kDauAbsDcaxyMin, ///< Min. |DCAxy| of the daughters from primary vertex
  kDauTpcClsMin,   ///< Min. number of TPC clusters of daughter

  // pid selection for daughters
  kPosDaughTpcPion,   ///< TPC Pion PID for positive daughter
  kPosDaughTpcProton, ///< TPC Proton PID for positive daughter
  kNegDaughTpcPion,   ///< TPC Pion PID for negative daughter
  kNegDaughTpcProton, ///< TPC Proton PID for negative daughter

  kPosDaughTofPion,   ///< TOF Pion PID for positive daughter
  kPosDaughTofProton, ///< TOF Proton PID for positive daughter
  kNegDaughTofPion,   ///< TOF Pion PID for negative daughter
  kNegDaughTofProton, ///< TOF Proton PID for negative daughter

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
  {kDauAbsDcaxyMin, "Min. |DCAxy| of the daughters from primary vertex"},
  {kDauTpcClsMin, "Min. number of TPC clusters of daughters"},

  {kPosDaughTpcPion, "TPC Pion PID for positive daughter"},
  {kPosDaughTpcProton, "TPC Proton PID for positive daughter"},
  {kNegDaughTpcPion, "TPC Pion PID for negative daughter"},
  {kNegDaughTpcProton, "TPC Proton PID for negative daughter"},
  {kPosDaughTofPion, "TOF Pion PID for positive daughter"},
  {kPosDaughTofProton, "TOF Proton PID for positive daughter"},
  {kNegDaughTofPion, "TOF Pion PID for negative daughter"},
  {kNegDaughTofProton, "TOF Proton PID for negative daughter"}};

// enum for all track filters (loose pre-selection, applied before quality/PID cuts)
enum V0Filters {
  kPtMin,
  kPtMax,
  kEtaMin,
  kEtaMax,
  kPhiMin,
  kPhiMax,
  kRejectionK0shortMass,
  kK0shortMassMin,
  kK0shortMassMax,
  kRejectionLambdaMass,
  kLambdaMassMin,
  kLambdaMassMax,
  kV0FiltersMax
};

constexpr char LambdaFilterHistName[] = "hLambdaFilters";
constexpr char AntiLambdaFilterHistName[] = "hAntiLambdaFilters";
constexpr char K0shortFilterHistName[] = "hK0shortFilters";
const std::unordered_map<V0Filters, std::string> v0FilterNames = {
  {kPtMin, "ptMin"},
  {kPtMax, "ptMax"},
  {kEtaMin, "etaMin"},
  {kEtaMax, "etaMax"},
  {kPhiMin, "phiMin"},
  {kPhiMax, "phiMax"},
  {kRejectionK0shortMass, "rejectK0short"},
  {kK0shortMassMin, "k0shortMassMin"},
  {kK0shortMassMax, "k0shortMassMax"},
  {kRejectionLambdaMass, "rejectLambda"},
  {kLambdaMassMin, "lambdaMassMin"},
  {kLambdaMassMax, "lambdaMassMax"},
};

/// \brief Cut class to contain and execute all cuts applied to v0s
template <modes::V0 v0Type, auto& SelectionHistName, auto& FilterHistName>
class V0Selection : public baseselection::BaseSelection<float, datatypes::V0MaskType, kV0SelsMax>
{
 public:
  V0Selection() = default;
  ~V0Selection() override = default;

  template <typename T1, typename T2>
  void configure(o2::framework::HistogramRegistry* registry, T1& config, T2& filter)
  {
    this->init(config.passThrough.value);

    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;
    mRequireTof = config.requireTof.value;
    mKeepTracksWithoutTof = config.keepTracksWithoutTof.value;

    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda) || modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      mMassLambdaLowerLimit = filter.massMinLambda.value;
      mMassLambdaUpperLimit = filter.massMaxLambda.value;
      mRejectK0shortHypothesis = filter.rejectHypothesisK0short.value;
      mMassK0shortLowerLimit = filter.rejectMassMinK0short.value;
      mMassK0shortUpperLimit = filter.rejectMassMaxK0short.value;

      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        this->addSelection(kPosDaughTpcProton, v0SelectionNames.at(kPosDaughTpcProton), config.posDauTpcProton.value, limits::kAbsUpperLimit, true, true, false);
        this->addSelection(kNegDaughTpcPion, v0SelectionNames.at(kNegDaughTpcPion), config.negDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
        this->addSelection(kPosDaughTofProton, v0SelectionNames.at(kPosDaughTofProton), config.posDauTofProton.value, limits::kAbsUpperLimit, true, mRequireTof, false);
        this->addSelection(kNegDaughTofPion, v0SelectionNames.at(kNegDaughTofPion), config.negDauTofPion.value, limits::kAbsUpperLimit, true, mRequireTof, false);
      }

      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        this->addSelection(kPosDaughTpcPion, v0SelectionNames.at(kPosDaughTpcPion), config.posDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
        this->addSelection(kNegDaughTpcProton, v0SelectionNames.at(kNegDaughTpcProton), config.negDauTpcProton.value, limits::kAbsUpperLimit, true, true, false);
        this->addSelection(kPosDaughTofPion, v0SelectionNames.at(kPosDaughTofPion), config.posDauTofPion.value, limits::kAbsUpperLimit, true, mRequireTof, false);
        this->addSelection(kNegDaughTofProton, v0SelectionNames.at(kNegDaughTofProton), config.negDauTofProton.value, limits::kAbsUpperLimit, true, mRequireTof, false);
      }
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      mMassK0shortLowerLimit = filter.massMinK0short.value;
      mMassK0shortUpperLimit = filter.massMaxK0short.value;
      mRejectLambdaHypothesis = filter.rejectHypothesisLambda.value;
      mMassLambdaLowerLimit = filter.rejectMassMinLambda.value;
      mMassLambdaUpperLimit = filter.rejectMassMaxLambda.value;

      this->addSelection(kPosDaughTpcPion, v0SelectionNames.at(kPosDaughTpcPion), config.posDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
      this->addSelection(kNegDaughTpcPion, v0SelectionNames.at(kNegDaughTpcPion), config.negDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
      this->addSelection(kPosDaughTofPion, v0SelectionNames.at(kPosDaughTofPion), config.posDauTofPion.value, limits::kAbsUpperLimit, true, mRequireTof, false);
      this->addSelection(kNegDaughTofPion, v0SelectionNames.at(kNegDaughTofPion), config.negDauTofPion.value, limits::kAbsUpperLimit, true, mRequireTof, false);
    }

    this->addSelection(kDcaDaughMax, v0SelectionNames.at(kDcaDaughMax), config.dcaDauMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kDecayVtxMax, v0SelectionNames.at(kDecayVtxMax), config.decayVtxMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kCpaMin, v0SelectionNames.at(kCpaMin), config.cpaMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTransRadMin, v0SelectionNames.at(kTransRadMin), config.transRadMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTransRadMax, v0SelectionNames.at(kTransRadMax), config.transRadMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kDauAbsEtaMax, v0SelectionNames.at(kDauAbsEtaMax), config.dauAbsEtaMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kDauAbsDcaxyMin, v0SelectionNames.at(kDauAbsDcaxyMin), config.dauAbsDcaxyMin.value, limits::kAbsLowerLimit, true, true, false);
    this->addSelection(kDauTpcClsMin, v0SelectionNames.at(kDauTpcClsMin), config.dauTpcClustersMin.value, limits::kLowerLimit, true, true, false);

    this->setupSelectionHistogram<SelectionHistName>(registry);
    this->template setupFilterHistogram<FilterHistName>(
      registry,
      {
        {v0FilterNames.at(kPtMin), mPtMin},
        {v0FilterNames.at(kPtMax), mPtMax},
        {v0FilterNames.at(kEtaMin), mEtaMin},
        {v0FilterNames.at(kEtaMax), mEtaMax},
        {v0FilterNames.at(kPhiMin), mPhiMin},
        {v0FilterNames.at(kPhiMax), mPhiMax},
        {v0FilterNames.at(kRejectionK0shortMass), mRejectK0shortHypothesis ? 1 : 0},
        {v0FilterNames.at(kK0shortMassMin), mMassK0shortLowerLimit},
        {v0FilterNames.at(kK0shortMassMax), mMassK0shortUpperLimit},
        {v0FilterNames.at(kRejectionLambdaMass), mRejectLambdaHypothesis ? 1 : 0},
        {v0FilterNames.at(kLambdaMassMin), mMassLambdaLowerLimit},
        {v0FilterNames.at(kLambdaMassMax), mMassLambdaUpperLimit},
      });
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

    std::array<float, 2> etaAbsDaughters = {std::fabs(v0candidate.positiveeta()), std::fabs(v0candidate.negativeeta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaAbsDaughters.begin(), etaAbsDaughters.end()));

    std::array<float, 2> dcaxyAbsDaughters = {std::fabs(v0candidate.dcapostopv()), std::fabs(v0candidate.dcanegtopv())};
    this->evaluateObservable(kDauAbsDcaxyMin, *std::min_element(dcaxyAbsDaughters.begin(), dcaxyAbsDaughters.end()));

    std::array<float, 2> clustersDaughters = {1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauTpcClsMin, *std::min_element(clustersDaughters.begin(), clustersDaughters.end()));

    // PID of the daughters under the mass hypothesis of this v0 type
    // TPC nSigma comes from the daughter track, TOF nSigma from the strangeness-tagged v0 candidate
    // if the daughter has no TOF signal, feed 0 so the bit passes any limit (opt-in via keepTracksWithoutTof)
    auto evaluateDaughterPid = [this](V0Sels tpcBit, float tpcNSigma,
                                      V0Sels tofBit, float tofNSigma, bool hasTof) {
      this->evaluateObservable(tpcBit, tpcNSigma);
      if (hasTof) {
        this->evaluateObservable(tofBit, tofNSigma);
      } else if (mKeepTracksWithoutTof) {
        this->evaluateObservable(tofBit, 0.f);
      }
    };

    const bool posHasTof = v0candidate.positiveHasTOF();
    const bool negHasTof = v0candidate.negativeHasTOF();

    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
      // Lambda -> p pi-
      evaluateDaughterPid(kPosDaughTpcProton, posDaughter.tpcNSigmaPr(),
                          kPosDaughTofProton, v0candidate.tofNSigmaLaPr(), posHasTof);
      evaluateDaughterPid(kNegDaughTpcPion, negDaughter.tpcNSigmaPi(),
                          kNegDaughTofPion, v0candidate.tofNSigmaLaPi(), negHasTof);
    } else if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      // AntiLambda -> pbar pi+
      evaluateDaughterPid(kPosDaughTpcPion, posDaughter.tpcNSigmaPi(),
                          kPosDaughTofPion, v0candidate.tofNSigmaALaPi(), posHasTof);
      evaluateDaughterPid(kNegDaughTpcProton, negDaughter.tpcNSigmaPr(),
                          kNegDaughTofProton, v0candidate.tofNSigmaALaPr(), negHasTof);
    } else if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      // K0short -> pi+ pi-
      evaluateDaughterPid(kPosDaughTpcPion, posDaughter.tpcNSigmaPi(),
                          kPosDaughTofPion, v0candidate.tofNSigmaK0PiPlus(), posHasTof);
      evaluateDaughterPid(kNegDaughTpcPion, negDaughter.tpcNSigmaPi(),
                          kNegDaughTofPion, v0candidate.tofNSigmaK0PiMinus(), negHasTof);
    }

    this->assembleBitmask<SelectionHistName>();
  }

  template <typename T>
  bool checkFilters(const T& v0) const
  {
    bool pass = true;
    bool p = false;

    // kinematics
    p = v0.pt() > mPtMin;
    this->template fillFilter<FilterHistName>(kPtMin, p);
    pass &= p;

    p = v0.pt() < mPtMax;
    this->template fillFilter<FilterHistName>(kPtMax, p);
    pass &= p;

    p = v0.eta() > mEtaMin;
    this->template fillFilter<FilterHistName>(kEtaMin, p);
    pass &= p;

    p = v0.eta() < mEtaMax;
    this->template fillFilter<FilterHistName>(kEtaMax, p);
    pass &= p;

    p = v0.phi() > mPhiMin;
    this->template fillFilter<FilterHistName>(kPhiMin, p);
    pass &= p;

    p = v0.phi() < mPhiMax;
    this->template fillFilter<FilterHistName>(kPhiMax, p);
    pass &= p;

    // mass hypothesis: signal window is a gating AND-cut (two bounds, both must pass);
    // the competing hypothesis is a rejection OR-cut (pass if outside the window, or if
    // rejection is disabled) — its own bin records whether the candidate cleared the
    // rejection, while the min/max bins are diagnostic bounds on where the mass sits.
    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda) || modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      float const mass = modes::isEqual(v0Type, modes::V0::kLambda) ? v0.mLambda() : v0.mAntiLambda();

      p = mass > mMassLambdaLowerLimit;
      this->template fillFilter<FilterHistName>(kLambdaMassMin, p);
      pass &= p;

      p = mass < mMassLambdaUpperLimit;
      this->template fillFilter<FilterHistName>(kLambdaMassMax, p);
      pass &= p;

      bool const belowK0s = v0.mK0Short() < mMassK0shortLowerLimit;
      this->template fillFilter<FilterHistName>(kK0shortMassMin, belowK0s);
      bool const aboveK0s = v0.mK0Short() > mMassK0shortUpperLimit;
      this->template fillFilter<FilterHistName>(kK0shortMassMax, aboveK0s);

      p = !mRejectK0shortHypothesis || belowK0s || aboveK0s;
      this->template fillFilter<FilterHistName>(kRejectionK0shortMass, p);
      pass &= p;
    }

    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      p = v0.mK0Short() > mMassK0shortLowerLimit;
      this->template fillFilter<FilterHistName>(kK0shortMassMin, p);
      pass &= p;

      p = v0.mK0Short() < mMassK0shortUpperLimit;
      this->template fillFilter<FilterHistName>(kK0shortMassMax, p);
      pass &= p;

      bool const belowLambda = v0.mLambda() < mMassLambdaLowerLimit;
      bool const aboveLambda = v0.mLambda() > mMassLambdaUpperLimit;
      bool const belowAntiLambda = v0.mAntiLambda() < mMassLambdaLowerLimit;
      bool const aboveAntiLambda = v0.mAntiLambda() > mMassLambdaUpperLimit;
      this->template fillFilter<FilterHistName>(kLambdaMassMin, belowLambda && belowAntiLambda);
      this->template fillFilter<FilterHistName>(kLambdaMassMax, aboveLambda && aboveAntiLambda);

      p = !mRejectLambdaHypothesis || ((belowLambda || aboveLambda) && (belowAntiLambda || aboveAntiLambda));
      this->template fillFilter<FilterHistName>(kRejectionLambdaMass, p);
      pass &= p;
    }

    this->template fillFilterSummary<FilterHistName>(pass);

    return this->isPassThrough() || pass;
  }

 protected:
  float mMassK0shortLowerLimit = 0.483f;
  float mMassK0shortUpperLimit = 0.503f;
  bool mRejectK0shortHypothesis = false;

  float mMassLambdaLowerLimit = 1.105f;
  float mMassLambdaUpperLimit = 1.125f;
  bool mRejectLambdaHypothesis = false;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -1.f;
  float mEtaMax = 1.f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;

  bool mRequireTof = false;
  bool mKeepTracksWithoutTof = false;
};

struct V0BuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FLambdas> producedLambdas;
  o2::framework::Produces<o2::aod::FLiteLambdas> producedLiteLambdas;
  o2::framework::Produces<o2::aod::FLambdaMasks> producedLambdaMasks;
  o2::framework::Produces<o2::aod::FLambdaExtras> producedLambdaExtras;
  o2::framework::Produces<o2::aod::FK0shorts> producedK0shorts;
  o2::framework::Produces<o2::aod::FLiteK0shorts> producedLiteK0shorts;
  o2::framework::Produces<o2::aod::FK0shortMasks> producedK0shortMasks;
  o2::framework::Produces<o2::aod::FK0shortExtras> producedK0shortExtras;
};

struct ConfV0Tables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("V0Tables");
  o2::framework::Configurable<int> produceLambdas{"produceLambdas", -1, "Produce Lambdas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteLambdas{"produceLiteLambdas", -1, "Produce LiteLambdas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLambdaMasks{"produceLambdaMasks", -1, "Produce LambdaMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLambdaExtras{"produceLambdaExtras", -1, "Produce LambdaExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceK0shorts{"produceK0shorts", -1, "Produce K0shorts (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteK0shorts{"produceLiteK0shorts", -1, "Produce LiteK0shorts (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceK0shortMasks{"produceK0shortMasks", -1, "Produce K0shortMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceK0shortExtras{"produceK0shortExtras", -1, "Produce K0shortExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::V0 v0Type, auto& SelectionHistName, auto& FilterHistName>
class V0Builder
{
 public:
  V0Builder() = default;
  ~V0Builder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext, T5& trackBuilder)
  {
    if constexpr (modes::isEqual(v0Type, modes::V0::kLambda) || modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        LOG(info) << "Initialize femto Lambda builder...";
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        LOG(info) << "Initialize femto AntiLambda builder...";
      }
      mProduceLambdas = utils::enableTable("FLambdas_001", table.produceLambdas.value, initContext);
      mProduceLiteLambdas = utils::enableTable("FLiteLambdas_001", table.produceLiteLambdas.value, initContext);
      mProduceLambdaMasks = utils::enableTable("FLambdaMasks_002", table.produceLambdaMasks.value, initContext);
      mProduceLambdaExtras = utils::enableTable("FLambdaExtras_001", table.produceLambdaExtras.value, initContext);

      if (mProduceLambdas && mProduceLiteLambdas) {
        LOG(fatal) << "FLambdas and FLiteLambdas are mutually exclusive -- enable only one. "
                   << "FLiteLambdas is meant to replace FLambdas at the producer stage (for better compression in derived data); "
                   << "use the dedicated converter task to reconstruct FLambdas from FLiteLambdas downstream.";
      }
      if (mProduceLambdas && !trackBuilder.producingTracks()) {
        LOG(fatal) << "FLambdas is enabled, but the track builder is not producing FTracks (full precision). "
                   << "FLambdas stores daughter indices into FTracks -- enable TrackTables.produceTracks, "
                   << "or switch to FLiteLambdas if TrackTables.produceLiteTracks is enabled instead.";
      }
      if (mProduceLiteLambdas && !trackBuilder.producingLiteTracks()) {
        LOG(fatal) << "FLiteLambdas is enabled, but the track builder is not producing FLiteTracks. "
                   << "FLiteLambdas stores daughter indices into FLiteTracks -- enable TrackTables.produceLiteTracks, "
                   << "or switch to FLambdas if TrackTables.produceTracks is enabled instead.";
      }
    }
    if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
      LOG(info) << "Initialize femto K0short builder...";
      mProduceK0shorts = utils::enableTable("FK0shorts_001", table.produceK0shorts.value, initContext);
      mProduceLiteK0shorts = utils::enableTable("FLiteK0shorts_001", table.produceLiteK0shorts.value, initContext);
      mProduceK0shortMasks = utils::enableTable("FK0shortMasks_002", table.produceK0shortMasks.value, initContext);
      mProduceK0shortExtras = utils::enableTable("FK0shortExtras_001", table.produceK0shortExtras.value, initContext);

      if (mProduceK0shorts && mProduceLiteK0shorts) {
        LOG(fatal) << "FK0shorts and FLiteK0shorts are mutually exclusive -- enable only one. "
                   << "FLiteK0shorts is meant to replace FK0shorts at the producer stage (for better compression in derived data); "
                   << "use the dedicated converter task to reconstruct FK0shorts from FLiteK0shorts downstream.";
      }
      if (mProduceK0shorts && !trackBuilder.producingTracks()) {
        LOG(fatal) << "FK0shorts is enabled, but the track builder is not producing FTracks (full precision). "
                   << "FK0shorts stores daughter indices into FTracks -- enable TrackTables.produceTracks, "
                   << "or switch to FLiteK0shorts if TrackTables.produceLiteTracks is enabled instead.";
      }
      if (mProduceLiteK0shorts && !trackBuilder.producingLiteTracks()) {
        LOG(fatal) << "FLiteK0shorts is enabled, but the track builder is not producing FLiteTracks. "
                   << "FLiteK0shorts stores daughter indices into FLiteTracks -- enable TrackTables.produceLiteTracks, "
                   << "or switch to FK0shorts if TrackTables.produceTracks is enabled instead.";
      }
    }
    if (mProduceLambdas || mProduceLiteLambdas || mProduceLambdaMasks || mProduceLambdaExtras ||
        mProduceK0shorts || mProduceLiteK0shorts || mProduceK0shortMasks || mProduceK0shortExtras) {
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

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillV0s(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts, T5& v0Products, T6 const& v0s, T7 const& tracks, T8& trackBuilder)
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

      auto posDaughter = v0.template posTrack_as<T7>();
      auto negDaughter = v0.template negTrack_as<T7>();

      posDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(posDaughter, trackProducts, collisionBuilder);
      negDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(negDaughter, trackProducts, collisionBuilder);

      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        fillLambda(collisionBuilder, v0Products, v0, 1.f, posDaughterIndex, negDaughterIndex);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        fillLambda(collisionBuilder, v0Products, v0, -1.f, posDaughterIndex, negDaughterIndex);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
        fillK0short(collisionBuilder, v0Products, v0, posDaughterIndex, negDaughterIndex);
      }
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
  void fillMcV0s(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts, T6& v0Products, T7 const& v0s, T8 const& tracks, T9& trackBuilder, T10 const& mcParticles, T11& mcBuilder, T12& mcProducts)
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

      auto posDaughter = v0.template posTrack_as<T8>();
      posDaughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kV0Daughter>(posDaughter, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

      auto negDaughter = v0.template negTrack_as<T8>();
      negDaughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kV0Daughter>(negDaughter, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

      if constexpr (modes::isEqual(v0Type, modes::V0::kLambda)) {
        fillLambda(collisionBuilder, v0Products, v0, 1.f, posDaughterIndex, negDaughterIndex);
        mcBuilder.template fillMcLambdaWithLabel<system>(v0, mcParticles, mcCols, mcProducts);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kAntiLambda)) {
        fillLambda(collisionBuilder, v0Products, v0, -1.f, posDaughterIndex, negDaughterIndex);
        mcBuilder.template fillMcLambdaWithLabel<system>(v0, mcParticles, mcCols, mcProducts);
      }
      if constexpr (modes::isEqual(v0Type, modes::V0::kK0short)) {
        fillK0short(collisionBuilder, v0Products, v0, posDaughterIndex, negDaughterIndex);
        mcBuilder.template fillMcK0shortWithLabel<system>(v0, mcParticles, mcCols, mcProducts);
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillLambda(T1& collisionBuilder, T2& v0Products, T3 const& v0, float sign, int64_t posDaughterIndex, int64_t negDaughterIndex)
  {
    float mass = 0;
    float massAnti = 0;
    float strangeTofPosDau = 0;
    float strangeTofNegDau = 0;
    if (sign > 0.f) {
      mass = v0.mLambda();
      massAnti = v0.mAntiLambda();
      strangeTofPosDau = v0.tofNSigmaLaPr();
      strangeTofNegDau = v0.tofNSigmaLaPi();
    } else {
      mass = v0.mAntiLambda();
      massAnti = v0.mLambda();
      strangeTofPosDau = v0.tofNSigmaALaPi();
      strangeTofNegDau = v0.tofNSigmaALaPr();
    }
    if (mProduceLambdas) {
      v0Products.producedLambdas(collisionBuilder.collisionIndex(),
                                 sign * v0.pt(),
                                 v0.eta(),
                                 v0.phi(),
                                 mass,
                                 posDaughterIndex,
                                 negDaughterIndex);
    }
    if (mProduceLiteLambdas) {
      v0Products.producedLiteLambdas(collisionBuilder.collisionIndex(),
                                     o2::aod::femtobase::lite::binSignedPt(sign * v0.pt()),
                                     o2::aod::femtobase::lite::binEta(v0.eta()),
                                     o2::aod::femtobase::lite::binPhi(v0.phi()),
                                     o2::aod::femtov0s::lite::binLambdaMass(mass),
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
        strangeTofPosDau,
        strangeTofNegDau,
        v0.v0radius(),
        v0.x(),
        v0.y(),
        v0.z());
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillK0short(T1& collisionBuilder, T2& v0Products, T3 const& v0, int64_t posDaughterIndex, int64_t negDaughterIndex)
  {
    if (mProduceK0shorts) {
      v0Products.producedK0shorts(collisionBuilder.collisionIndex(),
                                  v0.pt(),
                                  v0.eta(),
                                  v0.phi(),
                                  v0.mK0Short(),
                                  posDaughterIndex,
                                  negDaughterIndex);
    }
    if (mProduceLiteK0shorts) {
      v0Products.producedLiteK0shorts(collisionBuilder.collisionIndex(),
                                      o2::aod::femtobase::lite::binUnsignedPt(v0.pt()),
                                      o2::aod::femtobase::lite::binEta(v0.eta()),
                                      o2::aod::femtobase::lite::binPhi(v0.phi()),
                                      o2::aod::femtov0s::lite::binK0shortMass(v0.mK0Short()),
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
        v0.tofNSigmaK0PiPlus(),
        v0.tofNSigmaK0PiMinus(),
        v0.v0radius(),
        v0.x(),
        v0.y(),
        v0.z());
    }
  }

  [[nodiscard]] bool fillAnyTable() const { return mFillAnyTable; }
  [[nodiscard]] bool isPassThrough() const { return mV0Selection.isPassThrough(); }

 private:
  V0Selection<v0Type, SelectionHistName, FilterHistName> mV0Selection;
  bool mFillAnyTable = false;
  bool mProduceLambdas = false;
  bool mProduceLiteLambdas = false;
  bool mProduceLambdaMasks = false;
  bool mProduceLambdaExtras = false;
  bool mProduceK0shorts = false;
  bool mProduceLiteK0shorts = false;
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

    if (mLimitLambda == 0 && mLimitK0short == 0) {
      LOG(fatal) << "Both lambda limit and k0short limit are 0. Breaking...";
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewLambdas(T1 const& col, T2 const& /*lambdaTable*/, T3& partitionLambda, T4& cache)
  {
    auto lambdaSlice = partitionLambda->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return lambdaSlice.size() < mLimitLambda;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewK0shorts(T1 const& col, T2 const& /*k0shortTable*/, T3& partitionK0short, T4& cache)
  {
    auto k0shortSlice = partitionK0short->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return k0shortSlice.size() < mLimitK0short;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processLambdas(T1 const& col, T2 const& /*lambdaTable*/, T3 const& oldTrackTable, T4& partitionLambda, T5& trackBuilder, T6& cache, T7& newLambdaTable, T8& newTrackTable, T9& newCollisionTable)
  {
    auto lambdaSlice = partitionLambda->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& lambda : lambdaSlice) {

      // auto posDaughter = lambda.template posDau_as<T3>();
      // auto negDaughter = lambda.template negDau_as<T3>();
      auto posDaughter = oldTrackTable.rawIteratorAt(lambda.posDauId() - oldTrackTable.offset());
      auto negDaughter = oldTrackTable.rawIteratorAt(lambda.negDauId() - oldTrackTable.offset());

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
  void processK0shorts(T1 const& col, T2 const& /*k0shortTable*/, T3 const& oldTrackTable, T4& partitionK0short, T5& trackBuilder, T6& cache, T7& newK0shortTable, T8& newTrackTable, T9& newCollisionTable)
  {
    auto k0shortSlice = partitionK0short->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& k0short : k0shortSlice) {

      // auto posDaughter = k0short.template posDau_as<T3>();
      // auto negDaughter = k0short.template negDau_as<T3>();
      auto posDaughter = oldTrackTable.rawIteratorAt(k0short.posDauId() - oldTrackTable.offset());
      auto negDaughter = oldTrackTable.rawIteratorAt(k0short.negDauId() - oldTrackTable.offset());

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

} // namespace o2::analysis::femto::v0builder
#endif // PWGCF_FEMTO_CORE_V0BUILDER_H_
