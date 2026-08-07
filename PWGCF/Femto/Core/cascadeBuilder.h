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

namespace o2::analysis::femto::cascadebuilder
{
struct ConfCascadeFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massLambdaMin{"massLambdaMin", 1.0f, "Minimum mass of Lambda daughters in cascade decay"};
  o2::framework::Configurable<float> massLambdaMax{"massLambdaMax", 1.2f, "Maximum mass of Lambda daughters in cascade decay"};
  o2::framework::Configurable<float> massXiMin{"massXiMin", 1.2f, "Minimum Xi mass"};
  o2::framework::Configurable<float> massXiMax{"massXiMax", 1.4f, "Maximum Xi mass"};
  o2::framework::Configurable<bool> rejectHypothesisOmega{"rejectHypothesisOmega", false, "Rejection of Omega hypothesis for Xi candidates"};
  o2::framework::Configurable<float> rejectMassOmegaMin{"rejectMassOmegaMin", 1.668f, "Reject minimum Omega mass for Xi hypothesis"};
  o2::framework::Configurable<float> rejectMassOmegaMax{"rejectMassOmegaMax", 1.676f, "Reject maximum Omega mass for Xi hypothesis"};

  o2::framework::Configurable<float> massOmegaMin{"massOmegaMin", 1.5f, "Minimum Omega mass"};
  o2::framework::Configurable<float> massOmegaMax{"massOmegaMax", 1.9f, "Maximum Omega mass"};
  o2::framework::Configurable<bool> rejectHypothesisXi{"rejectHypothesisXi", false, "Rejection of Xi hypothesis for Omega candidates"};
  o2::framework::Configurable<float> rejectMassXiMin{"rejectMassXiMin", 1.317f, "Reject Minimum Xi mass for Omega hypothesis"};
  o2::framework::Configurable<float> rejectMassXiMax{"rejectMassXiMax", 1.325f, "Reject Maximum Xi mass for Omega hypothesis"};
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CASCADE_DEFAULT_BITS                                                                                                                                           \
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, all Cascades are passed through. Bits for all selections are stored."};                \
  o2::framework::Configurable<std::vector<float>> cascadeCpaMin{"cascadeCpaMin", {0.95f}, "Minimum cosine of pointing angle"};                                         \
  o2::framework::Configurable<std::vector<float>> cascadeTransRadMin{"cascadeTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                                  \
  o2::framework::Configurable<std::vector<float>> cascadeDcaDauMax{"cascadeDcaDauMax", {0.25f}, "Maximum DCA between the daughters at decay vertex (cm)"};             \
  o2::framework::Configurable<std::vector<float>> lambdaCpaMin{"lambdaCpaMin", {0.78f}, "Minimum cosine of pointing angle"};                                           \
  o2::framework::Configurable<std::vector<float>> lambdaTransRadMin{"lambdaTransRadMin", {0.9f}, "Minimum transverse radius (cm)"};                                    \
  o2::framework::Configurable<std::vector<float>> lambdaDcaDauMax{"lambdaDcaDauMax", {0.5f}, "Maximum DCA between the daughters at decay vertex (cm)"};                \
  o2::framework::Configurable<std::vector<float>> lambdaDcaToPvMin{"lambdaDcaToPvMin", {0.3f}, "Minimum DCA between the lambda and primary vertex"};                   \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.8f}, "Maximum |eta| of all daughters"};                                              \
  o2::framework::Configurable<std::vector<float>> dauAbsDcaxyMin{"dauAbsDcaxyMin", {0.05f}, "Minimum |DCAxy| of the daughters and bachelor from primary vertex (cm)"}; \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {80.f}, "Minimum number of TPC clusters for daughter tracks"};                \
  o2::framework::Configurable<std::vector<float>> posDauTpc{"posDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for positive daughter tracks"};                      \
  o2::framework::Configurable<std::vector<float>> negDauTpc{"negDauTpc", {5.f}, "Maximum |nsimga_Pion/Proton| TPC for negative daughter tracks"};                      \
  o2::framework::Configurable<std::vector<float>> posDauTof{"posDauTof", {}, "Maximum |nsimga_Pion/Proton| TOF for positive daughter tracks"};                         \
  o2::framework::Configurable<std::vector<float>> negDauTof{"negDauTof", {}, "Maximum |nsigma_Pion/Proton| TOF for negative daughter tracks"};                         \
  o2::framework::Configurable<bool> requireTof{"requireTof", false, "If true, TOF PID is a mandatory selection"};                                                      \
  o2::framework::Configurable<bool> keepTracksWithoutTof{"keepTracksWithoutTof", true, "If true, candidates whose daughters have no TOF signal are kept"};

struct ConfXiBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("XiBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcPion{"bachelorTpcPion", {5.f}, "Maximum |nsimga_Pion| TPC for bachelor tracks"};
  o2::framework::Configurable<std::vector<float>> bachelorTofPion{"bachelorTofPion", {}, "Maximum |nsimga_Pion| TOF for bachelor tracks"};
};

struct ConfOmegaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaBits");
  CASCADE_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> bachelorTpcKaon{"bachelorTpcKaon", {5.f}, "Maximum |nsimga_Kaon| TPC for bachelor tracks"};
  o2::framework::Configurable<std::vector<float>> bachelorTofKaon{"bachelorTofKaon", {}, "Maximum |nsimga_Kaon| TOF for bachelor tracks"};
};

#undef CASCADE_DEFAULT_BITS

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define CASCADE_DEFAULT_SELECTION(defaultMassMin, defaultMassMax, defaultPdgCode)                                                         \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", (defaultPdgCode), "Cascade PDG code. Set sign to +1 to select antiparticle"}; \
  o2::framework::Configurable<int> sign{"sign", -1, "Sign of the charge of the Cascade"};                                                 \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                                   \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                                                 \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                                              \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                                               \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                                \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                                   \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass for Cascade"};                          \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass for Cascade"};                          \
  o2::framework::Configurable<o2::analysis::femto::datatypes::CascadeMaskType> mask{"mask", 0x0, "Bitmask for cascade selection"};

struct ConfXiSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("XiSelection");
  CASCADE_DEFAULT_SELECTION(1.22, 1.42, 3312)
};

struct ConfOmegaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaSelection");
  CASCADE_DEFAULT_SELECTION(1.57, 1.77, 3334)
};

#undef CASCADE_DEFAULT_SELECTION

/// The different selections this task is capable of doing
enum CascadeSels {
  // selections for cascades
  kCascadeCpaMin,      ///< Min. CPA (cosine pointing angle)
  kCascadeDcaDaughMax, ///< Max. DCA of the daughers at decay vertex
  kCascadeTransRadMin, ///< max. transverse radius

  // selection for lambda daughter
  kLambdaCpaMin,      ///< Min. CPA of the lambda
  kLambdaDcaDauMax,   ///< Max. DCA between the lambda daughters at lambda decay vertex
  kLambdaTransRadMin, ///< Min. tranverse radius of the lambda
  kLambdaDcaToPvMin,  ///< Min. DCA of the lambda to the primary vertex

  // selection for bachelor/daugthers
  kDauAbsEtaMax,   ///< Max. |eta| of daughter tracks
  kDauTpcClsMin,   ///< Min. number of TPC clusters of daughters/bachelor
  kDauAbsDcaxyMin, ///< Min. |DCAxy| of the daughers and bachelor from primary vertex

  // PID selection for cascade bachelor
  kBachelorTpcPion, ///< TPC Pion PID for bachelor
  kBachelorTpcKaon, ///< TPC Kaon PID for bachelor
  kBachelorTofPion, ///< TOF Pion PID for bachelor
  kBachelorTofKaon, ///< TOF Kaon PID for bachelor
                    ///
  // PID selection for lambda daughers
  kPosDauTpc, ///< TPC PID for positive daughter
  kNegDauTpc, ///< TPC PID for negative daughter
  kPosDauTof, ///< TOF PID for positive daughter
  kNegDauTof, ///< TOF PID for negative daughter

  kCascadeSelsMax
};

constexpr char XiSelHistName[] = "hXiSelection";
constexpr char OmegaSelHistName[] = "hOmegaSelection";
constexpr char CascadeSelsName[] = "Cascade Selection Object";
const std::unordered_map<CascadeSels, std::string> cascadeSelectionNames = {
  {kCascadeCpaMin, "Cascade CPA Min"},
  {kCascadeDcaDaughMax, "Cascade DCA Daughters Max"},
  {kCascadeTransRadMin, "Cascade Transverse Radius Min"},

  {kLambdaCpaMin, "Lambda CPA Min"},
  {kLambdaDcaDauMax, "Lambda DCA Daughter Max"},
  {kLambdaTransRadMin, "Lambda Transverse Radius Min"},
  {kLambdaDcaToPvMin, "Lambda DCA to PV Min"},

  {kDauAbsEtaMax, "Daughter Abs Eta Max"},
  {kDauTpcClsMin, "Daughter TPC Clusters Min"},
  {kDauAbsDcaxyMin, "Daughter |DCAxy| Min"},

  {kBachelorTpcPion, "Bachelor TPC Pion PID"},
  {kBachelorTpcKaon, "Bachelor TPC Kaon PID"},
  {kBachelorTofPion, "Bachelor TOF Pion PID"},
  {kBachelorTofKaon, "Bachelor TOF Kaon PID"},
  {kPosDauTpc, "Positive Daughter TPC PID"},
  {kNegDauTpc, "Negative Daughter TPC PID"},
  {kPosDauTof, "Positive Daughter TOF PID"},
  {kNegDauTof, "Negative Daughter TOF PID"}};

/// enum for all cascade pre-filters (evaluated in checkFilters, before the selection bitmask)
enum CascadeFilters {
  kPtMin,
  kPtMax,
  kEtaMin,
  kEtaMax,
  kPhiMin,
  kPhiMax,
  kLambdaMassMin,
  kLambdaMassMax,
  kXiMassMin,
  kXiMassMax,
  kRejectionOmegaMass,
  kOmegaMassMin,
  kOmegaMassMax,
  kRejectionXiMass,
  kCascadeFiltersMax
};

constexpr char XiFilterHistName[] = "hXiFilters";
constexpr char OmegaFilterHistName[] = "hOmegaFilters";
const std::unordered_map<CascadeFilters, std::string> cascadeFilterNames = {
  {kPtMin, "ptMin"},
  {kPtMax, "ptMax"},
  {kEtaMin, "etaMin"},
  {kEtaMax, "etaMax"},
  {kPhiMin, "phiMin"},
  {kPhiMax, "phiMax"},
  {kLambdaMassMin, "lambdaMassMin"},
  {kLambdaMassMax, "lambdaMassMax"},
  {kXiMassMin, "xiMassMin"},
  {kXiMassMax, "xiMassMax"},
  {kRejectionOmegaMass, "rejectOmega"},
  {kOmegaMassMin, "omegaMassMin"},
  {kOmegaMassMax, "omegaMassMax"},
  {kRejectionXiMass, "rejectXi"},
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <modes::Cascade cascadeType, auto& SelectionHistName, auto& FilterHistName>
class CascadeSelection : public baseselection::BaseSelection<float, o2::analysis::femto::datatypes::CascadeMaskType, kCascadeSelsMax>
{
 public:
  CascadeSelection() = default;
  ~CascadeSelection() override = default;

  template <typename T1, typename T2>
  void configure(o2::framework::HistogramRegistry* registry, T1 const& config, T2 const& filter)
  {
    this->init(config.passThrough.value);

    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;
    mLambdaMassMin = filter.massLambdaMin.value;
    mLambdaMassMax = filter.massLambdaMax.value;
    mRequireTof = config.requireTof.value;
    mKeepTracksWithoutTof = config.keepTracksWithoutTof.value;

    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      mXiMassLowerLimit = filter.massXiMin.value;
      mXiMassUpperLimit = filter.massXiMax.value;
      mRejectOmegaHypothesis = filter.rejectHypothesisOmega.value;
      mOmegaMassLowerLimit = filter.rejectMassOmegaMin.value;
      mOmegaMassUpperLimit = filter.rejectMassOmegaMax.value;
      this->addSelection(kBachelorTpcPion, cascadeSelectionNames.at(kBachelorTpcPion), config.bachelorTpcPion.value, limits::kAbsUpperLimit, true, true, false);
      this->addSelection(kBachelorTofPion, cascadeSelectionNames.at(kBachelorTofPion), config.bachelorTofPion.value, limits::kAbsUpperLimit, true, mRequireTof, false);
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      mOmegaMassLowerLimit = filter.massOmegaMin.value;
      mOmegaMassUpperLimit = filter.massOmegaMax.value;
      mRejectXiHypothesis = filter.rejectHypothesisXi.value;
      mXiMassLowerLimit = filter.rejectMassXiMin.value;
      mXiMassUpperLimit = filter.rejectMassXiMax.value;
      this->addSelection(kBachelorTpcKaon, cascadeSelectionNames.at(kBachelorTpcKaon), config.bachelorTpcKaon.value, limits::kAbsUpperLimit, true, true, false);
      this->addSelection(kBachelorTofKaon, cascadeSelectionNames.at(kBachelorTofKaon), config.bachelorTofKaon.value, limits::kAbsUpperLimit, true, mRequireTof, false);
    }

    this->addSelection(kPosDauTpc, cascadeSelectionNames.at(kPosDauTpc), config.posDauTpc.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kNegDauTpc, cascadeSelectionNames.at(kNegDauTpc), config.negDauTpc.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kPosDauTof, cascadeSelectionNames.at(kPosDauTof), config.posDauTof.value, limits::kAbsUpperLimit, true, mRequireTof, false);
    this->addSelection(kNegDauTof, cascadeSelectionNames.at(kNegDauTof), config.negDauTof.value, limits::kAbsUpperLimit, true, mRequireTof, false);

    this->addSelection(kCascadeCpaMin, cascadeSelectionNames.at(kCascadeCpaMin), config.cascadeCpaMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kCascadeTransRadMin, cascadeSelectionNames.at(kCascadeTransRadMin), config.cascadeTransRadMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kCascadeDcaDaughMax, cascadeSelectionNames.at(kCascadeDcaDaughMax), config.cascadeDcaDauMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kLambdaCpaMin, cascadeSelectionNames.at(kLambdaCpaMin), config.lambdaCpaMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kLambdaTransRadMin, cascadeSelectionNames.at(kLambdaTransRadMin), config.lambdaTransRadMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kLambdaDcaDauMax, cascadeSelectionNames.at(kLambdaDcaDauMax), config.lambdaDcaDauMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kLambdaDcaToPvMin, cascadeSelectionNames.at(kLambdaDcaToPvMin), config.lambdaDcaToPvMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kDauAbsEtaMax, cascadeSelectionNames.at(kDauAbsEtaMax), config.dauAbsEtaMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kDauAbsDcaxyMin, cascadeSelectionNames.at(kDauAbsDcaxyMin), config.dauAbsDcaxyMin.value, limits::kAbsLowerLimit, true, true, false);
    this->addSelection(kDauTpcClsMin, cascadeSelectionNames.at(kDauTpcClsMin), config.dauTpcClustersMin.value, limits::kLowerLimit, true, true, false);

    this->setupSelectionHistogram<SelectionHistName>(registry);
    this->template setupFilterHistogram<FilterHistName>(
      registry,
      {
        {cascadeFilterNames.at(kPtMin), mPtMin},
        {cascadeFilterNames.at(kPtMax), mPtMax},
        {cascadeFilterNames.at(kEtaMin), mEtaMin},
        {cascadeFilterNames.at(kEtaMax), mEtaMax},
        {cascadeFilterNames.at(kPhiMin), mPhiMin},
        {cascadeFilterNames.at(kPhiMax), mPhiMax},
        {cascadeFilterNames.at(kLambdaMassMin), mLambdaMassMin},
        {cascadeFilterNames.at(kLambdaMassMax), mLambdaMassMax},
        {cascadeFilterNames.at(kXiMassMin), mXiMassLowerLimit},
        {cascadeFilterNames.at(kXiMassMax), mXiMassUpperLimit},
        {cascadeFilterNames.at(kRejectionOmegaMass), mRejectOmegaHypothesis ? 1.f : 0.f},
        {cascadeFilterNames.at(kOmegaMassMin), mOmegaMassLowerLimit},
        {cascadeFilterNames.at(kOmegaMassMax), mOmegaMassUpperLimit},
        {cascadeFilterNames.at(kRejectionXiMass), mRejectXiHypothesis ? 1.f : 0.f},
      });
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
    std::array<float, 3> etaDaughters = {std::fabs(cascade.bacheloreta()), std::fabs(cascade.positiveeta()), std::fabs(cascade.negativeeta())};
    this->evaluateObservable(kDauAbsEtaMax, *std::max_element(etaDaughters.begin(), etaDaughters.end()));

    std::array<float, 3> dcaDaughters = {std::fabs(cascade.dcabachtopv()), std::fabs(cascade.dcapostopv()), std::fabs(cascade.dcanegtopv())};
    this->evaluateObservable(kDauAbsDcaxyMin, *std::min_element(dcaDaughters.begin(), dcaDaughters.end()));

    std::array<float, 3> clustersDaughters = {1.f * bachelor.tpcNClsFound(), 1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauTpcClsMin, *std::min_element(clustersDaughters.begin(), clustersDaughters.end()));

    // pid selections
    // TPC nSigma comes from the daughter track, TOF nSigma and the has-TOF flags from the cascade candidate
    // if a daughter has no TOF signal, feed 0 so the bit passes any limit (opt-in via keepTracksWithoutTof)
    auto evaluatePid = [this](CascadeSels tpcBit, float tpcNSigma,
                              CascadeSels tofBit, float tofNSigma, bool hasTof) {
      this->evaluateObservable(tpcBit, tpcNSigma);
      if (hasTof) {
        this->evaluateObservable(tofBit, tofNSigma);
      } else if (mKeepTracksWithoutTof) {
        this->evaluateObservable(tofBit, 0.f);
      }
    };

    const bool bachHasTof = cascade.bachelorHasTOF();
    const bool posHasTof = cascade.positiveHasTOF();
    const bool negHasTof = cascade.negativeHasTOF();

    // bachelor: pion for Xi, kaon for Omega
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      evaluatePid(kBachelorTpcPion, bachelor.tpcNSigmaPi(),
                  kBachelorTofPion, cascade.tofNSigmaXiPi(), bachHasTof);
    } else if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      evaluatePid(kBachelorTpcKaon, bachelor.tpcNSigmaKa(),
                  kBachelorTofKaon, cascade.tofNSigmaOmKa(), bachHasTof);
    }

    // v0 daughters: charge of the cascade fixes the Lambda vs. AntiLambda hypothesis
    if (cascade.sign() == 0) {
      LOG(warn) << "Encountered Cascade candidate with 0 charge";
    } else {
      // sign < 0: Xi-/Omega- -> Lambda -> p pi-   (pos = proton, neg = pion)
      // sign > 0: Xi+/Omega+ -> AntiLambda        (pos = pion,   neg = antiproton)
      const bool isMatter = cascade.sign() < 0;

      const float tpcPosDau = isMatter ? posDaughter.tpcNSigmaPr() : posDaughter.tpcNSigmaPi();
      const float tpcNegDau = isMatter ? negDaughter.tpcNSigmaPi() : negDaughter.tpcNSigmaPr();

      float tofPosDau = 0.f;
      float tofNegDau = 0.f;
      if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
        tofPosDau = isMatter ? cascade.tofNSigmaXiLaPr() : cascade.tofNSigmaXiLaPi();
        tofNegDau = isMatter ? cascade.tofNSigmaXiLaPi() : cascade.tofNSigmaXiLaPr();
      } else if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
        tofPosDau = isMatter ? cascade.tofNSigmaOmLaPr() : cascade.tofNSigmaOmLaPi();
        tofNegDau = isMatter ? cascade.tofNSigmaOmLaPi() : cascade.tofNSigmaOmLaPr();
      }

      evaluatePid(kPosDauTpc, tpcPosDau, kPosDauTof, tofPosDau, posHasTof);
      evaluatePid(kNegDauTpc, tpcNegDau, kNegDauTof, tofNegDau, negHasTof);
    }

    this->assembleBitmask<SelectionHistName>();
  };

  template <typename T>
  bool checkFilters(const T& cascade) const
  {
    bool pass = true;
    bool p = false;

    // kinematics
    p = cascade.pt() > mPtMin;
    this->template fillFilter<FilterHistName>(kPtMin, p);
    pass &= p;

    p = cascade.pt() < mPtMax;
    this->template fillFilter<FilterHistName>(kPtMax, p);
    pass &= p;

    p = cascade.eta() > mEtaMin;
    this->template fillFilter<FilterHistName>(kEtaMin, p);
    pass &= p;

    p = cascade.eta() < mEtaMax;
    this->template fillFilter<FilterHistName>(kEtaMax, p);
    pass &= p;

    p = cascade.phi() > mPhiMin;
    this->template fillFilter<FilterHistName>(kPhiMin, p);
    pass &= p;

    p = cascade.phi() < mPhiMax;
    this->template fillFilter<FilterHistName>(kPhiMax, p);
    pass &= p;

    // mass of daughter lambda (gating AND-cut)
    p = cascade.mLambda() > mLambdaMassMin;
    this->template fillFilter<FilterHistName>(kLambdaMassMin, p);
    pass &= p;

    p = cascade.mLambda() < mLambdaMassMax;
    this->template fillFilter<FilterHistName>(kLambdaMassMax, p);
    pass &= p;

    // mass hypothesis: signal window is a gating AND-cut, competing hypothesis is a
    // rejection OR-cut (pass if outside its window, or if rejection is disabled);
    // the competing window's own min/max bins are diagnostic only.
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      p = cascade.mXi() > mXiMassLowerLimit;
      this->template fillFilter<FilterHistName>(kXiMassMin, p);
      pass &= p;

      p = cascade.mXi() < mXiMassUpperLimit;
      this->template fillFilter<FilterHistName>(kXiMassMax, p);
      pass &= p;

      bool const belowOmega = cascade.mOmega() < mOmegaMassLowerLimit;
      this->template fillFilter<FilterHistName>(kOmegaMassMin, belowOmega);
      bool const aboveOmega = cascade.mOmega() > mOmegaMassUpperLimit;
      this->template fillFilter<FilterHistName>(kOmegaMassMax, aboveOmega);

      p = !mRejectOmegaHypothesis || belowOmega || aboveOmega;
      this->template fillFilter<FilterHistName>(kRejectionOmegaMass, p);
      pass &= p;
    }

    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      p = cascade.mOmega() > mOmegaMassLowerLimit;
      this->template fillFilter<FilterHistName>(kOmegaMassMin, p);
      pass &= p;

      p = cascade.mOmega() < mOmegaMassUpperLimit;
      this->template fillFilter<FilterHistName>(kOmegaMassMax, p);
      pass &= p;

      bool const belowXi = cascade.mXi() < mXiMassLowerLimit;
      this->template fillFilter<FilterHistName>(kXiMassMin, belowXi);
      bool const aboveXi = cascade.mXi() > mXiMassUpperLimit;
      this->template fillFilter<FilterHistName>(kXiMassMax, aboveXi);

      p = !mRejectXiHypothesis || belowXi || aboveXi;
      this->template fillFilter<FilterHistName>(kRejectionXiMass, p);
      pass &= p;
    }

    this->template fillFilterSummary<FilterHistName>(pass);

    return this->isPassThrough() || pass;
  }

 protected:
  bool mRejectOmegaHypothesis = false;
  float mOmegaMassLowerLimit = 0.f;
  float mOmegaMassUpperLimit = 999.f;

  bool mRejectXiHypothesis = false;
  float mXiMassLowerLimit = 0.f;
  float mXiMassUpperLimit = 999.f;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -0.9f;
  float mEtaMax = 0.9f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;
  float mLambdaMassMin = 1.f;
  float mLambdaMassMax = 1.2f;
  bool mRequireTof = false;
  bool mKeepTracksWithoutTof = false;
};

struct CascadeBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FXis> producedXis;
  o2::framework::Produces<o2::aod::FLiteXis> producedLiteXis;
  o2::framework::Produces<o2::aod::FXiMasks> producedXiMasks;
  o2::framework::Produces<o2::aod::FXiExtras> producedXiExtras;
  o2::framework::Produces<o2::aod::FOmegas> producedOmegas;
  o2::framework::Produces<o2::aod::FLiteOmegas> producedLiteOmegas;
  o2::framework::Produces<o2::aod::FOmegaMasks> producedOmegaMasks;
  o2::framework::Produces<o2::aod::FOmegaExtras> producedOmegaExtras;
};

struct ConfCascadeTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeTables");
  o2::framework::Configurable<int> produceXis{"produceXis", -1, "Produce Xis (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteXis{"produceLiteXis", -1, "Produce LiteXis (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceXiMasks{"produceXiMasks", -1, "Produce XiMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceXiExtras{"produceXiExtras", -1, "Produce XiExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOmegas{"produceOmegas", -1, "Produce Omegas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteOmegas{"produceLiteOmegas", -1, "Produce LiteOmegas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOmegaMasks{"produceOmegaMasks", -1, "Produce OmegaMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOmegaExtras{"produceOmegaExtras", -1, "Produce OmegaExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::Cascade cascadeType, auto& SelectionHistName, auto& FilterHistName>
class CascadeBuilder
{
 public:
  CascadeBuilder() = default;
  ~CascadeBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext, T5& trackBuilder)
  {
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      LOG(info) << "Initialize femto Xi builder...";
      mProduceLiteXis = utils::enableTable("FLiteXis_001", table.produceLiteXis.value, initContext);
      mProduceXis = utils::enableTable("FXis_001", table.produceXis.value, initContext);
      mProduceXiMasks = utils::enableTable("FXiMasks_002", table.produceXiMasks.value, initContext);
      mProduceXiExtras = utils::enableTable("FXiExtras_001", table.produceXiExtras.value, initContext);

      if (mProduceXis && mProduceLiteXis) {
        LOG(fatal) << "FXis and FLiteXis are mutually exclusive -- enable only one. "
                   << "FLiteXis is meant to replace FXis at the producer stage (for better compression in derived data); "
                   << "use the dedicated converter task to reconstruct FXis from FLiteXis downstream.";
      }
      if (mProduceXis && !trackBuilder.producingTracks()) {
        LOG(fatal) << "FXis is enabled, but the track builder is not producing FTracks (full precision). "
                   << "FXis stores daughter indices into FTracks -- enable TrackTables.produceTracks, "
                   << "or switch to FLiteXis if TrackTables.produceLiteTracks is enabled instead.";
      }
      if (mProduceLiteXis && !trackBuilder.producingLiteTracks()) {
        LOG(fatal) << "FLiteXis is enabled, but the track builder is not producing FLiteTracks. "
                   << "FLiteXis stores daughter indices into FLiteTracks -- enable TrackTables.produceLiteTracks, "
                   << "or switch to FXis if TrackTables.produceTracks is enabled instead.";
      }
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      LOG(info) << "Initialize femto Omega builder...";
      mProduceOmegas = utils::enableTable("FOmegas_001", table.produceOmegas.value, initContext);
      mProduceLiteOmegas = utils::enableTable("FLiteOmegas_001", table.produceLiteOmegas.value, initContext);
      mProduceOmegaMasks = utils::enableTable("FOmegaMasks_002", table.produceOmegaMasks.value, initContext);
      mProduceOmegaExtras = utils::enableTable("FOmegaExtras_001", table.produceOmegaExtras.value, initContext);

      if (mProduceOmegas && mProduceLiteOmegas) {
        LOG(fatal) << "FOmegas and FLiteOmegas are mutually exclusive -- enable only one. "
                   << "FLiteOmegas is meant to replace FOmegas at the producer stage (for better compression in derived data); "
                   << "use the dedicated converter task to reconstruct FOmegas from FLiteOmegas downstream.";
      }
      if (mProduceOmegas && !trackBuilder.producingTracks()) {
        LOG(fatal) << "FOmegas is enabled, but the track builder is not producing FTracks (full precision). "
                   << "FOmegas stores daughter indices into FTracks -- enable TrackTables.produceTracks, "
                   << "or switch to FLiteOmegas if TrackTables.produceLiteTracks is enabled instead.";
      }
      if (mProduceLiteOmegas && !trackBuilder.producingLiteTracks()) {
        LOG(fatal) << "FLiteOmegas is enabled, but the track builder is not producing FLiteTracks. "
                   << "FLiteOmegas stores daughter indices into FLiteTracks -- enable TrackTables.produceLiteTracks, "
                   << "or switch to FOmegas if TrackTables.produceTracks is enabled instead.";
      }
    }

    if (mProduceXis || mProduceLiteXis || mProduceXiExtras || mProduceXiMasks ||
        mProduceOmegas || mProduceLiteOmegas || mProduceOmegaMasks || mProduceOmegaExtras) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mCascadeSelection.configure(registry, config, filter);
    mCascadeSelection.printSelections(CascadeSelsName);
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillCascades(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts, T5& cascadeProducts, T6 const& cascades, T7 const& tracks, T8& trackBuilder)
  {
    if (!mFillAnyTable) {
      return;
    }

    int64_t bachelorIndex = 0;
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& cascade : cascades) {
      if (!mCascadeSelection.checkFilters(cascade)) {
        continue;
      }
      mCascadeSelection.applySelections(cascade, tracks, col);
      if (!mCascadeSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillCollision<system>(collisionProducts, col);

      auto bachelor = cascade.template bachelor_as<T7>();
      bachelorIndex = trackBuilder.template getDaughterIndex<modes::Track::kCascadeBachelor>(bachelor, trackProducts, collisionBuilder);

      auto posDaughter = cascade.template posTrack_as<T7>();
      posDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(posDaughter, trackProducts, collisionBuilder);

      auto negDaughter = cascade.template negTrack_as<T7>();
      negDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(negDaughter, trackProducts, collisionBuilder);

      fillCascade(collisionBuilder, cascadeProducts, cascade, col, bachelorIndex, posDaughterIndex, negDaughterIndex);
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
  void fillMcCascades(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts, T6& cascadeProducts, T7 const& cascades, T8 const& tracks, T9& trackBuilder, T10 const& mcParticles, T11& mcBuilder, T12& mcProducts)
  {
    if (!mFillAnyTable) {
      return;
    }

    int64_t bachelorIndex = 0;
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& cascade : cascades) {
      if (!mCascadeSelection.checkFilters(cascade)) {
        continue;
      }
      mCascadeSelection.applySelections(cascade, tracks, col);
      if (!mCascadeSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);

      auto bachelor = cascade.template bachelor_as<T8>();
      bachelorIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCascadeBachelor>(bachelor, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

      auto posDaughter = cascade.template posTrack_as<T8>();
      posDaughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kV0Daughter>(posDaughter, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

      auto negDaughter = cascade.template negTrack_as<T8>();
      negDaughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kV0Daughter>(negDaughter, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

      fillCascade(collisionBuilder, cascadeProducts, cascade, col, bachelorIndex, posDaughterIndex, negDaughterIndex);
      if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
        mcBuilder.template fillMcXiWithLabel<system>(cascade, mcParticles, mcCols, mcProducts);
      }
      if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
        mcBuilder.template fillMcOmegaWithLabel<system>(cascade, mcParticles, mcCols, mcProducts);
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillCascade(T1& collisionBuilder, T2& cascadeProducts, T3 const& cascade, T4 const& col, int bachelorIndex, int posDaughterIndex, int negDaughterIndex)
  {
    float strangeTofBachelor = 0.f;
    float strangeTofPosDau = 0.f;
    float strangeTofNegDau = 0.f;
    const bool isMatter = cascade.sign() < 0; // Xi-/Omega- -> Lambda -> p pi- (pos=proton, neg=pion)

    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      strangeTofBachelor = cascade.tofNSigmaXiPi();
      strangeTofPosDau = isMatter ? cascade.tofNSigmaXiLaPr() : cascade.tofNSigmaXiLaPi();
      strangeTofNegDau = isMatter ? cascade.tofNSigmaXiLaPi() : cascade.tofNSigmaXiLaPr();
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      strangeTofBachelor = cascade.tofNSigmaOmKa();
      strangeTofPosDau = isMatter ? cascade.tofNSigmaOmLaPr() : cascade.tofNSigmaOmLaPi();
      strangeTofNegDau = isMatter ? cascade.tofNSigmaOmLaPi() : cascade.tofNSigmaOmLaPr();
    }

    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kXi)) {
      if (mProduceXis) {
        cascadeProducts.producedXis(collisionBuilder.collisionIndex(),
                                    cascade.sign() * cascade.pt(),
                                    cascade.eta(),
                                    cascade.phi(),
                                    cascade.mXi(),
                                    bachelorIndex,
                                    posDaughterIndex,
                                    negDaughterIndex);
      }
      if (mProduceLiteXis) {
        cascadeProducts.producedLiteXis(collisionBuilder.collisionIndex(),
                                        o2::aod::femtobase::lite::binSignedPt(cascade.sign() * cascade.pt()),
                                        o2::aod::femtobase::lite::binEta(cascade.eta()),
                                        o2::aod::femtobase::lite::binPhi(cascade.phi()),
                                        o2::aod::femtocascades::lite::binXiMass(cascade.mXi()),
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
          cascade.mLambda(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posX(), col.posY(), col.posZ()),
          strangeTofBachelor,
          strangeTofPosDau,
          strangeTofNegDau);
      }
    }
    if constexpr (modes::isEqual(cascadeType, modes::Cascade::kOmega)) {
      if (mProduceOmegas) {
        cascadeProducts.producedOmegas(collisionBuilder.collisionIndex(),
                                       cascade.sign() * cascade.pt(),
                                       cascade.eta(),
                                       cascade.phi(),
                                       cascade.mOmega(),
                                       bachelorIndex,
                                       posDaughterIndex,
                                       negDaughterIndex);
      }
      if (mProduceLiteOmegas) {
        cascadeProducts.producedLiteOmegas(collisionBuilder.collisionIndex(),
                                           o2::aod::femtobase::lite::binSignedPt(cascade.sign() * cascade.pt()),
                                           o2::aod::femtobase::lite::binEta(cascade.eta()),
                                           o2::aod::femtobase::lite::binPhi(cascade.phi()),
                                           o2::aod::femtocascades::lite::binOmegaMass(cascade.mOmega()),
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
          cascade.mLambda(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posX(), col.posY(), col.posZ()),
          strangeTofBachelor,
          strangeTofPosDau,
          strangeTofNegDau);
      }
    }
  }

  [[nodiscard]] bool fillAnyTable() const { return mFillAnyTable; }
  [[nodiscard]] bool isPassThrough() const { return mCascadeSelection.isPassThrough(); }

 private:
  CascadeSelection<cascadeType, SelectionHistName, FilterHistName> mCascadeSelection;
  bool mFillAnyTable = false;
  bool mProduceXis = false;
  bool mProduceLiteXis = false;
  bool mProduceXiMasks = false;
  bool mProduceXiExtras = false;
  bool mProduceOmegas = false;
  bool mProduceLiteOmegas = false;
  bool mProduceOmegaMasks = false;
  bool mProduceOmegaExtras = false;
};

struct ConfCascadeTablesDerivedToDerived : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CascadeTables");
  o2::framework::Configurable<int> limitXi{"limitXi", 1, "At least this many xi need to be in the collision"};
  o2::framework::Configurable<int> limitOmega{"limitOmega", 0, "At least this many omega need to be in the collision"};
};

struct CascadeBuilderDerivedToDerivedProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::StoredFXis> producedXis;
  o2::framework::Produces<o2::aod::StoredFXiMasks> producedXiMasks;
  o2::framework::Produces<o2::aod::StoredFOmegas> producedOmegas;
  o2::framework::Produces<o2::aod::StoredFOmegaMasks> producedOmegaMasks;
};

class CascadeBuilderDerivedToDerived
{
 public:
  CascadeBuilderDerivedToDerived() = default;
  ~CascadeBuilderDerivedToDerived() = default;

  template <typename T>
  void init(T& config)
  {
    mLimitXi = config.limitXi.value;
    mLimitOmega = config.limitOmega.value;

    if (mLimitXi == 0 && mLimitOmega == 0) {
      LOG(fatal) << "Both xi limit and omega limit are 0. Breaking...";
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewXis(T1 const& col, T2 const& /*xiTable*/, T3& partitionXi, T4& cache)
  {
    auto xiSlice = partitionXi->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return xiSlice.size() < mLimitXi;
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewOmegas(T1 const& col, T2 const& /*omegaTable*/, T3& partitionOmega, T4& cache)
  {
    auto omegaSlice = partitionOmega->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return omegaSlice.size() < mLimitOmega;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processXis(T1 const& col, T2 const& /*xiTable*/, T3 const& oldTrackTable, T4& partitionXi, T5& trackBuilder, T6& cache, T7& newXiTable, T8& newTrackTable, T9& newCollisionTable)
  {
    auto xiSlice = partitionXi->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& xi : xiSlice) {

      // auto bachelor = xi.template bachelor_as<T3>();
      // auto posDaughter = xi.template posDau_as<T3>();
      // auto negDaughter = xi.template negDau_as<T3>();
      auto bachelor = oldTrackTable.rawIteratorAt(xi.bachelorId() - oldTrackTable.offset());
      auto posDaughter = oldTrackTable.rawIteratorAt(xi.posDauId() - oldTrackTable.offset());
      auto negDaughter = oldTrackTable.rawIteratorAt(xi.negDauId() - oldTrackTable.offset());

      int bachelorIndex = trackBuilder.getDaughterIndex(bachelor, newTrackTable, newCollisionTable);
      int posDaughterIndex = trackBuilder.getDaughterIndex(posDaughter, newTrackTable, newCollisionTable);
      int negDaughterIndex = trackBuilder.getDaughterIndex(negDaughter, newTrackTable, newCollisionTable);

      newXiTable.producedXis(newCollisionTable.producedCollision.lastIndex(),
                             xi.signedPt(),
                             xi.eta(),
                             xi.phi(),
                             xi.mass(),
                             bachelorIndex,
                             posDaughterIndex,
                             negDaughterIndex);
      newXiTable.producedXiMasks(xi.mask());
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processOmegas(T1 const& col, T2 const& /*omegaTable*/, T3 const& oldTrackTable, T4& partitionOmega, T5& trackBuilder, T6& cache, T7& newOmegaTable, T8& newTrackTable, T9& newCollisionTable)
  {
    auto omegaSlice = partitionOmega->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& omega : omegaSlice) {

      // auto bachelor = omega.template bachelor_as<T3>();
      // auto posDaughter = omega.template posDau_as<T3>();
      // auto negDaughter = omega.template negDau_as<T3>();
      auto bachelor = oldTrackTable.rawIteratorAt(omega.bachelorId() - oldTrackTable.offset());
      auto posDaughter = oldTrackTable.rawIteratorAt(omega.posDauId() - oldTrackTable.offset());
      auto negDaughter = oldTrackTable.rawIteratorAt(omega.negDauId() - oldTrackTable.offset());

      int bachelorIndex = trackBuilder.getDaughterIndex(bachelor, newTrackTable, newCollisionTable);
      int posDaughterIndex = trackBuilder.getDaughterIndex(posDaughter, newTrackTable, newCollisionTable);
      int negDaughterIndex = trackBuilder.getDaughterIndex(negDaughter, newTrackTable, newCollisionTable);

      newOmegaTable.producedOmegas(newCollisionTable.producedCollision.lastIndex(),
                                   omega.signedPt(),
                                   omega.eta(),
                                   omega.phi(),
                                   omega.mass(),
                                   bachelorIndex,
                                   posDaughterIndex,
                                   negDaughterIndex);
      newOmegaTable.producedOmegaMasks(omega.mask());
    }
  }

 private:
  int mLimitXi = 0;
  int mLimitOmega = 0;
};

} // namespace o2::analysis::femto::cascadebuilder
#endif // PWGCF_FEMTO_CORE_CASCADEBUILDER_H_
