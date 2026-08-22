// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file charmHadronBuilder.h
/// \brief charm hadron builder
/// \author Igor Ptak, WUT, igor.tomasz.ptak@cern.ch

#ifndef PWGCF_FEMTO_CORE_CHARMHADRONBUILDER_H_
#define PWGCF_FEMTO_CORE_CHARMHADRONBUILDER_H_

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto::charmhadronbuilder
{
// number of ML probability classes (background, prompt, non-prompt) provided by the HF ML selector
constexpr std::size_t NSizeMLScore{3u};

// filter applied in the producer task
struct ConfCharmHadronFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CharmHadronFilters");
  // kinematic windows
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 24.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.8f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.8f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  // rapidity acceptance (HF convention: y, not eta)
  o2::framework::Configurable<bool> useYCut{"useYCut", true, "cut on y (true) or eta (false)"};
  o2::framework::Configurable<float> yMin{"yMin", -0.8f, "Minimum rapidity"};
  o2::framework::Configurable<float> yMax{"yMax", 0.8f, "Maximum rapidity"};
  // invariant-mass window
  o2::framework::Configurable<float> massMin{"massMin", 1.7f, "Minimum invariant mass"};
  o2::framework::Configurable<float> massMax{"massMax", 2.0f, "Maximum invariant mass"};
  o2::framework::Configurable<bool> rejectAmbiguousHypothesis{"rejectAmbiguousHypothesis", false, "Reject candidates selected under both particle and antiparticle hypotheses"};
};

// derived selection bits for D0s
struct ConfD0Bits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("D0Bits");
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, all D0s are passed through. Bits for all selections are stored."};
  o2::framework::Configurable<std::vector<float>> cpaMin{"cpaMin", {0.9f}, "Minimum cosine of pointing angle"};
  o2::framework::Configurable<std::vector<float>> decayLengthMin{"decayLengthMin", {0.02f}, "Minimum decay length (cm)"};
  o2::framework::Configurable<std::vector<float>> impactParameterProductMax{"impactParameterProductMax", {0.f}, "Maximum product of prong impact parameters d0*d0 (cm^2)"};
  o2::framework::Configurable<std::vector<float>> cosThetaStarMax{"cosThetaStarMax", {1.f}, "Maximum |cos(theta*)| of the decay"};
};

// derived selection bits for Lc
struct ConfLcBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("LcBits");
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, all Lcs are passed through. Bits for all selections are stored."};
  o2::framework::Configurable<std::vector<float>> cpaMin{"cpaMin", {0.9f}, "Minimum cosine of pointing angle"};
  o2::framework::Configurable<std::vector<float>> decayLengthMin{"decayLengthMin", {0.02f}, "Minimum decay length (cm)"};
  o2::framework::Configurable<std::vector<float>> chi2PcaMax{"chi2PcaMax", {1.f}, "Maximum sum of distances of the secondary vertex to its prongs (cm)"};
  o2::framework::Configurable<std::vector<float>> impactParameterProngSqSumMin{"impactParameterProngSqSumMin", {0.f}, "Minimum sum of squared prong impact parameters (cm^2)"};
};

// base selection for analysis task for D0s
template <auto& Prefix>
struct ConfD0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", o2::constants::physics::Pdg::kD0, "PDG code (D0)"};
  o2::framework::Configurable<int> sign{"sign", 0, "Particle sign (+1: D0; -1: D0bar; 0: both)"};
  o2::framework::Configurable<float> ptMin{"ptMin", 1.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 3.f, "Maximum pT"};
  // acceptance is applied as a rapidity cut in the builder; the eta/phi windows
  // are kept open and exist only to satisfy macro
  o2::framework::Configurable<float> etaMin{"etaMin", -0.8f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.8f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  // signal region; side-bands remain available in the derived data
  o2::framework::Configurable<float> massMin{"massMin", 1.81f, "Minimum invariant mass for D0"};
  o2::framework::Configurable<float> massMax{"massMax", 1.922f, "Maximum invariant mass for D0"};
  o2::framework::Configurable<datatypes::CharmHadronMaskType> mask{"mask", 0, "Bitmask for D0 selection"};
};

// base selection for analysis task for Lcs
template <auto& Prefix>
struct ConfLcSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", o2::constants::physics::Pdg::kLambdaCPlus, "PDG code (Lc)"};
  o2::framework::Configurable<int> sign{"sign", 0, "Particle sign (+1: Lc; -1: Lcbar; 0: both)"};
  o2::framework::Configurable<float> ptMin{"ptMin", 1.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 3.f, "Maximum pT"};
  // acceptance is applied as a rapidity cut in the builder; the eta/phi windows
  // are kept open and exist only to satisfy the partition macro
  o2::framework::Configurable<float> etaMin{"etaMin", -0.8f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.8f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  // signal region; side-bands remain available in the derived data
  o2::framework::Configurable<float> massMin{"massMin", 2.24f, "Minimum invariant mass for Lc"};
  o2::framework::Configurable<float> massMax{"massMax", 2.33f, "Maximum invariant mass for Lc"};
  o2::framework::Configurable<datatypes::CharmHadronMaskType> mask{"mask", 0, "Bitmask for Lc selection"};
};

constexpr const char PrefixD0Selection1[] = "D0Selection1";
constexpr const char PrefixD0Selection2[] = "D0Selection2";
using ConfD0Selection1 = ConfD0Selection<PrefixD0Selection1>;
using ConfD0Selection2 = ConfD0Selection<PrefixD0Selection2>;

constexpr const char PrefixLcSelection1[] = "LcSelection1";
constexpr const char PrefixLcSelection2[] = "LcSelection2";
using ConfLcSelection1 = ConfLcSelection<PrefixLcSelection1>;
using ConfLcSelection2 = ConfLcSelection<PrefixLcSelection2>;

/// The different selections for charm hadrons.
/// The enum is shared by all species in the family; which bits are actually
/// registered is decided per species in CharmHadronSelection::configure.
enum CharmHadronSels {
  // topological selections
  kCpaMin,         ///< Min. CPA (cosine pointing angle)
  kDecayLengthMin, ///< Min. decay length

  // topological selections, 2-prong only (D0)
  kImpactParameterProductMax, ///< Max. product of prong impact parameters (d0*d0)
  kCosThetaStarMax,           ///< Max. |cos(theta*)| of the decay

  // topological selections, 3-prong only (Lc)
  kChi2PcaMax,                   ///< Max. sum of distances of the secondary vertex to its prongs
  kImpactParameterProngSqSumMin, ///< Min. sum of squared prong impact parameters

  kCharmHadronSelsMax
};

constexpr char D0SelHistName[] = "hD0Selection";
constexpr char D0barSelHistName[] = "hD0barSelection";
constexpr char CharmHadronSelsName[] = "Charm hadron selection object";
constexpr char LcSelHistName[] = "hLcSelection";
constexpr char LcBarSelHistName[] = "hLcBarSelection";
const std::unordered_map<CharmHadronSels, std::string> charmHadronSelectionNames = {
  {kCpaMin, "Min. CPA (cosine pointing angle)"},
  {kDecayLengthMin, "Min. decay length"},

  {kImpactParameterProductMax, "Max. product of prong impact parameters (d0*d0)"},
  {kCosThetaStarMax, "Max. |cos(theta*)| of the decay"},

  {kChi2PcaMax, "Max. sum of distances of the secondary vertex to its prongs (cm)"},
  {kImpactParameterProngSqSumMin, "Min. sum of squared prong impact parameters (cm^2)"}};

/// enum for all charm hadron filters (loose kinematic pre-selection, applied before the bit selections).
enum CharmHadronFilters {
  kDecayChannel,
  kHypothesis,
  kPtMin,
  kPtMax,
  kEtaMin,
  kEtaMax,
  kPhiMin,
  kPhiMax,
  kYMin,
  kYMax,
  kMassMin,
  kMassMax,
  kRejectAmbiguous,
  kCharmHadronFiltersMax
};

constexpr char D0FilterHistName[] = "hD0Filters";
constexpr char D0barFilterHistName[] = "hD0barFilters";
constexpr char LcFilterHistName[] = "hLcFilters";
constexpr char LcBarFilterHistName[] = "hLcBarFilters";
const std::unordered_map<CharmHadronFilters, std::string> charmHadronFilterNames = {
  {kDecayChannel, "Decay channel"},
  {kHypothesis, "Mass hypothesis"},
  {kPtMin, "Minimum pT"},
  {kPtMax, "Maximum pT"},
  {kEtaMin, "Minimum eta"},
  {kEtaMax, "Maximum eta"},
  {kPhiMin, "Minimum phi"},
  {kPhiMax, "Maximum phi"},
  {kYMin, "Minimum rapidity"},
  {kYMax, "Maximum rapidity"},
  {kMassMin, "Minimum invariant mass"},
  {kMassMax, "Maximum invariant mass"},
  {kRejectAmbiguous, "Reject ambiguous mass hypothesis"}};

template <modes::CharmHadron hadronType>
constexpr bool isThreeProng()
{
  return modes::isEqual(hadronType, modes::CharmHadron::kLc) ||
         modes::isEqual(hadronType, modes::CharmHadron::kLcBar);
}

template <modes::CharmHadron hadronType>
constexpr bool isParticle()
{
  return modes::isEqual(hadronType, modes::CharmHadron::kD0) ||
         modes::isEqual(hadronType, modes::CharmHadron::kLc);
}

template <modes::CharmHadron hadronType, auto& SelectionHistName, auto& FilterHistName>
class CharmHadronSelection : public baseselection::BaseSelection<float, o2::analysis::femto::datatypes::CharmHadronMaskType, kCharmHadronSelsMax>
{
 public:
  CharmHadronSelection() = default;
  ~CharmHadronSelection() override = default;

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
    mMassMin = filter.massMin.value;
    mMassMax = filter.massMax.value;
    mUseYCut = filter.useYCut.value;
    mYMin = filter.yMin.value;
    mYMax = filter.yMax.value;
    mRejectAmbiguousHypothesis = filter.rejectAmbiguousHypothesis.value;

    this->addSelection(kCpaMin, charmHadronSelectionNames.at(kCpaMin), config.cpaMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kDecayLengthMin, charmHadronSelectionNames.at(kDecayLengthMin), config.decayLengthMin.value, limits::kLowerLimit, true, true, false);

    if constexpr (isThreeProng<hadronType>()) {
      this->addSelection(kChi2PcaMax, charmHadronSelectionNames.at(kChi2PcaMax), config.chi2PcaMax.value, limits::kUpperLimit, true, true, false);
      this->addSelection(kImpactParameterProngSqSumMin, charmHadronSelectionNames.at(kImpactParameterProngSqSumMin), config.impactParameterProngSqSumMin.value, limits::kLowerLimit, true, true, false);
    } else {
      this->addSelection(kImpactParameterProductMax, charmHadronSelectionNames.at(kImpactParameterProductMax), config.impactParameterProductMax.value, limits::kUpperLimit, true, true, false);
      this->addSelection(kCosThetaStarMax, charmHadronSelectionNames.at(kCosThetaStarMax), config.cosThetaStarMax.value, limits::kAbsUpperLimit, true, true, false);
    }

    this->template setupSelectionHistogram<SelectionHistName>(registry);
    this->template setupFilterHistogram<FilterHistName>(
      registry,
      {
        {charmHadronFilterNames.at(kDecayChannel), 1},
        {charmHadronFilterNames.at(kHypothesis), 1},
        {charmHadronFilterNames.at(kPtMin), mPtMin},
        {charmHadronFilterNames.at(kPtMax), mPtMax},
        {charmHadronFilterNames.at(kEtaMin), mEtaMin},
        {charmHadronFilterNames.at(kEtaMax), mEtaMax},
        {charmHadronFilterNames.at(kPhiMin), mPhiMin},
        {charmHadronFilterNames.at(kPhiMax), mPhiMax},
        {charmHadronFilterNames.at(kYMin), mYMin},
        {charmHadronFilterNames.at(kYMax), mYMax},
        {charmHadronFilterNames.at(kMassMin), mMassMin},
        {charmHadronFilterNames.at(kMassMax), mMassMax},
        {charmHadronFilterNames.at(kRejectAmbiguous), static_cast<float>(mRejectAmbiguousHypothesis)},
      });
  }

  template <typename T1>
  void applySelections(T1 const& candidate)
  {
    this->reset();
    this->evaluateObservable(kCpaMin, candidate.cpa());
    this->evaluateObservable(kDecayLengthMin, candidate.decayLength());

    if constexpr (isThreeProng<hadronType>()) {
      this->evaluateObservable(kChi2PcaMax, candidate.chi2PCA());
      this->evaluateObservable(kImpactParameterProngSqSumMin, candidate.impactParameterProngSqSum());
    } else {
      this->evaluateObservable(kImpactParameterProductMax, candidate.impactParameter0() * candidate.impactParameter1());
      this->evaluateObservable(kCosThetaStarMax, mHfHelper.cosThetaStarD0(candidate));
    }
    this->template assembleBitmask<SelectionHistName>();
  }

  template <typename T>
  bool checkFilters(const T& candidate) const
  {
    bool pass = true;
    bool p = false;

    // decay channel
    if constexpr (isThreeProng<hadronType>()) {
      p = (candidate.hfflag() & (1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi)) != 0;
    } else {
      p = (candidate.hfflag() & (1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) != 0;
    }
    this->template fillFilter<FilterHistName>(kDecayChannel, p);
    pass &= p;

    if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
      p = candidate.isSelD0();
    } else if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0Bar)) {
      p = candidate.isSelD0bar();
    } else { // kLc / kLcBar
      p = candidate.isSelLcToPKPi() || candidate.isSelLcToPiKP();
    }
    this->template fillFilter<FilterHistName>(kHypothesis, p);
    pass &= p;

    bool competing = false;
    if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
      competing = candidate.isSelD0bar();
    } else if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0Bar)) {
      competing = candidate.isSelD0();
    } else { // Lc / Lcbar: ambiguous means both mass hypotheses passed
      competing = candidate.isSelLcToPKPi() && candidate.isSelLcToPiKP();
    }
    p = !mRejectAmbiguousHypothesis || !competing;
    this->template fillFilter<FilterHistName>(kRejectAmbiguous, p);
    pass &= p;

    p = candidate.pt() > mPtMin;
    this->template fillFilter<FilterHistName>(kPtMin, p);
    pass &= p;

    p = candidate.pt() < mPtMax;
    this->template fillFilter<FilterHistName>(kPtMax, p);
    pass &= p;

    p = candidate.eta() > mEtaMin;
    this->template fillFilter<FilterHistName>(kEtaMin, p);
    pass &= p;

    p = candidate.eta() < mEtaMax;
    this->template fillFilter<FilterHistName>(kEtaMax, p);
    pass &= p;

    p = candidate.phi() > mPhiMin;
    this->template fillFilter<FilterHistName>(kPhiMin, p);
    pass &= p;

    p = candidate.phi() < mPhiMax;
    this->template fillFilter<FilterHistName>(kPhiMax, p);
    pass &= p;

    // rapidity
    if (mUseYCut) {
      float y = 0.f;
      if constexpr (isThreeProng<hadronType>()) {
        y = mHfHelper.yLc(candidate);
      } else {
        y = mHfHelper.yD0(candidate);
      }
      p = y > mYMin;
      this->template fillFilter<FilterHistName>(kYMin, p);
      pass &= p;

      p = y < mYMax;
      this->template fillFilter<FilterHistName>(kYMax, p);
      pass &= p;
    }

    this->template fillFilterSummary<FilterHistName>(pass);
    return this->isPassThrough() || pass;
  }

  [[nodiscard]] float getMassMin() const
  {
    return mMassMin;
  }

  [[nodiscard]] float getMassMax() const
  {
    return mMassMax;
  }

 private:
  HfHelper mHfHelper;
  float mPtMin = 0.f;
  float mPtMax = 24.f;
  float mEtaMin = -0.8f;
  float mEtaMax = 0.8f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;
  float mMassMin = 1.7f;
  float mMassMax = 2.0f;
  bool mUseYCut = true;
  float mYMin = -0.8f;
  float mYMax = 0.8f;
  bool mRejectAmbiguousHypothesis = false;
};

// tables produced by the charm hadron builder: kinematics, bitmask and QA
struct CharmHadronBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FD0s> producedD0s;
  o2::framework::Produces<o2::aod::FD0Masks> producedD0Masks;
  o2::framework::Produces<o2::aod::FD0Extras> producedD0Extras;
  o2::framework::Produces<o2::aod::FLcs> producedLcs;
  o2::framework::Produces<o2::aod::FLcMasks> producedLcMasks;
  o2::framework::Produces<o2::aod::FLcExtras> producedLcExtras;
};

// per-table produce switches (-1: auto = produce only if a downstream device subscribes; 0 off; 1 on)
struct ConfCharmHadronTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CharmHadronTables");
  o2::framework::Configurable<int> produceD0s{"produceD0s", -1, "Produce D0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceD0Masks{"produceD0Masks", -1, "Produce D0Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceD0Extras{"produceD0Extras", -1, "Produce D0Extras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLcs{"produceLcs", -1, "Produce Lcs (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLcMasks{"produceLcMasks", -1, "Produce LcMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLcExtras{"produceLcExtras", -1, "Produce LcExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::CharmHadron hadronType, auto& SelectionHistName, auto& FilterHistName>
class CharmHadronBuilder
{
 public:
  CharmHadronBuilder() = default;
  ~CharmHadronBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext)
  {
    if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0) || modes::isEqual(hadronType, modes::CharmHadron::kD0Bar)) {
      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
        LOG(info) << "Initialize femto D0 builder...";
      }
      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0Bar)) {
        LOG(info) << "Initialize femto D0bar builder...";
      }
      mProduceD0s = utils::enableTable("FD0s_001", table.produceD0s.value, initContext);
      mProduceD0Masks = utils::enableTable("FD0Masks_001", table.produceD0Masks.value, initContext);
      mProduceD0Extras = utils::enableTable("FD0Extras_001", table.produceD0Extras.value, initContext);
    }
    if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kLc) || modes::isEqual(hadronType, modes::CharmHadron::kLcBar)) {
      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kLc)) {
        LOG(info) << "Initialize femto Lc builder...";
      }
      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kLcBar)) {
        LOG(info) << "Initialize femto Lcbar builder...";
      }
      mProduceLcs = utils::enableTable("FLcs_001", table.produceLcs.value, initContext);
      mProduceLcMasks = utils::enableTable("FLcMasks_001", table.produceLcMasks.value, initContext);
      mProduceLcExtras = utils::enableTable("FLcExtras_001", table.produceLcExtras.value, initContext);
    }
    if (mProduceD0s || mProduceD0Masks || mProduceD0Extras ||
        mProduceLcs || mProduceLcMasks || mProduceLcExtras) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mCharmHadronSelection.configure(registry, config, filter);
    mCharmHadronSelection.printSelections(CharmHadronSelsName);
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3>
  void fillD0Tables(T1& collisionProducts, T2& charmHadronProducts, T3 const& candidate,
                    float signedPt, float mass, int64_t posDauIndex, int64_t negDauIndex)
  {
    if (mProduceD0s) {
      charmHadronProducts.producedD0s(collisionProducts.producedCollision.lastIndex(),
                                      signedPt,
                                      candidate.eta(),
                                      candidate.phi(),
                                      mass,
                                      posDauIndex,
                                      negDauIndex);
    }
    if (mProduceD0Masks) {
      charmHadronProducts.producedD0Masks(mCharmHadronSelection.getBitmask());
    }
    if (mProduceD0Extras) {
      charmHadronProducts.producedD0Extras(
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.impactParameter0() * candidate.impactParameter1(),
        mHfHelper.cosThetaStarD0(candidate),
        candidate.mlProbD0().size() < NSizeMLScore ? -1.f : candidate.mlProbD0()[0],
        candidate.mlProbD0().size() < NSizeMLScore ? -1.f : candidate.mlProbD0()[1],
        candidate.mlProbD0().size() < NSizeMLScore ? -1.f : candidate.mlProbD0()[2],
        candidate.mlProbD0bar().size() < NSizeMLScore ? -1.f : candidate.mlProbD0bar()[0],
        candidate.mlProbD0bar().size() < NSizeMLScore ? -1.f : candidate.mlProbD0bar()[1],
        candidate.mlProbD0bar().size() < NSizeMLScore ? -1.f : candidate.mlProbD0bar()[2],
        static_cast<int8_t>(candidate.isSelD0()),
        static_cast<int8_t>(candidate.isSelD0bar()));
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillLcTables(T1& collisionProducts, T2& charmHadronProducts, T3 const& candidate,
                    float signedPt, float mass, float massCompeting,
                    int64_t protonDauIndex, int64_t kaonDauIndex, int64_t pionDauIndex)
  {
    if (mProduceLcs) {
      charmHadronProducts.producedLcs(collisionProducts.producedCollision.lastIndex(),
                                      signedPt,
                                      candidate.eta(),
                                      candidate.phi(),
                                      mass,
                                      protonDauIndex,
                                      kaonDauIndex,
                                      pionDauIndex);
    }
    if (mProduceLcMasks) {
      charmHadronProducts.producedLcMasks(mCharmHadronSelection.getBitmask());
    }
    if (mProduceLcExtras) {
      charmHadronProducts.producedLcExtras(
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.chi2PCA(),
        candidate.impactParameterProngSqSum(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        massCompeting,
        candidate.mlProbLcToPKPi().size() < NSizeMLScore ? -1.f : candidate.mlProbLcToPKPi()[0],
        candidate.mlProbLcToPKPi().size() < NSizeMLScore ? -1.f : candidate.mlProbLcToPKPi()[1],
        candidate.mlProbLcToPKPi().size() < NSizeMLScore ? -1.f : candidate.mlProbLcToPKPi()[2],
        candidate.mlProbLcToPiKP().size() < NSizeMLScore ? -1.f : candidate.mlProbLcToPiKP()[0],
        candidate.mlProbLcToPiKP().size() < NSizeMLScore ? -1.f : candidate.mlProbLcToPiKP()[1],
        candidate.mlProbLcToPiKP().size() < NSizeMLScore ? -1.f : candidate.mlProbLcToPiKP()[2],
        static_cast<int8_t>(candidate.isSelLcToPKPi()),
        static_cast<int8_t>(candidate.isSelLcToPiKP()));
    }
  }

  /// pass-through accepts every candidate once, under the pKPi prong mapping
  template <bool isPKPi, typename T>
  bool acceptsHypothesis(T const& candidate) const
  {
    if (mCharmHadronSelection.isPassThrough()) {
      return isPKPi;
    }
    return isPKPi ? candidate.isSelLcToPKPi() : candidate.isSelLcToPiKP();
  }

  /// a candidate passing both hypotheses is stored twice, once per hypothesis
  template <bool isPKPi, modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillLcHypothesis(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts,
                        T5& charmHadronProducts, T6 const& candidate, T7& trackBuilder,
                        T8 const& prong0, T8 const& prong1, T8 const& prong2)
  {
    if (!this->acceptsHypothesis<isPKPi>(candidate)) {
      return;
    }

    // remap the prongs onto the accepted hypothesis so that the proton is always first
    auto const& protonProng = isPKPi ? prong0 : prong2;
    auto const& pionProng = isPKPi ? prong2 : prong0;
    float const mass = isPKPi ? mHfHelper.invMassLcToPKPi(candidate) : mHfHelper.invMassLcToPiKP(candidate);
    float const massCompeting = isPKPi ? mHfHelper.invMassLcToPiKP(candidate) : mHfHelper.invMassLcToPKPi(candidate);

    collisionBuilder.template fillCollision<system>(collisionProducts, col);
    int64_t const protonDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kCharmDaughter>(protonProng, trackProducts, collisionBuilder);
    int64_t const kaonDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kCharmDaughter>(prong1, trackProducts, collisionBuilder);
    int64_t const pionDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kCharmDaughter>(pionProng, trackProducts, collisionBuilder);

    float const signedPt = isParticle<hadronType>() ? candidate.pt() : -candidate.pt();
    this->fillLcTables(collisionProducts, charmHadronProducts, candidate, signedPt, mass, massCompeting, protonDauIndex, kaonDauIndex, pionDauIndex);
  }

  /// Monte Carlo counterpart of fillLcHypothesis. The label is written inside this method,
  /// so FLcs and FLcLabels stay in lockstep also for candidates passing both hypotheses.
  template <bool isPKPi, modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13>
  void fillMcLcHypothesis(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts,
                          T6& charmHadronProducts, T7 const& candidate, T8 const& tracks, T9& trackBuilder,
                          T10 const& mcParticles, T11& mcBuilder, T12& mcProducts,
                          T13 const& prong0, T13 const& prong1, T13 const& prong2)
  {
    if (!this->acceptsHypothesis<isPKPi>(candidate)) {
      return;
    }

    // remap the prongs onto the accepted hypothesis so that the proton is always first
    auto const& protonProng = isPKPi ? prong0 : prong2;
    auto const& pionProng = isPKPi ? prong2 : prong0;
    float const mass = isPKPi ? mHfHelper.invMassLcToPKPi(candidate) : mHfHelper.invMassLcToPiKP(candidate);
    float const massCompeting = isPKPi ? mHfHelper.invMassLcToPiKP(candidate) : mHfHelper.invMassLcToPKPi(candidate);

    collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);

    int64_t const protonDauIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCharmDaughter>(protonProng, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);
    int64_t const kaonDauIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCharmDaughter>(prong1, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);
    int64_t const pionDauIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCharmDaughter>(pionProng, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

    float const signedPt = isParticle<hadronType>() ? candidate.pt() : -candidate.pt();
    this->fillLcTables(collisionProducts, charmHadronProducts, candidate, signedPt, mass, massCompeting, protonDauIndex, kaonDauIndex, pionDauIndex);
    mcBuilder.template fillMcLcWithLabel<system>(candidate, tracks, mcParticles, mcCols, mcProducts);
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillCharmHadrons(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts,
                        T5& charmHadronProducts, T6 const& candidates, T7 const& /*tracks*/, T8& trackBuilder)
  {
    if (!mFillAnyTable) {
      return;
    }
    for (const auto& candidate : candidates) {
      if (!mCharmHadronSelection.checkFilters(candidate)) {
        continue;
      }
      mCharmHadronSelection.applySelections(candidate);
      if (!mCharmHadronSelection.passesAllRequiredSelections()) {
        continue;
      }

      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0) || modes::isEqual(hadronType, modes::CharmHadron::kD0Bar)) {
        auto prong0 = candidate.template prong0_as<T7>();
        auto prong1 = candidate.template prong1_as<T7>();

        collisionBuilder.template fillCollision<system>(collisionProducts, col);

        int64_t posDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kCharmDaughter>(prong0, trackProducts, collisionBuilder);
        int64_t negDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kCharmDaughter>(prong1, trackProducts, collisionBuilder);

        if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
          this->fillD0Tables(collisionProducts, charmHadronProducts, candidate, candidate.pt(), mHfHelper.invMassD0ToPiK(candidate), posDauIndex, negDauIndex);
        } else {
          this->fillD0Tables(collisionProducts, charmHadronProducts, candidate, -candidate.pt(), mHfHelper.invMassD0barToKPi(candidate), posDauIndex, negDauIndex);
        }
      }

      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kLc) || modes::isEqual(hadronType, modes::CharmHadron::kLcBar)) {
        auto prong0 = candidate.template prong0_as<T7>();
        auto prong1 = candidate.template prong1_as<T7>();
        auto prong2 = candidate.template prong2_as<T7>();

        // the charge of prong 0 is the charge of the Lc, so the species is decided here
        if ((prong0.sign() > 0) != isParticle<hadronType>()) {
          continue;
        }

        this->fillLcHypothesis<true, system>(col, collisionBuilder, collisionProducts, trackProducts, charmHadronProducts, candidate, trackBuilder, prong0, prong1, prong2);
        this->fillLcHypothesis<false, system>(col, collisionBuilder, collisionProducts, trackProducts, charmHadronProducts, candidate, trackBuilder, prong0, prong1, prong2);
      }
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
  void fillMcCharmHadrons(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts,
                          T6& charmHadronProducts, T7 const& candidates, T8 const& tracks, T9& trackBuilder, T10 const& mcParticles, T11& mcBuilder, T12& mcProducts)
  {
    if (!mFillAnyTable) {
      return;
    }
    for (const auto& candidate : candidates) {
      if (!mCharmHadronSelection.checkFilters(candidate)) {
        continue;
      }

      mCharmHadronSelection.applySelections(candidate);
      if (!mCharmHadronSelection.passesAllRequiredSelections()) {
        continue;
      }

      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0) || modes::isEqual(hadronType, modes::CharmHadron::kD0Bar)) {
        auto prong0 = candidate.template prong0_as<T8>();
        auto prong1 = candidate.template prong1_as<T8>();

        collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);

        int64_t posDauIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCharmDaughter>(prong0, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);
        int64_t negDauIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCharmDaughter>(prong1, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

        if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
          this->fillD0Tables(collisionProducts, charmHadronProducts, candidate, candidate.pt(), mHfHelper.invMassD0ToPiK(candidate), posDauIndex, negDauIndex);
        } else {
          this->fillD0Tables(collisionProducts, charmHadronProducts, candidate, -candidate.pt(), mHfHelper.invMassD0barToKPi(candidate), posDauIndex, negDauIndex);
        }
        mcBuilder.template fillMcD0WithLabel<system>(candidate, tracks, mcParticles, mcCols, mcProducts);
      }

      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kLc) || modes::isEqual(hadronType, modes::CharmHadron::kLcBar)) {
        auto prong0 = candidate.template prong0_as<T8>();
        auto prong1 = candidate.template prong1_as<T8>();
        auto prong2 = candidate.template prong2_as<T8>();

        // the charge of prong 0 is the charge of the Lc, so the species is decided here
        if ((prong0.sign() > 0) != isParticle<hadronType>()) {
          continue;
        }

        this->fillMcLcHypothesis<true, system>(col, collisionBuilder, collisionProducts, mcCols, trackProducts, charmHadronProducts, candidate, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts, prong0, prong1, prong2);
        this->fillMcLcHypothesis<false, system>(col, collisionBuilder, collisionProducts, mcCols, trackProducts, charmHadronProducts, candidate, tracks, trackBuilder, mcParticles, mcBuilder, mcProducts, prong0, prong1, prong2);
      }
    }
  }

  [[nodiscard]] bool fillAnyTable() const { return mFillAnyTable; }
  [[nodiscard]] bool isPassThrough() const { return mCharmHadronSelection.isPassThrough(); }

 private:
  CharmHadronSelection<hadronType, SelectionHistName, FilterHistName> mCharmHadronSelection;
  HfHelper mHfHelper;

  bool mProduceD0s = false;
  bool mProduceD0Masks = false;
  bool mProduceD0Extras = false;
  bool mProduceLcs = false;
  bool mProduceLcMasks = false;
  bool mProduceLcExtras = false;
  bool mFillAnyTable = false;
};
} // namespace o2::analysis::femto::charmhadronbuilder

#endif // PWGCF_FEMTO_CORE_CHARMHADRONBUILDER_H_
