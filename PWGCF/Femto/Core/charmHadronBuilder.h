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
struct ConfD0Filters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("D0Filters");
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
  o2::framework::Configurable<float> massMin{"massMin", 1.7f, "Minimum invariant mass for D0"};
  o2::framework::Configurable<float> massMax{"massMax", 2.0f, "Maximum invariant mass for D0"};
  o2::framework::Configurable<bool> rejectAmbiguousHypothesis{"rejectAmbiguousHypothesis", false, "Reject candidates selected under both D0 and D0bar hypotheses"};
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

// base selection for analysis task for D0s
template <auto& Prefix>
struct ConfD0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", o2::constants::physics::Pdg::kD0, "PDG code (D0)"};
  o2::framework::Configurable<int> sign{"sign", 0, "Particle sign (+1: D0; -1: D0bar; 0: both)"};
  o2::framework::Configurable<float> ptMin{"ptMin", 1.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 3.f, "Maximum pT"};
  // acceptance is applied as a rapidity cut in the builder; the eta/phi windows
  // are kept open and exist only to satisfy MAKE_D0_PARTITION
  o2::framework::Configurable<float> etaMin{"etaMin", -0.8f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.8f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  // signal region; side-bands remain available in the derived data
  o2::framework::Configurable<float> massMin{"massMin", 1.81f, "Minimum invariant mass for D0"};
  o2::framework::Configurable<float> massMax{"massMax", 1.922f, "Maximum invariant mass for D0"};
  o2::framework::Configurable<datatypes::CharmHadronMaskType> mask{"mask", 0, "Bitmask for D0 selection"};
};

constexpr const char PrefixD0Selection1[] = "D0Selection1";
constexpr const char PrefixD0Selection2[] = "D0Selection2";
using ConfD0Selection1 = ConfD0Selection<PrefixD0Selection1>;
using ConfD0Selection2 = ConfD0Selection<PrefixD0Selection2>;

/// The different selections for D0s
enum D0Sels {
  // topological selections
  kCpaMin,                    ///< Min. CPA (cosine pointing angle)
  kDecayLengthMin,            ///< Min. decay length
  kImpactParameterProductMax, ///< Max. product of prong impact parameters (d0*d0)
  kCosThetaStarMax,           ///< Max. |cos(theta*)| of the decay

  kD0SelsMax
};

constexpr char D0SelHistName[] = "hD0Selection";
constexpr char D0barSelHistName[] = "hD0barSelection";
constexpr char D0SelsName[] = "D0 selection object";
const std::unordered_map<D0Sels, std::string> d0SelectionNames = {
  {kCpaMin, "Min. CPA (cosine pointing angle)"},
  {kDecayLengthMin, "Min. decay length"},
  {kImpactParameterProductMax, "Max. product of prong impact parameters (d0*d0)"},
  {kCosThetaStarMax, "Max. |cos(theta*)| of the decay"}};

/// enum for all D0 filters (loose kinematic pre-selection, applied before the bit selections)
enum D0Filters {
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
  kD0FiltersMax
};

constexpr char D0FilterHistName[] = "hD0Filters";
constexpr char D0barFilterHistName[] = "hD0barFilters";
const std::unordered_map<D0Filters, std::string> d0FilterNames = {
  {kDecayChannel, "D0 -> K pi decay channel"},
  {kHypothesis, "D0/D0bar hypothesis"},
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
  {kRejectAmbiguous, "Reject ambiguous D0/D0bar hypothesis"}};

template <modes::CharmHadron hadronType, auto& SelectionHistName, auto& FilterHistName>
class D0Selection : public baseselection::BaseSelection<float, o2::analysis::femto::datatypes::CharmHadronMaskType, kD0SelsMax>
{
 public:
  D0Selection() = default;
  ~D0Selection() override = default;

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

    this->addSelection(kCpaMin, d0SelectionNames.at(kCpaMin), config.cpaMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kDecayLengthMin, d0SelectionNames.at(kDecayLengthMin), config.decayLengthMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kImpactParameterProductMax, d0SelectionNames.at(kImpactParameterProductMax), config.impactParameterProductMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kCosThetaStarMax, d0SelectionNames.at(kCosThetaStarMax), config.cosThetaStarMax.value, limits::kAbsUpperLimit, true, true, false);

    this->template setupSelectionHistogram<SelectionHistName>(registry);
    this->template setupFilterHistogram<FilterHistName>(
      registry,
      {
        {d0FilterNames.at(kDecayChannel), 1},
        {d0FilterNames.at(kHypothesis), 1},
        {d0FilterNames.at(kPtMin), mPtMin},
        {d0FilterNames.at(kPtMax), mPtMax},
        {d0FilterNames.at(kEtaMin), mEtaMin},
        {d0FilterNames.at(kEtaMax), mEtaMax},
        {d0FilterNames.at(kPhiMin), mPhiMin},
        {d0FilterNames.at(kPhiMax), mPhiMax},
        {d0FilterNames.at(kYMin), mYMin},
        {d0FilterNames.at(kYMax), mYMax},
        {d0FilterNames.at(kMassMin), mMassMin},
        {d0FilterNames.at(kMassMax), mMassMax},
        {d0FilterNames.at(kRejectAmbiguous), static_cast<float>(mRejectAmbiguousHypothesis)},
      });
  }

  template <typename T1>
  void applySelections(T1 const& d0candidate)
  {
    this->reset();
    this->evaluateObservable(kCpaMin, d0candidate.cpa());
    this->evaluateObservable(kDecayLengthMin, d0candidate.decayLength());
    this->evaluateObservable(kImpactParameterProductMax, d0candidate.impactParameter0() * d0candidate.impactParameter1());
    this->evaluateObservable(kCosThetaStarMax, mHfHelper.cosThetaStarD0(d0candidate));
    this->template assembleBitmask<SelectionHistName>();
  }

  template <typename T>
  bool checkFilters(const T& d0candidate) const
  {
    bool pass = true;
    bool p = false;

    // decay channel
    p = (d0candidate.hfflag() & (1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK)) != 0;
    this->template fillFilter<FilterHistName>(kDecayChannel, p);
    pass &= p;

    if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
      p = d0candidate.isSelD0();
    } else {
      p = d0candidate.isSelD0bar();
    }

    this->template fillFilter<FilterHistName>(kHypothesis, p);
    pass &= p;

    bool competing = false;
    if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
      competing = d0candidate.isSelD0bar();
    } else {
      competing = d0candidate.isSelD0();
    }
    p = !mRejectAmbiguousHypothesis || !competing;
    this->template fillFilter<FilterHistName>(kRejectAmbiguous, p);
    pass &= p;

    p = d0candidate.pt() > mPtMin;
    this->template fillFilter<FilterHistName>(kPtMin, p);
    pass &= p;

    p = d0candidate.pt() < mPtMax;
    this->template fillFilter<FilterHistName>(kPtMax, p);
    pass &= p;

    p = d0candidate.eta() > mEtaMin;
    this->template fillFilter<FilterHistName>(kEtaMin, p);
    pass &= p;

    p = d0candidate.eta() < mEtaMax;
    this->template fillFilter<FilterHistName>(kEtaMax, p);
    pass &= p;

    p = d0candidate.phi() > mPhiMin;
    this->template fillFilter<FilterHistName>(kPhiMin, p);
    pass &= p;

    p = d0candidate.phi() < mPhiMax;
    this->template fillFilter<FilterHistName>(kPhiMax, p);
    pass &= p;

    // rapidity
    if (mUseYCut) {
      const float y = mHfHelper.yD0(d0candidate);
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

// tables produced by the D0 builder: kinematics, bitmask and QA
struct CharmHadronBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FD0s> producedD0s;
  o2::framework::Produces<o2::aod::FD0Masks> producedD0Masks;
  o2::framework::Produces<o2::aod::FD0Extras> producedD0Extras;
};

// per-table produce switches (-1: auto = produce only if a downstream device subscribes; 0 off; 1 on)
struct ConfD0Tables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("D0Tables");
  o2::framework::Configurable<int> produceD0s{"produceD0s", -1, "Produce D0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceD0Masks{"produceD0Masks", -1, "Produce D0Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceD0Extras{"produceD0Extras", -1, "Produce D0Extras (-1: auto; 0 off; 1 on)"};
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
    LOG(info) << "Initialize femto D0 builder...";
    mProduceD0s = utils::enableTable("FD0s_001", table.produceD0s.value, initContext);
    mProduceD0Masks = utils::enableTable("FD0Masks_001", table.produceD0Masks.value, initContext);
    mProduceD0Extras = utils::enableTable("FD0Extras_001", table.produceD0Extras.value, initContext);

    if (mProduceD0s || mProduceD0Masks || mProduceD0Extras) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No D0 tables configured, selection object will not be configured...";
      return;
    }

    mD0Selection.configure(registry, config, filter);
    mD0Selection.printSelections(D0SelsName);
  }

  template <typename T1, typename T2, typename T3>
  void fillD0Tables(T1& collisionProducts, T2& d0Products, T3 const& candidate,
                    float signedPt, float mass, int64_t posDauIndex, int64_t negDauIndex)
  {
    if (mProduceD0s) {
      d0Products.producedD0s(collisionProducts.producedCollision.lastIndex(),
                             signedPt,
                             candidate.eta(),
                             candidate.phi(),
                             mass,
                             posDauIndex,
                             negDauIndex);
    }
    if (mProduceD0Masks) {
      d0Products.producedD0Masks(mD0Selection.getBitmask());
    }
    if (mProduceD0Extras) {
      d0Products.producedD0Extras(
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

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillD0s(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts,
               T5& d0Products, T6 const& candidates, T7 const& /*tracks*/, T8& trackBuilder)
  {
    if (!mFillAnyTable) {
      return;
    }
    int64_t posDauIndex = 0;
    int64_t negDauIndex = 0;
    for (const auto& candidate : candidates) {
      if (!mD0Selection.checkFilters(candidate)) {
        continue;
      }
      mD0Selection.applySelections(candidate);
      if (!mD0Selection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillCollision<system>(collisionProducts, col);

      auto prong0 = candidate.template prong0_as<T7>();
      auto prong1 = candidate.template prong1_as<T7>();
      posDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kCharmDaughter>(prong0, trackProducts, collisionBuilder);
      negDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kCharmDaughter>(prong1, trackProducts, collisionBuilder);

      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
        this->fillD0Tables(collisionProducts, d0Products, candidate, candidate.pt(), mHfHelper.invMassD0ToPiK(candidate), posDauIndex, negDauIndex);
      } else {
        this->fillD0Tables(collisionProducts, d0Products, candidate, -candidate.pt(), mHfHelper.invMassD0barToKPi(candidate), posDauIndex, negDauIndex);
      }
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
  void fillMcD0s(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts,
                 T6& d0Products, T7 const& candidates, T8 const& tracks, T9& trackBuilder, T10 const& mcParticles, T11& mcBuilder, T12& mcProducts)
  {
    if (!mFillAnyTable) {
      return;
    }
    int64_t posDauIndex = 0;
    int64_t negDauIndex = 0;
    for (const auto& candidate : candidates) {
      if (!mD0Selection.checkFilters(candidate)) {
        continue;
      }

      mD0Selection.applySelections(candidate);
      if (!mD0Selection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);

      auto prong0 = candidate.template prong0_as<T8>();
      auto prong1 = candidate.template prong1_as<T8>();
      posDauIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCharmDaughter>(prong0, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);
      negDauIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kCharmDaughter>(prong1, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

      if constexpr (modes::isEqual(hadronType, modes::CharmHadron::kD0)) {
        this->fillD0Tables(collisionProducts, d0Products, candidate, candidate.pt(), mHfHelper.invMassD0ToPiK(candidate), posDauIndex, negDauIndex);
      } else {
        this->fillD0Tables(collisionProducts, d0Products, candidate, -candidate.pt(), mHfHelper.invMassD0barToKPi(candidate), posDauIndex, negDauIndex);
      }
      mcBuilder.template fillMcD0WithLabel<system>(candidate, tracks, mcParticles, mcCols, mcProducts);
    }
  }

  [[nodiscard]] bool fillAnyTable() const { return mFillAnyTable; }
  [[nodiscard]] bool isPassThrough() const { return mD0Selection.isPassThrough(); }

 private:
  D0Selection<hadronType, SelectionHistName, FilterHistName> mD0Selection;
  HfHelper mHfHelper;

  bool mProduceD0s = false;
  bool mProduceD0Masks = false;
  bool mProduceD0Extras = false;
  bool mFillAnyTable = false;
};
} // namespace o2::analysis::femto::charmhadronbuilder

#endif // PWGCF_FEMTO_CORE_CHARMHADRONBUILDER_H_
