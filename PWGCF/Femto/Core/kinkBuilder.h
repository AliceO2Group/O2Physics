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

/// \file kinkBuilder.h
/// \brief kink builder
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#ifndef PWGCF_FEMTO_CORE_KINKBUILDER_H_
#define PWGCF_FEMTO_CORE_KINKBUILDER_H_

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/Core/RecoDecay.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{
namespace kinkbuilder
{

// filters applied in the producer task
struct ConfKinkFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("KinkFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 99.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMinSigma{"massMinSigma", 1.1f, "Minimum mass for Sigma hypothesis"};
  o2::framework::Configurable<float> massMaxSigma{"massMaxSigma", 1.3f, "Maximum mass for Sigma hypothesis"};
  o2::framework::Configurable<float> massMinSigmaPlus{"massMinSigmaPlus", 1.1f, "Minimum mass for SigmaPlus hypothesis"};
  o2::framework::Configurable<float> massMaxSigmaPlus{"massMaxSigmaPlus", 1.3f, "Maximum mass for SigmaPlus hypothesis"};
};

// selections bits for all kinks
#define KINK_DEFAULT_BITS                                                                                                                     \
  o2::framework::Configurable<std::vector<float>> kinkTopoDcaMax{"kinkTopoDcaMax", {2.0f}, "Maximum kink topological DCA"};                   \
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2f}, "Minimum transverse radius (cm)"};                       \
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};                      \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {0.8f}, "Maximum absolute pseudorapidity for daughter track"}; \
  o2::framework::Configurable<std::vector<float>> dauDcaPvMin{"dauDcaPvMin", {0.0f}, "Minimum DCA of daughter from primary vertex (cm)"};     \
  o2::framework::Configurable<std::vector<float>> mothDcaPvMax{"mothDcaPvMax", {1.0f}, "Maximum DCA of mother from primary vertex (cm)"};     \
  o2::framework::Configurable<std::vector<float>> alphaAPMin{"alphaAPMin", {-1.0f}, "Minimum Alpha_AP for Sigma candidates"};                 \
  o2::framework::Configurable<std::vector<float>> alphaAPMax{"alphaAPMax", {0.0f}, "Maximum Alpha_AP for Sigma candidates"};                  \
  o2::framework::Configurable<std::vector<float>> qtAPMin{"qtAPMin", {0.15f}, "Minimum qT_AP for Sigma candidates"};                          \
  o2::framework::Configurable<std::vector<float>> qtAPMax{"qtAPMax", {0.2f}, "Maximum qT_AP for Sigma candidates"};                           \
  o2::framework::Configurable<std::vector<float>> cosPointingAngleMin{"cosPointingAngleMin", {0.0f}, "Minimum cosine of pointing angle"};

// derived selection bits for sigma
struct ConfSigmaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("SigmaBits");
  KINK_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> chaDauTpcPion{"chaDauTpcPion", {5.f}, "Maximum |nsigma_Pion| TPC for charged daughter tracks"};
};

// derived selection bits for sigma plus
struct ConfSigmaPlusBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("SigmaPlusBits");
  KINK_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> chaDauTpcProton{"chaDauTpcProton", {5.f}, "Maximum |nsigma_Proton| TPC for charged daughter tracks"};
  o2::framework::Configurable<std::vector<float>> chaDauTofProton{"chaDauTofProton", {5.f}, "Maximum combined |nsigma_Proton| (TPC+TOF) for charged daughter tracks"};
  o2::framework::Configurable<float> pidThres{"pidThres", 0.75f, "Momentum threshold for using TOF/combined pid for daughter tracks (GeV/c)"};
};

#undef KINK_DEFAULT_BITS

// base selection for analysis task for kinks
#define KINK_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                              \
  o2::framework::Configurable<int> pdgCode{"pdgCode", defaultPdgCode, "Kink PDG code"};                      \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                      \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                    \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                 \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                  \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                   \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};      \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Sigma"}; \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Sigma"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::KinkMaskType> mask{"mask", 0x0, "Bitmask for kink selection"};

// base selection for analysis task for sigmas
template <const char* Prefix>
struct ConfSigmaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_SELECTIONS(1.1, 1.3, 3112)
  o2::framework::Configurable<int> sign{"sign", -1, "Sign of the Sigma mother track (e.g. -1 for Sigma- or +1 for AntiSigma-)"};
};

// base selection for analysis task for sigma plus
template <const char* Prefix>
struct ConfSigmaPlusSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_SELECTIONS(1.1, 1.3, 3222)
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Sigma mother track (e.g. +1 for Sigma+ or -1 for AntiSigma+)"};
};

#undef KINK_DEFAULT_SELECTIONS

constexpr const char PrefixSigmaSelection1[] = "SigmaSelection1";
constexpr const char PrefixSigmaSelection2[] = "SigmaSelection2";
using ConfSigmaSelection1 = ConfSigmaSelection<PrefixSigmaSelection1>;
using ConfSigmaSelection2 = ConfSigmaSelection<PrefixSigmaSelection2>;
constexpr const char PrefixSigmaPlusSelection1[] = "SigmaPlusSelection1";
constexpr const char PrefixSigmaPlusSelection2[] = "SigmaPlusSelection2";
using ConfSigmaPlusSelection1 = ConfSigmaPlusSelection<PrefixSigmaPlusSelection1>;
using ConfSigmaPlusSelection2 = ConfSigmaPlusSelection<PrefixSigmaPlusSelection2>;

/// The different selections for kinks
enum KinkSeles {
  kKinkTopoDcaMax,
  kTransRadMin,
  kTransRadMax,

  kDauAbsEtaMax,
  kDauDcaPvMin,
  kMothDcaPvMax,

  kChaDaughTpcPion,
  kChaDaughTpcProton,
  kChaDaughTofProton,

  kAlphaAPMin,
  kAlphaAPMax,
  kQtAPMin,
  kQtAPMax,
  kCosPointingAngleMin,

  kKinkSelsMax
};

constexpr char SigmaSelHistName[] = "hSigmaSelection";
constexpr char SigmaPlusSelHistName[] = "hSigmaPlusSelection";
const char kinkSelsName[] = "Kink selection object";
const std::unordered_map<KinkSeles, std::string> kinkSelectionNames = {
  {kKinkTopoDcaMax, "kinkTopoDcaMax"},
  {kTransRadMin, "transRadMin"},
  {kTransRadMax, "transRadMax"},
  {kDauAbsEtaMax, "dauAbsEtaMax"},
  {kDauDcaPvMin, "dauDcaPvMin"},
  {kMothDcaPvMax, "mothDcaPvMax"},
  {kChaDaughTpcPion, "chaDauTpcPion"},
  {kChaDaughTpcProton, "chaDauTpcProton"},
  {kChaDaughTofProton, "chaDauTofProton"},
  {kAlphaAPMin, "alphaAPMin"},
  {kAlphaAPMax, "alphaAPMax"},
  {kQtAPMin, "qtAPMin"},
  {kQtAPMax, "qtAPMax"},
  {kCosPointingAngleMin, "cosPointingAngleMin"}};

/// \class KinkCuts
/// \brief Cut class to contain and execute all cuts applied to kinks
template <modes::Kink kinkType, const char* HistName>
class KinkSelection : public BaseSelection<float, o2::aod::femtodatatypes::KinkMaskType, kKinkSelsMax>
{
 public:
  KinkSelection() = default;
  ~KinkSelection() = default;

  template <typename T1, typename T2>
  void configure(o2::framework::HistogramRegistry* registry, T1& config, T2& filter)
  {
    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      mMassSigmaLowerLimit = filter.massMinSigma.value;
      mMassSigmaUpperLimit = filter.massMaxSigma.value;
      // Only add PID selection if we need it - will be checked at runtime
      this->addSelection(kChaDaughTpcPion, kinkSelectionNames.at(kChaDaughTpcPion), config.chaDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
    }

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      mMassSigmaPlusLowerLimit = filter.massMinSigmaPlus.value;
      mMassSigmaPlusUpperLimit = filter.massMaxSigmaPlus.value;
      mPidThreshold = config.pidThres.value;
      this->addSelection(kChaDaughTpcProton, kinkSelectionNames.at(kChaDaughTpcProton), config.chaDauTpcProton.value, limits::kAbsUpperLimit, false, false, true);
      this->addSelection(kChaDaughTofProton, kinkSelectionNames.at(kChaDaughTofProton), config.chaDauTofProton.value, limits::kUpperLimit, false, false, true);
    }

    this->addSelection(kKinkTopoDcaMax, kinkSelectionNames.at(kKinkTopoDcaMax), config.kinkTopoDcaMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kTransRadMin, kinkSelectionNames.at(kTransRadMin), config.transRadMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTransRadMax, kinkSelectionNames.at(kTransRadMax), config.transRadMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kDauAbsEtaMax, kinkSelectionNames.at(kDauAbsEtaMax), config.dauAbsEtaMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kDauDcaPvMin, kinkSelectionNames.at(kDauDcaPvMin), config.dauDcaPvMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kMothDcaPvMax, kinkSelectionNames.at(kMothDcaPvMax), config.mothDcaPvMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kAlphaAPMin, kinkSelectionNames.at(kAlphaAPMin), config.alphaAPMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kAlphaAPMax, kinkSelectionNames.at(kAlphaAPMax), config.alphaAPMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kQtAPMin, kinkSelectionNames.at(kQtAPMin), config.qtAPMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kQtAPMax, kinkSelectionNames.at(kQtAPMax), config.qtAPMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kCosPointingAngleMin, kinkSelectionNames.at(kCosPointingAngleMin), config.cosPointingAngleMin.value, limits::kLowerLimit, true, true, false);

    this->setupContainers<HistName>(registry);
  };

  template <typename T1, typename T2>
  void applySelections(T1 const& kinkCand, T2 const& /*tracks*/)
  {
    this->reset();
    // kink selections
    std::array<float, 3> momMother = {kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
    std::array<float, 3> momDaughter = {kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};

    // Alpha_AP
    std::array<float, 3> momMissing = {momMother[0] - momDaughter[0], momMother[1] - momDaughter[1], momMother[2] - momDaughter[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momDaughter.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    float alphaAP = (lQlP + lQlN != 0.f) ? (lQlP - lQlN) / (lQlP + lQlN) : 0.f;
    this->evaluateObservable(kAlphaAPMin, alphaAP);
    this->evaluateObservable(kAlphaAPMax, alphaAP);

    // qT_AP
    float dp = lQlP;
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momDaughter.begin(), momDaughter.end(), momDaughter.begin(), 0.f);
    float qtAP = std::sqrt(std::max(0.f, p2A - dp * dp / p2V0));
    this->evaluateObservable(kQtAPMin, qtAP);
    this->evaluateObservable(kQtAPMax, qtAP);

    std::array<float, 3> vMother = {kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx()};
    float pMother = std::sqrt(std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f));
    float vMotherNorm = std::sqrt(std::inner_product(vMother.begin(), vMother.end(), vMother.begin(), 0.f));
    float cosPointingAngle = (vMotherNorm > 0.f && pMother > 0.f) ? (std::inner_product(momMother.begin(), momMother.end(), vMother.begin(), 0.f)) / (pMother * vMotherNorm) : 0.f;
    this->evaluateObservable(kCosPointingAngleMin, cosPointingAngle);

    this->evaluateObservable(kKinkTopoDcaMax, kinkCand.dcaKinkTopo());

    // Compute transRadius
    float transRadius = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx());
    this->evaluateObservable(kTransRadMin, transRadius);
    this->evaluateObservable(kTransRadMax, transRadius);

    // Compute daughter eta
    float pxDaug = kinkCand.pxDaug();
    float pyDaug = kinkCand.pyDaug();
    float pzDaug = kinkCand.pzDaug();
    float pDaug = std::sqrt(pxDaug * pxDaug + pyDaug * pyDaug + pzDaug * pzDaug);
    float etaDaug = (pDaug > 0.f) ? 0.5f * std::log((pDaug + pzDaug) / (pDaug - pzDaug)) : 0.f;
    this->evaluateObservable(kDauAbsEtaMax, std::fabs(etaDaug));

    this->evaluateObservable(kDauDcaPvMin, std::abs(kinkCand.dcaDaugPv()));
    this->evaluateObservable(kMothDcaPvMax, std::abs(kinkCand.dcaMothPv()));

    auto chaDaughter = kinkCand.template trackDaug_as<T2>();

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      this->evaluateObservable(kChaDaughTpcPion, chaDaughter.tpcNSigmaPi());
    }
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      if (pDaug < mPidThreshold) {
        this->evaluateObservable(kChaDaughTpcProton, chaDaughter.tpcNSigmaPr());
      } else {
        if (chaDaughter.hasTOF()) {
          this->evaluateObservable(kChaDaughTofProton, std::abs(chaDaughter.tofNSigmaPr()));
        } else {
          this->evaluateObservable(kChaDaughTofProton, 999.f);
        }
      }
    }

    this->assembleBitmask<HistName>();
  };

  template <typename T>
  bool checkFilters(const T& kink) const
  {
    float pt = kink.ptMoth();
    // Compute mother eta and phi
    float px = kink.pxMoth();
    float py = kink.pyMoth();
    float pz = kink.pzMoth();
    float p = std::sqrt(px * px + py * py + pz * pz);
    float eta = (p > 0.f) ? 0.5f * std::log((p + pz) / (p - pz)) : 0.f;
    float phi = RecoDecay::constrainAngle(std::atan2(py, px));

    return ((pt > mPtMin && pt < mPtMax) &&
            (eta > mEtaMin && eta < mEtaMax) &&
            (phi > mPhiMin && phi < mPhiMax));
  }

  template <typename T>
  bool checkMass(T const& kinkCand) const
  {
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      float sigmaMass = kinkCand.mSigmaMinus();
      return (sigmaMass > mMassSigmaLowerLimit && sigmaMass < mMassSigmaUpperLimit);
    }

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      float sigmaMass = kinkCand.mSigmaPlus();
      return (sigmaMass > mMassSigmaPlusLowerLimit && sigmaMass < mMassSigmaPlusUpperLimit);
    }
    return false;
  }

 public:
  float mMassSigmaLowerLimit = 1.15f;
  float mMassSigmaUpperLimit = 1.25f;
  float mMassSigmaPlusLowerLimit = 1.15f;
  float mMassSigmaPlusUpperLimit = 1.25f;
  float mPidThreshold = 0.75f;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -1.f;
  float mEtaMax = 1.f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;
};

struct KinkBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FSigmas> producedSigmas;
  o2::framework::Produces<o2::aod::FSigmaMasks> producedSigmaMasks;
  o2::framework::Produces<o2::aod::FSigmaExtras> producedSigmaExtras;
  o2::framework::Produces<o2::aod::FSigmaPlus> producedSigmaPlus;
  o2::framework::Produces<o2::aod::FSigmaPlusMasks> producedSigmaPlusMasks;
  o2::framework::Produces<o2::aod::FSigmaPlusExtras> producedSigmaPlusExtras;
};

struct ConfKinkTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("KinkTables");
  o2::framework::Configurable<int> produceSigmas{"produceSigmas", -1, "Produce Sigmas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaMasks{"produceSigmaMasks", -1, "Produce SigmaMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaExtras{"produceSigmaExtras", -1, "Produce SigmaExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlus{"produceSigmaPlus", -1, "Produce SigmaPlus (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusMasks{"produceSigmaPlusMasks", -1, "Produce SigmaPlusMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusExtras{"produceSigmaPlusExtras", -1, "Produce SigmaPlusExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::Kink kinkType, char const* HistName>
class KinkBuilder
{
 public:
  KinkBuilder() = default;
  ~KinkBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext)
  {
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      LOG(info) << "Initialize femto Sigma builder...";
      mProduceSigmas = utils::enableTable("FSigmas_001", table.produceSigmas.value, initContext);
      mProduceSigmaMasks = utils::enableTable("FSigmaMasks_001", table.produceSigmaMasks.value, initContext);
      mProduceSigmaExtras = utils::enableTable("FSigmaExtras_001", table.produceSigmaExtras.value, initContext);
    }

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      LOG(info) << "Initialize femto SigmaPlus builder...";
      mProduceSigmaPlus = utils::enableTable("FSigmaPlus_001", table.produceSigmaPlus.value, initContext);
      mProduceSigmaPlusMasks = utils::enableTable("FSigmaPlusMasks_001", table.produceSigmaPlusMasks.value, initContext);
      mProduceSigmaPlusExtras = utils::enableTable("FSigmaPlusExtras_001", table.produceSigmaPlusExtras.value, initContext);
    }

    if (mProduceSigmas || mProduceSigmaMasks || mProduceSigmaExtras || mProduceSigmaPlus || mProduceSigmaPlusMasks || mProduceSigmaPlusExtras) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mKinkSelection.configure(registry, config, filter);
    mKinkSelection.printSelections(kinkSelsName);
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillKinks(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts, T5& kinkProducts, T6 const& kinks, T7 const& tracks, T8& trackBuilder)
  {
    if (!mFillAnyTable) {
      return;
    }
    int64_t daughterIndex = 0;

    for (const auto& kink : kinks) {
      if (!mKinkSelection.checkFilters(kink)) {
        continue;
      }

      if (!mKinkSelection.checkMass(kink)) {
        continue;
      }

      // Apply selections
      mKinkSelection.applySelections(kink, tracks);

      if (!mKinkSelection.passesAllRequiredSelections()) {
        continue;
      }

      auto daughter = kink.template trackDaug_as<T7>();
      collisionBuilder.template fillCollision<system>(collisionProducts, col);
      daughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kKinkDaughter>(daughter, trackProducts, collisionProducts);
      if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
        fillSigma(collisionProducts, kinkProducts, kink, daughterIndex);
      }
      if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
        fillSigmaPlus(collisionProducts, kinkProducts, kink, daughterIndex);
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillSigma(T1& collisionProducts, T2& kinkProducts, T3 const& kink, int64_t daughterIndex)
  {
    float mass = kink.mSigmaMinus();

    if (mProduceSigmas) {
      float pt = kink.ptMoth();
      // Compute mother eta and phi
      float px = kink.pxMoth();
      float py = kink.pyMoth();
      float pz = kink.pzMoth();
      float p = std::sqrt(px * px + py * py + pz * pz);
      float eta = (p > 0.f) ? 0.5f * std::log((p + pz) / (p - pz)) : 0.f;
      float phi = RecoDecay::constrainAngle(std::atan2(py, px));

      kinkProducts.producedSigmas(collisionProducts.producedCollision.lastIndex(),
                                  kink.mothSign() * pt,
                                  eta,
                                  phi,
                                  mass,
                                  daughterIndex);
    }
    if (mProduceSigmaMasks) {
      kinkProducts.producedSigmaMasks(mKinkSelection.getBitmask());
    }
    if (mProduceSigmaExtras) {
      // Compute kink angle using momentum components
      float pxmoth = kink.pxMoth();
      float pymoth = kink.pyMoth();
      float pzmoth = kink.pzMoth();
      float pxch = kink.pxDaug();
      float pych = kink.pyDaug();
      float pzch = kink.pzDaug();

      float pMoth = std::sqrt(pxmoth * pxmoth + pymoth * pymoth + pzmoth * pzmoth);
      float pDaug = std::sqrt(pxch * pxch + pych * pych + pzch * pzch);
      float kinkAngle = 0.f;
      if (pMoth > 0.f && pDaug > 0.f) {
        float dotProduct = pxmoth * pxch + pymoth * pych + pzmoth * pzch;
        float cosAngle = dotProduct / (pMoth * pDaug);
        cosAngle = std::max(-1.0f, std::min(1.0f, cosAngle)); // Clamp
        kinkAngle = std::acos(cosAngle);
      }

      float transRadius = std::hypot(kink.xDecVtx(), kink.yDecVtx());

      kinkProducts.producedSigmaExtras(
        kinkAngle,
        kink.dcaDaugPv(),
        kink.dcaMothPv(),
        kink.xDecVtx(),
        kink.yDecVtx(),
        kink.zDecVtx(),
        transRadius);
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillSigmaPlus(T1& collisionProducts, T2& kinkProducts, T3 const& kink, int64_t daughterIndex)
  {
    float mass = kink.mSigmaPlus();

    if (mProduceSigmaPlus) {
      float pt = kink.ptMoth();
      // Compute mother eta and phi
      float px = kink.pxMoth();
      float py = kink.pyMoth();
      float pz = kink.pzMoth();
      float p = std::sqrt(px * px + py * py + pz * pz);
      float eta = (p > 0.f) ? 0.5f * std::log((p + pz) / (p - pz)) : 0.f;
      float phi = RecoDecay::constrainAngle(std::atan2(py, px));

      kinkProducts.producedSigmaPlus(collisionProducts.producedCollision.lastIndex(),
                                     kink.mothSign() * pt,
                                     eta,
                                     phi,
                                     mass,
                                     daughterIndex);
    }
    if (mProduceSigmaPlusMasks) {
      kinkProducts.producedSigmaPlusMasks(mKinkSelection.getBitmask());
    }
    if (mProduceSigmaPlusExtras) {
      // Compute kink angle using momentum components
      float pxmoth = kink.pxMoth();
      float pymoth = kink.pyMoth();
      float pzmoth = kink.pzMoth();
      float pxch = kink.pxDaug();
      float pych = kink.pyDaug();
      float pzch = kink.pzDaug();

      float pMoth = std::sqrt(pxmoth * pxmoth + pymoth * pymoth + pzmoth * pzmoth);
      float pDaug = std::sqrt(pxch * pxch + pych * pych + pzch * pzch);
      float kinkAngle = 0.f;
      if (pMoth > 0.f && pDaug > 0.f) {
        float dotProduct = pxmoth * pxch + pymoth * pych + pzmoth * pzch;
        float cosAngle = dotProduct / (pMoth * pDaug);
        cosAngle = std::max(-1.0f, std::min(1.0f, cosAngle)); // Clamp
        kinkAngle = std::acos(cosAngle);
      }

      float transRadius = std::hypot(kink.xDecVtx(), kink.yDecVtx());

      kinkProducts.producedSigmaPlusExtras(
        kinkAngle,
        kink.dcaDaugPv(),
        kink.dcaMothPv(),
        kink.xDecVtx(),
        kink.yDecVtx(),
        kink.zDecVtx(),
        transRadius);
    }
  }

  bool fillAnyTable() { return mFillAnyTable; }

 private:
  KinkSelection<kinkType, HistName> mKinkSelection;
  bool mFillAnyTable = false;
  bool mProduceSigmas = false;
  bool mProduceSigmaMasks = false;
  bool mProduceSigmaExtras = false;
  bool mProduceSigmaPlus = false;
  bool mProduceSigmaPlusMasks = false;
  bool mProduceSigmaPlusExtras = false;
};

} // namespace kinkbuilder
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_KINKBUILDER_H_
