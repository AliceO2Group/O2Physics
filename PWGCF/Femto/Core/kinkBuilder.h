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

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto::kinkbuilder
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
  o2::framework::Configurable<float> massMinSigmaPlusPcm{"massMinSigmaPlusPcm", 1.1f, "Minimum mass for SigmaPlus (PCM) hypothesis"};
  o2::framework::Configurable<float> massMaxSigmaPlusPcm{"massMaxSigmaPlusPcm", 1.3f, "Maximum mass for SigmaPlus (PCM) hypothesis"};
  o2::framework::Configurable<float> mlScoreMinSigmaPlusPcm{"mlScoreMinSigmaPlusPcm", -1.f, "Minimum BDT score for SigmaPlus (PCM) (-1: no cut)"};
};

// selections bits for all kinks
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define KINK_DEFAULT_BITS                                                                                                                                 \
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, all Kinks are passed through. Bits for all selections are stored."};      \
  o2::framework::Configurable<std::vector<float>> kinkTopoDcaMax{"kinkTopoDcaMax", {2.0f}, "Maximum kink topological DCA"};                               \
  o2::framework::Configurable<std::vector<float>> transRadMin{"transRadMin", {19.6f}, "Minimum transverse radius (cm)"};                                  \
  o2::framework::Configurable<std::vector<float>> transRadMax{"transRadMax", {100.f}, "Maximum transverse radius (cm)"};                                  \
  o2::framework::Configurable<std::vector<float>> dauAbsEtaMax{"dauAbsEtaMax", {1.0f}, "Maximum absolute pseudorapidity for daughter track"};             \
  o2::framework::Configurable<std::vector<float>> dauDcaPvMin{"dauDcaPvMin", {0.1f}, "Minimum DCA of daughter from primary vertex (cm)"};                 \
  o2::framework::Configurable<std::vector<float>> mothDcaPvMax{"mothDcaPvMax", {1.0f}, "Maximum DCA of mother from primary vertex (cm)"};                 \
  o2::framework::Configurable<std::vector<float>> qtAPMin{"qtAPMin", {0.15f}, "Minimum qT_AP for Sigma candidates"};                                      \
  o2::framework::Configurable<std::vector<float>> qtAPMax{"qtAPMax", {0.2f}, "Maximum qT_AP for Sigma candidates"};                                       \
  o2::framework::Configurable<std::vector<float>> cosPointingAngleMin{"cosPointingAngleMin", {0.0f}, "Minimum cosine of pointing angle"};                 \
  o2::framework::Configurable<std::vector<float>> ptOriginalMin{"ptOriginalMin", {1.2f}, "Minimum original (not recalculated) pT of the mother (GeV/c)"}; \
  o2::framework::Configurable<std::vector<float>> ptOriginalMax{"ptOriginalMax", {10.f}, "Maximum original (not recalculated) pT of the mother (GeV/c)"};

// derived selection bits for sigma
struct ConfSigmaBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("SigmaBits");
  KINK_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> alphaAPMin{"alphaAPMin", {-1.0f}, "Minimum Alpha_AP for Sigma candidates"};
  o2::framework::Configurable<std::vector<float>> alphaAPMax{"alphaAPMax", {0.0f}, "Maximum Alpha_AP for Sigma candidates"};
  o2::framework::Configurable<std::vector<float>> chaDauTpcPion{"chaDauTpcPion", {5.f}, "Maximum |nsigma_Pion| TPC for charged daughter tracks"};
};

// derived selection bits for sigma plus
struct ConfSigmaPlusBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("SigmaPlusBits");
  KINK_DEFAULT_BITS
  o2::framework::Configurable<std::vector<float>> alphaAPMin{"alphaAPMin", {0.0f}, "Minimum Alpha_AP for SigmaPlus candidates"};
  o2::framework::Configurable<std::vector<float>> alphaAPMax{"alphaAPMax", {1.0f}, "Maximum Alpha_AP for SigmaPlus candidates"};
  o2::framework::Configurable<std::vector<float>> chaDauTpcProton{"chaDauTpcProton", {5.f}, "Maximum |nsigma_Proton| TPC for charged daughter tracks"};
  o2::framework::Configurable<std::vector<float>> chaDauTofProton{"chaDauTofProton", {}, "Maximum |nsigma_Proton| TOF for charged daughter tracks"};
  o2::framework::Configurable<bool> requireTof{"requireTof", false, "If true, TOF PID is a minimal selection. If false, TOF PID is optional"};
  o2::framework::Configurable<bool> keepTracksWithoutTof{"keepTracksWithoutTof", true, "If true, the bit mask for the TOF selection will be true for all limits if the daughter track has no TOF"};
};

#undef KINK_DEFAULT_BITS

// selection bits for sigma plus reconstructed with a PCM photon
struct ConfSigmaPlusPcmBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("SigmaPlusPcmBits");
  o2::framework::Configurable<bool> passThrough{"passThrough", false, "If true, all SigmaPlusPcm are passed through. Bits for all selections are stored."};
  o2::framework::Configurable<std::vector<float>> flightDirTiltMax{"flightDirTiltMax", {0.05f}, "Maximum tilt of the flight direction (rad)"};
  o2::framework::Configurable<std::vector<float>> mothDcaPvMax{"mothDcaPvMax", {1.0f}, "Maximum DCA of mother from primary vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> chaDauDcaZMax{"chaDauDcaZMax", {0.3f}, "Maximum |DCAz| of charged daughter from primary vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> photonDcaZMax{"photonDcaZMax", {2.0f}, "Maximum |DCAz| of the photon from primary vertex (cm)"};
  o2::framework::Configurable<std::vector<float>> chaDauTpcProton{"chaDauTpcProton", {3.f}, "Maximum |nsigma_Proton| TPC for charged daughter tracks"};
  o2::framework::Configurable<std::vector<float>> chaDauTofProton{"chaDauTofProton", {}, "Maximum |nsigma_Proton| TOF for charged daughter tracks"};
  o2::framework::Configurable<bool> requireTof{"requireTof", false, "If true, TOF PID is a minimal selection. If false, TOF PID is optional"};
  o2::framework::Configurable<bool> keepTracksWithoutTof{"keepTracksWithoutTof", true, "If true, the bit mask for the TOF selection will be true for all limits if the daughter track has no TOF"};
};

// base selection for analysis task for kinks
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define KINK_DEFAULT_SELECTIONS(defaultMassMin, defaultMassMax, defaultPdgCode)                                          \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", (defaultPdgCode), "PDG code. Select antiparticle via sign"}; \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                  \
  o2::framework::Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};                                                \
  o2::framework::Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};                                             \
  o2::framework::Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};                                              \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                               \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                  \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass"};                     \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass"};                     \
  o2::framework::Configurable<o2::analysis::femto::datatypes::KinkMaskType> mask{"mask", 0x0, "Bitmask for kink selection"};

// base selection for analysis task for sigmas
template <auto& Prefix>
struct ConfSigmaSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_SELECTIONS(1.1, 1.3, 3112)
  o2::framework::Configurable<int> sign{"sign", -1, "Sign of the Sigma mother track (e.g. -1 for Sigma- or +1 for AntiSigma-)"};
};

// base selection for analysis task for sigma plus
template <auto& Prefix>
struct ConfSigmaPlusSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_SELECTIONS(1.1, 1.3, 3222)
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Sigma mother track (e.g. +1 for Sigma+ or -1 for AntiSigma+)"};
};

// base selection for analysis task for sigma plus reconstructed with a PCM photon
template <auto& Prefix>
struct ConfSigmaPlusPcmSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  KINK_DEFAULT_SELECTIONS(1.1, 1.3, 3222)
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the Sigma candidate (e.g. +1 for Sigma+ or -1 for AntiSigma+)"};
  o2::framework::Configurable<float> mlScoreMin{"mlScoreMin", -1.f, "Minimum BDT score (-1: no cut)"};
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
constexpr const char PrefixSigmaPlusPcmSelection1[] = "SigmaPlusPcmSelection1";
constexpr const char PrefixSigmaPlusPcmSelection2[] = "SigmaPlusPcmSelection2";
using ConfSigmaPlusPcmSelection1 = ConfSigmaPlusPcmSelection<PrefixSigmaPlusPcmSelection1>;
using ConfSigmaPlusPcmSelection2 = ConfSigmaPlusPcmSelection<PrefixSigmaPlusPcmSelection2>;

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

  kPtOriginalMin,
  kPtOriginalMax,

  kFlightDirTiltMax,
  kChaDauDcaZMax,
  kPhotonDcaZMax,

  kKinkSelsMax
};

constexpr char SigmaSelHistName[] = "hSigmaSelection";
constexpr char SigmaPlusSelHistName[] = "hSigmaPlusSelection";
constexpr char SigmaPlusPcmSelHistName[] = "hSigmaPlusPcmSelection";
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
  {kCosPointingAngleMin, "cosPointingAngleMin"},
  {kPtOriginalMin, "ptOriginalMin"},
  {kPtOriginalMax, "ptOriginalMax"},
  {kFlightDirTiltMax, "flightDirTiltMax"},
  {kChaDauDcaZMax, "chaDauDcaZMax"},
  {kPhotonDcaZMax, "photonDcaZMax"}};

/// enum for all kink pre-filters (evaluated in checkFilters, before the selection bitmask)
enum KinkFilters {
  kPtMin,
  kPtMax,
  kEtaMin,
  kEtaMax,
  kPhiMin,
  kPhiMax,
  kMassMin,
  kMassMax,
  kMlScoreMin,
  kKinkFiltersMax
};

constexpr char SigmaFilterHistName[] = "hSigmaFilters";
constexpr char SigmaPlusFilterHistName[] = "hSigmaPlusFilters";
constexpr char SigmaPlusPcmFilterHistName[] = "hSigmaPlusPcmFilters";
const std::unordered_map<KinkFilters, std::string> kinkFilterNames = {
  {kPtMin, "ptMin"},
  {kPtMax, "ptMax"},
  {kEtaMin, "etaMin"},
  {kEtaMax, "etaMax"},
  {kPhiMin, "phiMin"},
  {kPhiMax, "phiMax"},
  {kMassMin, "massMin"},
  {kMassMax, "massMax"},
  {kMlScoreMin, "mlScoreMin"}};

/// \class KinkCuts
/// \brief Cut class to contain and execute all cuts applied to kinks
template <modes::Kink kinkType, auto& SelectionHistName, auto& FilterHistName>
class KinkSelection : public baseselection::BaseSelection<float, o2::analysis::femto::datatypes::KinkMaskType, kKinkSelsMax>
{
 public:
  KinkSelection() = default;
  ~KinkSelection() override = default;

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

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      mMassSigmaLowerLimit = filter.massMinSigma.value;
      mMassSigmaUpperLimit = filter.massMaxSigma.value;
      // Only add PID selection if we need it - will be checked at runtime
      this->addSelection(kChaDaughTpcPion, kinkSelectionNames.at(kChaDaughTpcPion), config.chaDauTpcPion.value, limits::kAbsUpperLimit, true, true, false);
    }

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      mMassSigmaPlusLowerLimit = filter.massMinSigmaPlus.value;
      mMassSigmaPlusUpperLimit = filter.massMaxSigmaPlus.value;
      mKeepTracksWithoutTof = config.keepTracksWithoutTof.value;
      this->addSelection(kChaDaughTpcProton, kinkSelectionNames.at(kChaDaughTpcProton), config.chaDauTpcProton.value, limits::kAbsUpperLimit, true, true, false);
      this->addSelection(kChaDaughTofProton, kinkSelectionNames.at(kChaDaughTofProton), config.chaDauTofProton.value, limits::kAbsUpperLimit, true, config.requireTof.value, false);
      this->addComments(kChaDaughTofProton, "requireTof = " + std::to_string(config.requireTof.value) + "; keepTracksWithoutTof = " + std::to_string(mKeepTracksWithoutTof));
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
    this->addSelection(kPtOriginalMin, kinkSelectionNames.at(kPtOriginalMin), config.ptOriginalMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kPtOriginalMax, kinkSelectionNames.at(kPtOriginalMax), config.ptOriginalMax.value, limits::kUpperLimit, true, true, false);

    this->setupSelectionHistogram<SelectionHistName>(registry);

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      this->template setupFilterHistogram<FilterHistName>(
        registry,
        {
          {kinkFilterNames.at(kPtMin), mPtMin},
          {kinkFilterNames.at(kPtMax), mPtMax},
          {kinkFilterNames.at(kEtaMin), mEtaMin},
          {kinkFilterNames.at(kEtaMax), mEtaMax},
          {kinkFilterNames.at(kPhiMin), mPhiMin},
          {kinkFilterNames.at(kPhiMax), mPhiMax},
          {kinkFilterNames.at(kMassMin), mMassSigmaLowerLimit},
          {kinkFilterNames.at(kMassMax), mMassSigmaUpperLimit},
        });
    }
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      this->template setupFilterHistogram<FilterHistName>(
        registry,
        {
          {kinkFilterNames.at(kPtMin), mPtMin},
          {kinkFilterNames.at(kPtMax), mPtMax},
          {kinkFilterNames.at(kEtaMin), mEtaMin},
          {kinkFilterNames.at(kEtaMax), mEtaMax},
          {kinkFilterNames.at(kPhiMin), mPhiMin},
          {kinkFilterNames.at(kPhiMax), mPhiMax},
          {kinkFilterNames.at(kMassMin), mMassSigmaPlusLowerLimit},
          {kinkFilterNames.at(kMassMax), mMassSigmaPlusUpperLimit},
        });
    }
  };

  template <typename T1, typename T2, typename T3>
  void computeKinkKinematics(T1 const& kinkCand, T2 const& /*tracks*/, T3 const& /*col*/)
  {
    std::array<float, 3> momMother = {kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
    float kinkMomP = RecoDecay::p(momMother);
    std::array<float, 3> momDaughter = {kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};
    float kinkDauP = RecoDecay::p(momDaughter);

    // Alpha_AP
    std::array<float, 3> momMissing = {momMother[0] - momDaughter[0], momMother[1] - momDaughter[1], momMother[2] - momDaughter[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momDaughter.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    mAlphaAp = (lQlP + lQlN != 0.f) ? (lQlP - lQlN) / (lQlP + lQlN) : 0.f;

    // qT_AP
    float dp = lQlP;
    float p2V0 = kinkMomP * kinkMomP;
    float p2A = kinkDauP * kinkDauP;
    mQtAp = (p2V0 > 0.f) ? std::sqrt(std::max(0.f, p2A - dp * dp / p2V0)) : 0.f;

    std::array<float, 3> vMother = {kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx()};
    float vMotherNorm = std::sqrt(std::inner_product(vMother.begin(), vMother.end(), vMother.begin(), 0.f));
    mCosPointingAngle = (vMotherNorm > 0.f && kinkMomP > 0.f) ? (std::inner_product(momMother.begin(), momMother.end(), vMother.begin(), 0.f)) / (kinkMomP * vMotherNorm) : 0.f;
    mTransRadius = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx());

    mKinkDauEta = RecoDecay::eta(momDaughter);

    mKinkAngle = 0.f;
    if (kinkMomP > 0.f && kinkDauP > 0.f) {
      float dotProduct = lQlP;
      float cosAngle = dotProduct / (kinkMomP * kinkDauP);
      cosAngle = std::max(-1.0f, std::min(1.0f, cosAngle)); // Clamp
      mKinkAngle = std::acos(cosAngle);
    }
  }

  template <typename T1, typename T2>
  void applySelections(T1 const& kinkCand, T2 const& /*tracks*/)
  {
    this->reset();
    // kink selections
    this->evaluateObservable(kAlphaAPMin, mAlphaAp);
    this->evaluateObservable(kAlphaAPMax, mAlphaAp);
    this->evaluateObservable(kQtAPMin, mQtAp);
    this->evaluateObservable(kQtAPMax, mQtAp);
    this->evaluateObservable(kCosPointingAngleMin, mCosPointingAngle);
    // the stored pT is the recalculated one, the original pT is only available as selection bits
    this->evaluateObservable(kPtOriginalMin, mKinkMotherPtOriginal);
    this->evaluateObservable(kPtOriginalMax, mKinkMotherPtOriginal);
    this->evaluateObservable(kKinkTopoDcaMax, kinkCand.dcaKinkTopo());

    // Compute transRadius
    this->evaluateObservable(kTransRadMin, mTransRadius);
    this->evaluateObservable(kTransRadMax, mTransRadius);

    // Compute daughter eta
    this->evaluateObservable(kDauAbsEtaMax, std::fabs(mKinkDauEta));

    this->evaluateObservable(kDauDcaPvMin, std::abs(kinkCand.dcaDaugPv()));
    this->evaluateObservable(kMothDcaPvMax, std::abs(kinkCand.dcaMothPv()));

    auto chaDaughter = kinkCand.template trackDaug_as<T2>();

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      this->evaluateObservable(kChaDaughTpcPion, chaDaughter.tpcNSigmaPi());
    }
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      this->evaluateObservable(kChaDaughTpcProton, chaDaughter.tpcNSigmaPr());
      if (chaDaughter.hasTOF()) {
        this->evaluateObservable(kChaDaughTofProton, chaDaughter.tofNSigmaPr());
      } else {
        this->evaluateObservable(kChaDaughTofProton, mKeepTracksWithoutTof ? 0.f : 999.f);
      }
    }

    this->assembleBitmask<SelectionHistName>();
  };

  template <typename T>
  void computeKinkMotherKinematics(const T& kinkCand)
  {
    std::array<float, 3> momMother = {kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
    std::array<float, 3> momDaughter = {kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};

    mKinkMotherEta = RecoDecay::eta(momMother);
    mKinkMotherPhi = RecoDecay::phi(momMother);

    // Recalculate pT using kinematic constraints of the decay channel
    float ptRecalc = -999.f;
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      // Sigma- -> pi- n
      ptRecalc = utils::calcPtnew(momMother[0], momMother[1], momMother[2], momDaughter[0], momDaughter[1], momDaughter[2],
                                  o2::constants::physics::MassSigmaMinus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassNeutron);
    }
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      // Sigma+ -> p pi0
      ptRecalc = utils::calcPtnew(momMother[0], momMother[1], momMother[2], momDaughter[0], momDaughter[1], momDaughter[2],
                                  o2::constants::physics::MassSigmaPlus, o2::constants::physics::MassProton, o2::constants::physics::MassPionNeutral);
    }
    mKinkMotherPtOriginal = std::hypot(momMother[0], momMother[1]);
    mKinkMotherPt = (ptRecalc > 0.f) ? ptRecalc : mKinkMotherPtOriginal;
  }

  template <typename T>
  bool checkFilters(const T& kinkCand) const
  {
    bool pass = true;
    bool p = false;

    p = mKinkMotherPt > mPtMin;
    this->template fillFilter<FilterHistName>(kPtMin, p);
    pass &= p;

    p = mKinkMotherPt < mPtMax;
    this->template fillFilter<FilterHistName>(kPtMax, p);
    pass &= p;

    p = mKinkMotherEta > mEtaMin;
    this->template fillFilter<FilterHistName>(kEtaMin, p);
    pass &= p;

    p = mKinkMotherEta < mEtaMax;
    this->template fillFilter<FilterHistName>(kEtaMax, p);
    pass &= p;

    p = mKinkMotherPhi > mPhiMin;
    this->template fillFilter<FilterHistName>(kPhiMin, p);
    pass &= p;

    p = mKinkMotherPhi < mPhiMax;
    this->template fillFilter<FilterHistName>(kPhiMax, p);
    pass &= p;

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      float sigmaMass = kinkCand.mSigmaMinus();

      p = sigmaMass > mMassSigmaLowerLimit;
      this->template fillFilter<FilterHistName>(kMassMin, p);
      pass &= p;

      p = sigmaMass < mMassSigmaUpperLimit;
      this->template fillFilter<FilterHistName>(kMassMax, p);
      pass &= p;
    }

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      float sigmaMass = kinkCand.mSigmaPlus();

      p = sigmaMass > mMassSigmaPlusLowerLimit;
      this->template fillFilter<FilterHistName>(kMassMin, p);
      pass &= p;

      p = sigmaMass < mMassSigmaPlusUpperLimit;
      this->template fillFilter<FilterHistName>(kMassMax, p);
      pass &= p;
    }

    this->template fillFilterSummary<FilterHistName>(pass);

    return this->isPassThrough() || pass;
  }

  [[nodiscard]] float getKinkMotherPt() const { return mKinkMotherPt; }
  [[nodiscard]] float getKinkMotherEta() const { return mKinkMotherEta; }
  [[nodiscard]] float getKinkMotherPhi() const { return mKinkMotherPhi; }
  [[nodiscard]] float getKinkTransRadius() const { return mTransRadius; }
  [[nodiscard]] float getKinkAngle() const { return mKinkAngle; }

 private:
  float mMassSigmaLowerLimit = 1.15f;
  float mMassSigmaUpperLimit = 1.25f;
  float mMassSigmaPlusLowerLimit = 1.15f;
  float mMassSigmaPlusUpperLimit = 1.25f;
  bool mKeepTracksWithoutTof = true;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -1.f;
  float mEtaMax = 1.f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;

  // mother kinematic
  float mKinkMotherPt = 0.f;         // recalculated pT (falls back to the original pT if the recalculation fails)
  float mKinkMotherPtOriginal = 0.f; // pT of the mother track as reconstructed
  float mKinkMotherEta = 0.f;
  float mKinkMotherPhi = 0.f;

  // qa variables
  float mAlphaAp = 0.f;
  float mQtAp = 0.f;
  float mCosPointingAngle = 0.f;
  float mTransRadius = 0.f;
  float mKinkDauEta = 0.f;
  float mKinkAngle = 0.f;
};

struct KinkBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FSigmas> producedSigmas;
  o2::framework::Produces<o2::aod::FLiteSigmas> producedLiteSigmas;
  o2::framework::Produces<o2::aod::FSigmaMasks> producedSigmaMasks;
  o2::framework::Produces<o2::aod::FSigmaExtras> producedSigmaExtras;
  o2::framework::Produces<o2::aod::FSigmaPlus> producedSigmaPlus;
  o2::framework::Produces<o2::aod::FLiteSigmaPlus> producedLiteSigmaPlus;
  o2::framework::Produces<o2::aod::FSigmaPlusMasks> producedSigmaPlusMasks;
  o2::framework::Produces<o2::aod::FSigmaPlusExtras> producedSigmaPlusExtras;
  o2::framework::Produces<o2::aod::FSigmaPlusPcm> producedSigmaPlusPcm;
  o2::framework::Produces<o2::aod::FLiteSigmaPlusPcm> producedLiteSigmaPlusPcm;
  o2::framework::Produces<o2::aod::FSigmaPlusPcmMasks> producedSigmaPlusPcmMasks;
  o2::framework::Produces<o2::aod::FSigmaPlusPcmMlScores> producedSigmaPlusPcmMlScores;
  o2::framework::Produces<o2::aod::FSigmaPlusPcmExtras> producedSigmaPlusPcmExtras;
};

struct ConfKinkTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("KinkTables");
  o2::framework::Configurable<int> produceSigmas{"produceSigmas", -1, "Produce Sigmas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteSigmas{"produceLiteSigmas", -1, "Produce LiteSigmas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaMasks{"produceSigmaMasks", -1, "Produce SigmaMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaExtras{"produceSigmaExtras", -1, "Produce SigmaExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlus{"produceSigmaPlus", -1, "Produce SigmaPlus (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteSigmaPlus{"produceLiteSigmaPlus", -1, "Produce LiteSigmaPlus (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusMasks{"produceSigmaPlusMasks", -1, "Produce SigmaPlusMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusExtras{"produceSigmaPlusExtras", -1, "Produce SigmaPlusExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusPcm{"produceSigmaPlusPcm", -1, "Produce SigmaPlusPcm (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceLiteSigmaPlusPcm{"produceLiteSigmaPlusPcm", -1, "Produce LiteSigmaPlusPcm (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusPcmMasks{"produceSigmaPlusPcmMasks", -1, "Produce SigmaPlusPcmMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusPcmMlScores{"produceSigmaPlusPcmMlScores", -1, "Produce SigmaPlusPcmMlScores (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSigmaPlusPcmExtras{"produceSigmaPlusPcmExtras", -1, "Produce SigmaPlusPcmExtras (-1: auto; 0 off; 1 on)"};
};

template <modes::Kink kinkType, auto& SelectionHistName, auto& FilterHistName>
class KinkBuilder
{
 public:
  KinkBuilder() = default;
  ~KinkBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext, T5& trackBuilder)
  {
    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
      LOG(info) << "Initialize femto Sigma builder...";
      mProduceSigmas = utils::enableTable("FSigmas_002", table.produceSigmas.value, initContext);
      mProduceLiteSigmas = utils::enableTable("FLiteSigmas_001", table.produceLiteSigmas.value, initContext);
      mProduceSigmaMasks = utils::enableTable("FSigmaMasks_001", table.produceSigmaMasks.value, initContext);
      mProduceSigmaExtras = utils::enableTable("FSigmaExtras_001", table.produceSigmaExtras.value, initContext);

      if (mProduceSigmas && mProduceLiteSigmas) {
        LOG(fatal) << "FSigmas and FLiteSigmas are mutually exclusive -- enable only one. "
                   << "FLiteSigmas is meant to replace FSigmas at the producer stage (for better compression in derived data); "
                   << "use the dedicated converter task to reconstruct FSigmas from FLiteSigmas downstream.";
      }
      if (mProduceSigmas && !trackBuilder.producingTracks()) {
        LOG(fatal) << "FSigmas is enabled, but the track builder is not producing FTracks (full precision). "
                   << "FSigmas stores the daughter index into FTracks -- enable TrackTables.produceTracks, "
                   << "or switch to FLiteSigmas if TrackTables.produceLiteTracks is enabled instead.";
      }
      if (mProduceLiteSigmas && !trackBuilder.producingLiteTracks()) {
        LOG(fatal) << "FLiteSigmas is enabled, but the track builder is not producing FLiteTracks. "
                   << "FLiteSigmas stores the daughter index into FLiteTracks -- enable TrackTables.produceLiteTracks, "
                   << "or switch to FSigmas if TrackTables.produceTracks is enabled instead.";
      }
    }

    if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
      LOG(info) << "Initialize femto SigmaPlus builder...";
      mProduceSigmaPlus = utils::enableTable("FSigmaPlus_001", table.produceSigmaPlus.value, initContext);
      mProduceLiteSigmaPlus = utils::enableTable("FLiteSigmaPlus_001", table.produceLiteSigmaPlus.value, initContext);
      mProduceSigmaPlusMasks = utils::enableTable("FSigmaPlusMasks_001", table.produceSigmaPlusMasks.value, initContext);
      mProduceSigmaPlusExtras = utils::enableTable("FSigmaPlusExtras_002", table.produceSigmaPlusExtras.value, initContext);

      if (mProduceSigmaPlus && mProduceLiteSigmaPlus) {
        LOG(fatal) << "FSigmaPlus and FLiteSigmaPlus are mutually exclusive -- enable only one. "
                   << "FLiteSigmaPlus is meant to replace FSigmaPlus at the producer stage (for better compression in derived data); "
                   << "use the dedicated converter task to reconstruct FSigmaPlus from FLiteSigmaPlus downstream.";
      }
      if (mProduceSigmaPlus && !trackBuilder.producingTracks()) {
        LOG(fatal) << "FSigmaPlus is enabled, but the track builder is not producing FTracks (full precision). "
                   << "FSigmaPlus stores the daughter index into FTracks -- enable TrackTables.produceTracks, "
                   << "or switch to FLiteSigmaPlus if TrackTables.produceLiteTracks is enabled instead.";
      }
      if (mProduceLiteSigmaPlus && !trackBuilder.producingLiteTracks()) {
        LOG(fatal) << "FLiteSigmaPlus is enabled, but the track builder is not producing FLiteTracks. "
                   << "FLiteSigmaPlus stores the daughter index into FLiteTracks -- enable TrackTables.produceLiteTracks, "
                   << "or switch to FSigmaPlus if TrackTables.produceTracks is enabled instead.";
      }
    }

    if (mProduceSigmas || mProduceLiteSigmas || mProduceSigmaMasks || mProduceSigmaExtras ||
        mProduceSigmaPlus || mProduceLiteSigmaPlus || mProduceSigmaPlusMasks || mProduceSigmaPlusExtras) {
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
      // compute mother kinematics before checking filters
      mKinkSelection.computeKinkMotherKinematics(kink);
      if (!mKinkSelection.checkFilters(kink)) {
        continue;
      }
      // compute qa variables before applying selections
      mKinkSelection.computeKinkKinematics(kink, tracks, col);
      mKinkSelection.applySelections(kink, tracks);
      if (!mKinkSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillCollision<system>(collisionProducts, col);

      auto daughter = kink.template trackDaug_as<T7>();
      daughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kKinkDaughter>(daughter, trackProducts, collisionBuilder);
      if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
        fillSigma(collisionBuilder, kinkProducts, kink, daughterIndex);
      }
      if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
        auto mother = kink.template trackMoth_as<T7>();
        fillSigmaPlus(collisionBuilder, kinkProducts, kink, daughter, mother, daughterIndex);
      }
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
  void fillMcKinks(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts, T6& kinkProducts, T7 const& kinks, T8 const& tracks, T9& trackBuilder, T10 const& mcParticles, T11& mcBuilder, T12& mcProducts)
  {

    if (!mFillAnyTable) {
      return;
    }
    int64_t daughterIndex = 0;
    for (const auto& kink : kinks) {
      // compute mother kinematics before checking filters
      mKinkSelection.computeKinkMotherKinematics(kink);
      if (!mKinkSelection.checkFilters(kink)) {
        continue;
      }
      // compute qa variables before applying selections
      mKinkSelection.computeKinkKinematics(kink, tracks, col);
      mKinkSelection.applySelections(kink, tracks);
      if (!mKinkSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);

      auto daughter = kink.template trackDaug_as<T8>();
      daughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kKinkDaughter>(daughter, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);

      if constexpr (modes::isEqual(kinkType, modes::Kink::kSigma)) {
        fillSigma(collisionBuilder, kinkProducts, kink, daughterIndex);
        mcBuilder.template fillMcSigmaWithLabel<system>(daughter, mcParticles, mcCols, mcProducts);
      }
      if constexpr (modes::isEqual(kinkType, modes::Kink::kSigmaPlus)) {
        auto mother = kink.template trackMoth_as<T8>();
        fillSigmaPlus(collisionBuilder, kinkProducts, kink, daughter, mother, daughterIndex);
        mcBuilder.template fillMcSigmaPlusWithLabel<system>(daughter, mcParticles, mcCols, mcProducts);
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillSigma(T1& collisionBuilder, T2& kinkProducts, T3 const& kink, int64_t daughterIndex)
  {
    // Mass is calculated from original momentum components stored in kink table
    float mass = kink.mSigmaMinus();
    float signedPt = kink.mothSign() * mKinkSelection.getKinkMotherPt(); // Recalculated pT

    if (mProduceSigmas) {
      kinkProducts.producedSigmas(collisionBuilder.collisionIndex(),
                                  signedPt,
                                  mKinkSelection.getKinkMotherEta(),
                                  mKinkSelection.getKinkMotherPhi(),
                                  mass,
                                  daughterIndex);
    }
    if (mProduceLiteSigmas) {
      kinkProducts.producedLiteSigmas(collisionBuilder.collisionIndex(),
                                      o2::aod::femtobase::lite::binSignedPt(signedPt),
                                      o2::aod::femtobase::lite::binEta(mKinkSelection.getKinkMotherEta()),
                                      o2::aod::femtobase::lite::binPhi(mKinkSelection.getKinkMotherPhi()),
                                      o2::aod::femtokinks::lite::binSigmaMass(mass),
                                      daughterIndex);
    }
    if (mProduceSigmaMasks) {
      kinkProducts.producedSigmaMasks(mKinkSelection.getBitmask());
    }
    if (mProduceSigmaExtras) {
      kinkProducts.producedSigmaExtras(
        mKinkSelection.getKinkAngle(),
        kink.dcaDaugPv(),
        kink.dcaMothPv(),
        kink.xDecVtx(),
        kink.yDecVtx(),
        kink.zDecVtx(),
        mKinkSelection.getKinkTransRadius());
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillSigmaPlus(T1& collisionBuilder, T2& kinkProducts, T3 const& kink, T4 const& daughter, T5 const& mother, int64_t daughterIndex)
  {
    // Mass is calculated from original momentum components stored in kink table
    float mass = kink.mSigmaPlus();
    float signedPt = kink.mothSign() * mKinkSelection.getKinkMotherPt(); // Recalculated pT
    if (mProduceSigmaPlus) {
      kinkProducts.producedSigmaPlus(collisionBuilder.collisionIndex(),
                                     signedPt,
                                     mKinkSelection.getKinkMotherEta(),
                                     mKinkSelection.getKinkMotherPhi(),
                                     mass,
                                     daughterIndex);
    }
    if (mProduceLiteSigmaPlus) {
      kinkProducts.producedLiteSigmaPlus(collisionBuilder.collisionIndex(),
                                         o2::aod::femtobase::lite::binSignedPt(signedPt),
                                         o2::aod::femtobase::lite::binEta(mKinkSelection.getKinkMotherEta()),
                                         o2::aod::femtobase::lite::binPhi(mKinkSelection.getKinkMotherPhi()),
                                         o2::aod::femtokinks::lite::binSigmaMass(mass),
                                         daughterIndex);
    }
    if (mProduceSigmaPlusMasks) {
      kinkProducts.producedSigmaPlusMasks(mKinkSelection.getBitmask());
    }
    if (mProduceSigmaPlusExtras) {
      kinkProducts.producedSigmaPlusExtras(
        mKinkSelection.getKinkAngle(),
        kink.dcaDaugPv(),
        kink.dcaMothPv(),
        kink.xDecVtx(),
        kink.yDecVtx(),
        kink.zDecVtx(),
        mKinkSelection.getKinkTransRadius(),
        kink.dcaKinkTopo(),
        daughter.tpcNSigmaPr(),
        daughter.tpcNSigmaPi(),
        daughter.tofNSigmaPr(),
        kink.ptMoth(),
        mother.itsClusterSizes());
    }
  }

  [[nodiscard]] bool fillAnyTable() const { return mFillAnyTable; }
  [[nodiscard]] bool isPassThrough() const { return mKinkSelection.isPassThrough(); }

 private:
  KinkSelection<kinkType, SelectionHistName, FilterHistName> mKinkSelection;
  bool mFillAnyTable = false;
  bool mProduceSigmas = false;
  bool mProduceLiteSigmas = false;
  bool mProduceSigmaMasks = false;
  bool mProduceSigmaExtras = false;
  bool mProduceSigmaPlus = false;
  bool mProduceLiteSigmaPlus = false;
  bool mProduceSigmaPlusMasks = false;
  bool mProduceSigmaPlusExtras = false;
};

/// \class SigmaPlusPcmSelection
/// \brief Cut class to contain and execute all cuts applied to sigma plus reconstructed with a PCM photon
template <auto& SelectionHistName, auto& FilterHistName>
class SigmaPlusPcmSelection : public baseselection::BaseSelection<float, o2::analysis::femto::datatypes::KinkMaskType, kKinkSelsMax>
{
 public:
  SigmaPlusPcmSelection() = default;
  ~SigmaPlusPcmSelection() override = default;

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
    mMassLowerLimit = filter.massMinSigmaPlusPcm.value;
    mMassUpperLimit = filter.massMaxSigmaPlusPcm.value;
    mMlScoreMin = filter.mlScoreMinSigmaPlusPcm.value;
    mKeepTracksWithoutTof = config.keepTracksWithoutTof.value;

    this->addSelection(kFlightDirTiltMax, kinkSelectionNames.at(kFlightDirTiltMax), config.flightDirTiltMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kMothDcaPvMax, kinkSelectionNames.at(kMothDcaPvMax), config.mothDcaPvMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kChaDauDcaZMax, kinkSelectionNames.at(kChaDauDcaZMax), config.chaDauDcaZMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kPhotonDcaZMax, kinkSelectionNames.at(kPhotonDcaZMax), config.photonDcaZMax.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kChaDaughTpcProton, kinkSelectionNames.at(kChaDaughTpcProton), config.chaDauTpcProton.value, limits::kAbsUpperLimit, true, true, false);
    this->addSelection(kChaDaughTofProton, kinkSelectionNames.at(kChaDaughTofProton), config.chaDauTofProton.value, limits::kAbsUpperLimit, true, config.requireTof.value, false);
    this->addComments(kChaDaughTofProton, "requireTof = " + std::to_string(config.requireTof.value) + "; keepTracksWithoutTof = " + std::to_string(mKeepTracksWithoutTof));

    this->setupSelectionHistogram<SelectionHistName>(registry);
    this->template setupFilterHistogram<FilterHistName>(
      registry,
      {
        {kinkFilterNames.at(kPtMin), mPtMin},
        {kinkFilterNames.at(kPtMax), mPtMax},
        {kinkFilterNames.at(kEtaMin), mEtaMin},
        {kinkFilterNames.at(kEtaMax), mEtaMax},
        {kinkFilterNames.at(kPhiMin), mPhiMin},
        {kinkFilterNames.at(kPhiMax), mPhiMax},
        {kinkFilterNames.at(kMassMin), mMassLowerLimit},
        {kinkFilterNames.at(kMassMax), mMassUpperLimit},
        {kinkFilterNames.at(kMlScoreMin), mMlScoreMin},
      });
  };

  template <typename T>
  void computeKinematics(T const& sigmaPlusPcmCand)
  {
    std::array<float, 3> momMother = {sigmaPlusPcmCand.pxSigmaPlus(), sigmaPlusPcmCand.pySigmaPlus(), sigmaPlusPcmCand.pzSigmaPlus()};
    mPt = RecoDecay::pt(momMother);
    mEta = RecoDecay::eta(momMother);
    mPhi = RecoDecay::phi(momMother);
    mMass = sigmaPlusPcmCand.massSigmaPlus();
    mMlScore = sigmaPlusPcmCand.mlScore();
  }

  template <typename T>
  bool checkFilters(T const& /*sigmaPlusPcmCand*/) const
  {
    bool pass = true;
    bool p = false;

    p = mPt > mPtMin;
    this->template fillFilter<FilterHistName>(kPtMin, p);
    pass &= p;

    p = mPt < mPtMax;
    this->template fillFilter<FilterHistName>(kPtMax, p);
    pass &= p;

    p = mEta > mEtaMin;
    this->template fillFilter<FilterHistName>(kEtaMin, p);
    pass &= p;

    p = mEta < mEtaMax;
    this->template fillFilter<FilterHistName>(kEtaMax, p);
    pass &= p;

    p = mPhi > mPhiMin;
    this->template fillFilter<FilterHistName>(kPhiMin, p);
    pass &= p;

    p = mPhi < mPhiMax;
    this->template fillFilter<FilterHistName>(kPhiMax, p);
    pass &= p;

    p = mMass > mMassLowerLimit;
    this->template fillFilter<FilterHistName>(kMassMin, p);
    pass &= p;

    p = mMass < mMassUpperLimit;
    this->template fillFilter<FilterHistName>(kMassMax, p);
    pass &= p;

    p = mMlScore >= mMlScoreMin;
    this->template fillFilter<FilterHistName>(kMlScoreMin, p);
    pass &= p;

    this->template fillFilterSummary<FilterHistName>(pass);

    return this->isPassThrough() || pass;
  }

  template <typename T1, typename T2>
  void applySelections(T1 const& sigmaPlusPcmCand, T2 const& /*tracks*/)
  {
    this->reset();
    this->evaluateObservable(kFlightDirTiltMax, sigmaPlusPcmCand.flightDirTilt());
    this->evaluateObservable(kMothDcaPvMax, sigmaPlusPcmCand.candDcaToPV());
    this->evaluateObservable(kChaDauDcaZMax, sigmaPlusPcmCand.protonDcaZ());
    this->evaluateObservable(kPhotonDcaZMax, sigmaPlusPcmCand.photonDcaZToPV());

    auto chaDaughter = sigmaPlusPcmCand.template protonTrack_as<T2>();
    this->evaluateObservable(kChaDaughTpcProton, chaDaughter.tpcNSigmaPr());
    if (chaDaughter.hasTOF()) {
      this->evaluateObservable(kChaDaughTofProton, chaDaughter.tofNSigmaPr());
    } else {
      this->evaluateObservable(kChaDaughTofProton, mKeepTracksWithoutTof ? 0.f : 999.f);
    }

    this->assembleBitmask<SelectionHistName>();
  };

  [[nodiscard]] float getPt() const { return mPt; }
  [[nodiscard]] float getEta() const { return mEta; }
  [[nodiscard]] float getPhi() const { return mPhi; }
  [[nodiscard]] float getMass() const { return mMass; }

 private:
  float mMassLowerLimit = 1.1f;
  float mMassUpperLimit = 1.3f;
  float mMlScoreMin = -1.f;
  bool mKeepTracksWithoutTof = true;

  // kinematic filters
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -1.f;
  float mEtaMax = 1.f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;

  // mother kinematic
  float mPt = 0.f;
  float mEta = 0.f;
  float mPhi = 0.f;
  float mMass = 0.f;
  float mMlScore = 0.f;
};

template <auto& SelectionHistName, auto& FilterHistName>
class SigmaPlusPcmBuilder
{
 public:
  SigmaPlusPcmBuilder() = default;
  ~SigmaPlusPcmBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext, T5& trackBuilder)
  {
    LOG(info) << "Initialize femto SigmaPlusPcm builder...";
    mProduceSigmaPlusPcm = utils::enableTable("FSigmaPlusPcm_001", table.produceSigmaPlusPcm.value, initContext);
    mProduceLiteSigmaPlusPcm = utils::enableTable("FLiteSigmaPlusPcm_001", table.produceLiteSigmaPlusPcm.value, initContext);
    mProduceSigmaPlusPcmMasks = utils::enableTable("FSigmaPlusPcmMasks_001", table.produceSigmaPlusPcmMasks.value, initContext);
    mProduceSigmaPlusPcmMlScores = utils::enableTable("FSigmaPlusPcmMlScores_001", table.produceSigmaPlusPcmMlScores.value, initContext);
    mProduceSigmaPlusPcmExtras = utils::enableTable("FSigmaPlusPcmExtras_001", table.produceSigmaPlusPcmExtras.value, initContext);

    if (mProduceSigmaPlusPcm && mProduceLiteSigmaPlusPcm) {
      LOG(fatal) << "FSigmaPlusPcm and FLiteSigmaPlusPcm are mutually exclusive -- enable only one. "
                 << "FLiteSigmaPlusPcm is meant to replace FSigmaPlusPcm at the producer stage (for better compression in derived data); "
                 << "use the dedicated converter task to reconstruct FSigmaPlusPcm from FLiteSigmaPlusPcm downstream.";
    }
    if (mProduceSigmaPlusPcm && !trackBuilder.producingTracks()) {
      LOG(fatal) << "FSigmaPlusPcm is enabled, but the track builder is not producing FTracks (full precision). "
                 << "FSigmaPlusPcm stores the daughter index into FTracks -- enable TrackTables.produceTracks, "
                 << "or switch to FLiteSigmaPlusPcm if TrackTables.produceLiteTracks is enabled instead.";
    }
    if (mProduceLiteSigmaPlusPcm && !trackBuilder.producingLiteTracks()) {
      LOG(fatal) << "FLiteSigmaPlusPcm is enabled, but the track builder is not producing FLiteTracks. "
                 << "FLiteSigmaPlusPcm stores the daughter index into FLiteTracks -- enable TrackTables.produceLiteTracks, "
                 << "or switch to FSigmaPlusPcm if TrackTables.produceTracks is enabled instead.";
    }

    if (mProduceSigmaPlusPcm || mProduceLiteSigmaPlusPcm || mProduceSigmaPlusPcmMasks || mProduceSigmaPlusPcmMlScores || mProduceSigmaPlusPcmExtras) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mSigmaPlusPcmSelection.configure(registry, config, filter);
    mSigmaPlusPcmSelection.printSelections(kinkSelsName);
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void fillSigmaPlusPcms(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4& trackProducts, T5& kinkProducts, T6 const& sigmaPlusPcms, T7 const& tracks, T8& trackBuilder)
  {
    if (!mFillAnyTable) {
      return;
    }
    int64_t daughterIndex = 0;

    for (const auto& sigmaPlusPcm : sigmaPlusPcms) {
      mSigmaPlusPcmSelection.computeKinematics(sigmaPlusPcm);
      if (!mSigmaPlusPcmSelection.checkFilters(sigmaPlusPcm)) {
        continue;
      }
      mSigmaPlusPcmSelection.applySelections(sigmaPlusPcm, tracks);
      if (!mSigmaPlusPcmSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillCollision<system>(collisionProducts, col);

      auto daughter = sigmaPlusPcm.template protonTrack_as<T7>();
      daughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kKinkDaughter>(daughter, trackProducts, collisionBuilder);
      fillSigmaPlusPcm(collisionBuilder, kinkProducts, sigmaPlusPcm, daughterIndex);
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
  void fillMcSigmaPlusPcms(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5& trackProducts, T6& kinkProducts, T7 const& sigmaPlusPcms, T8 const& tracks, T9& trackBuilder, T10 const& mcParticles, T11& mcBuilder, T12& mcProducts)
  {
    if (!mFillAnyTable) {
      return;
    }
    int64_t daughterIndex = 0;

    for (const auto& sigmaPlusPcm : sigmaPlusPcms) {
      mSigmaPlusPcmSelection.computeKinematics(sigmaPlusPcm);
      if (!mSigmaPlusPcmSelection.checkFilters(sigmaPlusPcm)) {
        continue;
      }
      mSigmaPlusPcmSelection.applySelections(sigmaPlusPcm, tracks);
      if (!mSigmaPlusPcmSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);

      auto daughter = sigmaPlusPcm.template protonTrack_as<T8>();
      daughterIndex = trackBuilder.template getDaughterIndex<system, modes::Track::kKinkDaughter>(daughter, trackProducts, mcCols, collisionBuilder, mcParticles, mcBuilder, mcProducts);
      fillSigmaPlusPcm(collisionBuilder, kinkProducts, sigmaPlusPcm, daughterIndex);
      mcBuilder.template fillMcSigmaPlusPcmWithLabel<system>(sigmaPlusPcm, daughter, mcParticles, mcCols, mcProducts);
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillSigmaPlusPcm(T1& collisionBuilder, T2& kinkProducts, T3 const& sigmaPlusPcm, int64_t daughterIndex)
  {
    float signedPt = sigmaPlusPcm.protonSign() * mSigmaPlusPcmSelection.getPt();
    if (mProduceSigmaPlusPcm) {
      kinkProducts.producedSigmaPlusPcm(collisionBuilder.collisionIndex(),
                                        signedPt,
                                        mSigmaPlusPcmSelection.getEta(),
                                        mSigmaPlusPcmSelection.getPhi(),
                                        mSigmaPlusPcmSelection.getMass(),
                                        daughterIndex);
    }
    if (mProduceLiteSigmaPlusPcm) {
      kinkProducts.producedLiteSigmaPlusPcm(collisionBuilder.collisionIndex(),
                                            o2::aod::femtobase::lite::binSignedPt(signedPt),
                                            o2::aod::femtobase::lite::binEta(mSigmaPlusPcmSelection.getEta()),
                                            o2::aod::femtobase::lite::binPhi(mSigmaPlusPcmSelection.getPhi()),
                                            o2::aod::femtokinks::lite::binSigmaMass(mSigmaPlusPcmSelection.getMass()),
                                            daughterIndex);
    }
    if (mProduceSigmaPlusPcmMasks) {
      kinkProducts.producedSigmaPlusPcmMasks(mSigmaPlusPcmSelection.getBitmask());
    }
    if (mProduceSigmaPlusPcmMlScores) {
      kinkProducts.producedSigmaPlusPcmMlScores(sigmaPlusPcm.mlScore());
    }
    if (mProduceSigmaPlusPcmExtras) {
      kinkProducts.producedSigmaPlusPcmExtras(
        sigmaPlusPcm.flightDirTilt(),
        sigmaPlusPcm.candDcaToPV(),
        sigmaPlusPcm.protonDcaXY(),
        sigmaPlusPcm.protonDcaZ(),
        sigmaPlusPcm.nSigmaTPCProton(),
        sigmaPlusPcm.nSigmaTOFProton(),
        sigmaPlusPcm.photonDcaXYToPV(),
        sigmaPlusPcm.photonDcaZToPV(),
        sigmaPlusPcm.photonCosPAToPV(),
        sigmaPlusPcm.photonConvRadius(),
        sigmaPlusPcm.antiSigmaPointingAngle(),
        sigmaPlusPcm.dcaProtonGamma(),
        sigmaPlusPcm.radius());
    }
  }

  [[nodiscard]] bool fillAnyTable() const { return mFillAnyTable; }
  [[nodiscard]] bool isPassThrough() const { return mSigmaPlusPcmSelection.isPassThrough(); }

 private:
  SigmaPlusPcmSelection<SelectionHistName, FilterHistName> mSigmaPlusPcmSelection;
  bool mFillAnyTable = false;
  bool mProduceSigmaPlusPcm = false;
  bool mProduceLiteSigmaPlusPcm = false;
  bool mProduceSigmaPlusPcmMasks = false;
  bool mProduceSigmaPlusPcmMlScores = false;
  bool mProduceSigmaPlusPcmExtras = false;
};

struct ConfKinkTablesDerivedToDerived : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("KinkTables");
  o2::framework::Configurable<int> limitSigma{"limitSigma", 1, "Require at least this many Sigmas in the collision. Set to 0 to skip Sigmas entirely (not written to output)."};
  o2::framework::Configurable<int> limitSigmaPlus{"limitSigmaPlus", 0, "Require at least this many SigmaPlus in the collision. Set to 0 to skip SigmaPlus entirely (not written to output)."};
  o2::framework::Configurable<int> limitSigmaPlusPcm{"limitSigmaPlusPcm", 0, "Require at least this many SigmaPlusPcm in the collision. Set to 0 to skip SigmaPlusPcm entirely (not written to output)."};
};

struct KinkBuilderDerivedToDerivedProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::StoredFSigmas> producedSigmas;
  o2::framework::Produces<o2::aod::StoredFSigmaMasks> producedSigmaMasks;
  o2::framework::Produces<o2::aod::StoredFSigmaPlus> producedSigmaPluses;
  o2::framework::Produces<o2::aod::StoredFSigmaPlusMasks> producedSigmaPlusMasks;
  o2::framework::Produces<o2::aod::StoredFSigmaPlusExtras> producedSigmaPlusExtras;
  o2::framework::Produces<o2::aod::StoredFSigmaPlusPcm> producedSigmaPlusPcms;
  o2::framework::Produces<o2::aod::StoredFSigmaPlusPcmMasks> producedSigmaPlusPcmMasks;
  o2::framework::Produces<o2::aod::StoredFSigmaPlusPcmMlScores> producedSigmaPlusPcmMlScores;
  o2::framework::Produces<o2::aod::StoredFSigmaPlusPcmExtras> producedSigmaPlusPcmExtras;
};

class KinkBuilderDerivedToDerived
{
 public:
  template <typename T>
  void init(T& config)
  {
    mLimitSigma = config.limitSigma.value;
    mLimitSigmaPlus = config.limitSigmaPlus.value;
    mLimitSigmaPlusPcm = config.limitSigmaPlusPcm.value;

    if (mLimitSigma < 0 || mLimitSigmaPlus < 0 || mLimitSigmaPlusPcm < 0) {
      LOG(fatal) << "Kink limits must be non-negative (got " << mLimitSigma << ", " << mLimitSigmaPlus << " and " << mLimitSigmaPlusPcm << "). Breaking...";
    }
    if (mLimitSigma == 0 && mLimitSigmaPlus == 0 && mLimitSigmaPlusPcm == 0) {
      LOG(fatal) << "Sigma, sigmaplus and sigmapluspcm limits are all 0. Breaking...";
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewSigma(T1 const& col, T2 const& /*sigmaTable*/, T3& partitionSigma, T4& cache)
  {
    if (mLimitSigma == 0) { // sigmas disabled, cannot reject on them
      return false;
    }
    auto sigmaSlice = partitionSigma->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return sigmaSlice.size() < static_cast<int64_t>(mLimitSigma);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewSigmaPlus(T1 const& col, T2 const& /*sigmaPlusTable*/, T3& partitionSigmaPlus, T4& cache)
  {
    if (mLimitSigmaPlus == 0) { // sigma pluses disabled, cannot reject on them
      return false;
    }
    auto sigmaPlusSlice = partitionSigmaPlus->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return sigmaPlusSlice.size() < static_cast<int64_t>(mLimitSigmaPlus);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewSigmaPlusPcm(T1 const& col, T2 const& /*sigmaPlusPcmTable*/, T3& partitionSigmaPlusPcm, T4& cache)
  {
    if (mLimitSigmaPlusPcm == 0) { // sigma pluses (PCM) disabled, cannot reject on them
      return false;
    }
    auto sigmaPlusPcmSlice = partitionSigmaPlusPcm->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return sigmaPlusPcmSlice.size() < static_cast<int64_t>(mLimitSigmaPlusPcm);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processSigma(T1 const& col, T2 const& /*sigmaTable*/, T3 const& oldTrackTable, T4& partitionSigma, T5& trackBuilder, T6& cache, T7& newSigmaTable, T8& newTrackTable, T9& newCollisionTable)
  {
    if (mLimitSigma == 0) { // sigmas disabled
      return;
    }

    auto sigmaSlice = partitionSigma->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sigma : sigmaSlice) {
      auto chaDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(sigma.chaDauId(), oldTrackTable));

      int64_t chaDaughterIndex = trackBuilder.getDaughterIndex(chaDaughter, newTrackTable, newCollisionTable);

      newSigmaTable.producedSigmas(newCollisionTable.producedCollision.lastIndex(),
                                   sigma.signedPt(),
                                   sigma.eta(),
                                   sigma.phi(),
                                   sigma.mass(),
                                   chaDaughterIndex);
      newSigmaTable.producedSigmaMasks(sigma.mask());
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processSigmaPlus(T1 const& col, T2 const& /*sigmaPlusTable*/, T3 const& oldTrackTable, T4& partitionSigmaPlus, T5& trackBuilder, T6& cache, T7& newSigmaPlusTable, T8& newTrackTable, T9& newCollisionTable)
  {
    if (mLimitSigmaPlus == 0) { // sigma pluses disabled
      return;
    }

    auto sigmaPlusSlice = partitionSigmaPlus->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sigmaPlus : sigmaPlusSlice) {
      auto chaDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(sigmaPlus.chaDauId(), oldTrackTable));

      int64_t chaDaughterIndex = trackBuilder.getDaughterIndex(chaDaughter, newTrackTable, newCollisionTable);

      newSigmaPlusTable.producedSigmaPluses(newCollisionTable.producedCollision.lastIndex(),
                                            sigmaPlus.signedPt(),
                                            sigmaPlus.eta(),
                                            sigmaPlus.phi(),
                                            sigmaPlus.mass(),
                                            chaDaughterIndex);
      newSigmaPlusTable.producedSigmaPlusMasks(sigmaPlus.mask());
      fillSigmaPlusExtensions(sigmaPlus, newSigmaPlusTable);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15>
  void processSigmaWithMc(T1 const& col, T2 const& /*sigmaTable*/, T3 const& oldTrackTable, T4& partitionSigma, T5& trackBuilder, T6& cache, T7& newSigmaTable, T8& newTrackTable, T9& newCollisionTable,
                          T10& mcBuilder, T11 const& mcCols, T12 const& mcParticles, T13 const& mcMothers, T14 const& mcPartonicMothers, T15& mcProducts)
  {
    if (mLimitSigma == 0) {
      return;
    }

    auto sigmaSlice = partitionSigma->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sigma : sigmaSlice) {
      auto chaDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(sigma.chaDauId(), oldTrackTable));

      int64_t chaDaughterIndex = trackBuilder.getDaughterIndex(chaDaughter, newTrackTable, newCollisionTable, mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcProducts);

      newSigmaTable.producedSigmas(newCollisionTable.producedCollision.lastIndex(),
                                   sigma.signedPt(),
                                   sigma.eta(),
                                   sigma.phi(),
                                   sigma.mass(),
                                   chaDaughterIndex);
      newSigmaTable.producedSigmaMasks(sigma.mask());
      mcBuilder.fillSigmaWithLabel(sigma, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcProducts);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15>
  void processSigmaPlusWithMc(T1 const& col, T2 const& /*sigmaPlusTable*/, T3 const& oldTrackTable, T4& partitionSigmaPlus, T5& trackBuilder, T6& cache, T7& newSigmaPlusTable, T8& newTrackTable, T9& newCollisionTable,
                              T10& mcBuilder, T11 const& mcCols, T12 const& mcParticles, T13 const& mcMothers, T14 const& mcPartonicMothers, T15& mcProducts)
  {
    if (mLimitSigmaPlus == 0) {
      return;
    }

    auto sigmaPlusSlice = partitionSigmaPlus->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sigmaPlus : sigmaPlusSlice) {
      auto chaDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(sigmaPlus.chaDauId(), oldTrackTable));

      int64_t chaDaughterIndex = trackBuilder.getDaughterIndex(chaDaughter, newTrackTable, newCollisionTable, mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcProducts);

      newSigmaPlusTable.producedSigmaPluses(newCollisionTable.producedCollision.lastIndex(),
                                            sigmaPlus.signedPt(),
                                            sigmaPlus.eta(),
                                            sigmaPlus.phi(),
                                            sigmaPlus.mass(),
                                            chaDaughterIndex);
      newSigmaPlusTable.producedSigmaPlusMasks(sigmaPlus.mask());
      fillSigmaPlusExtensions(sigmaPlus, newSigmaPlusTable);
      mcBuilder.fillSigmaPlusWithLabel(sigmaPlus, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcProducts);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processSigmaPlusPcm(T1 const& col, T2 const& /*sigmaPlusPcmTable*/, T3 const& oldTrackTable, T4& partitionSigmaPlusPcm, T5& trackBuilder, T6& cache, T7& newSigmaPlusPcmTable, T8& newTrackTable, T9& newCollisionTable)
  {
    if (mLimitSigmaPlusPcm == 0) { // sigma pluses (PCM) disabled
      return;
    }

    auto sigmaPlusPcmSlice = partitionSigmaPlusPcm->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sigmaPlusPcm : sigmaPlusPcmSlice) {
      auto chaDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(sigmaPlusPcm.chaDauId(), oldTrackTable));

      int64_t chaDaughterIndex = trackBuilder.getDaughterIndex(chaDaughter, newTrackTable, newCollisionTable);

      newSigmaPlusPcmTable.producedSigmaPlusPcms(newCollisionTable.producedCollision.lastIndex(),
                                                 sigmaPlusPcm.signedPt(),
                                                 sigmaPlusPcm.eta(),
                                                 sigmaPlusPcm.phi(),
                                                 sigmaPlusPcm.mass(),
                                                 chaDaughterIndex);
      newSigmaPlusPcmTable.producedSigmaPlusPcmMasks(sigmaPlusPcm.mask());
      fillSigmaPlusPcmExtensions(sigmaPlusPcm, newSigmaPlusPcmTable);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15>
  void processSigmaPlusPcmWithMc(T1 const& col, T2 const& /*sigmaPlusPcmTable*/, T3 const& oldTrackTable, T4& partitionSigmaPlusPcm, T5& trackBuilder, T6& cache, T7& newSigmaPlusPcmTable, T8& newTrackTable, T9& newCollisionTable,
                                 T10& mcBuilder, T11 const& mcCols, T12 const& mcParticles, T13 const& mcMothers, T14 const& mcPartonicMothers, T15& mcProducts)
  {
    if (mLimitSigmaPlusPcm == 0) {
      return;
    }

    auto sigmaPlusPcmSlice = partitionSigmaPlusPcm->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& sigmaPlusPcm : sigmaPlusPcmSlice) {
      auto chaDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(sigmaPlusPcm.chaDauId(), oldTrackTable));

      int64_t chaDaughterIndex = trackBuilder.getDaughterIndex(chaDaughter, newTrackTable, newCollisionTable, mcBuilder, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcProducts);

      newSigmaPlusPcmTable.producedSigmaPlusPcms(newCollisionTable.producedCollision.lastIndex(),
                                                 sigmaPlusPcm.signedPt(),
                                                 sigmaPlusPcm.eta(),
                                                 sigmaPlusPcm.phi(),
                                                 sigmaPlusPcm.mass(),
                                                 chaDaughterIndex);
      newSigmaPlusPcmTable.producedSigmaPlusPcmMasks(sigmaPlusPcm.mask());
      fillSigmaPlusPcmExtensions(sigmaPlusPcm, newSigmaPlusPcmTable);
      mcBuilder.fillSigmaPlusPcmWithLabel(sigmaPlusPcm, mcCols, mcParticles, mcMothers, mcPartonicMothers, mcProducts);
    }
  }

 private:
  template <typename T1, typename T2>
  void fillSigmaPlusExtensions(T1 const& sigmaPlus, T2& newSigmaPlusTable)
  {
    if constexpr (requires(T1 t) { t.transRadius(); }) {
      newSigmaPlusTable.producedSigmaPlusExtras(sigmaPlus.kinkAngle(),
                                                sigmaPlus.dcaDaugToPV(),
                                                sigmaPlus.dcaMothToPV(),
                                                sigmaPlus.decayVtxX(),
                                                sigmaPlus.decayVtxY(),
                                                sigmaPlus.decayVtxZ(),
                                                sigmaPlus.transRadius(),
                                                sigmaPlus.dcaKinkTopo(),
                                                sigmaPlus.dauTpcNSigmaPr(),
                                                sigmaPlus.dauTpcNSigmaPi(),
                                                sigmaPlus.dauTofNSigmaPr(),
                                                sigmaPlus.ptOriginal(),
                                                sigmaPlus.itsClusterSizes());
    }
  }

  template <typename T1, typename T2>
  void fillSigmaPlusPcmExtensions(T1 const& sigmaPlusPcm, T2& newSigmaPlusPcmTable)
  {
    if constexpr (requires(T1 t) { t.mlScore(); }) {
      newSigmaPlusPcmTable.producedSigmaPlusPcmMlScores(sigmaPlusPcm.mlScore());
    }
    if constexpr (requires(T1 t) { t.photonConvRadius(); }) {
      newSigmaPlusPcmTable.producedSigmaPlusPcmExtras(sigmaPlusPcm.flightDirTilt(),
                                                      sigmaPlusPcm.dcaMothToPV(),
                                                      sigmaPlusPcm.dauDcaXY(),
                                                      sigmaPlusPcm.dauDcaZ(),
                                                      sigmaPlusPcm.dauTpcNSigmaPr(),
                                                      sigmaPlusPcm.dauTofNSigmaPr(),
                                                      sigmaPlusPcm.photonDcaXYToPV(),
                                                      sigmaPlusPcm.photonDcaZToPV(),
                                                      sigmaPlusPcm.photonCosPAToPV(),
                                                      sigmaPlusPcm.photonConvRadius(),
                                                      sigmaPlusPcm.antiSigmaPointingAngle(),
                                                      sigmaPlusPcm.dcaDauPhoton(),
                                                      sigmaPlusPcm.transRadius());
    }
  }

  int mLimitSigma = 0;
  int mLimitSigmaPlus = 0;
  int mLimitSigmaPlusPcm = 0;
};
} // namespace o2::analysis::femto::kinkbuilder
#endif // PWGCF_FEMTO_CORE_KINKBUILDER_H_
