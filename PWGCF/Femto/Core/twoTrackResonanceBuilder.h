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

/// \file twoTrackResonanceBuilder.h
/// \brief two track resonance builder
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TWOTRACKRESONANCEBUILDER_H_
#define PWGCF_FEMTO_CORE_TWOTRACKRESONANCEBUILDER_H_

#include "RecoDecay.h"

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "Math/Vector4D.h"

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
namespace twotrackresonancebuilder
{

struct ConfTwoTrackResonanceDaughterFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TwoTrackResonanceDaughterFilter");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.2f, "Minimum pT of daughters"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT of daughters"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta of daughters"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta of daughters"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi of daughters"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi of daughters"};
};

template <const char* Prefix>
struct ConfTwoTrackResonanceFilters : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.2f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.f, "Minimum invariant mass for Resonance"};
  o2::framework::Configurable<float> massMax{"massMax", 6.f, "Maximum invariant mass for Resonance"};
};
constexpr const char PrefixRhoFilters[] = "Rho0Filters1";
constexpr const char PrefixPhiFilters[] = "PhiFilters1";
constexpr const char PrefixKstarFilters[] = "Kstar0Filters1";
using ConfRhoFilters = ConfTwoTrackResonanceFilters<PrefixRhoFilters>;
using ConfPhiFilters = ConfTwoTrackResonanceFilters<PrefixPhiFilters>;
using ConfKstarFilters = ConfTwoTrackResonanceFilters<PrefixKstarFilters>;

#define TWOTRACKRESONANCE_DEFAULT_BITS(posThres, negThres)                                                                                                          \
  o2::framework::Configurable<std::vector<float>> dauEtaMax{"dauEtaMax", {0.8f}, "Maximum |eta| "};                                                                 \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {90.f}, "Minimum number of clusters in TPC"};                              \
  o2::framework::Configurable<std::vector<std::string>> dauDcaxyMax{"dauDcaxyMax", {"0.004 + 0.013*pow(x, -1)"}, "Maximum |dca_xy| as a function of pT"};           \
  o2::framework::Configurable<std::vector<std::string>> dauDcazMax{"dauDcazMax", {"0.004 + 0.013*pow(x, -1)"}, "Maximum |dca_z| as a function of pT"};              \
  o2::framework::Configurable<std::vector<float>> posDauPtMin{"posDauPtMin", {0.2f}, "Minimum pT of positive daughter "};                                           \
  o2::framework::Configurable<std::vector<float>> posDauPtMax{"posDauPtMax", {6.f}, "Maximum pT of the positive daughter"};                                         \
  o2::framework::Configurable<std::vector<float>> negDauPtMin{"negDauPtMin", {0.2f}, "Minimum pT of negative daughter "};                                           \
  o2::framework::Configurable<std::vector<float>> negDauPtMax{"negDauPtMax", {6.f}, "Maximum pT of the negative daughter"};                                         \
  o2::framework::Configurable<std::vector<float>> posDauMinMomForTof{"posDauMinMomForTof", {posThres}, "Minimum momentum to require TOF PID (positive daughters)"}; \
  o2::framework::Configurable<std::vector<float>> negDauMinMomForTof{"negDauMinMomForTof", {negThres}, "Minimum momentum to require TOF PID (negative daughters)"};

#define TWOTRACKRESONANCE_PIONPID_BITS                                                                                                                      \
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDauTpcPion", {3.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> posDauTofPion{"posDauTofPion", {}, "Maximum |nsimga_Pion| TOF for positive daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> posDauTpctofPion{"posDauTpctofPion", {3.f}, "Maximum |nsimga_Pion| TPCTOF for positive daughter tracks"}; \
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDauTpcPion", {3.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> negDauTofPion{"negDauTofPion", {}, "Maximum |nsimga_Pion| TOF for negative daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> negDauTpctofPion{"negDauTpctofPion", {3.f}, "Maximum |nsimga_Pion| TPCTOF for negative daughter tracks"};

#define TWOTRACKRESONANCE_KAONPID_BITS                                                                                                                      \
  o2::framework::Configurable<std::vector<float>> posDauTpcKaon{"posDauTpcKaon", {3.f}, "Maximum |nsimga_Kaon| TPC for positive daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> posDauTofKaon{"posDauTofKaon", {}, "Maximum |nsimga_Kaon| TOF for positive daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> posDauTpctofKaon{"posDauTpctofKaon", {3.f}, "Maximum |nsimga_Kaon| TPCTOF for positive daughter tracks"}; \
  o2::framework::Configurable<std::vector<float>> negDauTpcKaon{"negDauTpcKaon", {3.f}, "Maximum |nsimga_Kaon| TPC for negative daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> negDauTofKaon{"negDauTofKaon", {}, "Maximum |nsimga_Kaon| TOF for negative daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> negDauTpctofKaon{"negDauTpctofKaon", {3.f}, "Maximum |nsimga_Kaon| TPCTOF for negative daughter tracks"};

struct ConfPhiBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhiBits");
  TWOTRACKRESONANCE_DEFAULT_BITS(0.4f, 0.4f)
  TWOTRACKRESONANCE_KAONPID_BITS
};

struct ConfRho0Bits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Rho0Bits");
  TWOTRACKRESONANCE_DEFAULT_BITS(0.5f, 0.5f)
  TWOTRACKRESONANCE_PIONPID_BITS
};

struct ConfKstar0Bits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Kstar0Bits");
  TWOTRACKRESONANCE_DEFAULT_BITS(0.5f, 0.4f)
  TWOTRACKRESONANCE_PIONPID_BITS
  TWOTRACKRESONANCE_KAONPID_BITS
};

#undef TWOTRACKRESONANCE_DEFAULT_BITS
#undef TWOTRACKRESONANCE_KAONPID_BITS
#undef TWOTRACKRESONANCE_PIONPID_BITS

#define TWOTRACKRESONANCE_DEFAULT_SELECTION(defaultMassMin, defaultMassMax)                                                                                                            \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                                                                                \
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};                                                                                                                \
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};                                                                                                           \
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};                                                                                                            \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                                                                             \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                                                                                \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Resonance"};                                                                       \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Resonance"};                                                                       \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> posDauBitForThres{"posDauBitForThres", 32u, "Bit marking momentum threshold for positive daughter"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> posDauMaskBelowThres{"posDauMaskBelowThres", 16u, "Bitmask for resonance selection"};                \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> posDauMaskAboveThres{"posDauMaskAboveThres", 8u, "Bitmask for resonance selection"};                 \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> negDauBitForThres{"negDauBitForThres", 4u, "Bit marking resonance selection for negative daughter"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> negDauMaskBelowThres{"negDauMaskBelowThres", 2u, "Bitmask for resonance selection"};                 \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> negDauMaskAboveThres{"negDauMaskAboveThres", 1u, "Bitmask for resonance selection"};

struct ConfPhiSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhiSelection");
  TWOTRACKRESONANCE_DEFAULT_SELECTION(0.95f, 1.05f)
};

struct ConfRho0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("RhoSelection");
  TWOTRACKRESONANCE_DEFAULT_SELECTION(0.7f, 0.84f)
};

struct ConfKstar0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Kstar0Selection");
  o2::framework::Configurable<int> sign{"sign", 1, "Sign (+1 for Kstar0 and -1 for Kstar0Bar) "};
  TWOTRACKRESONANCE_DEFAULT_SELECTION(0.8f, 1.0f)
};

#undef TWOTRACKRESONANCE_DEFAULT_SELECTION

/// The different selections this task is capable of doing
enum TwoTrackResonanceSels {

  // common selections for both daughters
  kDauEtaAbsMax,     ///< max |eta|
  kDauTpcClusterMin, ///< min number of TPC cluster
  kDauDcaxyAbsMax,   ///< max |DCA_xy|
  kDauDcazAbsMax,    ///< max |DCA_z|

  // selection for positive daughter
  // add one bit for the momentum threshold
  // when the partition for a resonance is build, we do not have information about the daughter tracks so have to store everything needed for the selection here
  kPosDauMinMomForTof, ///< min p for TOF
  kPosDauPtMin,        ///< min pt
  kPosDauPtMax,        ///< max pt
  kPosDauTpcPion,      /// < max |nsigma_TPC| for pion
  kPosDauTofPion,      /// < max |nsigma_TOF| for pion
  kPosDauTpctofPion,   /// < max |nsigma_TPC+TOF| for pion
  kPosDauTpcKaon,      /// < max |nsigma_TPC| for kaon
  kPosDauTofKaon,      /// < max |nsigma_TOF| for kaon
  kPosDauTpctofKaon,   /// < max |nsigma_TPC+TOF| for kaon

  // selection for negative daughter
  kNegDauMinMomForTof, ///< min p for TOF
  kNegDauPtMin,        ///< min pt
  kNegDauPtMax,        ///< max pt
  kNegDauTpcPion,      /// < max |nsigma_TPC| for pion
  kNegDauTofPion,      /// < max |nsigma_TOF| for pion
  kNegDauTpctofPion,   /// < max |nsigma_TPC+TOF| for pion
  kNegDauTpcKaon,      /// < max |nsigma_TPC| for kaon
  kNegDauTofKaon,      /// < max |nsigma_TOF| for kaon
  kNegDauTpctofKaon,   /// < max |nsigma_TPC+TOF| for kaon

  kResonanceSelsMax
};

const char twoTrackResonanceSelsName[] = "TwoTrackResonance Selection Object";
const std::unordered_map<TwoTrackResonanceSels, std::string> twoTrackResonanceSelsToString = {
  {kDauEtaAbsMax, "Max. |eta| of daughters"},
  {kDauTpcClusterMin, "Min. number of TPC clusters of daughters"},
  {kDauDcaxyAbsMax, "Max. |DCA_xy| of daughters"},
  {kDauDcazAbsMax, "Max. |DCA_z| of the daughters"},
  {kPosDauMinMomForTof, "Min. p of TOF PID of positive daughter"},
  {kPosDauPtMin, "Min. pt of positive daughter"},
  {kPosDauPtMax, "Max. pt of positive daughter"},
  {kPosDauTpcPion, "Max. |sigma_TPC| for pion of positive daughter"},
  {kPosDauTofPion, "Max. |sigma_TOF| for pion of positive daughter"},
  {kPosDauTpctofPion, "Max. |sigma_TPCTOF| for pion of positive daughter"},
  {kPosDauTpcKaon, "Max. |sigma_TPC| for kaon of positive daughter"},
  {kPosDauTofKaon, "Max. |sigma_TOF| for kaon of positive daughter"},
  {kPosDauTpctofKaon, "Max. |sigma_TPCTOF| for kaon of positive daughter"},
  {kNegDauMinMomForTof, "Min. p for TOF PID of negative daughter"},
  {kNegDauPtMin, "Min. pt of negative daughter"},
  {kNegDauPtMax, "Max. pt of negative daughter"},
  {kNegDauTpcPion, "Max. |sigma_TPC| for pion of negative daughter"},
  {kNegDauTofPion, "Max. |sigma_TOF| for pion of negative daughter"},
  {kNegDauTpctofPion, "Max. |sigma_TPCTOF| for pion of negative daughter"},
  {kNegDauTpcKaon, "Max. |sigma_TPC| for kaon of negative daughter"},
  {kNegDauTofKaon, "Max. |sigma_TOF| for kaon of negative daughter"},
  {kNegDauTpctofKaon, "Max. |sigma_TPCTOF| for kaon of negative daughter"}};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <modes::TwoTrackResonance resoType>
class TwoTrackResonanceSelection : public BaseSelection<float, o2::aod::femtodatatypes::TwoTrackResonanceMaskType, kResonanceSelsMax>
{
 public:
  TwoTrackResonanceSelection() {}
  virtual ~TwoTrackResonanceSelection() = default;

  template <typename T1, typename T2, typename T3>
  void configure(T1& config, T2& filter, T3& daughterFilter)
  {
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kPhi)) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mNegDaughterMass = o2::constants::physics::MassKMinus;
      this->addSelection(config.posDauTpcKaon.value, kPosDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofKaon.value, kPosDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofKaon.value, kPosDauTpctofKaon, limits::kUpperLimit, false, false);
      this->addSelection(config.negDauTpcKaon.value, kNegDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofKaon.value, kNegDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofKaon.value, kNegDauTpctofKaon, limits::kUpperLimit, false, false);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kRho0)) {
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
      mNegDaughterMass = o2::constants::physics::MassPiMinus;
      this->addSelection(config.posDauTpcPion.value, kPosDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofPion.value, kPosDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofPion.value, kPosDauTpctofPion, limits::kUpperLimit, false, false);
      this->addSelection(config.negDauTpcPion.value, kNegDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofPion.value, kNegDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofPion.value, kNegDauTpctofPion, limits::kUpperLimit, false, false);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0)) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mNegDaughterMass = o2::constants::physics::MassPiMinus;
      this->addSelection(config.posDauTpcKaon.value, kPosDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofKaon.value, kPosDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofKaon.value, kPosDauTpctofKaon, limits::kUpperLimit, false, false);
      this->addSelection(config.negDauTpcPion.value, kNegDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofPion.value, kNegDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofPion.value, kNegDauTpctofPion, limits::kUpperLimit, false, false);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
      mNegDaughterMass = o2::constants::physics::MassKMinus;
      this->addSelection(config.posDauTpcPion.value, kPosDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofPion.value, kPosDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofPion.value, kPosDauTpctofPion, limits::kUpperLimit, false, false);
      this->addSelection(config.negDauTpcKaon.value, kNegDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofKaon.value, kNegDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofKaon.value, kNegDauTpctofKaon, limits::kUpperLimit, false, false);
    }

    mMassMin = filter.massMin.value;
    mMassMax = filter.massMax.value;
    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;

    this->addSelection(config.dauEtaMax.value, kDauEtaAbsMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, kDauTpcClusterMin, limits::kLowerLimit, true, true);
    this->addSelection(config.dauDcaxyMax.name, daughterFilter.ptMin.value, daughterFilter.ptMax.value, config.dauDcaxyMax.value, kDauDcaxyAbsMax, limits::kAbsUpperFunctionLimit, true, true);
    this->addSelection(config.dauDcazMax.name, daughterFilter.ptMin.value, daughterFilter.ptMax.value, config.dauDcazMax.value, kDauDcazAbsMax, limits::kAbsUpperFunctionLimit, true, true);

    this->addSelection(config.posDauMinMomForTof.value, kPosDauMinMomForTof, limits::kLowerLimit, false, false); // momentum threshold for TOF is no minimal cut
    this->addSelection(config.posDauPtMin.value, kPosDauPtMin, limits::kLowerLimit, true, true);
    this->addSelection(config.posDauPtMax.value, kPosDauPtMax, limits::kUpperLimit, true, true);

    this->addSelection(config.negDauMinMomForTof.value, kNegDauMinMomForTof, limits::kLowerLimit, false, false); // momentum threshold for TOF is no minimal cut
    this->addSelection(config.negDauPtMin.value, kNegDauPtMin, limits::kLowerLimit, true, true);
    this->addSelection(config.negDauPtMax.value, kNegDauPtMax, limits::kUpperLimit, true, true);
  };

  template <typename Tracks>
  void reconstructResonance(Tracks const& posDaughter, Tracks const& negDaughter)
  {

    ROOT::Math::PtEtaPhiMVector vecPosDaughter{posDaughter.pt(), posDaughter.eta(), posDaughter.phi(), mPosDaughterMass};
    ROOT::Math::PtEtaPhiMVector vecNegDaughter{negDaughter.pt(), negDaughter.eta(), negDaughter.phi(), mNegDaughterMass};
    ROOT::Math::PtEtaPhiMVector vecResonance = vecPosDaughter + vecNegDaughter;

    // cache kinematics
    mMass = vecResonance.M();
    mPt = vecResonance.Pt();
    mEta = vecResonance.Eta();
    mPhi = RecoDecay::constrainAngle(vecResonance.Phi());
  }

  bool checkFilters() const
  {
    return ((mMass > mMassMin && mMass < mMassMax) &&
            (mPt > mPtMin && mPt < mPtMax) &&
            (mEta > mEtaMin && mEta < mEtaMax) &&
            (mPhi > mPhiMin && mPhi < mPhiMax));
  }

  template <typename T>
  bool hasTofAboveThreshold(T const& positiveDaughter, T const& negativeDaughter)
  {
    bool posDauHasTofAboveThreshold = !this->passesOptionalSelection(kPosDauMinMomForTof) || positiveDaughter.hasTOF(); // is always true if momentum is below threshold
    bool negDauHasTofAboveThreshold = !this->passesOptionalSelection(kNegDauMinMomForTof) || negativeDaughter.hasTOF(); // is always true if momentum is below threshold
    return posDauHasTofAboveThreshold && negDauHasTofAboveThreshold;
  }

  float getPt() const { return mPt; }
  float getEta() const { return mEta; }
  float getPhi() const { return mPhi; }
  float getMass() const { return mMass; }

  template <typename Tracks>
  void applySelections(Tracks const& posDaughter, Tracks const& negDaughter)
  {
    this->reset();
    // for resoanace topological selectsion are in general not possible, so only selections on the daughters are performed

    // common daugher selections
    std::array<float, 2> etaDaughters = {std::fabs(posDaughter.eta()), std::fabs(negDaughter.eta())};
    this->evaluateObservable(kDauEtaAbsMax, *std::max_element(etaDaughters.begin(), etaDaughters.end()));
    std::array<float, 2> tpcClusterDaughters = {1.f * posDaughter.tpcNClsFound(), 1.f * negDaughter.tpcNClsFound()};
    this->evaluateObservable(kDauTpcClusterMin, *std::min_element(tpcClusterDaughters.begin(), tpcClusterDaughters.end()));

    // check pt dependend dca cut on both daughters
    // we apply the same cut to both daughters so we only want to store the result were both daughters survive the cut
    // since momenta of daughters are different, we compute the bitmask for both, combine them with logical AND and keep the result
    uint64_t bitmaskDcaPos, bitmaskDcaNeg, bitmaskDca;
    this->updateLimits(kDauDcaxyAbsMax, posDaughter.pt());
    this->evaluateObservable(kDauDcaxyAbsMax, posDaughter.dcaXY());
    bitmaskDcaPos = this->getBitmask(kDauDcaxyAbsMax);
    this->updateLimits(kDauDcaxyAbsMax, negDaughter.pt());
    this->evaluateObservable(kDauDcaxyAbsMax, negDaughter.dcaXY());
    bitmaskDcaNeg = this->getBitmask(kDauDcaxyAbsMax);
    bitmaskDca = bitmaskDcaPos & bitmaskDcaNeg;
    this->setBitmask(kDauDcaxyAbsMax, bitmaskDca);

    this->updateLimits(kDauDcazAbsMax, posDaughter.pt());
    this->evaluateObservable(kDauDcazAbsMax, posDaughter.dcaZ());
    bitmaskDcaPos = this->getBitmask(kDauDcazAbsMax);
    this->updateLimits(kDauDcazAbsMax, negDaughter.pt());
    this->evaluateObservable(kDauDcazAbsMax, negDaughter.dcaZ());
    bitmaskDcaNeg = this->getBitmask(kDauDcazAbsMax);
    bitmaskDca = bitmaskDcaPos & bitmaskDcaNeg;
    this->setBitmask(kDauDcazAbsMax, bitmaskDca);

    // positive daughter selections
    this->evaluateObservable(kPosDauMinMomForTof, posDaughter.p());
    this->evaluateObservable(kPosDauPtMin, posDaughter.pt());
    this->evaluateObservable(kPosDauPtMax, posDaughter.pt());

    this->evaluateObservable(kPosDauTpcPion, posDaughter.tpcNSigmaPi());
    this->evaluateObservable(kPosDauTofPion, posDaughter.tofNSigmaPi());
    this->evaluateObservable(kPosDauTpctofPion, std::hypot(posDaughter.tpcNSigmaPi(), posDaughter.tofNSigmaPi()));

    this->evaluateObservable(kPosDauTpcKaon, posDaughter.tpcNSigmaKa());
    this->evaluateObservable(kPosDauTofKaon, posDaughter.tofNSigmaKa());
    this->evaluateObservable(kPosDauTpctofKaon, std::hypot(posDaughter.tpcNSigmaKa(), posDaughter.tofNSigmaKa()));

    // negative daughter selections
    this->evaluateObservable(kNegDauMinMomForTof, negDaughter.p());
    this->evaluateObservable(kNegDauPtMin, negDaughter.pt());
    this->evaluateObservable(kNegDauPtMax, negDaughter.pt());

    this->evaluateObservable(kNegDauTpcPion, negDaughter.tpcNSigmaPi());
    this->evaluateObservable(kNegDauTofPion, negDaughter.tofNSigmaPi());
    this->evaluateObservable(kNegDauTpctofPion, std::hypot(negDaughter.tpcNSigmaPi(), negDaughter.tofNSigmaPi()));

    this->evaluateObservable(kNegDauTpcKaon, negDaughter.tpcNSigmaKa());
    this->evaluateObservable(kNegDauTofKaon, negDaughter.tofNSigmaKa());
    this->evaluateObservable(kNegDauTpctofKaon, std::hypot(negDaughter.tpcNSigmaKa(), negDaughter.tofNSigmaKa()));

    this->assembleBitmask();
  };

  bool checkHypothesis()
  {
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kRho0)) {
      return (this->passesOptionalSelection(kPosDauTpcPion) || this->passesOptionalSelection(kPosDauTofPion) || this->passesOptionalSelection(kPosDauTpctofPion)) &&
             (this->passesOptionalSelection(kNegDauTpcPion) || this->passesOptionalSelection(kNegDauTofPion) || this->passesOptionalSelection(kNegDauTpctofPion));
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kPhi)) {
      return (this->passesOptionalSelection(kPosDauTpcKaon) || this->passesOptionalSelection(kPosDauTofKaon) || this->passesOptionalSelection(kPosDauTpctofKaon)) &&
             (this->passesOptionalSelection(kNegDauTpcKaon) || this->passesOptionalSelection(kNegDauTofKaon) || this->passesOptionalSelection(kNegDauTpctofKaon));
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0)) {
      return (this->passesOptionalSelection(kPosDauTpcKaon) || this->passesOptionalSelection(kPosDauTofKaon) || this->passesOptionalSelection(kPosDauTpctofKaon)) &&
             (this->passesOptionalSelection(kNegDauTpcPion) || this->passesOptionalSelection(kNegDauTofPion) || this->passesOptionalSelection(kNegDauTpctofPion));
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
      return (this->passesOptionalSelection(kPosDauTpcPion) || this->passesOptionalSelection(kPosDauTofPion) || this->passesOptionalSelection(kPosDauTpctofPion)) &&
             (this->passesOptionalSelection(kNegDauTpcKaon) || this->passesOptionalSelection(kNegDauTofKaon) || this->passesOptionalSelection(kNegDauTpctofKaon));
    }
    return false;
  }

 protected:
  // (cached) kinematic variables of the resonance
  float mPt = 0.f;
  float mEta = 0.f;
  float mPhi = 0.f;
  float mMass = 0.f;

  // kinematic selections of the resonance
  float mMassMin = 0.f;
  float mMassMax = 6.f;
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -0.9f;
  float mEtaMax = 0.9f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;

  // daughter masses
  float mPosDaughterMass = 0.f;
  float mNegDaughterMass = 0.f;
};

struct TwoTrackResonanceBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FPhis> producedPhis;
  o2::framework::Produces<o2::aod::FPhiMasks> producedPhiMasks;
  o2::framework::Produces<o2::aod::FKstar0s> producedKstars;
  o2::framework::Produces<o2::aod::FKstar0Masks> producedKstarMasks;
  o2::framework::Produces<o2::aod::FRho0s> producedRhos;
  o2::framework::Produces<o2::aod::FRho0Masks> producedRhoMasks;
};

struct ConfTwoTrackResonanceTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TwoTrackResonanceTables");
  o2::framework::Configurable<int> producePhis{"producePhis", -1, "Produce Phis (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producePhiMasks{"producePhiMasks", -1, "Produce PhiMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceKstar0s{"produceKstar0s", -1, "Produce K0stars (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceKstar0Masks{"produceKstar0Masks", -1, "Produce Kstar0Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceRho0s{"produceRho0s", -1, "Produce Rho0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceRho0Masks{"produceRho0Masks", -1, "Produce Rho0Masks (-1: auto; 0 off; 1 on)"};
};

template <modes::TwoTrackResonance resoType>
class TwoTrackResonanceBuilder
{
 public:
  TwoTrackResonanceBuilder() {}
  virtual ~TwoTrackResonanceBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void init(T1& config, T2& filter, T3& daughterFilter, T4& table, T5 initContext)
  {
    twoTrackResonanceSelection.configure(config, filter, daughterFilter);
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kPhi)) {
      LOG(info) << "Initialize femto Phi builder...";
      producePhis = utils::enableTable("FPhis_001", table.producePhis.value, initContext);
      producePhiMasks = utils::enableTable("FPhiMasks_001", table.producePhiMasks.value, initContext);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0) || modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
      if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0)) {
        LOG(info) << "Initialize femto Kstar0 builder...";
      }
      if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
        LOG(info) << "Initialize femto Kstar0Bar builder...";
      }
      produceKstar0s = utils::enableTable("FKstar0s_001", table.produceKstar0s.value, initContext);
      produceKstar0Masks = utils::enableTable("FKstar0Masks_001", table.produceKstar0Masks.value, initContext);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kRho0)) {
      LOG(info) << "Initialize femto Rho0 builder...";
      produceRho0s = utils::enableTable("FRho0s_001", table.produceRho0s.value, initContext);
      produceRho0Masks = utils::enableTable("FRho0Masks_001", table.produceRho0Masks.value, initContext);
    }

    if (producePhis || producePhiMasks || produceKstar0s || produceKstar0Masks || produceRho0s || produceRho0Masks) {
      fillAnyTable = true;
      twoTrackResonanceSelection.printSelections(twoTrackResonanceSelsName, twoTrackResonanceSelsToString);
    } else {
      LOG(info) << "No tables configured";
    }
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void fillResonances(T1& collisionProducts, T2& trackProducts, T3& resonanceProducts, T4& groupPositiveTracks, T5& groupNegativeTracks, T6& trackBuilder, T7& indexMap)
  {
    if (!fillAnyTable) {
      return;
    }
    for (auto const& positiveTrack : groupPositiveTracks) {
      for (auto const& negativeTrack : groupNegativeTracks) {
        this->fillResonance(collisionProducts, trackProducts, resonanceProducts, positiveTrack, negativeTrack, trackBuilder, indexMap);
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fillResonance(T1& collisionProducts, T2& trackProducts, T3& resonanceProducts, T4 const& posDaughter, T4 const& negDaughter, T5& trackBuilder, T6& indexMap)
  {

    twoTrackResonanceSelection.applySelections(posDaughter, negDaughter); // for resonances selection are only applied to daughter tracks
    if (!twoTrackResonanceSelection.hasTofAboveThreshold(posDaughter, negDaughter) || !twoTrackResonanceSelection.passesAllRequiredSelections()) {
      return;
    }
    twoTrackResonanceSelection.reconstructResonance(posDaughter, negDaughter);
    if (!twoTrackResonanceSelection.checkFilters() || !twoTrackResonanceSelection.checkHypothesis()) {
      return;
    }
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    posDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kResonanceDaughter>(posDaughter, trackProducts, collisionProducts, indexMap);
    negDaughterIndex = trackBuilder.template getDaughterIndex<modes::Track::kResonanceDaughter>(negDaughter, trackProducts, collisionProducts, indexMap);

    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kRho0)) {
      if (produceRho0s) {
        resonanceProducts.producedRhos(
          collisionProducts.producedCollision.lastIndex(),
          twoTrackResonanceSelection.getPt(),
          twoTrackResonanceSelection.getEta(),
          twoTrackResonanceSelection.getPhi(),
          twoTrackResonanceSelection.getMass(),
          posDaughterIndex,
          negDaughterIndex);
      }
      if (produceRho0Masks) {
        resonanceProducts.producedRhoMasks(twoTrackResonanceSelection.getBitmask());
      }
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kPhi)) {
      if (producePhis) {
        resonanceProducts.producedPhis(
          collisionProducts.producedCollision.lastIndex(),
          twoTrackResonanceSelection.getPt(),
          twoTrackResonanceSelection.getEta(),
          twoTrackResonanceSelection.getPhi(),
          twoTrackResonanceSelection.getMass(),
          posDaughterIndex,
          negDaughterIndex);
      }
      if (producePhiMasks) {
        resonanceProducts.producedPhiMasks(twoTrackResonanceSelection.getBitmask());
      }
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0)) {
      if (produceKstar0s) {
        resonanceProducts.producedKstars(
          collisionProducts.producedCollision.lastIndex(),
          twoTrackResonanceSelection.getPt(),
          twoTrackResonanceSelection.getEta(),
          twoTrackResonanceSelection.getPhi(),
          twoTrackResonanceSelection.getMass(),
          posDaughterIndex,
          negDaughterIndex);
      }
      if (produceKstar0Masks) {
        resonanceProducts.producedKstarMasks(twoTrackResonanceSelection.getBitmask());
      }
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
      if (produceKstar0s) {
        resonanceProducts.producedKstars(
          collisionProducts.producedCollision.lastIndex(),
          -1.f * twoTrackResonanceSelection.getPt(),
          twoTrackResonanceSelection.getEta(),
          twoTrackResonanceSelection.getPhi(),
          twoTrackResonanceSelection.getMass(),
          posDaughterIndex,
          negDaughterIndex);
      }
      if (produceKstar0Masks) {
        resonanceProducts.producedKstarMasks(twoTrackResonanceSelection.getBitmask());
      }
    }
  }

 private:
  TwoTrackResonanceSelection<resoType> twoTrackResonanceSelection;
  bool fillAnyTable = false;
  bool producePhis = false;
  bool producePhiMasks = false;
  bool produceKstar0s = false;
  bool produceKstar0Masks = false;
  bool produceRho0s = false;
  bool produceRho0Masks = false;
}; // namespace twotrackresonancebuilder

} // namespace twotrackresonancebuilder
} // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_TWOTRACKRESONANCEBUILDER_H_
