// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file twoTrackResonanceSelection.h
/// \brief two track resonance selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCESELECTION_H_

#include "RecoDecay.h"

#include "PWGCF/FemtoUnited/Core/baseSelection.h"
#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/Core/modes.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/Configurable.h"
#include <CommonConstants/PhysicsConstants.h>

#include "Math/Vector4D.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace twotrackresonanceselection
{

struct ConfTwoTrackResonanceDaughterFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TwoTrackResonanceDaughterFilter");
  o2::framework::Configurable<float> posDauPtMin{"posDauPtMin", 0.f, "Minimum pT for positive Daughter"};
  o2::framework::Configurable<float> posDauPtMax{"posDauPtMax", 6.f, "Maximum pT for positive Daughter"};
  o2::framework::Configurable<float> negDauPtMin{"negDauPtMin", 0.f, "Minimum pT for negative Daughter"};
  o2::framework::Configurable<float> negDauPtMax{"negDauPtMax", 6.f, "Maximum pT for negative Daughter"};
  o2::framework::Configurable<float> posDauEtaMin{"posDauEtaMin", -0.9f, "Minimum eta for positive Daughter"};
  o2::framework::Configurable<float> posDauEtaMax{"posDauEtaMax", 0.9f, "Maximum eta for positive Daughter"};
  o2::framework::Configurable<float> negDauEtaMin{"negDauEtaMin", -0.9f, "Minimum eta for negative Daughter"};
  o2::framework::Configurable<float> negDauEtaMax{"negDauEtaMax", 0.9f, "Maximum eta for negative Daughter"};
  o2::framework::Configurable<float> posDauPhiMin{"posDauPhiMin", 0.f, "Minimum phi for positive Daughter"};
  o2::framework::Configurable<float> posDauPhiMax{"posDauPhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi for positive Daughter"};
  o2::framework::Configurable<float> negDauPhiMin{"negDauPhiMin", 0.f, "Minimum phi for negative Daughter"};
  o2::framework::Configurable<float> negDauPhiMax{"negDauPhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi for negative Daughter"};
};

template <const char* Prefix>
struct ConfTwoTrackResonanceFilters : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.f, "Minimum invariant mass for Resonance"};
  o2::framework::Configurable<float> massMax{"massMax", 6.f, "Maximum invariant mass for Resonance"};
};
constexpr const char PrefixRhoFilters[] = "RhoFilters1";
constexpr const char PrefixPhiFilters[] = "PhiFilters1";
constexpr const char PrefixKstarFilters[] = "KstarFilters1";
using ConfRhoFilters = ConfTwoTrackResonanceFilters<PrefixRhoFilters>;
using ConfPhiFilters = ConfTwoTrackResonanceFilters<PrefixPhiFilters>;
using ConfKstarFilters = ConfTwoTrackResonanceFilters<PrefixKstarFilters>;

#define TWOTRACKRESONANCE_DEFAULT_BITS                                                                                                                        \
  o2::framework::Configurable<std::vector<float>> dauEtaMax{"dauEtaMax", {0.8f}, "Maximum |eta| "};                                                           \
  o2::framework::Configurable<std::vector<float>> dauTpcClustersMin{"dauTpcClustersMin", {90.f}, "Minimum number of clusters in TPC"};                        \
  o2::framework::Configurable<std::vector<float>> posDauPtMin{"posDauPtMin", {0.f}, "Minimum pT of positive daughter "};                                      \
  o2::framework::Configurable<std::vector<float>> posDauPtMax{"posDauPtMax", {6.f}, "Maximum pT of the positive daughter"};                                   \
  o2::framework::Configurable<std::vector<std::string>> posDauDcaxyMax{"posDauDcaxyMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_xy| as a function of pT"}; \
  o2::framework::Configurable<std::vector<std::string>> posDauDcazMax{"posDauDcazMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_z| as a function of pT"};    \
  o2::framework::Configurable<std::vector<float>> negDauPtMin{"negDauPtMin", {0.f}, "Minimum pT of negative daughter "};                                      \
  o2::framework::Configurable<std::vector<float>> negDauPtMax{"negDauPtMax", {6.f}, "Maximum pT of the negative daughter"};                                   \
  o2::framework::Configurable<std::vector<std::string>> negDauDcaxyMax{"negDauDcaxyMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_xy| as a function of pT"}; \
  o2::framework::Configurable<std::vector<std::string>> negDauDcazMax{"negDauDcazMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_z| as a function of pT"};    \
  o2::framework::Configurable<float> posDauMinMomentumForTof{"posDauMinMomentumForTof", 0.4f, "Minimum momentum to required TOF PID (positive daughers)"};    \
  o2::framework::Configurable<float> negDauMinMomentumForTof{"negDauMinMomentumForTof", 0.4f, "Minimum momentum to required TOF PID (negative daughers)"};

#define TWOTRACKRESONANCE_PIONPID_BITS                                                                                                                      \
  o2::framework::Configurable<std::vector<float>> posDauItsPion{"posDauItsPion", {}, "Maximum |nsimga_Pion| ITS for positive daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> posDauTpcPion{"posDauTpcPion", {3.f}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> posDauTofPion{"posDauTofPion", {}, "Maximum |nsimga_Pion| TOF for positive daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> posDauTpctofPion{"posDauTpctofPion", {3.f}, "Maximum |nsimga_Pion| TPCTOF for positive daughter tracks"}; \
  o2::framework::Configurable<std::vector<float>> negDauItsPion{"negDauItsPion", {}, "Maximum |nsimga_Pion| ITS for negative daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> negDauTpcPion{"negDauTpcPion", {3.f}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> negDauTofPion{"negDauTofPion", {}, "Maximum |nsimga_Pion| TOF for negative daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> negDauTpctofPion{"negDauTpctofPion", {3.f}, "Maximum |nsimga_Pion| TPCTOF for negative daughter tracks"};

#define TWOTRACKRESONANCE_KAONPID_BITS                                                                                                                      \
  o2::framework::Configurable<std::vector<float>> posDauItsKaon{"posDauItsKaon", {}, "Maximum |nsimga_Kaon| ITS for positive daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> posDauTpcKaon{"posDauTpcKaon", {3.f}, "Maximum |nsimga_Kaon| TPC for positive daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> posDauTofKaon{"posDauTofKaon", {}, "Maximum |nsimga_Kaon| TOF for positive daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> posDauTpctofKaon{"posDauTpctofKaon", {3.f}, "Maximum |nsimga_Kaon| TPCTOF for positive daughter tracks"}; \
  o2::framework::Configurable<std::vector<float>> negDauItsKaon{"negDauItsKaon", {}, "Maximum |nsimga_Kaon| ITS for negative daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> negDauTpcKaon{"negDauTpcKaon", {3.f}, "Maximum |nsimga_Kaon| TPC for negative daughter tracks"};          \
  o2::framework::Configurable<std::vector<float>> negDauTofKaon{"negDauTofKaon", {}, "Maximum |nsimga_Kaon| TOF for negative daughter tracks"};             \
  o2::framework::Configurable<std::vector<float>> negDauTpctofKaon{"negDauTpctofKaon", {3.f}, "Maximum |nsimga_Kaon| TPCTOF for negative daughter tracks"};

struct ConfPhiBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhiBits");
  TWOTRACKRESONANCE_DEFAULT_BITS
  TWOTRACKRESONANCE_KAONPID_BITS
};

struct ConfRho0Bits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Rho0Bits");
  TWOTRACKRESONANCE_DEFAULT_BITS
  TWOTRACKRESONANCE_PIONPID_BITS
};

struct ConfKstar0Bits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Kstar0Bits");
  TWOTRACKRESONANCE_DEFAULT_BITS
  TWOTRACKRESONANCE_PIONPID_BITS
  TWOTRACKRESONANCE_KAONPID_BITS
};

#undef TWOTRACKRESONANCE_DEFAULT_BITS
#undef TWOTRACKRESONANCE_KAONPID_BITS
#undef TWOTRACKRESONANCE_PIONPID_BITS

#define TWOTRACKRESONANCE_DEFAULT_SELECTION(defaultMassMin, defaultMassMax)                                                                                            \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                                                                \
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};                                                                                                \
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};                                                                                           \
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};                                                                                            \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                                                             \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                                                                \
  o2::framework::Configurable<float> massMin{"massMin", defaultMassMin, "Minimum invariant mass for Resonance"};                                                       \
  o2::framework::Configurable<float> massMax{"massMax", defaultMassMax, "Maximum invariant mass for Resonance"};                                                       \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> posDauMaskBelowThres{"posDauMaskBelowThres", 8u, "Bitmask for resonance selection"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> posDauMaskAboveThres{"posDauMaskAboveThres", 4u, "Bitmask for resonance selection"}; \
  o2::framework::Configurable<o2::aod::femtodatatypes::TwoTrackResonanceMaskType> negDauMaskBelowThres{"negDauMaskBelowThres", 2u, "Bitmask for resonance selection"}; \
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

  // selection for positive daughter
  kPosDauPtMin,       ///< min pt
  kPosDauPtMax,       ///< max pt
  kPosDauDcaxyAbsMax, ///< max |DCA_xy|
  kPosDauDcazAbsMax,  ///< max |DCA_z|
  kPosDauTpcPion,     /// < max |nsigma_TPC| for pion
  kPosDauTofPion,     /// < max |nsigma_TOF| for pion
  kPosDauTpctofPion,  /// < max |nsigma_TPC+TOF| for pion
  kPosDauTpcKaon,     /// < max |nsigma_TPC| for kaon
  kPosDauTofKaon,     /// < max |nsigma_TOF| for kaon
  kPosDauTpctofKaon,  /// < max |nsigma_TPC+TOF| for kaon

  // selection for negative daughter
  kNegDauPtMin,       ///< min pt
  kNegDauPtMax,       ///< max pt
  kNegDauDcaxyAbsMax, ///< max |DCA_xy|
  kNegDauDcazAbsMax,  ///< max |DCA_z|
  kNegDauTpcPion,     /// < max |nsigma_TPC| for pion
  kNegDauTofPion,     /// < max |nsigma_TOF| for pion
  kNegDauTpctofPion,  /// < max |nsigma_TPC+TOF| for pion
  kNegDauTpcKaon,     /// < max |nsigma_TPC| for kaon
  kNegDauTofKaon,     /// < max |nsigma_TOF| for kaon
  kNegDauTpctofKaon,  /// < max |nsigma_TPC+TOF| for kaon

  kResonanceSelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <o2::analysis::femtounited::modes::TwoTrackResonance reso>
class TwoTrackResonanceSelection : public BaseSelection<float, o2::aod::femtodatatypes::TwoTrackResonanceMaskType, twotrackresonanceselection::kResonanceSelsMax>
{
 public:
  TwoTrackResonanceSelection() {}
  virtual ~TwoTrackResonanceSelection() = default;

  template <typename T1, typename T2, typename T3>
  void configure(T1 const& config, T2 const& filter, T3 const& daughterFilter)
  {
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kPhi)) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mNegDaughterMass = o2::constants::physics::MassKMinus;
      this->addSelection(config.posDauTpcKaon.value, kPosDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofKaon.value, kPosDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofKaon.value, kPosDauTpctofKaon, limits::kUpperLimit, false, false);
      this->addSelection(config.negDauTpcKaon.value, kNegDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofKaon.value, kNegDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofKaon.value, kNegDauTpctofKaon, limits::kUpperLimit, false, false);
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kRho0)) {
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
      mNegDaughterMass = o2::constants::physics::MassPiMinus;
      this->addSelection(config.posDauTpcPion.value, kPosDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofPion.value, kPosDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofPion.value, kPosDauTpctofPion, limits::kUpperLimit, false, false);
      this->addSelection(config.negDauTpcPion.value, kNegDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofPion.value, kNegDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofPion.value, kNegDauTpctofPion, limits::kUpperLimit, false, false);
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kKstar0)) {
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mNegDaughterMass = o2::constants::physics::MassPiMinus;
      this->addSelection(config.posDauTpcKaon.value, kPosDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofKaon.value, kPosDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofKaon.value, kPosDauTpctofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpcPion.value, kNegDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofPion.value, kNegDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofPion.value, kNegDauTpctofPion, limits::kAbsUpperLimit, false, false);
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kKstarBar0)) {
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
      mNegDaughterMass = o2::constants::physics::MassKMinus;
      this->addSelection(config.posDauTpcPion.value, kPosDauTpcPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTofPion.value, kPosDauTofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.posDauTpctofPion.value, kPosDauTpctofPion, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpcKaon.value, kNegDauTpcKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTofKaon.value, kNegDauTofKaon, limits::kAbsUpperLimit, false, false);
      this->addSelection(config.negDauTpctofKaon.value, kNegDauTpctofKaon, limits::kAbsUpperLimit, false, false);
    }

    mMassMin = filter.massMin.value;
    mMassMax = filter.massMax.value;
    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;

    mPosDauMinimalMomentumForTof = config.posDauMinMomentumForTof.value;
    mNegDauMinimalMomentumForTof = config.negDauMinMomentumForTof.value;

    this->addSelection(config.dauEtaMax.value, kDauEtaAbsMax, limits::kAbsUpperLimit, true, true);
    this->addSelection(config.dauTpcClustersMin.value, kDauTpcClusterMin, limits::kLowerLimit, true, true);

    this->addSelection(config.posDauPtMin.value, kPosDauPtMin, limits::kLowerLimit, true, true);
    this->addSelection(config.posDauPtMax.value, kPosDauPtMax, limits::kUpperLimit, true, true);
    this->addSelection(config.negDauPtMin.value, kNegDauPtMin, limits::kLowerLimit, true, true);
    this->addSelection(config.negDauPtMax.value, kNegDauPtMax, limits::kUpperLimit, true, true);

    this->addSelection(config.posDauDcaxyMax.name, daughterFilter.posDauPtMin.value, daughterFilter.posDauPtMax.value, config.posDauDcaxyMax.value, kPosDauDcaxyAbsMax, limits::kAbsUpperFunctionLimit, true, true);
    this->addSelection(config.posDauDcazMax.name, daughterFilter.posDauPtMin.value, daughterFilter.posDauPtMax.value, config.posDauDcazMax.value, kPosDauDcazAbsMax, limits::kAbsUpperFunctionLimit, true, true);

    this->addSelection(config.negDauDcaxyMax.name, daughterFilter.negDauPtMin.value, daughterFilter.negDauPtMax.value, config.negDauDcaxyMax.value, kNegDauDcaxyAbsMax, limits::kAbsUpperFunctionLimit, true, true);
    this->addSelection(config.negDauDcazMax.name, daughterFilter.negDauPtMin.value, daughterFilter.negDauPtMax.value, config.negDauDcazMax.value, kNegDauDcazAbsMax, limits::kAbsUpperFunctionLimit, true, true);
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
    mPosDauMomAboveThres = positiveDaughter.p() > mPosDauMomAboveThres;
    mNegDauMomAboveThres = negativeDaughter.p() > mNegDauMomAboveThres;
    return !(mPosDauMomAboveThres && !positiveDaughter.hasTOF()) &&
           !(mNegDauMomAboveThres && !negativeDaughter.hasTOF());
  }

  float getPt() const { return mPt; }
  float getEta() const { return mEta; }
  float getPhi() const { return mPhi; }
  float getMass() const { return mMass; }
  bool getPosDauMomAboveThres() { return mPosDauMomAboveThres; }
  bool getNegDauMomAboveThres() { return mNegDauMomAboveThres; }

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

    // positive daughter selections
    this->evaluateObservable(kPosDauPtMin, posDaughter.pt());
    this->evaluateObservable(kPosDauPtMax, posDaughter.pt());

    this->updateLimits(kPosDauDcaxyAbsMax, posDaughter.pt());
    this->evaluateObservable(kPosDauDcaxyAbsMax, posDaughter.dcaXY());
    this->updateLimits(kPosDauDcazAbsMax, posDaughter.pt());
    this->evaluateObservable(kPosDauDcazAbsMax, posDaughter.dcaZ());

    this->evaluateObservable(kPosDauTpcPion, posDaughter.tpcNSigmaPi());
    this->evaluateObservable(kPosDauTofPion, posDaughter.tofNSigmaPi());
    this->evaluateObservable(kPosDauTpctofPion, std::hypot(posDaughter.tpcNSigmaPi(), posDaughter.tofNSigmaPi()));

    this->evaluateObservable(kPosDauTpcKaon, posDaughter.tpcNSigmaKa());
    this->evaluateObservable(kPosDauTofKaon, posDaughter.tofNSigmaKa());
    this->evaluateObservable(kPosDauTpctofKaon, std::hypot(posDaughter.tpcNSigmaKa(), posDaughter.tofNSigmaKa()));

    // negative daughter selections
    this->evaluateObservable(kNegDauPtMin, negDaughter.pt());
    this->evaluateObservable(kNegDauPtMax, negDaughter.pt());

    this->updateLimits(kNegDauDcaxyAbsMax, negDaughter.pt());
    this->evaluateObservable(kNegDauDcaxyAbsMax, negDaughter.dcaXY());
    this->updateLimits(kNegDauDcazAbsMax, negDaughter.pt());
    this->evaluateObservable(kNegDauDcazAbsMax, negDaughter.dcaZ());

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
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kRho0)) {
      return (this->passesOptionalCut(kPosDauTpcPion) || this->passesOptionalCut(kPosDauTofPion) || this->passesOptionalCut(kPosDauTpctofPion)) &&
             (this->passesOptionalCut(kNegDauTpcPion) || this->passesOptionalCut(kNegDauTofPion) || this->passesOptionalCut(kNegDauTpctofPion));
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kPhi)) {
      return (this->passesOptionalCut(kPosDauTpcKaon) || this->passesOptionalCut(kPosDauTofKaon) || this->passesOptionalCut(kPosDauTpctofKaon)) &&
             (this->passesOptionalCut(kNegDauTpcKaon) || this->passesOptionalCut(kNegDauTofKaon) || this->passesOptionalCut(kNegDauTpctofKaon));
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kKstar0)) {
      return (this->passesOptionalCut(kPosDauTpcKaon) || this->passesOptionalCut(kPosDauTofKaon) || this->passesOptionalCut(kPosDauTpctofKaon)) &&
             (this->passesOptionalCut(kNegDauTpcPion) || this->passesOptionalCut(kNegDauTofPion) || this->passesOptionalCut(kNegDauTpctofPion));
    }
    if constexpr (o2::analysis::femtounited::modes::isEqual(reso, o2::analysis::femtounited::modes::TwoTrackResonance::kKstarBar0)) {
      return (this->passesOptionalCut(kPosDauTpcPion) || this->passesOptionalCut(kPosDauTofPion) || this->passesOptionalCut(kPosDauTpctofPion)) &&
             (this->passesOptionalCut(kNegDauTpcKaon) || this->passesOptionalCut(kNegDauTofKaon) || this->passesOptionalCut(kNegDauTpctofKaon));
    }
    return false;
  }

 protected:
  // (cached) kinematic variables of the resonance
  float mPt = 0.f;
  float mEta = 0.f;
  float mPhi = 0.f;
  float mMass = 0.f;

  float mAntiPt = 0.f;
  float mAntiEta = 0.f;
  float mAntiPhi = 0.f;
  float mAntiMass = 0.f;

  // kinematic selections of the resonance
  float mMassMin = 0.f;
  float mMassMax = 6.f;
  float mPtMin = 0.f;
  float mPtMax = 6.f;
  float mEtaMin = -0.9f;
  float mEtaMax = 0.9f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;

  // minimum momentum of the daughers to ask for tof information
  float mPosDauMinimalMomentumForTof = 0.5f;
  float mNegDauMinimalMomentumForTof = 0.5f;
  bool mPosDauMomAboveThres = false;
  bool mNegDauMomAboveThres = false;

  // daughter masses
  float mPosDaughterMass = 0.f;
  float mNegDaughterMass = 0.f;
};

} // namespace twotrackresonanceselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_TWOTRACKRESONANCESELECTION_H_
