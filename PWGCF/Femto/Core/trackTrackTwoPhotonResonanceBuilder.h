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

/// \file trackTrackTwoPhotonResonanceBuilder.h
/// \brief builder for resonances reconstructed from two tracks and a two-photon resonance
///        (e.g. omega(782) -> pi+ pi- pi0)
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TRACKTRACKTWOPHOTONRESONANCEBUILDER_H_
#define PWGCF_FEMTO_CORE_TRACKTRACKTWOPHOTONRESONANCEBUILDER_H_

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <string>

namespace o2::analysis::femto::trackstwophotonresonancebuilder
{
// filters applied in the producer task (loose pre-selection, applied before the mass window itself)
struct ConfOmegaMesonFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("OmegaMesonFilters");
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 10.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.6f, "Minimum invariant mass"};
  o2::framework::Configurable<float> massMax{"massMax", 0.95f, "Maximum invariant mass"};
};

// selection used downstream (QA/pairing) to re-select already produced resonances by mass window
// and daughter mask. both charged daughters are pions (just opposite charge, no species swap the
// way kstar's kaon/pion does), so -- same as TWOTRACKRESONANCE_SAME_SPECIES_MASK for phi/rho0 --
// a single pair of mask thresholds is applied to both. plus a 3rd mask for the neutral daughter (no
// momentum-threshold split there, same as femtotwophotonresonances' single dauMask)
template <auto& Prefix>
struct ConfOmegaMesonSelection : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", o2::constants::physics::Pdg::kOmega, "Absolute PDG code of the resonance (pdg[omega(782)]=223)"};
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 10.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  o2::framework::Configurable<float> massMin{"massMin", 0.6f, "Minimum invariant mass"};
  o2::framework::Configurable<float> massMax{"massMax", 0.95f, "Maximum invariant mass"};
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> pionDauMaskBelowThres{"pionDauMaskBelowThres", 0x10u, "Bitmask for both pion daughters below threshold"};
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> pionDauMaskAboveThres{"pionDauMaskAboveThres", 0x8u, "Bitmask for both pion daughters above threshold"};
  o2::framework::Configurable<o2::analysis::femto::datatypes::PhotonMaskType> neutralDauMask{"neutralDauMask", 0, "Bitmask required for the neutral (two-photon-resonance) daughter"};
};

constexpr const char PrefixOmegaMesonSelection1[] = "OmegaMesonSelection1";
using ConfOmegaMesonSelection1 = ConfOmegaMesonSelection<PrefixOmegaMesonSelection1>;

struct TrackTrackTwoPhotonResonanceBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FOmegaMesons> producedOmegaMesons;
  o2::framework::Produces<o2::aod::FOmegaMesonMasks> producedOmegaMesonMasks;
};

struct ConfTrackTrackTwoPhotonResonanceTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TrackTrackTwoPhotonResonanceTables");
  o2::framework::Configurable<int> produceOmegaMesons{"produceOmegaMesons", -1, "Produce OmegaMesons (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOmegaMesonMasks{"produceOmegaMesonMasks", -1, "Produce OmegaMesonMasks (-1: auto; 0 off; 1 on)"};
};

/// \brief Builder for resonances reconstructed from two tracks and a two-photon resonance
/// (e.g. omega(782) -> pi+ pi- pi0).
///
/// posDau/negDau are drawn from already-selected FTracks partitions (sign and PID handled by the
/// normal track selection, same as any other track daughter use case), and neutralDau from an
/// already-built FPi0s/FEtas partition. This builder does the 3-body combinatorics, the
/// resonance's own mass window, and -- like TwoTrackResonanceBuilder -- copies the daughters'
/// bitmasks onto the candidate so systematic variations can be reselected downstream from the
/// candidate table alone, without rejoining the daughter tables.
template <modes::TrackTrackTwoPhotonResonance resoType>
class TrackTrackTwoPhotonResonanceBuilder
{
 public:
  TrackTrackTwoPhotonResonanceBuilder() = default;
  ~TrackTrackTwoPhotonResonanceBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void init(T1& confFilter, T2& confPosDauSelection, T3& confNegDauSelection, T4& confTable, T5& initContext)
  {
    mMassMin = confFilter.massMin.value;
    mMassMax = confFilter.massMax.value;
    mPtMin = confFilter.ptMin.value;
    mPtMax = confFilter.ptMax.value;
    mEtaMin = confFilter.etaMin.value;
    mEtaMax = confFilter.etaMax.value;
    mPhiMin = confFilter.phiMin.value;
    mPhiMax = confFilter.phiMax.value;

    mPosDauThreshold = confPosDauSelection.pidThres.value;
    mNegDauThreshold = confNegDauSelection.pidThres.value;

    if constexpr (modes::isEqual(resoType, modes::TrackTrackTwoPhotonResonance::kOmegaMeson)) {
      LOG(info) << "Initialize femto OmegaMeson builder...";
      mPosDauMass = o2::constants::physics::MassPiPlus;
      mNegDauMass = o2::constants::physics::MassPiMinus;
      mProduceOmegaMesons = utils::enableTable("FOmegaMesons_001", confTable.produceOmegaMesons.value, initContext);
      mProduceOmegaMesonMasks = utils::enableTable("FOmegaMesonMasks_001", confTable.produceOmegaMesonMasks.value, initContext);
    }

    if (mProduceOmegaMesons || mProduceOmegaMesonMasks) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fillResonances(T1 const& col, T2& resonanceProducts, T3& posDauPartition, T4& negDauPartition, T5& neutralDauPartition, T6& cache)
  {
    if (!mFillAnyTable) {
      return;
    }
    auto posDauSlice = posDauPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto negDauSlice = negDauPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto neutralDauSlice = neutralDauPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& [posDau, negDau, neutralDau] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posDauSlice, negDauSlice, neutralDauSlice))) {
      this->fillResonance(col, posDau, negDau, neutralDau, resonanceProducts);
    }
  }

 private:
  template <typename T1, typename T2, typename T3>
  void reconstructResonance(T1 const& posDau, T2 const& negDau, T3 const& neutralDau)
  {
    ROOT::Math::PtEtaPhiMVector vecPosDau{posDau.pt(), posDau.eta(), posDau.phi(), mPosDauMass};
    ROOT::Math::PtEtaPhiMVector vecNegDau{negDau.pt(), negDau.eta(), negDau.phi(), mNegDauMass};
    ROOT::Math::PtEtaPhiMVector vecNeutralDau{neutralDau.pt(), neutralDau.eta(), neutralDau.phi(), neutralDau.mass()};
    ROOT::Math::PtEtaPhiMVector vecResonance = vecPosDau + vecNegDau + vecNeutralDau;

    mPt = vecResonance.Pt();
    mEta = vecResonance.Eta();
    mPhi = RecoDecay::constrainAngle(vecResonance.Phi());
    mMass = vecResonance.M();
  }

  [[nodiscard]] bool checkFilters() const
  {
    return ((mMass > mMassMin && mMass < mMassMax) &&
            (mPt > mPtMin && mPt < mPtMax) &&
            (mEta > mEtaMin && mEta < mEtaMax) &&
            (mPhi >= mPhiMin && mPhi < mPhiMax));
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillResonance(T1 const& col, T2 const& posDau, T3 const& negDau, T4 const& neutralDau, T5& resonanceProducts)
  {
    reconstructResonance(posDau, negDau, neutralDau);
    if (!checkFilters()) {
      return;
    }

    bool posDauHasHighMomentum = posDau.p() > mPosDauThreshold;
    bool negDauHasHighMomentum = negDau.p() > mNegDauThreshold;

    if constexpr (modes::isEqual(resoType, modes::TrackTrackTwoPhotonResonance::kOmegaMeson)) {
      if (mProduceOmegaMesons) {
        resonanceProducts.producedOmegaMesons(col.globalIndex(),
                                              mPt,
                                              mEta,
                                              mPhi,
                                              mMass,
                                              posDau.globalIndex(),
                                              negDau.globalIndex(),
                                              neutralDau.globalIndex());
      }
      if (mProduceOmegaMesonMasks) {
        resonanceProducts.producedOmegaMesonMasks(posDau.mask(),
                                                  posDauHasHighMomentum,
                                                  negDau.mask(),
                                                  negDauHasHighMomentum,
                                                  neutralDau.maskDau());
      }
    }
  }

  // cached kinematics of the resonance
  float mPt = 0.f;
  float mEta = 0.f;
  float mPhi = 0.f;
  float mMass = 0.f;

  float mMassMin = 0.f;
  float mMassMax = 1.f;
  float mPtMin = 0.f;
  float mPtMax = 10.f;
  float mEtaMin = -0.9f;
  float mEtaMax = 0.9f;
  float mPhiMin = 0.f;
  float mPhiMax = o2::constants::math::TwoPI;

  float mPosDauMass = 0.f;
  float mNegDauMass = 0.f;

  float mPosDauThreshold = 0.f;
  float mNegDauThreshold = 0.f;

  bool mFillAnyTable = false;
  bool mProduceOmegaMesons = false;
  bool mProduceOmegaMesonMasks = false;
};

} // namespace o2::analysis::femto::trackstwophotonresonancebuilder
#endif // PWGCF_FEMTO_CORE_TRACKTRACKTWOPHOTONRESONANCEBUILDER_H_
