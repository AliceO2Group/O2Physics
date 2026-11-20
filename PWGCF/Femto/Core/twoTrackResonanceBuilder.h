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

#include <cstdint>
#include <string>

namespace o2::analysis::femto::twotrackresonancebuilder
{
template <auto& Prefix>
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

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TWOTRACKRESONANCE_DEFAULT_SELECTION(defaultPdgCode, defaultMassMin, defaultMassMax)                                                \
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", (defaultPdgCode), "Resonance PDG code. Set sign to minus 1 for antiparticle"}; \
  o2::framework::Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};                                                                    \
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};                                                                    \
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};                                                               \
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};                                                                \
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};                                                                 \
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};                                    \
  o2::framework::Configurable<float> massMin{"massMin", (defaultMassMin), "Minimum invariant mass for Resonance"};                         \
  o2::framework::Configurable<float> massMax{"massMax", (defaultMassMax), "Maximum invariant mass for Resonance"};

// phi -> K+ K-, rho0 -> pi+ pi-: both daughters are the same species, so a single pair of mask
// thresholds is enough, applied to both the positive and negative daughter alike
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TWOTRACKRESONANCE_SAME_SPECIES_MASK                                                                                                                               \
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> dauMaskBelowThres{"dauMaskBelowThres", 0x10u, "Bitmask for both daughters below threshold"}; \
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> dauMaskAboveThres{"dauMaskAboveThres", 0x8u, "Bitmask for both daughters above threshold"};

struct ConfPhiSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("PhiSelection");
  TWOTRACKRESONANCE_DEFAULT_SELECTION(333, 0.95f, 1.05f)
  TWOTRACKRESONANCE_SAME_SPECIES_MASK
  o2::framework::Configurable<int> sign{"sign", 1, "Dummy value for compatibility"};
};

struct ConfRho0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Rho0Selection");
  TWOTRACKRESONANCE_DEFAULT_SELECTION(113, 0.7f, 0.84f)
  TWOTRACKRESONANCE_SAME_SPECIES_MASK
  o2::framework::Configurable<int> sign{"sign", 1, "Dummy value for compatibility"};
};

#undef TWOTRACKRESONANCE_SAME_SPECIES_MASK

// K*0 -> K+ pi-, K*0bar -> K- pi+: the physical pos/neg daughter's species flips with the
// candidate's charge sign, so the masks are named by species (kaon/pion) rather than by charge --
// MAKE_RESONANCE_1_PARTITON resolves which stored column (maskPosDau/maskNegDau) each applies to
// from the candidate's own signedPt, so the user never has to reason about pos/neg <-> kaon/pion
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TWOTRACKRESONANCE_KSTAR_MASK                                                                                                                                                 \
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> kaonDauMaskBelowThres{"kaonDauMaskBelowThres", 0x10u, "Bitmask for the kaon daughter below threshold"}; \
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> kaonDauMaskAboveThres{"kaonDauMaskAboveThres", 0x8u, "Bitmask for the kaon daughter above threshold"};  \
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> pionDauMaskBelowThres{"pionDauMaskBelowThres", 0x2u, "Bitmask for the pion daughter below threshold"};  \
  o2::framework::Configurable<o2::analysis::femto::datatypes::TrackMaskType> pionDauMaskAboveThres{"pionDauMaskAboveThres", 0x1u, "Bitmask for the pion daughter above threshold"};

struct ConfKstar0Selection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Kstar0Selection");
  o2::framework::Configurable<int> sign{"sign", 1, "Sign (+1 for Kstar0 and -1 for Kstar0Bar) "};
  TWOTRACKRESONANCE_DEFAULT_SELECTION(313, 0.8f, 1.0f)
  TWOTRACKRESONANCE_KSTAR_MASK
};

#undef TWOTRACKRESONANCE_DEFAULT_SELECTION
#undef TWOTRACKRESONANCE_KSTAR_MASK

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
  o2::framework::Configurable<int> produceKstar0s{"produceKstar0s", -1, "Produce Kstar0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceKstar0Masks{"produceKstar0Masks", -1, "Produce Kstar0Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceRho0s{"produceRho0s", -1, "Produce Rho0s (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceRho0Masks{"produceRho0Masks", -1, "Produce Rho0Masks (-1: auto; 0 off; 1 on)"};
};

template <modes::TwoTrackResonance resoType>
class TwoTrackResonanceBuilder
{
 public:
  TwoTrackResonanceBuilder() = default;
  ~TwoTrackResonanceBuilder() = default;

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

    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kPhi)) {
      LOG(info) << "Initialize femto Phi builder...";
      mPosDaughterMass = o2::constants::physics::MassKPlus;
      mNegDaughterMass = o2::constants::physics::MassKMinus;
      mProducePhis = utils::enableTable("FPhis_001", confTable.producePhis.value, initContext);
      mProducePhiMasks = utils::enableTable("FPhiMasks_001", confTable.producePhiMasks.value, initContext);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0) || modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
      if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0)) {
        LOG(info) << "Initialize femto Kstar0 builder...";
        mPosDaughterMass = o2::constants::physics::MassKPlus;
        mNegDaughterMass = o2::constants::physics::MassPiMinus;
      }
      if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
        LOG(info) << "Initialize femto Kstar0Bar builder...";
        mPosDaughterMass = o2::constants::physics::MassPiPlus;
        mNegDaughterMass = o2::constants::physics::MassKMinus;
      }
      mProduceKstar0s = utils::enableTable("FKstar0s_001", confTable.produceKstar0s.value, initContext);
      mProduceKstar0Masks = utils::enableTable("FKstar0Masks_001", confTable.produceKstar0Masks.value, initContext);
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kRho0)) {
      LOG(info) << "Initialize femto Rho0 builder...";
      mPosDaughterMass = o2::constants::physics::MassPiPlus;
      mNegDaughterMass = o2::constants::physics::MassPiMinus;
      mProduceRho0s = utils::enableTable("FRho0s_001", confTable.produceRho0s.value, initContext);
      mProduceRho0Masks = utils::enableTable("FRho0Masks_001", confTable.produceRho0Masks.value, initContext);
    }

    if (mProducePhis || mProducePhiMasks || mProduceKstar0s || mProduceKstar0Masks || mProduceRho0s || mProduceRho0Masks) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void fillResonances(T1 const& col, T2& resonanceProducts, T3& posDaughterPartition, T4& negDaughterPartition, T5 const& /*trackTable*/, T6& cache)
  {
    if (!mFillAnyTable) {
      return;
    }
    auto posDaughterSlice = posDaughterPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto negDaughterSlice = negDaughterPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    for (auto const& [posDaughter, negDaughter] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posDaughterSlice, negDaughterSlice))) {
      this->fillResonance(col, posDaughter, negDaughter, resonanceProducts);
    }
  }

 private:
  template <typename T1, typename T2>
  void reconstructResonance(T1 const& posDaughter, T2 const& negDaughter)
  {
    ROOT::Math::PtEtaPhiMVector vecPosDaughter{posDaughter.pt(), posDaughter.eta(), posDaughter.phi(), mPosDaughterMass};
    ROOT::Math::PtEtaPhiMVector vecNegDaughter{negDaughter.pt(), negDaughter.eta(), negDaughter.phi(), mNegDaughterMass};
    ROOT::Math::PtEtaPhiMVector vecResonance = vecPosDaughter + vecNegDaughter;

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

  template <typename T1, typename T2, typename T3, typename T4>
  void fillResonance(T1 const& col, T2 const& posDaughter, T3 const& negDaughter, T4& resonanceProducts)
  {
    reconstructResonance(posDaughter, negDaughter);
    if (!checkFilters()) {
      return;
    }

    bool posDauHasHighMomentum = posDaughter.p() > mPosDauThreshold;
    bool negDauHasHighMomentum = negDaughter.p() > mNegDauThreshold;

    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kRho0)) {
      if (mProduceRho0s) {
        resonanceProducts.producedRhos(col.globalIndex(),
                                       mPt,
                                       mEta,
                                       mPhi,
                                       mMass,
                                       posDaughter.globalIndex(),
                                       negDaughter.globalIndex());
      }
      if (mProduceRho0Masks) {
        resonanceProducts.producedRhoMasks(posDaughter.mask(),
                                           posDauHasHighMomentum,
                                           negDaughter.mask(),
                                           negDauHasHighMomentum);
      }
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kPhi)) {
      if (mProducePhis) {
        resonanceProducts.producedPhis(col.globalIndex(),
                                       mPt,
                                       mEta,
                                       mPhi,
                                       mMass,
                                       posDaughter.globalIndex(),
                                       negDaughter.globalIndex());
      }
      if (mProducePhiMasks) {
        resonanceProducts.producedPhiMasks(posDaughter.mask(),
                                           posDauHasHighMomentum,
                                           negDaughter.mask(),
                                           negDauHasHighMomentum);
      }
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0)) {
      if (mProduceKstar0s) {
        resonanceProducts.producedKstars(col.globalIndex(),
                                         mPt,
                                         mEta,
                                         mPhi,
                                         mMass,
                                         posDaughter.globalIndex(),
                                         negDaughter.globalIndex());
      }
      if (mProduceKstar0Masks) {
        resonanceProducts.producedKstarMasks(posDaughter.mask(),
                                             posDauHasHighMomentum,
                                             negDaughter.mask(),
                                             negDauHasHighMomentum);
      }
    }
    if constexpr (modes::isEqual(resoType, modes::TwoTrackResonance::kKstar0Bar)) {
      if (mProduceKstar0s) {
        resonanceProducts.producedKstars(col.globalIndex(),
                                         -1.f * mPt,
                                         mEta,
                                         mPhi,
                                         mMass,
                                         posDaughter.globalIndex(),
                                         negDaughter.globalIndex());
      }
      if (mProduceKstar0Masks) {
        resonanceProducts.producedKstarMasks(posDaughter.mask(),
                                             posDauHasHighMomentum,
                                             negDaughter.mask(),
                                             negDauHasHighMomentum);
      }
    }
  }

  // cached kinamtics of the resonance
  float mPt = 0;
  float mEta = 0;
  float mPhi = 0;
  float mMass = 0;

  float mPosDaughterMass = 0.f;
  float mNegDaughterMass = 0.f;

  float mPosDauThreshold = 0.f;
  float mNegDauThreshold = 0.f;

  float mMassMin = 0.f;
  float mMassMax = 0.f;
  float mPtMin = 0.f;
  float mPtMax = 0.f;
  float mEtaMin = 0.f;
  float mEtaMax = 0.f;
  float mPhiMin = 0.f;
  float mPhiMax = 0.f;

  bool mFillAnyTable = false;
  bool mProducePhis = false;
  bool mProducePhiMasks = false;
  bool mProduceKstar0s = false;
  bool mProduceKstar0Masks = false;
  bool mProduceRho0s = false;
  bool mProduceRho0Masks = false;

}; // namespace twotrackresonancebuilder

struct ConfTwoTrackResonanceTablesDerivedToDerived : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TwoTrackResonanceTables");
  o2::framework::Configurable<int> limitRho0{"limitRho0", 1, "Require at least this many Rho0s in the collision. Set to 0 to skip Rho0s entirely (not written to output)."};
  o2::framework::Configurable<int> limitPhi{"limitPhi", 0, "Require at least this many Phis in the collision. Set to 0 to skip Phis entirely (not written to output)."};
};

struct TwoTrackResonanceBuilderDerivedToDerivedProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::StoredFRho0s> producedRho0s;
  o2::framework::Produces<o2::aod::StoredFRho0Masks> producedRho0Masks;
  o2::framework::Produces<o2::aod::StoredFPhis> producedPhis;
  o2::framework::Produces<o2::aod::StoredFPhiMasks> producedPhiMasks;
};

class TwoTrackResonanceBuilderDerivedToDerived
{
 public:
  template <typename T>
  void init(T& config)
  {
    mLimitRho0 = config.limitRho0.value;
    mLimitPhi = config.limitPhi.value;

    if (mLimitRho0 < 0 || mLimitPhi < 0) {
      LOG(fatal) << "Two-track resonance limits must be non-negative (got " << mLimitRho0 << " and " << mLimitPhi << "). Breaking...";
    }
    if (mLimitRho0 == 0 && mLimitPhi == 0) {
      LOG(fatal) << "Both Rho0 limit and Phi limit are 0. Breaking...";
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewRho0s(T1 const& col, T2 const& /*rho0Table*/, T3& partitionRho0, T4& cache)
  {
    if (mLimitRho0 == 0) { // Rho0s disabled, cannot reject on them
      return false;
    }
    auto rho0Slice = partitionRho0->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return rho0Slice.size() < static_cast<int64_t>(mLimitRho0);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  bool collisionHasTooFewPhis(T1 const& col, T2 const& /*phiTable*/, T3& partitionPhi, T4& cache)
  {
    if (mLimitPhi == 0) { // Phis disabled, cannot reject on them
      return false;
    }
    auto phiSlice = partitionPhi->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return phiSlice.size() < static_cast<int64_t>(mLimitPhi);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processRho0s(T1 const& col, T2 const& /*rho0Table*/, T3 const& oldTrackTable, T4& partitionRho0, T5& trackBuilder, T6& cache, T7& newRho0Table, T8& newTrackTable, T9& newCollisionTable)
  {
    if (mLimitRho0 == 0) { // Rho0s disabled
      return;
    }

    auto rho0Slice = partitionRho0->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& rho0 : rho0Slice) {
      auto posDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(rho0.posDauId(), oldTrackTable));
      auto negDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(rho0.negDauId(), oldTrackTable));

      int64_t posDaughterIndex = trackBuilder.getDaughterIndex(posDaughter, newTrackTable, newCollisionTable);
      int64_t negDaughterIndex = trackBuilder.getDaughterIndex(negDaughter, newTrackTable, newCollisionTable);

      newRho0Table.producedRho0s(newCollisionTable.producedCollision.lastIndex(),
                                 rho0.pt(),
                                 rho0.eta(),
                                 rho0.phi(),
                                 rho0.mass(),
                                 posDaughterIndex,
                                 negDaughterIndex);
      newRho0Table.producedRho0Masks(rho0.maskPosDau(),
                                     rho0.posDauHasHighMomentum(),
                                     rho0.maskNegDau(),
                                     rho0.negDauHasHighMomentum());
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processPhis(T1 const& col, T2 const& /*phiTable*/, T3 const& oldTrackTable, T4& partitionPhi, T5& trackBuilder, T6& cache, T7& newPhiTable, T8& newTrackTable, T9& newCollisionTable)
  {
    if (mLimitPhi == 0) { // Phis disabled
      return;
    }

    auto phiSlice = partitionPhi->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& phi : phiSlice) {
      auto posDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(phi.posDauId(), oldTrackTable));
      auto negDaughter = oldTrackTable.rawIteratorAt(utils::daughterRow(phi.negDauId(), oldTrackTable));

      int64_t posDaughterIndex = trackBuilder.getDaughterIndex(posDaughter, newTrackTable, newCollisionTable);
      int64_t negDaughterIndex = trackBuilder.getDaughterIndex(negDaughter, newTrackTable, newCollisionTable);

      newPhiTable.producedPhis(newCollisionTable.producedCollision.lastIndex(),
                               phi.pt(),
                               phi.eta(),
                               phi.phi(),
                               phi.mass(),
                               posDaughterIndex,
                               negDaughterIndex);
      newPhiTable.producedPhiMasks(phi.maskPosDau(),
                                   phi.posDauHasHighMomentum(),
                                   phi.maskNegDau(),
                                   phi.negDauHasHighMomentum());
    }
  }

 private:
  int mLimitRho0 = 0;
  int mLimitPhi = 0;
};

} // namespace o2::analysis::femto::twotrackresonancebuilder
#endif // PWGCF_FEMTO_CORE_TWOTRACKRESONANCEBUILDER_H_
