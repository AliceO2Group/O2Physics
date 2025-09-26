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

/// \file pairBuilder.h
/// \brief histogram manager for pair tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PAIRBUILDER_H_
#define PWGCF_FEMTO_CORE_PAIRBUILDER_H_

#include "PWGCF/Femto/Core/cascadeHistManager.h"
#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/kinkHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/pairProcessHelpers.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/twoTrackResonanceHistManager.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/HistogramRegistry.h"

#include <map>
#include <random>
#include <vector>

namespace o2::analysis::femto
{
namespace pairbuilder
{

template <
  const char* prefixTrack1,
  const char* prefixTrack2,
  const char* prefixSe,
  const char* prefixMe,
  const char* prefixCprSe,
  const char* prefixCprMe,
  modes::Mode mode>
class PairTrackTrackBuilder
{
 public:
  PairTrackTrackBuilder() = default;

  template <typename T1,
            typename T2,
            typename T3,
            typename T4,
            typename T5,
            typename T6,
            typename T7,
            typename T8,
            typename T9>
  void init(o2::framework::HistogramRegistry* registry,
            T1& confTrackSelection1,
            T2& confTrackSelection2,
            T3& confCpr,
            T4& confMixing,
            std::map<T5, std::vector<framework::AxisSpec>>& colHistSpec,
            std::map<T6, std::vector<framework::AxisSpec>>& trackHistSpec1,
            std::map<T7, std::vector<framework::AxisSpec>>& trackHistSpec2,
            std::map<T8, std::vector<framework::AxisSpec>>& pairHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>>& cprHistSpec)
  {

    // check if correlate the same tracks or not
    mSameSpecies = confMixing.sameSpecies.value;

    mColHistManager.init(registry, colHistSpec);
    mPairHistManagerSe.init(registry, pairHistSpec);
    mPairHistManagerMe.init(registry, pairHistSpec);

    if (mSameSpecies) {
      mTrackHistManager1.init(registry, trackHistSpec1);

      mPairHistManagerSe.setMass(confTrackSelection1.pdgCode.value, confTrackSelection1.pdgCode.value);
      mPairHistManagerSe.setCharge(confTrackSelection1.charge.value, confTrackSelection1.charge.value);
      mCprSe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection1.charge.value, confTrackSelection1.charge.value, confCpr.on.value);

      mPairHistManagerMe.setMass(confTrackSelection1.pdgCode.value, confTrackSelection1.pdgCode.value);
      mPairHistManagerMe.setCharge(confTrackSelection1.charge.value, confTrackSelection1.charge.value);
      mCprMe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection1.charge.value, confTrackSelection1.charge.value, confCpr.on.value);
    } else {
      mTrackHistManager1.init(registry, trackHistSpec1);
      mTrackHistManager2.init(registry, trackHistSpec2);

      mPairHistManagerSe.setMass(confTrackSelection1.pdgCode.value, confTrackSelection2.pdgCode.value);
      mPairHistManagerSe.setCharge(confTrackSelection1.charge.value, confTrackSelection2.charge.value);
      mCprSe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection1.charge.value, confTrackSelection2.charge.value, confCpr.on.value);

      mPairHistManagerMe.setMass(confTrackSelection1.pdgCode.value, confTrackSelection2.pdgCode.value);
      mPairHistManagerMe.setCharge(confTrackSelection1.charge.value, confTrackSelection2.charge.value);
      mCprMe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection1.charge.value, confTrackSelection2.charge.value, confCpr.on.value);
    }

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPoliciy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;

    // setup rng if necessary
    if (confMixing.seed.value >= 0) {
      uint64_t randomSeed = 0;
      mMixIdenticalParticles = true;
      if (confMixing.seed.value == 0) {
        randomSeed = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      } else {
        randomSeed = static_cast<uint64_t>(confMixing.seed.value);
      }
      mRng = std::mt19937(randomSeed);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void processSameEvent(T1 const& col, T2& /*trackTable*/, T3& partition1, T4& partition2, T5& cache)
  {
    if (mSameSpecies) {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0) {
        return;
      }
      mColHistManager.fill(col);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent(trackSlice1, mTrackHistManager1, mPairHistManagerSe, mCprSe, mRng, mMixIdenticalParticles);
    } else {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
      auto trackSlice2 = partition2->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice2.size() == 0) {
        return;
      }
      mColHistManager.fill(col);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent(trackSlice1, trackSlice2, mTrackHistManager1, mTrackHistManager2, mPairHistManagerSe, mCprSe, mPc);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& /*trackTable*/, T3& partition1, T4& partition2, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {

    if (mSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent(cols, partition1, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent(cols, partition1, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent(cols, partition1, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent(cols, partition1, partition2, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent(cols, partition1, partition2, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent(cols, partition1, partition2, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }

 private:
  colhistmanager::CollisionHistManager<mode> mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack1, mode> mTrackHistManager1;
  trackhistmanager::TrackHistManager<prefixTrack2, mode> mTrackHistManager2;
  pairhistmanager::PairHistManager<prefixSe, mode> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, mode> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackTrack<prefixCprSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackTrack<prefixCprMe> mCprMe;
  paircleaner::TrackTrackPairCleaner mPc;
  std::mt19937 mRng;
  pairhistmanager::MixingPoliciy mMixingPolicy = pairhistmanager::MixingPoliciy::kVtxMult;
  bool mSameSpecies = false;
  int mMixingDepth = 5;
  bool mMixIdenticalParticles = false;
};

template <
  const char* prefixTrack,
  const char* prefixV0,
  const char* prefixPosDau,
  const char* prefixNegDau,
  const char* prefixSe,
  const char* prefixMe,
  const char* prefixCprSe,
  const char* prefixCprMe,
  modes::Mode mode,
  modes::V0 v0Type>
class PairTrackV0Builder
{
 public:
  PairTrackV0Builder() = default;

  template <typename T1,
            typename T2,
            typename T3,
            typename T4,
            typename T5,
            typename T6,
            typename T7,
            typename T8,
            typename T9,
            typename T10,
            typename T11>
  void init(o2::framework::HistogramRegistry* registry,
            T1& confTrackSelection,
            T2& confV0Selection,
            T3& confCpr,
            T4& confMixing,
            std::map<T5, std::vector<framework::AxisSpec>>& colHistSpec,
            std::map<T6, std::vector<framework::AxisSpec>>& trackHistSpec,
            std::map<T7, std::vector<framework::AxisSpec>>& v0HistSpec,
            std::map<T8, std::vector<framework::AxisSpec>>& posDauHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>>& negDauHistSpec,
            std::map<T10, std::vector<framework::AxisSpec>>& pairHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>>& cprHistSpec)
  {
    mColHistManager.init(registry, colHistSpec);

    mTrackHistManager.init(registry, trackHistSpec);
    mV0HistManager.init(registry, v0HistSpec, posDauHistSpec, negDauHistSpec);

    mPairHistManagerSe.init(registry, pairHistSpec);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCode.value, confV0Selection.pdgCode.value);
    mPairHistManagerSe.setCharge(confTrackSelection.charge.value, 1);
    mCprSe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confCpr.on.value);

    mPairHistManagerMe.init(registry, pairHistSpec);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCode.value, confV0Selection.pdgCode.value);
    mPairHistManagerMe.setCharge(confTrackSelection.charge.value, 1);
    mCprMe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confCpr.on.value);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPoliciy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*v0table*/, T5& v0Partition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    auto v0Slice = v0Partition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || v0Slice.size() == 0) {
      return;
    }
    mColHistManager.fill(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent(trackSlice, v0Slice, trackTable, mTrackHistManager, mV0HistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& v0Partition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, v0Partition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager<mode> mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack, mode> mTrackHistManager;
  v0histmanager::V0HistManager<prefixV0, prefixPosDau, prefixNegDau, mode, v0Type> mV0HistManager;
  pairhistmanager::PairHistManager<prefixSe, mode> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, mode> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackV0<prefixCprSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackV0<prefixCprMe> mCprMe;
  paircleaner::TrackV0PairCleaner mPc;
  pairhistmanager::MixingPoliciy mMixingPolicy = pairhistmanager::MixingPoliciy::kVtxMult;
  int mMixingDepth = 5;
};

template <
  const char* prefixTrack,
  const char* prefixResonance,
  const char* prefixPosDau,
  const char* prefixNegDau,
  const char* prefixSe,
  const char* prefixMe,
  const char* prefixCprSe,
  const char* prefixCprMe,
  modes::Mode mode,
  modes::TwoTrackResonance resonanceType>
class PairTrackTwoTrackResonanceBuilder
{
 public:
  PairTrackTwoTrackResonanceBuilder() = default;

  template <typename T1,
            typename T2,
            typename T3,
            typename T4,
            typename T5,
            typename T6,
            typename T7,
            typename T8,
            typename T9,
            typename T10,
            typename T11>
  void init(o2::framework::HistogramRegistry* registry,
            T1& confTrackSelection,
            T2& confResonanceSelection,
            T3& confCpr,
            T4& confMixing,
            std::map<T5, std::vector<framework::AxisSpec>>& colHistSpec,
            std::map<T6, std::vector<framework::AxisSpec>>& trackHistSpec,
            std::map<T7, std::vector<framework::AxisSpec>>& resonanceHistSpec,
            std::map<T8, std::vector<framework::AxisSpec>>& posDauHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>>& negDauHistSpec,
            std::map<T10, std::vector<framework::AxisSpec>>& pairHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>>& cprHistSpec)
  {
    mColHistManager.init(registry, colHistSpec);

    mTrackHistManager.init(registry, trackHistSpec);
    mResonanceHistManager.init(registry, resonanceHistSpec, posDauHistSpec, negDauHistSpec);

    mPairHistManagerSe.init(registry, pairHistSpec);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCode.value, confResonanceSelection.pdgCode.value);
    mPairHistManagerSe.setCharge(confTrackSelection.charge.value, 1);
    mCprSe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confCpr.on.value);

    mPairHistManagerMe.init(registry, pairHistSpec);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCode.value, confResonanceSelection.pdgCode.value);
    mPairHistManagerMe.setCharge(confTrackSelection.charge.value, 1);
    mCprMe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confCpr.on.value);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPoliciy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*resonanceTable*/, T5& resonancePartition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    auto v0Slice = resonancePartition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || v0Slice.size() == 0) {
      return;
    }
    mColHistManager.fill(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent(trackSlice, v0Slice, trackTable, mTrackHistManager, mResonanceHistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& resonancePartition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, resonancePartition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, resonancePartition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, resonancePartition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager<mode> mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack, mode> mTrackHistManager;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<prefixResonance, prefixPosDau, prefixNegDau, mode, resonanceType> mResonanceHistManager;
  pairhistmanager::PairHistManager<prefixSe, mode> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, mode> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackV0<prefixCprSe> mCprSe; // cpr for twotrackresonances and v0 work the same way
  closepairrejection::ClosePairRejectionTrackV0<prefixCprMe> mCprMe; // cpr for twotrackresonances and v0 work the same way
  paircleaner::TrackV0PairCleaner mPc;                               // pc for twotrackresonances and v0 work the same way
  pairhistmanager::MixingPoliciy mMixingPolicy = pairhistmanager::MixingPoliciy::kVtxMult;
  int mMixingDepth = 5;
};

template <
  const char* prefixTrack,
  const char* prefixKink,
  const char* prefixChaDau,
  const char* prefixSe,
  const char* prefixMe,
  const char* prefixCprSe,
  const char* prefixCprMe,
  modes::Mode mode,
  modes::Kink kinkType>
class PairTrackKinkBuilder
{
 public:
  PairTrackKinkBuilder() = default;

  template <typename T1,
            typename T2,
            typename T3,
            typename T4,
            typename T5,
            typename T6,
            typename T7,
            typename T8,
            typename T9,
            typename T10>
  void init(o2::framework::HistogramRegistry* registry,
            T1& confTrackSelection,
            T2& confKinkSelection,
            T3& confCpr,
            T4& confMixing,
            std::map<T5, std::vector<framework::AxisSpec>>& colHistSpec,
            std::map<T6, std::vector<framework::AxisSpec>>& trackHistSpec,
            std::map<T7, std::vector<framework::AxisSpec>>& kinkHistSpec,
            std::map<T8, std::vector<framework::AxisSpec>>& chaDauHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>>& pairHistSpec,
            std::map<T10, std::vector<framework::AxisSpec>>& cprHistSpec)
  {
    mColHistManager.init(registry, colHistSpec);

    mTrackHistManager.init(registry, trackHistSpec);
    mKinkHistManager.init(registry, kinkHistSpec, chaDauHistSpec);

    mPairHistManagerSe.init(registry, pairHistSpec);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCode.value, confKinkSelection.pdgCode.value);
    mPairHistManagerSe.setCharge(confTrackSelection.charge.value, confKinkSelection.sign.value);
    mCprSe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confKinkSelection.sign.value, confCpr.on.value);

    mPairHistManagerMe.init(registry, pairHistSpec);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCode.value, confKinkSelection.pdgCode.value);
    mPairHistManagerMe.setCharge(confTrackSelection.charge.value, confKinkSelection.sign.value);
    mCprMe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confKinkSelection.sign.value, confCpr.on.value);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPoliciy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*kinktable*/, T5& kinkPartition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    auto kinkSlice = kinkPartition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || kinkSlice.size() == 0) {
      return;
    }
    mColHistManager.fill(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent(trackSlice, kinkSlice, trackTable, mTrackHistManager, mKinkHistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& kinkPartition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, kinkPartition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, kinkPartition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, kinkPartition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager<mode> mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack, mode> mTrackHistManager;
  kinkhistmanager::KinkHistManager<prefixKink, prefixChaDau, mode, kinkType> mKinkHistManager;
  pairhistmanager::PairHistManager<prefixSe, mode> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, mode> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackKink<prefixCprSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackKink<prefixCprMe> mCprMe;
  paircleaner::TrackKinkPairCleaner mPc;
  pairhistmanager::MixingPoliciy mMixingPolicy = pairhistmanager::MixingPoliciy::kVtxMult;
  int mMixingDepth = 5;
};

template <
  const char* prefixTrack,
  const char* prefixCascade,
  const char* prefixBachelor,
  const char* prefixPosDau,
  const char* prefixNegDau,
  const char* prefixSe,
  const char* prefixMe,
  const char* prefixCprSe,
  const char* prefixCprMe,
  modes::Mode mode,
  modes::Cascade cascadeType>
class PairTrackCascadeBuilder
{
 public:
  PairTrackCascadeBuilder() = default;

  template <typename T1,
            typename T2,
            typename T3,
            typename T4,
            typename T5,
            typename T6,
            typename T7,
            typename T8,
            typename T9,
            typename T10,
            typename T11,
            typename T12>
  void init(o2::framework::HistogramRegistry* registry,
            T1& confTrackSelection,
            T2& confCascadeSelection,
            T3& confCpr,
            T4& confMixing,
            std::map<T5, std::vector<framework::AxisSpec>>& colHistSpec,
            std::map<T6, std::vector<framework::AxisSpec>>& trackHistSpec,
            std::map<T7, std::vector<framework::AxisSpec>>& cascadeHistSpec,
            std::map<T8, std::vector<framework::AxisSpec>>& bachelorHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>>& posDauHistSpec,
            std::map<T10, std::vector<framework::AxisSpec>>& negDauHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>>& pairHistSpec,
            std::map<T12, std::vector<framework::AxisSpec>>& cprHistSpec)
  {
    mColHistManager.init(registry, colHistSpec);

    mTrackHistManager.init(registry, trackHistSpec);
    mCascadeHistManager.init(registry, cascadeHistSpec, bachelorHistSpec, posDauHistSpec, negDauHistSpec);

    mPairHistManagerSe.init(registry, pairHistSpec);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCode.value, confCascadeSelection.pdgCode.value);
    mPairHistManagerSe.setCharge(confTrackSelection.charge.value, confCascadeSelection.sign.value);
    mCprSe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confCascadeSelection.sign.value, confCpr.on.value);

    mPairHistManagerMe.init(registry, pairHistSpec);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCode.value, confCascadeSelection.pdgCode.value);
    mPairHistManagerMe.setCharge(confTrackSelection.charge.value, confCascadeSelection.sign.value);
    mCprMe.init(registry, cprHistSpec, confCpr.detaMax.value, confCpr.dphistarMax.value, confTrackSelection.charge.value, confCascadeSelection.sign.value, confCpr.on.value);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPoliciy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*cascadeTable*/, T5& v0Partition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    auto v0Slice = v0Partition->sliceByCached(o2::aod::femtobase::stored::collisionId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || v0Slice.size() == 0) {
      return;
    }
    mColHistManager.fill(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent(trackSlice, v0Slice, trackTable, mTrackHistManager, mCascadeHistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& v0Partition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, v0Partition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager<mode> mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack, mode> mTrackHistManager;
  cascadehistmanager::CascadeHistManager<prefixCascade, prefixBachelor, prefixPosDau, prefixNegDau, mode, cascadeType> mCascadeHistManager;
  pairhistmanager::PairHistManager<prefixSe, mode> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, mode> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackCascade<prefixCprSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackCascade<prefixCprMe> mCprMe;
  paircleaner::TrackCascadePairCleaner mPc;
  pairhistmanager::MixingPoliciy mMixingPolicy = pairhistmanager::MixingPoliciy::kVtxMult;
  int mMixingDepth = 5;
};

} // namespace pairbuilder
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_PAIRBUILDER_H_
