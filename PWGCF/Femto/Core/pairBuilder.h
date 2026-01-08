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
#include "Framework/HistogramSpec.h"

#include "fairlogger/Logger.h"

#include <chrono>
#include <cstdint>
#include <map>
#include <random>
#include <vector>

namespace o2::analysis::femto
{
namespace pairbuilder
{

template <const char* prefixTrack1,
          const char* prefixTrack2,
          const char* prefixSe,
          const char* prefixMe,
          const char* prefixCprSe,
          const char* prefixCprMe>
class PairTrackTrackBuilder
{
 public:
  PairTrackTrackBuilder() = default;
  ~PairTrackTrackBuilder() = default;

  template <modes::Mode mode,
            typename T1,
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
            T1 const& confTrackSelection1,
            T2 const& confTrackSelection2,
            T3 const& confCpr,
            T4 const& confMixing,
            T5 const& confPairBinning,
            T6 const& confPairCuts,
            std::map<T7, std::vector<framework::AxisSpec>> const& colHistSpec,
            std::map<T8, std::vector<framework::AxisSpec>> const& trackHistSpec1,
            std::map<T9, std::vector<framework::AxisSpec>> const& trackHistSpec2,
            std::map<T10, std::vector<framework::AxisSpec>> const& pairHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>> const& cprHistSpec)
  {

    // check if correlate the same tracks or not
    mSameSpecies = confMixing.sameSpecies.value;

    mColHistManager.template init<mode>(registry, colHistSpec);
    mPairHistManagerSe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerMe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPc.template init<mode>(confPairCuts);

    if (mSameSpecies) {
      mTrackHistManager1.template init<mode>(registry, trackHistSpec1, confTrackSelection1);

      mPairHistManagerSe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value);
      mPairHistManagerSe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);
      mCprSe.init(registry, cprHistSpec, confCpr, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);

      mPairHistManagerMe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value);
      mPairHistManagerMe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);
      mCprMe.init(registry, cprHistSpec, confCpr, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);
    } else {
      mTrackHistManager1.template init<mode>(registry, trackHistSpec1, confTrackSelection1);
      mTrackHistManager2.template init<mode>(registry, trackHistSpec2, confTrackSelection2);

      mPairHistManagerSe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection2.pdgCodeAbs.value);
      mPairHistManagerSe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value);
      mCprSe.init(registry, cprHistSpec, confCpr, confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value);

      mPairHistManagerMe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection2.pdgCodeAbs.value);
      mPairHistManagerMe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value);
      mCprMe.init(registry, cprHistSpec, confCpr, confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value);
    }

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPolicy>(confMixing.policy.value);
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

  // data
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processSameEvent(T1 const& col, T2& trackTable, T3& partition1, T4& partition2, T5& cache)
  {
    if (mSameSpecies) {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent<mode>(trackSlice1, trackTable, col, mTrackHistManager1, mPairHistManagerSe, mCprSe, mPc, mRng, mMixIdenticalParticles);
    } else {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice2 = partition2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice2.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent<mode>(trackSlice1, trackSlice2, trackTable, col, mTrackHistManager1, mTrackHistManager2, mPairHistManagerSe, mCprSe, mPc);
    }
  }

  // mc
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processSameEvent(T1 const& col, T2 const& mcCols, T3& trackTable, T4& partition1, T5& partition2, T6 const& mcParticles, T7 const& mcMothers, T8 const& mcPartonicMothers, T9& cache)
  {
    if (mSameSpecies) {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col, mcCols);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent<mode>(trackSlice1, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mTrackHistManager1, mPairHistManagerSe, mCprSe, mPc, mRng, mMixIdenticalParticles);
    } else {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice2 = partition2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice2.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col, mcCols);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent<mode>(trackSlice1, trackSlice2, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mTrackHistManager1, mTrackHistManager2, mPairHistManagerSe, mCprSe, mPc);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& partition1, T4& partition2, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {

    if (mSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void processMixedEvent(T1 const& cols, T2 const& mcCols, T3& trackTable, T4& partition1, T5& partition2, T6 const& mcParticles, T7& cache, T8& binsVtxMult, T9& binsVtxCent, T10& binsVtxMultCent)
  {
    if (mSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, trackTable, mcParticles, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, trackTable, mcParticles, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, trackTable, mcParticles, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition2, trackTable, mcParticles, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition2, trackTable, mcParticles, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition2, trackTable, mcParticles, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }

 private:
  colhistmanager::CollisionHistManager mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack1> mTrackHistManager1;
  trackhistmanager::TrackHistManager<prefixTrack2> mTrackHistManager2;
  pairhistmanager::PairHistManager<prefixSe, modes::Particle::kTrack, modes::Particle::kTrack> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, modes::Particle::kTrack, modes::Particle::kTrack> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackTrack<prefixCprSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackTrack<prefixCprMe> mCprMe;
  paircleaner::TrackTrackPairCleaner mPc;
  std::mt19937 mRng;
  pairhistmanager::MixingPolicy mMixingPolicy = pairhistmanager::MixingPolicy::kVtxMult;
  bool mSameSpecies = false;
  int mMixingDepth = 5;
  bool mMixIdenticalParticles = false;
};

template <const char* prefixV01,
          const char* prefixPosDau1,
          const char* prefixNegDau1,
          const char* prefixV02,
          const char* prefixPosDau2,
          const char* prefixNegDau2,
          const char* prefixSe,
          const char* prefixMe,
          const char* prefixCprPosSe,
          const char* prefixCprNegSe,
          const char* prefixCprPosMe,
          const char* prefixCprNegMe,
          modes::V0 v0Type1,
          modes::V0 v0Type2>
class PairV0V0Builder
{
 public:
  PairV0V0Builder() = default;
  ~PairV0V0Builder() = default;

  template <modes::Mode mode,
            typename T1,
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
            typename T12,
            typename T13,
            typename T14,
            typename T15>
  void init(o2::framework::HistogramRegistry* registry,
            T1 const& confV0Selection1,
            T2 const& confV0Selection2,
            T3 const& confCprPos,
            T4 const& confCprNeg,
            T5 const& confMixing,
            T6 const& confPairBinning,
            T7 const& confPairCuts,
            std::map<T8, std::vector<framework::AxisSpec>> const& colHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>> const& V0HistSpec1,
            std::map<T10, std::vector<framework::AxisSpec>> const& V0HistSpec2,
            std::map<T11, std::vector<framework::AxisSpec>> const& PosDauHistSpec,
            std::map<T12, std::vector<framework::AxisSpec>> const& NegDauHistSpec,
            std::map<T13, std::vector<framework::AxisSpec>> const& pairHistSpec,
            std::map<T14, std::vector<framework::AxisSpec>> const& cprHistSpecPos,
            std::map<T15, std::vector<framework::AxisSpec>> const& cprHistSpecNeg)
  {

    // check if correlate the same tracks or not
    mSameSpecies = confMixing.sameSpecies.value;

    mColHistManager.template init<mode>(registry, colHistSpec);
    mPairHistManagerSe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerMe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);

    if (mSameSpecies) {
      mV0HistManager1.template init<mode>(registry, V0HistSpec1, confV0Selection1, PosDauHistSpec, NegDauHistSpec);

      mPairHistManagerSe.setMass(confV0Selection1.pdgCodeAbs.value, confV0Selection1.pdgCodeAbs.value);
      mPairHistManagerSe.setCharge(1, 1);
      mCprSe.init(registry, cprHistSpecPos, cprHistSpecNeg, confCprPos, confCprPos);

      mPairHistManagerMe.setMass(confV0Selection1.pdgCodeAbs.value, confV0Selection1.pdgCodeAbs.value);
      mPairHistManagerMe.setCharge(1, 1);
      mCprMe.init(registry, cprHistSpecPos, cprHistSpecNeg, confCprPos, confCprNeg);
    } else {
      mV0HistManager1.template init<mode>(registry, V0HistSpec1, confV0Selection1, PosDauHistSpec, NegDauHistSpec);
      mV0HistManager2.template init<mode>(registry, V0HistSpec2, confV0Selection2, PosDauHistSpec, NegDauHistSpec);

      mPairHistManagerSe.setMass(confV0Selection1.pdgCodeAbs.value, confV0Selection2.pdgCodeAbs.value);
      mPairHistManagerSe.setCharge(1, 1);
      mCprSe.init(registry, cprHistSpecPos, cprHistSpecNeg, confCprPos, confCprNeg);

      mPairHistManagerMe.setMass(confV0Selection1.pdgCodeAbs.value, confV0Selection2.pdgCodeAbs.value);
      mPairHistManagerMe.setCharge(1, 1);
      mCprMe.init(registry, cprHistSpecPos, cprHistSpecNeg, confCprPos, confCprNeg);
    }

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPolicy>(confMixing.policy.value);
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

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& /*lambdaTable*/, T4& partition1, T5& partition2, T6& cache)
  {
    if (mSameSpecies) {
      auto v0Slice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (v0Slice1.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent<mode>(v0Slice1, trackTable, col, mV0HistManager1, mPairHistManagerSe, mCprSe, mPc, mRng, mMixIdenticalParticles);
    } else {
      auto v0Slice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto v0Slice2 = partition2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (v0Slice1.size() == 0 || v0Slice2.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCprSe.setMagField(col.magField());
      pairprocesshelpers::processSameEvent<mode>(v0Slice1, v0Slice2, trackTable, col, mV0HistManager1, mV0HistManager2, mPairHistManagerSe, mCprSe, mPc);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& partition1, T4& partition2, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {

    if (mSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          pairprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }

 private:
  colhistmanager::CollisionHistManager mColHistManager;
  v0histmanager::V0HistManager<prefixV01, prefixPosDau1, prefixNegDau1, v0Type1> mV0HistManager1;
  v0histmanager::V0HistManager<prefixV02, prefixPosDau2, prefixNegDau2, v0Type2> mV0HistManager2;
  pairhistmanager::PairHistManager<prefixSe, modes::Particle::kV0, modes::Particle::kV0> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, modes::Particle::kV0, modes::Particle::kV0> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionV0V0<prefixCprPosSe, prefixCprNegSe> mCprSe;
  closepairrejection::ClosePairRejectionV0V0<prefixCprPosMe, prefixCprNegMe> mCprMe;
  paircleaner::V0V0PairCleaner mPc;
  std::mt19937 mRng;
  pairhistmanager::MixingPolicy mMixingPolicy = pairhistmanager::MixingPolicy::kVtxMult;
  bool mSameSpecies = false;
  int mMixingDepth = 5;
  bool mMixIdenticalParticles = false;
};

template <const char* prefixTrack,
          const char* prefixV0,
          const char* prefixPosDau,
          const char* prefixNegDau,
          const char* prefixSe,
          const char* prefixMe,
          const char* prefixCprSe,
          const char* prefixCprMe,
          modes::V0 v0Type>
class PairTrackV0Builder
{
 public:
  PairTrackV0Builder() = default;
  ~PairTrackV0Builder() = default;

  template <modes::Mode mode,
            typename T1,
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
            typename T12,
            typename T13>
  void init(o2::framework::HistogramRegistry* registry,
            T1 const& confTrackSelection,
            T2 const& confV0Selection,
            T3 const& confCpr,
            T4 const& confMixing,
            T5 const& confPairBinning,
            T6 const& confPairCuts,
            std::map<T7, std::vector<framework::AxisSpec>>& colHistSpec,
            std::map<T8, std::vector<framework::AxisSpec>>& trackHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>>& v0HistSpec,
            std::map<T10, std::vector<framework::AxisSpec>>& posDauHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>>& negDauHistSpec,
            std::map<T12, std::vector<framework::AxisSpec>>& pairHistSpec,
            std::map<T13, std::vector<framework::AxisSpec>>& cprHistSpec)
  {
    mColHistManager.template init<mode>(registry, colHistSpec);

    mTrackHistManager.template init<mode>(registry, trackHistSpec, confTrackSelection);
    mV0HistManager.template init<mode>(registry, v0HistSpec, confV0Selection, posDauHistSpec, negDauHistSpec);

    mPairHistManagerSe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCodeAbs.value, confV0Selection.pdgCodeAbs.value);
    mPairHistManagerSe.setCharge(confTrackSelection.chargeAbs.value, 1);
    mCprSe.init(registry, cprHistSpec, confCpr, confTrackSelection.chargeAbs.value);

    mPairHistManagerMe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCodeAbs.value, confV0Selection.pdgCodeAbs.value);
    mPairHistManagerMe.setCharge(confTrackSelection.chargeAbs.value, 1);
    mCprMe.init(registry, cprHistSpec, confCpr, confTrackSelection.chargeAbs.value);
    mPc.template init<mode>(confPairCuts);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPolicy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*v0table*/, T5& v0Partition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto v0Slice = v0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || v0Slice.size() == 0) {
      return;
    }
    mColHistManager.template fill<mode>(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent<mode>(trackSlice, v0Slice, trackTable, col, mTrackHistManager, mV0HistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void processSameEvent(T1 const& col, T2 const& mcCols, T3 const& trackTable, T4& trackPartition, T5 const& /*v0table*/, T6& v0Partition, T7 const& mcParticles, T8 const& mcMothers, T9 const& mcPartonicMothers, T10& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto v0Slice = v0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || v0Slice.size() == 0) {
      return;
    }
    mColHistManager.template fill<mode>(col, mcCols);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent<mode>(trackSlice, v0Slice, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mTrackHistManager, mV0HistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& v0Partition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, v0Partition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void processMixedEvent(T1 const& cols, T2 const& mcCols, T3& trackTable, T4& trackPartition, T5& v0Partition, T6 const& mcParticles, T7& cache, T8& binsVtxMult, T9& binsVtxCent, T10& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, trackPartition, v0Partition, trackTable, mcParticles, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, trackPartition, v0Partition, trackTable, mcParticles, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, trackPartition, v0Partition, trackTable, mcParticles, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack> mTrackHistManager;
  v0histmanager::V0HistManager<prefixV0, prefixPosDau, prefixNegDau, v0Type> mV0HistManager;
  pairhistmanager::PairHistManager<prefixSe, modes::Particle::kTrack, modes::Particle::kV0> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, modes::Particle::kTrack, modes::Particle::kV0> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackV0<prefixCprSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackV0<prefixCprMe> mCprMe;
  paircleaner::TrackV0PairCleaner mPc;
  pairhistmanager::MixingPolicy mMixingPolicy = pairhistmanager::MixingPolicy::kVtxMult;
  int mMixingDepth = 5;
};

template <const char* prefixTrack,
          const char* prefixResonance,
          const char* prefixPosDau,
          const char* prefixNegDau,
          const char* prefixSe,
          const char* prefixMe,
          const char* prefixCprSe,
          const char* prefixCprMe,
          modes::TwoTrackResonance resonanceType>
class PairTrackTwoTrackResonanceBuilder
{
 public:
  PairTrackTwoTrackResonanceBuilder() = default;
  ~PairTrackTwoTrackResonanceBuilder() = default;

  template <modes::Mode mode,
            typename T1,
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
            typename T12,
            typename T13>
  void init(o2::framework::HistogramRegistry* registry,
            T1 const& confTrackSelection,
            T2 const& confResonanceSelection,
            T3 const& confCpr,
            T4 const& confMixing,
            T5 const& confPairBinning,
            T6 const& confPairCuts,
            std::map<T7, std::vector<framework::AxisSpec>> const& colHistSpec,
            std::map<T8, std::vector<framework::AxisSpec>> const& trackHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>> const& resonanceHistSpec,
            std::map<T10, std::vector<framework::AxisSpec>> const& posDauHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>> const& negDauHistSpec,
            std::map<T12, std::vector<framework::AxisSpec>> const& pairHistSpec,
            std::map<T13, std::vector<framework::AxisSpec>> const& cprHistSpec)
  {
    mColHistManager.template init<mode>(registry, colHistSpec);

    mTrackHistManager.template init<mode>(registry, trackHistSpec, confTrackSelection);
    mResonanceHistManager.template init<mode>(registry, resonanceHistSpec, confResonanceSelection, posDauHistSpec, negDauHistSpec);

    mPairHistManagerSe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCodeAbs.value, confResonanceSelection.pdgCodeAbs.value);
    mPairHistManagerSe.setCharge(confTrackSelection.chargeAbs.value, 1);
    mCprSe.init(registry, cprHistSpec, confCpr, confTrackSelection.chargeAbs.value);

    mPairHistManagerMe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCodeAbs.value, confResonanceSelection.pdgCodeAbs.value);
    mPairHistManagerMe.setCharge(confTrackSelection.chargeAbs.value, 1);
    mCprMe.init(registry, cprHistSpec, confCpr, confTrackSelection.chargeAbs.value);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPolicy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*resonanceTable*/, T5& resonancePartition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto v0Slice = resonancePartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || v0Slice.size() == 0) {
      return;
    }
    mColHistManager.template fill<mode>(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent<mode>(trackSlice, v0Slice, trackTable, col, mTrackHistManager, mResonanceHistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& resonancePartition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, resonancePartition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, resonancePartition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, resonancePartition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack> mTrackHistManager;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<prefixResonance, prefixPosDau, prefixNegDau, resonanceType> mResonanceHistManager;
  pairhistmanager::PairHistManager<prefixSe, modes::Particle::kTrack, modes::Particle::kTwoTrackResonance> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, modes::Particle::kTrack, modes::Particle::kTwoTrackResonance> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackV0<prefixCprSe> mCprSe; // cpr for twotrackresonances and v0 work the same way
  closepairrejection::ClosePairRejectionTrackV0<prefixCprMe> mCprMe; // cpr for twotrackresonances and v0 work the same way
  paircleaner::TrackV0PairCleaner mPc;                               // pc for twotrackresonances and v0 work the same way
  pairhistmanager::MixingPolicy mMixingPolicy = pairhistmanager::MixingPolicy::kVtxMult;
  int mMixingDepth = 5;
};

template <const char* prefixTrack,
          const char* prefixKink,
          const char* prefixChaDau,
          const char* prefixSe,
          const char* prefixMe,
          const char* prefixCprSe,
          const char* prefixCprMe,
          modes::Kink kinkType>
class PairTrackKinkBuilder
{
 public:
  PairTrackKinkBuilder() = default;
  ~PairTrackKinkBuilder() = default;

  template <modes::Mode mode,
            typename T1,
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
            T1 const& confTrackSelection,
            T2 const& confKinkSelection,
            T3 const& confCpr,
            T4 const& confMixing,
            T5 const& confPairBinning,
            T6 const& confPairCuts,
            std::map<T7, std::vector<framework::AxisSpec>> const& colHistSpec,
            std::map<T8, std::vector<framework::AxisSpec>> const& trackHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>> const& kinkHistSpec,
            std::map<T10, std::vector<framework::AxisSpec>> const& chaDauHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>> const& pairHistSpec,
            std::map<T12, std::vector<framework::AxisSpec>> const& cprHistSpec)
  {
    mColHistManager.template init<mode>(registry, colHistSpec);

    mTrackHistManager.template init<mode>(registry, trackHistSpec, confTrackSelection);
    mKinkHistManager.template init<mode>(registry, kinkHistSpec, confKinkSelection, chaDauHistSpec);

    mPairHistManagerSe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCodeAbs.value, confKinkSelection.pdgCodeAbs.value);
    mPairHistManagerSe.setCharge(confTrackSelection.chargeAbs.value, 1); // abs charge of kink daughter is always 1
    mCprSe.init(registry, cprHistSpec, confCpr, confTrackSelection.chargeAbs.value);

    mPairHistManagerMe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCodeAbs.value, confKinkSelection.pdgCodeAbs.value);
    mPairHistManagerMe.setCharge(confTrackSelection.chargeAbs.value, 1); // abs charge of kink daughter is always 1
    mCprMe.init(registry, cprHistSpec, confCpr, confTrackSelection.chargeAbs.value);
    mPc.template init<mode>(confPairCuts);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPolicy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*kinktable*/, T5& kinkPartition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto kinkSlice = kinkPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || kinkSlice.size() == 0) {
      return;
    }
    mColHistManager.fill<mode>(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent<mode>(trackSlice, kinkSlice, trackTable, col, mTrackHistManager, mKinkHistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void processSameEvent(T1 const& col, T2 const& mcCols, T3& trackTable, T4& trackPartition, T5& /*kinktable*/, T6& kinkPartition, T7 const& mcParticles, T8 const& mcMothers, T9 const& mcPartonicMothers, T10& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto kinkSlice = kinkPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || kinkSlice.size() == 0) {
      return;
    }
    mColHistManager.fill<mode>(col, mcCols);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent<mode>(trackSlice, kinkSlice, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mTrackHistManager, mKinkHistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& kinkPartition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, kinkPartition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, kinkPartition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, kinkPartition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void processMixedEvent(T1 const& cols, T2 const& mcCols, T3& trackTable, T4& trackPartition, T5& kinkPartition, T6 const& mcParticles, T7& cache, T8& binsVtxMult, T9& binsVtxCent, T10& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, trackPartition, kinkPartition, trackTable, mcParticles, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, trackPartition, kinkPartition, trackTable, mcParticles, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, trackPartition, kinkPartition, trackTable, mcParticles, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack> mTrackHistManager;
  kinkhistmanager::KinkHistManager<prefixKink, prefixChaDau, kinkType> mKinkHistManager;
  pairhistmanager::PairHistManager<prefixSe, modes::Particle::kTrack, modes::Particle::kKink> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, modes::Particle::kTrack, modes::Particle::kKink> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackKink<prefixCprSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackKink<prefixCprMe> mCprMe;
  paircleaner::TrackKinkPairCleaner mPc;
  pairhistmanager::MixingPolicy mMixingPolicy = pairhistmanager::MixingPolicy::kVtxMult;
  int mMixingDepth = 5;
};

template <const char* prefixTrack,
          const char* prefixCascade,
          const char* prefixBachelor,
          const char* prefixPosDau,
          const char* prefixNegDau,
          const char* prefixSe,
          const char* prefixMe,
          const char* prefixCprBachelorSe,
          const char* prefixCprV0DaughterSe,
          const char* prefixCprBachelorMe,
          const char* prefixCprV0DaughterMe,
          modes::Cascade cascadeType>
class PairTrackCascadeBuilder
{
 public:
  PairTrackCascadeBuilder() = default;
  ~PairTrackCascadeBuilder() = default;

  template <modes::Mode mode,
            typename T1,
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
            typename T12,
            typename T13,
            typename T14,
            typename T15,
            typename T16>
  void init(o2::framework::HistogramRegistry* registry,
            T1 const& confTrackSelection,
            T2 const& confCascadeSelection,
            T3 const& confCprBachelor,
            T4 const& confCprV0Daughter,
            T5 const& confMixing,
            T6 const& confPairBinning,
            T7 const& confPairCuts,
            std::map<T8, std::vector<framework::AxisSpec>> const& colHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>> const& trackHistSpec,
            std::map<T10, std::vector<framework::AxisSpec>> const& cascadeHistSpec,
            std::map<T11, std::vector<framework::AxisSpec>> const& bachelorHistSpec,
            std::map<T12, std::vector<framework::AxisSpec>> const& posDauHistSpec,
            std::map<T13, std::vector<framework::AxisSpec>> const& negDauHistSpec,
            std::map<T14, std::vector<framework::AxisSpec>> const& pairHistSpec,
            std::map<T15, std::vector<framework::AxisSpec>> const& cprHistSpecBachelor,
            std::map<T16, std::vector<framework::AxisSpec>> const& cprHistSpecV0Daughter)
  {
    mColHistManager.template init<mode>(registry, colHistSpec);

    mTrackHistManager.template init<mode>(registry, trackHistSpec, confTrackSelection);
    mCascadeHistManager.template init<mode>(registry, cascadeHistSpec, confCascadeSelection, bachelorHistSpec, posDauHistSpec, negDauHistSpec);

    mPairHistManagerSe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerSe.setMass(confTrackSelection.pdgCodeAbs.value, confCascadeSelection.pdgCodeAbs.value);
    mPairHistManagerSe.setCharge(confTrackSelection.chargeAbs.value, 1);
    mCprSe.init(registry, cprHistSpecBachelor, cprHistSpecV0Daughter, confCprBachelor, confCprV0Daughter, confTrackSelection.chargeAbs.value);

    mPairHistManagerMe.template init<mode>(registry, pairHistSpec, confPairBinning, confPairCuts);
    mPairHistManagerMe.setMass(confTrackSelection.pdgCodeAbs.value, confCascadeSelection.pdgCodeAbs.value);
    mPairHistManagerMe.setCharge(confTrackSelection.chargeAbs.value, 1);
    mCprMe.init(registry, cprHistSpecBachelor, cprHistSpecV0Daughter, confCprBachelor, confCprV0Daughter, confTrackSelection.chargeAbs.value);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPolicy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& trackPartition, T4& /*cascadeTable*/, T5& v0Partition, T6& cache)
  {
    auto trackSlice = trackPartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto v0Slice = v0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (trackSlice.size() == 0 || v0Slice.size() == 0) {
      return;
    }
    mColHistManager.template fill<mode>(col);
    mCprSe.setMagField(col.magField());
    pairprocesshelpers::processSameEvent<mode>(trackSlice, v0Slice, trackTable, col, mTrackHistManager, mCascadeHistManager, mPairHistManagerSe, mCprSe, mPc);
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& trackPartition, T4& v0Partition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, v0Partition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, trackPartition, v0Partition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPc);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager mColHistManager;
  trackhistmanager::TrackHistManager<prefixTrack> mTrackHistManager;
  cascadehistmanager::CascadeHistManager<prefixCascade, prefixBachelor, prefixPosDau, prefixNegDau, cascadeType> mCascadeHistManager;
  pairhistmanager::PairHistManager<prefixSe, modes::Particle::kTrack, modes::Particle::kCascade> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, modes::Particle::kTrack, modes::Particle::kCascade> mPairHistManagerMe;
  closepairrejection::ClosePairRejectionTrackCascade<prefixCprBachelorSe, prefixCprV0DaughterSe> mCprSe;
  closepairrejection::ClosePairRejectionTrackCascade<prefixCprBachelorMe, prefixCprV0DaughterMe> mCprMe;
  paircleaner::TrackCascadePairCleaner mPc;
  pairhistmanager::MixingPolicy mMixingPolicy = pairhistmanager::MixingPolicy::kVtxMult;
  int mMixingDepth = 5;
};

} // namespace pairbuilder
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_PAIRBUILDER_H_
