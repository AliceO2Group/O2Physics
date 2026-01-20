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

/// \file tripletBuilder.h
/// \brief histogram manager for pair tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_TRIPLETBUILDER_H_
#define PWGCF_FEMTO_CORE_TRIPLETBUILDER_H_

#include "PWGCF/Femto/Core/closeTripletRejection.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/tripletCleaner.h"
#include "PWGCF/Femto/Core/tripletHistManager.h"
#include "PWGCF/Femto/Core/tripletProcessHelpers.h"
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
namespace tripletbuilder
{

template <const char* prefixTrack1,
          const char* prefixTrack2,
          const char* prefixTrack3,
          const char* prefixSe,
          const char* prefixMe,
          const char* prefixCtrSeTrack1Track2,
          const char* prefixCtrSeTrack2Track3,
          const char* prefixCtrSeTrack1Track3,
          const char* prefixCtrMeTrack1Track2,
          const char* prefixCtrMeTrack2Track3,
          const char* prefixCtrMeTrack1Track3>
class TripletTrackTrackTrackBuilder
{
 public:
  TripletTrackTrackTrackBuilder() = default;
  ~TripletTrackTrackTrackBuilder() = default;

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
            T1 const& confTrackSelection1,
            T2 const& confTrackSelection2,
            T3 const& confTrackSelection3,
            T4 const& confCtr,
            T5 const& confMixing,
            T6 const& confTripletBinning,
            T7 const& confTripletCuts,
            std::map<T8, std::vector<framework::AxisSpec>> const& colHistSpec,
            std::map<T9, std::vector<framework::AxisSpec>> const& trackHistSpec1,
            std::map<T10, std::vector<framework::AxisSpec>> const& trackHistSpec2,
            std::map<T11, std::vector<framework::AxisSpec>> const& trackHistSpec3,
            std::map<T12, std::vector<framework::AxisSpec>> const& pairHistSpec,
            std::map<T13, std::vector<framework::AxisSpec>> const& cprHistSpec)
  {
    // check if correlate the same tracks or not
    mTrack1Track2Track3AreSameSpecies = confMixing.particle123AreSameSpecies.value;
    mTrack1Track2AreSameSpecies = confMixing.particle12AreSameSpecies.value;

    if (mTrack1Track2Track3AreSameSpecies && mTrack1Track2AreSameSpecies) {
      LOG(fatal) << "Option Track 1&2 are identical and Option Track 1&2&3 are identical is activated. Breaking...";
    }

    mColHistManager.template init<mode>(registry, colHistSpec);
    mTripletHistManagerSe.template init<mode>(registry, pairHistSpec, confTripletBinning, confTripletCuts);
    mTripletHistManagerMe.template init<mode>(registry, pairHistSpec, confTripletBinning, confTripletCuts);

    mTc.template init<mode>(confTripletCuts);

    if (mTrack1Track2Track3AreSameSpecies) {
      // Track1 & Track2 & Track3 are the same particle species
      mTrackHistManager1.template init<mode>(registry, trackHistSpec1, confTrackSelection1);

      mTripletHistManagerSe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value);
      mTripletHistManagerSe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);
      mCtrSe.init(registry, cprHistSpec, confCtr, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);

      mTripletHistManagerMe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value);
      mTripletHistManagerMe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);
      mCtrMe.init(registry, cprHistSpec, confCtr, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value);
    } else if (mTrack1Track2AreSameSpecies) {
      // Track1 & Track2 & are the same particle species and track 3 is something else
      mTrackHistManager1.template init<mode>(registry, trackHistSpec1, confTrackSelection1);
      mTrackHistManager3.template init<mode>(registry, trackHistSpec3, confTrackSelection2);

      mTripletHistManagerSe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value, confTrackSelection3.pdgCodeAbs.value);
      mTripletHistManagerSe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection3.chargeAbs.value);
      mCtrSe.init(registry, cprHistSpec, confCtr, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection3.chargeAbs.value);

      mTripletHistManagerMe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection1.pdgCodeAbs.value, confTrackSelection3.pdgCodeAbs.value);
      mTripletHistManagerMe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection3.chargeAbs.value);
      mCtrMe.init(registry, cprHistSpec, confCtr, confTrackSelection1.chargeAbs.value, confTrackSelection1.chargeAbs.value, confTrackSelection3.chargeAbs.value);
    } else {
      // all three tracks are different
      mTrackHistManager1.template init<mode>(registry, trackHistSpec1, confTrackSelection1);
      mTrackHistManager2.template init<mode>(registry, trackHistSpec2, confTrackSelection2);
      mTrackHistManager3.template init<mode>(registry, trackHistSpec3, confTrackSelection2);

      mTripletHistManagerSe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection2.pdgCodeAbs.value, confTrackSelection3.pdgCodeAbs.value);
      mTripletHistManagerSe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value, confTrackSelection3.chargeAbs.value);
      mCtrSe.init(registry, cprHistSpec, confCtr, confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value, confTrackSelection3.chargeAbs.value);

      mTripletHistManagerMe.setMass(confTrackSelection1.pdgCodeAbs.value, confTrackSelection2.pdgCodeAbs.value, confTrackSelection3.pdgCodeAbs.value);
      mTripletHistManagerMe.setCharge(confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value, confTrackSelection3.chargeAbs.value);
      mCtrMe.init(registry, cprHistSpec, confCtr, confTrackSelection1.chargeAbs.value, confTrackSelection2.chargeAbs.value, confTrackSelection3.chargeAbs.value);
    }

    // setup mixing
    mMixingPolicy = static_cast<triplethistmanager::MixingPolicy>(confMixing.policy.value);
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
      if (mTrack1Track2Track3AreSameSpecies) {
        mDist = std::uniform_int_distribution<>(tripletprocesshelpers::kOrder123, tripletprocesshelpers::kOrder321);
      }
      if (mTrack1Track2AreSameSpecies) {
        mDist = std::uniform_int_distribution<>(tripletprocesshelpers::kOrder123, tripletprocesshelpers::kOrder213);
      }
    }
  }

  // data
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processSameEvent(T1 const& col, T2& trackTable, T3& partition1, T4& partition2, T5& partition3, T6& cache)
  {
    tripletprocesshelpers::TripletOrder tripletOrder = tripletprocesshelpers::kOrder123;
    if (mTrack1Track2Track3AreSameSpecies) {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCtrSe.setMagField(col.magField());
      if (mMixIdenticalParticles) {
        tripletOrder = static_cast<tripletprocesshelpers::TripletOrder>(mDist(mRng));
      }
      tripletprocesshelpers::processSameEvent<mode>(trackSlice1, trackTable, col, mTrackHistManager1, mTripletHistManagerSe, mCtrSe, mTc, tripletOrder);
    } else if (mTrack1Track2AreSameSpecies) {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice3 = partition3->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice3.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCtrSe.setMagField(col.magField());
      if (mMixIdenticalParticles) {
        tripletOrder = static_cast<tripletprocesshelpers::TripletOrder>(mDist(mRng));
      }
      tripletprocesshelpers::processSameEvent<mode>(trackSlice1, trackSlice3, trackTable, col, mTrackHistManager1, mTrackHistManager3, mTripletHistManagerSe, mCtrSe, mTc, tripletOrder);
    } else {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice2 = partition2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice3 = partition3->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice2.size() == 0 || trackSlice3.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCtrSe.setMagField(col.magField());
      tripletprocesshelpers::processSameEvent<mode>(trackSlice1, trackSlice2, trackSlice3, trackTable, col, mTrackHistManager1, mTrackHistManager2, mTrackHistManager3, mTripletHistManagerSe, mCtrSe, mTc);
    }
  }

  // mc
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void processSameEvent(T1 const& col, T2 const& mcCols, T3& trackTable, T4& partition1, T5& partition2, T6& partition3, T7 const& mcParticles, T8 const& mcMothers, T9 const& mcPartonicMothers, T10& cache)
  {
    tripletprocesshelpers::TripletOrder tripletOrder = tripletprocesshelpers::kOrder123;
    if (mTrack1Track2Track3AreSameSpecies) {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCtrSe.setMagField(col.magField());
      if (mMixIdenticalParticles) {
        tripletOrder = static_cast<tripletprocesshelpers::TripletOrder>(mDist(mRng));
      }
      tripletprocesshelpers::processSameEvent<mode>(trackSlice1, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mTrackHistManager1, mTripletHistManagerSe, mCtrSe, mTc, tripletOrder);
    } else if (mTrack1Track2AreSameSpecies) {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice3 = partition3->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice3.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCtrSe.setMagField(col.magField());
      if (mMixIdenticalParticles) {
        tripletOrder = static_cast<tripletprocesshelpers::TripletOrder>(mDist(mRng));
      }
      tripletprocesshelpers::processSameEvent<mode>(trackSlice1, trackSlice3, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mTrackHistManager1, mTrackHistManager3, mTripletHistManagerSe, mCtrSe, mTc, tripletOrder);
    } else {
      auto trackSlice1 = partition1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice2 = partition2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      auto trackSlice3 = partition3->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      if (trackSlice1.size() == 0 || trackSlice2.size() == 0 || trackSlice3.size() == 0) {
        return;
      }
      mColHistManager.template fill<mode>(col);
      mCtrSe.setMagField(col.magField());
      tripletprocesshelpers::processSameEvent<mode>(trackSlice1, trackSlice2, trackSlice3, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mTrackHistManager1, mTrackHistManager2, mTrackHistManager3, mTripletHistManagerSe, mCtrSe, mTc);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void processMixedEvent(T1 const& cols, T2& trackTable, T3& partition1, T4& partition2, T5& partition3, T6& cache, T7& binsVtxMult, T8& binsVtxCent, T9& binsVtxMultCent)
  {
    if (mTrack1Track2Track3AreSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition1, partition1, trackTable, cache, binsVtxMult, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition1, partition1, trackTable, cache, binsVtxCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition1, partition1, trackTable, cache, binsVtxMultCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else if (mTrack1Track2AreSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition1, partition3, trackTable, cache, binsVtxMult, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition1, partition3, trackTable, cache, binsVtxCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition1, partition3, trackTable, cache, binsVtxMultCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, partition3, trackTable, cache, binsVtxMult, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, partition3, trackTable, cache, binsVtxCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, partition1, partition2, partition3, trackTable, cache, binsVtxMultCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
  void processMixedEvent(T1 const& cols, T2 const& mcCols, T3& trackTable, T4& partition1, T5& partition2, T6& partition3, T7 const& mcParticles, T8& cache, T9& binsVtxMult, T10& binsVtxCent, T11& binsVtxMultCent)
  {
    if (mTrack1Track2Track3AreSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition1, partition1, trackTable, mcParticles, cache, binsVtxMult, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition1, partition1, trackTable, mcParticles, cache, binsVtxCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition1, partition1, trackTable, mcParticles, cache, binsVtxMultCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else if (mTrack1Track2AreSameSpecies) {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition1, partition3, trackTable, mcParticles, cache, binsVtxMult, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition1, partition3, trackTable, mcParticles, cache, binsVtxCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition1, partition3, trackTable, mcParticles, cache, binsVtxMultCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        default:
          LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
      }
    } else {
      switch (mMixingPolicy) {
        case static_cast<int>(pairhistmanager::kVtxMult):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition2, partition3, trackTable, mcParticles, cache, binsVtxMult, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition2, partition3, trackTable, mcParticles, cache, binsVtxCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
          break;
        case static_cast<int>(pairhistmanager::kVtxMultCent):
          tripletprocesshelpers::processMixedEvent<mode>(cols, mcCols, partition1, partition2, partition3, trackTable, mcParticles, cache, binsVtxMultCent, mMixingDepth, mTripletHistManagerMe, mCtrMe, mTc);
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
  trackhistmanager::TrackHistManager<prefixTrack3> mTrackHistManager3;
  triplethistmanager::TripletHistManager<prefixSe, modes::Particle::kTrack, modes::Particle::kTrack, modes::Particle::kTrack> mTripletHistManagerSe;
  triplethistmanager::TripletHistManager<prefixMe, modes::Particle::kTrack, modes::Particle::kTrack, modes::Particle::kTrack> mTripletHistManagerMe;

  closetripletrejection::CloseTripletRejectionTrackTrackTrack<prefixCtrSeTrack1Track2, prefixCtrSeTrack2Track3, prefixCtrSeTrack1Track3> mCtrSe;
  closetripletrejection::CloseTripletRejectionTrackTrackTrack<prefixCtrMeTrack1Track2, prefixCtrMeTrack2Track3, prefixCtrMeTrack1Track3> mCtrMe;
  tripletcleaner::TrackTrackTrackTripletCleaner mTc;
  triplethistmanager::MixingPolicy mMixingPolicy = triplethistmanager::MixingPolicy::kVtxMult;
  bool mTrack1Track2Track3AreSameSpecies = false;
  bool mTrack1Track2AreSameSpecies = false;
  int mMixingDepth = 5;
  bool mMixIdenticalParticles = false;
  std::mt19937 mRng;
  std::uniform_int_distribution<> mDist;
};

} // namespace tripletbuilder
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_TRIPLETBUILDER_H_
