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

/// \file collisionBuilder.h
/// \brief collision builder
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_COLLISIONBUILDER_H_
#define PWGCF_FEMTO_CORE_COLLISIONBUILDER_H_

#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/CCDB/EventSelectionParams.h"

#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "fairlogger/Logger.h"

#include <cmath>
#include <string>

namespace o2::analysis::femto
{
namespace collisionbuilder
{

// configurables for collision selection
struct ConfCollisionFilter : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionFilter");
  o2::framework::Configurable<float> vtxZMin{"vtxZMin", -10.f, "Minimum vertex Z position (cm)"};
  o2::framework::Configurable<float> vtxZMax{"vtxZMax", 10.f, "Maximum vertex Z position (cm)"};
  o2::framework::Configurable<float> multMin{"multMin", 0.f, "Minimum multiplicity"};
  o2::framework::Configurable<float> multMax{"multMax", 999.f, "Maximum multiplicity"};
  o2::framework::Configurable<float> centMin{"centMin", 0.f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> centMax{"centMax", 999.f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> spherMin{"spherMin", 0.f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> spherMax{"spherMax", 2.f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> occupancyMin{"occupancyMin", 0.f, "Minimum occupancy (Cut works at producer level and if occupancy is stored also at consumer level)"};
  o2::framework::Configurable<float> occupancyMax{"occupancyMax", 1e6f, "Maximum occupancy (Cut works at producer level and if occupancy is stored also at consumer level)"};
  o2::framework::Configurable<float> magFieldMin{"magFieldMin", -1.f, "Minimum magnetic field strength (T)"};
  o2::framework::Configurable<float> magFieldMax{"magFieldMax", 1.f, "Maximum magnetic field strength (T)"};
};

struct ConfCollisionFlags : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<bool> sel8{"sel8", true, "Use sel8 (Should suffice for pp, for PbPb the other flags are more important)"};
  o2::framework::Configurable<bool> noSameBunchPileup{"noSameBunchPileup", false, "Reject collisions in case of pileup with another collision in the same foundBC"};
  o2::framework::Configurable<bool> isVertexItsTpc{"isVertexItsTpc", false, "At least one ITS-TPC track found for the vertex"};
  o2::framework::Configurable<bool> isGoodZvtxFt0VsPv{"isGoodZvtxFt0VsPv", false, "small difference between z-vertex from PV and from FT0"};
  o2::framework::Configurable<bool> noCollInTimeRangeNarrow{"noCollInTimeRangeNarrow", false, "no other collisions in specified time range (narrower than Strict)"};
  o2::framework::Configurable<bool> noCollInTimeRangeStrict{"noCollInTimeRangeStrict", false, "no other collisions in specified time range"};
  o2::framework::Configurable<bool> noCollInTimeRangeStandard{"noCollInTimeRangeStandard", false, "no other collisions in specified time range with per-collision multiplicity above threshold"};
  o2::framework::Configurable<bool> noCollInRofStrict{"noCollInRofStrict", false, "no other collisions in this Readout Frame"};
  o2::framework::Configurable<bool> noCollInRofStandard{"noCollInRofStandard", false, "no other collisions in this Readout Frame with per-collision multiplicity above threshold"};
  o2::framework::Configurable<bool> noHighMultCollInPrevRof{"noHighMultCollInPrevRof", false, "veto an event if FT0C amplitude in previous ITS ROF is above threshold"};
  o2::framework::Configurable<bool> isGoodItsLayer3{"isGoodItsLayer3", false, "number of inactive chips on ITS layer 3 is below maximum allowed value"};
  o2::framework::Configurable<bool> isGoodItsLayer0123{"isGoodItsLayer0123", false, "numbers of inactive chips on ITS layers 0-3 are below maximum allowed values"};
  o2::framework::Configurable<bool> isGoodItsLayersAll{"isGoodItsLayersAll", false, "numbers of inactive chips on all ITS layers are below maximum allowed values"};
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class CollisionSelection
{
 public:
  CollisionSelection() {}
  virtual ~CollisionSelection() = default;

  template <typename T1, typename T2>
  void configure(T1 const& filter, T2 const& flags)
  {
    // flags
    mSel8 = flags.sel8.value;
    mNoSameBunchPileup = flags.noSameBunchPileup.value;
    mIsVertexItsTpc = flags.isVertexItsTpc.value;
    mIsGoodZvtxFt0VsPv = flags.isGoodZvtxFt0VsPv.value;
    mNoCollInTimeRangeNarrow = flags.noCollInTimeRangeNarrow.value;
    mNoCollInTimeRangeStrict = flags.noCollInTimeRangeStrict.value;
    mNoCollInTimeRangeStandard = flags.noCollInTimeRangeStandard.value;
    mNoCollInRofStrict = flags.noCollInRofStrict.value;
    mNoCollInRofStandard = flags.noCollInRofStandard.value;
    mNoHighMultCollInPrevRof = flags.noHighMultCollInPrevRof.value;
    mIsGoodItsLayer3 = flags.isGoodItsLayer3.value;
    mIsGoodItsLayer0123 = flags.isGoodItsLayer0123.value;
    mIsGoodItsLayersAll = flags.isGoodItsLayersAll.value;

    // cuts
    mVtxZMin = filter.vtxZMin.value;
    mVtxZMax = filter.vtxZMax.value;
    mMagFieldMin = filter.magFieldMin.value;
    mMagFieldMax = filter.magFieldMax.value;
    mMultMin = filter.multMin.value;
    mMultMax = filter.multMax.value;
    mCentMin = filter.centMin.value;
    mCentMax = filter.centMax.value;
    mSphericityMin = filter.spherMin.value;
    mSphericityMax = filter.spherMax.value;
    mOccupancyMin = filter.occupancyMin.value;
    mOccupancyMax = filter.occupancyMax.value;
  };

  void setMagneticField(float MagField)
  {
    mMagField = MagField;
  }

  float getMagneticField()
  {
    return mMagField;
  }

  template <typename T>
  void setSphericity(T tracks)
  {
    mSphericity = utils::sphericity(tracks);
  }

  float getSphericity() const { return mSphericity; }

  template <modes::System system, typename T>
  void setCentrality(const T& col)
  {
    if constexpr (modes::isFlagSet(system, modes::System::kPP)) {
      mCentrality = col.centFT0M();
    }
    if constexpr (modes::isFlagSet(system, modes::System::kPbPb)) {
      mCentrality = col.centFT0C();
    }
  }
  float getCentrality() const { return mCentrality; }

  template <typename T>
  bool checkCuts(T const& col)
  {

    // flags
    if (mSel8 && !col.sel8()) { // might change in the future, but sel8 checks TVXTrigger and cuts out time frame border and ITS readout frame border
      return false;
    }
    if (mNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (mIsVertexItsTpc && !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (mIsGoodZvtxFt0VsPv && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (mNoCollInTimeRangeNarrow && !col.selection_bit(aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (mNoCollInTimeRangeStrict && !col.selection_bit(aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (mNoCollInTimeRangeStandard && !col.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (mNoCollInRofStrict && !col.selection_bit(aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (mNoCollInRofStandard && !col.selection_bit(aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (mNoHighMultCollInPrevRof && !col.selection_bit(aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if (mIsGoodItsLayer3 && !col.selection_bit(aod::evsel::kIsGoodITSLayer3)) {
      return false;
    }
    if (mIsGoodItsLayer0123 && !col.selection_bit(aod::evsel::kIsGoodITSLayer0123)) {
      return false;
    }
    if (mIsGoodItsLayersAll && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }

    // cuts
    if (col.posZ() < mVtxZMin || col.posZ() > mVtxZMax) {
      return false;
    }
    if (col.multNTracksPV() < mMultMin || col.multNTracksPV() > mMultMax) {
      return false;
    }
    if (mCentrality < mCentMin || mCentrality > mCentMax) {
      return false;
    }
    if (mMagField < mMagFieldMin || mMagField > mMagFieldMax) {
      return false;
    }
    if (mSphericity < mSphericityMin || mSphericity > mSphericityMax) {
      return false;
    }
    return true;
  }

 private:
  float mCentrality = 0.f;

  // flags
  bool mSel8 = false;
  bool mNoSameBunchPileup = false;
  bool mIsVertexItsTpc = false;
  bool mIsGoodZvtxFt0VsPv = false;
  bool mNoCollInTimeRangeNarrow = false;
  bool mNoCollInTimeRangeStrict = false;
  bool mNoCollInTimeRangeStandard = false;
  bool mNoCollInRofStrict = false;
  bool mNoCollInRofStandard = false;
  bool mNoHighMultCollInPrevRof = false;
  bool mIsGoodItsLayer3 = false;
  bool mIsGoodItsLayer0123 = false;
  bool mIsGoodItsLayersAll = false;

  // cuts
  float mVtxZMin = -12.f;
  float mVtxZMax = -12.f;
  float mSphericityMin = 0.f;
  float mSphericityMax = 2.f;
  float mMagFieldMin = -1.f;
  float mMagFieldMax = 1.f;
  float mMultMin = 0.f;
  float mMultMax = 999.f;
  float mCentMin = 0.f;
  float mCentMax = 999.f;
  float mOccupancyMin = 0.;
  float mOccupancyMax = 1e6f;

  float mMagField = 0.f;
  float mSphericity = 0.f;
};

struct CollisionBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FCols> producedCollision;
  o2::framework::Produces<o2::aod::FColOccs> producedOccupancy;
  o2::framework::Produces<o2::aod::FColQns> producedQns;
  o2::framework::Produces<o2::aod::FColPos> producedPositions;
  o2::framework::Produces<o2::aod::FColMults> producedMultiplicityEstimators;
  o2::framework::Produces<o2::aod::FColCents> producedCentralityEstimators;
};

struct ConfCollisionTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionTables");
  o2::framework::Configurable<int> produceCollisions{"produceCollisions", -1, "Produce Collisions (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceOccupancy{"produceOccupancy", -1, "Produce Occupancy (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceQns{"produceQns", -1, "Produce Qn (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producePositions{"producePositions", -1, "Produce Positions (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMults{"produceMults", -1, "Produce Multiplicities (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceCents{"produceCents", -1, "Produce Centralities (-1: auto; 0 off; 1 on)"};
};

class CollisionBuilder
{
 public:
  CollisionBuilder() {}
  virtual ~CollisionBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(T1& filter, T2& flags, T3& table, T4& initContext)
  {
    collisionSelection.configure(filter, flags);
    LOG(info) << "Initialize femto collision builder...";
    producedCollisions = utils::enableTable("FCols_001", table.produceCollisions.value, initContext);
    produceOccupancy = utils::enableTable("FColOccs_001", table.produceOccupancy.value, initContext);
    produceQns = utils::enableTable("FColQnBins_001", table.produceQns.value, initContext);
    producedPositions = utils::enableTable("FColPos_001", table.producePositions.value, initContext);
    producedMultiplicities = utils::enableTable("FColMults_001", table.produceMults.value, initContext);
    producedCentralities = utils::enableTable("FColCents_001", table.produceCents.value, initContext);
    if (producedCollisions || producedPositions || producedMultiplicities || producedCentralities) {
      fillAnyTable = true;
    } else {
      LOG(info) << "No tables configured";
    }
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2>
  void buildCollision(T1& col, T2 tracks, float magField)
  {
    collisionSelection.setMagneticField(magField);
    collisionSelection.setSphericity(tracks);
    collisionSelection.setCentrality<system>(col);
  }

  template <typename T>
  bool checkCuts(T const& col)
  {
    return collisionSelection.checkCuts(col);
  }

  template <modes::System system, typename T1, typename T2>
  void fillCollision(T1& collisionProducts, T2 const& col)
  {
    if (!fillAnyTable) {
      return;
    }
    if (producedCollisions) {
      collisionProducts.producedCollision(col.posZ(),
                                          col.multNTracksPV(),
                                          collisionSelection.getCentrality(),
                                          collisionSelection.getSphericity(),
                                          collisionSelection.getMagneticField());
    }

    if (produceOccupancy) {
      collisionProducts.producedOccupancy(col.trackOccupancyInTimeRange());
    }

    if (producedPositions) {
      collisionProducts.producedPositions(col.posX(),
                                          col.posY());
    }
    if (producedMultiplicities) {
      collisionProducts.producedMultiplicityEstimators(
        col.multFT0A(),
        col.multFT0C(),
        col.multNTracksPVeta1(),
        col.multNTracksPVetaHalf(),
        col.trackOccupancyInTimeRange(),
        col.ft0cOccupancyInTimeRange());
    }
    if (producedCentralities) {
      collisionProducts.producedCentralityEstimators(
        col.centFT0A(),
        col.centFT0C());
    }

    // PbPb specific columns
    if constexpr (modes::isFlagSet(system, modes::System::kPbPb)) {
      if (produceQns) {
        collisionProducts.producedQns(utils::qn(col));
      }
    }
  }

 private:
  CollisionSelection collisionSelection;
  bool fillAnyTable = false;
  bool producedCollisions = false;
  bool produceOccupancy = false;
  bool produceQns = false;
  bool producedPositions = false;
  bool producedMultiplicities = false;
  bool producedCentralities = false;
};
}; // namespace collisionbuilder
}; // namespace o2::analysis::femto
;
#endif // PWGCF_FEMTO_CORE_COLLISIONBUILDER_H_
