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

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/Zorro.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "fairlogger/Logger.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{
namespace collisionbuilder
{

// configurables for collision selection
struct ConfCollisionFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionFilter");
  o2::framework::Configurable<float> vtxZMin{"vtxZMin", -10.f, "Minimum vertex Z position (cm)"};
  o2::framework::Configurable<float> vtxZMax{"vtxZMax", 10.f, "Maximum vertex Z position (cm)"};
  o2::framework::Configurable<float> multMin{"multMin", 0.f, "Minimum multiplicity"};
  o2::framework::Configurable<float> multMax{"multMax", 5000.f, "Maximum multiplicity"};
  o2::framework::Configurable<float> centMin{"centMin", 0.f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> centMax{"centMax", 100.f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> sphericityMin{"sphericityMin", 0.f, "Minimum sphericity"};
  o2::framework::Configurable<float> sphericityMax{"sphericityMax", 1.f, "Maximum sphericity"};
  o2::framework::Configurable<int> magFieldMin{"magFieldMin", -5, "Minimum magnetic field strength (kG)"};
  o2::framework::Configurable<int> magFieldMax{"magFieldMax", 5, "Maximum magnetic field strength (kG)"};
};

struct ConfCollisionBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionBits");
  o2::framework::Configurable<int> sel8{"sel8", 1, "Use sel8 (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> noSameBunchPileup{"noSameBunchPileup", 0, "Reject collisions in case of pileup with another collision in the same foundBC (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> isVertexItsTpc{"isVertexItsTpc", 0, "At least one ITS-TPC track found for the vertex (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> isGoodZvtxFt0VsPv{"isGoodZvtxFt0VsPv", 0, "small difference between z-vertex from PV and from FT0 (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> noCollInTimeRangeNarrow{"noCollInTimeRangeNarrow", 0, "no other collisions in specified time range (narrower than Strict)(-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> noCollInTimeRangeStrict{"noCollInTimeRangeStrict", 0, "no other collisions in specified time range strict (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> noCollInTimeRangeStandard{"noCollInTimeRangeStandard", 0, "no other collisions in specified time range with per-collision multiplicity above threshold (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> noCollInRofStrict{"noCollInRofStrict", 0, "no other collisions in this Readout Frame strict (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> noCollInRofStandard{"noCollInRofStandard", 0, "no other collisions in this Readout Frame with per-collision multiplicity above threshold (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> noHighMultCollInPrevRof{"noHighMultCollInPrevRof", 0, "veto an event if FT0C amplitude in previous ITS ROF is above threshold (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> isGoodItsLayer3{"isGoodItsLayer3", 0, "number of inactive chips on ITS layer 3 is below maximum allowed value (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> isGoodItsLayer0123{"isGoodItsLayer0123", 0, "numbers of inactive chips on ITS layers 0-3 are below maximum allowed values (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<int> isGoodItsLayersAll{"isGoodItsLayersAll", 0, "numbers of inactive chips on all ITS layers are below maximum allowed values (-1: stored in bitmaks; 0 off; 1 on)"};
  o2::framework::Configurable<std::vector<float>> occupancyMin{"occupancyMin", {}, "Minimum occpancy"};
  o2::framework::Configurable<std::vector<float>> occupancyMax{"occupancyMax", {}, "Maximum occpancy"};
  o2::framework::Configurable<std::vector<float>> sphericityMin{"sphericityMin", {}, "Minimum sphericity"};
  o2::framework::Configurable<std::vector<float>> sphericityMax{"sphericityMax", {}, "Maximum sphericity"};
  o2::framework::Configurable<std::vector<std::string>> triggers{"triggers", {}, "List of all triggers to be used"};
};

struct ConfCcdb : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("ConfCcdb");
  o2::framework::Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL to ccdb"};
  o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/Config/GRPMagField", "Path to GRP object (Run3 -> GLO/Config/GRPMagField/Run2 -> GLO/GRP/GRP"};
  o2::framework::Configurable<std::string> triggerPath{"triggerPath", "EventFiltering/Zorro/", "CCDB path for trigger information"};
};

struct ConfCollisionRctFlags : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionRctFlags");
  o2::framework::Configurable<bool> useRctFlags{"useRctFlags", true, "Set to true to use RCT flags"};
  o2::framework::Configurable<std::string> label{"label", std::string("CBT_hadronPID"), "Which RCT flag to check"};
  o2::framework::Configurable<bool> useZdc{"useZdc", false, "Whether to use ZDC (only use for PbPb)"};
  o2::framework::Configurable<bool> treatLimitedAcceptanceAsBad{"treatLimitedAcceptanceAsBad", false, "Whether to treat limited acceptance as bad or not"};
};

// configurables for collision selection
struct ConfCollisionSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionSelection");
  o2::framework::Configurable<float> vtxZMin{"vtxZMin", -10.f, "Minimum vertex Z position (cm)"};
  o2::framework::Configurable<float> vtxZMax{"vtxZMax", 10.f, "Maximum vertex Z position (cm)"};
  o2::framework::Configurable<float> multMin{"multMin", 0.f, "Minimum multiplicity"};
  o2::framework::Configurable<float> multMax{"multMax", 5000.f, "Maximum multiplicity"};
  o2::framework::Configurable<float> centMin{"centMin", 0.f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> centMax{"centMax", 100.f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<int> magFieldMin{"magFieldMin", -5, "Minimum magnetic field strength (kG)"};
  o2::framework::Configurable<int> magFieldMax{"magFieldMax", 5, "Maximum magnetic field strength (kG)"};
  o2::framework::Configurable<aod::femtodatatypes::CollisionMaskType> collisionMask{"collisionMask", 0x0, "Bitmask for collision"};
};

/// enum for all collision selections
enum CollisionSels {
  // collsion selection flags
  kSel8,                      ///< Sel8
  kNoSameBunchPileUp,         ///< Reject collisions in case of pileup with another collision in the same foundBC
  kIsVertexItsTpc,            ///< At least one ITS-TPC track found for the vertex
  kIsGoodZvtxFt0VsPv,         ///< small difference between z-vertex from PV and from FT0
  kNoCollInTimeRangeNarrow,   ///< no other collisions in specified time range (narrower than Strict)
  kNoCollInTimeRangeStrict,   ///< no other collisions in specified time range strict
  kNoCollInTimeRangeStandard, ///< no other collisions in specified time range
  kNoCollInRofStrict,         ///< no other collisions in this Readout Frame strict
  kNoCollInRofStandard,       ///< no other collisions in this Readout Frame
  kNoHighMultCollInPrevRof,   ///< veto an event if FT0C amplitude in previous ITS ROF is above threshold
  kIsGoodItsLayer3,           ///< number of inactive chips on ITS layer 3 is below maximum allowed value
  kIsGoodItsLayer0123,        ///< numbers of inactive chips on ITS layers 0-3 are below maximum allowed values
  kIsGoodItsLayersAll,        ///< numbers of inactive chips on all ITS layers are below maximum allowed values
  kOccupancyMin,              ///< Min. occupancy
  kOccupancyMax,              ///< Max. occupancy
  kSphericityMin,             ///< Min. sphericity
  kSphericityMax,             ///< Max. sphericity

  kTriggers,

  kCollisionSelsMax
};

const char colSelsName[] = "Collision Selection Object";
const std::unordered_map<CollisionSels, std::string> colSelsToString = {
  {kSel8, "Sel8"},
  {kNoSameBunchPileUp, "No same bunch pileup"},
  {kIsVertexItsTpc, "Is vertex ITS TPC"},
  {kIsGoodZvtxFt0VsPv, "Is good zvtx FT0 vs PV"},
  {kNoCollInTimeRangeNarrow, "No collision in time range narrow"},
  {kNoCollInTimeRangeStrict, "No collision in time range strict"},
  {kNoCollInTimeRangeStandard, "No collission in time range standard"},
  {kNoCollInRofStrict, "No collsion in ROF strict"},
  {kNoCollInRofStandard, "No collision in ROF standard"},
  {kNoHighMultCollInPrevRof, "No high mult collsions in previous ROF"},
  {kIsGoodItsLayer3, "Is good ITS layer 3"},
  {kIsGoodItsLayer0123, "Is good ITS layer 0-3"},
  {kIsGoodItsLayersAll, "Is good ITS layer all"},
  {kOccupancyMin, "Minimum Occupancy"},
  {kOccupancyMax, "Maximum Occupancy"},
  {kSphericityMin, "Minimum Sphericity"},
  {kSphericityMax, "Maximum Sphericity"},

  {kTriggers, "Triggers"}

};

class CollisionSelection : public BaseSelection<float, o2::aod::femtodatatypes::CollisionMaskType, kCollisionSelsMax>
{
 public:
  CollisionSelection() = default;
  ~CollisionSelection() = default;

  template <typename T1, typename T2>
  void configure(T1 const& filter, T2 const& config)
  {
    // cuts
    mVtxZMin = filter.vtxZMin.value;
    mVtxZMax = filter.vtxZMax.value;
    mMagFieldMin = filter.magFieldMin.value;
    mMagFieldMax = filter.magFieldMax.value;
    mMultMin = filter.multMin.value;
    mMultMax = filter.multMax.value;
    mCentMin = filter.centMin.value;
    mCentMax = filter.centMax.value;
    mSphericityMin = filter.sphericityMin.value;
    mSphericityMax = filter.sphericityMax.value;

    // flags
    this->addSelection(config.sel8.value, kSel8);
    this->addSelection(config.noSameBunchPileup.value, kNoSameBunchPileUp);
    this->addSelection(config.isGoodZvtxFt0VsPv.value, kIsGoodZvtxFt0VsPv);
    this->addSelection(config.noCollInTimeRangeNarrow.value, kNoCollInTimeRangeNarrow);
    this->addSelection(config.noCollInTimeRangeStrict.value, kNoCollInTimeRangeStrict);
    this->addSelection(config.noCollInTimeRangeStandard.value, kNoCollInTimeRangeStandard);
    this->addSelection(config.noCollInRofStrict.value, kNoCollInRofStrict);
    this->addSelection(config.noCollInRofStandard.value, kNoCollInRofStandard);
    this->addSelection(config.noHighMultCollInPrevRof.value, kNoHighMultCollInPrevRof);
    this->addSelection(config.isGoodItsLayer3.value, kIsGoodItsLayer3);
    this->addSelection(config.isGoodItsLayer0123.value, kIsGoodItsLayer0123);
    this->addSelection(config.isGoodItsLayersAll.value, kIsGoodItsLayersAll);
    this->addSelection(config.occupancyMin.value, kOccupancyMin, limits::kLowerLimit, true, true);
    this->addSelection(config.occupancyMax.value, kOccupancyMax, limits::kUpperLimit, true, true);
    this->addSelection(config.sphericityMin.value, kSphericityMin, limits::kLowerLimit, true, true);
    this->addSelection(config.sphericityMax.value, kSphericityMax, limits::kUpperLimit, true, true);

    std::vector<float> triggerValues(config.triggers.value.size(), 1.f);
    this->addSelection(triggerValues, kTriggers, limits::kEqualArray, false, true);
    this->addComments(kTriggers, config.triggers.value);
  };

  void setMagneticField(int MagField)
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

  template <modes::System system, typename T>
  void setMultiplicity(const T& col)
  {
    if constexpr (modes::isFlagSet(system, modes::System::kPP)) {
      mMultiplicity = col.multNTracksPV();
    }
    if constexpr (modes::isFlagSet(system, modes::System::kPbPb)) {
      // change multiplicity estimator for PbPb?
      mMultiplicity = col.multNTracksPV();
    }
  }
  float getMultiplicity() const { return mMultiplicity; }

  template <typename T>
  bool checkFilters(T const& col) const
  {
    if (col.posZ() < mVtxZMin || col.posZ() > mVtxZMax) {
      return false;
    }
    if (mMultiplicity < mMultMin || mMultiplicity > mMultMax) {
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

  template <typename T>
  void applySelections(T const& col, std::vector<bool> const& triggerDecisions)
  {
    this->reset();

    // casting bool to float gurantees false -> 0 and true -> 1
    // and we check for equality to 1, so evaluation succeeds if the selection bit is true
    this->evaluateObservable(kSel8, static_cast<float>(col.sel8()));
    this->evaluateObservable(kNoSameBunchPileUp, static_cast<float>(col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)));
    this->evaluateObservable(kIsVertexItsTpc, static_cast<float>(col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)));
    this->evaluateObservable(kIsGoodZvtxFt0VsPv, static_cast<float>(col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)));
    this->evaluateObservable(kNoCollInTimeRangeNarrow, static_cast<float>(col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)));
    this->evaluateObservable(kNoCollInTimeRangeStrict, static_cast<float>(col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)));
    this->evaluateObservable(kNoCollInRofStrict, static_cast<float>(col.selection_bit(o2::aod::evsel::kNoCollInRofStrict)));
    this->evaluateObservable(kNoCollInRofStandard, static_cast<float>(col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)));
    this->evaluateObservable(kNoHighMultCollInPrevRof, static_cast<float>(col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)));
    this->evaluateObservable(kIsGoodItsLayer3, static_cast<float>(col.selection_bit(o2::aod::evsel::kIsGoodITSLayer3)));
    this->evaluateObservable(kIsGoodItsLayer0123, static_cast<float>(col.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)));
    this->evaluateObservable(kIsGoodItsLayersAll, static_cast<float>(col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)));

    this->evaluateObservable(kOccupancyMin, col.trackOccupancyInTimeRange());
    this->evaluateObservable(kOccupancyMax, col.trackOccupancyInTimeRange());
    this->evaluateObservable(kSphericityMin, mSphericity);
    this->evaluateObservable(kSphericityMax, mSphericity);

    // for the trigger we need to pass an vector of 0 (false) and 1 (true) for all configured trigger selections
    if (!triggerDecisions.empty()) {
      std::vector<float> trigger(triggerDecisions.size());
      std::transform(triggerDecisions.begin(), triggerDecisions.end(), trigger.begin(), [](bool b) { return b ? 1.0f : 0.0f; });
      this->evaluateObservable(kTriggers, trigger);
    }

    this->assembleBitmask();
  };

 protected:
  // filter cuts
  float mVtxZMin = -12.f;
  float mVtxZMax = -12.f;
  float mSphericityMin = 0.f;
  float mSphericityMax = 1.f;
  float mMagFieldMin = -5.f;
  float mMagFieldMax = 5.f;
  float mMultMin = 0.f;
  float mMultMax = 5000.f;
  float mCentMin = 0.f;
  float mCentMax = 100.f;

  int mMagField = 0.f;
  float mSphericity = 0.f;
  float mCentrality = 0.f;
  float mMultiplicity = 0.f;
};

struct CollisionBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FCols> producedCollision;
  o2::framework::Produces<o2::aod::FColMasks> producedCollisionMask;
  o2::framework::Produces<o2::aod::FColPos> producedPositions;
  o2::framework::Produces<o2::aod::FColSphericities> producedSphericities;
  o2::framework::Produces<o2::aod::FColMults> producedMultiplicityEstimators;
  o2::framework::Produces<o2::aod::FColCents> producedCentralityEstimators;
  o2::framework::Produces<o2::aod::FColQns> producedQns;
};

struct ConfCollisionTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionTables");
  o2::framework::Configurable<int> produceCollisions{"produceCollisions", -1, "Produce Collisions (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceCollisionMasks{"produceCollisionMasks", -1, "Produce Collision Masks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producePositions{"producePositions", -1, "Produce Positions (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceSphericities{"produceSphericities", -1, "Produce Sphericity (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceMults{"produceMults", -1, "Produce Multiplicities (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceCents{"produceCents", -1, "Produce Centralities (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceQns{"produceQns", -1, "Produce Qn (-1: auto; 0 off; 1 on)"};
};

class CollisionBuilder
{
 public:
  CollisionBuilder() = default;
  ~CollisionBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void init(T1& confFilter, T2& confBits, T3& confRct, T4& confCcdb, T5& confTable, T6& initContext)
  {
    mCollisionSelection.configure(confFilter, confBits);
    if (!confBits.triggers.value.empty()) {
      mUseTrigger = true;
      for (size_t i = 0; i < confBits.triggers.value.size(); ++i) {
        mTriggerNames += confBits.triggers.value[i];
        if (i != confBits.triggers.value.size() - 1) {
          mTriggerNames += ",";
        }
      }
      mZorro.setBaseCCDBPath(confCcdb.triggerPath.value);
    }
    if (confRct.useRctFlags.value) {
      mUseRctFlags = true;
      mRctFlagsChecker.init(confRct.label.value, confRct.useZdc.value, confRct.treatLimitedAcceptanceAsBad.value);
    }
    mGrpPath = confCcdb.grpPath.value;

    LOG(info) << "Initialize femto collision builder...";
    mProducedCollisions = utils::enableTable("FCols_001", confTable.produceCollisions.value, initContext);
    mProducedCollisionMasks = utils::enableTable("FColMasks_001", confTable.produceCollisionMasks.value, initContext);
    mProducedPositions = utils::enableTable("FColPos_001", confTable.producePositions.value, initContext);
    mProducedSphericities = utils::enableTable("FColSphericities_001", confTable.produceSphericities.value, initContext);
    mProducedMultiplicities = utils::enableTable("FColMults_001", confTable.produceMults.value, initContext);
    mProducedCentralities = utils::enableTable("FColCents_001", confTable.produceCents.value, initContext);
    mProduceQns = utils::enableTable("FColQnBins_001", confTable.produceQns.value, initContext);
    if (mProducedCollisions || mProducedCollisionMasks || mProducedPositions || mProducedSphericities || mProducedMultiplicities || mProducedCentralities) {
      mFillAnyTable = true;
      mCollisionSelection.printSelections(colSelsName, colSelsToString);
    } else {
      LOG(info) << "No tables configured";
    }
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void initCollision(T1& bc, T2& col, T3& tracks, T4& ccdb, T5& histRegistry)
  {
    if (mRunNumber != bc.runNumber()) {
      mRunNumber = bc.runNumber();
      static o2::parameters::GRPMagField* grpo = nullptr;
      grpo = ccdb->template getForRun<o2::parameters::GRPMagField>(mGrpPath, mRunNumber);
      if (grpo == nullptr) {
        LOG(fatal) << "GRP object not found for Run " << mRunNumber;
      }
      mMagField = static_cast<int>(grpo->getNominalL3Field()); // get magnetic field in kG

      if (mUseTrigger) {
        mZorro.initCCDB(ccdb.service, mRunNumber, bc.timestamp(), mTriggerNames);
        mZorro.populateHistRegistry(histRegistry, mRunNumber);
      }
    }

    mCollisionSelection.setMagneticField(mMagField);
    mCollisionSelection.setSphericity(tracks);
    mCollisionSelection.setMultiplicity<system>(col);
    mCollisionSelection.setCentrality<system>(col);

    std::vector<bool> triggerDecisions = {};
    if (mUseTrigger) {
      triggerDecisions = mZorro.getTriggerOfInterestResults(bc.globalBC());
    }

    mCollisionSelection.applySelections(col, triggerDecisions);
  }

  template <typename T1>
  bool checkCollision(T1 const& col)
  {
    // check RCT flags first
    if (mUseRctFlags && !mRctFlagsChecker(col)) {
      return false;
    }
    // make other checks
    return mCollisionSelection.checkFilters(col) &&
           mCollisionSelection.passesAllRequiredSelections();
  }

  template <modes::System system, typename T1, typename T2>
  void fillCollision(T1& collisionProducts, T2 const& col)
  {
    if (!mFillAnyTable) {
      return;
    }
    if (mProducedCollisions) {
      collisionProducts.producedCollision(col.posZ(),
                                          col.multNTracksPV(),
                                          mCollisionSelection.getCentrality(),
                                          static_cast<int8_t>(mCollisionSelection.getMagneticField()));
    }
    if (mProducedCollisionMasks) {
      collisionProducts.producedCollisionMask(mCollisionSelection.getBitmask());
    }
    if (mProducedPositions) {
      collisionProducts.producedPositions(col.posX(),
                                          col.posY());
    }
    if (mProducedSphericities) {
      collisionProducts.producedSphericities(mCollisionSelection.getSphericity());
    }
    if (mProducedMultiplicities) {
      collisionProducts.producedMultiplicityEstimators(
        col.multFT0A(),
        col.multFT0C(),
        col.multNTracksPVeta1(),
        col.multNTracksPVetaHalf(),
        col.trackOccupancyInTimeRange(),
        col.ft0cOccupancyInTimeRange());
    }
    if (mProducedCentralities) {
      collisionProducts.producedCentralityEstimators(
        col.centFT0A(),
        col.centFT0C());
    }

    // PbPb specific columns
    if constexpr (modes::isFlagSet(system, modes::System::kPbPb)) {
      if (mProduceQns) {
        collisionProducts.producedQns(utils::qn(col));
      }
    }
  }

 private:
  CollisionSelection mCollisionSelection;
  Zorro mZorro;
  bool mUseTrigger = false;
  int mRunNumber = -1;
  std::string mGrpPath = std::string("");
  int mMagField = 0;
  aod::rctsel::RCTFlagsChecker mRctFlagsChecker;
  bool mUseRctFlags = false;
  std::string mTriggerNames = std::string("");
  bool mFillAnyTable = false;
  bool mProducedCollisions = false;
  bool mProducedCollisionMasks = false;
  bool mProducedPositions = false;
  bool mProducedSphericities = false;
  bool mProducedMultiplicities = false;
  bool mProducedCentralities = false;
  bool mProduceQns = false;
};

struct CollisionBuilderDerivedToDerivedProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::StoredFCols> producedCollision;
  o2::framework::Produces<o2::aod::StoredFColMasks> producedCollisionMask;
};

class CollisionBuilderDerivedToDerived
{
 public:
  CollisionBuilderDerivedToDerived() = default;
  ~CollisionBuilderDerivedToDerived() = default;

  template <typename T1, typename T2>
  void processCollision(T1& col, T2& newCollisionTable)
  {
    newCollisionTable.producedCollision(col.posZ(),
                                        col.mult(),
                                        col.cent(),
                                        col.magField());
    newCollisionTable.producedCollisionMask(col.mask());
  }
};

} // namespace collisionbuilder
} // namespace o2::analysis::femto
;
#endif // PWGCF_FEMTO_CORE_COLLISIONBUILDER_H_
