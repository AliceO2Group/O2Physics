// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
/// \file lumiVertex.cxx
/// \brief Codes for the vdM scan analysis based on the Primary Vertex reconstruction. The code is based on the lumiTask in O2Physics/PWGMM/Lumi/Tasks/lumi.cxx.
///
/// \author Minjae Kim <minjae.kim@cern.ch>
/// \since Aug.06. 2026

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DataFormatsFT0/Digit.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PrimaryVertex.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <vector>

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, uint64_t);
DECLARE_SOA_COLUMN(VertexX, vertexX, double);
DECLARE_SOA_COLUMN(VertexY, vertexY, double);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, double);

DECLARE_SOA_COLUMN(VertexXX, vertexXX, double);
DECLARE_SOA_COLUMN(VertexYY, vertexYY, double);
DECLARE_SOA_COLUMN(VertexXY, vertexXY, double);

DECLARE_SOA_COLUMN(VertexChi2, vertexChi2, double);
DECLARE_SOA_COLUMN(NContrib, nContrib, int);
DECLARE_SOA_COLUMN(Bcid, bcid, int32_t);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(LocalBC, localBC, int32_t);
DECLARE_SOA_COLUMN(InputMask, inputMask, uint64_t);
DECLARE_SOA_COLUMN(HasFT0, hasFT0, bool);
DECLARE_SOA_COLUMN(HasFoundFT0, hasFoundFT0, bool);
DECLARE_SOA_COLUMN(HasFT0TVX, hasFT0TVX, bool);
DECLARE_SOA_COLUMN(IsCollisionTVX, isCollisionTVX, bool);
DECLARE_SOA_COLUMN(HasAnyFT0Trigger, hasAnyFT0Trigger, bool);
DECLARE_SOA_COLUMN(HasFT0AC, hasFT0AC, bool);
DECLARE_SOA_COLUMN(Ft0TriggerMask, ft0TriggerMask, uint8_t);
DECLARE_SOA_COLUMN(Ft0Rate, ft0Rate, float);
DECLARE_SOA_COLUMN(NCollisions, nCollisions, int32_t);
DECLARE_SOA_COLUMN(NSelectedCollisions, nSelectedCollisions, int32_t);
DECLARE_SOA_COLUMN(IsCollidingBC, isCollidingBC, bool);
DECLARE_SOA_COLUMN(IsFT0SuperLeading, isFT0SuperLeading, bool);
DECLARE_SOA_COLUMN(HasPreviousFT0Activity, hasPreviousFT0Activity, bool);      //! Whether a preceding FT0 CTP activity was observed
DECLARE_SOA_COLUMN(BcsSinceLastFT0Activity, bcsSinceLastFT0Activity, int64_t); //! -1 when the preceding activity is unknown
DECLARE_SOA_COLUMN(PvRefitValid, pvRefitValid, bool);
} // namespace full
DECLARE_SOA_TABLE(EventInfo, "AOD", "EventInfo", full::TimeStamp, full::VertexX,
                  full::VertexY, full::VertexZ,

                  full::VertexXX, full::VertexYY, full::VertexXY,

                  full::VertexChi2, full::NContrib, full::Bcid, full::RunNumber,
                  full::GlobalBC, full::HasFT0TVX, full::IsCollisionTVX,
                  full::HasAnyFT0Trigger,
                  full::IsCollidingBC, full::IsFT0SuperLeading,
                  full::HasPreviousFT0Activity, full::BcsSinceLastFT0Activity,
                  full::PvRefitValid);
DECLARE_SOA_TABLE(EventInfoBC, "AOD", "EventInfoBC", full::RunNumber,
                  full::TimeStamp, full::GlobalBC, full::LocalBC,
                  full::InputMask, full::HasFT0, full::HasFoundFT0,
                  full::HasFT0TVX, full::HasAnyFT0Trigger, full::HasFT0AC,
                  full::Ft0TriggerMask, full::Ft0Rate, full::NCollisions,
                  full::NSelectedCollisions, full::IsCollidingBC,
                  full::IsFT0SuperLeading, full::HasPreviousFT0Activity,
                  full::BcsSinceLastFT0Activity);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using CollisionsWithEvSels = soa::Join<aod::Collisions, aod::EvSels>;
using UnfilteredTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra>;

struct LumiVertex {
  Produces<o2::aod::EventInfo> rowEventInfo;
  Produces<o2::aod::EventInfoBC> rowEventInfoBC;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  const char* ccdbPathGrp = "GLO/Config/GRPMagField";
  const char* ccdbPathLut = "GLO/Param/MatLUT";
  const char* ccdburl = "http://alice-ccdb.cern.ch";
  int mFieldRunNumber = -1;
  o2::base::MatLayerCylSet* mMatLUT = nullptr;
  int mLHCIFRunNumber = -1;
  o2::parameters::GRPLHCIFData* mLHCIFData = nullptr;
  std::bitset<o2::constants::lhc::LHCMaxBunches> mCollidingBCPattern;

  bool mHasPreviousFT0Activity = false;
  uint64_t mGlobalBCOfLastFT0Activity = 0;
  bool mHasLastProcessedBC = false;
  uint64_t mLastProcessedGlobalBC = 0;
  uint64_t mSkippedZeroTimestampBCs = 0;

  ctpRateFetcher mRateFetcher;
  SliceCache cache;
  PresliceUnsorted<CollisionsWithEvSels> collisionsPerFoundBC =
    aod::evsel::foundBCId;
  Preslice<UnfilteredTracks> tracksPerCollision =
    aod::track::collisionId;

  struct : ConfigurableGroup {
    Configurable<uint64_t> ftts{"ftts", 0,
                                "Timestamp reference in ms; zero uses the fill start from GRPLHCIF"};
    Configurable<int> nContribMax{"nContribMax", 2500,
                                  "Maximum number of contributors"};
    Configurable<int> nContribMin{"nContribMin", 1,
                                  "Minimum number of contributors"};
    Configurable<float> minVertexChi2{"minVertexChi2", 0,
                                      "Minimum chi2 for PV selection"};
    Configurable<float> maxVertexChi2PerContributor{"maxVertexChi2PerContributor", 999,
                                                    "Maximum refitted-vertex chi2 per contributor"};
    Configurable<bool> doBasicCheck{"doBasicCheck", true,
                                    "Perform basic checks on the refitted vertex"};
    Configurable<bool> doCheckTVXCollision{"doCheckTVXCollision", true,
                                           "Require the collision-level kIsTriggerTVX selection bit"};
    Configurable<bool> doPVrefit{"doPVrefit", true,
                                 "Refit the primary vertex"};
    Configurable<bool> useMatCorrLUT{"useMatCorrLUT", true,
                                     "Use the CCDB material LUT during PV refit; false disables material correction"};
    Configurable<bool> doHistogramAnalyisis{"doHistogramAnalyisis", true,
                                            "Fill histograms for PV refit analysis"};
    Configurable<bool> doFetchFT0Rate{"doFetchFT0Rate", false,
                                      "Fetch T0VTX rate from CTP scalers and store it in the BC table"};
    Configurable<int> nBCsBeforeFT0SuperLeading{"nBCsBeforeFT0SuperLeading", 1,
                                                "Minimum global-BC distance from the previous FT0 CTP activity"};

  } basicConfig;

  struct : ConfigurableGroup {
    Configurable<bool> requirePVContributor{"requirePVContributor", true,
                                            "Use only tracks marked as PV contributors in the refit"};
    Configurable<float> itsTrackPt{"itsTrackPt", 0.5,
                                   "Minimum pt of tracks used for PV refit"};
    Configurable<int> itsNCls{"itsNCls", 5,
                              "Minimum number of ITS clusters for tracks used for PV refit"};
  } trackConfig;

  HistogramRegistry histos{
    "histos",
    {{"vertexx", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}},     //
     {"vertexy", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}},     //
     {"timestamp", "", {HistType::kTH1F, {{20000, 0, 2e7, "t"}}}}, //
     {"vertexx_timestamp",
      "",
      {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "x"}}}}, //
     {"vertexy_timestamp",
      "",
      {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "y"}}}},     //
     {"chisquare", "", {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}}, //

     {"vertexx_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}}, //
     {"vertexy_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}}, //
     {"vertexx_Refitted_timestamp",
      "",
      {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "x"}}}}, //
     {"vertexy_Refitted_timestamp",
      "",
      {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "y"}}}}, //
     {"chisquare_Refitted",
      "",
      {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}}, //

     {"vertexx_Refitted_vertexx",
      "",
      {HistType::kTH2F, {{1000, -1, 1, "x"}, {1000, -1, 1, "rx"}}}}, //
     {"vertexy_Refitted_vertexy",
      "",
      {HistType::kTH2F, {{1000, -1, 1, "y"}, {1000, -1, 1, "ry"}}}}, //

     {"bc/nStoredBCsVsBCID",
      "Stored BC rows;BC ID in orbit;Entries",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"bc/nSkippedZeroTimestamp",
      "BC rows skipped because timestamp is zero;;Entries",
      {HistType::kTH1D, {{1, 0.5, 1.5}}}},
     {"bc/nCollidingBCsVsBCID",
      "Stored colliding BC rows;BC ID in orbit;Entries",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"bc/nFT0ActivityBCsVsBCID",
      "Stored BC rows with any FT0 CTP activity;BC ID in orbit;Entries",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"bc/nSuperLeadingFT0BCsVsBCID",
      "TVX super-leading colliding BC rows;BC ID in orbit;Entries",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"bc/nNonSuperLeadingFT0BCsVsBCID",
      "TVX non-super-leading colliding BC rows;BC ID in orbit;Entries",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"bc/nNoPreviousFT0ActivityBCsVsBCID",
      "Colliding TVX BC rows without an observed preceding FT0 activity;BC ID in orbit;Entries",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},

     {"collision/nCollisionsVsBCID",
      "Collisions after the optional collision-TVX selection;found BC ID in orbit;Collisions",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/nSelectedCollisionsVsBCID",
      "Collisions passing the basic PV-refit selection;BC ID in orbit;Collisions",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/nSuperLeadingFT0CollisionsVsBCID",
      "Collisions in TVX super-leading BCs;BC ID in orbit;Collisions",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/nNonSuperLeadingFT0CollisionsVsBCID",
      "Collisions in TVX non-super-leading BCs;BC ID in orbit;Collisions",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/nNoPreviousFT0ActivityCollisionsVsBCID",
      "Collisions in TVX BCs without an observed preceding FT0 activity;BC ID in orbit;Collisions",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/nSelectedSuperLeadingFT0CollisionsVsBCID",
      "Selected collisions in TVX super-leading BCs;BC ID in orbit;Collisions",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/nSelectedNonSuperLeadingFT0CollisionsVsBCID",
      "Selected collisions in TVX non-super-leading BCs;BC ID in orbit;Collisions",
      {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanVertexXVsBCID",
      "Mean original vertex x by BC ID;BC ID in orbit;#LT x #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanVertexYVsBCID",
      "Mean original vertex y by BC ID;BC ID in orbit;#LT y #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanRefittedVertexXVsBCID",
      "Mean refitted vertex x by BC ID;BC ID in orbit;#LT x_{refit} #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanRefittedVertexYVsBCID",
      "Mean refitted vertex y by BC ID;BC ID in orbit;#LT y_{refit} #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanRefittedVertexXSuperLeadingVsBCID",
      "Mean refitted vertex x in TVX super-leading BCs;BC ID in orbit;#LT x_{refit} #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanRefittedVertexYSuperLeadingVsBCID",
      "Mean refitted vertex y in TVX super-leading BCs;BC ID in orbit;#LT y_{refit} #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanRefittedVertexXNonSuperLeadingVsBCID",
      "Mean refitted vertex x in TVX non-super-leading BCs;BC ID in orbit;#LT x_{refit} #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}},
     {"collision/meanRefittedVertexYNonSuperLeadingVsBCID",
      "Mean refitted vertex y in TVX non-super-leading BCs;BC ID in orbit;#LT y_{refit} #GT (cm)",
      {HistType::kTProfile, {{o2::constants::lhc::LHCMaxBunches, -0.5, o2::constants::lhc::LHCMaxBunches - 0.5, "BC ID in orbit"}}}}}};
  // CTP input 3 is MTVX; std::bitset uses a zero-based bit index.
  static constexpr uint8_t FT0TVXInputBit = 2;

  static bool hasFT0ACCoincidence(uint8_t triggerMask)
  {
    return (triggerMask & (1U << o2::ft0::Triggers::bitA)) != 0U &&
           (triggerMask & (1U << o2::ft0::Triggers::bitC)) != 0U;
  }

  static bool hasAnyFT0Trigger(const std::bitset<64>& ctpInputMask)
  {
    return ctpInputMask.test(0) || ctpInputMask.test(1) ||
           ctpInputMask.test(FT0TVXInputBit) || ctpInputMask.test(3) ||
           ctpInputMask.test(4);
  }

  struct FT0SuperLeadingInfo {
    bool isCollidingBC = false;
    bool isSuperLeading = false;
    bool hasPreviousFT0Activity = false;
    bool hasFT0Activity = false;
    bool hasFT0TVX = false;
    int64_t bcsSinceLastActivity = -1;
  };

  void resetFT0SuperLeadingHistory()
  {
    mHasPreviousFT0Activity = false;
    mGlobalBCOfLastFT0Activity = 0;
    mHasLastProcessedBC = false;
    mLastProcessedGlobalBC = 0;
  }

  void setLHCIFData(const auto& bc)
  {
    if (mLHCIFRunNumber == bc.runNumber()) {
      return;
    }
    if (mLHCIFRunNumber >= 0 && basicConfig.doHistogramAnalyisis) {
      LOGF(warn,
           "Run changed from %d to %d: BCID histograms aggregate runs. "
           "Use EventInfo tables or one run per job for run-separated results.",
           mLHCIFRunNumber, bc.runNumber());
    }

    mLHCIFData = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>(
      "GLO/Config/GRPLHCIF", bc.timestamp());
    if (mLHCIFData == nullptr) {
      LOGF(fatal,
           "GRPLHCIFData is not available for run=%d at timestamp=%llu",
           bc.runNumber(), bc.timestamp());
      return;
    }

    const auto beamPatternA =
      mLHCIFData->getBunchFilling().getBeamPattern(0);
    const auto beamPatternC =
      mLHCIFData->getBunchFilling().getBeamPattern(1);
    mCollidingBCPattern = beamPatternA & beamPatternC;
    mLHCIFRunNumber = bc.runNumber();
    resetFT0SuperLeadingHistory();

    LOGF(info, "Loaded filling scheme for run %d: %zu colliding BC IDs",
         mLHCIFRunNumber, mCollidingBCPattern.count());
  }

  uint64_t getRelativeTimestamp(const auto& bc) const
  {
    const uint64_t referenceTimestamp =
      basicConfig.ftts.value > 0
        ? basicConfig.ftts.value
        : static_cast<uint64_t>(mLHCIFData->getFillNumberTime());
    if (bc.timestamp() < referenceTimestamp) {
      return 0;
    }
    return bc.timestamp() - referenceTimestamp;
  }

  FT0SuperLeadingInfo classifyFT0SuperLeading(
    uint64_t globalBC,
    int32_t localBC,
    const std::bitset<64>& ctpInputMask)
  {
    if (mHasLastProcessedBC && globalBC <= mLastProcessedGlobalBC) {
      LOGF(warn,
           "Non-increasing global BC sequence (%llu after %llu); "
           "resetting FT0 super-leading history",
           globalBC, mLastProcessedGlobalBC);
      resetFT0SuperLeadingHistory();
    }

    FT0SuperLeadingInfo result;
    result.isCollidingBC = mCollidingBCPattern.test(localBC);
    result.hasFT0Activity = hasAnyFT0Trigger(ctpInputMask);
    result.hasFT0TVX = ctpInputMask.test(FT0TVXInputBit);
    result.hasPreviousFT0Activity =
      mHasPreviousFT0Activity && globalBC > mGlobalBCOfLastFT0Activity;

    if (result.hasPreviousFT0Activity) {
      result.bcsSinceLastActivity =
        static_cast<int64_t>(globalBC - mGlobalBCOfLastFT0Activity);
    }
    const bool hasSufficientFT0Gap =
      !result.hasPreviousFT0Activity ||
      result.bcsSinceLastActivity >=
        static_cast<int64_t>(basicConfig.nBCsBeforeFT0SuperLeading.value);
    result.isSuperLeading = result.hasFT0TVX &&
                            result.isCollidingBC &&
                            hasSufficientFT0Gap;

    // The current BC must be classified against preceding activity, so update
    // the state only after the classification above.
    if (result.hasFT0Activity) {
      mGlobalBCOfLastFT0Activity = globalBC;
      mHasPreviousFT0Activity = true;
    }
    mLastProcessedGlobalBC = globalBC;
    mHasLastProcessedBC = true;
    return result;
  }

  void updateFieldAndMaterial(const auto& bc)
  {
    if (mFieldRunNumber == bc.runNumber()) {
      return;
    }

    auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(
      ccdbPathGrp, bc.timestamp());
    if (grpo == nullptr) {
      LOGF(fatal,
           "GRP object is not available in CCDB for run=%d at timestamp=%llu",
           bc.runNumber(), bc.timestamp());
    }
    o2::base::Propagator::initFieldFromGRP(grpo);

    if (basicConfig.useMatCorrLUT) {
      auto* rawLUT = ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(
        ccdbPathLut, bc.timestamp());
      if (rawLUT == nullptr) {
        LOGF(fatal,
             "Material LUT is not available in CCDB for run=%d at timestamp=%llu",
             bc.runNumber(), bc.timestamp());
      }
      mMatLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(rawLUT);
      o2::base::Propagator::Instance()->setMatLUT(mMatLUT);
    }

    mFieldRunNumber = bc.runNumber();
  }

  void init(InitContext&)
  {
    if (basicConfig.nBCsBeforeFT0SuperLeading < 1) {
      LOGF(fatal,
           "nBCsBeforeFT0SuperLeading must be at least 1 (received %d)",
           basicConfig.nBCsBeforeFT0SuperLeading.value);
    }
    if (basicConfig.nContribMin < 0 ||
        basicConfig.nContribMax <= basicConfig.nContribMin) {
      LOG(fatal) << "Require 0 <= nContribMin < nContribMax";
    }
    if (basicConfig.maxVertexChi2PerContributor <= 0) {
      LOG(fatal) << "maxVertexChi2PerContributor must be positive";
    }
    if (!basicConfig.doPVrefit && basicConfig.doBasicCheck) {
      LOG(warn) << "doPVrefit=false with doBasicCheck=true will reject all "
                   "collision rows; set doBasicCheck=false to store them";
    }

    // A mean-vertex constraint requires loading and passing the CCDB mean
    // vertex object to PVertexer. This standalone task intentionally refits
    // without that external constraint.
    o2::conf::ConfigurableParam::updateFromString(
      "pvertexer.useMeanVertexConstraint=true");

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  bool processCollision(auto const& collision,
                        auto const& bc,
                        auto const& unfilteredTracks,
                        uint64_t relativeTimestamp,
                        uint64_t globalBC,
                        int32_t localBC,
                        const FT0SuperLeadingInfo& ft0SLInfo)
  {
    // Check only the collision-level TVX bit. EvSels::sel8() would also apply
    // the time-frame and ITS readout-frame border selections.
    const bool isCollisionTVX =
      collision.selection_bit(o2::aod::evsel::kIsTriggerTVX);
    if (basicConfig.doCheckTVXCollision && !isCollisionTVX) {
      return false;
    }

    std::vector<o2::track::TrackParCov> pvContributorTracks;
    int nContrib = 0;
    for (const auto& track : unfilteredTracks) {
      if (trackConfig.requirePVContributor && !track.isPVContributor()) {
        continue;
      }
      if (!track.hasITS()) {
        continue;
      }
      if (track.pt() < trackConfig.itsTrackPt ||
          track.itsNCls() < trackConfig.itsNCls) {
        continue;
      }
      pvContributorTracks.push_back(getTrackParCov(track));
      nContrib++;
    }

    o2::dataformats::VertexBase primaryVertex;
    primaryVertex.setX(collision.posX());
    primaryVertex.setY(collision.posY());
    primaryVertex.setZ(collision.posZ());
    primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(),
                         collision.covXZ(), collision.covYZ(), collision.covZZ());

    double chi2 = -1.;
    double refitX = -9999.;
    double refitY = -9999.;
    double refitZ = -9999.;
    double refitXX = -9999.;
    double refitYY = -9999.;
    double refitXY = -9999.;

    bool pvRefitValid = false;
    if (basicConfig.doPVrefit && !pvContributorTracks.empty()) {
      updateFieldAndMaterial(bc);
      o2::vertexing::PVertexer vertexer;
      vertexer.init();
      vertexer.setMatCorrType(
        basicConfig.useMatCorrLUT
          ? o2::base::Propagator::MatCorrType::USEMatCorrLUT
          : o2::base::Propagator::MatCorrType::USEMatCorrNONE);

      const bool pvRefitDoable =
        vertexer.prepareVertexRefit(pvContributorTracks, primaryVertex);
      nContrib = static_cast<int>(vertexer.getTracksPool().size());
      if (pvRefitDoable) {
        std::vector<bool> useTrackForPVRefit(pvContributorTracks.size(), true);
        const auto refittedVertex =
          vertexer.refitVertex(useTrackForPVRefit, primaryVertex);
        const double candidateChi2 = refittedVertex.getChi2();
        if (candidateChi2 >= 0. &&
            std::isfinite(candidateChi2) &&
            std::isfinite(refittedVertex.getX()) &&
            std::isfinite(refittedVertex.getY()) &&
            std::isfinite(refittedVertex.getZ())) {
          pvRefitValid = true;
          nContrib = refittedVertex.getNContributors();
          chi2 = candidateChi2;
          refitX = refittedVertex.getX();
          refitY = refittedVertex.getY();
          refitZ = refittedVertex.getZ();
          refitXX = refittedVertex.getSigmaX2();
          refitYY = refittedVertex.getSigmaY2();
          refitXY = refittedVertex.getSigmaXY();
        }
      }
    }

    const bool passesBasicCheck =
      pvRefitValid &&
      chi2 >= basicConfig.minVertexChi2 &&
      nContrib >= basicConfig.nContribMin &&
      nContrib <= basicConfig.nContribMax &&
      (chi2 / nContrib) <=
        basicConfig.maxVertexChi2PerContributor;
    const bool passesOutputSelection =
      !basicConfig.doBasicCheck || passesBasicCheck;
    if (passesOutputSelection) {
      rowEventInfo(relativeTimestamp, refitX, refitY, refitZ, refitXX,
                   refitYY, refitXY, chi2, nContrib, localBC,
                   static_cast<int32_t>(bc.runNumber()), globalBC,
                   ft0SLInfo.hasFT0TVX, isCollisionTVX,
                   ft0SLInfo.hasFT0Activity,
                   ft0SLInfo.isCollidingBC, ft0SLInfo.isSuperLeading,
                   ft0SLInfo.hasPreviousFT0Activity,
                   ft0SLInfo.bcsSinceLastActivity, pvRefitValid);
    }

    if (!basicConfig.doHistogramAnalyisis) {
      return passesOutputSelection;
    }

    histos.fill(HIST("collision/nCollisionsVsBCID"), localBC);
    const bool isCollidingTVX =
      ft0SLInfo.isCollidingBC && ft0SLInfo.hasFT0TVX && isCollisionTVX;
    if (isCollidingTVX) {
      if (!ft0SLInfo.hasPreviousFT0Activity) {
        histos.fill(
          HIST("collision/nNoPreviousFT0ActivityCollisionsVsBCID"),
          localBC);
      }
      if (ft0SLInfo.isSuperLeading) {
        histos.fill(
          HIST("collision/nSuperLeadingFT0CollisionsVsBCID"), localBC);
      } else {
        histos.fill(
          HIST("collision/nNonSuperLeadingFT0CollisionsVsBCID"),
          localBC);
      }
    }
    if (passesBasicCheck) {
      histos.fill(HIST("collision/nSelectedCollisionsVsBCID"), localBC);
      if (isCollidingTVX) {
        if (ft0SLInfo.isSuperLeading) {
          histos.fill(
            HIST("collision/nSelectedSuperLeadingFT0CollisionsVsBCID"),
            localBC);
        } else {
          histos.fill(
            HIST("collision/nSelectedNonSuperLeadingFT0CollisionsVsBCID"),
            localBC);
        }
      }
    }

    if (pvRefitValid) {
      histos.fill(HIST("chisquare_Refitted"), chi2);
    }
    const bool passesRefitHistogramSelection = passesBasicCheck;
    if (passesRefitHistogramSelection) {
      histos.fill(HIST("vertexx_Refitted"), refitX);
      histos.fill(HIST("vertexy_Refitted"), refitY);
      histos.fill(HIST("vertexx_Refitted_timestamp"), relativeTimestamp,
                  refitX);
      histos.fill(HIST("vertexy_Refitted_timestamp"), relativeTimestamp,
                  refitY);
      histos.fill(HIST("collision/meanRefittedVertexXVsBCID"), localBC,
                  refitX);
      histos.fill(HIST("collision/meanRefittedVertexYVsBCID"), localBC,
                  refitY);

      if (isCollidingTVX) {
        if (ft0SLInfo.isSuperLeading) {
          histos.fill(
            HIST("collision/meanRefittedVertexXSuperLeadingVsBCID"),
            localBC, refitX);
          histos.fill(
            HIST("collision/meanRefittedVertexYSuperLeadingVsBCID"),
            localBC, refitY);
        } else {
          histos.fill(
            HIST("collision/meanRefittedVertexXNonSuperLeadingVsBCID"),
            localBC, refitX);
          histos.fill(
            HIST("collision/meanRefittedVertexYNonSuperLeadingVsBCID"),
            localBC, refitY);
        }
      }
    }

    histos.fill(HIST("chisquare"), collision.chi2());
    const int originalNContrib = collision.numContrib();
    const bool passesOriginalHistogramSelection =
      originalNContrib > 0 &&
      originalNContrib >= basicConfig.nContribMin &&
      originalNContrib <= basicConfig.nContribMax &&
      (collision.chi2() / originalNContrib) <=
        basicConfig.maxVertexChi2PerContributor;
    if (passesOriginalHistogramSelection) {
      histos.fill(HIST("vertexx"), collision.posX());
      histos.fill(HIST("vertexy"), collision.posY());
      histos.fill(HIST("timestamp"), relativeTimestamp);
      histos.fill(HIST("vertexx_timestamp"), relativeTimestamp,
                  collision.posX());
      histos.fill(HIST("vertexy_timestamp"), relativeTimestamp,
                  collision.posY());
      histos.fill(HIST("collision/meanVertexXVsBCID"), localBC,
                  collision.posX());
      histos.fill(HIST("collision/meanVertexYVsBCID"), localBC,
                  collision.posY());
    }

    if (passesOriginalHistogramSelection &&
        passesRefitHistogramSelection) {
      histos.fill(HIST("vertexx_Refitted_vertexx"), collision.posX(),
                  refitX);
      histos.fill(HIST("vertexy_Refitted_vertexy"), collision.posY(),
                  refitY);
    }
    return passesOutputSelection;
  }

  void process(MyBCs const& bcs, CollisionsWithEvSels const& collisions,
               UnfilteredTracks const& unfilteredTracks, aod::FT0s const&)
  {
    for (const auto& bc : bcs) {
      if (bc.timestamp() == 0) {
        mSkippedZeroTimestampBCs++;
        resetFT0SuperLeadingHistory();
        if (basicConfig.doHistogramAnalyisis) {
          histos.fill(HIST("bc/nSkippedZeroTimestamp"), 1.);
        }
        if (mSkippedZeroTimestampBCs == 1) {
          LOG(warn) << "Skipping BC rows with timestamp zero and resetting "
                       "the FT0 super-leading history";
        }
        continue;
      }

      setLHCIFData(bc);

      const uint64_t relativeTimestamp = getRelativeTimestamp(bc);
      const uint64_t globalBC = bc.globalBC();
      const auto localBC = static_cast<int32_t>(globalBC % o2::constants::lhc::LHCMaxBunches);
      const std::bitset<64> ctpInputMask(bc.inputMask());
      const bool hasFT0 = bc.has_ft0();
      const uint8_t ft0TriggerMask =
        hasFT0 ? bc.ft0().triggerMask() : 0;
      const auto ft0SLInfo =
        classifyFT0SuperLeading(globalBC, localBC, ctpInputMask);

      const bool hasFoundFT0 = bc.has_foundFT0();
      const bool hasFT0TVX = ft0SLInfo.hasFT0TVX;
      const bool hasFT0AC = hasFT0ACCoincidence(ft0TriggerMask);
      const float ft0Rate = basicConfig.doFetchFT0Rate ? static_cast<float>(mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), "T0VTX", true) * 1.e-3) : -1.f;

      // EvSel TVX is evaluated on the collision's found BC, so group
      // collisions by foundBCId rather than the nominal collision::bcId.
      auto collisionsInBC =
        collisions.sliceBy(collisionsPerFoundBC, bc.globalIndex());

      if (basicConfig.doHistogramAnalyisis) {
        histos.fill(HIST("bc/nStoredBCsVsBCID"), localBC);
        if (ft0SLInfo.hasFT0Activity) {
          histos.fill(HIST("bc/nFT0ActivityBCsVsBCID"), localBC);
        }
        if (ft0SLInfo.isCollidingBC) {
          histos.fill(HIST("bc/nCollidingBCsVsBCID"), localBC);
          if (ft0SLInfo.hasFT0TVX) {
            if (!ft0SLInfo.hasPreviousFT0Activity) {
              histos.fill(
                HIST("bc/nNoPreviousFT0ActivityBCsVsBCID"), localBC);
            }
            if (ft0SLInfo.isSuperLeading) {
              histos.fill(
                HIST("bc/nSuperLeadingFT0BCsVsBCID"), localBC);
            } else {
              histos.fill(
                HIST("bc/nNonSuperLeadingFT0BCsVsBCID"), localBC);
            }
          }
        }
      }

      int32_t nSelectedCollisions = 0;
      for (const auto& collision : collisionsInBC) {
        auto tracksInCollision =
          unfilteredTracks.sliceBy(
            tracksPerCollision,
            collision.globalIndex());

        if (processCollision(collision, bc, tracksInCollision,
                             relativeTimestamp, globalBC, localBC,
                             ft0SLInfo)) {
          nSelectedCollisions++;
        }
      }

      rowEventInfoBC(
        static_cast<int32_t>(bc.runNumber()), relativeTimestamp, globalBC,
        localBC, bc.inputMask(), hasFT0, hasFoundFT0, hasFT0TVX,
        ft0SLInfo.hasFT0Activity, hasFT0AC, ft0TriggerMask, ft0Rate,
        static_cast<int32_t>(collisionsInBC.size()), nSelectedCollisions,
        ft0SLInfo.isCollidingBC, ft0SLInfo.isSuperLeading,
        ft0SLInfo.hasPreviousFT0Activity,
        ft0SLInfo.bcsSinceLastActivity);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w{adaptAnalysisTask<LumiVertex>(cfgc)};
  return w;
}
