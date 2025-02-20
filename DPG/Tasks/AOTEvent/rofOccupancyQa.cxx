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

/// \file rofOccupancyQa.cxx
/// \brief ROF occupancy QA task
///
/// \author Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include <vector>

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "Framework/HistogramRegistry.h"
// #include "DataFormatsParameters/GRPLHCIFData.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

using BCsWithBcSelsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;

struct RofOccupancyQaTask {
  // configurables for occupancy-based event selection
  Configurable<float> confTimeIntervalForOccupancyCalculationMin{"TimeIntervalForOccupancyCalculationMin", -40, "Min time diff window for TPC occupancy calculation, us"};                      // o2-linter: disable=name/configurable
  Configurable<float> confTimeIntervalForOccupancyCalculationMax{"TimeIntervalForOccupancyCalculationMax", 100, "Max time diff window for TPC occupancy calculation, us"};                      // o2-linter: disable=name/configurable
  Configurable<float> confTimeRangeVetoOnCollStandard{"TimeRangeVetoOnCollStandard", 10.0, "Exclusion of a collision if there are other collisions nearby, +/- us"};                            // o2-linter: disable=name/configurable
  Configurable<float> confTimeRangeVetoOnCollNarrow{"TimeRangeVetoOnCollNarrow", 2.0, "Exclusion of a collision if there are other collisions nearby, +/- us"};                                 // o2-linter: disable=name/configurable
  Configurable<int> confNtracksCutVetoOnCollInTimeRange{"NtracksCutVetoOnCollInTimeRange", 800, "Max allowed N tracks (PV contributors) for each nearby collision in +/- time range"};          // o2-linter: disable=name/configurable
  Configurable<float> confEpsilonDistanceForVzDependentVetoTPC{"EpsilonDistanceForVzDependentVetoTPC", 2.5, "Epsilon for vZ-dependent veto on drifting TPC tracks from nearby collisions, cm"}; // o2-linter: disable=name/configurable
  Configurable<float> confFT0CamplCutVetoOnCollInROF{"FT0CamplPerCollCutVetoOnCollInROF", 5000, "Max allowed FT0C amplitude for each nearby collision inside this ITS ROF"};                    // o2-linter: disable=name/configurable
  Configurable<float> confEpsilonVzDiffVetoInROF{"EpsilonVzDiffVetoInROF", 0.3, "Minumum distance to nearby collisions along z inside this ITS ROF, cm"};                                       // o2-linter: disable=name/configurable
  Configurable<bool> confUseWeightsForOccupancyVariable{"UseWeightsForOccupancyEstimator", 1, "Use or not the delta-time weights for the occupancy estimator"};                                 // o2-linter: disable=name/configurable
  Configurable<float> confFactorForHistRange{"kFactorForHistRange", 1.0, "To change axes b/n pp and Pb-Pb"};                                                                                    // o2-linter: disable=name/configurable

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int lastRun = -1; // last run number (needed to access ccdb only if run!=lastRun)
  // std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB; // bc pattern of colliding bunches

  int64_t bcSOR = -1;     // global bc of the start of the first orbit
  int64_t nBCsPerTF = -1; // duration of TF in bcs, should be 128*3564 or 32*3564
  int rofOffset = -1;     // ITS ROF offset, in bc
  int rofLength = -1;     // ITS ROF length, in bc

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    float k = confFactorForHistRange;
    histos.add("hDeltaTime", "", kTH1D, {{1500, -50, 100}});
    histos.add("hDeltaTimeAboveNtracksCut", "", kTH1D, {{1500, -50, 100}});
    histos.add("hDeltaTime_vZ10cm", "", kTH1D, {{1500, -50, 100}});
    histos.add("hDeltaTime_sel8", "", kTH1D, {{1500, -50, 100}});
    histos.add("hDeltaTime_sel8_vZ10cm", "", kTH1D, {{1500, -50, 100}});
    histos.add("hDeltaTimeAboveNtracksCut_sel8_vZ10cm", "", kTH1D, {{1500, -50, 100}});

    histos.add("hOccupancyWeights", "", kTH1D, {{150, -50, 100}});
    histos.add("hOccupancyByTracks", "", kTH1D, {{250, 0., 25000 * k}});
    histos.add("hOccupancyByFT0C", "", kTH1D, {{250, 0., 2.5e5 * k}});
    histos.add("hOccupancyByTrInROF", "", kTH1D, {{250, 0., 25000 * k}});
    histos.add("hOccupancyByFT0C_vs_ByTracks", "", kTH2D, {{500, 0., 25000 * k}, {500, 0., 2.5e5 * k}});
    histos.add("hOccupancyByFT0C_vs_ByTracks_vZ_TF_ROF_border_cuts", "", kTH2D, {{500, 0., 25000 * k}, {500, 0., 2.5e5 * k}});
    histos.add("hOccupancyByFT0C_vs_ByTracks_afterNarrowDeltaTimeCut", "", kTH2D, {{500, 0., 25000 * k}, {500, 0., 2.5e5 * k}});
    histos.add("hOccupancyByFT0C_vs_ByTracks_afterStrictDeltaTimeCut", "", kTH2D, {{500, 0., 25000 * k}, {500, 0., 2.5e5 * k}});
    histos.add("hOccupancyByFT0C_vs_ByTracks_afterStandardDeltaTimeCut", "", kTH2D, {{500, 0., 25000 * k}, {500, 0., 2.5e5 * k}});
    histos.add("hOccupancyByFT0C_vs_ByTracks_afterVzDependentDeltaTimeCut", "", kTH2D, {{500, 0., 25000 * k}, {500, 0., 2.5e5 * k}});

    histos.add("hOccupancyByTracks_CROSSCHECK", "", kTH1D, {{250, 0., 25000 * k}});
    histos.add("hOccupancyByFT0C_CROSSCHECK", "", kTH1D, {{250, 0., 2.5e5 * k}});

    // this ev nITStr vs FT0C
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/all", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/vZ_TF_ROF_border_cuts", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/afterNarrowDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/afterStrictDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/afterStandardDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/afterVzDependentDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C/kNoCollInRofStrict", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/kNoCollInRofStandard", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/kNoCollInRofWithCloseVz", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C/NarrowDeltaCut_StdTimeAndRofCuts", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C/occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/hThisEvITSTPCTr_vs_ThisEvFT0C_occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/hThisEvITSTPCTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    // CROSS-CHECK SEL BITS:
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_afterNarrowDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_afterStrictDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_afterStandardDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_kNoCollInRofStrict", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_kNoCollInRofStandard", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    // histos.add("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_kNoCollInRofWithCloseVz", "", kTH2D, {{250, 0., 1e5*k}, {250, 0., 10000*k}});

    // this ev nITSTPCtr vs nITStr
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/all", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/vZ_TF_ROF_border_cuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/afterNarrowDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/afterStrictDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/afterStandardDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/afterVzDependentDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/kNoCollInRofStrict", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/kNoCollInRofStandard", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/kNoCollInRofWithCloseVz", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/NarrowDeltaCut_StdTimeAndRofCuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/occupBelow2000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/occupBelow2000_NarrowDeltaCut_StdTimeAndRofCuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr/occupBelow2000_StrictDeltaTimeCutAndRofCuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITStr_vs_vZcut", "", kTH2D, {{200, 0., 10.}, {200, 0., 8000 * k}});
    histos.add("hThisEvITSTPCtr_vs_vZcut", "", kTH2D, {{200, 0., 10.}, {200, 0., 5000 * k}});

    histos.add("hThisEvITStr_vs_vZ", "", kTH2D, {{400, -20, 20.}, {200, 0., 8000 * k}});
    histos.add("hThisEvITSTPCtr_vs_vZ", "", kTH2D, {{400, -20, 20.}, {200, 0., 5000 * k}});

    histos.add("hVz", "", kTH1D, {{1600, -40, 40.}});
    histos.add("hDeltaVz", "", kTH1D, {{1600, -40, 40.}});
    histos.add("hDeltaVzAfterCuts", "", kTH1D, {{1600, -40, 40.}});
    histos.add("hDeltaVzAfterTFandROFborderCuts", "", kTH1D, {{1600, -40, 40.}});
    histos.add("hDeltaVzGivenCollAbove100NearbyBelow100", "", kTH1D, {{1600, -40, 40.}});
    histos.add("hDeltaVzGivenCollBelow100NearbyAbove100", "", kTH1D, {{1600, -40, 40.}});
    histos.add("hDeltaVzAfterAllCuts", "", kTH1D, {{1600, -40, 40.}});

    //
    histos.add("hDeltaVzVsDeltaTime1", "", kTH2D, {{400, -25, 25}, {100, -20, 20}});
    histos.add("hDeltaVzVsDeltaTime2", "", kTH2D, {{400, -25, 25}, {100, -20, 20}});
    histos.add("hDeltaVzVsDeltaTime3", "", kTH2D, {{400, -25, 25}, {100, -20, 20}});
    histos.add("hDeltaVzVsDeltaTime4", "", kTH2D, {{400, -25, 25}, {100, -20, 20}});

    histos.add("hEtaVz02", "", kTH1D, {{500, -2.5, 2.5}});
    histos.add("hEtaVzPlus10", "", kTH1D, {{500, -2.5, 2.5}});
    histos.add("hEtaVzMinus10", "", kTH1D, {{500, -2.5, 2.5}});
    histos.add("hEtaVzPlus15", "", kTH1D, {{500, -2.5, 2.5}});
    histos.add("hEtaVzMinus15", "", kTH1D, {{500, -2.5, 2.5}});
    histos.add("hEtaVsVz", "", kTH2D, {{250, -25, 25}, {250, -2.5, 2.5}});
    histos.add("hNPVcontribVsVz", "", kTH2D, {{250, -25, 25}, {500, 0., 8000}});
    histos.add("hNPVcontribVsVz_eta08", "", kTH2D, {{250, -25, 25}, {500, 0., 8000}});

    //
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterStrictDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterStrictDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterStrictDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterStrictDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterStandardDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterStandardDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterStandardDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterStandardDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterVzDependentDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterVzDependentDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("afterVzDependentDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("afterVzDependentDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStrict/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStrict/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStrict/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStrict/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStandard/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStandard/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStandard/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofStandard/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTPCTr_vs_occupancyByTracks", "", kTH2D, {{250, 0., 25000 * k}, {250, 0., 8000}});
    histos.add("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});
    histos.add("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTPCTr_vs_occupancyByFT0C", "", kTH2D, {{250, 0., 2.5e5 * k}, {250, 0., 8000}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF", "", kTH2D, {{250, 0., 12000 * k}, {250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyInROF", "", kTH2D, {{250, 0., 12000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF", "", kTH2D, {{250, 0., 12000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyInROF", "", kTH2D, {{250, 0., 12000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_HasNeighbours", "", kTH2D, {{250, 0., 12000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_HasNeighbours", "", kTH2D, {{250, 0., 12000 * k}, {250, 0., 8000 * k}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyInROF_HasNeighbours", "", kTH2D, {{250, 0., 12000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyFT0CInROF", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyFT0CInROF", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyFT0CInROF", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyFT0CInROF", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 8000 * k}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 100000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 100000 * k}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 100000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyInROF_2coll", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvFT0C_vs_occupancyFT0CInROF_2coll", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 100000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvFT0C_vs_occupancyFT0CInROF_2coll", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 100000 * k}});
    histos.add("kNoCollInRofWithCloseVz/hThisEvFT0C_vs_occupancyFT0CInROF_2coll", "", kTH2D, {{250, 0., 100000 * k}, {250, 0., 100000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInAnotherEarlierROF_1collPerROF", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInAnotherEarlierROF_1collPerROF", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    // histos.add("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyInAnotherEarlierROF_1collPerROF", "", kTH2D, {{250, 0., 8000*k}, {250, 0., 8000*k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPrevPrevROF_1collPerROF", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    // coll on x axis always has more tracks:
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    // 2,3,4 colls in ROF
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_noVzCutOnOtherVertices", "", kTH2D, {{500, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_3coll_noVzCutOnOtherVertices", "", kTH2D, {{500, 0., 20000 * k}, {250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_4coll_noVzCutOnOtherVertices", "", kTH2D, {{500, 0., 20000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_1coll_in_ROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_2coll_in_ROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_3coll_in_ROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_4coll_in_ROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_5collOrMore_in_ROF", "", kTH1D, {{250, 0., 8000 * k}});

    // 1D
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_allOccup_2coll_inROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_lowOccup_2coll_inROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_highOccup_2coll_inROF", "", kTH1D, {{250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_allOccup_1collPerROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_lowOccup_1collPerROF", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_highOccup_1collPerROF", "", kTH1D, {{250, 0., 8000 * k}});

    // now with the ratio on y-axis:
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});

    histos.add("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_inROF", "", kTH1D, {{500, 0., 15000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_thisROFprevROF", "", kTH1D, {{500, 0., 15000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_thisROFprevPrevROF", "", kTH1D, {{500, 0., 15000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_thisROFearlierThanPrevPrevROF", "", kTH1D, {{500, 0., 15000 * k}});

    //
    histos.add("hNcollPerROF", "", kTH1D, {{16, -0.5, 15.5}});

    // ROF-by-ROF study:
    histos.add("ROFbyROF/nPV_vs_ROFid", "", kTH2D, {{800, 0., 8000 * k}, {10, -0.5, 9.5}});
    histos.add("ROFbyROF/nPV_vs_subROFid", "", kTH2D, {{800, 0., 8000 * k}, {20, -0.5, 19.5}});

    histos.add("ROFbyROF/FT0C_vs_ROFid", "", kTH2D, {{800, 0., 80000 * k}, {10, -0.5, 9.5}});
    histos.add("ROFbyROF/FT0C_vs_subROFid", "", kTH2D, {{800, 0., 80000 * k}, {20, -0.5, 19.5}});

    histos.add("ROFbyROF/nPV_00x00", "", kTH1D, {{250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_10x00", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_01x00", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_00x10", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_00x01", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_11x00", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_01x10", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_00x11", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_10x00_nearbyByFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_01x00_nearbyByFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_00x10_nearbyByFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_00x01_nearbyByFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_10x00_thisFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_01x00_thisFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_00x10_thisFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_00x01_thisFT0C", "", kTH2D, {{250, 0., 80000 * k}, {250, 0., 8000 * k}});

    // histos.add("ROFbyROF/nPV_11x11", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    // ### sub-ROFs:
    histos.add("ROFbyROF/nPV_0_x00_0", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_0x0_0", "", kTH1D, {{250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_00x_0", "", kTH1D, {{250, 0., 8000 * k}});

    // corr with prev subROFs:
    histos.add("ROFbyROF/nPV_0_00x_y00_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_0x0_y00_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_x00_y00_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_0_00x_0y0_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_0x0_0y0_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_x00_0y0_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_0_00x_00y_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_0x0_00y_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_0_x00_00y_000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    // corr with next subROFs:
    histos.add("ROFbyROF/nPV_000_y00_00x_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_000_y00_0x0_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_000_y00_x00_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_000_0y0_00x_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_000_0y0_0x0_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_000_0y0_x00_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("ROFbyROF/nPV_000_00y_00x_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_000_00y_0x0_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("ROFbyROF/nPV_000_00y_x00_0", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    // #### new occupancy studies
    histos.add("nPV_vs_occupancyByTracks/sel8", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInTimeRangeNarrow", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInTimeRangeStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInTimeRangeStandard", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInRofStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInRofStandard", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInTimeAndRofStandard", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInTimeAndRofStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/NoCollInTimeAndRofStrict_vZ_5cm", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/kNoHighMultCollInPrevRof", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});
    histos.add("nPV_vs_occupancyByTracks/kNoHighMultCollInPrevRofAndRofStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 25000 * k}});

    histos.add("nPV_vs_occupancyByFT0C/sel8", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInTimeRangeNarrow", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInTimeRangeStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInTimeRangeStandard", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInRofStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInRofStandard", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInTimeAndRofStandard", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInTimeAndRofStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/NoCollInTimeAndRofStrict_vZ_5cm", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/kNoHighMultCollInPrevRof", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
    histos.add("nPV_vs_occupancyByFT0C/kNoHighMultCollInPrevRofAndRofStrict", "", kTH2D, {{125, 0., 8000 * k}, {100, 0., 2.5e5 * k}});
  }

  Partition<FullTracksIU> pvTracks = ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  // Partition<FullTracksIU> pvTracks = ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Preslice<FullTracksIU> perCollision = aod::track::collisionId;

  using ColEvSels = soa::Join<aod::Collisions, aod::EvSels>; //, aod::Mults, aod::CentFT0Cs>;
  void processRun3(ColEvSels const& cols, FullTracksIU const&, BCsWithBcSelsRun3 const& bcs, aod::FT0s const&)
  {
    int run = bcs.iteratorAt(0).runNumber();
    // extract bc pattern from CCDB for data or anchored MC only
    if (run != lastRun && run >= 500000) {
      lastRun = run;
      auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run);
      // first bc of the first orbit
      bcSOR = runInfo.orbitSOR * o2::constants::lhc::LHCMaxBunches;
      // duration of TF in bcs
      nBCsPerTF = runInfo.orbitsPerTF * o2::constants::lhc::LHCMaxBunches;

      // extract ITS ROF parameters
      int64_t ts = bcs.iteratorAt(0).timestamp();
      auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
      rofOffset = alppar->roFrameBiasInBC;
      rofLength = alppar->roFrameLengthInBC;
      LOGP(info, "rofOffset={} rofLength={}", rofOffset, rofLength);
    }

    std::vector<int> vTracksITS567perColl(cols.size(), 0);                                    // counter of tracks per collision for occupancy studies
    std::vector<int> vTracksITSTPCperColl(cols.size(), 0);                                    // counter of tracks per collision for occupancy studies
    std::vector<int> vTracksITS567eta08perColl(cols.size(), 0);                               // counter of tracks per collision for occupancy studies
    std::vector<float> vAmpFT0CperColl(cols.size(), 0);                                       // amplitude FT0C per collision
    std::vector<bool> vIsFullInfoForOccupancy(cols.size(), 0);                                // info for occupancy in +/- windows is available (i.e. a given coll is not too close to the TF borders)
    const float timeWinOccupancyCalcMinNS = confTimeIntervalForOccupancyCalculationMin * 1e3; // ns
    const float timeWinOccupancyCalcMaxNS = confTimeIntervalForOccupancyCalculationMax * 1e3; // ns
    std::vector<bool> vIsVertexITSTPC(cols.size(), 0);                                        // at least one of vertex contributors is ITS-TPC track
    std::vector<bool> vIsVertexTOFmatched(cols.size(), 0);                                    // at least one of vertex contributors is matched to TOF
    std::vector<bool> vIsVertexTRDmatched(cols.size(), 0);                                    // at least one of vertex contributors is matched to TRD

    // std::vector<int> vCollisionsPerBc(bcs.size(), 0);    // counter of collisions per found bc for pileup checks
    std::vector<int> vFoundBCindex(cols.size(), -1);     // indices of found bcs
    std::vector<int64_t> vFoundGlobalBC(cols.size(), 0); // global BCs for collisions

    std::vector<float> vCollVz(cols.size(), 0); // vector with vZ positions for each collision
    std::vector<bool> vIsSel8(cols.size(), 0);
    std::vector<bool> vCombCond(cols.size(), 0);

    std::vector<int> vCollRofId(cols.size(), 0);            // rof Id for each collision
    std::vector<int> vCollRofIdPerOrbit(cols.size(), 0);    // rof Id for each collision, per orbit
    std::vector<int> vCollRofSubId(cols.size(), 0);         // rof sub-Id for each collision
    std::vector<int> vCollRofSubIdPerOrbit(cols.size(), 0); // rof sub-Id for each collision, per orbit

    // first loop over collisions - collecting info
    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      // auto bc = col.bc_as<BCsWithBcSelsRun3>();
      const auto& bc = col.foundBC_as<BCsWithBcSelsRun3>();
      int64_t globalBC = bc.globalBC();

      int32_t foundBC = bc.globalIndex();
      vFoundBCindex[colIndex] = foundBC;
      vFoundGlobalBC[colIndex] = globalBC; // bc.globalBC();

      if (bc.has_foundFT0())
        vAmpFT0CperColl[colIndex] = bc.foundFT0().sumAmpC();

      vCollVz[colIndex] = col.posZ();
      vIsSel8[colIndex] = col.sel8();
      vCombCond[colIndex] = vIsSel8[colIndex] && (std::fabs(vCollVz[colIndex]) < 8) && (vAmpFT0CperColl[colIndex] > 500 /* a.u.*/);

      int bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
      vIsFullInfoForOccupancy[colIndex] = ((bcInTF - 300) * bcNS > -timeWinOccupancyCalcMinNS) && ((nBCsPerTF - 4000 - bcInTF) * bcNS > timeWinOccupancyCalcMaxNS) ? true : false;

      // int64_t rofId = (globalBC + 3564 - rofOffset) / rofLength;
      int rofId = (bcInTF - rofOffset) / rofLength;
      vCollRofId[colIndex] = rofId;

      int rofIdPerOrbit = rofId % (3564 / rofLength);
      vCollRofIdPerOrbit[colIndex] = rofIdPerOrbit;

      int bcInITSROF = (globalBC + 3564 - rofOffset) % rofLength;
      int subRofId = bcInITSROF / (rofLength / 3);
      vCollRofSubId[colIndex] = subRofId;
      vCollRofSubIdPerOrbit[colIndex] = 3 * rofIdPerOrbit + subRofId;
      // LOGP(info, ">> rofId={} rofIdPerOrbit={} subRofId={} vCollRofSubId={}", rofId, rofIdPerOrbit, subRofId, vCollRofSubId[colIndex]);

      auto colPvTracks = pvTracks.sliceBy(perCollision, col.globalIndex());

      for (const auto& track : colPvTracks) {
        if (track.itsNCls() >= 5) {
          vTracksITS567perColl[colIndex]++;
          if (std::fabs(track.eta()) < 0.8)
            vTracksITS567eta08perColl[colIndex]++;
          if (track.tpcNClsFound() > 70)
            vTracksITSTPCperColl[colIndex]++;
          if (std::fabs(col.posZ()) < 1)
            histos.fill(HIST("hEtaVz02"), track.eta());
          else if (col.posZ() > 9.5 && col.posZ() < 10.5)
            histos.fill(HIST("hEtaVzPlus10"), track.eta());
          else if (col.posZ() > -10.5 && col.posZ() < -9.5)
            histos.fill(HIST("hEtaVzMinus10"), track.eta());
          else if (col.posZ() > 14.5 && col.posZ() < 15.5)
            histos.fill(HIST("hEtaVzPlus15"), track.eta());
          else if (col.posZ() > -15.5 && col.posZ() < -14.5)
            histos.fill(HIST("hEtaVzMinus15"), track.eta());

          histos.fill(HIST("hEtaVsVz"), col.posZ(), track.eta());
        }
        if (track.hasTRD())
          vIsVertexTRDmatched[colIndex] = 1;
        if (track.hasTPC())
          vIsVertexITSTPC[colIndex] = 1;
        if (track.hasTOF()) {
          vIsVertexTOFmatched[colIndex] = 1;
        }
      }

      if (col.sel8()) {
        histos.fill(HIST("hNPVcontribVsVz"), col.posZ(), vTracksITS567perColl[colIndex]);
        histos.fill(HIST("hNPVcontribVsVz_eta08"), col.posZ(), vTracksITS567eta08perColl[colIndex]);
      }
    }

    // ROF-by-ROF study:
    int nColls = vCombCond.size();
    for (const auto& col : cols) {
      int32_t k = col.globalIndex();

      if (k - 2 < 0 || k + 2 > nColls - 1)
        continue;

      // the "purest case"
      if (vCombCond[k]) {
        if (vCollRofId[k - 1] < vCollRofId[k] - 2 && /* next coll is far */ vCollRofId[k + 1] > vCollRofId[k] + 2) // 00x00
          histos.fill(HIST("ROFbyROF/nPV_00x00"), vTracksITS567perColl[k]);

        if (vCollRofId[k - 1] < vCollRofId[k] - 1 && /* next coll is far */ vCollRofId[k + 1] > vCollRofId[k] + 1) // 0x0
        {
          if (vCollRofSubId[k] == 0)
            histos.fill(HIST("ROFbyROF/nPV_0_x00_0"), vTracksITS567perColl[k]);
          if (vCollRofSubId[k] == 1)
            histos.fill(HIST("ROFbyROF/nPV_0_0x0_0"), vTracksITS567perColl[k]);
          if (vCollRofSubId[k] == 2)
            histos.fill(HIST("ROFbyROF/nPV_0_00x_0"), vTracksITS567perColl[k]);
        }
      }

      // prev 1 coll
      if (vCombCond[k] && vCombCond[k - 1]) {
        if (vCollRofId[k - 2] < vCollRofId[k] - 2 && vCollRofId[k - 1] == vCollRofId[k] - 2 && /* next coll is far */ vCollRofId[k + 1] > vCollRofId[k] + 2) // 10x00
        {
          histos.fill(HIST("ROFbyROF/nPV_10x00"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_10x00_nearbyByFT0C"), vAmpFT0CperColl[k - 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_10x00_thisFT0C"), vAmpFT0CperColl[k], vTracksITS567perColl[k]);
        }
        if (vCollRofId[k - 2] < vCollRofId[k] - 2 && vCollRofId[k - 1] == vCollRofId[k] - 1 && /* next coll is far */ vCollRofId[k + 1] > vCollRofId[k] + 2) // 01x00
        {
          histos.fill(HIST("ROFbyROF/nPV_01x00"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_01x00_nearbyByFT0C"), vAmpFT0CperColl[k - 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_01x00_thisFT0C"), vAmpFT0CperColl[k], vTracksITS567perColl[k]);
        }

        if (vCollRofId[k - 2] < vCollRofId[k] - 2 && vCollRofId[k - 1] == vCollRofId[k] - 1 && /* next coll is far */ vCollRofId[k + 1] > vCollRofId[k] + 1) // 01x0
        {
          // sub-ROFs:
          if (vCollRofSubId[k - 1] == 2 && vCollRofSubId[k] == 0)
            histos.fill(HIST("ROFbyROF/nPV_0_00x_y00_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k - 1] == 1 && vCollRofSubId[k] == 0)
            histos.fill(HIST("ROFbyROF/nPV_0_0x0_y00_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k - 1] == 0 && vCollRofSubId[k] == 0)
            histos.fill(HIST("ROFbyROF/nPV_0_x00_y00_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);

          if (vCollRofSubId[k - 1] == 2 && vCollRofSubId[k] == 1)
            histos.fill(HIST("ROFbyROF/nPV_0_00x_0y0_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k - 1] == 1 && vCollRofSubId[k] == 1)
            histos.fill(HIST("ROFbyROF/nPV_0_0x0_0y0_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k - 1] == 0 && vCollRofSubId[k] == 1)
            histos.fill(HIST("ROFbyROF/nPV_0_x00_0y0_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);

          if (vCollRofSubId[k - 1] == 2 && vCollRofSubId[k] == 2)
            histos.fill(HIST("ROFbyROF/nPV_0_00x_00y_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k - 1] == 1 && vCollRofSubId[k] == 2)
            histos.fill(HIST("ROFbyROF/nPV_0_0x0_00y_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k - 1] == 0 && vCollRofSubId[k] == 2)
            histos.fill(HIST("ROFbyROF/nPV_0_x00_00y_000"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
        }
      }
      // next 1 coll
      if (vCombCond[k] && vCombCond[k + 1]) {
        if (vCollRofId[k - 1] < vCollRofId[k] - 2 /* prev coll is far */ && vCollRofId[k + 1] == vCollRofId[k] + 1 && vCollRofId[k + 2] > vCollRofId[k] + 2) // 00x10
        {
          histos.fill(HIST("ROFbyROF/nPV_00x10"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_00x10_nearbyByFT0C"), vAmpFT0CperColl[k + 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_00x10_thisFT0C"), vAmpFT0CperColl[k], vTracksITS567perColl[k]);
        }

        if (vCollRofId[k - 1] < vCollRofId[k] - 1 /* prev coll is far */ && vCollRofId[k + 1] == vCollRofId[k] + 1 && vCollRofId[k + 2] > vCollRofId[k] + 2) // 0x10
        {
          // sub-ROFs:
          if (vCollRofSubId[k + 1] == 2 && vCollRofSubId[k] == 0)
            histos.fill(HIST("ROFbyROF/nPV_000_y00_00x_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k + 1] == 1 && vCollRofSubId[k] == 0)
            histos.fill(HIST("ROFbyROF/nPV_000_y00_0x0_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k + 1] == 0 && vCollRofSubId[k] == 0)
            histos.fill(HIST("ROFbyROF/nPV_000_y00_x00_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);

          if (vCollRofSubId[k + 1] == 2 && vCollRofSubId[k] == 1)
            histos.fill(HIST("ROFbyROF/nPV_000_0y0_00x_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k + 1] == 1 && vCollRofSubId[k] == 1)
            histos.fill(HIST("ROFbyROF/nPV_000_0y0_0x0_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k + 1] == 0 && vCollRofSubId[k] == 1)
            histos.fill(HIST("ROFbyROF/nPV_000_0y0_x00_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);

          if (vCollRofSubId[k + 1] == 2 && vCollRofSubId[k] == 2)
            histos.fill(HIST("ROFbyROF/nPV_000_00y_00x_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k + 1] == 1 && vCollRofSubId[k] == 2)
            histos.fill(HIST("ROFbyROF/nPV_000_00y_0x0_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          if (vCollRofSubId[k + 1] == 0 && vCollRofSubId[k] == 2)
            histos.fill(HIST("ROFbyROF/nPV_000_00y_x00_0"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
        }

        if (vCollRofId[k - 1] < vCollRofId[k] - 2 /* prev coll is far */ && vCollRofId[k + 1] == vCollRofId[k] + 2 && vCollRofId[k + 2] > vCollRofId[k] + 2) // 00x01
        {
          histos.fill(HIST("ROFbyROF/nPV_00x01"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_00x01_nearbyByFT0C"), vAmpFT0CperColl[k + 1], vTracksITS567perColl[k]);
          histos.fill(HIST("ROFbyROF/nPV_00x01_thisFT0C"), vAmpFT0CperColl[k], vTracksITS567perColl[k]);
        }
      }

      // 2 colls
      if (vCombCond[k] && vCombCond[k - 1] && vCombCond[k - 2]) {
        if (vCollRofId[k - 2] == vCollRofId[k] - 2 && vCollRofId[k - 1] == vCollRofId[k] - 1 && vCollRofId[k + 1] > vCollRofId[k] + 2) // 11x00
          histos.fill(HIST("ROFbyROF/nPV_11x00"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
      }

      if (vCombCond[k] && vCombCond[k - 1] && vCombCond[k + 1]) {
        if (vCollRofId[k - 2] < vCollRofId[k] - 2 && vCollRofId[k - 1] == vCollRofId[k] - 1 && vCollRofId[k + 1] == vCollRofId[k] + 1 && vCollRofId[k + 2] > vCollRofId[k] + 2) // 01x10
          histos.fill(HIST("ROFbyROF/nPV_01x10"), vTracksITS567perColl[k - 1], vTracksITS567perColl[k]);
      }

      if (vCombCond[k] && vCombCond[k + 1] && vCombCond[k + 2]) {
        if (vCollRofId[k - 1] < vCollRofId[k] - 2 && vCollRofId[k + 1] == vCollRofId[k] + 1 && vCollRofId[k + 2] == vCollRofId[k] + 2) // 00x11
          histos.fill(HIST("ROFbyROF/nPV_00x11"), vTracksITS567perColl[k + 1], vTracksITS567perColl[k]);
      }

      // many colls around
      // histos.add("ROFbyROF/nPV_11x11", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    }

    // save indices of collisions in time range for occupancy calculation
    std::vector<std::vector<int>> vCollsInTimeWin;
    std::vector<std::vector<int>> vCollsInSameITSROF;
    std::vector<std::vector<float>> vTimeDeltaForColls; // delta time wrt a given collision
    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];

      // int bcInTF = (foundGlobalBC - bcSOR) % nBCsPerTF;
      // int bcInITSROF = (foundGlobalBC + 3564 - rofOffset) % rofLength;
      int64_t tfId = (foundGlobalBC - bcSOR) / nBCsPerTF;
      int64_t rofId = (foundGlobalBC + 3564 - rofOffset) / rofLength;
      // int rofIdInTF = (bcInTF - rofOffset) / rofLength;

      // ### in-ROF occupancy
      std::vector<int> vAssocToSameROF;
      // find all collisions in the same ROF before a given collision
      int32_t minColIndex = colIndex - 1;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];
        // check if this is still the same TF
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        // int thisRofIdInTF = (thisBC - rofOffset) / rofLength;
        int64_t thisRofId = (thisBC + 3564 - rofOffset) / rofLength;

        // check if we are within the same ROF
        if (thisRofId != rofId)
          break;
        vAssocToSameROF.push_back(minColIndex);
        minColIndex--;
      }
      // find all collisions in the same ROF after the current one
      int32_t maxColIndex = colIndex + 1;
      while (maxColIndex < cols.size()) {
        int64_t thisBC = vFoundGlobalBC[maxColIndex];
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        // int thisRofIdInTF = (thisBC - rofOffset) / rofLength;
        int64_t thisRofId = (thisBC + 3564 - rofOffset) / rofLength;
        if (thisRofId != rofId)
          break;
        vAssocToSameROF.push_back(maxColIndex);
        maxColIndex++;
      }
      vCollsInSameITSROF.push_back(vAssocToSameROF);

      // ### occupancy in time windows
      std::vector<int> vAssocToThisCol;
      std::vector<float> vCollsTimeDeltaWrtGivenColl;
      // protection against the TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) {
        vCollsInTimeWin.push_back(vAssocToThisCol);
        vTimeDeltaForColls.push_back(vCollsTimeDeltaWrtGivenColl);
        continue;
      }
      // find all collisions in time window before the current one
      minColIndex = colIndex - 1;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];
        // check if this is still the same TF
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        float dt = (thisBC - foundGlobalBC) * bcNS; // ns
        // check if we are within the chosen time range
        if (dt < timeWinOccupancyCalcMinNS)
          break;
        vAssocToThisCol.push_back(minColIndex);
        vCollsTimeDeltaWrtGivenColl.push_back(dt);
        minColIndex--;
      }
      // find all collisions in time window after the current one
      maxColIndex = colIndex + 1;
      while (maxColIndex < cols.size()) {
        int64_t thisBC = vFoundGlobalBC[maxColIndex];
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != tfId)
          break;
        float dt = (thisBC - foundGlobalBC) * bcNS; // ns
        if (dt > timeWinOccupancyCalcMaxNS)
          break;
        vAssocToThisCol.push_back(maxColIndex);
        vCollsTimeDeltaWrtGivenColl.push_back(dt);
        maxColIndex++;
      }
      vCollsInTimeWin.push_back(vAssocToThisCol);
      vTimeDeltaForColls.push_back(vCollsTimeDeltaWrtGivenColl);
    }

    // perform the occupancy calculation per ITS ROF and also in the pre-defined time window
    std::vector<int> vNumTracksITS567inFullTimeWin(cols.size(), 0); // counter of tracks in full time window for occupancy studies (excluding given event)
    std::vector<float> vSumAmpFT0CinFullTimeWin(cols.size(), 0);    // sum of FT0C of tracks in full time window for occupancy studies (excluding given event)
    std::vector<int> vNumTracksITS567inROF(cols.size(), 0);         // counter of tracks in given ROF (excluding given event)
    std::vector<int> vSumAmpFT0CinROF(cols.size(), 0);              // counter of tracks in given ROF (excluding given event)
    std::vector<int> vNumCollinROF(cols.size(), 0);                 // counter of tracks in given ROF (excluding given event)
    std::vector<int> vNumCollinROFinVz10(cols.size(), 0);           // counter of tracks in given ROF (excluding given event)
    std::vector<int> vInROFcollIndex(cols.size(), 0);               // counter of tracks in given ROF (excluding given event)
    std::vector<int> vROFidThisColl(cols.size(), 0);                // counter of tracks in given ROF (excluding given event)

    std::vector<bool> vNoCollInTimeRangeStrict(cols.size(), 0);      // no collisions in a specified time range
    std::vector<bool> vNoCollInTimeRangeNarrow(cols.size(), 0);      // no collisions in a specified time range (narrow)
    std::vector<bool> vNoHighMultCollInTimeRange(cols.size(), 0);    // no high-mult collisions in a specified time range
    std::vector<bool> vNoCollInVzDependentTimeRange(cols.size(), 0); // no collisions in a vZ-dependent time range

    std::vector<bool> vNoCollInSameRofStrict(cols.size(), 0);      // to veto events with other collisions in the same ITS ROF
    std::vector<bool> vNoCollInSameRofStandard(cols.size(), 0);    // to veto events with other collisions in the same ITS ROF, with per-collision multiplicity above threshold
    std::vector<bool> vNoCollInSameRofWithCloseVz(cols.size(), 0); // to veto events with nearby collisions with close vZ
    std::vector<std::vector<bool>> vArrNoCollInSameRofWithCloseVz; //(cols.size(), 0); // to veto events with nearby collisions with close vZ

    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      float vZ = col.posZ();

      // QA:
      if (vAmpFT0CperColl[colIndex] > 5000) {
        histos.fill(HIST("hThisEvITStr_vs_vZ"), vZ, vTracksITS567perColl[colIndex]);
        histos.fill(HIST("hThisEvITSTPCtr_vs_vZ"), vZ, vTracksITSTPCperColl[colIndex]);
      }

      // ### in-ROF occupancy
      // int64_t rofId = (vFoundGlobalBC[colIndex] + 3564 - rofOffset) / rofLength;
      // int bcInTF = (vFoundGlobalBC[colIndex] - bcSOR) % nBCsPerTF;
      int bcInITSROF = (vFoundGlobalBC[colIndex] + 3564 - rofOffset) % rofLength;
      int rofIdInTF = (vFoundGlobalBC[colIndex] + 3564 - rofOffset) / rofLength;
      // auto bc = bcs.iteratorAt(vFoundBCindex[colIndex]);
      // LOGP(info, "#### starting new coll: bc={} bcInTF={} bcInITSROF={} rofId={};  noROFborder={};   rofOffset={} rofLength={}", vFoundGlobalBC[colIndex], bcInTF, bcInITSROF, rofId, bc.selection_bit(kNoITSROFrameBorder), rofOffset, rofLength);
      // LOGP(info, "#### starting new coll: bcInTF={} bcInITSROF={} rofIdInTF={};  noROFborder={},  vZ={} mult={};   rofOffset={} rofLength={}", bcInTF, bcInITSROF, rofIdInTF, bc.selection_bit(kNoITSROFrameBorder), vZ, vTracksITS567perColl[colIndex], rofOffset, rofLength);

      std::vector<int> vAssocToSameROF = vCollsInSameITSROF[colIndex];
      int nITS567tracksForRofVetoStrict = 0;  // to veto events with other collisions in the same ITS ROF
      float nSumAmplFT0CforRofVetoStrict = 0; // to veto events with other collisions in the same ITS ROF
      // int nITS567tracksForRofVetoStandard = 0;           // to veto events with other collisions in the same ITS ROF, with per-collision multiplicity above threshold
      int nCollsInRofWithFT0CAboveVetoStandard = 0;      // to veto events with other collisions in the same ITS ROF, with per-collision multiplicity above threshold
      int nITS567tracksForRofVetoOnCloseVz = 0;          // to veto events with nearby collisions with close vZ
      int nArrITS567tracksForRofVetoOnCloseVz[200] = {}; // to veto events with nearby collisions with close vZ
      vNumCollinROF[colIndex] = 1;
      vInROFcollIndex[colIndex] = 0;
      vROFidThisColl[colIndex] = rofIdInTF;

      if (std::fabs(vZ) < 10)
        vNumCollinROFinVz10[colIndex] = 1;
      for (uint32_t iCol = 0; iCol < vAssocToSameROF.size(); iCol++) {
        int thisColIndex = vAssocToSameROF[iCol];
        // int64_t thisRofId = (vFoundGlobalBC[thisColIndex] + 3564 - rofOffset) / rofLength;
        // int thisBcInTF = (vFoundGlobalBC[thisColIndex] - bcSOR) % nBCsPerTF;
        int thisBcInITSROF = (vFoundGlobalBC[thisColIndex] + 3564 - rofOffset) % rofLength;
        // int thisRofIdInTF = (vFoundGlobalBC[thisColIndex] + 3564 - rofOffset) / rofLength;
        // auto bcAssoc = bcs.iteratorAt(vFoundBCindex[thisColIndex]);
        // LOGP(info, ">> assoc: bc={} bcInTF={} bcInITSROF={} rofId={} noROFborder={}", vFoundGlobalBC[thisColIndex], thisBcInTF, thisBcInITSROF, thisRofId, bcAssoc.selection_bit(kNoITSROFrameBorder));
        // LOGP(info, ">> assoc: bcInTF={} bcInITSROF={} rofIdInTF={} noROFborder={} vZ={} mult={}", thisBcInTF, thisBcInITSROF, thisRofIdInTF, bcAssoc.selection_bit(kNoITSROFrameBorder), vCollVz[thisColIndex], vTracksITS567perColl[thisColIndex]);

        // if (std::fabs(vTracksITS567perColl[thisColIndex]) > confNtracksCutVetoOnCollInROF)
        nITS567tracksForRofVetoStrict += vTracksITS567perColl[thisColIndex];
        nSumAmplFT0CforRofVetoStrict += vAmpFT0CperColl[thisColIndex];
        vNumCollinROF[colIndex]++;
        if (std::fabs(vCollVz[thisColIndex]) < 10)
          vNumCollinROFinVz10[colIndex]++;
        vInROFcollIndex[colIndex] = thisBcInITSROF > bcInITSROF ? 0 : 1; // if colIndex is for the first coll in ROF => inROFindex=0, otherwise =1

        // if (vTracksITS567perColl[thisColIndex] > confNtracksCutVetoOnCollInROF)
        // nITS567tracksForRofVetoStandard += vTracksITS567perColl[thisColIndex];

        if (vAmpFT0CperColl[thisColIndex] > confFT0CamplCutVetoOnCollInROF)
          nCollsInRofWithFT0CAboveVetoStandard++;

        if (std::fabs(vCollVz[thisColIndex] - vZ) < confEpsilonVzDiffVetoInROF)
          nITS567tracksForRofVetoOnCloseVz += vTracksITS567perColl[thisColIndex];
        for (int i = 0; i < 200; i++) {
          // if (std::fabs(vCollVz[thisColIndex] - vZ) < 0.05 * i && vTracksITS567perColl[thisColIndex] > 50)
          // if (vTracksITS567perColl[colIndex]>100 && vTracksITS567perColl[colIndex]<1000 &&
          if (vAmpFT0CperColl[colIndex] > 4000 && vAmpFT0CperColl[colIndex] < 15000 &&
              (vCollVz[thisColIndex] - vZ) > 0.05 * i && std::fabs(vCollVz[thisColIndex] - vZ) < (0.1 + 0.05) * i && vTracksITS567perColl[thisColIndex] > 20) // 0.05 * (i + 1))
            // std::fabs(vCollVz[thisColIndex] - vZ) < 0.05 * i && vTracksITS567perColl[thisColIndex] > 30) // 0.05 * (i + 1))
            nArrITS567tracksForRofVetoOnCloseVz[i]++;
        }

        if (std::fabs(vZ) < 10) {
          histos.fill(HIST("hDeltaVz"), vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] >= 100 && vTracksITS567perColl[thisColIndex] < 100)
            histos.fill(HIST("hDeltaVzGivenCollAbove100NearbyBelow100"), vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] <= 100 && vTracksITS567perColl[thisColIndex] > 100)
            histos.fill(HIST("hDeltaVzGivenCollBelow100NearbyAbove100"), vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] > 20 && vTracksITS567perColl[thisColIndex] > 20)
            histos.fill(HIST("hDeltaVzAfterCuts"), vCollVz[thisColIndex] - vZ);
          if (col.sel8()) // bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder) && bc.selection_bit(kNoITSROFrameBorder))
            histos.fill(HIST("hDeltaVzAfterTFandROFborderCuts"), vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] > 20 && vTracksITS567perColl[thisColIndex] > 20 && col.sel8()) // bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder) && bc.selection_bit(kNoITSROFrameBorder))
            histos.fill(HIST("hDeltaVzAfterAllCuts"), vCollVz[thisColIndex] - vZ);
        }
      }
      vNumTracksITS567inROF[colIndex] = nITS567tracksForRofVetoStrict; // occupancy in ROF (excluding a given collision)
      vSumAmpFT0CinROF[colIndex] = nSumAmplFT0CforRofVetoStrict;       // occupancy in ROF (excluding a given collision)

      // in-ROF occupancy flags
      vNoCollInSameRofStrict[colIndex] = (nITS567tracksForRofVetoStrict == 0);
      // vNoCollInSameRofStandard[colIndex] = (nITS567tracksForRofVetoStandard == 0);
      vNoCollInSameRofStandard[colIndex] = (nCollsInRofWithFT0CAboveVetoStandard == 0);
      vNoCollInSameRofWithCloseVz[colIndex] = (nITS567tracksForRofVetoOnCloseVz == 0);

      std::vector<bool> vVzCutThisColl;

      // ### occupancy in time windows
      // protection against TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) { // occupancy in undefined (too close to TF borders)
        vNumTracksITS567inFullTimeWin[colIndex] = -1;
        vSumAmpFT0CinFullTimeWin[colIndex] = -1;
        // vNumTracksITS567inROF[colIndex] = -1;
        vArrNoCollInSameRofWithCloseVz.push_back(vVzCutThisColl);
        continue;
      }
      std::vector<int> vAssocToThisCol = vCollsInTimeWin[colIndex];
      std::vector<float> vCollsTimeDeltaWrtGivenColl = vTimeDeltaForColls[colIndex];
      int nITS567tracksInFullTimeWindow = 0;
      int sumAmpFT0CInFullTimeWindow = 0;
      int nITS567tracksForVetoNarrow = 0;      // to veto events with nearby collisions (narrower range)
      int nITS567tracksForVetoStrict = 0;      // to veto events with nearby collisions
      int nITS567tracksForVetoStandard = 0;    // to veto events with per-collision multiplicity above threshold
      int nITS567tracksForVetoVzDependent = 0; // to veto events with nearby collisions, vZ-dependent time cut
      for (uint32_t iCol = 0; iCol < vAssocToThisCol.size(); iCol++) {
        int thisColIndex = vAssocToThisCol[iCol];
        float dt = vCollsTimeDeltaWrtGivenColl[iCol] / 1e3; // ns -> us
        histos.fill(HIST("hDeltaTime"), dt);
        if (vTracksITS567perColl[colIndex] > 50 && vTracksITS567perColl[thisColIndex] > 50)
          histos.fill(HIST("hDeltaTimeAboveNtracksCut"), dt);

        if (std::fabs(vCollVz[colIndex]) < 10 && std::fabs(vCollVz[thisColIndex]) < 10)
          histos.fill(HIST("hDeltaTime_vZ10cm"), dt);

        if (vIsSel8[colIndex] && vIsSel8[thisColIndex])
          histos.fill(HIST("hDeltaTime_sel8"), dt);

        if (std::fabs(vCollVz[colIndex]) < 10 && vIsSel8[colIndex] && std::fabs(vCollVz[thisColIndex]) < 10 && vIsSel8[thisColIndex]) {
          histos.fill(HIST("hDeltaTime_sel8_vZ10cm"), dt);
          if (vTracksITS567perColl[colIndex] > 50 && vTracksITS567perColl[thisColIndex] > 50)
            histos.fill(HIST("hDeltaTimeAboveNtracksCut_sel8_vZ10cm"), dt);
        }

        float wOccup = 1.;
        if (confUseWeightsForOccupancyVariable) {
          // weighted occupancy
          wOccup = 0;
          if (dt >= -40 && dt < -5) // collisions in the past
            wOccup = 1. / 1225 * (dt + 40) * (dt + 40);
          else if (dt >= -5 && dt < 15) // collisions near a given one
            wOccup = 1;
          // else if (dt >= 15 && dt < 100) // collisions from the future
          //   wOccup = -1. / 85 * dt + 20. / 17;
          else if (dt >= 15 && dt < 40) // collisions from the future
            wOccup = -0.4 / 25 * dt + 1.24;
          else if (dt >= 40 && dt < 100) // collisions from the distant future
            wOccup = -0.4 / 60 * dt + 0.6 + 0.8 / 3;
          if (wOccup > 0)
            histos.fill(HIST("hOccupancyWeights"), dt, wOccup);
        }
        nITS567tracksInFullTimeWindow += wOccup * vTracksITS567perColl[thisColIndex];
        sumAmpFT0CInFullTimeWindow += wOccup * vAmpFT0CperColl[thisColIndex];

        // counting tracks from other collisions in fixed time windows
        if (std::fabs(dt) < confTimeRangeVetoOnCollNarrow)
          nITS567tracksForVetoNarrow += vTracksITS567perColl[thisColIndex];
        if (std::fabs(dt) < confTimeRangeVetoOnCollStandard)
          nITS567tracksForVetoStrict += vTracksITS567perColl[thisColIndex];

        // if (std::fabs(dt) < confTimeRangeVetoOnCollStandard + 0.5) { // add 0.5 us safety margin
        // standard cut on other collisions vs delta-times
        const float driftV = 2.5;  // drift velocity in cm/us, TPC drift_length / drift_time = 250 cm / 100 us
        if (std::fabs(dt) < 2.0) { // us, complete veto on other collisions
          nITS567tracksForVetoStandard += vTracksITS567perColl[thisColIndex];
        } else if (dt > -4.0 && dt <= -2.0) { // us, strict veto to suppress fake ITS-TPC matches more
          if (vTracksITS567perColl[thisColIndex] > confNtracksCutVetoOnCollInTimeRange / 5)
            nITS567tracksForVetoStandard += vTracksITS567perColl[thisColIndex];
        } else if (std::fabs(dt) < 8 + std::fabs(vZ) / driftV) { // loose veto, 8 us corresponds to maximum possible |vZ|, which is ~20 cm
          // counting number of other collisions with mult above threshold
          if (vTracksITS567perColl[thisColIndex] > confNtracksCutVetoOnCollInTimeRange)
            nITS567tracksForVetoStandard += vTracksITS567perColl[thisColIndex];
        }
        // vZ-dependent time cut to avoid collinear tracks from other collisions (experimental)
        if (std::fabs(dt) < 8 + std::fabs(vZ) / driftV) {
          if (dt < 0) {
            // check distance between given vZ and (moving in two directions) vZ of drifting tracks from past collisions
            if ((std::fabs(vCollVz[thisColIndex] - std::fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC) ||
                (std::fabs(vCollVz[thisColIndex] + std::fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC))
              nITS567tracksForVetoVzDependent += vTracksITS567perColl[thisColIndex];

            // FOR QA:
            if (std::fabs(vCollVz[thisColIndex] - std::fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.fill(HIST("hDeltaVzVsDeltaTime1"), vCollVz[thisColIndex] - std::fabs(dt) * driftV - vZ, dt);
            if (std::fabs(vCollVz[thisColIndex] + std::fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.fill(HIST("hDeltaVzVsDeltaTime2"), vCollVz[thisColIndex] + std::fabs(dt) * driftV - vZ, dt);

          } else { // dt>0
            // check distance between drifted vZ of given collision (in two directions) and vZ of future collisions
            if ((std::fabs(vZ - dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC) ||
                (std::fabs(vZ + dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC))
              nITS567tracksForVetoVzDependent += vTracksITS567perColl[thisColIndex];

            // FOR QA:
            if (std::fabs(vZ - dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.fill(HIST("hDeltaVzVsDeltaTime3"), vZ - dt * driftV - vCollVz[thisColIndex], dt);
            if (std::fabs(vZ + dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.fill(HIST("hDeltaVzVsDeltaTime4"), vZ + dt * driftV - vCollVz[thisColIndex], dt);
          }
        }
      }
      vNumTracksITS567inFullTimeWin[colIndex] = nITS567tracksInFullTimeWindow; // occupancy by ITS tracks (without a current collision)
      vSumAmpFT0CinFullTimeWin[colIndex] = sumAmpFT0CInFullTimeWindow;         // occupancy by FT0C (without a current collision)
      // occupancy flags based on nearby collisions
      vNoCollInTimeRangeNarrow[colIndex] = (nITS567tracksForVetoNarrow == 0);
      vNoCollInTimeRangeStrict[colIndex] = (nITS567tracksForVetoStrict == 0);
      vNoHighMultCollInTimeRange[colIndex] = (nITS567tracksForVetoStandard == 0);
      vNoCollInVzDependentTimeRange[colIndex] = (nITS567tracksForVetoVzDependent == 0);

      for (int i = 0; i < 200; i++)
        vVzCutThisColl.push_back(nArrITS567tracksForRofVetoOnCloseVz[i] == 0);
      vArrNoCollInSameRofWithCloseVz.push_back(vVzCutThisColl);
    }

    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      bool sel8 = col.sel8(); // bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder) && bc.selection_bit(kNoITSROFrameBorder);

      float vZ = vCollVz[colIndex];
      int occTracks = col.trackOccupancyInTimeRange();
      float occFT0C = col.ft0cOccupancyInTimeRange();

      if (sel8) {
        // histos.fill(HIST("hColCounterAcc"),Form("%d", bc.runNumber()), 1);
        if (vInROFcollIndex[colIndex] == 0)
          histos.fill(HIST("hNcollPerROF"), vNumCollinROF[colIndex]);

        histos.fill(HIST("hVz"), vZ);

        int nPV = vTracksITS567perColl[colIndex];
        float ft0C = vAmpFT0CperColl[colIndex];

        // ROF-by-ROF
        if (std::fabs(vZ) < 8) {
          histos.fill(HIST("ROFbyROF/nPV_vs_ROFid"), nPV, vCollRofIdPerOrbit[colIndex]);
          histos.fill(HIST("ROFbyROF/nPV_vs_subROFid"), nPV, vCollRofSubIdPerOrbit[colIndex]);

          histos.fill(HIST("ROFbyROF/FT0C_vs_ROFid"), ft0C, vCollRofIdPerOrbit[colIndex]);
          histos.fill(HIST("ROFbyROF/FT0C_vs_subROFid"), ft0C, vCollRofSubIdPerOrbit[colIndex]);
        }
        // vs occupancy
        if (occTracks >= 0 && std::fabs(vZ) < 8) {
          histos.fill(HIST("nPV_vs_occupancyByTracks/sel8"), nPV, occTracks);
          if (col.selection_bit(kNoCollInTimeRangeNarrow))
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInTimeRangeNarrow"), nPV, occTracks);
          if (col.selection_bit(kNoCollInTimeRangeStrict))
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInTimeRangeStrict"), nPV, occTracks);
          if (col.selection_bit(kNoCollInTimeRangeStandard))
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInTimeRangeStandard"), nPV, occTracks);
          if (col.selection_bit(kNoCollInRofStrict))
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInRofStrict"), nPV, occTracks);
          if (col.selection_bit(kNoCollInRofStandard))
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInRofStandard"), nPV, occTracks);
          if (col.selection_bit(kNoCollInTimeRangeStandard) && col.selection_bit(kNoCollInRofStandard))
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInTimeAndRofStandard"), nPV, occTracks);
          if (col.selection_bit(kNoCollInTimeRangeStrict) && col.selection_bit(kNoCollInRofStrict))
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInTimeAndRofStrict"), nPV, occTracks);
          if (col.selection_bit(kNoCollInTimeRangeStrict) && col.selection_bit(kNoCollInRofStrict) && std::fabs(vZ) < 5)
            histos.fill(HIST("nPV_vs_occupancyByTracks/NoCollInTimeAndRofStrict_vZ_5cm"), nPV, occTracks);
          if (col.selection_bit(kNoHighMultCollInPrevRof))
            histos.fill(HIST("nPV_vs_occupancyByTracks/kNoHighMultCollInPrevRof"), nPV, occTracks);
          if (col.selection_bit(kNoHighMultCollInPrevRof) && col.selection_bit(kNoCollInRofStrict))
            histos.fill(HIST("nPV_vs_occupancyByTracks/kNoHighMultCollInPrevRofAndRofStrict"), nPV, occTracks);
        }
        if (occFT0C >= 0 && std::fabs(vZ) < 8) {
          histos.fill(HIST("nPV_vs_occupancyByFT0C/sel8"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInTimeRangeNarrow))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInTimeRangeNarrow"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInTimeRangeStrict))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInTimeRangeStrict"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInTimeRangeStandard))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInTimeRangeStandard"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInRofStrict))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInRofStrict"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInRofStandard))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInRofStandard"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInTimeRangeStandard) && col.selection_bit(kNoCollInRofStandard))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInTimeAndRofStandard"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInTimeRangeStrict) && col.selection_bit(kNoCollInRofStrict))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInTimeAndRofStrict"), nPV, occFT0C);
          if (col.selection_bit(kNoCollInTimeRangeStrict) && col.selection_bit(kNoCollInRofStrict) && std::fabs(vZ) < 5)
            histos.fill(HIST("nPV_vs_occupancyByFT0C/NoCollInTimeAndRofStrict_vZ_5cm"), nPV, occFT0C);
          if (col.selection_bit(kNoHighMultCollInPrevRof))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/kNoHighMultCollInPrevRof"), nPV, occFT0C);
          if (col.selection_bit(kNoHighMultCollInPrevRof) && col.selection_bit(kNoCollInRofStrict))
            histos.fill(HIST("nPV_vs_occupancyByFT0C/kNoHighMultCollInPrevRofAndRofStrict"), nPV, occFT0C);
        }
      }

      if (occTracks >= 0)
        histos.fill(HIST("hOccupancyByTracks_CROSSCHECK"), occTracks);
      if (occFT0C >= 0)
        histos.fill(HIST("hOccupancyByFT0C_CROSSCHECK"), occFT0C);

      if (vNumTracksITS567inFullTimeWin[colIndex] >= 0) {
        histos.fill(HIST("hOccupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex]);
        histos.fill(HIST("hOccupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex]);

        histos.fill(HIST("hOccupancyByTrInROF"), vNumTracksITS567inROF[colIndex]);
        histos.fill(HIST("hOccupancyByFT0C_vs_ByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
        histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/all"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
        histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/all"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

        if (sel8 && std::fabs(col.posZ()) < 8) {
          histos.fill(HIST("hOccupancyByFT0C_vs_ByTracks_vZ_TF_ROF_border_cuts"), vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);

          // if (vAmpFT0CperColl[colIndex] > 5000 && vAmpFT0CperColl[colIndex] < 10000) {
          // if (vAmpFT0CperColl[colIndex] > 500) {
          if (vAmpFT0CperColl[colIndex] > 0) { // 50) {

            histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
            histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
            histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);

            histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
            histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyInROF"), vNumTracksITS567inROF[colIndex], vTracksITSTPCperColl[colIndex]);

            if (vNumTracksITS567inROF[colIndex] > 0) {
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_HasNeighbours"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours"), vSumAmpFT0CinROF[colIndex], vTracksITS567perColl[colIndex]);

              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours"), vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
            }

            // two collisions in one ROF (both with |vZ|<10 cm)
            if (vNumCollinROF[colIndex] == 2 && vNumCollinROFinVz10[colIndex] == 2 && vInROFcollIndex[colIndex] == 1) {
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvFT0C_vs_occupancyFT0CInROF_2coll"), vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);

              if (vNumTracksITS567inROF[colIndex] > vTracksITS567perColl[colIndex]) {
                histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                if (vNumTracksITS567inROF[colIndex] > 0)
                  histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"), vNumTracksITS567inROF[colIndex], 1.0 * vTracksITS567perColl[colIndex] / vNumTracksITS567inROF[colIndex]);

              } else {
                histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"), vTracksITS567perColl[colIndex], vNumTracksITS567inROF[colIndex]);

                histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_allOccup_2coll_inROF"), vTracksITS567perColl[colIndex]);
                if (occTracks < 500)
                  histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_lowOccup_2coll_inROF"), vTracksITS567perColl[colIndex]);
                else if (occTracks > 1000)
                  histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_highOccup_2coll_inROF"), vTracksITS567perColl[colIndex]);

                if (vTracksITS567perColl[colIndex] > 0)
                  histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"), vTracksITS567perColl[colIndex], 1.0 * vNumTracksITS567inROF[colIndex] / vTracksITS567perColl[colIndex]);
              }
            }

            // 3 or 4 collisions in one ROF
            if (vNumCollinROF[colIndex] == 2 && vInROFcollIndex[colIndex] == 1)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_noVzCutOnOtherVertices"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
            if (vNumCollinROF[colIndex] == 3 && vInROFcollIndex[colIndex] == 1)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_3coll_noVzCutOnOtherVertices"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
            if (vNumCollinROF[colIndex] == 4 && vInROFcollIndex[colIndex] == 1)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_4coll_noVzCutOnOtherVertices"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);

            // now 1D histograms vs nCollInROF
            if (vNumCollinROF[colIndex] == 1)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_1coll_in_ROF"), vTracksITS567perColl[colIndex]);
            if (vNumCollinROF[colIndex] == 2)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_2coll_in_ROF"), vTracksITS567perColl[colIndex]);
            if (vNumCollinROF[colIndex] == 3)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_3coll_in_ROF"), vTracksITS567perColl[colIndex]);
            if (vNumCollinROF[colIndex] == 4)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_4coll_in_ROF"), vTracksITS567perColl[colIndex]);
            if (vNumCollinROF[colIndex] >= 5)
              histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_5collOrMore_in_ROF"), vTracksITS567perColl[colIndex]);

            // compare with previous ROF
            if (colIndex - 1 >= 0) {
              if (vNumCollinROF[colIndex] == 1 && vNumCollinROFinVz10[colIndex] == 1 && vInROFcollIndex[colIndex] == 0 && vNumCollinROF[colIndex - 1] == 1 && vNumCollinROFinVz10[colIndex - 1] == 1 && vInROFcollIndex[colIndex - 1] == 0) {
                histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInAnotherEarlierROF_1collPerROF"), vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                if (vROFidThisColl[colIndex] == vROFidThisColl[colIndex - 1] + 1) // one ROF right after the previous
                {
                  histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF"), vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                  if (vTracksITS567perColl[colIndex - 1] > vTracksITS567perColl[colIndex]) {
                    histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"), vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                    if (vTracksITS567perColl[colIndex - 1] > 0)
                      histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"), vTracksITS567perColl[colIndex - 1], 1.0 * vTracksITS567perColl[colIndex] / vTracksITS567perColl[colIndex - 1]);
                  } else {
                    histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"), vTracksITS567perColl[colIndex], vTracksITS567perColl[colIndex - 1]);

                    histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_allOccup_1collPerROF"), vTracksITS567perColl[colIndex]);
                    if (occTracks < 500)
                      histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_lowOccup_1collPerROF"), vTracksITS567perColl[colIndex]);
                    else if (occTracks > 1000)
                      histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_highOccup_1collPerROF"), vTracksITS567perColl[colIndex]);

                    if (vTracksITS567perColl[colIndex] > 0)
                      histos.fill(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"), vTracksITS567perColl[colIndex], 1.0 * vTracksITS567perColl[colIndex - 1] / vTracksITS567perColl[colIndex]);
                  }
                }
              }
            }
          }

          if (vNoCollInTimeRangeNarrow[colIndex])
            histos.fill(HIST("hOccupancyByFT0C_vs_ByTracks_afterNarrowDeltaTimeCut"), vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
          if (vNoCollInTimeRangeStrict[colIndex])
            histos.fill(HIST("hOccupancyByFT0C_vs_ByTracks_afterStrictDeltaTimeCut"), vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
          if (vNoHighMultCollInTimeRange[colIndex])
            histos.fill(HIST("hOccupancyByFT0C_vs_ByTracks_afterStandardDeltaTimeCut"), vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
          if (vNoCollInVzDependentTimeRange[colIndex])
            histos.fill(HIST("hOccupancyByFT0C_vs_ByTracks_afterVzDependentDeltaTimeCut"), vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);

          // same-event 2D correlations:
          histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/vZ_TF_ROF_border_cuts"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (vNoCollInTimeRangeNarrow[colIndex])
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/afterNarrowDeltaTimeCut"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInTimeRangeStrict[colIndex])
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/afterStrictDeltaTimeCut"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoHighMultCollInTimeRange[colIndex])
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/afterStandardDeltaTimeCut"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInVzDependentTimeRange[colIndex])
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/afterVzDependentDeltaTimeCut"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (vNoCollInSameRofStrict[colIndex])
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/kNoCollInRofStrict"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInSameRofStandard[colIndex])
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/kNoCollInRofStandard"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInSameRofWithCloseVz[colIndex])
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/kNoCollInRofWithCloseVz"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          // CROSS CHECK WITH SEL BITS:
          if (col.selection_bit(kNoCollInTimeRangeNarrow))
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_afterNarrowDeltaTimeCut"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (col.selection_bit(kNoCollInTimeRangeStrict))
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_afterStrictDeltaTimeCut"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (col.selection_bit(kNoCollInTimeRangeStandard))
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_afterStandardDeltaTimeCut"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (col.selection_bit(kNoCollInRofStrict))
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_kNoCollInRofStrict"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (col.selection_bit(kNoCollInRofStandard))
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/CROSSCHECK_kNoCollInRofStandard"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (vNumTracksITS567inFullTimeWin[colIndex] < 2000) {
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/occupBelow2000"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/hThisEvITSTPCTr_vs_ThisEvFT0C_occupBelow2000"), vAmpFT0CperColl[colIndex], vTracksITSTPCperColl[colIndex]);
          }

          if (vNoCollInTimeRangeNarrow[colIndex] && vNoHighMultCollInTimeRange[colIndex] && vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
            histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/NarrowDeltaCut_StdTimeAndRofCuts"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

            if (vNumTracksITS567inFullTimeWin[colIndex] < 2000) {
              histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000"), vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("hThisEvITSTr_vs_ThisEvFT0C/hThisEvITSTPCTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000"), vAmpFT0CperColl[colIndex], vTracksITSTPCperColl[colIndex]);
            }
          }

          // now ITSTPC vs ITS tr (this event)
          histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/vZ_TF_ROF_border_cuts"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNoCollInTimeRangeNarrow[colIndex])
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/afterNarrowDeltaTimeCut"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInTimeRangeStrict[colIndex])
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/afterStrictDeltaTimeCut"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoHighMultCollInTimeRange[colIndex])
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/afterStandardDeltaTimeCut"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInVzDependentTimeRange[colIndex])
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/afterVzDependentDeltaTimeCut"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNoCollInSameRofStrict[colIndex])
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/kNoCollInRofStrict"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInSameRofStandard[colIndex])
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/kNoCollInRofStandard"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInSameRofWithCloseVz[colIndex])
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/kNoCollInRofWithCloseVz"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNumTracksITS567inFullTimeWin[colIndex] < 2000)
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/occupBelow2000"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNoCollInTimeRangeNarrow[colIndex] && vNoHighMultCollInTimeRange[colIndex] && vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
            histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/NarrowDeltaCut_StdTimeAndRofCuts"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
            if (vNumTracksITS567inFullTimeWin[colIndex] < 2000)
              histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/occupBelow2000_NarrowDeltaCut_StdTimeAndRofCuts"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          }

          if (vNoCollInTimeRangeStrict[colIndex] && vNoCollInSameRofStrict[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
            if (vNumTracksITS567inFullTimeWin[colIndex] < 2000)
              histos.fill(HIST("hThisEvITSTPCTr_vs_ThisEvITStr/occupBelow2000_StrictDeltaTimeCutAndRofCuts"), vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          }

          // vZ bins to tune vZthresh cut
          if (vNoCollInTimeRangeNarrow[colIndex]) {

            for (int i = 0; i < 200; i++) {
              if (std::fabs(col.posZ()) < 8 && !vArrNoCollInSameRofWithCloseVz[colIndex][i]) {
                histos.fill(HIST("hThisEvITStr_vs_vZcut"), 0.025 + 0.05 * i, vTracksITS567perColl[colIndex]);
                histos.fill(HIST("hThisEvITSTPCtr_vs_vZcut"), 0.025 + 0.05 * i, vTracksITSTPCperColl[colIndex]);
              }
            }
          }

          // ### this event vs Occupancy 2D histos
          // if (vAmpFT0CperColl[colIndex] > 5000 && vAmpFT0CperColl[colIndex] < 10000) {
          // if (vAmpFT0CperColl[colIndex] > 500) {
          if (vAmpFT0CperColl[colIndex] > 0) { // 100) {

            if (vNoCollInTimeRangeNarrow[colIndex]) {
              histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);

              histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyInROF"), vNumTracksITS567inROF[colIndex], vTracksITSTPCperColl[colIndex]);

              if (vNumTracksITS567inROF[colIndex] > 0) {
                histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_HasNeighbours"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours"), vSumAmpFT0CinROF[colIndex], vTracksITS567perColl[colIndex]);

                histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours"), vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
              }

              if (vNumCollinROF[colIndex] == 2 && vNumCollinROFinVz10[colIndex] == 2 && vInROFcollIndex[colIndex] == 1) {
                histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvFT0C_vs_occupancyFT0CInROF_2coll"), vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);

                if (vNumTracksITS567inROF[colIndex] > vTracksITS567perColl[colIndex]) {
                  histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                  if (vNumTracksITS567inROF[colIndex] > 0)
                    histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"), vNumTracksITS567inROF[colIndex], 1.0 * vTracksITS567perColl[colIndex] / vNumTracksITS567inROF[colIndex]);
                } else {
                  histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"), vTracksITS567perColl[colIndex], vNumTracksITS567inROF[colIndex]);
                  if (vTracksITS567perColl[colIndex] > 0)
                    histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"), vTracksITS567perColl[colIndex], 1.0 * vNumTracksITS567inROF[colIndex] / vTracksITS567perColl[colIndex]);
                }

                // the sum of v1 and v2:
                if (vSumAmpFT0CinROF[colIndex] > 4000 && vAmpFT0CperColl[colIndex] > 4000)
                  histos.fill(HIST("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_inROF"), vTracksITS567perColl[colIndex] + vNumTracksITS567inROF[colIndex]);
              }
              // compare with previous ROF
              if (colIndex - 1 >= 0) {
                if (vNumCollinROF[colIndex] == 1 && vNumCollinROFinVz10[colIndex] == 1 && vInROFcollIndex[colIndex] == 0 && vNumCollinROF[colIndex - 1] == 1 && vNumCollinROFinVz10[colIndex - 1] == 1 && vInROFcollIndex[colIndex - 1] == 0) {
                  histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInAnotherEarlierROF_1collPerROF"), vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                  if (vROFidThisColl[colIndex] == vROFidThisColl[colIndex - 1] + 1) // one ROF right after the previous
                  {
                    histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF"), vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);

                    if (vTracksITS567perColl[colIndex - 1] > vTracksITS567perColl[colIndex]) {
                      histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"), vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                      if (vTracksITS567perColl[colIndex - 1] > 0)
                        histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"), vTracksITS567perColl[colIndex - 1], 1.0 * vTracksITS567perColl[colIndex] / vTracksITS567perColl[colIndex - 1]);
                    } else {
                      histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"), vTracksITS567perColl[colIndex], vTracksITS567perColl[colIndex - 1]);
                      if (vTracksITS567perColl[colIndex] > 0)
                        histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"), vTracksITS567perColl[colIndex], 1.0 * vTracksITS567perColl[colIndex - 1] / vTracksITS567perColl[colIndex]);
                    }
                    // the sum of v1 and v2:
                    if (vAmpFT0CperColl[colIndex] > 4000 && vAmpFT0CperColl[colIndex - 1] > 4000)
                      histos.fill(HIST("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_thisROFprevROF"), vTracksITS567perColl[colIndex] + vTracksITS567perColl[colIndex - 1]);
                  } else if (vROFidThisColl[colIndex] == vROFidThisColl[colIndex - 1] + 2) {
                    // ROF vs ROF-2
                    histos.fill(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPrevPrevROF_1collPerROF"), vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                    if (vAmpFT0CperColl[colIndex] > 4000 && vAmpFT0CperColl[colIndex - 1] > 4000)
                      histos.fill(HIST("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_thisROFprevPrevROF"), vTracksITS567perColl[colIndex] + vTracksITS567perColl[colIndex - 1]);
                  } else {
                    // ROF is earlier than previous
                    // the sum of v1 and v2:
                    if (vAmpFT0CperColl[colIndex] > 4000 && vAmpFT0CperColl[colIndex - 1] > 4000)
                      histos.fill(HIST("afterNarrowDeltaTimeCut/hSum_2coll_withFT0above4000_thisROFearlierThanPrevPrevROF"), vTracksITS567perColl[colIndex] + vTracksITS567perColl[colIndex - 1]);
                  }
                }
              }
            }
            if (vNoCollInTimeRangeStrict[colIndex]) {
              histos.fill(HIST("afterStrictDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterStrictDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("afterStrictDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterStrictDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }
            if (vNoHighMultCollInTimeRange[colIndex]) {
              histos.fill(HIST("afterStandardDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterStandardDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("afterStandardDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterStandardDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInVzDependentTimeRange[colIndex]) {
              histos.fill(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInSameRofStrict[colIndex]) {
              histos.fill(HIST("kNoCollInRofStrict/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("kNoCollInRofStrict/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("kNoCollInRofStrict/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("kNoCollInRofStrict/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInSameRofStandard[colIndex]) {
              histos.fill(HIST("kNoCollInRofStandard/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("kNoCollInRofStandard/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("kNoCollInRofStandard/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("kNoCollInRofStandard/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInSameRofWithCloseVz[colIndex]) {
              histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);

              if (vNumTracksITS567inROF[colIndex] > 0) {
                histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyInROF_HasNeighbours"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours"), vSumAmpFT0CinROF[colIndex], vTracksITS567perColl[colIndex]);

                histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours"), vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
              }

              if (vNumCollinROF[colIndex] == 2 && vNumCollinROFinVz10[colIndex] == 2 && vInROFcollIndex[colIndex] == 1) {
                histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyInROF_2coll"), vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.fill(HIST("kNoCollInRofWithCloseVz/hThisEvFT0C_vs_occupancyFT0CInROF_2coll"), vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
              }
            }

            if (vNoCollInTimeRangeNarrow[colIndex] && vNoHighMultCollInTimeRange[colIndex] && vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
              histos.fill(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTPCTr_vs_occupancyByFT0C"), vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.fill(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.fill(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTPCTr_vs_occupancyByTracks"), vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(RofOccupancyQaTask, processRun3, "Process Run3 ROF occupancy QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<RofOccupancyQaTask>(cfgc)};
}
