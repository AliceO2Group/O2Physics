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

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "Framework/HistogramRegistry.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "MetadataHelper.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"

#include "TH1D.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

MetadataHelper metadataInfo; // Metadata helper

using BCsWithRun2InfosTimestampsAndMatches = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::Run2MatchedToBCSparse>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using BCsWithBcSelsRun2 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run2BCInfos, aod::Run2MatchedToBCSparse>;
using BCsWithBcSelsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;

struct RofOccupancyQaTask {
  SliceCache cache;
  // Produces<aod::EvSels> evsel;
  Configurable<int> muonSelection{"muonSelection", 0, "0 - barrel, 1 - muon selection with pileup cuts, 2 - muon selection without pileup cuts"};
  Configurable<float> maxDiffZvtxFT0vsPV{"maxDiffZvtxFT0vsPV", 1., "maximum difference (in cm) between z-vertex from FT0 and PV"};
  Configurable<int> isMC{"isMC", 0, "-1 - autoset, 0 - data, 1 - MC"};
  Configurable<int> confSigmaBCforHighPtTracks{"confSigmaBCforHighPtTracks", 4, "Custom sigma (in bcs) for collisions with high-pt tracks"};

  // configurables for occupancy-based event selection
  Configurable<float> confTimeIntervalForOccupancyCalculationMin{"TimeIntervalForOccupancyCalculationMin", -40, "Min time diff window for TPC occupancy calculation, us"};
  Configurable<float> confTimeIntervalForOccupancyCalculationMax{"TimeIntervalForOccupancyCalculationMax", 100, "Max time diff window for TPC occupancy calculation, us"};
  Configurable<float> confTimeRangeVetoOnCollStandard{"TimeRangeVetoOnCollStandard", 10.0, "Exclusion of a collision if there are other collisions nearby, +/- us"};
  Configurable<float> confTimeRangeVetoOnCollNarrow{"TimeRangeVetoOnCollNarrow", 2.0, "Exclusion of a collision if there are other collisions nearby, +/- us"};
  Configurable<int> confNtracksCutVetoOnCollInTimeRange{"NtracksCutVetoOnCollInTimeRange", 800, "Max allowed N tracks (PV contributors) for each nearby collision in +/- time range"};
  Configurable<float> confEpsilonDistanceForVzDependentVetoTPC{"EpsilonDistanceForVzDependentVetoTPC", 2.5, "Epsilon for vZ-dependent veto on drifting TPC tracks from nearby collisions, cm"};
  // Configurable<int> confNtracksCutVetoOnCollInROF{"NtracksCutVetoOnCollInROF", 500, "Max allowed N tracks (PV contributors) for each nearby collision inside this ITS ROF"};
  Configurable<float> confFT0CamplCutVetoOnCollInROF{"FT0CamplPerCollCutVetoOnCollInROF", 5000, "Max allowed FT0C amplitude for each nearby collision inside this ITS ROF"};
  Configurable<float> confEpsilonVzDiffVetoInROF{"EpsilonVzDiffVetoInROF", 0.3, "Minumum distance to nearby collisions along z inside this ITS ROF, cm"};
  Configurable<bool> confUseWeightsForOccupancyVariable{"UseWeightsForOccupancyEstimator", 1, "Use or not the delta-time weights for the occupancy estimator"};

  Configurable<float> confFactorForHistRange{"kFactorForHistRange", 1.0, "To change axes b/n pp and Pb-Pb"};

  Partition<aod::Tracks> tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int lastRun = -1;                                          // last run number (needed to access ccdb only if run!=lastRun)
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB; // bc pattern of colliding bunches

  int64_t bcSOR = -1;     // global bc of the start of the first orbit
  int64_t nBCsPerTF = -1; // duration of TF in bcs, should be 128*3564 or 32*3564
  int rofOffset = -1;     // ITS ROF offset, in bc
  int rofLength = -1;     // ITS ROF length, in bc

  int32_t findClosest(int64_t globalBC, std::map<int64_t, int32_t>& bcs)
  {
    auto it = bcs.lower_bound(globalBC);
    int64_t bc1 = it->first;
    int32_t index1 = it->second;
    if (it != bcs.begin())
      --it;
    int64_t bc2 = it->first;
    int32_t index2 = it->second;
    int64_t dbc1 = std::abs(bc1 - globalBC);
    int64_t dbc2 = std::abs(bc2 - globalBC);
    return (dbc1 <= dbc2) ? index1 : index2;
  }

  // helper function to find median time in the vector of TOF or TRD-track times
  float getMedian(std::vector<float> v)
  {
    int medianIndex = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + medianIndex, v.end());
    return v[medianIndex];
  }

  // helper function to find closest TVX signal in time and in zVtx
  int64_t findBestGlobalBC(int64_t meanBC, int64_t sigmaBC, int32_t nContrib, float zVtxCol, std::map<int64_t, float>& mapGlobalBcVtxZ)
  {
    int64_t minBC = meanBC - 3 * sigmaBC;
    int64_t maxBC = meanBC + 3 * sigmaBC;
    // TODO: use ITS ROF bounds to reduce the search range?

    float zVtxSigma = 2.7 * pow(nContrib, -0.466) + 0.024;
    zVtxSigma += 1.0; // additional uncertainty due to imperfectections of FT0 time calibration

    auto itMin = mapGlobalBcVtxZ.lower_bound(minBC);
    auto itMax = mapGlobalBcVtxZ.upper_bound(maxBC);

    float bestChi2 = 1e+10;
    int64_t bestGlobalBC = 0;
    for (std::map<int64_t, float>::iterator it = itMin; it != itMax; ++it) {
      float chi2 = pow((it->second - zVtxCol) / zVtxSigma, 2) + pow((it->first - meanBC) / sigmaBC, 2.);
      if (chi2 < bestChi2) {
        bestChi2 = chi2;
        bestGlobalBC = it->first;
      }
    }

    return bestGlobalBC;
  }

  void init(InitContext&)
  {
    if (metadataInfo.isFullyDefined()) { // Check if the metadata is initialized (only if not forced from the workflow configuration)
                                         // if (!doprocessRun2 && !doprocessRun3) {
                                         // if (!doprocessRun3) {
      LOG(info) << "Autosetting the processing mode (Run2 or Run3) based on metadata";
      if (metadataInfo.isRun3()) {
        doprocessRun3.value = true;
      } else {
        // doprocessRun2.value = false;
      }
      // }
      if (isMC == -1) {
        LOG(info) << "Autosetting the MC mode based on metadata";
        if (metadataInfo.isMC()) {
          isMC.value = 1;
        } else {
          isMC.value = 0;
        }
      }
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    histos.add("hColCounterAll", "", kTH1D, {{1, 0., 1.}});
    histos.add("hColCounterTVX", "", kTH1D, {{1, 0., 1.}});
    histos.add("hColCounterAcc", "", kTH1D, {{1, 0., 1.}});

    float k = confFactorForHistRange;
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
    histos.add("hThisEvITSTr_vs_ThisEvFT0C", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_vZ_TF_ROF_border_cuts", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_afterNarrowDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_afterStrictDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_afterStandardDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_afterVzDependentDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C_kNoCollInRofStrict", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_kNoCollInRofStandard", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_kNoCollInRofWithCloseVz", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C_occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvFT0C_occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    // CROSS-CHECK SEL BITS:
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterNarrowDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterStrictDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterStandardDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterVzDependentDeltaTimeCut", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});

    histos.add("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_kNoCollInRofStrict", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    histos.add("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_kNoCollInRofStandard", "", kTH2D, {{250, 0., 1e5 * k}, {250, 0., 10000 * k}});
    // histos.add("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_kNoCollInRofWithCloseVz", "", kTH2D, {{250, 0., 1e5*k}, {250, 0., 10000*k}});

    // this ev nITSTPCtr vs nITStr
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_vZ_TF_ROF_border_cuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_afterNarrowDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_afterStrictDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_afterStandardDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_afterVzDependentDeltaTimeCut", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_kNoCollInRofStrict", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_kNoCollInRofStandard", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_kNoCollInRofWithCloseVz", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_NarrowDeltaCut_StdTimeAndRofCuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_occupBelow2000", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_occupBelow2000_NarrowDeltaCut_StdTimeAndRofCuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});
    histos.add("hThisEvITSTPCTr_vs_ThisEvITStr_occupBelow2000_StrictDeltaTimeCutAndRofCuts", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 5000 * k}});

    histos.add("hThisEvITStr_vs_vZcut", "", kTH2D, {{200, 0., 10.}, {200, 0., 8000 * k}});
    histos.add("hThisEvITSTPCtr_vs_vZcut", "", kTH2D, {{200, 0., 10.}, {200, 0., 5000 * k}});

    histos.add("hThisEvITStr_vs_vZ", "", kTH2D, {{400, -20, 20.}, {200, 0., 8000 * k}});
    histos.add("hThisEvITSTPCtr_vs_vZ", "", kTH2D, {{400, -20, 20.}, {200, 0., 5000 * k}});

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

    // coll on x axis always has more tracks:
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins", "", kTH2D, {{250, 0., 8000 * k}, {250, 0., 8000 * k}});

    // now with the ratio on y-axis:
    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});

    histos.add("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});
    histos.add("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1", "", kTH2D, {{250, 0., 8000 * k}, {220, 0., 1.1}});

    //
    histos.add("hNcollPerROF", "", kTH1D, {{16, -0.5, 15.5}});
  }

  Partition<FullTracksIU> pvTracks = ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Preslice<FullTracksIU> perCollision = aod::track::collisionId;

  using ColEvSels = soa::Join<aod::Collisions, aod::EvSels>; //, aod::Mults, aod::CentFT0Cs>;
  // void processRun3(aod::Collisions const& cols, FullTracksIU const&, BCsWithBcSelsRun3 const& bcs, aod::FT0s const&)
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
      // colliding bc pattern
      int64_t ts = bcs.iteratorAt(0).timestamp();
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", ts);
      bcPatternB = grplhcif->getBunchFilling().getBCPattern();

      // extract ITS ROF parameters
      // auto timestamps = ccdb->getRunDuration(run, true); /// fatalise if timestamps are not found
      // int64_t sorTimestamp = timestamps.first;           // timestamp of the SOR/SOX/STF in ms
      // int64_t eorTimestamp = timestamps.second;          // timestamp of the EOR/EOX/ETF in ms
      // int64_t ts = eorTimestamp / 2 + sorTimestamp / 2;  // timestamp of the middle of the run
      auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
      rofOffset = alppar->roFrameBiasInBC;
      rofLength = alppar->roFrameLengthInBC;
      LOGP(info, "rofOffset={} rofLength={}", rofOffset, rofLength);
    }

    // create maps from globalBC to bc index for TVX-fired bcs
    // to be used for closest TVX searches
    std::map<int64_t, int32_t> mapGlobalBcWithTVX;
    std::map<int64_t, float> mapGlobalBcVtxZ;
    for (auto& bc : bcs) {
      int64_t globalBC = bc.globalBC();
      // skip non-colliding bcs for data and anchored runs
      if (run >= 500000 && bcPatternB[globalBC % o2::constants::lhc::LHCMaxBunches] == 0) {
        continue;
      }
      if (bc.selection_bit(kIsTriggerTVX)) {
        mapGlobalBcWithTVX[globalBC] = bc.globalIndex();
        mapGlobalBcVtxZ[globalBC] = bc.has_ft0() ? bc.ft0().posZ() : 0;
      }
    }

    // protection against empty FT0 maps
    if (mapGlobalBcWithTVX.size() == 0) {
      LOGP(error, "FT0 table is empty or corrupted.");
      // for (auto& col : cols) {
      //   auto bc = col.bc_as<BCsWithBcSelsRun3>();
      //   int32_t foundBC = bc.globalIndex();
      //   int32_t foundFT0 = bc.foundFT0Id();
      //   int32_t foundFV0 = bc.foundFV0Id();
      //   int32_t foundFDD = bc.foundFDDId();
      //   int32_t foundZDC = bc.foundZDCId();
      //   int bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
      //   // evsel(bc.alias_raw(), bc.selection_raw(), kFALSE, kFALSE, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, bcInTF, -1, -1, -1);
      // }
      return;
    }
    std::vector<int> vTracksITS567perColl(cols.size(), 0);                                    // counter of tracks per collision for occupancy studies
    std::vector<int> vTracksITSTPCperColl(cols.size(), 0);                                    // counter of tracks per collision for occupancy studies
    std::vector<float> vAmpFT0CperColl(cols.size(), 0);                                       // amplitude FT0C per collision
    std::vector<float> vCollVz(cols.size(), 0);                                               // vector with vZ positions for each collision
    std::vector<bool> vIsFullInfoForOccupancy(cols.size(), 0);                                // info for occupancy in +/- windows is available (i.e. a given coll is not too close to the TF borders)
    const float timeWinOccupancyCalcMinNS = confTimeIntervalForOccupancyCalculationMin * 1e3; // ns
    const float timeWinOccupancyCalcMaxNS = confTimeIntervalForOccupancyCalculationMax * 1e3; // ns
    std::vector<bool> vIsVertexITSTPC(cols.size(), 0);                                        // at least one of vertex contributors is ITS-TPC track
    std::vector<bool> vIsVertexTOFmatched(cols.size(), 0);                                    // at least one of vertex contributors is matched to TOF
    std::vector<bool> vIsVertexTRDmatched(cols.size(), 0);                                    // at least one of vertex contributors is matched to TRD

    std::vector<int> vCollisionsPerBc(bcs.size(), 0);    // counter of collisions per found bc for pileup checks
    std::vector<int> vFoundBCindex(cols.size(), -1);     // indices of found bcs
    std::vector<int64_t> vFoundGlobalBC(cols.size(), 0); // global BCs for collisions

    std::vector<bool> vIsVertexTOF(cols.size(), 0);
    std::vector<bool> vIsVertexTRD(cols.size(), 0);
    std::vector<bool> vIsVertexTPC(cols.size(), 0);
    std::vector<bool> vIsVertexHighPtTPC(cols.size(), 0);
    std::vector<int> vNcontributors(cols.size(), 0);
    std::vector<float> vWeightedTimesTPCnoTOFnoTRD(cols.size(), 0);
    std::vector<float> vWeightedSigmaTPCnoTOFnoTRD(cols.size(), 0);

    // temporary vectors to find tracks with median time
    std::vector<float> vTrackTimesTOF;
    std::vector<float> vTrackTimesTRDnoTOF;

    // first loop to match collisions to TVX, also extract other per-collision information for further use
    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      auto bc = col.bc_as<BCsWithBcSelsRun3>();

      vCollVz[colIndex] = col.posZ();

      int64_t globalBC = bc.globalBC();
      int bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
      vIsFullInfoForOccupancy[colIndex] = ((bcInTF - 300) * bcNS > -timeWinOccupancyCalcMinNS) && ((nBCsPerTF - 4000 - bcInTF) * bcNS > timeWinOccupancyCalcMaxNS) ? true : false;

      // const auto& colPvTracks = pvTracks.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
      // const auto& colPvTracks = pvTracks.sliceBy(aod::track::collisionId, col.globalIndex(), cache);
      auto colPvTracks = pvTracks.sliceBy(perCollision, col.globalIndex());

      vTrackTimesTOF.clear();
      vTrackTimesTRDnoTOF.clear();
      int nPvTracksTPCnoTOFnoTRD = 0;
      int nPvTracksHighPtTPCnoTOFnoTRD = 0;
      float sumTime = 0, sumW = 0, sumHighPtTime = 0, sumHighPtW = 0;
      for (auto& track : colPvTracks) {
        float trackTime = track.trackTime();
        if (track.itsNCls() >= 5) {
          vTracksITS567perColl[colIndex]++;
          if (track.tpcNClsFound() > 70)
            vTracksITSTPCperColl[colIndex]++;
          if (fabs(col.posZ()) < 1)
            histos.get<TH1>(HIST("hEtaVz02"))->Fill(track.eta());
          else if (col.posZ() > 8 && col.posZ() < 10)
            histos.get<TH1>(HIST("hEtaVzPlus10"))->Fill(track.eta());
          else if (col.posZ() > -10 && col.posZ() < -8)
            histos.get<TH1>(HIST("hEtaVzMinus10"))->Fill(track.eta());
          else if (col.posZ() > 14 && col.posZ() < 16)
            histos.get<TH1>(HIST("hEtaVzPlus15"))->Fill(track.eta());
          else if (col.posZ() > -16 && col.posZ() < -14)
            histos.get<TH1>(HIST("hEtaVzMinus15"))->Fill(track.eta());
        }
        if (track.hasTRD())
          vIsVertexTRDmatched[colIndex] = 1;
        if (track.hasTPC())
          vIsVertexITSTPC[colIndex] = 1;
        if (track.hasTOF()) {
          vTrackTimesTOF.push_back(trackTime);
          vIsVertexTOFmatched[colIndex] = 1;
        } else if (track.hasTRD()) {
          vTrackTimesTRDnoTOF.push_back(trackTime);
        } else if (track.hasTPC()) {
          float trackTimeRes = track.trackTimeRes();
          float trackPt = track.pt();
          float w = 1. / (trackTimeRes * trackTimeRes);
          sumTime += trackTime * w;
          sumW += w;
          nPvTracksTPCnoTOFnoTRD++;
          if (trackPt > 1) {
            sumHighPtTime += trackTime * w;
            sumHighPtW += w;
            nPvTracksHighPtTPCnoTOFnoTRD++;
          }
        }
      }
      vWeightedTimesTPCnoTOFnoTRD[colIndex] = sumW > 0 ? sumTime / sumW : 0;
      vWeightedSigmaTPCnoTOFnoTRD[colIndex] = sumW > 0 ? sqrt(1. / sumW) : 0;
      vNcontributors[colIndex] = colPvTracks.size();
      int nPvTracksTOF = vTrackTimesTOF.size();
      int nPvTracksTRDnoTOF = vTrackTimesTRDnoTOF.size();
      // collision type
      vIsVertexTOF[colIndex] = nPvTracksTOF > 0;
      vIsVertexTRD[colIndex] = nPvTracksTRDnoTOF > 0;
      vIsVertexTPC[colIndex] = nPvTracksTPCnoTOFnoTRD > 0;
      vIsVertexHighPtTPC[colIndex] = nPvTracksHighPtTPCnoTOFnoTRD > 0;

      int64_t foundGlobalBC = 0;
      int32_t foundBCindex = -1;

      if (nPvTracksTOF > 0) {
        // for collisions with TOF tracks:
        // take bc corresponding to TOF track with median time
        int64_t tofGlobalBC = globalBC + TMath::Nint(getMedian(vTrackTimesTOF) / bcNS);
        std::map<int64_t, int32_t>::iterator it = mapGlobalBcWithTVX.find(tofGlobalBC);
        if (it != mapGlobalBcWithTVX.end()) {
          foundGlobalBC = it->first;
          foundBCindex = it->second;
        }
      } else if (nPvTracksTPCnoTOFnoTRD == 0 && nPvTracksTRDnoTOF > 0) {
        // for collisions with TRD tracks but without TOF or ITSTPC-only tracks:
        // take bc corresponding to TRD track with median time
        int64_t trdGlobalBC = globalBC + TMath::Nint(getMedian(vTrackTimesTRDnoTOF) / bcNS);
        std::map<int64_t, int32_t>::iterator it = mapGlobalBcWithTVX.find(trdGlobalBC);
        if (it != mapGlobalBcWithTVX.end()) {
          foundGlobalBC = it->first;
          foundBCindex = it->second;
        }
      } else if (nPvTracksHighPtTPCnoTOFnoTRD > 0) {
        // for collisions with high-pt ITSTPC-nonTOF-nonTRD tracks
        // search in 3*confSigmaBCforHighPtTracks range (3*4 bcs by default)
        int64_t meanBC = globalBC + TMath::Nint(sumHighPtTime / sumHighPtW / bcNS);
        int64_t bestGlobalBC = findBestGlobalBC(meanBC, confSigmaBCforHighPtTracks, vNcontributors[colIndex], col.posZ(), mapGlobalBcVtxZ);
        if (bestGlobalBC > 0) {
          foundGlobalBC = bestGlobalBC;
          foundBCindex = mapGlobalBcWithTVX[bestGlobalBC];
        }
      }

      // fill foundBC indices and global BCs
      // keep current bc if TVX matching failed at this step
      vFoundBCindex[colIndex] = foundBCindex >= 0 ? foundBCindex : bc.globalIndex();
      vFoundGlobalBC[colIndex] = foundGlobalBC > 0 ? foundGlobalBC : globalBC;

      // erase found global BC with TVX from the pool of bcs for the next loop over low-pt TPCnoTOFnoTRD collisions
      if (foundBCindex >= 0)
        mapGlobalBcVtxZ.erase(foundGlobalBC);
    }

    // second loop to match remaining low-pt TPCnoTOFnoTRD collisions
    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      if (vIsVertexTPC[colIndex] > 0 && vIsVertexTOF[colIndex] == 0 && vIsVertexHighPtTPC[colIndex] == 0) {
        float weightedTime = vWeightedTimesTPCnoTOFnoTRD[colIndex];
        float weightedSigma = vWeightedSigmaTPCnoTOFnoTRD[colIndex];
        auto bc = col.bc_as<BCsWithBcSelsRun3>();
        int64_t globalBC = bc.globalBC();
        int64_t meanBC = globalBC + TMath::Nint(weightedTime / bcNS);
        int64_t bestGlobalBC = findBestGlobalBC(meanBC, weightedSigma / bcNS, vNcontributors[colIndex], col.posZ(), mapGlobalBcVtxZ);
        vFoundGlobalBC[colIndex] = bestGlobalBC > 0 ? bestGlobalBC : globalBC;
        vFoundBCindex[colIndex] = bestGlobalBC > 0 ? mapGlobalBcWithTVX[bestGlobalBC] : bc.globalIndex();
      }
      // fill pileup counter
      vCollisionsPerBc[vFoundBCindex[colIndex]]++;
    }

    // save indices of collisions in time range for occupancy calculation
    std::vector<std::vector<int>> vCollsInTimeWin;
    std::vector<std::vector<int>> vCollsInSameITSROF;
    std::vector<std::vector<float>> vTimeDeltaForColls; // delta time wrt a given collision
    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];

      auto bc = bcs.iteratorAt(vFoundBCindex[colIndex]);
      if (bc.has_foundFT0())
        vAmpFT0CperColl[colIndex] = bc.foundFT0().sumAmpC();

      // int bcInTF = (foundGlobalBC - bcSOR) % nBCsPerTF;
      // int bcInITSROF = (foundGlobalBC + 3564 - rofOffset) % rofLength;
      int64_t TFid = (foundGlobalBC - bcSOR) / nBCsPerTF;
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
        if (thisTFid != TFid)
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
        if (thisTFid != TFid)
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
        if (thisTFid != TFid)
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
        if (thisTFid != TFid)
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

    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      float vZ = col.posZ();

      // QA:
      if (vAmpFT0CperColl[colIndex] > 5000) {
        histos.get<TH2>(HIST("hThisEvITStr_vs_vZ"))->Fill(vZ, vTracksITS567perColl[colIndex]);
        histos.get<TH2>(HIST("hThisEvITSTPCtr_vs_vZ"))->Fill(vZ, vTracksITSTPCperColl[colIndex]);
      }

      // ### in-ROF occupancy
      // int64_t rofId = (vFoundGlobalBC[colIndex] + 3564 - rofOffset) / rofLength;
      // int bcInTF = (vFoundGlobalBC[colIndex] - bcSOR) % nBCsPerTF;
      int bcInITSROF = (vFoundGlobalBC[colIndex] + 3564 - rofOffset) % rofLength;
      int rofIdInTF = (vFoundGlobalBC[colIndex] + 3564 - rofOffset) / rofLength;
      auto bc = bcs.iteratorAt(vFoundBCindex[colIndex]);
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

      if (fabs(vZ) < 10)
        vNumCollinROFinVz10[colIndex]++;
      for (uint32_t iCol = 0; iCol < vAssocToSameROF.size(); iCol++) {
        int thisColIndex = vAssocToSameROF[iCol];
        // int64_t thisRofId = (vFoundGlobalBC[thisColIndex] + 3564 - rofOffset) / rofLength;
        // int thisBcInTF = (vFoundGlobalBC[thisColIndex] - bcSOR) % nBCsPerTF;
        int thisBcInITSROF = (vFoundGlobalBC[thisColIndex] + 3564 - rofOffset) % rofLength;
        // int thisRofIdInTF = (vFoundGlobalBC[thisColIndex] + 3564 - rofOffset) / rofLength;
        // auto bcAssoc = bcs.iteratorAt(vFoundBCindex[thisColIndex]);
        // LOGP(info, ">> assoc: bc={} bcInTF={} bcInITSROF={} rofId={} noROFborder={}", vFoundGlobalBC[thisColIndex], thisBcInTF, thisBcInITSROF, thisRofId, bcAssoc.selection_bit(kNoITSROFrameBorder));
        // LOGP(info, ">> assoc: bcInTF={} bcInITSROF={} rofIdInTF={} noROFborder={} vZ={} mult={}", thisBcInTF, thisBcInITSROF, thisRofIdInTF, bcAssoc.selection_bit(kNoITSROFrameBorder), vCollVz[thisColIndex], vTracksITS567perColl[thisColIndex]);

        // if (fabs(vTracksITS567perColl[thisColIndex]) > confNtracksCutVetoOnCollInROF)
        nITS567tracksForRofVetoStrict += vTracksITS567perColl[thisColIndex];
        nSumAmplFT0CforRofVetoStrict += vAmpFT0CperColl[thisColIndex];
        vNumCollinROF[colIndex]++;
        vInROFcollIndex[colIndex] = thisBcInITSROF > bcInITSROF ? 0 : 1; // if colIndex is for the first coll in ROF => inROFindex=0, otherwise =1
        if (fabs(vCollVz[thisColIndex]) < 10)
          vNumCollinROFinVz10[colIndex]++;

        // if (vTracksITS567perColl[thisColIndex] > confNtracksCutVetoOnCollInROF)
        // nITS567tracksForRofVetoStandard += vTracksITS567perColl[thisColIndex];

        if (vAmpFT0CperColl[thisColIndex] > confFT0CamplCutVetoOnCollInROF)
          nCollsInRofWithFT0CAboveVetoStandard++;

        if (fabs(vCollVz[thisColIndex] - vZ) < confEpsilonVzDiffVetoInROF)
          nITS567tracksForRofVetoOnCloseVz += vTracksITS567perColl[thisColIndex];
        for (int i = 0; i < 200; i++) {
          // if (fabs(vCollVz[thisColIndex] - vZ) < 0.05 * i && vTracksITS567perColl[thisColIndex] > 50)
          // if (vTracksITS567perColl[colIndex]>100 && vTracksITS567perColl[colIndex]<1000 &&
          if (vAmpFT0CperColl[colIndex] > 4000 && vAmpFT0CperColl[colIndex] < 15000 &&
              (vCollVz[thisColIndex] - vZ) > 0.05 * i && fabs(vCollVz[thisColIndex] - vZ) < (0.1 + 0.05) * i && vTracksITS567perColl[thisColIndex] > 20) // 0.05 * (i + 1))
            // fabs(vCollVz[thisColIndex] - vZ) < 0.05 * i && vTracksITS567perColl[thisColIndex] > 30) // 0.05 * (i + 1))
            nArrITS567tracksForRofVetoOnCloseVz[i]++;
        }

        if (fabs(vZ) < 8.) {
          histos.get<TH1>(HIST("hDeltaVz"))->Fill(vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] >= 100 && vTracksITS567perColl[thisColIndex] < 100)
            histos.get<TH1>(HIST("hDeltaVzGivenCollAbove100NearbyBelow100"))->Fill(vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] <= 100 && vTracksITS567perColl[thisColIndex] > 100)
            histos.get<TH1>(HIST("hDeltaVzGivenCollBelow100NearbyAbove100"))->Fill(vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] > 20 && vTracksITS567perColl[thisColIndex] > 20)
            histos.get<TH1>(HIST("hDeltaVzAfterCuts"))->Fill(vCollVz[thisColIndex] - vZ);
          if (bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder) && bc.selection_bit(kNoITSROFrameBorder))
            histos.get<TH1>(HIST("hDeltaVzAfterTFandROFborderCuts"))->Fill(vCollVz[thisColIndex] - vZ);
          if (vTracksITS567perColl[colIndex] > 20 && vTracksITS567perColl[thisColIndex] > 20 && bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder) && bc.selection_bit(kNoITSROFrameBorder))
            histos.get<TH1>(HIST("hDeltaVzAfterAllCuts"))->Fill(vCollVz[thisColIndex] - vZ);
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
            histos.get<TH1>(HIST("hOccupancyWeights"))->Fill(dt, wOccup);
        }
        nITS567tracksInFullTimeWindow += wOccup * vTracksITS567perColl[thisColIndex];
        sumAmpFT0CInFullTimeWindow += wOccup * vAmpFT0CperColl[thisColIndex];

        // counting tracks from other collisions in fixed time windows
        if (fabs(dt) < confTimeRangeVetoOnCollNarrow)
          nITS567tracksForVetoNarrow += vTracksITS567perColl[thisColIndex];
        if (fabs(dt) < confTimeRangeVetoOnCollStandard)
          nITS567tracksForVetoStrict += vTracksITS567perColl[thisColIndex];

        // if (fabs(dt) < confTimeRangeVetoOnCollStandard + 0.5) { // add 0.5 us safety margin
        // standard cut on other collisions vs delta-times
        const float driftV = 2.5; // drift velocity in cm/us, TPC drift_length / drift_time = 250 cm / 100 us
        if (fabs(dt) < 2.0) {     // us, complete veto on other collisions
          nITS567tracksForVetoStandard += vTracksITS567perColl[thisColIndex];
        } else if (dt > -4.0 && dt <= -2.0) { // us, strict veto to suppress fake ITS-TPC matches more
          if (vTracksITS567perColl[thisColIndex] > confNtracksCutVetoOnCollInTimeRange / 5)
            nITS567tracksForVetoStandard += vTracksITS567perColl[thisColIndex];
        } else if (fabs(dt) < 8 + fabs(vZ) / driftV) { // loose veto, 8 us corresponds to maximum possible |vZ|, which is ~20 cm
          // counting number of other collisions with mult above threshold
          if (vTracksITS567perColl[thisColIndex] > confNtracksCutVetoOnCollInTimeRange)
            nITS567tracksForVetoStandard += vTracksITS567perColl[thisColIndex];
        }
        // vZ-dependent time cut to avoid collinear tracks from other collisions (experimental)
        if (fabs(dt) < 8 + fabs(vZ) / driftV) {
          if (dt < 0) {
            // check distance between given vZ and (moving in two directions) vZ of drifting tracks from past collisions
            if ((fabs(vCollVz[thisColIndex] - fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC) ||
                (fabs(vCollVz[thisColIndex] + fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC))
              nITS567tracksForVetoVzDependent += vTracksITS567perColl[thisColIndex];

            // FOR QA:
            if (fabs(vCollVz[thisColIndex] - fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.get<TH2>(HIST("hDeltaVzVsDeltaTime1"))->Fill(vCollVz[thisColIndex] - fabs(dt) * driftV - vZ, dt);
            if (fabs(vCollVz[thisColIndex] + fabs(dt) * driftV - vZ) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.get<TH2>(HIST("hDeltaVzVsDeltaTime2"))->Fill(vCollVz[thisColIndex] + fabs(dt) * driftV - vZ, dt);

          } else { // dt>0
            // check distance between drifted vZ of given collision (in two directions) and vZ of future collisions
            if ((fabs(vZ - dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC) ||
                (fabs(vZ + dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC))
              nITS567tracksForVetoVzDependent += vTracksITS567perColl[thisColIndex];

            // FOR QA:
            if (fabs(vZ - dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.get<TH2>(HIST("hDeltaVzVsDeltaTime3"))->Fill(vZ - dt * driftV - vCollVz[thisColIndex], dt);
            if (fabs(vZ + dt * driftV - vCollVz[thisColIndex]) < confEpsilonDistanceForVzDependentVetoTPC)
              histos.get<TH2>(HIST("hDeltaVzVsDeltaTime4"))->Fill(vZ + dt * driftV - vCollVz[thisColIndex], dt);
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

    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      int32_t foundBC = vFoundBCindex[colIndex];
      auto bc = bcs.iteratorAt(foundBC);
      // int32_t foundFT0 = bc.foundFT0Id();
      // int32_t foundFV0 = bc.foundFV0Id();
      // int32_t foundFDD = bc.foundFDDId();
      // int32_t foundZDC = bc.foundZDCId();

      // // compare zVtx from FT0 and from PV
      // bool isGoodZvtxFT0vsPV = bc.has_foundFT0() ? fabs(bc.foundFT0().posZ() - col.posZ()) < maxDiffZvtxFT0vsPV : 0;

      // // copy alias decisions from bcsel table
      // uint32_t alias = bc.alias_raw();

      // // copy selection decisions from bcsel table
      // uint64_t selection = bc.selection_raw();
      // selection |= vCollisionsPerBc[foundBC] <= 1 ? BIT(kNoSameBunchPileup) : 0;
      // selection |= vIsVertexITSTPC[colIndex] ? BIT(kIsVertexITSTPC) : 0;
      // selection |= vIsVertexTOFmatched[colIndex] ? BIT(kIsVertexTOFmatched) : 0;
      // selection |= vIsVertexTRDmatched[colIndex] ? BIT(kIsVertexTRDmatched) : 0;
      // selection |= isGoodZvtxFT0vsPV ? BIT(kIsGoodZvtxFT0vsPV) : 0;

      // // selection bits based on occupancy time pattern
      // selection |= vNoCollInTimeRangeNarrow[colIndex] ? BIT(kNoCollInTimeRangeNarrow) : 0;
      // selection |= vNoCollInTimeRangeStrict[colIndex] ? BIT(kNoCollInTimeRangeStrict) : 0;
      // selection |= vNoHighMultCollInTimeRange[colIndex] ? BIT(kNoCollInTimeRangeStandard) : 0;
      // selection |= vNoCollInVzDependentTimeRange[colIndex] ? BIT(kNoCollInTimeRangeVzDependent) : 0;

      // selection |= vNoCollInSameRofStrict[colIndex] ? BIT(kNoCollInRofStrict) : 0;
      // selection |= (vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) ? BIT(kNoCollInRofStandard) : 0;
      // selection |= vNoCollInSameRofWithCloseVz[colIndex] ? BIT(kNoCollInRofWithCloseVz) : 0;

      // apply int7-like selections
      // bool sel7 = 0;

      // TODO apply other cuts for sel8
      // TODO introduce sel1 etc?
      // TODO introduce array of sel[0]... sel[8] or similar?
      bool sel8 = bc.selection_bit(kIsTriggerTVX) && bc.selection_bit(kNoTimeFrameBorder) && bc.selection_bit(kNoITSROFrameBorder);

      // fill counters
      histos.get<TH1>(HIST("hColCounterAll"))->Fill(Form("%d", bc.runNumber()), 1);
      if (bc.selection_bit(kIsTriggerTVX)) {
        histos.get<TH1>(HIST("hColCounterTVX"))->Fill(Form("%d", bc.runNumber()), 1);
      }
      if (sel8) {
        histos.get<TH1>(HIST("hColCounterAcc"))->Fill(Form("%d", bc.runNumber()), 1);
        if (vInROFcollIndex[colIndex] == 0)
          histos.get<TH1>(HIST("hNcollPerROF"))->Fill(vNumCollinROF[colIndex]);
      }

      if (col.trackOccupancyInTimeRange() >= 0)
        histos.get<TH1>(HIST("hOccupancyByTracks_CROSSCHECK"))->Fill(col.trackOccupancyInTimeRange());
      if (col.ft0cOccupancyInTimeRange() >= 0)
        histos.get<TH1>(HIST("hOccupancyByFT0C_CROSSCHECK"))->Fill(col.ft0cOccupancyInTimeRange());

      if (vNumTracksITS567inFullTimeWin[colIndex] >= 0) {
        histos.get<TH1>(HIST("hOccupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex]);
        histos.get<TH1>(HIST("hOccupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex]);

        histos.get<TH1>(HIST("hOccupancyByTrInROF"))->Fill(vNumTracksITS567inROF[colIndex]);
        histos.get<TH2>(HIST("hOccupancyByFT0C_vs_ByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
        histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
        histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

        if (sel8 && fabs(col.posZ()) < 10) {
          histos.get<TH2>(HIST("hOccupancyByFT0C_vs_ByTracks_vZ_TF_ROF_border_cuts"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);

          // if (vAmpFT0CperColl[colIndex] > 5000 && vAmpFT0CperColl[colIndex] < 10000) {
          // if (vAmpFT0CperColl[colIndex] > 500) {
          if (vAmpFT0CperColl[colIndex] > 0) { // 50) {

            histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
            histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
            histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);

            histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
            histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTPCTr_vs_occupancyInROF"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITSTPCperColl[colIndex]);

            if (vNumTracksITS567inROF[colIndex] > 0) {
              histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_HasNeighbours"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours"))->Fill(vSumAmpFT0CinROF[colIndex], vTracksITS567perColl[colIndex]);

              histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours"))->Fill(vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
            }

            if (vNumCollinROF[colIndex] == 2 && vNumCollinROFinVz10[colIndex] == 2 && vInROFcollIndex[colIndex] == 1) {
              histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvFT0C_vs_occupancyFT0CInROF_2coll"))->Fill(vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);

              if (vNumTracksITS567inROF[colIndex] > vTracksITS567perColl[colIndex]) {
                histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                if (vNumTracksITS567inROF[colIndex] > 0)
                  histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"))->Fill(vNumTracksITS567inROF[colIndex], 1.0 * vTracksITS567perColl[colIndex] / vNumTracksITS567inROF[colIndex]);
              } else {
                histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"))->Fill(vTracksITS567perColl[colIndex], vNumTracksITS567inROF[colIndex]);
                if (vTracksITS567perColl[colIndex] > 0)
                  histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"))->Fill(vTracksITS567perColl[colIndex], 1.0 * vNumTracksITS567inROF[colIndex] / vTracksITS567perColl[colIndex]);
              }
            }
            // compare with previous ROF
            if (colIndex - 1 >= 0) {
              if (vNumCollinROF[colIndex] == 1 && vNumCollinROFinVz10[colIndex] == 1 && vInROFcollIndex[colIndex] == 0 && vNumCollinROF[colIndex - 1] == 1 && vNumCollinROFinVz10[colIndex - 1] == 1 && vInROFcollIndex[colIndex - 1] == 0) {
                histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInAnotherEarlierROF_1collPerROF"))->Fill(vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                if (vROFidThisColl[colIndex] == vROFidThisColl[colIndex - 1] + 1) // one ROF right after the previous
                {
                  histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF"))->Fill(vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                  if (vTracksITS567perColl[colIndex - 1] > vTracksITS567perColl[colIndex]) {
                    histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"))->Fill(vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                    if (vTracksITS567perColl[colIndex - 1] > 0)
                      histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"))->Fill(vTracksITS567perColl[colIndex - 1], 1.0 * vTracksITS567perColl[colIndex] / vTracksITS567perColl[colIndex - 1]);
                  } else {
                    histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"))->Fill(vTracksITS567perColl[colIndex], vTracksITS567perColl[colIndex - 1]);
                    if (vTracksITS567perColl[colIndex] > 0)
                      histos.get<TH2>(HIST("vZ_TF_ROF_border_cuts/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"))->Fill(vTracksITS567perColl[colIndex], 1.0 * vTracksITS567perColl[colIndex - 1] / vTracksITS567perColl[colIndex]);
                  }
                }
              }
            }
          }

          if (vNoCollInTimeRangeNarrow[colIndex])
            histos.get<TH2>(HIST("hOccupancyByFT0C_vs_ByTracks_afterNarrowDeltaTimeCut"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
          if (vNoCollInTimeRangeStrict[colIndex])
            histos.get<TH2>(HIST("hOccupancyByFT0C_vs_ByTracks_afterStrictDeltaTimeCut"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
          if (vNoHighMultCollInTimeRange[colIndex])
            histos.get<TH2>(HIST("hOccupancyByFT0C_vs_ByTracks_afterStandardDeltaTimeCut"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);
          if (vNoCollInVzDependentTimeRange[colIndex])
            histos.get<TH2>(HIST("hOccupancyByFT0C_vs_ByTracks_afterVzDependentDeltaTimeCut"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex]);

          // same-event 2D correlations:
          histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_vZ_TF_ROF_border_cuts"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (vNoCollInTimeRangeNarrow[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_afterNarrowDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInTimeRangeStrict[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_afterStrictDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoHighMultCollInTimeRange[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_afterStandardDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInVzDependentTimeRange[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_afterVzDependentDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (vNoCollInSameRofStrict[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_kNoCollInRofStrict"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInSameRofStandard[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_kNoCollInRofStandard"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (vNoCollInSameRofWithCloseVz[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_kNoCollInRofWithCloseVz"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          // CROSS CHECK WITH SEL BITS:
          if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow))
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterNarrowDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict))
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterStrictDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterStandardDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeVzDependent))
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_afterVzDependentDeltaTimeCut"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (col.selection_bit(o2::aod::evsel::kNoCollInRofStrict))
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_kNoCollInRofStrict"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
          if (col.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_CROSSCHECK_kNoCollInRofStandard"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

          if (vNumTracksITS567inFullTimeWin[colIndex] < 2000) {
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_occupBelow2000"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvFT0C_occupBelow2000"))->Fill(vAmpFT0CperColl[colIndex], vTracksITSTPCperColl[colIndex]);
          }

          if (vNoCollInTimeRangeNarrow[colIndex] && vNoHighMultCollInTimeRange[colIndex] && vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
            histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);

            if (vNumTracksITS567inFullTimeWin[colIndex] < 2000) {
              histos.get<TH2>(HIST("hThisEvITSTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000"))->Fill(vAmpFT0CperColl[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvFT0C_NarrowDeltaCut_StdTimeAndRofCuts_occupBelow2000"))->Fill(vAmpFT0CperColl[colIndex], vTracksITSTPCperColl[colIndex]);
            }
          }

          // now ITSTPC vs ITS tr (this event)
          histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_vZ_TF_ROF_border_cuts"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNoCollInTimeRangeNarrow[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_afterNarrowDeltaTimeCut"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInTimeRangeStrict[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_afterStrictDeltaTimeCut"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoHighMultCollInTimeRange[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_afterStandardDeltaTimeCut"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInVzDependentTimeRange[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_afterVzDependentDeltaTimeCut"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNoCollInSameRofStrict[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_kNoCollInRofStrict"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInSameRofStandard[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_kNoCollInRofStandard"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          if (vNoCollInSameRofWithCloseVz[colIndex])
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_kNoCollInRofWithCloseVz"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNumTracksITS567inFullTimeWin[colIndex] < 2000)
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_occupBelow2000"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);

          if (vNoCollInTimeRangeNarrow[colIndex] && vNoHighMultCollInTimeRange[colIndex] && vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
            histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_NarrowDeltaCut_StdTimeAndRofCuts"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
            if (vNumTracksITS567inFullTimeWin[colIndex] < 2000)
              histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_occupBelow2000_NarrowDeltaCut_StdTimeAndRofCuts"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          }

          if (vNoCollInTimeRangeStrict[colIndex] && vNoCollInSameRofStrict[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
            if (vNumTracksITS567inFullTimeWin[colIndex] < 2000)
              histos.get<TH2>(HIST("hThisEvITSTPCTr_vs_ThisEvITStr_occupBelow2000_StrictDeltaTimeCutAndRofCuts"))->Fill(vTracksITS567perColl[colIndex], vTracksITSTPCperColl[colIndex]);
          }

          // vZ bins to tune vZthresh cut
          if (vNoCollInTimeRangeNarrow[colIndex]) {

            for (int i = 0; i < 200; i++) {
              if (fabs(col.posZ()) < 8 && !vArrNoCollInSameRofWithCloseVz[colIndex][i]) {
                histos.get<TH2>(HIST("hThisEvITStr_vs_vZcut"))->Fill(0.025 + 0.05 * i, vTracksITS567perColl[colIndex]);
                histos.get<TH2>(HIST("hThisEvITSTPCtr_vs_vZcut"))->Fill(0.025 + 0.05 * i, vTracksITSTPCperColl[colIndex]);
              }
            }
          }

          // ### this event vs Occupancy 2D histos
          // if (vAmpFT0CperColl[colIndex] > 5000 && vAmpFT0CperColl[colIndex] < 10000) {
          // if (vAmpFT0CperColl[colIndex] > 500) {
          if (vAmpFT0CperColl[colIndex] > 0) { // 100) {

            if (vNoCollInTimeRangeNarrow[colIndex]) {
              histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);

              histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyInROF"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITSTPCperColl[colIndex]);

              if (vNumTracksITS567inROF[colIndex] > 0) {
                histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_HasNeighbours"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours"))->Fill(vSumAmpFT0CinROF[colIndex], vTracksITS567perColl[colIndex]);

                histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours"))->Fill(vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
              }

              if (vNumCollinROF[colIndex] == 2 && vNumCollinROFinVz10[colIndex] == 2 && vInROFcollIndex[colIndex] == 1) {
                histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvFT0C_vs_occupancyFT0CInROF_2coll"))->Fill(vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);

                if (vNumTracksITS567inROF[colIndex] > vTracksITS567perColl[colIndex]) {
                  histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                  if (vNumTracksITS567inROF[colIndex] > 0)
                    histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"))->Fill(vNumTracksITS567inROF[colIndex], 1.0 * vTracksITS567perColl[colIndex] / vNumTracksITS567inROF[colIndex]);
                } else {
                  histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins"))->Fill(vTracksITS567perColl[colIndex], vNumTracksITS567inROF[colIndex]);
                  if (vTracksITS567perColl[colIndex] > 0)
                    histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInROF_2coll_XaxisWins_RatioV2toV1"))->Fill(vTracksITS567perColl[colIndex], 1.0 * vNumTracksITS567inROF[colIndex] / vTracksITS567perColl[colIndex]);
                }
              }
              // compare with previous ROF
              if (colIndex - 1 >= 0) {
                if (vNumCollinROF[colIndex] == 1 && vNumCollinROFinVz10[colIndex] == 1 && vInROFcollIndex[colIndex] == 0 && vNumCollinROF[colIndex - 1] == 1 && vNumCollinROFinVz10[colIndex - 1] == 1 && vInROFcollIndex[colIndex - 1] == 0) {
                  histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInAnotherEarlierROF_1collPerROF"))->Fill(vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                  if (vROFidThisColl[colIndex] == vROFidThisColl[colIndex - 1] + 1) // one ROF right after the previous
                  {
                    histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF"))->Fill(vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);

                    if (vTracksITS567perColl[colIndex - 1] > vTracksITS567perColl[colIndex]) {
                      histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"))->Fill(vTracksITS567perColl[colIndex - 1], vTracksITS567perColl[colIndex]);
                      if (vTracksITS567perColl[colIndex - 1] > 0)
                        histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"))->Fill(vTracksITS567perColl[colIndex - 1], 1.0 * vTracksITS567perColl[colIndex] / vTracksITS567perColl[colIndex - 1]);
                    } else {
                      histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins"))->Fill(vTracksITS567perColl[colIndex], vTracksITS567perColl[colIndex - 1]);
                      if (vTracksITS567perColl[colIndex] > 0)
                        histos.get<TH2>(HIST("afterNarrowDeltaTimeCut/hThisEvITSTr_vs_occupancyInPreviousROF_1collPerROF_XaxisWins_RatioV2toV1"))->Fill(vTracksITS567perColl[colIndex], 1.0 * vTracksITS567perColl[colIndex - 1] / vTracksITS567perColl[colIndex]);
                    }
                  }
                }
              }
            }
            if (vNoCollInTimeRangeStrict[colIndex]) {
              histos.get<TH2>(HIST("afterStrictDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterStrictDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("afterStrictDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterStrictDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }
            if (vNoHighMultCollInTimeRange[colIndex]) {
              histos.get<TH2>(HIST("afterStandardDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterStandardDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("afterStandardDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterStandardDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInVzDependentTimeRange[colIndex]) {
              histos.get<TH2>(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("afterVzDependentDeltaTimeCut/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInSameRofStrict[colIndex]) {
              histos.get<TH2>(HIST("kNoCollInRofStrict/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofStrict/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofStrict/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofStrict/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInSameRofStandard[colIndex]) {
              histos.get<TH2>(HIST("kNoCollInRofStandard/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofStandard/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofStandard/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofStandard/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }

            if (vNoCollInSameRofWithCloseVz[colIndex]) {
              histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);

              if (vNumTracksITS567inROF[colIndex] > 0) {
                histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyInROF_HasNeighbours"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyFT0CInROF_HasNeighbours"))->Fill(vSumAmpFT0CinROF[colIndex], vTracksITS567perColl[colIndex]);

                histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvFT0C_vs_occupancyFT0CInROF_HasNeighbours"))->Fill(vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
              }

              if (vNumCollinROF[colIndex] == 2 && vNumCollinROFinVz10[colIndex] == 2 && vInROFcollIndex[colIndex] == 1) {
                histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvITSTr_vs_occupancyInROF_2coll"))->Fill(vNumTracksITS567inROF[colIndex], vTracksITS567perColl[colIndex]);
                histos.get<TH2>(HIST("kNoCollInRofWithCloseVz/hThisEvFT0C_vs_occupancyFT0CInROF_2coll"))->Fill(vSumAmpFT0CinROF[colIndex], vAmpFT0CperColl[colIndex]);
              }
            }

            if (vNoCollInTimeRangeNarrow[colIndex] && vNoHighMultCollInTimeRange[colIndex] && vNoCollInSameRofStandard[colIndex] && vNoCollInSameRofWithCloseVz[colIndex]) {
              histos.get<TH2>(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTPCTr_vs_occupancyByFT0C"))->Fill(vSumAmpFT0CinFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
              histos.get<TH2>(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITS567perColl[colIndex]);
              histos.get<TH2>(HIST("NarrowDeltaCut_StdTimeAndRofCuts/hThisEvITSTPCTr_vs_occupancyByTracks"))->Fill(vNumTracksITS567inFullTimeWin[colIndex], vTracksITSTPCperColl[colIndex]);
            }
          }
        }
      }

      // int bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;

      // evsel(alias, selection, sel7, sel8, foundBC, foundFT0, foundFV0, foundFDD, foundZDC, bcInTF,
      //       vNumTracksITS567inFullTimeWin[colIndex], vSumAmpFT0CinFullTimeWin[colIndex], vNumTracksITS567inROF[colIndex]);
    }
  }

  PROCESS_SWITCH(RofOccupancyQaTask, processRun3, "Process Run3 event selection", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);

  return WorkflowSpec{
    adaptAnalysisTask<RofOccupancyQaTask>(cfgc)};
}
