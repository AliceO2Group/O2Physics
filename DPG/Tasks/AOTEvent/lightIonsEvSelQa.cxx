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

/// \file detectorOccupancyQa.cxx
/// \brief Occupancy QA task
///
/// \author Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
// #include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/BunchFilling.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsFT0/RecPoints.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3.h"

#include <map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>; //, aod::CentFT0Cs>;
// using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;

struct LightIonsEvSelQa {
  Configurable<int> confCutMinTPCcls{"MinNumTPCcls", 50, "min number of TPC clusters for a current event"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nBinsTracks{"nBinsTracks", 450, "N bins in n tracks histo"};                            // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nMaxTracks{"nMaxTracks", 450, "N max in n tracks histo"};                               // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nMaxGlobalTracks{"nMaxGlobalTracks", 450, "N max in n tracks histo"};                   // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nBinsMultFwd{"nBinsMultFwd", 1000, "N bins in mult fwd histo"};                         // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> nMaxMultFwd{"nMaxMultFwd", 100000, "N max in mult fwd histo"};                        // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> timeBinWidthInSec{"TimeBinWidthInSec", 10, "Width of time bins in seconds"};          // o2-linter: disable=name/configurable (temporary fix)

  Configurable<float> nSigmaForVzDiff{"nSigmaForVzDiff", 3.0, "n +/- sigma for diff vZ"};    // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> safetyDiffMargin{"SafetyDiffVzMargin", 0.4, "margin for diff vZ, cm"}; // o2-linter: disable=name/configurable (temporary fix)

  uint64_t minGlobalBC = 0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  bool* applySelection = NULL;
  int nBCsPerOrbit = 3564;
  int lastRunNumber = -1;
  double maxSec = 1;
  double minSec = 0;
  int64_t bcSOR = 0;                     // global bc of the start of the first orbit, setting 0 by default for unanchored MC
  int64_t nBCsPerTF = 32 * nBCsPerOrbit; // duration of TF in bcs, should be 128*3564 or 32*3564, setting 128 orbits by default sfor unanchored MC

  // save time "slices" for several collisions for QA
  bool flagFillQAtimeOccupHist = false;
  int nCollisionsForTimeBinQA = 40;
  int counterQAtimeOccupHistos = 0;

  void init(InitContext&)
  {
    // ccdb->setURL("http://ccdb-test.cern.ch:8080");
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // const AxisSpec axisBCinTF{static_cast<int>(nBCsPerTF), 0, static_cast<double>(nBCsPerTF), "bc in TF"};
    // histos.add("hNcolVsBcInTF/hNcolVsBcInTF", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
    // histos.add("hNcolVsBcInTF/hNcolVsBcInTF_vertexTOFmatched", ";bc in TF; n collisions", kTH1F, {axisBCinTF});

    // ##############
    const AxisSpec axisBCs{nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit), ""};

    histos.add("bcQA/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("bcQA/pastActivity/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("bcQA/futureActivity/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("bcQA/noPastActivity/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("bcQA/noFutureActivity/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("bcQA/noPastFutureActivity/hBcFV0", "", kTH1F, {axisBCs});

    histos.add("bcQA/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("bcQA/pastActivity/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("bcQA/futureActivity/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("bcQA/noPastActivity/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("bcQA/noFutureActivity/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("bcQA/noPastFutureActivity/hBcFT0", "", kTH1F, {axisBCs});

    histos.add("bcQA/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("bcQA/pastActivity/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("bcQA/futureActivity/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("bcQA/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("bcQA/pastActivity/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("bcQA/futureActivity/hBcZDC", "", kTH1F, {axisBCs});

    const AxisSpec axisNtracks{nBinsTracks, -0.5, nMaxTracks - 0.5, "n tracks"};
    const AxisSpec axisNtracksGlobal{nBinsTracks, -0.5, nMaxGlobalTracks - 0.5, "n tracks"};
    const AxisSpec axisMultV0A{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd), "mult V0A"};
    const AxisSpec axisMultFT0A{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd * 0.4), "mult FT0C"};
    const AxisSpec axisMultFT0C{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd * 0.15), "mult FT0C"};
    const AxisSpec axisMultT0M{nBinsMultFwd * 2, 0., static_cast<float>(nMaxMultFwd * 0.4), "mult FT0M"};

    const AxisSpec axisVtxZ{800, -20., 20., ""};
    const AxisSpec axisBcDiff{600, -300., 300., "bc difference"};

    const AxisSpec axisNcontrib{801, -0.5, 800.5, "n contributors"};
    const AxisSpec axisColTimeRes{1500, 0., 1500., "collision time resolution (ns)"};

    histos.add("noSpecSelections/hBcColNoSel8", "", kTH1F, {axisBCs});
    // histos.add("noSpecSelections/hBcColNoSel8TOF", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noSpecSelections/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noSpecSelections/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noSpecSelections/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noSpecSelections/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noSpecSelections/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noSpecSelections/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noSpecSelections/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noSpecSelections/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noSpecSelections/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noSpecSelections/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("noPileup/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPileup/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPileup/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPileup/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPileup/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPileup/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPileup/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPileup/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPileup/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPileup/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPileup/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPileup/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPileup/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPileup/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPileup/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noPileup/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPileup/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("pvTOFmatched/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("pvTOFmatched/hBcColNoSel8TOF", "", kTH1F, {axisBCs});
    histos.add("pvTOFmatched/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("pvTOFmatched/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("pvTOFmatched/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("pvTOFmatched/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("pvTOFmatched/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("pvTOFmatched/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("pvTOFmatched/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("pvTOFmatched/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("pvTOFmatched/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("pvTOFmatched/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("pvTOFmatched/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("pvTOFmatched/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("pvTOFmatched/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("pvTOFmatched/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("pvTOFmatched/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("pvTOFmatched/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("bcDiffCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("bcDiffCut/hBcColNoSel8TOF", "", kTH1F, {axisBCs});
    histos.add("bcDiffCut/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("bcDiffCut/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("bcDiffCut/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("bcDiffCut/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("bcDiffCut/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("bcDiffCut/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("bcDiffCut/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("bcDiffCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("bcDiffCut/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("bcDiffCut/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("bcDiffCut/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("bcDiffCut/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("bcDiffCut/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("bcDiffCut/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("bcDiffCut/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("bcDiffCut/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("narrowTimeVeto/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("narrowTimeVeto/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("narrowTimeVeto/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("narrowTimeVeto/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("narrowTimeVeto/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("narrowTimeVeto/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("narrowTimeVeto/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("narrowTimeVeto/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("narrowTimeVeto/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("narrowTimeVeto/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("narrowTimeVeto/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("noCollSameROF/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noCollSameROF/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noCollSameROF/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noCollSameROF/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noCollSameROF/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noCollSameROF/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noCollSameROF/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noCollSameROF/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noCollSameROF/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noCollSameROF/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noCollSameROF/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("noPileup_LowMultCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPileup_LowMultCut/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPileup_LowMultCut/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPileup_LowMultCut/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPileup_LowMultCut/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPileup_LowMultCut/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPileup_LowMultCut/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPileup_LowMultCut/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPileup_LowMultCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPileup_LowMultCut/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPileup_LowMultCut/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPileup_LowMultCut/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPileup_LowMultCut/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPileup_LowMultCut/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPileup_LowMultCut/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noPileup_LowMultCut/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPileup_LowMultCut/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("noPileup_HighMultCloudCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPileup_HighMultCloudCut/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPileup_HighMultCloudCut/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPileup_HighMultCloudCut/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPileup_HighMultCloudCut/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPileup_HighMultCloudCut/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPileup_HighMultCloudCut/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPileup_HighMultCloudCut/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPileup_HighMultCloudCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPileup_HighMultCloudCut/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPileup_HighMultCloudCut/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPileup_HighMultCloudCut/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPileup_HighMultCloudCut/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPileup_HighMultCloudCut/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPileup_HighMultCloudCut/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noPileup_HighMultCloudCut/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPileup_HighMultCloudCut/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("badVzDiff/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("badVzDiff/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("badVzDiff/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("badVzDiff/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("badVzDiff/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("badVzDiff/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("badVzDiff/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("badVzDiff/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("badVzDiff/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("badVzDiff/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("badVzDiff/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("badVzDiff/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("badVzDiff/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("badVzDiff/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("badVzDiff/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("badVzDiff/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("badVzDiff/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("goodVzDiff/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("goodVzDiff/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("goodVzDiff/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("goodVzDiff/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("goodVzDiff/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("goodVzDiff/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("goodVzDiff/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("goodVzDiff/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("goodVzDiff/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("goodVzDiff/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("goodVzDiff/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("goodVzDiff/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("goodVzDiff/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("goodVzDiff/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("goodVzDiff/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("goodVzDiff/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("goodVzDiff/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("noPastActivity/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPastActivity/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPastActivity/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPastActivity/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPastActivity/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPastActivity/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPastActivity/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPastActivity/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPastActivity/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPastActivity/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPastActivity/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPastActivity/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPastActivity/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPastActivity/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPastActivity/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noPastActivity/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPastActivity/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("noFT0activityNearby/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noFT0activityNearby/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noFT0activityNearby/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noFT0activityNearby/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noFT0activityNearby/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noFT0activityNearby/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noFT0activityNearby/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noFT0activityNearby/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noFT0activityNearby/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noFT0activityNearby/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noFT0activityNearby/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noFT0activityNearby/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noFT0activityNearby/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noFT0activityNearby/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noFT0activityNearby/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noFT0activityNearby/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noFT0activityNearby/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hTVXvsBcDiff", "", kTH1F, {axisBcDiff});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPileup_cutByVzDiff_pvTOF_noFT0act/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
  }

  Preslice<FullTracksIU> perCollision = aod::track::collisionId;

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

  // #####
  void processRun3(
    ColEvSels const& cols,
    FullTracksIU const& tracks,
    BCsRun3 const& bcs,
    aod::FT0s const&)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber; // do it only once

      int64_t tsSOR = 0; // dummy start-of-run timestamp
      int64_t tsEOR = 1; // dummy end-of-run timestamp

      if (runNumber >= 500000) {
        auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), runNumber);
        // first bc of the first orbit
        bcSOR = runInfo.orbitSOR * o2::constants::lhc::LHCMaxBunches;
        // duration of TF in bcs
        nBCsPerTF = runInfo.orbitsPerTF * o2::constants::lhc::LHCMaxBunches;

        // start-of-run timestamp
        tsSOR = runInfo.sor;
        // end-of-run timestamp
        tsEOR = runInfo.eor;

        LOGP(info, "bcSOR = {}, nBCsPerTF = {}", bcSOR, nBCsPerTF);
      }

      minSec = floor(tsSOR / 1000.);
      maxSec = ceil(tsEOR / 1000.);
      int nTimeBins = static_cast<int>((maxSec - minSec) / timeBinWidthInSec);
      double timeInterval = nTimeBins * timeBinWidthInSec;

      const AxisSpec axisSeconds{nTimeBins, 0, timeInterval, "seconds"};
      histos.add("hSecondsCollisions/noPU", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisions/noPU_underLine", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisions/noPU_grassOnTheRight", "", kTH1D, {axisSeconds});
      histos.add("hSecondsCollisions/noPU_good", "", kTH1D, {axisSeconds});

    } // end of runNumber check

    // vectors of TVX flags used for past-future studies
    int nBCs = bcs.size();
    std::vector<bool> vIsTVX(nBCs, 0);
    std::vector<uint64_t> vGlobalBCs(nBCs, 0);

    std::vector<bool> vPastActivity(nBCs, 0);
    std::vector<bool> vFutureActivity(nBCs, 0);
    std::vector<bool> vNearbyFT0activity(nBCs, 0);

    // create maps from globalBC to bc index for TVX or FT0-OR fired bcs
    // to be used for closest TVX (FT0-OR) searches
    std::map<int64_t, int32_t> mapGlobalBcWithTVX;
    std::map<int64_t, int32_t> mapGlobalBcWithTOR;

    // ### BC loop
    for (const auto& bc : bcs) {
      uint64_t globalBC = bc.globalBC();

      int indexBc = bc.globalIndex();

      if (bc.selection_bit(kIsBBT0A) || bc.selection_bit(kIsBBT0C)) {
        mapGlobalBcWithTOR[globalBC] = indexBc;
      }
      if (bc.selection_bit(kIsTriggerTVX)) {
        mapGlobalBcWithTVX[globalBC] = indexBc;
      }

      // fill TVX flags for past-future searches
      vIsTVX[indexBc] = bc.selection_bit(kIsTriggerTVX);
      vGlobalBCs[indexBc] = globalBC;

      // ### checking nearby activities
      int deltaIndex = 0;  // backward move counts
      int deltaBC = 0;     // current difference wrt globalBC
      int maxDeltaBC = 30; // maximum difference

      bool nearbyFT0activity = 0;

      // past
      bool pastActivityFT0 = 0;
      bool pastActivityFDD = 0;
      bool pastActivityFV0 = 0;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() - deltaIndex < 0) {
          break;
        }
        const auto& bcPast = bcs.iteratorAt(bc.globalIndex() - deltaIndex);
        deltaBC = globalBC - bcPast.globalBC();
        if (deltaBC < maxDeltaBC) {
          pastActivityFT0 |= bcPast.has_ft0();
          pastActivityFV0 |= bcPast.has_fv0a();
          pastActivityFDD |= bcPast.has_fdd();
        }
        if (deltaBC < 2) {
          if (bcPast.has_ft0()) {
            std::bitset<8> triggers = bcPast.ft0().triggerMask();
            nearbyFT0activity |= triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideA];
          }
        }
      }
      bool pastActivity = pastActivityFT0 | pastActivityFV0 | pastActivityFDD;
      vPastActivity[indexBc] = pastActivity;

      // future
      deltaIndex = 0;
      deltaBC = 0;
      bool futureActivityFT0 = 0;
      bool futureActivityFDD = 0;
      bool futureActivityFV0 = 0;
      while (deltaBC < maxDeltaBC) {
        deltaIndex++;
        if (bc.globalIndex() + deltaIndex >= bcs.size()) {
          break;
        }
        const auto& bcFuture = bcs.iteratorAt(bc.globalIndex() + deltaIndex);
        deltaBC = bcFuture.globalBC() - globalBC;
        if (deltaBC < maxDeltaBC) {
          futureActivityFT0 |= bcFuture.has_ft0();
          futureActivityFV0 |= bcFuture.has_fv0a();
          futureActivityFDD |= bcFuture.has_fdd();
        }
        if (deltaBC < 2) {
          if (bcFuture.has_ft0()) {
            std::bitset<8> triggers = bcFuture.ft0().triggerMask();
            nearbyFT0activity |= triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideA];
          }
        }
      }
      bool futureActivity = futureActivityFT0 | futureActivityFV0 | futureActivityFDD;
      vFutureActivity[indexBc] = futureActivity;
      vNearbyFT0activity[indexBc] = nearbyFT0activity;

      // monitor BCs with nearby activity:

      int localBC = globalBC % nBCsPerOrbit;
      if (bc.has_fv0a()) {
        histos.fill(HIST("bcQA/hBcFV0"), localBC);
        if (pastActivity)
          histos.fill(HIST("bcQA/pastActivity/hBcFV0"), localBC);
        if (futureActivity)
          histos.fill(HIST("bcQA/futureActivity/hBcFV0"), localBC);
        if (!pastActivity)
          histos.fill(HIST("bcQA/noPastActivity/hBcFV0"), localBC);
        if (!futureActivity)
          histos.fill(HIST("bcQA/noFutureActivity/hBcFV0"), localBC);
        if (!pastActivity && !futureActivity)
          histos.fill(HIST("bcQA/noPastFutureActivity/hBcFV0"), localBC);
      }
      if (bc.has_ft0()) {
        histos.fill(HIST("bcQA/hBcFT0"), localBC);
        if (pastActivity)
          histos.fill(HIST("bcQA/pastActivity/hBcFT0"), localBC);
        if (futureActivity)
          histos.fill(HIST("bcQA/futureActivity/hBcFT0"), localBC);
        if (!pastActivity)
          histos.fill(HIST("bcQA/noPastActivity/hBcFT0"), localBC);
        if (!futureActivity)
          histos.fill(HIST("bcQA/noFutureActivity/hBcFT0"), localBC);
        if (!pastActivity && !futureActivity)
          histos.fill(HIST("bcQA/noPastFutureActivity/hBcFT0"), localBC);
      }
      if (bc.has_fdd()) {
        histos.fill(HIST("bcQA/hBcFDD"), localBC);
        if (pastActivity)
          histos.fill(HIST("bcQA/pastActivity/hBcFDD"), localBC);
        if (futureActivity)
          histos.fill(HIST("bcQA/futureActivity/hBcFDD"), localBC);
      }
      if (bc.has_zdc()) {
        histos.fill(HIST("bcQA/hBcZDC"), localBC);
        if (pastActivity)
          histos.fill(HIST("bcQA/pastActivity/hBcZDC"), localBC);
        if (futureActivity)
          histos.fill(HIST("bcQA/futureActivity/hBcZDC"), localBC);
      }

    } // end of bc loop

    // ### collision loop
    for (const auto& col : cols) {
      if (std::abs(col.posZ()) > 10)
        continue;

      const auto& foundBC = col.foundBC_as<BCsRun3>();
      uint64_t globalBC = foundBC.globalBC();
      // uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;

      int64_t ts = foundBC.timestamp();
      double secFromSOR = ts / 1000. - minSec;

      // search for nearest ft0a&ft0c entry
      uint64_t globalOrigBC = col.bc_as<BCsRun3>().globalBC();
      int32_t indexClosestTVX = findClosest(globalOrigBC, mapGlobalBcWithTVX);
      int bcToClosestTVXdiff = static_cast<int>(globalOrigBC - vGlobalBCs[indexClosestTVX]);

      // selection decisions:
      bool noPU = col.selection_bit(kNoSameBunchPileup);
      bool pvTOFmatched = col.selection_bit(kIsVertexTOFmatched);
      bool narrowDeltaTimeVeto = col.selection_bit(kNoCollInTimeRangeNarrow);
      bool noCollSameROF = col.selection_bit(kNoCollInRofStrict);
      bool bcDiffCut = (bcToClosestTVXdiff == 0);

      auto bcIndex = foundBC.globalIndex();
      // bool noNearbyActivity = (vPastActivity[bcIndex] == 0 && vFutureActivity[bcIndex] == 0);
      bool noPastActivity = (vPastActivity[bcIndex] == 0);

      // ### count tracks of different types
      int nPVtracks = 0;
      int nGlobalTracks = 0;
      // int nTOFtracks = 0;
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (const auto& track : tracksGrouped) {
        if (!track.isPVContributor()) {
          continue;
        }
        if (track.itsNCls() < 5)
          continue;
        if (track.pt() < 0.2 || track.pt() > 10)
          continue;
        if (std::abs(track.eta()) > 0.8)
          continue;

        nPVtracks++;
        // nTOFtracks += track.hasTOF();

        if (track.hasITS() && track.hasTPC() && track.tpcNClsFound() > 50 && track.tpcNClsCrossedRows() > 50 && track.tpcChi2NCl() < 4)
          nGlobalTracks++;
      } // end of track loop

      bool hasFT0 = foundBC.has_ft0();
      bool hasFV0A = foundBC.has_fv0a();

      // bool noFT0activityNearby = false;
      bool noFT0activityNearby = (vNearbyFT0activity[bcIndex] == 0);
      // check kIsFlangeEvent
      if (hasFT0) {
        std::bitset<8> triggers = foundBC.ft0().triggerMask();
        if (triggers[o2::ft0::RecPoints::ETriggerBits::kIsFlangeEvent])
          noFT0activityNearby = false;
      }

      float vZ = col.posZ();
      float vZft0 = hasFT0 ? foundBC.ft0().posZ() : -1000;
      float diffVz = vZft0 - vZ;

      float multV0A = hasFV0A ? col.multFV0A() : 0;
      float multT0A = hasFT0 ? col.multFT0A() : 0;
      float multT0C = hasFT0 ? col.multFT0C() : 0;

      float multT0M = multT0A + multT0C;
      // bool badVzDiff = col.selection_bit(kIsGoodZvtxFT0vsPV);
      // bool badVzDiff = hasFT0 && (multT0M > 5000) && (diffVz < -2.5 || diffVz > 2.5);
      float meanDiff = 0.0; // cm
      // O-O
      if (lastRunNumber == 564356)
        meanDiff = -0.01;
      if (lastRunNumber == 564359)
        meanDiff = -0.17;
      if (lastRunNumber == 564373)
        meanDiff = 0.99;
      if (lastRunNumber == 564374)
        meanDiff = 0.57;
      // Ne-Ne
      if (lastRunNumber == 564468)
        meanDiff = -0.51;
      if (lastRunNumber == 564468)
        meanDiff = -0.60;

      // float stdDev = (multT0M > 20) ? 0.54 + 6.46 / sqrt(multT0M) : 2.0; // cm
      float stdDev = (multT0M > 10) ? 0.144723 + 13.5345 / sqrt(multT0M) : 1.5; // cm
      if (multT0M > 4000)
        stdDev = 0.35; // cm
      bool badVzDiff = diffVz < (meanDiff - stdDev * nSigmaForVzDiff - safetyDiffMargin) || diffVz > (meanDiff + stdDev * nSigmaForVzDiff + safetyDiffMargin);

      bool underLine = false;
      if (hasFT0 && nPVtracks < 45. / 40000 * multV0A - 4 && nPVtracks < 20) {
        underLine = true;
      }
      bool grassOnTheRight = false;
      if (hasFT0 && nPVtracks < 220. / 40000 * multV0A - 100 && nPVtracks >= 25) {
        grassOnTheRight = true;
      }

      // vs time
      if (noPU) {
        histos.fill(HIST("hSecondsCollisions/noPU"), secFromSOR);
        if (underLine)
          histos.fill(HIST("hSecondsCollisions/noPU_underLine"), secFromSOR);
        if (grassOnTheRight)
          histos.fill(HIST("hSecondsCollisions/noPU_grassOnTheRight"), secFromSOR);
        if (!underLine && !grassOnTheRight)
          histos.fill(HIST("hSecondsCollisions/noPU_good"), secFromSOR);
      }

      // study bc diff:
      int nContributors = col.numContrib();
      float timeRes = col.collisionTimeRes();
      // int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      // if (col.selection_bit(kIsVertexTOFmatched)) {
      //   histos.fill(HIST("hColBcDiffVsNcontribWithTOF"), nContributors, bcToClosestTVXdiff);
      //   histos.fill(HIST("hColTimeResVsNcontribWithTOF"), nContributors, timeRes);
      //   histos.fill(HIST("hBcColTOF"), localBC);
      // }

      histos.fill(HIST("noSpecSelections/hBcColNoSel8"), localBC);

      if (noPU) {
        histos.fill(HIST("noPileup/hBcColNoSel8"), localBC);
      }
      if (pvTOFmatched) {
        histos.fill(HIST("pvTOFmatched/hBcColNoSel8"), localBC);
      }
      if (bcDiffCut) {
        histos.fill(HIST("bcDiffCut/hBcColNoSel8"), localBC);
      }
      if (noPastActivity) {
        histos.fill(HIST("noPastActivity/hBcColNoSel8"), localBC);
      }
      if (noFT0activityNearby) {
        histos.fill(HIST("noFT0activityNearby/hBcColNoSel8"), localBC);
      }
      if (badVzDiff) {
        histos.fill(HIST("badVzDiff/hBcColNoSel8"), localBC);
      }
      if (!badVzDiff) {
        histos.fill(HIST("goodVzDiff/hBcColNoSel8"), localBC);
      }
      if (narrowDeltaTimeVeto) {
        histos.fill(HIST("narrowTimeVeto/hBcColNoSel8"), localBC);
      }
      if (noCollSameROF) {
        histos.fill(HIST("noCollSameROF/hBcColNoSel8"), localBC);
      }
      if (noPU && underLine) {
        histos.fill(HIST("noPileup_LowMultCut/hBcColNoSel8"), localBC);
      }
      if (noPU && grassOnTheRight) {
        histos.fill(HIST("noPileup_HighMultCloudCut/hBcColNoSel8"), localBC);
      }
      if (noPU && pvTOFmatched && !badVzDiff && noFT0activityNearby) { // noPileup_cutByVzDiff_pvTOF_noFT0act
        histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcColNoSel8"), localBC);
      }

      // only here cut on sel8:
      if (!col.sel8())
        continue;

      histos.fill(HIST("noSpecSelections/hBcTVX"), localBC);
      histos.fill(HIST("noSpecSelections/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
      histos.fill(HIST("noSpecSelections/hColTimeResVsNcontrib"), nContributors, timeRes);

      if (noPU) {
        histos.fill(HIST("noPileup/hBcTVX"), localBC);
        histos.fill(HIST("noPileup/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPileup/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (pvTOFmatched) {
        histos.fill(HIST("pvTOFmatched/hBcTVX"), localBC);
        histos.fill(HIST("pvTOFmatched/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("pvTOFmatched/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (bcDiffCut) {
        histos.fill(HIST("bcDiffCut/hBcTVX"), localBC);
        histos.fill(HIST("bcDiffCut/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("bcDiffCut/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (noPastActivity) {
        histos.fill(HIST("noPastActivity/hBcTVX"), localBC);
        histos.fill(HIST("noPastActivity/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPastActivity/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (noFT0activityNearby) {
        histos.fill(HIST("noFT0activityNearby/hBcTVX"), localBC);
        histos.fill(HIST("noFT0activityNearby/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noFT0activityNearby/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (badVzDiff) {
        histos.fill(HIST("badVzDiff/hBcTVX"), localBC);
        histos.fill(HIST("badVzDiff/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("badVzDiff/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (!badVzDiff) {
        histos.fill(HIST("goodVzDiff/hBcTVX"), localBC);
        histos.fill(HIST("goodVzDiff/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("goodVzDiff/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (narrowDeltaTimeVeto) {
        histos.fill(HIST("narrowTimeVeto/hBcTVX"), localBC);
        histos.fill(HIST("narrowTimeVeto/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("narrowTimeVeto/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (noCollSameROF) {
        histos.fill(HIST("noCollSameROF/hBcTVX"), localBC);
        histos.fill(HIST("noCollSameROF/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noCollSameROF/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (noPU && underLine) {
        histos.fill(HIST("noPileup_LowMultCut/hBcTVX"), localBC);
        histos.fill(HIST("noPileup_LowMultCut/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPileup_LowMultCut/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (noPU && grassOnTheRight) {
        histos.fill(HIST("noPileup_HighMultCloudCut/hBcTVX"), localBC);
        histos.fill(HIST("noPileup_HighMultCloudCut/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPileup_HighMultCloudCut/hColTimeResVsNcontrib"), nContributors, timeRes);
      }
      if (noPU && pvTOFmatched && !badVzDiff && noFT0activityNearby) { // noPileup_cutByVzDiff_pvTOF_noFT0act
        histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcTVX"), localBC);
        histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hColTimeResVsNcontrib"), nContributors, timeRes);
      }

      if (foundBC.has_ft0()) {
        // float multT0A = foundBC.has_ft0() ? foundBC.ft0().sumAmpA() : -999.f;
        // float multT0C = foundBC.has_ft0() ? foundBC.ft0().sumAmpC() : -999.f;

        histos.fill(HIST("noSpecSelections/hBcFT0"), localBC);
        histos.fill(HIST("noSpecSelections/hVtxFT0VsVtxCol"), vZft0, vZ);
        histos.fill(HIST("noSpecSelections/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
        histos.fill(HIST("noSpecSelections/nTracksPV_vs_T0A"), multT0A, nPVtracks);
        histos.fill(HIST("noSpecSelections/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
        histos.fill(HIST("noSpecSelections/nTracksPV_vs_T0C"), multT0C, nPVtracks);
        histos.fill(HIST("noSpecSelections/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);

        if (noPU) {
          histos.fill(HIST("noPileup/hBcFT0"), localBC);
          histos.fill(HIST("noPileup/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPileup/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPileup/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPileup/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("noPileup/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPileup/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (pvTOFmatched) {
          histos.fill(HIST("pvTOFmatched/hBcFT0"), localBC);
          histos.fill(HIST("pvTOFmatched/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("pvTOFmatched/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("pvTOFmatched/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("pvTOFmatched/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("pvTOFmatched/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("pvTOFmatched/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (bcDiffCut) {
          histos.fill(HIST("bcDiffCut/hBcFT0"), localBC);
          histos.fill(HIST("bcDiffCut/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("bcDiffCut/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("bcDiffCut/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("bcDiffCut/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("bcDiffCut/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("bcDiffCut/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (badVzDiff) {
          histos.fill(HIST("badVzDiff/hBcFT0"), localBC);
          histos.fill(HIST("badVzDiff/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("badVzDiff/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("badVzDiff/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("badVzDiff/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("badVzDiff/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("badVzDiff/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (!badVzDiff) {
          histos.fill(HIST("goodVzDiff/hBcFT0"), localBC);
          histos.fill(HIST("goodVzDiff/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("goodVzDiff/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("goodVzDiff/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("goodVzDiff/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("goodVzDiff/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("goodVzDiff/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (noPastActivity) {
          histos.fill(HIST("noPastActivity/hBcFT0"), localBC);
          histos.fill(HIST("noPastActivity/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPastActivity/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPastActivity/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPastActivity/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("noPastActivity/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPastActivity/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (noFT0activityNearby) {
          histos.fill(HIST("noFT0activityNearby/hBcFT0"), localBC);
          histos.fill(HIST("noFT0activityNearby/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noFT0activityNearby/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noFT0activityNearby/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noFT0activityNearby/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("noFT0activityNearby/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noFT0activityNearby/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (narrowDeltaTimeVeto) {
          histos.fill(HIST("narrowTimeVeto/hBcFT0"), localBC);
          histos.fill(HIST("narrowTimeVeto/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("narrowTimeVeto/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("narrowTimeVeto/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("narrowTimeVeto/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("narrowTimeVeto/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("narrowTimeVeto/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (noCollSameROF) {
          histos.fill(HIST("noCollSameROF/hBcFT0"), localBC);
          histos.fill(HIST("noCollSameROF/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noCollSameROF/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noCollSameROF/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noCollSameROF/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("noCollSameROF/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noCollSameROF/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (noPU && underLine) {
          histos.fill(HIST("noPileup_LowMultCut/hBcFT0"), localBC);
          histos.fill(HIST("noPileup_LowMultCut/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPileup_LowMultCut/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPileup_LowMultCut/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPileup_LowMultCut/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("noPileup_LowMultCut/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPileup_LowMultCut/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (noPU && grassOnTheRight) {
          histos.fill(HIST("noPileup_HighMultCloudCut/hBcFT0"), localBC);
          histos.fill(HIST("noPileup_HighMultCloudCut/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPileup_HighMultCloudCut/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPileup_HighMultCloudCut/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPileup_HighMultCloudCut/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("noPileup_HighMultCloudCut/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPileup_HighMultCloudCut/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
        if (noPU && pvTOFmatched && !badVzDiff && noFT0activityNearby) { // noPileup_cutByVzDiff_pvTOF_noFT0act
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcFT0"), localBC);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracks);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracks);
        }
      }

      if (foundBC.has_fv0a()) {
        histos.fill(HIST("noSpecSelections/hBcFV0"), localBC);
        histos.fill(HIST("noSpecSelections/nTracksPV_vs_V0A"), multV0A, nPVtracks);
        histos.fill(HIST("noSpecSelections/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        if (noPU) {
          histos.fill(HIST("noPileup/hBcFV0"), localBC);
          histos.fill(HIST("noPileup/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPileup/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (pvTOFmatched) {
          histos.fill(HIST("pvTOFmatched/hBcFV0"), localBC);
          histos.fill(HIST("pvTOFmatched/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("pvTOFmatched/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (bcDiffCut) {
          histos.fill(HIST("bcDiffCut/hBcFV0"), localBC);
          histos.fill(HIST("bcDiffCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("bcDiffCut/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (badVzDiff) {
          histos.fill(HIST("badVzDiff/hBcFV0"), localBC);
          histos.fill(HIST("badVzDiff/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("badVzDiff/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (!badVzDiff) {
          histos.fill(HIST("goodVzDiff/hBcFV0"), localBC);
          histos.fill(HIST("goodVzDiff/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("goodVzDiff/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (noPastActivity) {
          histos.fill(HIST("noPastActivity/hBcFV0"), localBC);
          histos.fill(HIST("noPastActivity/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPastActivity/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (noFT0activityNearby) {
          histos.fill(HIST("noFT0activityNearby/hBcFV0"), localBC);
          histos.fill(HIST("noFT0activityNearby/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noFT0activityNearby/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (narrowDeltaTimeVeto) {
          histos.fill(HIST("narrowTimeVeto/hBcFV0"), localBC);
          histos.fill(HIST("narrowTimeVeto/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("narrowTimeVeto/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (noCollSameROF) {
          histos.fill(HIST("noCollSameROF/hBcFV0"), localBC);
          histos.fill(HIST("noCollSameROF/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noCollSameROF/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (noPU && underLine) {
          histos.fill(HIST("noPileup_LowMultCut/hBcFV0"), localBC);
          histos.fill(HIST("noPileup_LowMultCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPileup_LowMultCut/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (noPU && grassOnTheRight) {
          histos.fill(HIST("noPileup_HighMultCloudCut/hBcFV0"), localBC);
          histos.fill(HIST("noPileup_HighMultCloudCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPileup_HighMultCloudCut/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
        if (noPU && pvTOFmatched && !badVzDiff && noFT0activityNearby) { // noPileup_cutByVzDiff_pvTOF_noFT0act
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcFV0"), localBC);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracks);
        }
      }
      if (foundBC.has_zdc()) {
        histos.fill(HIST("noSpecSelections/hBcZDC"), localBC);
        if (noPU) {
          histos.fill(HIST("noPileup/hBcZDC"), localBC);
        }
        if (pvTOFmatched) {
          histos.fill(HIST("pvTOFmatched/hBcZDC"), localBC);
        }
        if (bcDiffCut) {
          histos.fill(HIST("bcDiffCut/hBcZDC"), localBC);
        }
        if (badVzDiff) {
          histos.fill(HIST("badVzDiff/hBcZDC"), localBC);
        }
        if (!badVzDiff) {
          histos.fill(HIST("goodVzDiff/hBcZDC"), localBC);
        }
        if (noPastActivity) {
          histos.fill(HIST("noPastActivity/hBcZDC"), localBC);
        }
        if (noFT0activityNearby) {
          histos.fill(HIST("noFT0activityNearby/hBcZDC"), localBC);
        }
        if (narrowDeltaTimeVeto) {
          histos.fill(HIST("narrowTimeVeto/hBcZDC"), localBC);
        }
        if (noCollSameROF) {
          histos.fill(HIST("noCollSameROF/hBcZDC"), localBC);
        }
        if (noPU && underLine) {
          histos.fill(HIST("noPileup_LowMultCut/hBcZDC"), localBC);
        }
        if (noPU && grassOnTheRight) {
          histos.fill(HIST("noPileup_HighMultCloudCut/hBcZDC"), localBC);
        }
        if (noPU && pvTOFmatched && !badVzDiff && noFT0activityNearby) { // noPileup_cutByVzDiff_pvTOF_noFT0act
          histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hBcZDC"), localBC);
        }
      }

      // bc diff
      // auto bc = col.bc_as<BCsRun3>();
      auto bcOriginal = globalOrigBC % 3564;
      float bcDiff = bcOriginal - localBC;

      histos.fill(HIST("noSpecSelections/hTVXvsBcDiff"), bcDiff);
      if (noPU) {
        histos.fill(HIST("noPileup/hTVXvsBcDiff"), bcDiff);
      }
      if (pvTOFmatched) {
        histos.fill(HIST("pvTOFmatched/hTVXvsBcDiff"), bcDiff);
      }
      if (bcDiffCut) {
        histos.fill(HIST("bcDiffCut/hTVXvsBcDiff"), bcDiff);
      }
      if (badVzDiff) {
        histos.fill(HIST("badVzDiff/hTVXvsBcDiff"), bcDiff);
      }
      if (!badVzDiff) {
        histos.fill(HIST("goodVzDiff/hTVXvsBcDiff"), bcDiff);
      }
      if (noPastActivity) {
        histos.fill(HIST("noPastActivity/hTVXvsBcDiff"), bcDiff);
      }
      if (noFT0activityNearby) {
        histos.fill(HIST("noFT0activityNearby/hTVXvsBcDiff"), bcDiff);
      }
      if (narrowDeltaTimeVeto) {
        histos.fill(HIST("narrowTimeVeto/hTVXvsBcDiff"), bcDiff);
      }
      if (noCollSameROF) {
        histos.fill(HIST("noCollSameROF/hTVXvsBcDiff"), bcDiff);
      }
      if (noPU && underLine) {
        histos.fill(HIST("noPileup_LowMultCut/hTVXvsBcDiff"), bcDiff);
      }
      if (noPU && grassOnTheRight) {
        histos.fill(HIST("noPileup_HighMultCloudCut/hTVXvsBcDiff"), bcDiff);
      }
      if (noPU && pvTOFmatched && !badVzDiff && noFT0activityNearby) { // noPileup_cutByVzDiff_pvTOF_noFT0act
        histos.fill(HIST("noPileup_cutByVzDiff_pvTOF_noFT0act/hTVXvsBcDiff"), bcDiff);
      }

    } // end of collisions loop
  }
  PROCESS_SWITCH(LightIonsEvSelQa, processRun3, "Process Run3 tracking vs detector occupancy QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LightIonsEvSelQa>(cfgc)};
}
