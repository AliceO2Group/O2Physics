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
  Configurable<int> nBinsMultFwd{"nBinsMultFwd", 800, "N bins in mult fwd histo"};                          // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> nMaxMultFwd{"nMaxMultFwd", 100000, "N max in mult fwd histo"};                        // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> timeBinWidthInSec{"TimeBinWidthInSec", 10, "Width of time bins in seconds"};          // o2-linter: disable=name/configurable (temporary fix)

  Configurable<float> nSigmaForVzDiff{"nSigmaForVzDiff", 2.5, "n +/- sigma for diff vZ"};      // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> safetyDiffVzMargin{"SafetyDiffVzMargin", 0.5, "margin for diff vZ, cm"}; // o2-linter: disable=name/configurable (temporary fix)

  Configurable<int> confUseDiffVzCutFromEvSel{"UseDiffVzCutFromEvSel", 0, "0 - custom diffVz cut from this task, 1 - cut from event selection"}; // o2-linter: disable=name/configurable (temporary fix)

  Configurable<bool> isMC{"isMC", false, "Run MC"};

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

    histos.add("bcQA/specFT0bits/hBc_kIsActiveSideA", "", kTH1F, {axisBCs});
    histos.add("bcQA/specFT0bits/hBc_kIsActiveSideC", "", kTH1F, {axisBCs});
    histos.add("bcQA/specFT0bits/hBc_kIsFlangeEvent", "", kTH1F, {axisBCs});

    histos.add("bcQA/specFT0bits/hBc_kIsActiveSideA_inTVX", "", kTH1F, {axisBCs});
    histos.add("bcQA/specFT0bits/hBc_kIsActiveSideC_inTVX", "", kTH1F, {axisBCs});
    histos.add("bcQA/specFT0bits/hBc_kIsFlangeEvent_inTVX", "", kTH1F, {axisBCs});

    const AxisSpec axisNtracks{nBinsTracks, -0.5, nMaxTracks - 0.5, "n tracks"};
    const AxisSpec axisNtracksGlobal{nBinsTracks, -0.5, nMaxGlobalTracks - 0.5, "n tracks"};
    const AxisSpec axisMultV0A{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd), "mult V0A"};
    const AxisSpec axisMultFT0A{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd * 0.4), "mult FT0C"};
    const AxisSpec axisMultFT0C{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd * 0.15), "mult FT0C"};
    const AxisSpec axisMultT0M{nBinsMultFwd * 2, 0., static_cast<float>(nMaxMultFwd * 0.4), "mult FT0M"};

    const AxisSpec axisVtxZ{800, -20., 20., ""};
    const AxisSpec axisBcDiff{601, -300.5, 300.5, "bc difference"};

    const AxisSpec axisNcontrib{601, -0.5, 600.5, "n contributors"};
    const AxisSpec axisColTimeRes{1500, 0., 1500., "collision time resolution (ns)"};

    AxisSpec axisVertexChi2{100, 0, 500, "Chi2 of vertex fit"};
    AxisSpec axisVertexChi2perContrib{100, 0, 10, "Chi2 of vertex fit"};

    const AxisSpec axisTimeZN{800, -20., 20., ""};
    const AxisSpec axisTimeDiff{150, -10., 10., ""};
    const AxisSpec axisTimeSum{150, -10., 10., ""};
    const AxisSpec axisZNampl{200, 0, 10000, ""};

    histos.add("noSpecSelections/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcOrigNoSel8", "", kTH1F, {axisBCs});
    // histos.add("noSpecSelections/hBcColNoSel8TOF", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noSpecSelections/hBcOrig", "", kTH1F, {axisBCs});
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
    histos.add("noSpecSelections/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noSpecSelections/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noSpecSelections/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noSpecSelections/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noSpecSelections/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});
    histos.add("noSpecSelections/hTimeZN_AC_sum_vs_diff", ";ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDiff, axisTimeSum});
    histos.add("noSpecSelections/hTimeZN_A_vs_C", ";ZNA time (ns); ZNC time (ns)", kTH2F, {axisTimeZN, axisTimeZN});
    histos.add("noSpecSelections/hTimeZNA", ";ZNA time (ns)", kTH1F, {axisTimeZN});
    histos.add("noSpecSelections/hTimeZNC", ";ZNC time (ns)", kTH1F, {axisTimeZN});

    histos.add("noPU/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU/hBcOrigNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPU/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});
    histos.add("noPU/hTimeZN_AC_sum_vs_diff", ";ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDiff, axisTimeSum});
    histos.add("noPU/hTimeZN_A_vs_C", ";ZNA time (ns); ZNC time (ns)", kTH2F, {axisTimeZN, axisTimeZN});
    histos.add("noPU/hTimeZNA", ";ZNA time (ns)", kTH1F, {axisTimeZN});
    histos.add("noPU/hTimeZNC", ";ZNC time (ns)", kTH1F, {axisTimeZN});
    histos.add("noPU/hAmplZNAC", "ZNC vs ZNA", kTH2F, {axisZNampl, axisZNampl});

    histos.add("noPU_pvTOFmatched/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTOFmatched/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTOFmatched/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTOFmatched/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTOFmatched/hBcFV0", "", kTH1F, {axisBCs});
    // histos.add("noPU_pvTOFmatched/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTOFmatched/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTOFmatched/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_pvTOFmatched/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_pvTOFmatched/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_pvTOFmatched/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_pvTOFmatched/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_pvTOFmatched/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_pvTOFmatched/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_pvTOFmatched/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_pvTOFmatched/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_pvTOFmatched/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_pvTOFmatched/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_pvTOFmatched/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_pvTOFmatched/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_pvTRDmatched/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTRDmatched/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTRDmatched/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTRDmatched/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTRDmatched/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTRDmatched/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTRDmatched/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_pvTRDmatched/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_pvTRDmatched/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_pvTRDmatched/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_pvTRDmatched/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_pvTRDmatched/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_pvTRDmatched/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_pvTRDmatched/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_pvTRDmatched/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_pvTRDmatched/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_pvTRDmatched/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_pvTRDmatched/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_pvTRDmatched/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_notTRDmatched/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_notTRDmatched/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_notTRDmatched/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_notTRDmatched/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_notTRDmatched/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_notTRDmatched/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_notTRDmatched/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_notTRDmatched/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_notTRDmatched/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_notTRDmatched/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_notTRDmatched/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_notTRDmatched/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_notTRDmatched/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_notTRDmatched/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_notTRDmatched/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_notTRDmatched/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_notTRDmatched/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_notTRDmatched/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_pvTOFmatched_notTRDmatched/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_pvTOFmatched_notTRDmatched/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_pvTOFmatched_notTRDmatched/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_pvTOFmatched_notTRDmatched/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("bcDiffWrtClosestTVXCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("bcDiffWrtClosestTVXCut/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("bcDiffWrtClosestTVXCut/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("bcDiffWrtClosestTVXCut/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("bcDiffWrtClosestTVXCut/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("bcDiffWrtClosestTVXCut/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("bcDiffWrtClosestTVXCut/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("bcDiffWrtClosestTVXCut/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("bcDiffWrtClosestTVXCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("bcDiffWrtClosestTVXCut/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("bcDiffWrtClosestTVXCut/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    // histos.add("bcDiffWrtClosestTVXCut/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    // histos.add("bcDiffWrtClosestTVXCut/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    // histos.add("bcDiffWrtClosestTVXCut/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("bcDiffWrtClosestTVXCut/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("bcDiffWrtClosestTVXCut/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("bcDiffWrtClosestTVXCut/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("bcDiffWrtClosestTVXCut/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_bcDiffWrtOriginalBcCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_bcDiffWrtOriginalBcCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_bcDiffWrtOriginalBcCut/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_bcDiffWrtOriginalBcCut/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    // histos.add("noPU_bcDiffWrtOriginalBcCut/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    // histos.add("noPU_bcDiffWrtOriginalBcCut/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    // histos.add("noPU_bcDiffWrtOriginalBcCut/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_bcDiffWrtOriginalBcCut/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_goodVertexChi2/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVertexChi2/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVertexChi2/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVertexChi2/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVertexChi2/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVertexChi2/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVertexChi2/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_goodVertexChi2/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_goodVertexChi2/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_goodVertexChi2/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_goodVertexChi2/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    // histos.add("noPU_goodVertexChi2/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    // histos.add("noPU_goodVertexChi2/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    // histos.add("noPU_goodVertexChi2/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_goodVertexChi2/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_goodVertexChi2/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_goodVertexChi2/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_goodVertexChi2/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("narrowTimeVeto/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("narrowTimeVeto/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("narrowTimeVeto/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("narrowTimeVeto/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("narrowTimeVeto/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("narrowTimeVeto/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("narrowTimeVeto/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("narrowTimeVeto/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("narrowTimeVeto/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("narrowTimeVeto/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("narrowTimeVeto/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("narrowTimeVeto/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("narrowTimeVeto/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("strictTimeVeto/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("strictTimeVeto/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("strictTimeVeto/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("strictTimeVeto/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("strictTimeVeto/hBcFV0", "", kTH1F, {axisBCs});
    // histos.add("strictTimeVeto/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("strictTimeVeto/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("strictTimeVeto/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("strictTimeVeto/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("strictTimeVeto/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("strictTimeVeto/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("strictTimeVeto/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("strictTimeVeto/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("strictTimeVeto/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("strictTimeVeto/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("strictTimeVeto/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("strictTimeVeto/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("strictTimeVeto/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("strictTimeVeto/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noCollSameROF/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcFV0", "", kTH1F, {axisBCs});
    // histos.add("noCollSameROF/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noCollSameROF/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noCollSameROF/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noCollSameROF/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noCollSameROF/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noCollSameROF/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noCollSameROF/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noCollSameROF/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noCollSameROF/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noCollSameROF/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noCollSameROF/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noCollSameROF/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noCollSameROF/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("lowMultCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("lowMultCut/hBcOrigNoSel8", "", kTH1F, {axisBCs});
    histos.add("lowMultCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});

    histos.add("noPU_lowMultCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_lowMultCut/hBcOrigNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_lowMultCut/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_lowMultCut/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_lowMultCut/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_lowMultCut/hBcFV0", "", kTH1F, {axisBCs});
    // histos.add("noPU_lowMultCut/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPU_lowMultCut/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_lowMultCut/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_lowMultCut/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_lowMultCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_lowMultCut/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_lowMultCut/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_lowMultCut/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_lowMultCut/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_lowMultCut/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_lowMultCut/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_lowMultCut/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_lowMultCut/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_lowMultCut/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_lowMultCut/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});
    histos.add("noPU_lowMultCut/hTimeZN_AC_sum_vs_diff", ";ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDiff, axisTimeSum});
    histos.add("noPU_lowMultCut/hTimeZN_A_vs_C", ";ZNA time (ns); ZNC time (ns)", kTH2F, {axisTimeZN, axisTimeZN});
    histos.add("noPU_lowMultCut/hTimeZNA", ";ZNA time (ns)", kTH1F, {axisTimeZN});
    histos.add("noPU_lowMultCut/hTimeZNC", ";ZNC time (ns)", kTH1F, {axisTimeZN});

    histos.add("highMultCloudCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("highMultCloudCut/hBcOrigNoSel8", "", kTH1F, {axisBCs});
    histos.add("highMultCloudCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});

    histos.add("noPU_highMultCloudCut/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_highMultCloudCut/hBcOrigNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_highMultCloudCut/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_highMultCloudCut/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_highMultCloudCut/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_highMultCloudCut/hBcFV0", "", kTH1F, {axisBCs});
    // histos.add("noPU_highMultCloudCut/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPU_highMultCloudCut/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_highMultCloudCut/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_highMultCloudCut/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_highMultCloudCut/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_highMultCloudCut/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_highMultCloudCut/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_highMultCloudCut/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_highMultCloudCut/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_highMultCloudCut/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_highMultCloudCut/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_highMultCloudCut/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_highMultCloudCut/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_highMultCloudCut/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_highMultCloudCut/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});
    histos.add("noPU_highMultCloudCut/hTimeZN_AC_sum_vs_diff", ";ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDiff, axisTimeSum});
    histos.add("noPU_highMultCloudCut/hTimeZN_A_vs_C", ";ZNA time (ns); ZNC time (ns)", kTH2F, {axisTimeZN, axisTimeZN});
    histos.add("noPU_highMultCloudCut/hTimeZNA", ";ZNA time (ns)", kTH1F, {axisTimeZN});
    histos.add("noPU_highMultCloudCut/hTimeZNC", ";ZNC time (ns)", kTH1F, {axisTimeZN});

    histos.add("noPU_badVzDiff/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_badVzDiff/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_badVzDiff/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_badVzDiff/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_badVzDiff/hBcFV0", "", kTH1F, {axisBCs});
    // histos.add("noPU_badVzDiff/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPU_badVzDiff/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_badVzDiff/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_badVzDiff/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_badVzDiff/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_badVzDiff/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_badVzDiff/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_badVzDiff/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_badVzDiff/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_badVzDiff/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_badVzDiff/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_badVzDiff/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_badVzDiff/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_badVzDiff/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_badVzDiff/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_goodVzDiff/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_goodVzDiff/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_goodVzDiff/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_goodVzDiff/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_goodVzDiff/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_goodVzDiff/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_goodVzDiff/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_goodVzDiff/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_goodVzDiff/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_goodVzDiff/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_goodVzDiff_narrowTimeVeto/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_goodVzDiff_narrowTimeVeto/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_goodVzDiff_strictTimeVeto/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_goodVzDiff_strictTimeVeto/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_goodVzDiff_strictTimeVeto/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_goodVzDiff_strictTimeVeto/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_goodVzDiff_strictTimeVeto/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff_strictTimeVeto/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff_strictTimeVeto/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hVertexChi2VsNcontrib", "", kTH2F, {axisNcontrib, axisVertexChi2perContrib});
    histos.add("noPU_goodVzDiff_strictTimeVeto/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_noPastActivity/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_noPastActivity/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_noPastActivity/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_noPastActivity/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_noPastActivity/hBcFV0", "", kTH1F, {axisBCs});
    // histos.add("noPU_noPastActivity/hBcFDD", "", kTH1F, {axisBCs});
    histos.add("noPU_noPastActivity/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_noPastActivity/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_noPastActivity/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_noPastActivity/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_noPastActivity/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_noPastActivity/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    // histos.add("noPU_noPastActivity/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    // histos.add("noPU_noPastActivity/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    // histos.add("noPU_noPastActivity/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_noPastActivity/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_noPastActivity/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_noPastActivity/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_noPastActivity/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_noFT0activityNearby/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_noFT0activityNearby/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_noFT0activityNearby/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_noFT0activityNearby/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_noFT0activityNearby/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_noFT0activityNearby/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_noFT0activityNearby/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_noFT0activityNearby/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_noFT0activityNearby/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_noFT0activityNearby/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_noFT0activityNearby/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    // histos.add("noPU_noFT0activityNearby/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    // histos.add("noPU_noFT0activityNearby/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    // histos.add("noPU_noFT0activityNearby/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_noFT0activityNearby/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_noFT0activityNearby/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_noFT0activityNearby/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_noFT0activityNearby/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_cutByVzDiff_pvTOF/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_pvTOF/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_pvTOF/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_pvTOF/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_pvTOF/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_pvTOF/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_pvTOF/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_cutByVzDiff_pvTOF/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_cutByVzDiff_pvTOF/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_cutByVzDiff_pvTOF/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_cutByVzDiff_pvTOF/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_cutByVzDiff_pvTOF/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    histos.add("noPU_cutByVzDiff_pvTOF/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    histos.add("noPU_cutByVzDiff_pvTOF/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_cutByVzDiff_pvTOF/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_cutByVzDiff_pvTOF/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_cutByVzDiff_pvTOF/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_cutByVzDiff_pvTOF/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});
    histos.add("noPU_cutByVzDiff_pvTOF/hTimeZN_AC_sum_vs_diff", ";ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDiff, axisTimeSum});
    histos.add("noPU_cutByVzDiff_pvTOF/hTimeZN_A_vs_C", ";ZNA time (ns); ZNC time (ns)", kTH2F, {axisTimeZN, axisTimeZN});
    histos.add("noPU_cutByVzDiff_pvTOF/hTimeZNA", ";ZNA time (ns)", kTH1F, {axisTimeZN});
    histos.add("noPU_cutByVzDiff_pvTOF/hTimeZNC", ";ZNC time (ns)", kTH1F, {axisTimeZN});

    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    // histos.add("noPU_cutByVzDiff_noFT0activityNearby/nTracksGlobal_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracksGlobal});
    // histos.add("noPU_cutByVzDiff_noFT0activityNearby/nTracksGlobal_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracksGlobal});
    // histos.add("noPU_cutByVzDiff_noFT0activityNearby/nTracksGlobal_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracksGlobal});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_cutByVzDiff_noFT0activityNearby/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});

    histos.add("noPU_CutOnZNACtime/hBcColNoSel8", "", kTH1F, {axisBCs});
    histos.add("noPU_CutOnZNACtime/hBcTVX", "", kTH1F, {axisBCs});
    histos.add("noPU_CutOnZNACtime/hBcOrig", "", kTH1F, {axisBCs});
    histos.add("noPU_CutOnZNACtime/hBcFT0", "", kTH1F, {axisBCs});
    histos.add("noPU_CutOnZNACtime/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_CutOnZNACtime/hBcZDC", "", kTH1F, {axisBCs});
    histos.add("noPU_CutOnZNACtime/hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});
    histos.add("noPU_CutOnZNACtime/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_CutOnZNACtime/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_CutOnZNACtime/nTracksPV_vs_T0A", "", kTH2F, {axisMultFT0A, axisNtracks});
    histos.add("noPU_CutOnZNACtime/nTracksPV_vs_T0C", "", kTH2F, {axisMultFT0C, axisNtracks});
    histos.add("noPU_CutOnZNACtime/hTVXvsBcDiffwrtOrigBc", "", kTH1F, {axisBcDiff});
    histos.add("noPU_CutOnZNACtime/hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("noPU_CutOnZNACtime/hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("noPU_CutOnZNACtime/hNPVvsNch", "", kTH2F, {axisNcontrib, axisNcontrib});
    histos.add("noPU_CutOnZNACtime/hTimeZN_AC_sum_vs_diff", ";ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDiff, axisTimeSum});
    histos.add("noPU_CutOnZNACtime/hAmplZNAC", "ZNC vs ZNA", kTH2F, {axisZNampl, axisZNampl});
    histos.add("noPU_CutOnZNACtime/hTimeZN_A_vs_C", ";ZNA time (ns); ZNC time (ns)", kTH2F, {axisTimeZN, axisTimeZN});
    histos.add("noPU_CutOnZNACtime/hTimeZNA", ";ZNA time (ns)", kTH1F, {axisTimeZN});
    histos.add("noPU_CutOnZNACtime/hTimeZNC", ";ZNC time (ns)", kTH1F, {axisTimeZN});

    histos.add("noPU_AntiCutOnZNACtime/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_AntiCutOnZNACtime/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_AntiCutOnZNACtime/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_AntiCutOnZNACtime/hAmplZNAC", "ZNC vs ZNA", kTH2F, {axisZNampl, axisZNampl});

    histos.add("noPU_AntiCutOnZNAampl/hBcFV0", "", kTH1F, {axisBCs});
    histos.add("noPU_AntiCutOnZNAampl/hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M});
    histos.add("noPU_AntiCutOnZNAampl/nTracksPV_vs_V0A", "", kTH2F, {axisMultV0A, axisNtracks});
    histos.add("noPU_AntiCutOnZNAampl/hAmplZNAC", "ZNC vs ZNA", kTH2F, {axisZNampl, axisZNampl});

    //
    histos.add("hNcontribColFromData", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccFromData", "", kTH1F, {axisNcontrib});

    // MC QA
    const AxisSpec axisVtxZdiff{400, -10., 10., ""};
    histos.add("MC/hMCdataVzDiff", "", kTH2F, {axisNcontrib, axisVtxZdiff});
    histos.add("MC/hMCdataBcDiffVsMult", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("MC/hMCdataFoundBcDiffVsMult", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("MCsel8/hMCdataVzDiff", "", kTH2F, {axisNcontrib, axisVtxZdiff});
    histos.add("MCsel8/hMCdataVzDiff_vertTRDmatched", "", kTH2F, {axisNcontrib, axisVtxZdiff});
    histos.add("MCsel8/hMCdataVzDiff_vertTOFmatched", "", kTH2F, {axisNcontrib, axisVtxZdiff});
    histos.add("MCsel8/hMCdataBcDiffVsMult", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("MCsel8/hMCdataFoundBcDiffVsMult", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("MCsel8/hMCdataFoundBcDiffVsMult_vertTRDmatched", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("MCsel8/hMCdataFoundBcDiffVsMult_vertTOFmatched", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("MCnonTVX/hMCdataVzDiff", "", kTH2F, {axisNcontrib, axisVtxZdiff});
    histos.add("MCnonTVX/hMCdataBcDiffVsMult", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("MCnonTVX/hMCdataFoundBcDiffVsMult", "", kTH2F, {axisNcontrib, axisBcDiff});

    //
    histos.add("MC_not_TF_ROF_borders/hNcontribColFromData", "", kTH1F, {axisNcontrib});
    histos.add("MC_not_TF_ROF_borders/hNcontribAccFromData", "", kTH1F, {axisNcontrib});
    histos.add("MC_not_TF_ROF_borders/hNcontribColFromData_foundBcDiff0", "", kTH1F, {axisNcontrib});
    histos.add("MC_not_TF_ROF_borders/hNcontribAccFromData_foundBcDiff0", "", kTH1F, {axisNcontrib});
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
    aod::Zdcs const&,
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
      histos.add("hSecondsCollisions/sel8", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCollisions/noPU", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCollisions/noPU_underLine", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCollisions/noPU_grassOnTheRight", "", kTH1F, {axisSeconds});
      histos.add("hSecondsCollisions/noPU_good", "", kTH1F, {axisSeconds});

      const AxisSpec axisDiffMeanVz{80, -4, 4, ""};
      const AxisSpec axisVzNarrow{40, -10, 10, ""};
      histos.add("hSecondsCollisions/noPU_meanDiffVz", "", kTH2F, {axisSeconds, axisDiffMeanVz});
      histos.add("hSecondsCollisions/noPU_meanDiffVz_lowMult", "", kTH2F, {axisSeconds, axisDiffMeanVz});
      histos.add("hSecondsCollisions/noPU_meanDiffVz_highMult", "", kTH2F, {axisSeconds, axisDiffMeanVz});
      histos.add("hSecondsCollisions/noPU_Vz", "", kTH2F, {axisSeconds, axisVzNarrow});
      histos.add("hSecondsCollisions/noPU_VzByFT0", "", kTH2F, {axisSeconds, axisVzNarrow});

      const AxisSpec axisVz{200, -20, 20, ""};
      histos.add("noSpecSelections/Vz", "", kTH1F, {axisVz});
      histos.add("noPU/Vz", "", kTH1F, {axisVz});
      histos.add("noPU_badVzDiff/Vz", "", kTH1F, {axisVz});
      histos.add("noPU_goodVzDiff/Vz", "", kTH1F, {axisVz});

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
      // bool pastActivityFDD = 0;
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
          // pastActivityFDD |= bcPast.has_fdd();
        }
        if (deltaBC < 2) {
          if (bcPast.has_ft0()) {
            std::bitset<8> triggers = bcPast.ft0().triggerMask();
            nearbyFT0activity |= (triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideA] || triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideC]);
          }
        }
      }
      // bool pastActivity = pastActivityFT0 | pastActivityFV0 | pastActivityFDD;
      bool pastActivity = pastActivityFT0 | pastActivityFV0; // | pastActivityFDD;
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
            nearbyFT0activity |= (triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideA] || triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideC]);
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

        // spec bits:
        std::bitset<8> triggers = bc.ft0().triggerMask();
        bool isTVX = bc.selection_bit(kIsTriggerTVX);
        if (triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideA]) {
          histos.fill(HIST("bcQA/specFT0bits/hBc_kIsActiveSideA"), localBC);
          if (isTVX)
            histos.fill(HIST("bcQA/specFT0bits/hBc_kIsActiveSideA_inTVX"), localBC);
        }
        if (triggers[o2::ft0::RecPoints::ETriggerBits::kIsActiveSideC]) {
          histos.fill(HIST("bcQA/specFT0bits/hBc_kIsActiveSideC"), localBC);
          if (isTVX)
            histos.fill(HIST("bcQA/specFT0bits/hBc_kIsActiveSideC_inTVX"), localBC);
        }
        if (triggers[o2::ft0::RecPoints::ETriggerBits::kIsFlangeEvent]) {
          histos.fill(HIST("bcQA/specFT0bits/hBc_kIsFlangeEvent"), localBC);
          if (isTVX)
            histos.fill(HIST("bcQA/specFT0bits/hBc_kIsFlangeEvent_inTVX"), localBC);
        }
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
      bool pvTRDmatched = col.selection_bit(kIsVertexTRDmatched);
      bool narrowTimeVeto = col.selection_bit(kNoCollInTimeRangeNarrow);
      bool strictTimeVeto = col.selection_bit(kNoCollInTimeRangeStrict);
      bool noCollSameROF = col.selection_bit(kNoCollInRofStrict);
      bool bcDiffWrtClosestTVXCut = (bcToClosestTVXdiff == 0);

      auto bcIndex = foundBC.globalIndex();
      // bool noNearbyActivity = (vPastActivity[bcIndex] == 0 && vFutureActivity[bcIndex] == 0);
      bool noPastActivity = (vPastActivity[bcIndex] == 0);

      // ### count tracks of different types
      int nPVtracks = 0;
      int nGlobalTracksPV = 0;
      int nGlobalTracksAll = 0;
      // int nTOFtracks = 0;
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (const auto& track : tracksGrouped) {
        if (track.itsNCls() < 5)
          continue;
        if (track.pt() < 0.2 || track.pt() > 10)
          continue;
        if (std::abs(track.eta()) > 0.8)
          continue;

        if (track.hasITS() && track.hasTPC() && track.tpcNClsFound() > 50 && track.tpcNClsCrossedRows() > 50 && track.tpcChi2NCl() < 4)
          nGlobalTracksAll++;

        if (!track.isPVContributor()) {
          continue;
        }
        nPVtracks++;
        if (track.hasITS() && track.hasTPC() && track.tpcNClsFound() > 50 && track.tpcNClsCrossedRows() > 50 && track.tpcChi2NCl() < 4)
          nGlobalTracksPV++;
      } // end of track loop

      if (col.selection_bit(kNoTimeFrameBorder) && col.selection_bit(kNoITSROFrameBorder)) {
        histos.fill(HIST("hNcontribColFromData"), nPVtracks);
        if (col.selection_bit(kIsTriggerTVX))
          histos.fill(HIST("hNcontribAccFromData"), nPVtracks);
      }

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

      float timeZNA = foundBC.has_zdc() ? foundBC.zdc().timeZNA() : -999.f;
      float timeZNC = foundBC.has_zdc() ? foundBC.zdc().timeZNC() : -999.f;
      float znSum = timeZNA + timeZNC;
      float znDiff = timeZNA - timeZNC;
      // bool goodZNACtime = fabs(znSum) < 2 && fabs(znDiff) < 2;
      bool goodZNACtime = (timeZNA > -5 && timeZNA < 2) && (timeZNC > -5 && timeZNC < 2);

      float multZNA = foundBC.has_zdc() ? foundBC.zdc().energyCommonZNA() : -999;
      float multZNC = foundBC.has_zdc() ? foundBC.zdc().energyCommonZNC() : -999;
      bool cutZNACampl = multZNA < 400 && multZNC < 400;

      // vZ diff (FT0 vs by tracks)
      bool badVzDiff = 0;
      if (confUseDiffVzCutFromEvSel)
        badVzDiff = !col.selection_bit(kIsGoodZvtxFT0vsPV);
      else {                  // tune by hand
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

        float stdDev = (multT0M > 10) ? 0.144723 + 13.5345 / sqrt(multT0M) : 1.5; // cm
        if (multT0M > 5000)
          stdDev = stdDev > 0.2 ? stdDev : 0.2; // 0.35; // cm
        badVzDiff = diffVz < (meanDiff - stdDev * nSigmaForVzDiff - safetyDiffVzMargin) || diffVz > (meanDiff + stdDev * nSigmaForVzDiff + safetyDiffVzMargin);
      }

      bool underLine = false;
      if (hasFT0 && nPVtracks < 45. / 40000 * multV0A - 4 && nPVtracks < 20) {
        underLine = true;
      }
      bool grassOnTheRight = false;
      if (hasFT0 && nPVtracks < 220. / 40000 * multV0A - 100 && nPVtracks >= 25) {
        grassOnTheRight = true;
      }

      // study bc diff wrt original bc:
      // auto bc = col.bc_as<BCsRun3>();
      auto bcOriginal = globalOrigBC % 3564;
      float bcDiffWrtOriginal = bcOriginal - localBC;
      bool bcDiffWrtOriginalBcCut = (bcDiffWrtOriginal == 0);

      // SPEC REMOVAL OF BC RANGES
      // if (bcOriginal > 1200 && bcOriginal < 1300)
      //   continue;
      // if (bcOriginal > 1700 && bcOriginal < 1900)
      //   continue;

      int nContributors = col.numContrib();
      float timeRes = col.collisionTimeRes();
      // int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      // if (col.selection_bit(kIsVertexTOFmatched)) {
      //   histos.fill(HIST("hColBcDiffVsNcontribWithTOF"), nContributors, bcToClosestTVXdiff);
      //   histos.fill(HIST("hColTimeResVsNcontribWithTOF"), nContributors, timeRes);
      //   histos.fill(HIST("hBcColTOF"), localBC);
      // }

      float vChi2 = col.chi2();
      float vChi2perContrib = nContributors > 0 ? vChi2 / nContributors : 0;
      bool goodVertexChi2 = (vChi2perContrib < 3.5);

      histos.fill(HIST("noSpecSelections/hBcColNoSel8"), localBC);
      histos.fill(HIST("noSpecSelections/hBcOrigNoSel8"), bcOriginal);
      histos.fill(HIST("noSpecSelections/Vz"), vZ);
      histos.fill(HIST("noSpecSelections/hTimeZN_AC_sum_vs_diff"), znDiff, znSum);
      histos.fill(HIST("noSpecSelections/hTimeZN_A_vs_C"), timeZNA, timeZNC);
      histos.fill(HIST("noSpecSelections/hTimeZNA"), timeZNA);
      histos.fill(HIST("noSpecSelections/hTimeZNC"), timeZNC);

      if (noPU) {
        histos.fill(HIST("noPU/hBcColNoSel8"), localBC);
        histos.fill(HIST("noPU/hBcOrigNoSel8"), bcOriginal);
        histos.fill(HIST("noPU/Vz"), vZ);
        histos.fill(HIST("noPU/hTimeZN_AC_sum_vs_diff"), znDiff, znSum);
        histos.fill(HIST("noPU/hTimeZN_A_vs_C"), timeZNA, timeZNC);
        histos.fill(HIST("noPU/hAmplZNAC"), multZNA, multZNC);
        histos.fill(HIST("noPU/hTimeZNA"), timeZNA);
        histos.fill(HIST("noPU/hTimeZNC"), timeZNC);
      }
      if (noPU && pvTOFmatched) {
        histos.fill(HIST("noPU_pvTOFmatched/hBcColNoSel8"), localBC);
      }
      if (noPU && pvTRDmatched) {
        histos.fill(HIST("noPU_pvTRDmatched/hBcColNoSel8"), localBC);
      }
      if (noPU && !pvTRDmatched) {
        histos.fill(HIST("noPU_notTRDmatched/hBcColNoSel8"), localBC);
      }
      if (bcDiffWrtClosestTVXCut) {
        histos.fill(HIST("bcDiffWrtClosestTVXCut/hBcColNoSel8"), localBC);
      }
      if (noPU && bcDiffWrtOriginalBcCut) {
        histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hBcColNoSel8"), localBC);
      }
      if (noPU && noPastActivity) {
        histos.fill(HIST("noPU_noPastActivity/hBcColNoSel8"), localBC);
      }
      if (noPU && noFT0activityNearby) {
        histos.fill(HIST("noPU_noFT0activityNearby/hBcColNoSel8"), localBC);
      }
      if (noPU && badVzDiff) {
        histos.fill(HIST("noPU_badVzDiff/hBcColNoSel8"), localBC);
        histos.fill(HIST("noPU_badVzDiff/Vz"), vZ);
      }
      if (noPU && !badVzDiff) {
        histos.fill(HIST("noPU_goodVzDiff/hBcColNoSel8"), localBC);
        histos.fill(HIST("noPU_goodVzDiff/Vz"), vZ);
      }
      if (noPU && !badVzDiff && narrowTimeVeto) {
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hBcColNoSel8"), localBC);
      }
      if (noPU && !badVzDiff && strictTimeVeto) {
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hBcColNoSel8"), localBC);
      }
      if (noPU && goodVertexChi2) {
        histos.fill(HIST("noPU_goodVertexChi2/hBcColNoSel8"), localBC);
      }
      if (narrowTimeVeto) {
        histos.fill(HIST("narrowTimeVeto/hBcColNoSel8"), localBC);
      }
      if (strictTimeVeto) {
        histos.fill(HIST("strictTimeVeto/hBcColNoSel8"), localBC);
      }
      if (noCollSameROF) {
        histos.fill(HIST("noCollSameROF/hBcColNoSel8"), localBC);
      }
      if (underLine) {
        histos.fill(HIST("lowMultCut/hBcColNoSel8"), localBC);
        histos.fill(HIST("lowMultCut/hBcOrigNoSel8"), bcOriginal);
      }
      if (noPU && underLine) {
        histos.fill(HIST("noPU_lowMultCut/hBcColNoSel8"), localBC);
        histos.fill(HIST("noPU_lowMultCut/hBcOrigNoSel8"), bcOriginal);
        histos.fill(HIST("noPU_lowMultCut/hTimeZN_AC_sum_vs_diff"), znDiff, znSum);
        histos.fill(HIST("noPU_lowMultCut/hTimeZN_A_vs_C"), timeZNA, timeZNC);
        histos.fill(HIST("noPU_lowMultCut/hTimeZNA"), timeZNA);
        histos.fill(HIST("noPU_lowMultCut/hTimeZNC"), timeZNC);
      }
      if (grassOnTheRight) {
        histos.fill(HIST("highMultCloudCut/hBcColNoSel8"), localBC);
        histos.fill(HIST("highMultCloudCut/hBcOrigNoSel8"), bcOriginal);
      }
      if (noPU && grassOnTheRight) {
        histos.fill(HIST("noPU_highMultCloudCut/hBcColNoSel8"), localBC);
        histos.fill(HIST("noPU_highMultCloudCut/hBcOrigNoSel8"), bcOriginal);
        histos.fill(HIST("noPU_highMultCloudCut/hTimeZN_AC_sum_vs_diff"), znDiff, znSum);
        histos.fill(HIST("noPU_highMultCloudCut/hTimeZN_A_vs_C"), timeZNA, timeZNC);
        histos.fill(HIST("noPU_highMultCloudCut/hTimeZNA"), timeZNA);
        histos.fill(HIST("noPU_highMultCloudCut/hTimeZNC"), timeZNC);
      }
      if (noPU && !badVzDiff && pvTOFmatched) { // noPileup_cutByVzDiff_pvTOF_noFT0act
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hBcColNoSel8"), localBC);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hTimeZN_AC_sum_vs_diff"), znDiff, znSum);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hTimeZN_A_vs_C"), timeZNA, timeZNC);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hTimeZNA"), timeZNA);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hTimeZNC"), timeZNC);
      }

      // only here cut on sel8:
      if (!col.sel8())
        continue;

      // vs time
      histos.fill(HIST("hSecondsCollisions/sel8"), secFromSOR);
      if (noPU) {
        histos.fill(HIST("hSecondsCollisions/noPU"), secFromSOR);
        if (underLine)
          histos.fill(HIST("hSecondsCollisions/noPU_underLine"), secFromSOR);
        if (grassOnTheRight)
          histos.fill(HIST("hSecondsCollisions/noPU_grassOnTheRight"), secFromSOR);
        if (!underLine && !grassOnTheRight)
          histos.fill(HIST("hSecondsCollisions/noPU_good"), secFromSOR);

        histos.fill(HIST("hSecondsCollisions/noPU_Vz"), secFromSOR, vZ);
        histos.fill(HIST("hSecondsCollisions/noPU_VzByFT0"), secFromSOR, vZft0);

        if (std::abs(diffVz) < 4) {
          histos.fill(HIST("hSecondsCollisions/noPU_meanDiffVz"), secFromSOR, diffVz);
          if (multT0M < 1000)
            histos.fill(HIST("hSecondsCollisions/noPU_meanDiffVz_lowMult"), secFromSOR, diffVz);
          if (multT0M > 10000)
            histos.fill(HIST("hSecondsCollisions/noPU_meanDiffVz_highMult"), secFromSOR, diffVz);
        }
      }

      histos.fill(HIST("noSpecSelections/hBcTVX"), localBC);
      histos.fill(HIST("noSpecSelections/hBcOrig"), bcOriginal);
      histos.fill(HIST("noSpecSelections/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
      histos.fill(HIST("noSpecSelections/hColTimeResVsNcontrib"), nContributors, timeRes);
      histos.fill(HIST("noSpecSelections/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
      histos.fill(HIST("noSpecSelections/hNPVvsNch"), nPVtracks, nGlobalTracksAll);

      if (noPU) {
        histos.fill(HIST("noPU/hBcTVX"), localBC);
        histos.fill(HIST("noPU/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && pvTOFmatched) {
        histos.fill(HIST("noPU_pvTOFmatched/hBcTVX"), localBC);
        histos.fill(HIST("noPU_pvTOFmatched/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_pvTOFmatched/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_pvTOFmatched/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_pvTOFmatched/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_pvTOFmatched/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && pvTRDmatched) {
        histos.fill(HIST("noPU_pvTRDmatched/hBcTVX"), localBC);
        histos.fill(HIST("noPU_pvTRDmatched/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_pvTRDmatched/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_pvTRDmatched/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_pvTRDmatched/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_pvTRDmatched/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && !pvTRDmatched) {
        histos.fill(HIST("noPU_notTRDmatched/hBcTVX"), localBC);
        histos.fill(HIST("noPU_notTRDmatched/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_notTRDmatched/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_notTRDmatched/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_notTRDmatched/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (bcDiffWrtClosestTVXCut) {
        histos.fill(HIST("bcDiffWrtClosestTVXCut/hBcTVX"), localBC);
        histos.fill(HIST("bcDiffWrtClosestTVXCut/hBcOrig"), bcOriginal);
        histos.fill(HIST("bcDiffWrtClosestTVXCut/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("bcDiffWrtClosestTVXCut/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("bcDiffWrtClosestTVXCut/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && bcDiffWrtOriginalBcCut) {
        histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hBcTVX"), localBC);
        histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && noPastActivity) {
        histos.fill(HIST("noPU_noPastActivity/hBcTVX"), localBC);
        histos.fill(HIST("noPU_noPastActivity/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_noPastActivity/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_noPastActivity/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_noPastActivity/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && noFT0activityNearby) {
        histos.fill(HIST("noPU_noFT0activityNearby/hBcTVX"), localBC);
        histos.fill(HIST("noPU_noFT0activityNearby/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_noFT0activityNearby/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_noFT0activityNearby/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_noFT0activityNearby/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && badVzDiff) {
        histos.fill(HIST("noPU_badVzDiff/hBcTVX"), localBC);
        histos.fill(HIST("noPU_badVzDiff/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_badVzDiff/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_badVzDiff/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_badVzDiff/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_badVzDiff/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && !badVzDiff) {
        histos.fill(HIST("noPU_goodVzDiff/hBcTVX"), localBC);
        histos.fill(HIST("noPU_goodVzDiff/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_goodVzDiff/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_goodVzDiff/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_goodVzDiff/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_goodVzDiff/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && !badVzDiff && narrowTimeVeto) {
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hBcTVX"), localBC);
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && !badVzDiff && strictTimeVeto) {
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hBcTVX"), localBC);
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && goodVertexChi2) {
        histos.fill(HIST("noPU_goodVertexChi2/hBcTVX"), localBC);
        histos.fill(HIST("noPU_goodVertexChi2/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_goodVertexChi2/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_goodVertexChi2/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_goodVertexChi2/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (narrowTimeVeto) {
        histos.fill(HIST("narrowTimeVeto/hBcTVX"), localBC);
        histos.fill(HIST("narrowTimeVeto/hBcOrig"), bcOriginal);
        histos.fill(HIST("narrowTimeVeto/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("narrowTimeVeto/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("narrowTimeVeto/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (strictTimeVeto) {
        histos.fill(HIST("strictTimeVeto/hBcTVX"), localBC);
        histos.fill(HIST("strictTimeVeto/hBcOrig"), bcOriginal);
        histos.fill(HIST("strictTimeVeto/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("strictTimeVeto/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("strictTimeVeto/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noCollSameROF) {
        histos.fill(HIST("noCollSameROF/hBcTVX"), localBC);
        histos.fill(HIST("noCollSameROF/hBcOrig"), bcOriginal);
        histos.fill(HIST("noCollSameROF/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noCollSameROF/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noCollSameROF/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && underLine) {
        histos.fill(HIST("noPU_lowMultCut/hBcTVX"), localBC);
        histos.fill(HIST("noPU_lowMultCut/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_lowMultCut/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_lowMultCut/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_lowMultCut/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_lowMultCut/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && grassOnTheRight) {
        histos.fill(HIST("noPU_highMultCloudCut/hBcTVX"), localBC);
        histos.fill(HIST("noPU_highMultCloudCut/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_highMultCloudCut/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_highMultCloudCut/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_highMultCloudCut/hVertexChi2VsNcontrib"), nContributors, vChi2perContrib);
        histos.fill(HIST("noPU_highMultCloudCut/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && !badVzDiff && pvTOFmatched) { // noPileup_cutByVzDiff_pvTOF_noFT0act
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hBcTVX"), localBC);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && !badVzDiff && noFT0activityNearby) {
        histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hBcTVX"), localBC);
        histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
      }
      if (noPU && goodZNACtime) {
        histos.fill(HIST("noPU_CutOnZNACtime/hBcTVX"), localBC);
        histos.fill(HIST("noPU_CutOnZNACtime/hBcOrig"), bcOriginal);
        histos.fill(HIST("noPU_CutOnZNACtime/hColBcDiffVsNcontrib"), nContributors, bcToClosestTVXdiff);
        histos.fill(HIST("noPU_CutOnZNACtime/hColTimeResVsNcontrib"), nContributors, timeRes);
        histos.fill(HIST("noPU_CutOnZNACtime/hNPVvsNch"), nPVtracks, nGlobalTracksAll);
        histos.fill(HIST("noPU_CutOnZNACtime/hTimeZN_AC_sum_vs_diff"), znDiff, znSum);
        histos.fill(HIST("noPU_CutOnZNACtime/hTimeZN_A_vs_C"), timeZNA, timeZNC);
        histos.fill(HIST("noPU_CutOnZNACtime/hTimeZNA"), timeZNA);
        histos.fill(HIST("noPU_CutOnZNACtime/hTimeZNC"), timeZNC);
      }

      if (foundBC.has_ft0()) {
        // float multT0A = foundBC.has_ft0() ? foundBC.ft0().sumAmpA() : -999.f;
        // float multT0C = foundBC.has_ft0() ? foundBC.ft0().sumAmpC() : -999.f;

        histos.fill(HIST("noSpecSelections/hBcFT0"), localBC);
        histos.fill(HIST("noSpecSelections/hVtxFT0VsVtxCol"), vZft0, vZ);
        histos.fill(HIST("noSpecSelections/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
        histos.fill(HIST("noSpecSelections/nTracksPV_vs_T0A"), multT0A, nPVtracks);
        histos.fill(HIST("noSpecSelections/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
        histos.fill(HIST("noSpecSelections/nTracksPV_vs_T0C"), multT0C, nPVtracks);
        histos.fill(HIST("noSpecSelections/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);

        if (noPU) {
          histos.fill(HIST("noPU/hBcFT0"), localBC);
          histos.fill(HIST("noPU/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && pvTOFmatched) {
          histos.fill(HIST("noPU_pvTOFmatched/hBcFT0"), localBC);
          histos.fill(HIST("noPU_pvTOFmatched/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_pvTOFmatched/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_pvTOFmatched/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_pvTOFmatched/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_pvTOFmatched/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_pvTOFmatched/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && pvTRDmatched) {
          histos.fill(HIST("noPU_pvTRDmatched/hBcFT0"), localBC);
          histos.fill(HIST("noPU_pvTRDmatched/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_pvTRDmatched/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_pvTRDmatched/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_pvTRDmatched/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_pvTRDmatched/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_pvTRDmatched/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && !pvTRDmatched) {
          histos.fill(HIST("noPU_notTRDmatched/hBcFT0"), localBC);
          histos.fill(HIST("noPU_notTRDmatched/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_notTRDmatched/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_notTRDmatched/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_notTRDmatched/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_notTRDmatched/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_notTRDmatched/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (bcDiffWrtClosestTVXCut) {
          histos.fill(HIST("bcDiffWrtClosestTVXCut/hBcFT0"), localBC);
          histos.fill(HIST("bcDiffWrtClosestTVXCut/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("bcDiffWrtClosestTVXCut/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("bcDiffWrtClosestTVXCut/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          // histos.fill(HIST("bcDiffWrtClosestTVXCut/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("bcDiffWrtClosestTVXCut/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          // histos.fill(HIST("bcDiffWrtClosestTVXCut/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && bcDiffWrtOriginalBcCut) {
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hBcFT0"), localBC);
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          // histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          // histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && badVzDiff) {
          histos.fill(HIST("noPU_badVzDiff/hBcFT0"), localBC);
          histos.fill(HIST("noPU_badVzDiff/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_badVzDiff/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_badVzDiff/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_badVzDiff/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_badVzDiff/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_badVzDiff/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff) {
          histos.fill(HIST("noPU_goodVzDiff/hBcFT0"), localBC);
          histos.fill(HIST("noPU_goodVzDiff/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_goodVzDiff/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_goodVzDiff/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_goodVzDiff/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && narrowTimeVeto) {
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hBcFT0"), localBC);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && strictTimeVeto) {
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hBcFT0"), localBC);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && noPastActivity) {
          histos.fill(HIST("noPU_noPastActivity/hBcFT0"), localBC);
          histos.fill(HIST("noPU_noPastActivity/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_noPastActivity/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_noPastActivity/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          // histos.fill(HIST("noPU_noPastActivity/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_noPastActivity/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          // histos.fill(HIST("noPU_noPastActivity/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && noFT0activityNearby) {
          histos.fill(HIST("noPU_noFT0activityNearby/hBcFT0"), localBC);
          histos.fill(HIST("noPU_noFT0activityNearby/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_noFT0activityNearby/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_noFT0activityNearby/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          // histos.fill(HIST("noPU_noFT0activityNearby/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_noFT0activityNearby/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          // histos.fill(HIST("noPU_noFT0activityNearby/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && goodVertexChi2) {
          histos.fill(HIST("noPU_goodVertexChi2/hBcFT0"), localBC);
          histos.fill(HIST("noPU_goodVertexChi2/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_goodVertexChi2/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_goodVertexChi2/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          // histos.fill(HIST("noPU_goodVertexChi2/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_goodVertexChi2/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          // histos.fill(HIST("noPU_goodVertexChi2/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (narrowTimeVeto) {
          histos.fill(HIST("narrowTimeVeto/hBcFT0"), localBC);
          histos.fill(HIST("narrowTimeVeto/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("narrowTimeVeto/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("narrowTimeVeto/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("narrowTimeVeto/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("narrowTimeVeto/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("narrowTimeVeto/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (strictTimeVeto) {
          histos.fill(HIST("strictTimeVeto/hBcFT0"), localBC);
          histos.fill(HIST("strictTimeVeto/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("strictTimeVeto/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("strictTimeVeto/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("strictTimeVeto/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("strictTimeVeto/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("strictTimeVeto/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noCollSameROF) {
          histos.fill(HIST("noCollSameROF/hBcFT0"), localBC);
          histos.fill(HIST("noCollSameROF/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noCollSameROF/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noCollSameROF/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noCollSameROF/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noCollSameROF/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noCollSameROF/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && underLine) {
          histos.fill(HIST("noPU_lowMultCut/hBcFT0"), localBC);
          histos.fill(HIST("noPU_lowMultCut/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_lowMultCut/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_lowMultCut/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_lowMultCut/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_lowMultCut/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_lowMultCut/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && grassOnTheRight) {
          histos.fill(HIST("noPU_highMultCloudCut/hBcFT0"), localBC);
          histos.fill(HIST("noPU_highMultCloudCut/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_highMultCloudCut/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_highMultCloudCut/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_highMultCloudCut/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_highMultCloudCut/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_highMultCloudCut/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && pvTOFmatched) { // noPileup_cutByVzDiff_pvTOF_noFT0act
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hBcFT0"), localBC);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && noFT0activityNearby) {
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hBcFT0"), localBC);
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          // histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/nTracksGlobal_vs_T0A"), multT0A, nGlobalTracksPV);
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          // histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/nTracksGlobal_vs_T0C"), multT0C, nGlobalTracksPV);
        }
        if (noPU && goodZNACtime) {
          histos.fill(HIST("noPU_CutOnZNACtime/hBcFT0"), localBC);
          histos.fill(HIST("noPU_CutOnZNACtime/hVtxFT0VsVtxCol"), vZft0, vZ);
          histos.fill(HIST("noPU_CutOnZNACtime/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_CutOnZNACtime/nTracksPV_vs_T0A"), multT0A, nPVtracks);
          histos.fill(HIST("noPU_CutOnZNACtime/nTracksPV_vs_T0C"), multT0C, nPVtracks);
          histos.fill(HIST("noPU_CutOnZNACtime/hAmplZNAC"), multZNA, multZNC);
        }
      }

      if (foundBC.has_fv0a()) {
        histos.fill(HIST("noSpecSelections/hBcFV0"), localBC);
        histos.fill(HIST("noSpecSelections/nTracksPV_vs_V0A"), multV0A, nPVtracks);
        histos.fill(HIST("noSpecSelections/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        if (noPU) {
          histos.fill(HIST("noPU/hBcFV0"), localBC);
          histos.fill(HIST("noPU/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && pvTOFmatched) {
          histos.fill(HIST("noPU_pvTOFmatched/hBcFV0"), localBC);
          histos.fill(HIST("noPU_pvTOFmatched/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_pvTOFmatched/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && pvTRDmatched) {
          histos.fill(HIST("noPU_pvTRDmatched/hBcFV0"), localBC);
          histos.fill(HIST("noPU_pvTRDmatched/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_pvTRDmatched/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && !pvTRDmatched) {
          histos.fill(HIST("noPU_notTRDmatched/hBcFV0"), localBC);
          histos.fill(HIST("noPU_notTRDmatched/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_notTRDmatched/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && pvTOFmatched && !pvTRDmatched) { // SPEC CHECK!
          histos.fill(HIST("noPU_pvTOFmatched_notTRDmatched/hBcFV0"), localBC);
          histos.fill(HIST("noPU_pvTOFmatched_notTRDmatched/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_pvTOFmatched_notTRDmatched/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (bcDiffWrtClosestTVXCut) {
          histos.fill(HIST("bcDiffWrtClosestTVXCut/hBcFV0"), localBC);
          histos.fill(HIST("bcDiffWrtClosestTVXCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          // histos.fill(HIST("bcDiffWrtClosestTVXCut/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && bcDiffWrtOriginalBcCut) {
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hBcFV0"), localBC);
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          // histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && badVzDiff) {
          histos.fill(HIST("noPU_badVzDiff/hBcFV0"), localBC);
          histos.fill(HIST("noPU_badVzDiff/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_badVzDiff/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff) {
          histos.fill(HIST("noPU_goodVzDiff/hBcFV0"), localBC);
          histos.fill(HIST("noPU_goodVzDiff/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && narrowTimeVeto) {
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hBcFV0"), localBC);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && strictTimeVeto) {
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hBcFV0"), localBC);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && noPastActivity) {
          histos.fill(HIST("noPU_noPastActivity/hBcFV0"), localBC);
          histos.fill(HIST("noPU_noPastActivity/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          // histos.fill(HIST("noPU_noPastActivity/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && noFT0activityNearby) {
          histos.fill(HIST("noPU_noFT0activityNearby/hBcFV0"), localBC);
          histos.fill(HIST("noPU_noFT0activityNearby/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          // histos.fill(HIST("noPU_noFT0activityNearby/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && goodVertexChi2) {
          histos.fill(HIST("noPU_goodVertexChi2/hBcFV0"), localBC);
          histos.fill(HIST("noPU_goodVertexChi2/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          // histos.fill(HIST("noPU_goodVertexChi2/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (narrowTimeVeto) {
          histos.fill(HIST("narrowTimeVeto/hBcFV0"), localBC);
          histos.fill(HIST("narrowTimeVeto/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("narrowTimeVeto/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (strictTimeVeto) {
          histos.fill(HIST("strictTimeVeto/hBcFV0"), localBC);
          histos.fill(HIST("strictTimeVeto/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("strictTimeVeto/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noCollSameROF) {
          histos.fill(HIST("noCollSameROF/hBcFV0"), localBC);
          histos.fill(HIST("noCollSameROF/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noCollSameROF/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (underLine) {
          histos.fill(HIST("lowMultCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
        }
        if (noPU && underLine) {
          histos.fill(HIST("noPU_lowMultCut/hBcFV0"), localBC);
          histos.fill(HIST("noPU_lowMultCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_lowMultCut/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (grassOnTheRight) {
          histos.fill(HIST("highMultCloudCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
        }
        if (noPU && grassOnTheRight) {
          histos.fill(HIST("noPU_highMultCloudCut/hBcFV0"), localBC);
          histos.fill(HIST("noPU_highMultCloudCut/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_highMultCloudCut/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && pvTOFmatched) { // noPileup_cutByVzDiff_pvTOF_noFT0act
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hBcFV0"), localBC);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }
        if (noPU && !badVzDiff && noFT0activityNearby) {
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hBcFV0"), localBC);
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          // histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/nTracksGlobal_vs_V0A"), multV0A, nGlobalTracksPV);
        }

        if (noPU && goodZNACtime) {
          histos.fill(HIST("noPU_CutOnZNACtime/hBcFV0"), localBC);
          histos.fill(HIST("noPU_CutOnZNACtime/nTracksPV_vs_V0A"), multV0A, nPVtracks);
        }
        if (noPU && !goodZNACtime) {
          histos.fill(HIST("noPU_AntiCutOnZNACtime/hBcFV0"), localBC);
          histos.fill(HIST("noPU_AntiCutOnZNACtime/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_AntiCutOnZNACtime/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_AntiCutOnZNACtime/hAmplZNAC"), multZNA, multZNC);
        }
        if (noPU && !cutZNACampl) {
          histos.fill(HIST("noPU_AntiCutOnZNAampl/hBcFV0"), localBC);
          histos.fill(HIST("noPU_AntiCutOnZNAampl/nTracksPV_vs_V0A"), multV0A, nPVtracks);
          histos.fill(HIST("noPU_AntiCutOnZNAampl/hVtxFT0MinusVtxColVsMultT0M"), diffVz, multT0A + multT0C);
          histos.fill(HIST("noPU_AntiCutOnZNAampl/hAmplZNAC"), multZNA, multZNC);
        }
      }
      if (foundBC.has_zdc()) {
        histos.fill(HIST("noSpecSelections/hBcZDC"), localBC);
        if (noPU) {
          histos.fill(HIST("noPU/hBcZDC"), localBC);
        }
        if (noPU && pvTOFmatched) {
          histos.fill(HIST("noPU_pvTOFmatched/hBcZDC"), localBC);
        }
        if (noPU && pvTRDmatched) {
          histos.fill(HIST("noPU_pvTRDmatched/hBcZDC"), localBC);
        }
        if (noPU && !pvTRDmatched) {
          histos.fill(HIST("noPU_notTRDmatched/hBcZDC"), localBC);
        }
        if (bcDiffWrtClosestTVXCut) {
          histos.fill(HIST("bcDiffWrtClosestTVXCut/hBcZDC"), localBC);
        }
        if (noPU && bcDiffWrtOriginalBcCut) {
          histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hBcZDC"), localBC);
        }
        if (noPU && badVzDiff) {
          histos.fill(HIST("noPU_badVzDiff/hBcZDC"), localBC);
        }
        if (noPU && !badVzDiff) {
          histos.fill(HIST("noPU_goodVzDiff/hBcZDC"), localBC);
        }
        if (noPU && !badVzDiff && narrowTimeVeto) {
          histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hBcZDC"), localBC);
        }
        if (noPU && !badVzDiff && strictTimeVeto) {
          histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hBcZDC"), localBC);
        }
        if (noPU && noPastActivity) {
          histos.fill(HIST("noPU_noPastActivity/hBcZDC"), localBC);
        }
        if (noPU && noFT0activityNearby) {
          histos.fill(HIST("noPU_noFT0activityNearby/hBcZDC"), localBC);
        }
        if (noPU && goodVertexChi2) {
          histos.fill(HIST("noPU_goodVertexChi2/hBcZDC"), localBC);
        }
        if (narrowTimeVeto) {
          histos.fill(HIST("narrowTimeVeto/hBcZDC"), localBC);
        }
        if (strictTimeVeto) {
          histos.fill(HIST("strictTimeVeto/hBcZDC"), localBC);
        }
        if (noCollSameROF) {
          histos.fill(HIST("noCollSameROF/hBcZDC"), localBC);
        }
        if (noPU && underLine) {
          histos.fill(HIST("noPU_lowMultCut/hBcZDC"), localBC);
        }
        if (noPU && grassOnTheRight) {
          histos.fill(HIST("noPU_highMultCloudCut/hBcZDC"), localBC);
        }
        if (noPU && !badVzDiff && pvTOFmatched) { // noPileup_cutByVzDiff_pvTOF_noFT0act
          histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hBcZDC"), localBC);
        }
        if (noPU && !badVzDiff && noFT0activityNearby) {
          histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hBcZDC"), localBC);
        }
      }

      // bc diff wrt original bc
      histos.fill(HIST("noSpecSelections/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      if (noPU) {
        histos.fill(HIST("noPU/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && pvTOFmatched) {
        histos.fill(HIST("noPU_pvTOFmatched/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && pvTRDmatched) {
        histos.fill(HIST("noPU_pvTRDmatched/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && !pvTRDmatched) {
        histos.fill(HIST("noPU_notTRDmatched/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (bcDiffWrtClosestTVXCut) {
        histos.fill(HIST("bcDiffWrtClosestTVXCut/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && bcDiffWrtOriginalBcCut) {
        histos.fill(HIST("noPU_bcDiffWrtOriginalBcCut/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && badVzDiff) {
        histos.fill(HIST("noPU_badVzDiff/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && !badVzDiff) {
        histos.fill(HIST("noPU_goodVzDiff/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && !badVzDiff && narrowTimeVeto) {
        histos.fill(HIST("noPU_goodVzDiff_narrowTimeVeto/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && !badVzDiff && strictTimeVeto) {
        histos.fill(HIST("noPU_goodVzDiff_strictTimeVeto/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && noPastActivity) {
        histos.fill(HIST("noPU_noPastActivity/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && noFT0activityNearby) {
        histos.fill(HIST("noPU_noFT0activityNearby/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && goodVertexChi2) {
        histos.fill(HIST("noPU_goodVertexChi2/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (narrowTimeVeto) {
        histos.fill(HIST("narrowTimeVeto/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (strictTimeVeto) {
        histos.fill(HIST("strictTimeVeto/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noCollSameROF) {
        histos.fill(HIST("noCollSameROF/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && underLine) {
        histos.fill(HIST("noPU_lowMultCut/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && grassOnTheRight) {
        histos.fill(HIST("noPU_highMultCloudCut/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && !badVzDiff && pvTOFmatched) { // noPileup_cutByVzDiff_pvTOF_noFT0act
        histos.fill(HIST("noPU_cutByVzDiff_pvTOF/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && !badVzDiff && noFT0activityNearby) {
        histos.fill(HIST("noPU_cutByVzDiff_noFT0activityNearby/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
      if (noPU && goodZNACtime) {
        histos.fill(HIST("noPU_CutOnZNACtime/hTVXvsBcDiffwrtOrigBc"), bcDiffWrtOriginal);
      }
    } // end of collisions loop
  }
  PROCESS_SWITCH(LightIonsEvSelQa, processRun3, "Process Run3 tracking vs detector occupancy QA", true);

  // ### MC QA
  using ColEvSelsWithMCLabels = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>; //, aod::CentFT0Cs>;
  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  void processMC(ColEvSelsWithMCLabels const& collisions,
                 BCsInfo const&,
                 aod::McCollisions const&)
  {
    for (const auto& col : collisions) {
      if (fabs(col.posZ()) > 10)
        continue;
      bool isSel8 = col.sel8();
      if (col.has_mcCollision()) {
        const auto mcCollision = col.mcCollision();
        LOGP(debug, "col.posZ() = {}, mcCollision.posZ() = {}", col.posZ(), mcCollision.posZ());

        float diffVz = col.posZ() - mcCollision.posZ();
        histos.fill(HIST("MC/hMCdataVzDiff"), col.numContrib(), diffVz);

        auto bc = col.bc_as<BCsInfo>();
        auto mcBc = mcCollision.bc_as<BCsInfo>();
        auto foundBC = col.foundBC_as<BCsInfo>();

        uint64_t globalBC = bc.globalBC();
        uint64_t foundGlobalBC = foundBC.globalBC();
        uint64_t mcGlobalBC = mcBc.globalBC();

        int bcDiff = static_cast<int>(globalBC - mcGlobalBC);
        int foundBcDiff = static_cast<int>(foundGlobalBC - mcGlobalBC);

        // restrict bc diff range - to see values in the diff histograms
        if (bcDiff > 300)
          bcDiff = 300;
        if (bcDiff < -300)
          bcDiff = -300;
        if (foundBcDiff > 300)
          foundBcDiff = 300;
        if (foundBcDiff < -300)
          foundBcDiff = -300;

        LOGP(debug, "globalBC = {}, mcGlobalBC = {}", globalBC, mcGlobalBC);
        histos.fill(HIST("MC/hMCdataBcDiffVsMult"), col.numContrib(), bcDiff);
        histos.fill(HIST("MC/hMCdataFoundBcDiffVsMult"), col.numContrib(), foundBcDiff);

        if (isSel8) {
          histos.fill(HIST("MCsel8/hMCdataVzDiff"), col.numContrib(), diffVz);
          histos.fill(HIST("MCsel8/hMCdataBcDiffVsMult"), col.numContrib(), bcDiff);
          histos.fill(HIST("MCsel8/hMCdataFoundBcDiffVsMult"), col.numContrib(), foundBcDiff);

          if (col.selection_bit(kIsVertexTRDmatched)) {
            histos.fill(HIST("MCsel8/hMCdataVzDiff_vertTRDmatched"), col.numContrib(), diffVz);
            histos.fill(HIST("MCsel8/hMCdataFoundBcDiffVsMult_vertTRDmatched"), col.numContrib(), foundBcDiff);
          }
          if (col.selection_bit(kIsVertexTOFmatched)) {
            histos.fill(HIST("MCsel8/hMCdataVzDiff_vertTOFmatched"), col.numContrib(), diffVz);
            histos.fill(HIST("MCsel8/hMCdataFoundBcDiffVsMult_vertTOFmatched"), col.numContrib(), foundBcDiff);
          }
        }
        if (col.selection_bit(kNoTimeFrameBorder) && col.selection_bit(kNoITSROFrameBorder) && !col.selection_bit(kIsTriggerTVX)) {
          histos.fill(HIST("MCnonTVX/hMCdataVzDiff"), col.numContrib(), diffVz);
          histos.fill(HIST("MCnonTVX/hMCdataBcDiffVsMult"), col.numContrib(), bcDiff);
          histos.fill(HIST("MCnonTVX/hMCdataFoundBcDiffVsMult"), col.numContrib(), foundBcDiff);
        }

        if (col.selection_bit(kNoTimeFrameBorder) && col.selection_bit(kNoITSROFrameBorder)) {
          histos.fill(HIST("MC_not_TF_ROF_borders/hNcontribColFromData"), col.numContrib());
          if (col.selection_bit(kIsTriggerTVX))
            histos.fill(HIST("MC_not_TF_ROF_borders/hNcontribAccFromData"), col.numContrib());

          if (foundBcDiff == 0) {
            histos.fill(HIST("MC_not_TF_ROF_borders/hNcontribColFromData_foundBcDiff0"), col.numContrib());
            if (col.selection_bit(kIsTriggerTVX))
              histos.fill(HIST("MC_not_TF_ROF_borders/hNcontribAccFromData_foundBcDiff0"), col.numContrib());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(LightIonsEvSelQa, processMC, "Process MC", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LightIonsEvSelQa>(cfgc)};
}
