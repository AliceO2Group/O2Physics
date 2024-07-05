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

#include "map"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "CommonDataFormat/BunchFilling.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3.h"

using namespace o2::framework;
using namespace o2;
using namespace o2::aod::evsel;

using BCsRun2 = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::BcSels, aod::Run2MatchedToBCSparse>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection>;

struct DetectorOccupancyQaTask {
  // configurables for study of occupancy in time windows
  Configurable<bool> confAddBasicQAhistos{"AddBasicQAhistos", true, "0 - add basic histograms, 1 - skip"};
  Configurable<float> confTimeIntervalForOccupancyCalculation{"TimeIntervalForOccupancyCalculation", 100, "Time interval for TPC occupancy calculation, us"};
  Configurable<float> confOccupancyHistCoeffNtracksForOccupancy{"HistCoeffNtracksForOccupancy", 1., "Coefficient for max nTracks in occupancy histos"};
  Configurable<float> confOccupancyHistCoeffNbins2D{"HistCoeffNbins2D", 1., "Coefficient for nBins in occupancy 2D histos"};
  Configurable<float> confOccupancyHistCoeffNbins3D{"HistCoeffNbins3D", 1., "Coefficient for nBins in occupancy 3D histos"};
  Configurable<float> confCoeffMaxNtracksThisEvent{"CoeffMaxNtracksThisEvent", 1., "Coefficient for max nTracks or FT0 ampl in histos in a given event"};
  Configurable<bool> confFlagApplyROFborderCut{"UseROFborderCut", true, "Use ROF border cut for a current event"};
  Configurable<int> confFlagWhichTimeRange{"FlagWhichTimeRange", 0, "Whicn time range for occupancy calculation: 0 - symmetric, 1 - only past, 2 - only future"};
  Configurable<int> confFlagUseGlobalTracks{"FlagUseGlobalTracks", 0, "For small time bins, use global tracks counter instead of ITSTPC tracks"};

  // configuration for small time binning
  Configurable<float> confTimeIntervalForSmallBins{"TimeIntervalForSmallBins", 100, "Time interval for TPC occupancy calculation in small bins, +/-, us"};
  Configurable<int> confNumberOfSmallTimeBins{"nSmallTimeBins", 40, "Number of small time bins"};

  // event and track cuts for given event
  Configurable<float> confCutVertZMinThisEvent{"VzMinThisEvent", -10, "vZ cut for a current event"};
  Configurable<float> confCutVertZMaxThisEvent{"VzMaxThisEvent", 10, "vZ cut for a current event"};
  Configurable<float> confCutPtMinThisEvent{"PtMinThisEvent", 0.2, "pt cut for particles in a current event"};
  Configurable<float> confCutPtMaxThisEvent{"PtMaxThisEvent", 100., "pt cut for particles in a current event"};
  Configurable<float> confCutEtaMinTracksThisEvent{"EtaMinTracksThisEvent", -0.8, "eta cut for particles in a current event"};
  Configurable<float> confCutEtaMaxTracksThisEvent{"EtaMaxTracksThisEvent", 0.8, "eta cut for particles in a current event"};
  Configurable<int> confCutMinTPCcls{"MinNumTPCcls", 70, "min number of TPC clusters for a current event"};

  // config for QA histograms
  Configurable<bool> confAddTracksVsFwdHistos{"AddTracksVsFwdHistos", true, "0 - add histograms, 1 - skip"};
  Configurable<int> nBinsTracks{"nBinsTracks", 400, "N bins in n tracks histo"};
  Configurable<int> nMaxTracks{"nMaxTracks", 8000, "N max in n tracks histo"};
  Configurable<int> nMaxGlobalTracks{"nMaxGlobalTracks", 3000, "N max in n tracks histo"};
  Configurable<int> nBinsMultFwd{"nBinsMultFwd", 400, "N bins in mult fwd histo"};
  Configurable<float> nMaxMultFwd{"nMaxMultFwd", 200000, "N max in mult fwd histo"};

  Configurable<int> nBinsOccupancy{"nBinsOccupancy", 150, "N bins for occupancy axis"};
  Configurable<float> nMaxOccupancy{"nMaxOccupancy", 15000, "N for max of the occupancy axis"};

  uint64_t minGlobalBC = 0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  bool* applySelection = NULL;
  int nBCsPerOrbit = 3564;
  int lastRunNumber = -1;
  int nOrbits;
  double minOrbit;
  int64_t bcSOR = 0;                      // global bc of the start of the first orbit, setting 0 by default for unanchored MC
  int64_t nBCsPerTF = 128 * nBCsPerOrbit; // duration of TF in bcs, should be 128*3564 or 32*3564, setting 128 orbits by default sfor unanchored MC

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

    // histograms for occupancy-in-time-window study
    double kMaxOccup = confOccupancyHistCoeffNtracksForOccupancy;
    double kMaxThisEv = confCoeffMaxNtracksThisEvent;

    AxisSpec axisBC{3601, -0.5, 3600.5, "bc"};
    histos.add("h2D_diff_FoundBC_vs_BC", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_multAbove10", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_multAbove20", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_multAbove50", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_multAbove100", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_hasTOF", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_hasTRD", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_hasTOF_multAbove10", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});
    histos.add("h2D_diff_FoundBC_vs_BC_hasTRD_multAbove10", "", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});

    if (confAddBasicQAhistos) {
      int nMax1D = kMaxThisEv * 8000;
      histos.add("hNumITS567tracksPerCollision", ";n tracks;n events", kTH1D, {{nMax1D, -0.5, nMax1D - 0.5}});
      histos.add("hNumITS567tracksPerCollisionSel", ";n tracks;n events", kTH1D, {{nMax1D, -0.5, nMax1D - 0.5}});
      histos.add("hNumITSTPCtracksPerCollision", ";n tracks;n events", kTH1D, {{nMax1D, -0.5, nMax1D - 0.5}});
      histos.add("hNumITSTPCtracksPerCollisionSel", ";n tracks;n events", kTH1D, {{nMax1D, -0.5, nMax1D - 0.5}});

      histos.add("hNumITS567tracksInTimeWindow", ";n tracks;n events", kTH1D, {{2500, -0.5, 25000.5}});
      histos.add("hNumITSTPCtracksInTimeWindow", ";n tracks;n events", kTH1D, {{2500, -0.5, 25000.5}});
      histos.add("hNumITS567tracksInTimeWindowSel", ";n tracks;n events", kTH1D, {{2500, -0.5, 25000.5}});
      histos.add("hNumITSTPCtracksInTimeWindowSel", ";n tracks;n events", kTH1D, {{2500, -0.5, 25000.5}});

      histos.add("hNumCollInTimeWindow", ";n collisions;n events", kTH1D, {{201, -0.5, 200.5}});
      histos.add("hNumCollInTimeWindowSel", ";n collisions;n events", kTH1D, {{201, -0.5, 200.5}});
      histos.add("hNumCollInTimeWindowSelITSTPC", ";n collisions;n events", kTH1D, {{201, -0.5, 200.5}});
      histos.add("hNumCollInTimeWindowSelIfTOF", ";n collisions;n events", kTH1D, {{201, -0.5, 200.5}});

      histos.add("hNumUniqueBCInTimeWindow", ";n collisions;n events", kTH1D, {{201, -0.5, 200.5}});
    }
    // 2D
    int nBins3D = 80 * confOccupancyHistCoeffNbins3D;
    int nBins3DoccupancyAxis = 100 * confOccupancyHistCoeffNbins3D;

    if (confAddBasicQAhistos) {
      int nBins2D = 200 * confOccupancyHistCoeffNbins2D;

      histos.add("hNumITSTPCtracksInTimeWindow_vs_TracksPerColl", ";n tracks this collision;n tracks in time window", kTH2D, {{nBins2D, 0, kMaxThisEv * 8000}, {nBins2D, 0, kMaxOccup * 25000}});
      histos.add("hNumITSTPCtracksInTimeWindow_vs_TracksPerColl_withoutThisCol", ";n tracks this collision;n tracks in time window", kTH2D, {{nBins2D, 0, kMaxThisEv * 8000}, {nBins2D, 0, kMaxOccup * 25000}});
      histos.add("hNumITS567tracksInTimeWindow_vs_TracksPerColl", ";n tracks this collision;n tracks in time window", kTH2D, {{nBins2D, 0, kMaxThisEv * 8000}, {nBins2D, 0, kMaxOccup * 25000}});
      histos.add("hNumITS567tracksInTimeWindow_vs_TracksPerColl_withoutThisCol", ";n tracks this collision;n tracks in time window", kTH2D, {{nBins2D, 0, kMaxThisEv * 8000}, {nBins2D, 0, kMaxOccup * 25000}});

      histos.add("hNumITS567tracksInTimeWindow_vs_FT0Campl", ";FT0C ampl. sum in time window;n ITS tracks with 5,6,7 hits in time window", kTH2D, {{nBins2D, 0, kMaxOccup * 250000}, {nBins2D, 0, kMaxOccup * 25000}});
      histos.add("hNumITSTPCtracksInTimeWindow_vs_FT0Campl", ";FT0C ampl. sum in time window;n ITS-TPC tracks in time window", kTH2D, {{nBins2D, 0, kMaxOccup * 250000}, {nBins2D, 0, kMaxOccup * 25000}});
      histos.add("hNumITSTPCtracksInTimeWindow_vs_ITS567tracks", ";n ITS tracks with 5,6,7 hits in time window;n ITS-TPC tracks in time window", kTH2D, {{nBins2D, 0, kMaxOccup * 25000}, {nBins2D, 0, kMaxOccup * 25000}});

      histos.add("hNumITS567tracks_vs_FT0Campl_ThisEvent", ";FT0C ampl.;n ITS tracks with 5,6,7 hits", kTH2D, {{nBins2D, 0, kMaxThisEv * 100000}, {nBins2D, 0, kMaxThisEv * 8000}});
      histos.add("hNumITSTPCtracks_vs_FT0Campl_ThisEvent", ";FT0C ampl.;n ITS-TPC tracks", kTH2D, {{nBins2D, 0, kMaxThisEv * 100000}, {nBins2D, 0, kMaxThisEv * 8000}});
      histos.add("hNumITSTPCtracks_vs_ITS567tracks_ThisEvent", ";n ITS tracks with 5,6,7 hits;n ITS-TPC tracks", kTH2D, {{nBins2D, 0, kMaxThisEv * 8000}, {nBins2D, 0, kMaxThisEv * 8000}});

      // 3D
      histos.add("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow", ";n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;FT0C ampl. sum in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 250000}});
      histos.add("hNumITSTPC_vs_ITS567tracksThisCol_vs_ITS567tracksInTimeWindow", ";n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 25000}});
    }
    // nD, time bins to cover the range -confTimeIntervalForSmallBins... +confTimeIntervalForSmallBins (us)
    double timeBinSize = 2 * confTimeIntervalForSmallBins / confNumberOfSmallTimeBins;
    std::vector<double> arrTimeBins;
    for (int i = 0; i < confNumberOfSmallTimeBins + 1; i++)
      arrTimeBins.push_back(-confTimeIntervalForSmallBins + i * timeBinSize);
    const AxisSpec axisTimeBins{arrTimeBins, "#Delta t, #mus"};
    histos.add("occupancyInTimeBins", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
    histos.add("thisEventITStracksInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("thisEventITSTPCtracksInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});

    // save dt information for several first collisions for QA
    histos.add<TH2>("histOccupInTimeBinsQA", ";dt;this coll id", kTH2F, {axisTimeBins, {nCollisionsForTimeBinQA, -0.5, nCollisionsForTimeBinQA - 0.5}});

    // QA of occupancy-based event selection
    histos.add("hOccupancy", "", kTH1D, {{15002, -1.5, 15000.5}});

    if (confAddTracksVsFwdHistos) {
      AxisSpec axisNtracks{nBinsTracks, -0.5, nMaxTracks - 0.5, "n tracks"};
      AxisSpec axisNtracksGlobal{nBinsTracks, -0.5, nMaxGlobalTracks - 0.5, "n tracks"};
      AxisSpec axisMultFw{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd), "mult Fwd"};
      AxisSpec axisOccupancy{nBinsOccupancy, 0., nMaxOccupancy, "occupancy (n ITS tracks weighted)"};

      histos.add("nTracksPV_vs_V0A_kNoHighOccupancyAgressive", "nTracksPV_vs_V0A_kNoHighOccupancyAgressive", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_kNoHighOccupancyStrict", "nTracksPV_vs_V0A_kNoHighOccupancyStrict", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_kNoHighOccupancyMedium", "nTracksPV_vs_V0A_kNoHighOccupancyMedium", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_kNoHighOccupancyRelaxed", "nTracksPV_vs_V0A_kNoHighOccupancyRelaxed", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_kNoHighOccupancyGentle", "nTracksPV_vs_V0A_kNoHighOccupancyGentle", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_kNoCollInTimeRangeStandard", "nTracksPV_vs_V0A_kNoCollInTimeRangeStandard", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_kNoCollInTimeRangeNarrow", "nTracksPV_vs_V0A_kNoCollInTimeRangeNarrow", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_250", "nTracksPV_vs_V0A_occup_0_250", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500", "nTracksPV_vs_V0A_occup_0_500", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_750", "nTracksPV_vs_V0A_occup_0_750", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", "nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", "nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_noOccupSel", "nTracksPV_vs_V0A_noOccupSel", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", "nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_ABOVE_750", "nTracksPV_vs_V0A_occup_ABOVE_750", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_Minus1", "nTracksPV_vs_V0A_occup_Minus1", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_AntiNoCollInTimeRangeStandard", "nTracksPV_vs_V0A_AntiNoCollInTimeRangeStandard", kTH2F, {axisMultFw, axisNtracks});
      histos.add("nTracksPV_vs_V0A_AntiNoCollInTimeRangeNarrow", "nTracksPV_vs_V0A_AntiNoCollInTimeRangeNarrow", kTH2F, {axisMultFw, axisNtracks});

      histos.add("nTracksGlobal_vs_V0A_kNoHighOccupancyAgressive", "nTracksGlobal_vs_V0A_kNoHighOccupancyAgressive", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_kNoHighOccupancyStrict", "nTracksGlobal_vs_V0A_kNoHighOccupancyStrict", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_kNoHighOccupancyMedium", "nTracksGlobal_vs_V0A_kNoHighOccupancyMedium", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_kNoHighOccupancyRelaxed", "nTracksGlobal_vs_V0A_kNoHighOccupancyRelaxed", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_kNoHighOccupancyGentle", "nTracksGlobal_vs_V0A_kNoHighOccupancyGentle", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_V0A_kNoCollInTimeRangeStandard", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_V0A_kNoCollInTimeRangeNarrow", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_250", "nTracksGlobal_vs_V0A_occup_0_250", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500", "nTracksGlobal_vs_V0A_occup_0_500", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_750", "nTracksGlobal_vs_V0A_occup_0_750", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_noOccupSel", "nTracksGlobal_vs_V0A_noOccupSel", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", "nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_ABOVE_750", "nTracksGlobal_vs_V0A_occup_ABOVE_750", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_Minus1", "nTracksGlobal_vs_V0A_occup_Minus1", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeStandard", "nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeStandard", kTH2F, {axisMultFw, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeNarrow", "nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeNarrow", kTH2F, {axisMultFw, axisNtracksGlobal});

      histos.add("nTracksGlobal_vs_nPV_kNoHighOccupancyAgressive", "nTracksGlobal_vs_nPV_kNoHighOccupancyAgressive", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_kNoHighOccupancyStrict", "nTracksGlobal_vs_nPV_kNoHighOccupancyStrict", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_kNoHighOccupancyMedium", "nTracksGlobal_vs_nPV_kNoHighOccupancyMedium", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_kNoHighOccupancyRelaxed", "nTracksGlobal_vs_nPV_kNoHighOccupancyRelaxed", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_kNoHighOccupancyGentle", "nTracksGlobal_vs_nPV_kNoHighOccupancyGentle", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_kNoCollInTimeRangeStandard", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_kNoCollInTimeRangeNarrow", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_250", "nTracksGlobal_vs_nPV_occup_0_250", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500", "nTracksGlobal_vs_nPV_occup_0_500", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_750", "nTracksGlobal_vs_nPV_occup_0_750", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeNarrow", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_noOccupSel", "nTracksGlobal_vs_nPV_noOccupSel", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", "nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_ABOVE_750", "nTracksGlobal_vs_nPV_occup_ABOVE_750", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_Minus1", "nTracksGlobal_vs_nPV_occup_Minus1", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeStandard", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeNarrow", kTH2F, {axisNtracks, axisNtracksGlobal});

      histos.add("nTracksGlobal_vs_nPV_vs_occup_pure", "nTracksGlobal_vs_nPV_vs_occup_pure", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancy});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancy});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancy});
    }
  }

  Preslice<FullTracksIU> perCollision = aod::track::collisionId;

  void processRun3(
    ColEvSels const& cols,
    FullTracksIU const& tracks,
    BCsRun3 const& bcs,
    aod::FT0s const&)
  {
    int runNumber = bcs.iteratorAt(0).runNumber();
    uint32_t nOrbitsPerTF = 128; // 128 in 2022, 32 in 2023
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber; // do it only once
      int64_t tsSOR = 0;
      int64_t tsEOR = 1;

      if (runNumber >= 500000) { // access CCDB for data or anchored MC only
        int64_t ts = bcs.iteratorAt(0).timestamp();

        EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", ts);
        // access orbit-reset timestamp
        auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
        int64_t tsOrbitReset = (*ctpx)[0]; // us
        // access TF duration, start-of-run and end-of-run timestamps from ECS GRP
        std::map<std::string, std::string> metadata;
        metadata["runNumber"] = Form("%d", runNumber);
        auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", ts, metadata);
        nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF;  nOrbitsPerTF=128 in 2022, 32 in 2023
        tsSOR = grpecs->getTimeStart();        // ms
        tsEOR = grpecs->getTimeEnd();          // ms
        // calculate SOR and EOR orbits
        int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
        int64_t orbitEOR = (tsEOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
        // adjust to the nearest TF edge
        orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF + par->fTimeFrameOrbitShift;
        orbitEOR = orbitEOR / nOrbitsPerTF * nOrbitsPerTF + par->fTimeFrameOrbitShift;
        // set nOrbits and minOrbit used for orbit-axis binning
        nOrbits = orbitEOR - orbitSOR;
        minOrbit = orbitSOR;
        // first bc of the first orbit (should coincide with TF start)
        bcSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;
        // duration of TF in bcs
        nBCsPerTF = nOrbitsPerTF * o2::constants::lhc::LHCMaxBunches;
        LOGP(info, "tsOrbitReset={} us, SOR = {} ms, EOR = {} ms, orbitSOR = {}, nBCsPerTF = {}", tsOrbitReset, tsSOR, tsEOR, orbitSOR, nBCsPerTF);
      }
    }

    // vectors with info for the occupancy study
    std::vector<int> vFoundBCindex(cols.size(), -1);                // indices of found bcs
    std::vector<int64_t> vFoundGlobalBC(cols.size(), 0);            // global BCs for collisions
    std::vector<bool> vIsVertexTOFmatched(cols.size(), 0);          // at least one of vertex contributors is matched to TOF
    std::vector<int> vTracksITS567perColl(cols.size(), 0);          // counter of tracks per found bc for occupancy studies
    std::vector<int> vTracksITS567perCollPtEtaCuts(cols.size(), 0); // counter of tracks per found bc for occupancy studies
    std::vector<int> vTracksGlobalPerCollPtEtaCuts(cols.size(), 0); // counter of tracks per found bc for occupancy studies
    std::vector<int> vTracksITSTPCperColl(cols.size(), 0);          // counter of tracks per found bc for occupancy studies
    std::vector<int> vTracksITSTPCperCollPtEtaCuts(cols.size(), 0); // counter of tracks per found bc for occupancy studies
    std::vector<int> vTFids(cols.size(), 0);
    std::vector<bool> vIsFullInfoForOccupancy(cols.size(), 0);

    const double timeWinOccupancyCalcNS = confTimeIntervalForOccupancyCalculation * 1e3; // ns, to be compared with TPC drift time
    const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;

    for (auto& col : cols) {
      const auto& bc = col.foundBC_as<BCsRun3>();

      // count tracks of different types
      int nITS567cls = 0;
      int nITS567clsPtEtaCuts = 0;
      int nGlobalPtEtaCuts = 0;
      int nITSTPCtracks = 0;
      int nITSTPCtracksPtEtaCuts = 0;
      int nTOFtracks = 0;
      int nTRDtracks = 0;
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (auto& track : tracksGrouped) {
        if (!track.isPVContributor()) {
          continue;
        }
        if (track.itsNCls() >= 5)
          nITS567cls++;
        nITSTPCtracks += track.hasITS() && track.hasTPC();
        nTOFtracks += track.hasTOF();
        nTRDtracks += track.hasTRD();

        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;

        if (track.itsNCls() >= 5)
          nITS567clsPtEtaCuts++;

        if (track.tpcNClsFound() < confCutMinTPCcls)
          continue;
        nITSTPCtracksPtEtaCuts += track.hasITS() && track.hasTPC();

        if (track.itsNCls() >= 5)
          nGlobalPtEtaCuts += track.isGlobalTrack();
      }

      int32_t foundBC = bc.globalIndex();
      int32_t colIndex = col.globalIndex();

      vFoundBCindex[colIndex] = foundBC;
      vFoundGlobalBC[colIndex] = bc.globalBC();

      vIsVertexTOFmatched[colIndex] = nTOFtracks > 0;

      vTracksITS567perColl[colIndex] += nITS567cls;
      vTracksITS567perCollPtEtaCuts[colIndex] += nITS567clsPtEtaCuts;
      vTracksGlobalPerCollPtEtaCuts[colIndex] += nGlobalPtEtaCuts;

      vTracksITSTPCperColl[colIndex] += nITSTPCtracks;
      vTracksITSTPCperCollPtEtaCuts[colIndex] += nITSTPCtracksPtEtaCuts;

      // TF ids within a given cols table
      int TFid = (bc.globalBC() - bcSOR) / nBCsPerTF;
      vTFids[colIndex] = TFid;

      // check that this collision has full information inside the time window (taking into account TF borders)
      int64_t bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
      vIsFullInfoForOccupancy[colIndex] = ((bcInTF - 300) * bcNS > timeWinOccupancyCalcNS) && ((nBCsPerTF - 4000 - bcInTF) * bcNS > timeWinOccupancyCalcNS) ? true : false;

      LOGP(debug, "###  check bcInTF cut: colIndex={} bcInTF={} vIsFullInfoForOccupancy={}", colIndex, bcInTF, static_cast<int>(vIsFullInfoForOccupancy[colIndex]));

      // additional QA:
      if (col.selection_bit(kNoTimeFrameBorder) && col.selection_bit(kNoITSROFrameBorder)) {
        auto bcFoundId = bc.globalBC() % 3564;
        auto bcNonFound = col.bc_as<BCsRun3>();
        auto bcNonFoundId = bcNonFound.globalBC() % 3564;
        int64_t diffFoundBC_vs_BC = (int64_t)bcFoundId - (int64_t)bcNonFoundId;
        histos.fill(HIST("h2D_diff_FoundBC_vs_BC"), bcNonFoundId, diffFoundBC_vs_BC);
        if (nITS567cls > 10)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_multAbove10"), bcNonFoundId, diffFoundBC_vs_BC);
        if (nITS567cls > 20)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_multAbove20"), bcNonFoundId, diffFoundBC_vs_BC);
        if (nITS567cls > 50)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_multAbove50"), bcNonFoundId, diffFoundBC_vs_BC);
        if (nITS567cls > 100)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_multAbove100"), bcNonFoundId, diffFoundBC_vs_BC);

        if (nTOFtracks > 0)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_hasTOF"), bcNonFoundId, diffFoundBC_vs_BC);
        if (nTRDtracks > 0)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_hasTRD"), bcNonFoundId, diffFoundBC_vs_BC);

        if (nITS567cls > 10 && nTOFtracks > 0)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_hasTOF_multAbove10"), bcNonFoundId, diffFoundBC_vs_BC);
        if (nITS567cls > 10 && nTRDtracks > 0)
          histos.fill(HIST("h2D_diff_FoundBC_vs_BC_hasTRD_multAbove10"), bcNonFoundId, diffFoundBC_vs_BC);
      }
    }

    // find for each collision all collisions within the defined time window
    std::vector<std::vector<int>> vCollsInTimeWin;
    std::vector<std::vector<float>> vTimeDeltaForColls; // delta time wrt a given collision
    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();
      std::vector<int> vCollsAssocToGivenColl;
      std::vector<float> vCollsTimeDeltaWrtGivenColl;

      // protection against the TF borders
      if (!vIsFullInfoForOccupancy[colIndex]) {
        vCollsInTimeWin.push_back(vCollsAssocToGivenColl);
        vTimeDeltaForColls.push_back(vCollsTimeDeltaWrtGivenColl);
        continue;
      }

      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
      int64_t TFid = (foundGlobalBC - bcSOR) / nBCsPerTF;

      // find all collisions in time window before the current one (start with the current collision)
      int32_t minColIndex = colIndex;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];

        // check if this is still the same TF
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != TFid)
          break;

        float dt = (thisBC - foundGlobalBC) * bcNS; // ns

        if (confFlagWhichTimeRange == 2 && dt < 0)
          break;

        // check if we are within the chosen time range
        if (dt < -timeWinOccupancyCalcNS)
          break;
        vCollsAssocToGivenColl.push_back(minColIndex);
        vCollsTimeDeltaWrtGivenColl.push_back(dt);
        minColIndex--;
      }

      // find all collisions in time window after the current one
      int32_t maxColIndex = colIndex + 1;
      while (maxColIndex < cols.size() && confFlagWhichTimeRange != 1) {
        int64_t thisBC = vFoundGlobalBC[maxColIndex];
        int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTFid != TFid)
          break;

        float dt = (thisBC - foundGlobalBC) * bcNS; // ns

        if (dt > timeWinOccupancyCalcNS)
          break;
        vCollsAssocToGivenColl.push_back(maxColIndex);
        vCollsTimeDeltaWrtGivenColl.push_back(dt);
        maxColIndex++;
      }

      vCollsInTimeWin.push_back(vCollsAssocToGivenColl);
      vTimeDeltaForColls.push_back(vCollsTimeDeltaWrtGivenColl);
    }

    // perform the occupancy calculation in the pre-defined time window
    uint32_t orbitAtCollIndexZero = 0;
    for (auto& col : cols) {
      int32_t colIndex = col.globalIndex();

      // protection against the TF borders
      if (!vIsFullInfoForOccupancy[colIndex])
        continue;

      // cut on vZ for a given collision
      if (col.posZ() < confCutVertZMinThisEvent || col.posZ() > confCutVertZMaxThisEvent)
        continue;

      // skip if collision is close to ROF border
      if (confFlagApplyROFborderCut && !col.selection_bit(kNoITSROFrameBorder))
        continue;

      std::vector<int> vCollsAssocToGivenColl = vCollsInTimeWin[colIndex];
      std::vector<float> vCollsTimeDeltaWrtGivenColl = vTimeDeltaForColls[colIndex];

      LOGP(debug, "  >> vCollsAssocToGivenColl.size={}", vCollsAssocToGivenColl.size());

      int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
      uint32_t orbit = foundGlobalBC / o2::constants::lhc::LHCMaxBunches;
      if (colIndex == 0)
        orbitAtCollIndexZero = orbit;

      int nITS567tracksInTimeWindow = 0;
      int nITSTPCtracksInTimeWindow = 0;
      int nITS567tracksInTimeWindowSel = 0;
      int nITSTPCtracksInTimeWindowSel = 0;

      int nCollInTimeWindow = 0;
      int nCollInTimeWindowSel = 0;
      int nCollInTimeWindowSelITSTPC = 0;
      int nCollInTimeWindowSelIfTOF = 0;
      double multFT0CmainCollision = 0.f;
      double multFT0CInTimeWindow = 0.f;
      map<int64_t, int32_t> mUniqueBC;

      bool sel = col.selection_bit(kIsTriggerTVX);

      // loop over nearby collisions
      for (int iCol = 0; iCol < vCollsAssocToGivenColl.size(); iCol++) {
        int thisColIndex = vCollsAssocToGivenColl[iCol];
        int64_t thisGlobBC = vFoundGlobalBC[thisColIndex];
        float thisColTimeDiff = vCollsTimeDeltaWrtGivenColl[iCol] / 1e3; // ns -> us

        // fill this-event time bins
        if (thisColIndex != colIndex && fabs(thisColTimeDiff) < confTimeIntervalForSmallBins) {
          LOGP(debug, " iCol={}/{}, thisColIndex={}, colIndex={}, thisColTimeDiff={} nITS={}", iCol, vCollsAssocToGivenColl.size(), thisColIndex, colIndex, thisColTimeDiff, vTracksITS567perColl[thisColIndex]);
          histos.fill(HIST("thisEventITStracksInTimeBins"), thisColTimeDiff, vTracksITS567perColl[thisColIndex]);
          // histos.fill(HIST("thisEventITSTPCtracksInTimeBins"), thisColTimeDiff, vTracksITSTPCperColl[thisColIndex]);
        }
        nCollInTimeWindow++;
        nITS567tracksInTimeWindow += vTracksITS567perColl[thisColIndex];
        nITSTPCtracksInTimeWindow += vTracksITSTPCperColl[thisColIndex];

        auto thisBC = bcs.iteratorAt(vFoundBCindex[thisColIndex]);
        bool selThisBCsel = thisBC.selection_bit(kIsTriggerTVX);

        if (sel && selThisBCsel) {
          nCollInTimeWindowSel++;
          nITS567tracksInTimeWindowSel += vTracksITS567perColl[thisColIndex];
          nITSTPCtracksInTimeWindowSel += vTracksITSTPCperColl[thisColIndex];

          mUniqueBC[thisGlobBC] = thisColIndex;
          if (vTracksITSTPCperColl[thisColIndex] >= 2)
            nCollInTimeWindowSelITSTPC++;

          if (vIsVertexTOFmatched[thisColIndex])
            nCollInTimeWindowSelIfTOF++;
        }

        if (thisBC.has_foundFT0()) {
          float multFT0C = thisBC.foundFT0().sumAmpC();
          if (iCol == 0) // the "middle" collision we study
            multFT0CmainCollision = multFT0C;
          multFT0CInTimeWindow += multFT0C;
        }
        LOGP(debug, "### Occupancy in time window study: colIndex={} thisColIndex={} thisGlobBC={} nTrThisCol={} nITSTPCtracksInTimeWindow={} timeDiff={} deltaBC={}", colIndex, thisColIndex, thisGlobBC, vTracksITSTPCperColl[thisColIndex], nITSTPCtracksInTimeWindow, (foundGlobalBC - thisGlobBC) * bcNS, foundGlobalBC - thisGlobBC);
      }

      LOGP(debug, "   --> ### summary: colIndex={}/{} BC={} orbit={} nCollInTimeWindow={} nCollInTimeWindowSel={} nITSTPCtracksInTimeWindow={} ", colIndex, cols.size(), foundGlobalBC, orbit - orbitAtCollIndexZero, nCollInTimeWindow, nCollInTimeWindowSel, nITSTPCtracksInTimeWindow);

      if (confAddBasicQAhistos) {
        histos.get<TH1>(HIST("hNumITS567tracksInTimeWindow"))->Fill(nITS567tracksInTimeWindow);
        histos.get<TH1>(HIST("hNumITSTPCtracksInTimeWindow"))->Fill(nITSTPCtracksInTimeWindow);

        histos.get<TH1>(HIST("hNumITSTPCtracksPerCollision"))->Fill(vTracksITSTPCperColl[colIndex]);
        histos.get<TH1>(HIST("hNumITS567tracksPerCollision"))->Fill(vTracksITS567perColl[colIndex]);

        histos.get<TH2>(HIST("hNumITSTPCtracksInTimeWindow_vs_TracksPerColl"))->Fill(vTracksITSTPCperColl[colIndex], nITSTPCtracksInTimeWindow);
        histos.get<TH2>(HIST("hNumITSTPCtracksInTimeWindow_vs_TracksPerColl_withoutThisCol"))->Fill(vTracksITSTPCperColl[colIndex], nITSTPCtracksInTimeWindow - vTracksITSTPCperColl[colIndex]);

        histos.get<TH2>(HIST("hNumITS567tracksInTimeWindow_vs_TracksPerColl"))->Fill(vTracksITS567perColl[colIndex], nITS567tracksInTimeWindow);
        histos.get<TH2>(HIST("hNumITS567tracksInTimeWindow_vs_TracksPerColl_withoutThisCol"))->Fill(vTracksITS567perColl[colIndex], nITS567tracksInTimeWindow - vTracksITS567perColl[colIndex]);

        histos.get<TH1>(HIST("hNumCollInTimeWindow"))->Fill(nCollInTimeWindow);

        histos.get<TH1>(HIST("hNumUniqueBCInTimeWindow"))->Fill(mUniqueBC.size());

        if (sel) {
          histos.get<TH1>(HIST("hNumITS567tracksInTimeWindowSel"))->Fill(nITS567tracksInTimeWindowSel);
          histos.get<TH1>(HIST("hNumITSTPCtracksInTimeWindowSel"))->Fill(nITSTPCtracksInTimeWindowSel);

          histos.get<TH1>(HIST("hNumITS567tracksPerCollisionSel"))->Fill(vTracksITS567perColl[colIndex]);
          histos.get<TH1>(HIST("hNumITSTPCtracksPerCollisionSel"))->Fill(vTracksITSTPCperCollPtEtaCuts[colIndex]);

          histos.get<TH1>(HIST("hNumCollInTimeWindowSel"))->Fill(nCollInTimeWindowSel);
          histos.get<TH1>(HIST("hNumCollInTimeWindowSelITSTPC"))->Fill(nCollInTimeWindowSelITSTPC);
          histos.get<TH1>(HIST("hNumCollInTimeWindowSelIfTOF"))->Fill(nCollInTimeWindowSelIfTOF);
        }

        // 2D histograms
        histos.get<TH2>(HIST("hNumITS567tracksInTimeWindow_vs_FT0Campl"))->Fill(multFT0CInTimeWindow, nITS567tracksInTimeWindow);
        histos.get<TH2>(HIST("hNumITSTPCtracksInTimeWindow_vs_FT0Campl"))->Fill(multFT0CInTimeWindow, nITSTPCtracksInTimeWindow);
        histos.get<TH2>(HIST("hNumITSTPCtracksInTimeWindow_vs_ITS567tracks"))->Fill(nITS567tracksInTimeWindow, nITSTPCtracksInTimeWindow);

        histos.get<TH2>(HIST("hNumITS567tracks_vs_FT0Campl_ThisEvent"))->Fill(multFT0CmainCollision, vTracksITS567perCollPtEtaCuts[colIndex]);
        histos.get<TH2>(HIST("hNumITSTPCtracks_vs_FT0Campl_ThisEvent"))->Fill(multFT0CmainCollision, vTracksITSTPCperCollPtEtaCuts[colIndex]);
        histos.get<TH2>(HIST("hNumITSTPCtracks_vs_ITS567tracks_ThisEvent"))->Fill(vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex]);

        // 3D histograms: ITS vs ITSTPC in this event vs occupancy from other events
        histos.get<TH3>(HIST("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow"))->Fill(vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex], multFT0CInTimeWindow - multFT0CmainCollision);
        histos.get<TH3>(HIST("hNumITSTPC_vs_ITS567tracksThisCol_vs_ITS567tracksInTimeWindow"))->Fill(vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex], nITS567tracksInTimeWindow - vTracksITS567perColl[colIndex]);
      }
      // loop over time axis in nD histograms:
      for (int iT = 0; iT < histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetNbinsX(); iT++) {
        int nITStrInTimeBin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetBinContent(iT + 1);
        if (nITStrInTimeBin == 0) // no collisions in this dt bin
          continue;
        // int nITSTPCtInTimeBin = histos.get<TH1>(HIST("thisEventITSTPCtracksInTimeBins"))->GetBinContent(iT + 1);

        float dt = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetBinCenter(iT + 1);
        histos.get<THn>(HIST("occupancyInTimeBins"))->Fill(dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nITStrInTimeBin);

        if (counterQAtimeOccupHistos < nCollisionsForTimeBinQA)
          histos.fill(HIST("histOccupInTimeBinsQA"), dt, counterQAtimeOccupHistos + 1, nITStrInTimeBin);
      }
      histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Reset();
      // histos.get<TH1>(HIST("thisEventITSTPCtracksInTimeBins"))->Reset();
      counterQAtimeOccupHistos++;
    } // end of occupancy calculation

    // ### occupancy event selection QA
    for (auto& col : cols) {
      if (!col.sel8()) {
        continue;
      }

      int occupancy = col.trackOccupancyInTimeRange();
      histos.fill(HIST("hOccupancy"), occupancy);

      if (!confAddTracksVsFwdHistos) {
        continue;
      }

      auto multV0A = col.multFV0A();
      // auto multT0A = col.multFT0A();
      // auto multT0C = col.multFT0C();
      int nPV = 0; // col.multNTracksPV();
      int nGlobalTracks = 0;

      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (auto& track : tracksGrouped) {
        if (!track.isPVContributor())
          continue;
        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;
        if (track.itsNCls() < 5)
          continue;
        nPV++;

        if (track.isGlobalTrack() && track.tpcNClsFound() >= confCutMinTPCcls)
          nGlobalTracks++;
      }

      // nPV tracks vs fwd amplitude
      histos.fill(HIST("nTracksPV_vs_V0A_noOccupSel"), multV0A, nPV);
      histos.fill(HIST("nTracksGlobal_vs_V0A_noOccupSel"), multV0A, nGlobalTracks);
      histos.fill(HIST("nTracksGlobal_vs_nPV_noOccupSel"), nPV, nGlobalTracks);

      if (occupancy >= 0)
        histos.fill(HIST("nTracksGlobal_vs_nPV_vs_occup_pure"), nPV, nGlobalTracks, occupancy);

      if (col.selection_bit(o2::aod::evsel::kNoHighOccupancyAgressive)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoHighOccupancyAgressive"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoHighOccupancyAgressive"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoHighOccupancyAgressive"), nPV, nGlobalTracks);
      }
      if (col.selection_bit(o2::aod::evsel::kNoHighOccupancyStrict)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoHighOccupancyStrict"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoHighOccupancyStrict"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoHighOccupancyStrict"), nPV, nGlobalTracks);
      }
      if (col.selection_bit(o2::aod::evsel::kNoHighOccupancyMedium)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoHighOccupancyMedium"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoHighOccupancyMedium"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoHighOccupancyMedium"), nPV, nGlobalTracks);
      }
      if (col.selection_bit(o2::aod::evsel::kNoHighOccupancyRelaxed)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoHighOccupancyRelaxed"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoHighOccupancyRelaxed"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoHighOccupancyRelaxed"), nPV, nGlobalTracks);
      }
      if (col.selection_bit(o2::aod::evsel::kNoHighOccupancyGentle)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoHighOccupancyGentle"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoHighOccupancyGentle"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoHighOccupancyGentle"), nPV, nGlobalTracks);
      }
      if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoCollInTimeRangeStandard"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoCollInTimeRangeStandard"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoCollInTimeRangeStandard"), nPV, nGlobalTracks);
        if (occupancy >= 0)
          histos.fill(HIST("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard"), nPV, nGlobalTracks, occupancy);
      }
      if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoCollInTimeRangeNarrow"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoCollInTimeRangeNarrow"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoCollInTimeRangeNarrow"), nPV, nGlobalTracks);
        if (occupancy >= 0)
          histos.fill(HIST("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow"), nPV, nGlobalTracks, occupancy);
      }
      if (occupancy >= 0 && occupancy < 250) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_0_250"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_0_250"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_250"), nPV, nGlobalTracks);
      }
      if (occupancy >= 0 && occupancy < 500) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_0_500"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_0_500"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_500"), nPV, nGlobalTracks);
      }
      if (occupancy >= 0 && occupancy < 750) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_0_750"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_0_750"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_750"), nPV, nGlobalTracks);
      }
      if (occupancy >= 0 && occupancy < 500 && col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard"), nPV, nGlobalTracks);
      }
      if (occupancy >= 0 && occupancy < 500 && col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeNarrow"), nPV, nGlobalTracks);
      }
      if (occupancy >= 0 && occupancy < 500 && col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) && col.selection_bit(kNoSameBunchPileup) && col.selection_bit(kIsGoodZvtxFT0vsPV)) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard_extraCuts"), nPV, nGlobalTracks);
      }

      // more checks
      if (occupancy >= 750) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_ABOVE_750"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_ABOVE_750"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_ABOVE_750"), nPV, nGlobalTracks);
      }
      if (occupancy == -1) {
        histos.fill(HIST("nTracksPV_vs_V0A_occup_Minus1"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_occup_Minus1"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_Minus1"), nPV, nGlobalTracks);
      }
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        histos.fill(HIST("nTracksPV_vs_V0A_AntiNoCollInTimeRangeStandard"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeStandard"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeStandard"), nPV, nGlobalTracks);
      }
      if (!col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
        histos.fill(HIST("nTracksPV_vs_V0A_AntiNoCollInTimeRangeNarrow"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeNarrow"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeNarrow"), nPV, nGlobalTracks);
      }
    }
  }
  PROCESS_SWITCH(DetectorOccupancyQaTask, processRun3, "Process Run3 tracking vs detector occupancy QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DetectorOccupancyQaTask>(cfgc)};
}
