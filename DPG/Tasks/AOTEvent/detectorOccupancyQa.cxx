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

#include <vector>
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
#include "Common/DataModel/Centrality.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

using BCsRun2 = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::BcSels, aod::Run2MatchedToBCSparse>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
// using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;
// using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA>;

struct DetectorOccupancyQaTask {
  // configurables for study of occupancy in time windows
  Configurable<bool> confAddBasicQAhistos{"AddBasicQAhistos", true, "0 - add basic histograms, 1 - skip"};                                                              // o2-linter: disable=name/configurable
  Configurable<float> confTimeIntervalForOccupancyCalculation{"TimeIntervalForOccupancyCalculation", 100, "Time interval for TPC occupancy calculation, us"};           // o2-linter: disable=name/configurable
  Configurable<float> confOccupancyHistCoeffNtracksForOccupancy{"HistCoeffNtracksForOccupancy", 1., "Coefficient for max nTracks in occupancy histos"};                 // o2-linter: disable=name/configurable
  Configurable<float> confOccupancyHistCoeffNbins2D{"HistCoeffNbins2D", 1., "Coefficient for nBins in occupancy 2D histos"};                                            // o2-linter: disable=name/configurable
  Configurable<float> confOccupancyHistCoeffNbins3D{"HistCoeffNbins3D", 1., "Coefficient for nBins in occupancy 3D histos"};                                            // o2-linter: disable=name/configurable
  Configurable<float> confCoeffMaxNtracksThisEvent{"CoeffMaxNtracksThisEvent", 1., "Coefficient for max nTracks or FT0 ampl in histos in a given event"};               // o2-linter: disable=name/configurable
  Configurable<bool> confFlagApplyROFborderCut{"ApplyROFborderCut", true, "Use ROF border cut for a current event"};                                                    // o2-linter: disable=name/configurable
  Configurable<bool> confFlagApplyTFborderCut{"ApplyTFborderCut", true, "Use TF border cut for a current event"};                                                       // o2-linter: disable=name/configurable
  Configurable<int> confFlagWhichTimeRange{"FlagWhichTimeRange", 0, "Whicn time range for occupancy calculation: 0 - symmetric, 1 - only past, 2 - only future"};       // o2-linter: disable=name/configurable
  Configurable<bool> confFlagUseGlobalTracks{"FlagUseGlobalTracks", false, "For small time bins, use global tracks counter instead of ITSTPC tracks"};                  // o2-linter: disable=name/configurable
  Configurable<bool> confFlagUseNoCollInRofStrict{"FlagUseNoCollInRofStrict", false, "Suppress same-ROF events for occupancy historams"};                               // o2-linter: disable=name/configurable
  Configurable<bool> confFlagUseNoHighMultCollInPrevRof{"FlagUseNoHighMultCollInPrevRof", false, "Suppress high-multiplicity prev-ROF events for occupancy historams"}; // o2-linter: disable=name/configurable

  // configuration for small time binning
  Configurable<float> confTimeIntervalForSmallBins{"TimeIntervalForSmallBins", 100, "Time interval for TPC occupancy calculation in small bins, +/-, us"}; // o2-linter: disable=name/configurable
  Configurable<int> confNumberOfSmallTimeBins{"nSmallTimeBins", 40, "Number of small time bins"};                                                          // o2-linter: disable=name/configurable

  // event and track cuts for given event
  Configurable<float> confCutVertZMinThisEvent{"VzMinThisEvent", -10, "vZ cut for a current event"};                           // o2-linter: disable=name/configurable
  Configurable<float> confCutVertZMaxThisEvent{"VzMaxThisEvent", 10, "vZ cut for a current event"};                            // o2-linter: disable=name/configurable
  Configurable<float> confCutPtMinThisEvent{"PtMinThisEvent", 0.2, "pt cut for particles in a current event"};                 // o2-linter: disable=name/configurable
  Configurable<float> confCutPtMaxThisEvent{"PtMaxThisEvent", 100., "pt cut for particles in a current event"};                // o2-linter: disable=name/configurable
  Configurable<float> confCutEtaMinTracksThisEvent{"EtaMinTracksThisEvent", -0.8, "eta cut for particles in a current event"}; // o2-linter: disable=name/configurable
  Configurable<float> confCutEtaMaxTracksThisEvent{"EtaMaxTracksThisEvent", 0.8, "eta cut for particles in a current event"};  // o2-linter: disable=name/configurable
  Configurable<int> confCutMinTPCcls{"MinNumTPCcls", 70, "min number of TPC clusters for a current event"};                    // o2-linter: disable=name/configurable

  // config for QA histograms
  Configurable<bool> confAddTracksVsFwdHistos{"AddTracksVsFwdHistos", true, "0 - add histograms, 1 - skip"}; // o2-linter: disable=name/configurable
  Configurable<int> nBinsTracks{"nBinsTracks", 400, "N bins in n tracks histo"};                             // o2-linter: disable=name/configurable
  Configurable<int> nMaxTracks{"nMaxTracks", 8000, "N max in n tracks histo"};                               // o2-linter: disable=name/configurable
  Configurable<int> nMaxGlobalTracks{"nMaxGlobalTracks", 3000, "N max in n tracks histo"};                   // o2-linter: disable=name/configurable
  Configurable<int> nBinsMultFwd{"nBinsMultFwd", 400, "N bins in mult fwd histo"};                           // o2-linter: disable=name/configurable
  Configurable<float> nMaxMultFwd{"nMaxMultFwd", 200000, "N max in mult fwd histo"};                         // o2-linter: disable=name/configurable

  Configurable<int> nBinsOccupancy{"nBinsOccupancy", 150, "N bins for occupancy axis"};         // o2-linter: disable=name/configurable
  Configurable<float> nMaxOccupancy{"nMaxOccupancy", 15000, "N for max of the occupancy axis"}; // o2-linter: disable=name/configurable

  Configurable<int> nMaxBcInTFforAnalysis{"nMaxBcInTFforAnalysis", -1, "When to stop taking collisions in TF, if -1: take all collisions"}; // o2-linter: disable=name/configurable

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

    const AxisSpec axisBCinTF{static_cast<int>(nBCsPerTF), 0, static_cast<double>(nBCsPerTF), "bc in TF"};
    histos.add("hNcolVsBcInTF", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
    histos.add("hNcolVsBcInTFafterMaxBcCut", ";bc in TF; n collisions", kTH1F, {axisBCinTF});

    // histograms for occupancy-in-time-window study
    double kMaxOccup = confOccupancyHistCoeffNtracksForOccupancy;
    double kMaxThisEv = confCoeffMaxNtracksThisEvent;

    // 1D, dE/dx, etc.
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
      histos.add("hNumCollInTimeWindowVsOrbit", ";orbit id;n collisions;n events", kTH2F, {{128, -0.5, 127.5}, {201, -0.5, 200.5}});

      histos.add("hNumUniqueBCInTimeWindow", ";n collisions;n events", kTH1D, {{201, -0.5, 200.5}});

      // dE/dx
      histos.add("dEdx_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}});
      histos.add("dEdx_vs_Momentum_occupBelow200", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}});
      histos.add("dEdx_vs_Momentum_occupBelow200_kNoCollStd", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}});
      histos.add("dEdx_vs_Momentum_occupAbove4000", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}});
      AxisSpec axisBinsOccupStudydEdx{{0., 500, 1000, 2000, 4000, 6000, 8000, 15000}, "p_{T}"};
      histos.add("dEdx_vs_Momentum_vs_occup", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}, axisBinsOccupStudydEdx});
      histos.add("dEdx_vs_Momentum_vs_occup_eta_02_04", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}, axisBinsOccupStudydEdx});
      histos.add("dEdx_vs_Momentum_vs_occup_eta_04_02", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}, axisBinsOccupStudydEdx});

      histos.add("dEdx_vs_centr_vs_occup_narrow_p_win", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, {60, 0, 15000, "occupancy"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}});

      // ### kinematic distributions for events with high occupancy at specified dt ranges
      histos.add("track_distr_nITStrThisEv_10_200/hEventCount", ";delta-time bin id;n events", kTH1D, {{5, -0.5, 4.5}});
      histos.add("track_distr_nITStrThisEv_above_2000/hEventCount", ";delta-time bin id;n events", kTH1D, {{5, -0.5, 4.5}});

      const int nEtaBins = 800;
      histos.add("track_distr_nITStrThisEv_10_200/hEta_lowOccupInTPC", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInRecentPast", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInCloseFuture", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInDistantFuture", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInNeighbourEvents", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});

      histos.add("track_distr_nITStrThisEv_above_2000/hEta_lowOccupInTPC", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInRecentPast", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInCloseFuture", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInDistantFuture", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInNeighbourEvents", ";#eta;n tracks", kTH1D, {{nEtaBins, -1.0, 1.0}});

      const int nPhiBins = 800;
      AxisSpec axisPhi{nPhiBins, 0, TMath::TwoPi(), "#varphi"}; // o2-linter: disable=external-pi
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents", ";#varphi;n tracks", kTH1D, {axisPhi});

      histos.add("track_distr_nITStrThisEv_above_2000/hPhi_lowOccupInTPC", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInRecentPast", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInCloseFuture", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInDistantFuture", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInNeighbourEvents", ";#varphi;n tracks", kTH1D, {axisPhi});

      // const int nPtBins = 800;
      AxisSpec axisLogPt{200, 0.05, 40, "p_{T}"};
      axisLogPt.makeLogarithmic();
      histos.add("track_distr_nITStrThisEv_10_200/hPt_lowOccupInTPC", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_10_200/hPt_highOccupInRecentPast", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_10_200/hPt_highOccupInCloseFuture", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_10_200/hPt_highOccupInDistantFuture", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_10_200/hPt_highOccupInNeighbourEvents", ";p_{T};n tracks", kTH1D, {axisLogPt});

      histos.add("track_distr_nITStrThisEv_above_2000/hPt_lowOccupInTPC", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_above_2000/hPt_highOccupInRecentPast", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_above_2000/hPt_highOccupInCloseFuture", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_above_2000/hPt_highOccupInDistantFuture", ";p_{T};n tracks", kTH1D, {axisLogPt});
      histos.add("track_distr_nITStrThisEv_above_2000/hPt_highOccupInNeighbourEvents", ";p_{T};n tracks", kTH1D, {axisLogPt});

      // 3D: pt vs centr vs occup
      histos.add("ptGlobal_vs_centr_vs_occup", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, {60, 0, 15000, "occupancy"}, axisLogPt});
      histos.add("ptPV_vs_centr_vs_occup", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, {60, 0, 15000, "occupancy"}, axisLogPt});
      histos.add("ptGlobal_vs_centr_vs_occup_NoCollStd", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, {60, 0, 15000, "occupancy"}, axisLogPt});
      histos.add("ptPV_vs_centr_vs_occup_NoCollStd", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, {60, 0, 15000, "occupancy"}, axisLogPt});
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
      histos.add("hNumITSTPC_vs_ITS567tracksThisCol_vs_ITS567tracksInTimeWindow_BEFORE_sel", ";n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 25000}});
      histos.add("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow_BEFORE_sel", ";n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;FT0C ampl. sum in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 250000}});

      histos.add("hNumITSTPC_vs_ITS567tracksThisCol_vs_ITS567tracksInTimeWindow", ";n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 25000}});
      histos.add("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow", ";n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;FT0C ampl. sum in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 250000}});
      histos.add("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow_kNoCollInTimeRangeNarrow", ";n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;FT0C ampl. sum in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 8000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 250000}});

      histos.add("hNumITSTPC_vs_FT0CthisCol_vs_FT0CamplInTimeWindow_kNoCollInTimeRangeNarrow", ";FT0C this collision;n ITS-TPC tracks, this collision;FT0C ampl. sum in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 80000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 250000}});
      histos.add("hNumITS567_vs_FT0CthisCol_vs_FT0CamplInTimeWindow_kNoCollInTimeRangeNarrow", ";FT0C this collision;n ITS567cls tracks, this collision;FT0C ampl. sum in time window", kTH3D, {{nBins3D, 0, kMaxThisEv * 80000}, {nBins3D, 0, kMaxThisEv * 8000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 250000}});
    }
    // nD, time bins to cover the range -confTimeIntervalForSmallBins... +confTimeIntervalForSmallBins (us)
    double timeBinSize = 2 * confTimeIntervalForSmallBins / confNumberOfSmallTimeBins;
    std::vector<double> arrTimeBins;
    for (int i = 0; i < confNumberOfSmallTimeBins + 1; i++)
      arrTimeBins.push_back(-confTimeIntervalForSmallBins + i * timeBinSize);
    const AxisSpec axisTimeBins{arrTimeBins, "#Delta t, #mus"};
    int nBinsX = 20;
    int nBinsY = 40;
    histos.add("occupancyInTimeBins_BEFORE_sel", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
    histos.add("occupancyInTimeBins_occupByFT0_BEFORE_sel", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});

    histos.add("occupancyInTimeBins", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
    histos.add("occupancyInTimeBins_occupByFT0", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});
    histos.add("occupancyInTimeBins_occupByFT0_kNoCollInTimeRangeNarrow", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});

    histos.add("occupancyInTimeBins_vs_FT0thisCol_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
    histos.add("occupancyInTimeBins_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});

    histos.add("occupancyInTimeBins_nITS567_vs_FT0thisCol_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS567cls tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
    histos.add("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS567cls tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});
    histos.add("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow_NoCollInRofStrict", ";time bin (#mus);FT0C this collision, this collision;n ITS567cls tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});

    histos.add("thisEventITStracksInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("thisEventITSTPCtracksInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("thisEventFT0CInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});

    histos.add("qaForHighOccupITStracksInTimeBinPast", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("qaForHighOccupITStracksInTimeBinFuture1", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("qaForHighOccupITStracksInTimeBinFuture2", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("qaForHighOccupITStracksForNeighbourEvents", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});

    // save dt information for several first collisions for QA
    histos.add<TH2>("histOccupInTimeBinsQA", ";dt;this coll id", kTH2F, {axisTimeBins, {nCollisionsForTimeBinQA, -0.5, nCollisionsForTimeBinQA - 0.5}});

    // QA of occupancy-based event selection
    histos.add("hOccupancy", "", kTH1D, {{15002, -1.5, 15000.5}});
    histos.add("hOccupancyVsOrbit", ";orbit id;weighted occupancy;n events", kTH2F, {{128, -0.5, 127.5}, {600, 0, 15000}});

    AxisSpec axisOccupancyTracks{nBinsOccupancy, 0., nMaxOccupancy, "occupancy (n ITS tracks weighted)"};
    AxisSpec axisCentrality{100, 0, 100, "centrality, %"};
    histos.add("hCentrVsOccupancy", "hCentrVsOccupancy", kTH2F, {axisCentrality, axisOccupancyTracks});
    histos.add("hCentrVsOccupancyNoCollStd", "hCentrVsOccupancyNoCollStd", kTH2F, {axisCentrality, axisOccupancyTracks});

    if (confAddTracksVsFwdHistos) {
      AxisSpec axisNtracks{nBinsTracks, -0.5, nMaxTracks - 0.5, "n tracks"};
      AxisSpec axisNtracksGlobal{nBinsTracks, -0.5, nMaxGlobalTracks - 0.5, "n tracks"};
      AxisSpec axisMultV0A{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd), "mult V0A"};
      AxisSpec axisMultFT0C{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd * 0.4), "mult FT0C"};

      histos.add("nTracksPV_vs_V0A_kNoCollInTimeRangeStandard", "nTracksPV_vs_V0A_kNoCollInTimeRangeStandard", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_kNoCollInTimeRangeNarrow", "nTracksPV_vs_V0A_kNoCollInTimeRangeNarrow", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_250", "nTracksPV_vs_V0A_occup_0_250", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500", "nTracksPV_vs_V0A_occup_0_500", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_750", "nTracksPV_vs_V0A_occup_0_750", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", "nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", "nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_noOccupSel", "nTracksPV_vs_V0A_noOccupSel", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", "nTracksPV_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_ABOVE_750", "nTracksPV_vs_V0A_occup_ABOVE_750", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_occup_Minus1", "nTracksPV_vs_V0A_occup_Minus1", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_AntiNoCollInTimeRangeStandard", "nTracksPV_vs_V0A_AntiNoCollInTimeRangeStandard", kTH2F, {axisMultV0A, axisNtracks});
      histos.add("nTracksPV_vs_V0A_AntiNoCollInTimeRangeNarrow", "nTracksPV_vs_V0A_AntiNoCollInTimeRangeNarrow", kTH2F, {axisMultV0A, axisNtracks});

      histos.add("nTracksGlobal_vs_V0A_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_V0A_kNoCollInTimeRangeStandard", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_V0A_kNoCollInTimeRangeNarrow", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_250", "nTracksGlobal_vs_V0A_occup_0_250", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500", "nTracksGlobal_vs_V0A_occup_0_500", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_750", "nTracksGlobal_vs_V0A_occup_0_750", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeNarrow", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_noOccupSel", "nTracksGlobal_vs_V0A_noOccupSel", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", "nTracksGlobal_vs_V0A_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_ABOVE_750", "nTracksGlobal_vs_V0A_occup_ABOVE_750", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_occup_Minus1", "nTracksGlobal_vs_V0A_occup_Minus1", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeStandard", "nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeStandard", kTH2F, {axisMultV0A, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeNarrow", "nTracksGlobal_vs_V0A_AntiNoCollInTimeRangeNarrow", kTH2F, {axisMultV0A, axisNtracksGlobal});

      histos.add("nTracksGlobal_vs_nPV_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_kNoCollInTimeRangeStandard", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_kNoCollInTimeRangeNarrow", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_250", "nTracksGlobal_vs_nPV_occup_0_250", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500", "nTracksGlobal_vs_nPV_occup_0_500", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_750", "nTracksGlobal_vs_nPV_occup_0_750", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_2000", "nTracksGlobal_vs_nPV_occup_0_2000", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeNarrow", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_noOccupSel", "nTracksGlobal_vs_nPV_noOccupSel", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", "nTracksGlobal_vs_nPV_occup_0_500_kNoCollInTimeRangeStandard_extraCuts", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_0_2000_kNoCollInTimeRangeStandard_extraCuts", "nTracksGlobal_vs_nPV_occup_0_2000_kNoCollInTimeRangeStandard_extraCuts", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_ABOVE_750", "nTracksGlobal_vs_nPV_occup_ABOVE_750", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occup_Minus1", "nTracksGlobal_vs_nPV_occup_Minus1", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeStandard", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_AntiNoCollInTimeRangeNarrow", kTH2F, {axisNtracks, axisNtracksGlobal});

      histos.add("nTracksGlobal_vs_nPV_QA_onlyVzCut_noTFROFborderCuts", "nTracksGlobal_vs_nPV_QA_onlyVzCut_noTFROFborderCuts", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_QA_after_TFborderCut", "nTracksGlobal_vs_nPV_QA_after_TFborderCut", kTH2F, {axisNtracks, axisNtracksGlobal});

      histos.add("nTracksGlobal_vs_nPV_occupByFT0C_0_2500", "nTracksGlobal_vs_nPV_occupByFT0C_0_2500", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occupByFT0C_0_20000", "nTracksGlobal_vs_nPV_occupByFT0C_0_20000", kTH2F, {axisNtracks, axisNtracksGlobal});

      // 3D histograms with occupancy axis
      histos.add("nTracksGlobal_vs_nPV_vs_occup_pure", "nTracksGlobal_vs_nPV_vs_occup_pure", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard_extraCuts", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});

      // 3D histograms: nGlobalTracks with cls567 as y-axis, V0A as x-axis:
      histos.add("nTracksGlobal_vs_V0A_vs_occup_pure", "", kTH3F, {axisMultV0A, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_V0A_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "", kTH3F, {axisMultV0A, axisNtracksGlobal, axisOccupancyTracks});
      // FT0C as x-axis:
      histos.add("nTracksGlobal_vs_FT0C_vs_occup_pure", "", kTH3F, {axisMultFT0C, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_FT0C_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "", kTH3F, {axisMultFT0C, axisNtracksGlobal, axisOccupancyTracks});

      // 3D histograms: now - nITStracks with cls567 as y-axis, V0A as x-axis:
      histos.add("nPV_vs_V0A_vs_occup_pure", "", kTH3F, {axisMultV0A, axisNtracks, axisOccupancyTracks});
      histos.add("nPV_vs_V0A_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "", kTH3F, {axisMultV0A, axisNtracks, axisOccupancyTracks});
      // FT0C as x-axis:
      histos.add("nPV_vs_FT0C_vs_occup_pure", "", kTH3F, {axisMultFT0C, axisNtracks, axisOccupancyTracks});
      histos.add("nPV_vs_FT0C_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "", kTH3F, {axisMultFT0C, axisNtracks, axisOccupancyTracks});
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
    if (runNumber != lastRunNumber) {
      lastRunNumber = runNumber; // do it only once

      if (runNumber >= 500000) {
        auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), runNumber);
        // first bc of the first orbit
        bcSOR = runInfo.orbitSOR * o2::constants::lhc::LHCMaxBunches;
        // duration of TF in bcs
        nBCsPerTF = runInfo.orbitsPerTF * o2::constants::lhc::LHCMaxBunches;

        LOGP(info, "bcSOR = {}, nBCsPerTF = {}", bcSOR, nBCsPerTF);
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
    std::vector<float> vAmpFT0CperColl(cols.size(), 0);             // amplitude FT0C per collision

    std::vector<int> vTFids(cols.size(), 0);
    std::vector<bool> vIsFullInfoForOccupancy(cols.size(), 0);
    std::vector<bool> vIsMarkedCollForAnalysis(cols.size(), 0); // cut on the max bcId in the time frame

    std::vector<int> vFlagsForEtaQAvsOccupancyInDeltaTimeWins(cols.size(), 0);

    const double timeWinOccupancyCalcNS = confTimeIntervalForOccupancyCalculation * 1e3; // ns, to be compared with TPC drift time
    const double bcNS = o2::constants::lhc::LHCBunchSpacingNS;

    for (const auto& col : cols) {
      const auto& bc = col.foundBC_as<BCsRun3>();

      // count tracks of different types
      int nITS567cls = 0;
      int nITS567clsPtEtaCuts = 0;
      int nGlobalPtEtaCuts = 0;
      int nITSTPCtracks = 0;
      int nITSTPCtracksPtEtaCuts = 0;
      int nTOFtracks = 0;
      // int nTRDtracks = 0;
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (const auto& track : tracksGrouped) {
        if (!track.isPVContributor()) {
          continue;
        }
        if (track.itsNCls() >= 5)
          nITS567cls++;
        nITSTPCtracks += track.hasITS() && track.hasTPC();
        nTOFtracks += track.hasTOF();
        // nTRDtracks += track.hasTRD();

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

      if (bc.has_foundFT0())
        vAmpFT0CperColl[colIndex] = bc.foundFT0().sumAmpC();

      vIsVertexTOFmatched[colIndex] = nTOFtracks > 0;

      vTracksITS567perColl[colIndex] += nITS567cls;
      vTracksITS567perCollPtEtaCuts[colIndex] += nITS567clsPtEtaCuts;
      vTracksGlobalPerCollPtEtaCuts[colIndex] += nGlobalPtEtaCuts;

      vTracksITSTPCperColl[colIndex] += nITSTPCtracks;
      vTracksITSTPCperCollPtEtaCuts[colIndex] += nITSTPCtracksPtEtaCuts;

      // TF ids within a given cols table
      int tfId = (bc.globalBC() - bcSOR) / nBCsPerTF;
      vTFids[colIndex] = tfId;

      // check that this collision has full information inside the time window (taking into account TF borders)
      int64_t bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
      vIsFullInfoForOccupancy[colIndex] = ((bcInTF - 300) * bcNS > timeWinOccupancyCalcNS) && ((nBCsPerTF - 4000 - bcInTF) * bcNS > timeWinOccupancyCalcNS) ? true : false;

      // cut on the max bc in the time frame
      vIsMarkedCollForAnalysis[colIndex] = nMaxBcInTFforAnalysis == -1 ? 1 : (bcInTF >= 300 && bcInTF < nMaxBcInTFforAnalysis);
      LOGP(debug, "###  check bcInTF cut: colIndex={} bcInTF={} vIsFullInfoForOccupancy={}", colIndex, bcInTF, static_cast<int>(vIsFullInfoForOccupancy[colIndex]));
    }

    // find for each collision all collisions within the defined time window
    std::vector<std::vector<int>> vCollsInTimeWin;
    std::vector<std::vector<float>> vTimeDeltaForColls; // delta time wrt a given collision
    for (const auto& col : cols) {
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
      int64_t tfId = (foundGlobalBC - bcSOR) / nBCsPerTF;

      // find all collisions in time window before the current one (start with the current collision)
      int32_t minColIndex = colIndex;
      while (minColIndex >= 0) {
        int64_t thisBC = vFoundGlobalBC[minColIndex];

        // check if this is still the same TF
        int64_t thisTfId = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTfId != tfId)
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
        int64_t thisTfId = (thisBC - bcSOR) / nBCsPerTF;
        if (thisTfId != tfId)
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
    for (const auto& col : cols) {
      int32_t colIndex = col.globalIndex();

      // protection against the TF borders
      if (!vIsFullInfoForOccupancy[colIndex])
        continue;

      // cut on the max bcId in the time frame (to avoid the artificial fade-out tail in the MC productions)
      if (!vIsMarkedCollForAnalysis[colIndex])
        continue;

      // cut on vZ for a given collision
      if (col.posZ() < confCutVertZMinThisEvent || col.posZ() > confCutVertZMaxThisEvent)
        continue;

      // skip if collision is close to TF border
      if (confFlagApplyTFborderCut && !col.selection_bit(kNoTimeFrameBorder))
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
      for (unsigned int iCol = 0; iCol < vCollsAssocToGivenColl.size(); iCol++) {
        int thisColIndex = vCollsAssocToGivenColl[iCol];
        int64_t thisGlobBC = vFoundGlobalBC[thisColIndex];
        float thisColTimeDiff = vCollsTimeDeltaWrtGivenColl[iCol] / 1e3; // ns -> us

        // fill this-event time bins
        if (thisColIndex != colIndex && std::fabs(thisColTimeDiff) < confTimeIntervalForSmallBins) {
          LOGP(debug, " iCol={}/{}, thisColIndex={}, colIndex={}, thisColTimeDiff={} nITS={}", iCol, vCollsAssocToGivenColl.size(), thisColIndex, colIndex, thisColTimeDiff, vTracksITS567perColl[thisColIndex]);
          histos.fill(HIST("thisEventITStracksInTimeBins"), thisColTimeDiff, vTracksITS567perColl[thisColIndex]);
          histos.fill(HIST("thisEventFT0CInTimeBins"), thisColTimeDiff, vAmpFT0CperColl[thisColIndex]);
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
        histos.fill(HIST("hNumITS567tracksInTimeWindow"), nITS567tracksInTimeWindow);
        histos.fill(HIST("hNumITSTPCtracksInTimeWindow"), nITSTPCtracksInTimeWindow);

        histos.fill(HIST("hNumITSTPCtracksPerCollision"), vTracksITSTPCperColl[colIndex]);
        histos.fill(HIST("hNumITS567tracksPerCollision"), vTracksITS567perColl[colIndex]);

        histos.fill(HIST("hNumITSTPCtracksInTimeWindow_vs_TracksPerColl"), vTracksITSTPCperColl[colIndex], nITSTPCtracksInTimeWindow);
        histos.fill(HIST("hNumITSTPCtracksInTimeWindow_vs_TracksPerColl_withoutThisCol"), vTracksITSTPCperColl[colIndex], nITSTPCtracksInTimeWindow - vTracksITSTPCperColl[colIndex]);

        histos.fill(HIST("hNumITS567tracksInTimeWindow_vs_TracksPerColl"), vTracksITS567perColl[colIndex], nITS567tracksInTimeWindow);
        histos.fill(HIST("hNumITS567tracksInTimeWindow_vs_TracksPerColl_withoutThisCol"), vTracksITS567perColl[colIndex], nITS567tracksInTimeWindow - vTracksITS567perColl[colIndex]);

        histos.fill(HIST("hNumCollInTimeWindow"), nCollInTimeWindow);

        int64_t bcInTF = (vFoundGlobalBC[colIndex] - bcSOR) % nBCsPerTF;
        int orbitId = bcInTF / o2::constants::lhc::LHCMaxBunches;
        histos.fill(HIST("hNumCollInTimeWindowVsOrbit"), orbitId, nCollInTimeWindow);

        histos.fill(HIST("hNumUniqueBCInTimeWindow"), mUniqueBC.size());

        // 3D before ev quality cut:
        histos.fill(HIST("hNumITSTPC_vs_ITS567tracksThisCol_vs_ITS567tracksInTimeWindow_BEFORE_sel"), vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex], nITS567tracksInTimeWindow - vTracksITS567perColl[colIndex]);
        histos.fill(HIST("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow_BEFORE_sel"), vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex], multFT0CInTimeWindow - multFT0CmainCollision);

        if (sel && std::fabs(col.posZ()) < 10) {
          histos.fill(HIST("hNumITS567tracksInTimeWindowSel"), nITS567tracksInTimeWindowSel);
          histos.fill(HIST("hNumITSTPCtracksInTimeWindowSel"), nITSTPCtracksInTimeWindowSel);

          histos.fill(HIST("hNumITS567tracksPerCollisionSel"), vTracksITS567perColl[colIndex]);
          histos.fill(HIST("hNumITSTPCtracksPerCollisionSel"), vTracksITSTPCperCollPtEtaCuts[colIndex]);

          histos.fill(HIST("hNumCollInTimeWindowSel"), nCollInTimeWindowSel);
          histos.fill(HIST("hNumCollInTimeWindowSelITSTPC"), nCollInTimeWindowSelITSTPC);
          histos.fill(HIST("hNumCollInTimeWindowSelIfTOF"), nCollInTimeWindowSelIfTOF);

          // 3D histograms: ITS vs ITSTPC in this event vs occupancy from other events
          histos.fill(HIST("hNumITSTPC_vs_ITS567tracksThisCol_vs_ITS567tracksInTimeWindow"), vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex], nITS567tracksInTimeWindow - vTracksITS567perColl[colIndex]);
          histos.fill(HIST("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow"), vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex], multFT0CInTimeWindow - multFT0CmainCollision);
          if (col.selection_bit(kNoCollInTimeRangeNarrow)) {
            histos.fill(HIST("hNumITSTPC_vs_ITS567tracksThisCol_vs_FT0CamplInTimeWindow_kNoCollInTimeRangeNarrow"), vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex], multFT0CInTimeWindow - multFT0CmainCollision);

            histos.fill(HIST("hNumITSTPC_vs_FT0CthisCol_vs_FT0CamplInTimeWindow_kNoCollInTimeRangeNarrow"), multFT0CmainCollision, vTracksITSTPCperCollPtEtaCuts[colIndex], multFT0CInTimeWindow - multFT0CmainCollision);
            histos.fill(HIST("hNumITS567_vs_FT0CthisCol_vs_FT0CamplInTimeWindow_kNoCollInTimeRangeNarrow"), multFT0CmainCollision, vTracksITS567perCollPtEtaCuts[colIndex], multFT0CInTimeWindow - multFT0CmainCollision);
          }
        }

        // 2D histograms
        histos.fill(HIST("hNumITS567tracksInTimeWindow_vs_FT0Campl"), multFT0CInTimeWindow, nITS567tracksInTimeWindow);
        histos.fill(HIST("hNumITSTPCtracksInTimeWindow_vs_FT0Campl"), multFT0CInTimeWindow, nITSTPCtracksInTimeWindow);
        histos.fill(HIST("hNumITSTPCtracksInTimeWindow_vs_ITS567tracks"), nITS567tracksInTimeWindow, nITSTPCtracksInTimeWindow);

        histos.fill(HIST("hNumITS567tracks_vs_FT0Campl_ThisEvent"), multFT0CmainCollision, vTracksITS567perCollPtEtaCuts[colIndex]);
        histos.fill(HIST("hNumITSTPCtracks_vs_FT0Campl_ThisEvent"), multFT0CmainCollision, vTracksITSTPCperCollPtEtaCuts[colIndex]);
        histos.fill(HIST("hNumITSTPCtracks_vs_ITS567tracks_ThisEvent"), vTracksITS567perCollPtEtaCuts[colIndex], vTracksITSTPCperCollPtEtaCuts[colIndex]);
      }

      // counters of occupancy in specified delta-time ranges, to monitor eta, phi, pt distributions later
      float integralFullDeltaTime = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral();
      int binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(-39.5); // us
      int binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(-10.5);
      float integralPast = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);
      binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(20.5);
      binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(49.5);
      float integralFuture1 = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);
      binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(50.5);
      binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(79.5);
      float integralFuture2 = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);
      binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(-9.5);
      binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(19.5);
      float integralNeighbourEvents = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);

      // recent past
      if (integralFullDeltaTime < 150) // ~empty detector
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 1;
      // recent past
      if (integralPast > /*3000*/ 2500 && (integralFullDeltaTime - integralPast) < 120) // low occupancy outside the dt region of interest
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 2;
      // close future
      if (integralFuture1 > /*3000*/ 2500 && (integralFullDeltaTime - integralFuture1) < 120)
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 3;
      // distant future
      if (integralFuture2 > /*3000*/ 2500 && (integralFullDeltaTime - integralFuture2) < 120)
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 4;
      // neighbour events
      if (integralNeighbourEvents > /*3000*/ 2500 && (integralFullDeltaTime - integralNeighbourEvents) < 120)
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 5;

      // loop over time axis in nD histograms:
      for (int iT = 0; iT < histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetNbinsX(); iT++) {
        int nITStrInTimeBin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetBinContent(iT + 1);
        if (nITStrInTimeBin == 0) // no collisions in this dt bin
          continue;
        // int nITSTPCtInTimeBin = histos.get<TH1>(HIST("thisEventITSTPCtracksInTimeBins"))->GetBinContent(iT + 1);

        float dt = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetBinCenter(iT + 1);

        int nFT0CInTimeBin = histos.get<TH1>(HIST("thisEventFT0CInTimeBins"))->GetBinContent(iT + 1);

        histos.fill(HIST("occupancyInTimeBins_BEFORE_sel"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nITStrInTimeBin);
        histos.fill(HIST("occupancyInTimeBins_occupByFT0_BEFORE_sel"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);

        bool flagFillOccupVsDt = true;
        if (confFlagUseNoCollInRofStrict && !col.selection_bit(kNoCollInRofStrict))
          flagFillOccupVsDt = false;
        if (confFlagUseNoHighMultCollInPrevRof && !col.selection_bit(kNoHighMultCollInPrevRof))
          flagFillOccupVsDt = false;

        if (sel && std::fabs(col.posZ()) < 10 && flagFillOccupVsDt) {
          histos.fill(HIST("occupancyInTimeBins"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nITStrInTimeBin);
          histos.fill(HIST("occupancyInTimeBins_occupByFT0"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);

          if (col.selection_bit(kNoCollInTimeRangeNarrow)) {
            histos.fill(HIST("occupancyInTimeBins_occupByFT0_kNoCollInTimeRangeNarrow"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);
            histos.fill(HIST("occupancyInTimeBins_vs_FT0thisCol_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nITStrInTimeBin);
            histos.fill(HIST("occupancyInTimeBins_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);

            histos.fill(HIST("occupancyInTimeBins_nITS567_vs_FT0thisCol_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], vTracksITS567perCollPtEtaCuts[colIndex], nITStrInTimeBin);
            histos.fill(HIST("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], vTracksITS567perCollPtEtaCuts[colIndex], nFT0CInTimeBin);
            if (col.selection_bit(kNoCollInRofStrict))
              histos.fill(HIST("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow_NoCollInRofStrict"), dt, vAmpFT0CperColl[colIndex], vTracksITS567perCollPtEtaCuts[colIndex], nFT0CInTimeBin);
          }
        }

        //

        if (counterQAtimeOccupHistos < nCollisionsForTimeBinQA)
          histos.fill(HIST("histOccupInTimeBinsQA"), dt, counterQAtimeOccupHistos + 1, nITStrInTimeBin);

        // QA for high occup in time bins
        if (vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] == 2)
          histos.fill(HIST("qaForHighOccupITStracksInTimeBinPast"), dt, nITStrInTimeBin);
        if (vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] == 3)
          histos.fill(HIST("qaForHighOccupITStracksInTimeBinFuture1"), dt, nITStrInTimeBin);
        if (vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] == 4)
          histos.fill(HIST("qaForHighOccupITStracksInTimeBinFuture2"), dt, nITStrInTimeBin);
        if (vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] == 5)
          histos.fill(HIST("qaForHighOccupITStracksForNeighbourEvents"), dt, nITStrInTimeBin);
      }

      // reset delta time hist for this event
      histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Reset();
      // histos.get<TH1>(HIST("thisEventITSTPCtracksInTimeBins"))->Reset();
      histos.get<TH1>(HIST("thisEventFT0CInTimeBins"))->Reset();
      counterQAtimeOccupHistos++;
    } // end of occupancy calculation

    // ### occupancy event selection QA
    for (const auto& col : cols) {
      // if (!col.sel8()) {
      //   continue;
      // }
      if (!col.selection_bit(kIsTriggerTVX))
        continue;
      // cut on vZ for a given collision
      if (col.posZ() < confCutVertZMinThisEvent || col.posZ() > confCutVertZMaxThisEvent)
        continue;

      int32_t colIndex = col.globalIndex();
      int64_t bcInTF = (vFoundGlobalBC[colIndex] - bcSOR) % nBCsPerTF;
      histos.fill(HIST("hNcolVsBcInTF"), bcInTF);

      // cut on the max bcId in the time frame (to avoid the artificial fade-out tail in the MC productions)
      if (!vIsMarkedCollForAnalysis[colIndex])
        continue;

      histos.fill(HIST("hNcolVsBcInTFafterMaxBcCut"), bcInTF);

      auto multV0A = col.multFV0A();
      // auto multT0A = col.multFT0A();
      auto multT0C = col.multFT0C();
      int nPV = 0; // col.multNTracksPV();
      int nGlobalTracks = 0;

      int occupancy = col.trackOccupancyInTimeRange();

      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (const auto& track : tracksGrouped) {
        if (!track.isPVContributor())
          continue;
        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;
        if (track.itsNCls() < 5)
          continue;
        nPV++;

        if (track.isGlobalTrack() && track.tpcNClsFound() >= confCutMinTPCcls) {
          nGlobalTracks++;

          if (track.passedTPCRefit()) {
            float signedP = track.sign() * track.tpcInnerParam();
            histos.fill(HIST("dEdx_vs_Momentum"), signedP, track.tpcSignal());
            if (occupancy >= 0 && occupancy < 200) {
              histos.fill(HIST("dEdx_vs_Momentum_occupBelow200"), signedP, track.tpcSignal());
              if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
                histos.fill(HIST("dEdx_vs_Momentum_occupBelow200_kNoCollStd"), signedP, track.tpcSignal());
            }
            if (occupancy > 4000)
              histos.fill(HIST("dEdx_vs_Momentum_occupAbove4000"), signedP, track.tpcSignal());

            if (occupancy >= 0) {
              histos.fill(HIST("dEdx_vs_Momentum_vs_occup"), signedP, track.tpcSignal(), occupancy);

              if (track.eta() > 0.2 && track.eta() < 0.4)
                histos.fill(HIST("dEdx_vs_Momentum_vs_occup_eta_02_04"), signedP, track.tpcSignal(), occupancy);
              if (track.eta() > -0.4 && track.eta() < -0.2)
                histos.fill(HIST("dEdx_vs_Momentum_vs_occup_eta_04_02"), signedP, track.tpcSignal(), occupancy);

              // dE/dx in narrow mom bin vs centrality and occupancy
              if (signedP > 0.38 && signedP < 0.4)
                histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win"), nPV, occupancy, track.tpcSignal());
            }
          }
        }
      }

      if (confAddTracksVsFwdHistos)
        histos.fill(HIST("nTracksGlobal_vs_nPV_QA_onlyVzCut_noTFROFborderCuts"), nPV, nGlobalTracks);

      // skip if collision is close to TF border
      if (confFlagApplyTFborderCut && !col.selection_bit(kNoTimeFrameBorder))
        continue;

      if (confAddTracksVsFwdHistos)
        histos.fill(HIST("nTracksGlobal_vs_nPV_QA_after_TFborderCut"), nPV, nGlobalTracks);

      // skip if collision is close to ROF border
      if (confFlagApplyROFborderCut && !col.selection_bit(kNoITSROFrameBorder))
        continue;

      histos.fill(HIST("hOccupancy"), occupancy);
      if (occupancy >= 0) {
        int orbitId = bcInTF / o2::constants::lhc::LHCMaxBunches;
        histos.fill(HIST("hOccupancyVsOrbit"), orbitId, occupancy);
      }

      // another track loop to fill track-level histograms
      if (confAddBasicQAhistos) {
        int flagWhichDeltaTimeWin = vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex];
        bool flagNoCollNearby = col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard);

        if (occupancy >= 0) {
          if (nPV >= 10 && nPV < 200) {
            if (flagNoCollNearby && flagWhichDeltaTimeWin != 5)
              histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEventCount"), flagWhichDeltaTimeWin);
            if (flagWhichDeltaTimeWin == 5) // nearby collisions --> avoid checking the flagNoCollNearby flag
              histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEventCount"), flagWhichDeltaTimeWin);
          }
          if (nPV >= 2000) {
            if (flagNoCollNearby && flagWhichDeltaTimeWin != 5)
              histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEventCount"), flagWhichDeltaTimeWin);
            if (flagWhichDeltaTimeWin == 5) // nearby collisions --> avoid checking the flagNoCollNearby flag
              histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEventCount"), flagWhichDeltaTimeWin);
          }
        }

        for (const auto& track : tracksGrouped) {
          if (!track.isPVContributor())
            continue;
          if (track.itsNCls() < 5)
            continue;
          if (!(track.isGlobalTrack() && track.tpcNClsFound() >= confCutMinTPCcls))
            continue;

          // pt vs centr vs occup
          if (occupancy >= 0) {
            histos.fill(HIST("ptGlobal_vs_centr_vs_occup"), nPV, occupancy, track.pt());
            histos.fill(HIST("ptPV_vs_centr_vs_occup"), nPV, occupancy, track.pt());
            if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
              histos.fill(HIST("ptGlobal_vs_centr_vs_occup_NoCollStd"), nPV, occupancy, track.pt());
              histos.fill(HIST("ptPV_vs_centr_vs_occup_NoCollStd"), nPV, occupancy, track.pt());
            }

            if (nPV >= 10 && nPV < 200) {
              if (flagWhichDeltaTimeWin == 1 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_lowOccupInTPC"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_lowOccupInTPC"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 2 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInRecentPast"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInRecentPast"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 3 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInCloseFuture"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInCloseFuture"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 4 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInDistantFuture"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInDistantFuture"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 5) {
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInNeighbourEvents"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInNeighbourEvents"), track.pt());
              }
            } else if (nPV >= 2000) {
              if (flagWhichDeltaTimeWin == 1 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_lowOccupInTPC"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_lowOccupInTPC"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_lowOccupInTPC"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 2 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInRecentPast"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInRecentPast"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInRecentPast"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 3 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInCloseFuture"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInCloseFuture"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInCloseFuture"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 4 && flagNoCollNearby) {
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInDistantFuture"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInDistantFuture"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInDistantFuture"), track.pt());
              }
              if (flagWhichDeltaTimeWin == 5) {
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInNeighbourEvents"), track.eta());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInNeighbourEvents"), track.phi());
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInNeighbourEvents"), track.pt());
              }
            }
          } // end of if (occupancy >= 0)
        }
      } // end of spec track loop to fill track histograms

      // occupancy vs centrality
      auto t0cCentr = col.centFT0C();
      if (occupancy >= 0) {
        histos.fill(HIST("hCentrVsOccupancy"), t0cCentr, occupancy);
        if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
          histos.fill(HIST("hCentrVsOccupancyNoCollStd"), t0cCentr, occupancy);
      }

      if (!confAddTracksVsFwdHistos) {
        continue;
      }

      // nPV tracks vs fwd amplitude
      histos.fill(HIST("nTracksPV_vs_V0A_noOccupSel"), multV0A, nPV);
      histos.fill(HIST("nTracksGlobal_vs_V0A_noOccupSel"), multV0A, nGlobalTracks);
      histos.fill(HIST("nTracksGlobal_vs_nPV_noOccupSel"), nPV, nGlobalTracks);

      if (occupancy >= 0) {
        histos.fill(HIST("nTracksGlobal_vs_nPV_vs_occup_pure"), nPV, nGlobalTracks, occupancy);
        histos.fill(HIST("nTracksGlobal_vs_V0A_vs_occup_pure"), multV0A, nGlobalTracks, occupancy);
        histos.fill(HIST("nTracksGlobal_vs_FT0C_vs_occup_pure"), multT0C, nGlobalTracks, occupancy);

        histos.fill(HIST("nPV_vs_V0A_vs_occup_pure"), multV0A, nPV, occupancy);
        histos.fill(HIST("nPV_vs_FT0C_vs_occup_pure"), multT0C, nPV, occupancy);
      }

      if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        histos.fill(HIST("nTracksPV_vs_V0A_kNoCollInTimeRangeStandard"), multV0A, nPV);
        histos.fill(HIST("nTracksGlobal_vs_V0A_kNoCollInTimeRangeStandard"), multV0A, nGlobalTracks);
        histos.fill(HIST("nTracksGlobal_vs_nPV_kNoCollInTimeRangeStandard"), nPV, nGlobalTracks);
        if (occupancy >= 0)
          histos.fill(HIST("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard"), nPV, nGlobalTracks, occupancy);
        if (occupancy >= 0 && col.selection_bit(kNoSameBunchPileup) && col.selection_bit(kIsGoodZvtxFT0vsPV)) {
          histos.fill(HIST("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard_extraCuts"), nPV, nGlobalTracks, occupancy);
          histos.fill(HIST("nTracksGlobal_vs_V0A_vs_occup_kNoCollInTimeRangeStandard_extraCuts"), multV0A, nGlobalTracks, occupancy);
          histos.fill(HIST("nTracksGlobal_vs_FT0C_vs_occup_kNoCollInTimeRangeStandard_extraCuts"), multT0C, nGlobalTracks, occupancy);

          histos.fill(HIST("nPV_vs_V0A_vs_occup_kNoCollInTimeRangeStandard_extraCuts"), multV0A, nPV, occupancy);
          histos.fill(HIST("nPV_vs_FT0C_vs_occup_kNoCollInTimeRangeStandard_extraCuts"), multT0C, nPV, occupancy);
        }
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
      if (occupancy >= 0 && occupancy < 2000) {
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_2000"), nPV, nGlobalTracks);
      }
      // ### now vs FT0C occupancy:
      float occupByFT0C = col.ft0cOccupancyInTimeRange();
      if (occupByFT0C >= 0 && occupByFT0C < 2500) {
        histos.fill(HIST("nTracksGlobal_vs_nPV_occupByFT0C_0_2500"), nPV, nGlobalTracks);
      }
      if (occupByFT0C >= 0 && occupByFT0C < 20000) {
        histos.fill(HIST("nTracksGlobal_vs_nPV_occupByFT0C_0_20000"), nPV, nGlobalTracks);
      }
      //
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
      if (occupancy >= 0 && occupancy < 2000 && col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) && col.selection_bit(kNoSameBunchPileup) && col.selection_bit(kIsGoodZvtxFT0vsPV)) {
        histos.fill(HIST("nTracksGlobal_vs_nPV_occup_0_2000_kNoCollInTimeRangeStandard_extraCuts"), nPV, nGlobalTracks);
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
