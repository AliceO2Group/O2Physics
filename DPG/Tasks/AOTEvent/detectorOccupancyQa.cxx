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
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/BunchFilling.h"
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
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA>;

struct DetectorOccupancyQaTask {
  // configurables for study of occupancy in time windows
  Configurable<bool> confAddBasicQAhistos{"FlagAddBasicQAhistos", true, "0 - add basic histograms, 1 - skip"};                                                // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confAddTimeDependentHistos{"FlagAddTimeDependentHistos", true, "0 - add time-dependent histograms, 1 - skip"};                           // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confTimeIntervalForOccupancyCalculation{"TimeIntervalForOccupancyCalculation", 100, "Time interval for TPC occupancy calculation, us"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confOccupancyHistCoeffNtracksForOccupancy{"HistCoeffNtracksForOccupancy", 1., "Coefficient for max nTracks in occupancy histos"};       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confOccupancyHistCoeffNbins2D{"HistCoeffNbins2D", 1., "Coefficient for nBins in occupancy 2D histos"};                                  // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confOccupancyHistCoeffNbins3D{"HistCoeffNbins3D", 1., "Coefficient for nBins in occupancy 3D histos"};                                  // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCoeffMaxNtracksThisEvent{"CoeffMaxNtracksThisEvent", 1., "Coefficient for max nTracks or FT0 ampl in histos in a given event"};     // o2-linter: disable=name/configurable (temporary fix)
  // Configurable<bool> confFlagApplyROFborderCut{"ApplyROFborderCut", true, "Use ROF border cut for a current event"};                                                    // o2-linter: disable=name/configurable (temporary fix)
  // Configurable<bool> confFlagApplyTFborderCut{"ApplyTFborderCut", true, "Use TF border cut for a current event"};                                                       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confFlagWhichTimeRange{"FlagWhichTimeRange", 0, "Whicn time range for occupancy calculation: 0 - symmetric, 1 - only past, 2 - only future"};       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confFlagUseGlobalTracks{"FlagUseGlobalTracks", false, "For small time bins, use global tracks counter instead of ITSTPC tracks"};                  // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confFlagUseNoCollInRofStrict{"FlagUseNoCollInRofStrict", false, "Suppress same-ROF events for occupancy historams"};                               // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confFlagUseNoHighMultCollInPrevRof{"FlagUseNoHighMultCollInPrevRof", false, "Suppress high-multiplicity prev-ROF events for occupancy historams"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confFlagCentralityIsAvailable{"FlagCentralityIsAvailable", true, "Fill centrality-related historams"};                                             // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confFlagManyHeavyHistos{"FlagManyHeavyHistos", true, "Fill more TH2, TH3, THn historams"};                                                         // o2-linter: disable=name/configurable (temporary fix)
  Configurable<bool> confFlagIsTOFIsTRDdtStudy{"FlagIsTOFIsTRDdtStudy", false, "Fill THn dt historams with isTOF and isTRD condition"};                                 // o2-linter: disable=name/configurable (temporary fix)

  // configuration for small time binning
  Configurable<float> confTimeIntervalForSmallBins{"TimeIntervalForSmallBins", 100, "Time interval for TPC occupancy calculation in small bins, +/-, us"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confNumberOfSmallTimeBins{"nSmallTimeBins", 40, "Number of small time bins"};                                                          // o2-linter: disable=name/configurable (temporary fix)

  // event and track cuts for given event
  Configurable<float> confCutVertZMinThisEvent{"VzMinThisEvent", -10, "vZ cut for a current event"};                           // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutVertZMaxThisEvent{"VzMaxThisEvent", 10, "vZ cut for a current event"};                            // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutPtMinThisEvent{"PtMinThisEvent", 0.2, "pt cut for particles in a current event"};                 // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutPtMaxThisEvent{"PtMaxThisEvent", 100., "pt cut for particles in a current event"};                // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutEtaMinTracksThisEvent{"EtaMinTracksThisEvent", -0.8, "eta cut for particles in a current event"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confCutEtaMaxTracksThisEvent{"EtaMaxTracksThisEvent", 0.8, "eta cut for particles in a current event"};  // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confCutMinTPCcls{"MinNumTPCcls", 50, "min number of TPC clusters for a current event"};                    // o2-linter: disable=name/configurable (temporary fix)

  // config for QA histograms
  Configurable<bool> confAddTracksVsFwdHistos{"FlagAddTracksVsFwdHistos", true, "0 - add histograms, 1 - skip"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nBinsTracks{"nBinsTracks", 400, "N bins in n tracks histo"};                                 // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nMaxTracks{"nMaxTracks", 8000, "N max in n tracks histo"};                                   // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nMaxGlobalTracks{"nMaxGlobalTracks", 3000, "N max in n tracks histo"};                       // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> nBinsMultFwd{"nBinsMultFwd", 400, "N bins in mult fwd histo"};                               // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> nMaxMultFwd{"nMaxMultFwd", 200000, "N max in mult fwd histo"};                             // o2-linter: disable=name/configurable (temporary fix)

  Configurable<int> nBinsOccupancy{"nBinsOccupancy", 150, "N bins for occupancy axis"};         // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> nMaxOccupancy{"nMaxOccupancy", 15000, "N for max of the occupancy axis"}; // o2-linter: disable=name/configurable (temporary fix)

  Configurable<int> nMaxBcInTFforAnalysis{"nMaxBcInTFforAnalysis", -1, "When to stop taking collisions in TF, if -1: take all collisions"}; // o2-linter: disable=name/configurable (temporary fix)

  Configurable<int> confNPhiBins{"nPhiBins", 810, "N phi bits for histograms"}; // o2-linter: disable=name/configurable (temporary fix)

  ConfigurableAxis confAxisPtBinsForPhiStudy{"PtBinsForPhiStudy", {VARIABLE_WIDTH, 0.2, 0.6, 1.0, 2.0, 10}, "pt axis"};
  ConfigurableAxis confAxisOccupForKine{"AxisOccupForKine", {VARIABLE_WIDTH, 0, 500, 1000, 2000, 4000, 6000, 8000, 10000, 20000}, "weighted occupancy"};

  Configurable<bool> confUsePhiAtTPCinnerR{"UsePhiAtTPCinnerR", false, "0 - not use, 1 - use"};                              // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confUseAorCsideForPhiStudy{"UseAorCsideForPhiStudy", -1, "-1 - use full eta range, 0 - A, 1 - C sides"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<float> confRadiusForPhiCorrection{"RadiusForPhiCorrection", 0.8, "default: inner TPC radius, cm"};            // o2-linter: disable=name/configurable (temporary fix)

  Configurable<int> confApplyGoodITSstavesFlaginEvSel{"ApplyGoodITSstavesFlaginEvSel", 0, "0 - no, 1 - yes"}; // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confMinITSclsPerTrack{"MinITSclsPerTrack", 5, "should be in 4..7"};                       // o2-linter: disable=name/configurable (temporary fix)

  Configurable<std::vector<float>> confTimeSlicesForPastFutureStudies{"TimeSlicesForPastFutureStudies", {-40, -10, 20, 50, 80}, "Time slices for past/future studies, us"};

  // configuration for THnD multi-dim histo(s):
  Configurable<bool> confFlagFillTHn{"FlagFillTHn", false, "Fill THn historams for multi-dim QA"};                                 // o2-linter: disable=name/configurable (temporary fix)
  Configurable<int> confTHnAxis_nPhiBins{"THn_nPhiBins", 180, "nPhiBins"};                                                         // o2-linter: disable=name/configurable (temporary fix)
  ConfigurableAxis confTHnAxis_R{"THn_R", {8, -0.5f, 7.5f}, "ids of radii"};                                                       // o2-linter: disable=name/configurable (temporary fix)
  ConfigurableAxis confTHnAxis_qOp{"THn_qOp", {16, -4.f, 4.f}, "qOp"};                                                             // o2-linter: disable=name/configurable (temporary fix)
  ConfigurableAxis confTHnAxis_IR{"THn_IR", {VARIABLE_WIDTH, 0, 12, 25, 38, 50}, "IR, kHz"};                                       // o2-linter: disable=name/configurable (temporary fix)
  ConfigurableAxis confTHnAxis_occ{"THn_occupancy", {VARIABLE_WIDTH, 0, 500, 1000, 2000, 4000, 6000, 8000}, "weighted occupancy"}; // o2-linter: disable=name/configurable (temporary fix)
  ConfigurableAxis confTHnAxis_centr{"THn_centr", {VARIABLE_WIDTH, 0, 500, 1000, 2000, 4000}, "centrality by nPVtracks"};          // o2-linter: disable=name/configurable (temporary fix)
  ConfigurableAxis confTHnAxis_eta{"THn_eta", {8, -0.8f, 0.8f}, "eta"};                                                            // o2-linter: disable=name/configurable (temporary fix)

  uint64_t minGlobalBC = 0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  bool* applySelection = NULL;
  int nBCsPerOrbit = 3564;
  int lastRunNumber = -1;
  int nOrbits;
  double minOrbit;
  int64_t bcSOR = 0;                     // global bc of the start of the first orbit, setting 0 by default for unanchored MC
  int64_t nBCsPerTF = 32 * nBCsPerOrbit; // duration of TF in bcs, should be 128*3564 or 32*3564, setting 128 orbits by default sfor unanchored MC
  ctpRateFetcher mRateFetcher;

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
    histos.add("hNcolVsBcInTF/hNcolVsBcInTF", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
    histos.add("hNcolVsBcInTF/hNcolVsBcInTF_vertexTOFmatched", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
    histos.add("hNcolVsBcInTF/hNcolVsBcInTFafterMaxBcCut", ";bc in TF; n collisions", kTH1F, {axisBCinTF});

    // QA of occupancy-based event selection
    histos.add("hOccupancy", "", kTH1D, {{15002, -1.5, 15000.5}});

    AxisSpec axisOccupancyTracks{nBinsOccupancy, 0., nMaxOccupancy, "occupancy (n ITS tracks weighted)"};
    if (confFlagCentralityIsAvailable) {
      AxisSpec axisCentrality{100, 0, 100, "centrality, %"};
      histos.add("hCentrVsOccupancy", "hCentrVsOccupancy", kTH2F, {axisCentrality, axisOccupancyTracks});
      histos.add("hCentrVsOccupancyNoCollStd", "hCentrVsOccupancyNoCollStd", kTH2F, {axisCentrality, axisOccupancyTracks});
    }
    // track QA counters
    histos.add("nTrackCounter_after_cuts_QA", "", kTH1D, {{12, -0.5, 11.5, "track QA"}});
    TAxis* axTrackCounters = reinterpret_cast<TAxis*>(histos.get<TH1>(HIST("nTrackCounter_after_cuts_QA"))->GetXaxis());
    axTrackCounters->SetBinLabel(1, "all");
    axTrackCounters->SetBinLabel(2, "PVcontrib");
    axTrackCounters->SetBinLabel(3, "ptCut");
    axTrackCounters->SetBinLabel(4, "etaCut");
    axTrackCounters->SetBinLabel(5, "itsNCls>=5");
    axTrackCounters->SetBinLabel(6, "isGlobal,nTPCcls>=70");
    axTrackCounters->SetBinLabel(7, "passedTPCRefit");
    axTrackCounters->SetBinLabel(8, "occupancy>=0");
    axTrackCounters->SetBinLabel(9, "fracton nClsNoPID (0,0.8)");
    axTrackCounters->SetBinLabel(10, "pos");
    axTrackCounters->SetBinLabel(11, "neg");

    // histograms for occupancy-in-time-window study
    double kMaxOccup = confOccupancyHistCoeffNtracksForOccupancy;
    double kMaxThisEv = confCoeffMaxNtracksThisEvent;

    // 1D, dE/dx, etc.
    if (confAddBasicQAhistos) {
      histos.add("hOccupancyVsOrbit", ";orbit id;weighted occupancy;n events", kTH2F, {{128, -0.5, 127.5}, {600, 0, 15000}});

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
      AxisSpec axisDeDx{800, 0.0, 800.0, "dE/dx (a. u.)"};
      histos.add("dEdx_vs_Momentum", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
      histos.add("dEdx_vs_Momentum_occupBelow200", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
      histos.add("dEdx_vs_Momentum_occupBelow200_kNoCollStd", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
      histos.add("dEdx_vs_Momentum_occupAbove4000", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
      histos.add("dEdx_vs_Momentum_NegativeFractionNclsPID", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
      histos.add("dEdx_vs_Momentum_HighFractionNclsNonPID", "dE/dx", kTH2F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx});
      AxisSpec axisBinsOccupStudydEdx{{0., 500, 1000, 2000, 4000, 6000, 8000, 15000}, "p_{T}"};
      histos.add("dEdx_vs_Momentum_vs_occup", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
      if (confFlagManyHeavyHistos) {
        histos.add("dEdx_vs_Momentum_vs_occup_eta_02_04", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
        histos.add("dEdx_vs_Momentum_vs_occup_eta_04_02", "dE/dx", kTH3F, {{1000, -5.0, 5.0, "#it{p}/Z (GeV/c)"}, axisDeDx, axisBinsOccupStudydEdx});
      }

      AxisSpec axisOccupancyForDeDxStudies{60, 0, 15000, "occupancy"};
      histos.add("dEdx_vs_centr_vs_occup_narrow_p_win", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
      histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_pos", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
      histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_neg", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
      histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_pos_FractionPIDclsInRange", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});
      histos.add("dEdx_vs_centr_vs_occup_narrow_p_win_neg_FractionPIDclsInRange", "dE/dx", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisDeDx});

      AxisSpec axisNTPCcls{160, 0, 160, "n TPC clusters"};
      histos.add("tpcNClsFound_vs_centr_vs_occup", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisNTPCcls});
      histos.add("tpcNClsFindable_vs_centr_vs_occup", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisNTPCcls});
      histos.add("tpcNClsShared_vs_centr_vs_occup", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisNTPCcls});
      histos.add("tpcNClsShared_vs_centr_vs_occup_Aside", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisNTPCcls});
      histos.add("tpcNClsShared_vs_centr_vs_occup_Cside", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisNTPCcls});
      histos.add("tpcNClsShared_vs_centr_vs_occup_pos", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisNTPCcls});
      histos.add("tpcNClsShared_vs_centr_vs_occup_neg", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisNTPCcls});

      // July 2025: more for data vs MC
      AxisSpec axisChi2TPC{150, 0, 15, "chi2Ncl TPC"};
      histos.add("QA_noTPCcuts/nPV_10_200/tpcNClsFindable_vs_occup_pt_02_05", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcNClsFindable_vs_occup_pt_05_10", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcNClsFindable_vs_occup_pt_above1_0", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcNClsFound_vs_occup_pt_02_05", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcNClsFound_vs_occup_pt_05_10", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcNClsFound_vs_occup_pt_above1_0", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcChi2NCl_vs_occup_pt_02_05", "", kTH2F, {axisChi2TPC, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcChi2NCl_vs_occup_pt_05_10", "", kTH2F, {axisChi2TPC, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_10_200/tpcChi2NCl_vs_occup_pt_above1_0", "", kTH2F, {axisChi2TPC, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcNClsFindable_vs_occup_pt_02_05", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcNClsFindable_vs_occup_pt_05_10", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcNClsFindable_vs_occup_pt_above1_0", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcNClsFound_vs_occup_pt_02_05", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcNClsFound_vs_occup_pt_05_10", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcNClsFound_vs_occup_pt_above1_0", "", kTH2F, {axisNTPCcls, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcChi2NCl_vs_occup_pt_02_05", "", kTH2F, {axisChi2TPC, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcChi2NCl_vs_occup_pt_05_10", "", kTH2F, {axisChi2TPC, axisOccupancyForDeDxStudies});
      histos.add("QA_noTPCcuts/nPV_above2000/tpcChi2NCl_vs_occup_pt_above1_0", "", kTH2F, {axisChi2TPC, axisOccupancyForDeDxStudies});

      AxisSpec axisFractionNclsFindableMinusPID{110, -1.1, 1.1, "TPC nClsFindableMinusPID / nClsFindable"};
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_peripheralByV0A", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_centralByV0A", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_eta02", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_pos", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_neg", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_lowPt", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});
      histos.add("fraction_tpcNClsFindableMinusPID_vs_occup_highPt", "", kTH2D, {axisOccupancyForDeDxStudies, axisFractionNclsFindableMinusPID});

      // more QA for TPC cls counting
      histos.add("tpcNClsFindable", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsFindableMinusFound", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsFindableMinusCrossedRows", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsShared", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsFindableMinusPID", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClUsedForPID", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsFound", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsFoundAsDiffByHand", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsFindableMinusPID_CORRECTED", "", kTH1D, {{601, -300.5, 300.5}});
      histos.add("tpcNClsFoundMinusPID_BY_HAND", "", kTH1D, {{601, -300.5, 300.5}});

      histos.add("tpcNClsUsedForPID_vs_Findable", ";tpcNClsFindable;tpcNClUsedForPID", kTH2D, {{601, -300.5, 300.5}, {601, -300.5, 300.5}});
      histos.add("tpcNClsUsedForPID_vs_Findable_CORRECTED", ";tpcNClsFindable;tpcNClUsedForPID", kTH2D, {{601, -300.5, 300.5}, {601, -300.5, 300.5}});
      histos.add("tpcNClsShared_vs_Findable", ";tpcNClsFindable;tpcNClsShared", kTH2D, {{601, -300.5, 300.5}, {601, -300.5, 300.5}});
      histos.add("tpcNClsFound_vs_Findable", ";tpcNClsFindable;tpcNClsFound", kTH2D, {{601, -300.5, 300.5}, {601, -300.5, 300.5}});
      histos.add("tpcNClsUsedForPID_vs_Shared", ";tpcNClsShared;tpcNClUsedForPID", kTH2D, {{601, -300.5, 300.5}, {601, -300.5, 300.5}});
      histos.add("tpcNClsUsedForPID_vs_Found", ";tpcNClsFound;tpcNClUsedForPID", kTH2D, {{601, -300.5, 300.5}, {601, -300.5, 300.5}});

      // ### kinematic distributions for events with high occupancy at specified dt ranges
      histos.add("track_distr_nITStrThisEv_10_200/hEventCount", ";delta-time bin id;n events", kTH1D, {{5, -0.5, 4.5}});
      histos.add("track_distr_nITStrThisEv_above_2000/hEventCount", ";delta-time bin id;n events", kTH1D, {{5, -0.5, 4.5}});

      const int nEtaBins = 800;
      AxisSpec axisEta{nEtaBins, -1.0, 1.0, "#eta"}; // o2-linter: disable=external-pi (temporary fix)
      histos.add("track_distr_nITStrThisEv_10_200/hEta_lowOccupInTPC", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInRecentPast", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInCloseFuture", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInDistantFuture", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_10_200/hEta_highOccupInNeighbourEvents", ";#eta;n tracks", kTH1D, {axisEta});

      histos.add("track_distr_nITStrThisEv_above_2000/hEta_lowOccupInTPC", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInRecentPast", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInCloseFuture", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInDistantFuture", ";#eta;n tracks", kTH1D, {axisEta});
      histos.add("track_distr_nITStrThisEv_above_2000/hEta_highOccupInNeighbourEvents", ";#eta;n tracks", kTH1D, {axisEta});

      const int nPhiBins = confNPhiBins;                        // 810=18*45
      AxisSpec axisPhi{nPhiBins, 0, TMath::TwoPi(), "#varphi"}; // o2-linter: disable=external-pi (temporary fix)
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture", ";#varphi;n tracks", kTH1D, {axisPhi});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents", ";#varphi;n tracks", kTH1D, {axisPhi});

      histos.add("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC_pos_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast_pos_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture_pos_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture_pos_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents_pos_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});

      histos.add("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC_neg_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast_neg_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture_neg_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture_neg_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents_neg_vs_pt", ";#varphi;n tracks", kTH2D, {axisPhi, confAxisPtBinsForPhiStudy});

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

      // July 2025: to compare data and MC (pt, eta, phi)
      // AxisSpec confAxisOccupForKine{{0, 500, 1000, 2000, 4000, 6000, 20000}, "weighted occupancy"};
      // AxisSpec confAxisOccupForKine{{0, 500, 1000, 2000, 4000, 6000, 8000, 10000, 20000}, "weighted occupancy"};
      // AxisSpec confAxisPtBinsForPhiStudy{{0.2, 0.6, 1.0, 2.0, 10}, "pt bins for phi study"};
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPt_pos", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPt_neg", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hEta_pos", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hEta_neg", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_pos", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_neg", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});

      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFindable_pos", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFound_pos", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsCrossedRows_pos", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFindable_neg", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFound_neg", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsCrossedRows_neg", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});

      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPt_pos", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPt_neg", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hEta_pos", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hEta_neg", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPhi_pos", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPhi_neg", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});

      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPt_pos", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPt_neg", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hEta_pos", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hEta_neg", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_pos", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_neg", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_posInitialQA", ";#varphi;n tracks", kTH1D, {{3 * 810, -TMath::TwoPi(), 2 * TMath::TwoPi(), "#varphi"}});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_posModifiedQA", ";#varphi;n tracks", kTH1D, {{3 * 810, -TMath::TwoPi(), 2 * TMath::TwoPi(), "#varphi"}});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_negInitialQA", ";#varphi;n tracks", kTH1D, {{3 * 810, -TMath::TwoPi(), 2 * TMath::TwoPi(), "#varphi"}});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_negModifiedQA", ";#varphi;n tracks", kTH1D, {{3 * 810, -TMath::TwoPi(), 2 * TMath::TwoPi(), "#varphi"}});

      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPt_pos", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPt_neg", ";p_{T};weighted occupancy", kTH2D, {axisLogPt, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hEta_pos", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hEta_neg", ";#eta;weighted occupancy", kTH2D, {axisEta, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPhi_pos", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});
      histos.add("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPhi_neg", ";#varphi;n tracks", kTH3D, {axisPhi, confAxisOccupForKine, confAxisPtBinsForPhiStudy});

      // QA nTPCcls
      AxisSpec axisNTPCclsPlusMinusQA{521, -260, 260, "n TPC clusters"};
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/QA_tpcNClsFindable_pos", "", kTH1D, {axisNTPCclsPlusMinusQA});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/QA_tpcNClsFound_pos", "", kTH1D, {axisNTPCclsPlusMinusQA});
      histos.add("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/QA_tpcNClsCrossedRows_pos", "", kTH1D, {axisNTPCclsPlusMinusQA});

      AxisSpec axisLogPtFor2D{50, 0.05, 10, "p_{T}"};
      AxisSpec axisLogPtTpcFor2D{50, 0.05, 10, "p_{T} TPC inner"};
      histos.add("track_distr_nITStrThisEv_10_200/hPt_vs_tpcInnerPt_vs_occup", ";p_{T};p_{T} TPC inner;weighted occupancy", kTH3D, {axisLogPtFor2D, axisLogPtTpcFor2D, confAxisOccupForKine});
      histos.add("track_distr_nITStrThisEv_above_2000/hPt_vs_tpcInnerPt_vs_occup", ";p_{T};p_{T} TPC inner;weighted occupancy", kTH3D, {axisLogPtFor2D, axisLogPtTpcFor2D, confAxisOccupForKine});

      histos.add("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy", ";bc in TF;weighted occupancy", kTH2F, {axisBCinTF, confAxisOccupForKine});
      histos.add("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy_vertexTOFmatched", ";bc in TF;weighted occupancy", kTH2F, {axisBCinTF, confAxisOccupForKine});
      histos.add("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy_nPV_10_200", ";bc in TF;weighted occupancy", kTH2F, {axisBCinTF, confAxisOccupForKine});
      histos.add("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy_nPV_above2000", ";bc in TF;weighted occupancy", kTH2F, {axisBCinTF, confAxisOccupForKine});
      // end of July 2025: to compare data and MC (pt, eta, phi)

      // 3D: pt vs centr vs occup
      // if (confFlagManyHeavyHistos) {
      //   histos.add("ptGlobal_vs_centr_vs_occup", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisLogPt});
      //   histos.add("ptPV_vs_centr_vs_occup", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisLogPt});
      //   histos.add("ptGlobal_vs_centr_vs_occup_NoCollStd", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisLogPt});
      //   histos.add("ptPV_vs_centr_vs_occup_NoCollStd", "", kTH3F, {{20, 0, 4000, "nITStrk cls567"}, axisOccupancyForDeDxStudies, axisLogPt});
      // }
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

    // THnD for Marian:
    if (confFlagFillTHn) {
      AxisSpec axis_THnF_phi{confTHnAxis_nPhiBins, 0, TMath::TwoPi(), ""}; // φ at radius: 360 bins
      histos.add("THnD_histos/phi_R_qOp_IR_occ_centr_eta", ";phi;R;qOp;IR;occ;cent;eta", kTHnF, {axis_THnF_phi, confTHnAxis_R, confTHnAxis_qOp, confTHnAxis_IR, confTHnAxis_occ, confTHnAxis_centr, confTHnAxis_eta});

      histos.add("THnD_histos/QA_under_asin", "", kTH1F, {{200, -4, 4}});
      histos.add("THnD_histos/QA_asin", "", kTH1F, {{200, -8, 8}});
    }
    // nD, time bins to cover the range -confTimeIntervalForSmallBins... +confTimeIntervalForSmallBins (us)
    double timeBinSize = 2 * confTimeIntervalForSmallBins / confNumberOfSmallTimeBins;
    std::vector<double> arrTimeBins;
    for (int i = 0; i < confNumberOfSmallTimeBins + 1; i++)
      arrTimeBins.push_back(-confTimeIntervalForSmallBins + i * timeBinSize);
    const AxisSpec axisTimeBins{arrTimeBins, "#Delta t, #mus"};
    int nBinsX = 20;
    int nBinsY = 40;

    if (confAddTimeDependentHistos) {
      histos.add("occupancyInTimeBins", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
      histos.add("occupancyInTimeBins_vs_FT0thisCol_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
      histos.add("occupancyInTimeBins_nITS567_vs_FT0thisCol_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS567cls tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});

      if (confFlagManyHeavyHistos) {
        histos.add("occupancyInTimeBins_BEFORE_sel", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
        histos.add("occupancyInTimeBins_occupByFT0_BEFORE_sel", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});

        histos.add("occupancyInTimeBins_occupByFT0", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});
        histos.add("occupancyInTimeBins_occupByFT0_kNoCollInTimeRangeNarrow", ";time bin (#mus);n ITS tracks with 5,6,7 cls, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBinsX, 0, kMaxThisEv * 4000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});

        histos.add("occupancyInTimeBins_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS-TPC tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});

        histos.add("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITS567cls tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});
        histos.add("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow_NoCollInRofStrict", ";time bin (#mus);FT0C this collision, this collision;n ITS567cls tracks, this collision;sum FT0 in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 100000}});
      }
      if (confFlagIsTOFIsTRDdtStudy) {
        histos.add("occupancyInTimeBins_nITSTOF_vs_FT0thisCol_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITSTOF tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
        histos.add("occupancyInTimeBins_nITSTRD_vs_FT0thisCol_kNoCollInTimeRangeNarrow", ";time bin (#mus);FT0C this collision, this collision;n ITSTRD tracks, this collision;ITS tracks with 5,6,7 cls in time window", kTHnF, {axisTimeBins, {nBins3D, 0, kMaxThisEv * 100000}, {nBinsY, 0, kMaxThisEv * 4000}, {nBins3DoccupancyAxis, 0, kMaxOccup * 10000}});
      }

      histos.add("qaForHighOccupITStracksInTimeBinPast", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
      histos.add("qaForHighOccupITStracksInTimeBinFuture1", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
      histos.add("qaForHighOccupITStracksInTimeBinFuture2", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
      histos.add("qaForHighOccupITStracksForNeighbourEvents", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});

      // save dt information for several first collisions for QA
      histos.add<TH2>("histOccupInTimeBinsQA", ";dt;this coll id", kTH2F, {axisTimeBins, {nCollisionsForTimeBinQA, -0.5, nCollisionsForTimeBinQA - 0.5}});
    }

    histos.add("thisEventITStracksInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("thisEventITSTPCtracksInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});
    histos.add("thisEventFT0CInTimeBins", ";time bin (#mus);n tracks", kTH1F, {axisTimeBins});

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

      // histos.add("nTracksGlobal_vs_nPV_QA_onlyVzCut_noTFROFborderCuts", "nTracksGlobal_vs_nPV_QA_onlyVzCut_noTFROFborderCuts", kTH2F, {axisNtracks, axisNtracksGlobal});
      // histos.add("nTracksGlobal_vs_nPV_QA_after_TFborderCut", "nTracksGlobal_vs_nPV_QA_after_TFborderCut", kTH2F, {axisNtracks, axisNtracksGlobal});

      histos.add("nTracksGlobal_vs_nPV_occupByFT0C_0_2500", "nTracksGlobal_vs_nPV_occupByFT0C_0_2500", kTH2F, {axisNtracks, axisNtracksGlobal});
      histos.add("nTracksGlobal_vs_nPV_occupByFT0C_0_20000", "nTracksGlobal_vs_nPV_occupByFT0C_0_20000", kTH2F, {axisNtracks, axisNtracksGlobal});

      // 3D histograms with occupancy axis
      histos.add("nTracksGlobal_vs_nPV_vs_occup_pure", "nTracksGlobal_vs_nPV_vs_occup_pure", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeStandard_extraCuts", kTH3F, {axisNtracks, axisNtracksGlobal, axisOccupancyTracks});

      // 3D histograms: nGlobalTracks with cls567 as y-axis, V0A as x-axis:
      histos.add("nTracksGlobal_vs_V0A_vs_occup_pure", "", kTH3F, {axisMultV0A, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_V0A_vs_occup_kNoCollInTimeRangeNarrow", "", kTH3F, {axisMultV0A, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_V0A_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "", kTH3F, {axisMultV0A, axisNtracksGlobal, axisOccupancyTracks});
      // FT0C as x-axis:
      histos.add("nTracksGlobal_vs_FT0C_vs_occup_pure", "", kTH3F, {axisMultFT0C, axisNtracksGlobal, axisOccupancyTracks});
      histos.add("nTracksGlobal_vs_FT0C_vs_occup_kNoCollInTimeRangeStandard_extraCuts", "", kTH3F, {axisMultFT0C, axisNtracksGlobal, axisOccupancyTracks});

      // 3D histograms: now - nITStracks with cls567 as y-axis, V0A as x-axis:
      histos.add("nPV_vs_V0A_vs_occup_pure", "", kTH3F, {axisMultV0A, axisNtracks, axisOccupancyTracks});
      histos.add("nPV_vs_V0A_vs_occup_kNoCollInTimeRangeNarrow", "", kTH3F, {axisMultV0A, axisNtracks, axisOccupancyTracks});
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
    std::vector<int> vTracksITSTOFperCollPtEtaCuts(cols.size(), 0); // counter of tracks per found bc for occupancy studies
    std::vector<int> vTracksITSTRDperCollPtEtaCuts(cols.size(), 0); // counter of tracks per found bc for occupancy studies
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
      int nITSTOFtracksPtEtaCuts = 0;
      int nITSTRDtracksPtEtaCuts = 0;
      int nTOFtracks = 0;
      // int nTRDtracks = 0;
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (const auto& track : tracksGrouped) {
        if (!track.isPVContributor()) {
          continue;
        }
        if (track.itsNCls() < confMinITSclsPerTrack)
          continue;
        nITS567cls++;
        nITSTPCtracks += track.hasITS() && track.hasTPC();
        nTOFtracks += track.hasTOF();

        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;

        nITS567clsPtEtaCuts++;
        nITSTOFtracksPtEtaCuts += track.hasITS() && track.hasTOF();
        nITSTRDtracksPtEtaCuts += track.hasITS() && track.hasTRD();

        if (track.tpcNClsFound() < confCutMinTPCcls)
          continue;
        nITSTPCtracksPtEtaCuts += track.hasITS() && track.hasTPC();

        nGlobalPtEtaCuts += track.isGlobalTrack();
      }

      int32_t foundBC = bc.globalIndex();
      int32_t colIndex = col.globalIndex();

      vFoundBCindex[colIndex] = foundBC;
      vFoundGlobalBC[colIndex] = bc.globalBC();

      if (bc.has_foundFT0())
        vAmpFT0CperColl[colIndex] = bc.foundFT0().sumAmpC();

      vIsVertexTOFmatched[colIndex] = nTOFtracks > 0;

      vTracksITS567perColl[colIndex] = nITS567cls;
      vTracksITS567perCollPtEtaCuts[colIndex] = nITS567clsPtEtaCuts;
      vTracksGlobalPerCollPtEtaCuts[colIndex] = nGlobalPtEtaCuts;

      vTracksITSTPCperColl[colIndex] = nITSTPCtracks;
      vTracksITSTPCperCollPtEtaCuts[colIndex] = nITSTPCtracksPtEtaCuts;
      vTracksITSTOFperCollPtEtaCuts[colIndex] = nITSTOFtracksPtEtaCuts;
      vTracksITSTRDperCollPtEtaCuts[colIndex] = nITSTRDtracksPtEtaCuts;

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
      if (!col.selection_bit(kNoTimeFrameBorder))
        continue;

      // skip if collision is close to ROF border
      if (!col.selection_bit(kNoITSROFrameBorder))
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
      std::map<int64_t, int32_t> mUniqueBC;

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
      int binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(0) + 0.5); // default was: -39.5 us
      int binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(1) - 0.5); // -10.5
      float integralPast = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);
      binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(2) + 0.5); // 20.5
      binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(3) - 0.5); // 49.5
      float integralFuture1 = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);
      binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(3) + 0.5); // 50.5
      binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(4) - 0.5); // 79.5
      float integralFuture2 = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);
      binMin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(1) + 0.5); // -9.5
      binMax = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->FindBin(confTimeSlicesForPastFutureStudies->at(2) - 0.5); // 19.5
      float integralNeighbourEvents = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Integral(binMin, binMax);

      // recent past
      if (integralFullDeltaTime < 200) // ~empty detector
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 1;
      // recent past
      if (integralPast > /*3000*/ 2500 && (integralFullDeltaTime - integralPast) < 180) // low occupancy outside the dt region of interest
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 2;
      // close future
      if (integralFuture1 > /*3000*/ 2500 && (integralFullDeltaTime - integralFuture1) < 180)
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 3;
      // distant future
      if (integralFuture2 > /*3000*/ 2500 && (integralFullDeltaTime - integralFuture2) < 180)
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 4;
      // neighbour events
      if (integralNeighbourEvents > /*3000*/ 2500 && (integralFullDeltaTime - integralNeighbourEvents) < 180)
        vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex] = 5;

      // loop over time axis in nD histograms:
      for (int iT = 0; iT < histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetNbinsX(); iT++) {
        int nITStrInTimeBin = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetBinContent(iT + 1);
        if (nITStrInTimeBin == 0) // no collisions in this dt bin
          continue;
        // int nITSTPCtInTimeBin = histos.get<TH1>(HIST("thisEventITSTPCtracksInTimeBins"))->GetBinContent(iT + 1);

        float dt = histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->GetBinCenter(iT + 1);

        int nFT0CInTimeBin = histos.get<TH1>(HIST("thisEventFT0CInTimeBins"))->GetBinContent(iT + 1);

        if (confAddTimeDependentHistos && confFlagManyHeavyHistos) {
          histos.fill(HIST("occupancyInTimeBins_BEFORE_sel"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nITStrInTimeBin);
          histos.fill(HIST("occupancyInTimeBins_occupByFT0_BEFORE_sel"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);
        }
        bool flagFillOccupVsDt = true;
        if (confFlagUseNoCollInRofStrict && !col.selection_bit(kNoCollInRofStrict))
          flagFillOccupVsDt = false;
        if (confFlagUseNoHighMultCollInPrevRof && !col.selection_bit(kNoHighMultCollInPrevRof))
          flagFillOccupVsDt = false;

        if (confAddTimeDependentHistos) {
          if (sel && std::fabs(col.posZ()) < 10 && flagFillOccupVsDt) {
            histos.fill(HIST("occupancyInTimeBins"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nITStrInTimeBin);
            if (confFlagManyHeavyHistos)
              histos.fill(HIST("occupancyInTimeBins_occupByFT0"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);

            if (col.selection_bit(kNoCollInTimeRangeNarrow)) {
              histos.fill(HIST("occupancyInTimeBins_vs_FT0thisCol_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nITStrInTimeBin);
              histos.fill(HIST("occupancyInTimeBins_nITS567_vs_FT0thisCol_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], vTracksITS567perCollPtEtaCuts[colIndex], nITStrInTimeBin);
              if (confFlagIsTOFIsTRDdtStudy) {
                histos.fill(HIST("occupancyInTimeBins_nITSTOF_vs_FT0thisCol_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], vTracksITSTOFperCollPtEtaCuts[colIndex], nITStrInTimeBin);
                histos.fill(HIST("occupancyInTimeBins_nITSTRD_vs_FT0thisCol_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], vTracksITSTRDperCollPtEtaCuts[colIndex], nITStrInTimeBin);
              }
              if (confFlagManyHeavyHistos) {
                histos.fill(HIST("occupancyInTimeBins_occupByFT0_kNoCollInTimeRangeNarrow"), dt, vTracksITS567perCollPtEtaCuts[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);
                histos.fill(HIST("occupancyInTimeBins_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], confFlagUseGlobalTracks ? vTracksGlobalPerCollPtEtaCuts[colIndex] : vTracksITSTPCperCollPtEtaCuts[colIndex], nFT0CInTimeBin);

                histos.fill(HIST("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow"), dt, vAmpFT0CperColl[colIndex], vTracksITS567perCollPtEtaCuts[colIndex], nFT0CInTimeBin);
                if (col.selection_bit(kNoCollInRofStrict))
                  histos.fill(HIST("occupancyInTimeBins_nITS567_vs_FT0thisCol_occupByFT0_kNoCollInTimeRangeNarrow_NoCollInRofStrict"), dt, vAmpFT0CperColl[colIndex], vTracksITS567perCollPtEtaCuts[colIndex], nFT0CInTimeBin);
              }
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
        } // end of confAddTimeDependentHistos
      }

      // reset delta time hist for this event
      histos.get<TH1>(HIST("thisEventITStracksInTimeBins"))->Reset();
      // histos.get<TH1>(HIST("thisEventITSTPCtracksInTimeBins"))->Reset();
      histos.get<TH1>(HIST("thisEventFT0CInTimeBins"))->Reset();
      counterQAtimeOccupHistos++;
    } // end of occupancy calculation

    // ### occupancy event selection QA
    for (const auto& col : cols) {
      if (!col.sel8())
        continue;

      if (confApplyGoodITSstavesFlaginEvSel && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll))
        continue;

      // if (!col.selection_bit(kIsTriggerTVX))
      // continue;

      // cut on vZ for a given collision
      if (col.posZ() < confCutVertZMinThisEvent || col.posZ() > confCutVertZMaxThisEvent)
        continue;

      int32_t colIndex = col.globalIndex();
      int64_t bcInTF = (vFoundGlobalBC[colIndex] - bcSOR) % nBCsPerTF;
      histos.fill(HIST("hNcolVsBcInTF/hNcolVsBcInTF"), bcInTF);
      if (col.selection_bit(kIsVertexTOFmatched))
        histos.fill(HIST("hNcolVsBcInTF/hNcolVsBcInTF_vertexTOFmatched"), bcInTF);

      // cut on the max bcId in the time frame (to avoid the artificial fade-out tail in the MC productions)
      if (!vIsMarkedCollForAnalysis[colIndex])
        continue;

      histos.fill(HIST("hNcolVsBcInTF/hNcolVsBcInTFafterMaxBcCut"), bcInTF);

      auto multV0A = col.multFV0A();
      // auto multT0A = col.multFT0A();
      auto multT0C = col.multFT0C();
      int nPV = 0; // col.multNTracksPV();
      int nGlobalTracks = 0;

      int occupancy = col.trackOccupancyInTimeRange();
      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());

      const auto& bc = col.foundBC_as<BCsRun3>();
      int64_t ts = bc.timestamp();
      double IR = mRateFetcher.fetch(ccdb.service, ts, runNumber, "ZNC hadronic") * 1.e-3; // kHz

      // pre-calc nPV
      for (const auto& track : tracksGrouped) {
        if (!track.isPVContributor())
          continue;
        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;
        if (track.itsNCls() < confMinITSclsPerTrack)
          continue;
        nPV++;
      }
      if (occupancy >= 0 && confAddBasicQAhistos) {
        histos.fill(HIST("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy"), bcInTF, occupancy);
        if (col.selection_bit(kIsVertexTOFmatched))
          histos.fill(HIST("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy_vertexTOFmatched"), bcInTF, occupancy);
        if (nPV >= 10 && nPV < 200)
          histos.fill(HIST("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy_nPV_10_200"), bcInTF, occupancy);
        else if (nPV >= 2000)
          histos.fill(HIST("hNcolVsBcInTF/hNcolVsBcInTF_vs_occupancy_nPV_above2000"), bcInTF, occupancy);
      }

      // main loop for dE/dx
      for (const auto& track : tracksGrouped) {
        histos.fill(HIST("nTrackCounter_after_cuts_QA"), 0);
        if (!track.isPVContributor())
          continue;
        histos.fill(HIST("nTrackCounter_after_cuts_QA"), 1);
        if (track.pt() < confCutPtMinThisEvent || track.pt() > confCutPtMaxThisEvent)
          continue;
        histos.fill(HIST("nTrackCounter_after_cuts_QA"), 2);
        if (track.eta() < confCutEtaMinTracksThisEvent || track.eta() > confCutEtaMaxTracksThisEvent)
          continue;
        histos.fill(HIST("nTrackCounter_after_cuts_QA"), 3);
        if (track.itsNCls() < confMinITSclsPerTrack)
          continue;
        histos.fill(HIST("nTrackCounter_after_cuts_QA"), 4);
        // nPV++;

        // July 2025: more for data vs MC:
        if (track.hasTPC() && occupancy >= 0 && confAddBasicQAhistos) {
          float pt = track.pt();
          // pt 0.2-0.5
          if (pt > 0.2 && pt < 0.5) {
            if (nPV >= 10 && nPV < 400) {
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcNClsFindable_vs_occup_pt_02_05"), track.tpcNClsFindable(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcNClsFound_vs_occup_pt_02_05"), track.tpcNClsFound(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcChi2NCl_vs_occup_pt_02_05"), track.tpcChi2NCl(), occupancy);
            } else if (nPV >= 2000) {
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcNClsFindable_vs_occup_pt_02_05"), track.tpcNClsFindable(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcNClsFound_vs_occup_pt_02_05"), track.tpcNClsFound(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcChi2NCl_vs_occup_pt_02_05"), track.tpcChi2NCl(), occupancy);
            }
          }
          // pt 0.5-1.0
          else if (pt > 0.5 && pt < 1.0) {
            if (nPV >= 10 && nPV < 400) {
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcNClsFindable_vs_occup_pt_05_10"), track.tpcNClsFindable(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcNClsFound_vs_occup_pt_05_10"), track.tpcNClsFound(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcChi2NCl_vs_occup_pt_05_10"), track.tpcChi2NCl(), occupancy);
            } else if (nPV >= 2000) {
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcNClsFindable_vs_occup_pt_05_10"), track.tpcNClsFindable(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcNClsFound_vs_occup_pt_05_10"), track.tpcNClsFound(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcChi2NCl_vs_occup_pt_05_10"), track.tpcChi2NCl(), occupancy);
            }
          }
          // pt > 1.0
          else if (pt > 1.0) {
            if (nPV >= 10 && nPV < 400) {
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcNClsFindable_vs_occup_pt_above1_0"), track.tpcNClsFindable(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcNClsFound_vs_occup_pt_above1_0"), track.tpcNClsFound(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_10_200/tpcChi2NCl_vs_occup_pt_above1_0"), track.tpcChi2NCl(), occupancy);
            } else if (nPV >= 2000) {
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcNClsFindable_vs_occup_pt_above1_0"), track.tpcNClsFindable(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcNClsFound_vs_occup_pt_above1_0"), track.tpcNClsFound(), occupancy);
              histos.fill(HIST("QA_noTPCcuts/nPV_above2000/tpcChi2NCl_vs_occup_pt_above1_0"), track.tpcChi2NCl(), occupancy);
            }
          }
        }

        if (track.isGlobalTrack() && track.tpcNClsFound() >= confCutMinTPCcls) {
          nGlobalTracks++;
          histos.fill(HIST("nTrackCounter_after_cuts_QA"), 5);

          if (track.passedTPCRefit() && confAddBasicQAhistos) {
            histos.fill(HIST("nTrackCounter_after_cuts_QA"), 6);

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
              histos.fill(HIST("nTrackCounter_after_cuts_QA"), 7);

              histos.fill(HIST("dEdx_vs_Momentum_vs_occup"), signedP, track.tpcSignal(), occupancy);

              if (confFlagManyHeavyHistos) {
                if (track.eta() > 0.2 && track.eta() < 0.4)
                  histos.fill(HIST("dEdx_vs_Momentum_vs_occup_eta_02_04"), signedP, track.tpcSignal(), occupancy);
                if (track.eta() > -0.4 && track.eta() < -0.2)
                  histos.fill(HIST("dEdx_vs_Momentum_vs_occup_eta_04_02"), signedP, track.tpcSignal(), occupancy);
              }
              // more QA for TPC cls counting
              histos.fill(HIST("tpcNClsFindable"), track.tpcNClsFindable());
              histos.fill(HIST("tpcNClsFindableMinusFound"), track.tpcNClsFindableMinusFound());
              histos.fill(HIST("tpcNClsFindableMinusCrossedRows"), track.tpcNClsFindableMinusCrossedRows());
              histos.fill(HIST("tpcNClsShared"), track.tpcNClsShared());
              histos.fill(HIST("tpcNClsFindableMinusPID"), track.tpcNClsFindableMinusPID());
              int tpcNClUsedForPID = track.tpcNClsFindable() - track.tpcNClsFindableMinusPID();
              histos.fill(HIST("tpcNClUsedForPID"), tpcNClUsedForPID);

              histos.fill(HIST("tpcNClsFound"), track.tpcNClsFound());
              histos.fill(HIST("tpcNClsFoundAsDiffByHand"), track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());

              histos.fill(HIST("tpcNClsUsedForPID_vs_Findable"), track.tpcNClsFindable(), tpcNClUsedForPID);
              histos.fill(HIST("tpcNClsShared_vs_Findable"), track.tpcNClsFindable(), track.tpcNClsShared());
              histos.fill(HIST("tpcNClsUsedForPID_vs_Shared"), track.tpcNClsShared(), tpcNClUsedForPID);
              histos.fill(HIST("tpcNClsFound_vs_Findable"), track.tpcNClsFindable(), track.tpcNClsFound());
              histos.fill(HIST("tpcNClsUsedForPID_vs_Found"), track.tpcNClsFound(), tpcNClUsedForPID);

              int tpcNClsCorrectedFindableMinusPID = track.tpcNClsFindableMinusPID();
              // correct for a buggy behaviour due to int8 and uint8 difference:
              if (tpcNClsCorrectedFindableMinusPID < -70)
                tpcNClsCorrectedFindableMinusPID += 256;
              histos.fill(HIST("tpcNClsFindableMinusPID_CORRECTED"), tpcNClsCorrectedFindableMinusPID);
              histos.fill(HIST("tpcNClsUsedForPID_vs_Findable_CORRECTED"), track.tpcNClsFindable(), track.tpcNClsFindable() - tpcNClsCorrectedFindableMinusPID);

              histos.fill(HIST("tpcNClsFoundMinusPID_BY_HAND"), (track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) - (track.tpcNClsFindable() - tpcNClsCorrectedFindableMinusPID));

              // check ratio tpcNClsFindableMinusPID / tpcNClsFindable
              // https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h#L242
              // https://github.com/AliceO2Group/AliceO2/blob/dev/Detectors/AOD/src/AODProducerWorkflowSpec.cxx#L2553C21-L2553C44
              float fractionTPCcls = (1.0 * tpcNClsCorrectedFindableMinusPID) / track.tpcNClsFindable();
              histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup"), occupancy, fractionTPCcls);
              if (fractionTPCcls >= 0 && fractionTPCcls < 0.8)
                histos.fill(HIST("nTrackCounter_after_cuts_QA"), 8);
              if (fractionTPCcls < 0)
                histos.fill(HIST("dEdx_vs_Momentum_HighFractionNclsNonPID"), signedP, track.tpcSignal());
              if (fractionTPCcls > 0.8)
                histos.fill(HIST("dEdx_vs_Momentum_NegativeFractionNclsPID"), signedP, track.tpcSignal());

              if (multV0A < 6800)
                histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup_peripheralByV0A"), occupancy, fractionTPCcls);
              else if (multV0A > 82850)
                histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup_centralByV0A"), occupancy, fractionTPCcls);

              if (std::fabs(track.eta()) < 0.2)
                histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup_eta02"), occupancy, fractionTPCcls);

              // vs charge
              if (signedP > 0) {
                histos.fill(HIST("nTrackCounter_after_cuts_QA"), 9);
                histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup_pos"), occupancy, fractionTPCcls);
              } else {
                histos.fill(HIST("nTrackCounter_after_cuts_QA"), 10);
                histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup_neg"), occupancy, fractionTPCcls);
              }
              // vs pt
              if (track.pt() > 0.2 && track.pt() < 0.8)
                histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup_lowPt"), occupancy, fractionTPCcls);
              if (track.pt() > 0.8 && track.pt() < 10)
                histos.fill(HIST("fraction_tpcNClsFindableMinusPID_vs_occup_highPt"), occupancy, fractionTPCcls);

              // dE/dx in narrow mom bin vs centrality and occupancy
              if (std::fabs(signedP) > 0.38 && std::fabs(signedP) < 0.4)
                histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win"), nPV, occupancy, track.tpcSignal());
              // vs charge
              if (signedP > 0.38 && signedP < 0.4) {
                histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win_pos"), nPV, occupancy, track.tpcSignal());
                if (fractionTPCcls >= 0 && fractionTPCcls < 0.8)
                  histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win_pos_FractionPIDclsInRange"), nPV, occupancy, track.tpcSignal());
              } else if (signedP > -0.4 && signedP < -0.38) {
                histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win_neg"), nPV, occupancy, track.tpcSignal());
                if (fractionTPCcls >= 0 && fractionTPCcls < 0.8)
                  histos.fill(HIST("dEdx_vs_centr_vs_occup_narrow_p_win_neg_FractionPIDclsInRange"), nPV, occupancy, track.tpcSignal());
              }

              // nTPCcls vs nITStr vs occup
              histos.fill(HIST("tpcNClsFound_vs_centr_vs_occup"), nPV, occupancy, track.tpcNClsFound());
              histos.fill(HIST("tpcNClsFindable_vs_centr_vs_occup"), nPV, occupancy, track.tpcNClsFindable());
              histos.fill(HIST("tpcNClsShared_vs_centr_vs_occup"), nPV, occupancy, track.tpcNClsShared());

              // nTPCsharedCls for A and C separately
              if (track.tgl() > 0.) // A side
                histos.fill(HIST("tpcNClsShared_vs_centr_vs_occup_Aside"), nPV, occupancy, track.tpcNClsShared());
              else // C side
                histos.fill(HIST("tpcNClsShared_vs_centr_vs_occup_Cside"), nPV, occupancy, track.tpcNClsShared());

              // nTPCsharedCls for pos and neg
              if (signedP > 0)
                histos.fill(HIST("tpcNClsShared_vs_centr_vs_occup_pos"), nPV, occupancy, track.tpcNClsShared());
              else
                histos.fill(HIST("tpcNClsShared_vs_centr_vs_occup_neg"), nPV, occupancy, track.tpcNClsShared());
            }
          }
        }
      } // end of track loop

      // if (confAddTracksVsFwdHistos)
      // histos.fill(HIST("nTracksGlobal_vs_nPV_QA_onlyVzCut_noTFROFborderCuts"), nPV, nGlobalTracks);

      // skip if collision is close to TF border
      // if (confFlagApplyTFborderCut && !col.selection_bit(kNoTimeFrameBorder))
      // continue;

      // if (confAddTracksVsFwdHistos)
      // histos.fill(HIST("nTracksGlobal_vs_nPV_QA_after_TFborderCut"), nPV, nGlobalTracks);

      // skip if collision is close to ROF border
      // if (confFlagApplyROFborderCut && !col.selection_bit(kNoITSROFrameBorder))
      // continue;

      histos.fill(HIST("hOccupancy"), occupancy);
      if (occupancy >= 0 && confAddBasicQAhistos) {
        int orbitId = bcInTF / o2::constants::lhc::LHCMaxBunches;
        histos.fill(HIST("hOccupancyVsOrbit"), orbitId, occupancy);
      }

      // another track loop to fill track-level histograms
      if (confAddBasicQAhistos) {
        int flagWhichDeltaTimeWin = vFlagsForEtaQAvsOccupancyInDeltaTimeWins[colIndex];
        bool flagNoCollNearby = col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow);

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
          if (track.itsNCls() < confMinITSclsPerTrack)
            continue;
          // if (!(track.isGlobalTrack() && track.tpcNClsFound() >= confCutMinTPCcls))
          // continue;

          bool isGoodGlobal = (track.isGlobalTrack() && track.tpcNClsFound() >= confCutMinTPCcls);
          bool hasTPCspecCuts = (track.hasTPC() && track.tpcNClsFound() >= confCutMinTPCcls && track.tpcNClsCrossedRows() > 80 && track.tpcChi2NCl() < 4);

          // ### kine distr vs centr vs occup
          float sign = track.sign();
          float pt = track.pt();
          float eta = track.eta();
          float phi = track.phi();
          float phiInitial = phi;

          if (confUsePhiAtTPCinnerR) {
            phi -= asin(confRadiusForPhiCorrection /*inner TPC radius*/ / 2 * 0.3 * sign * 0.5 / pt);
            if (phi < 0)
              phi += TMath::TwoPi();
            else if (phi > TMath::TwoPi())
              phi -= TMath::TwoPi();
          }

          bool etaInRange = true;
          if (confUseAorCsideForPhiStudy == 0 && eta < 0.1) // check if we are in A side
            etaInRange = false;
          if (confUseAorCsideForPhiStudy == 1 && eta > -0.1) // check if we are in C side
            etaInRange = false;

          if (occupancy >= 0 && fabs(eta) < 0.8 && pt > 0.15 && etaInRange) {
            if (nPV >= 10 && nPV < 200) {
              if (isGoodGlobal) {
                if (flagWhichDeltaTimeWin == 1 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_lowOccupInTPC"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_lowOccupInTPC"), pt);

                  if (sign > 0)
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC_pos_vs_pt"), phi, pt);
                  else
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_lowOccupInTPC_neg_vs_pt"), phi, pt);
                }
                if (flagWhichDeltaTimeWin == 2 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInRecentPast"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInRecentPast"), pt);

                  if (sign > 0)
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast_pos_vs_pt"), phi, pt);
                  else
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInRecentPast_neg_vs_pt"), phi, pt);
                }
                if (flagWhichDeltaTimeWin == 3 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInCloseFuture"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInCloseFuture"), pt);

                  if (sign > 0)
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture_pos_vs_pt"), phi, pt);
                  else
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInCloseFuture_neg_vs_pt"), phi, pt);
                }
                if (flagWhichDeltaTimeWin == 4 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInDistantFuture"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInDistantFuture"), pt);

                  if (sign > 0)
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture_pos_vs_pt"), phi, pt);
                  else
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInDistantFuture_neg_vs_pt"), phi, pt);
                }
                if (flagWhichDeltaTimeWin == 5) {
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hEta_highOccupInNeighbourEvents"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_highOccupInNeighbourEvents"), pt);

                  if (sign > 0)
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents_pos_vs_pt"), phi, pt);
                  else
                    histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPhi_highOccupInNeighbourEvents_neg_vs_pt"), phi, pt);
                }
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/hPt_vs_tpcInnerPt_vs_occup"), pt, track.tpcInnerParam(), occupancy);
              } // end of TPC good global

              // July 2025: for data vs MC kine distr comparison

              int tpcNClsFindable = track.tpcNClsFindable();
              int tpcNClsFound = track.tpcNClsFound();
              int tpcNClsCrossedRows = track.tpcNClsCrossedRows();

              if (sign > 0) // positive tracks
              {
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPt_pos"), pt, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hEta_pos"), eta, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPhi_pos"), phi, occupancy, pt);
                if (hasTPCspecCuts) {
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPt_pos"), pt, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hEta_pos"), eta, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_pos"), phi, occupancy, pt);

                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFindable_pos"), phi, occupancy, pt, tpcNClsFindable);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFound_pos"), phi, occupancy, pt, tpcNClsFound);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsCrossedRows_pos"), phi, occupancy, pt, tpcNClsCrossedRows);

                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/QA_tpcNClsFindable_pos"), tpcNClsFindable);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/QA_tpcNClsFound_pos"), tpcNClsFound);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/QA_tpcNClsCrossedRows_pos"), tpcNClsCrossedRows);
                }
              } else // negative tracks
              {
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPt_neg"), pt, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hEta_neg"), eta, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/PV_hPhi_neg"), phi, occupancy, pt);
                if (hasTPCspecCuts) {
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPt_neg"), pt, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hEta_neg"), eta, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_neg"), phi, occupancy, pt);

                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFindable_neg"), phi, occupancy, pt, tpcNClsFindable);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsFound_neg"), phi, occupancy, pt, tpcNClsFound);
                  histos.fill(HIST("track_distr_nITStrThisEv_10_200/kine_vs_weighted_occup/hPhi_tpcNClsCrossedRows_neg"), phi, occupancy, pt, tpcNClsCrossedRows);
                }
              }
              // end of July 2025: for data vs MC kine distr comparison

            } else if (nPV >= 2000) {
              if (isGoodGlobal) {
                if (flagWhichDeltaTimeWin == 1 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_lowOccupInTPC"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_lowOccupInTPC"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_lowOccupInTPC"), pt);
                }
                if (flagWhichDeltaTimeWin == 2 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInRecentPast"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInRecentPast"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInRecentPast"), pt);
                }
                if (flagWhichDeltaTimeWin == 3 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInCloseFuture"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInCloseFuture"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInCloseFuture"), pt);
                }
                if (flagWhichDeltaTimeWin == 4 && flagNoCollNearby) {
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInDistantFuture"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInDistantFuture"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInDistantFuture"), pt);
                }
                if (flagWhichDeltaTimeWin == 5) {
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hEta_highOccupInNeighbourEvents"), eta);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPhi_highOccupInNeighbourEvents"), phi);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_highOccupInNeighbourEvents"), pt);
                }
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/hPt_vs_tpcInnerPt_vs_occup"), pt, track.tpcInnerParam(), occupancy);
              } // end of TPC good global

              // July 2025: for data vs MC kine distr comparison
              if (sign > 0) // positive tracks
              {
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPt_pos"), pt, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hEta_pos"), eta, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPhi_pos"), phi, occupancy, pt);
                if (hasTPCspecCuts) {
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPt_pos"), pt, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hEta_pos"), eta, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_pos"), phi, occupancy, pt);
                  if (pt > 0.7 && pt < 1.0) {
                    histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_posInitialQA"), phiInitial);
                    histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_posModifiedQA"), phi);
                  }
                }
              } else // negative tracks
              {
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPt_neg"), pt, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hEta_neg"), eta, occupancy);
                histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/PV_hPhi_neg"), phi, occupancy, pt);
                if (hasTPCspecCuts) {
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPt_neg"), pt, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hEta_neg"), eta, occupancy);
                  histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_neg"), phi, occupancy, pt);
                  if (pt > 0.7 && pt < 1.0) {
                    histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_negInitialQA"), phiInitial);
                    histos.fill(HIST("track_distr_nITStrThisEv_above_2000/kine_vs_weighted_occup/hPhi_negModifiedQA"), phi);
                  }
                }
              }
              // end of July 2025: for data vs MC kine distr comparison
            } // end of if (nPV >= 2000)
          } // end of if (occupancy >= 0) && kine cuts
        } // end of spec track loop to fill track histograms
      } // end of if (confAddBasicQAhistos)

      // special loop to fill THn histograms
      if (confFlagFillTHn && occupancy >= 0) {
        for (const auto& track : tracksGrouped) {
          if (!track.isPVContributor())
            continue;
          if (track.itsNCls() < confMinITSclsPerTrack)
            continue;

          float pt = track.pt();
          float eta = track.eta();

          if (fabs(eta) > 0.8)
            continue;
          if (pt < 0.15)
            continue;

          bool hasTPCspecCuts = (track.hasTPC() && track.tpcNClsFound() >= confCutMinTPCcls && track.tpcNClsCrossedRows() > 80 && track.tpcChi2NCl() < 4);
          if (!hasTPCspecCuts)
            continue;

          float sign = track.sign();
          // if (sign < 0)
          // continue;

          float qpt = track.signed1Pt();

          // fill THnF:
          for (int iRadius = 0; iRadius < 8; iRadius++) {
            float R = (iRadius == 0 ? 0 : 0.8 + iRadius * 0.2); // cm
            float phiAtR = track.phi();
            if (iRadius > 0) {
              histos.fill(HIST("THnD_histos/QA_under_asin"), R / 2 * 0.3 * sign * 0.5 / pt);
              histos.fill(HIST("THnD_histos/QA_asin"), asin(R / 2 * 0.3 * sign * 0.5 / pt));

              phiAtR -= asin(R / 2 * 0.3 * sign * 0.5 / pt);
              if (phiAtR < 0)
                phiAtR += TMath::TwoPi();
              else if (phiAtR > TMath::TwoPi())
                phiAtR -= TMath::TwoPi();
            }
            histos.fill(HIST("THnD_histos/phi_R_qOp_IR_occ_centr_eta"), phiAtR, iRadius, qpt, IR, occupancy, nPV, eta);
          }
        }
      } // end of confFlagFillTHn

      // occupancy vs centrality
      if (confFlagCentralityIsAvailable) {
        auto t0cCentr = col.centFT0C();
        if (occupancy >= 0) {
          histos.fill(HIST("hCentrVsOccupancy"), t0cCentr, occupancy);
          if (col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
            histos.fill(HIST("hCentrVsOccupancyNoCollStd"), t0cCentr, occupancy);
        }
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
        if (occupancy >= 0) {
          histos.fill(HIST("nTracksGlobal_vs_nPV_vs_occup_kNoCollInTimeRangeNarrow"), nPV, nGlobalTracks, occupancy);
          histos.fill(HIST("nTracksGlobal_vs_V0A_vs_occup_kNoCollInTimeRangeNarrow"), multV0A, nGlobalTracks, occupancy);
          histos.fill(HIST("nPV_vs_V0A_vs_occup_kNoCollInTimeRangeNarrow"), multV0A, nPV, occupancy);
        }
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
