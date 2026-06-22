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
/// \file kstarInOO.cxx
/// \brief the pT spectra of k*0(892) resonance analysis in OO collisions
/// \author Jimun Lee <jimun.lee@cern.ch>

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

// jet
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/TrackJetQa.h"

#include <CCDB/BasicCCDBManager.h>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TTree.h>
#include <TVector2.h>

#include <RtypesCore.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <stdlib.h>
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct kstarInOO {
  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //==================================
  //|
  //|         Selection
  //|
  //==================================
  // Event Selection
  Configurable<std::string> cfgEventSelections{"cfgEventSelections", "sel8", "Set event selection"};
  Configurable<float> cfgEventVtxCut{"cfgEventVtxCut", 10.0, "V_z cut selection"};
  Configurable<double> cfgEventMaxEta{"cfgEventMaxEta", 1.0, "set INEL event eta cut"};
  ConfigurableAxis cfgCentAxis{"cfgCentAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0}, "Binning of the centrality axis"};
  Configurable<bool> cfgOccupancySel{"cfgOccupancySel", false, "Occupancy selection"};
  Configurable<float> cfgOccupancyMax{"cfgOccupancyMax", 999999., "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<float> cfgOccupancyMin{"cfgOccupancyMin", -100., "minimum occupancy of tracks in neighbouring collisions in a given time range"};

  // Track Selection
  // General
  Configurable<double> cfgTrackMinPt{"cfgTrackMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgTrackMaxEta{"cfgTrackMaxEta", 0.8, "set track max Eta"};
  Configurable<double> cfgTrackMaxDCArToPVcut{"cfgTrackMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgTrackMaxDCAzToPVcut{"cfgTrackMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgTrackGlobalSel{"cfgTrackGlobalSel", true, "Global track selection"};
  Configurable<bool> cfgTrackPrimaryTrack{"cfgTrackPrimaryTrack", true, "Primary track selection"};          // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgTrackConnectedToPV{"cfgTrackConnectedToPV", true, "PV contributor track selection"}; // PV Contriuibutor
  Configurable<bool> cfgTrackGlobalWoDCATrack{"cfgTrackGlobalWoDCATrack", true, "Global track selection without DCA"};
  // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  // TPC
  Configurable<double> cfgTrackFindableTPCClusters{"cfgTrackFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgTrackTPCCrossedRows{"cfgTrackTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgTrackRowsOverFindable{"cfgTrackRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgTrackTPCChi2{"cfgTrackTPCChi2", 4.0, "nTPC Chi2 per Cluster"};

  // ITS
  Configurable<double> cfgTrackITSChi2{"cfgTrackITSChi2", 36.0, "nITS Chi2 per Cluster"};

  // PID
  Configurable<bool> cfgTrackTPCPID{"cfgTrackTPCPID", true, "Enables TPC PID"};
  Configurable<bool> cfgTrackTOFPID{"cfgTrackTOFPID", true, "Enables TOF PID"};
  Configurable<bool> cfgTrackSquarePIDCut{"cfgTrackSqurePIDCut", true, "Enables PID cut shape square switch"};
  Configurable<bool> cfgTrackCirclePIDCut{"cfgTrackCirclePIDCut", true, "Enables PID cut shape circle switch"};
  Configurable<int> cfgTrackCircleValue{"cfgTrackCircleValue", 2, "Enables TOF TPC PID circle cut value"};
  Configurable<bool> cfgTrackTOFHard{"cfgTrackTOFHard", false, "Enables TOF Hard"};

  Configurable<float> cfgTrackTPCPIDnSig{"cfgTrackTPCPIDnSig", 4.0, "nTPC PID sigma"};
  Configurable<float> cfgTrackTOFPIDnSig{"cfgTrackTOFPIDnSig", 4.0, "nTOF PID sigma"};
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};

  // Mixing
  ConfigurableAxis cfgBinsMixMult{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  ConfigurableAxis cfgBinsMixVtx{"cfgBinsMixVtx", {VARIABLE_WIDTH, -10.0f, -5.f, 0.f, 5.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> cfgMixNMixedEvents{"cfgMixNMixedEvents", 10, "Number of mixed events per event"};
  Configurable<int> cfgVtxMixCut{"cfgVtxMixCut", 10, "Vertex Mix Cut"};

  // MCGen
  Configurable<bool> cfgForceGenReco{"cfgForceGenReco", false, "Only consider events which are reconstructed (neglect event-loss)"};
  Configurable<bool> cfgReqMcEffPID{"cfgReqMcEffPID", false, "Request McEfficiency PID"};
  Configurable<bool> cfgReqMcEffTrackQA{"cfgRecMcEffTrackQA", false, "Enable Jet Cut on Trig"};

  // Pair
  Configurable<int> cfgMinvNBins{"cfgMinvNBins", 300, "Number of bins for Minv axis"};
  Configurable<float> cfgMinvMin{"cfgMinvMin", 0.60, "Minimum Minv value"};
  Configurable<float> cfgMinvMax{"cfgMinvMax", 1.20, "Maximum Minv value"};

  // Histogram
  ConfigurableAxis binsDCAz{"binsDCAz", {40, -0.2, 0.2}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {40, -0.2, 0.2}, ""};
  Configurable<bool> cfgEventCutQA{"cfgEventCutsQA", false, "Enable Event QA Hists"};
  Configurable<bool> cfgTrackCutQA{"cfgTrackCutQA", false, "Enable Track QA Hists"};
  Configurable<bool> cfgJetQAHistos{"cfgJetQAHistos", false, "Enable Jet QA Histos"};
  Configurable<bool> cfgWoJetQA{"cfgWoJetQA", false, "Enable Without Jet QA Histos"};

  Configurable<bool> cfgMCHistos{"cfgMCHistos", false, "Enable MC Hists"};
  Configurable<bool> cfgMixedHistos{"cfgMixedHistos", false, "Enable Mixed Histos"};
  Configurable<bool> cfgJetHistos{"cfgJetHistos", false, "Enable Jet Histos"};
  Configurable<bool> cfgJetMCHistos{"cfgJetMCHistos", false, "Enable Jet MC Histos"};
  Configurable<bool> cfgJetdRHistos{"cfgJetdRHistos", false, "Enable Jet dR QA Histos"};
  Configurable<bool> cfgCutonTrig{"cfgCutonTrig", false, "Enable Jet Cut on Trig"};

  //======================
  //|
  //|        JET
  //|
  //======================
  Configurable<float> cfgJetpT{"cfgJetpT", 8.0, "Set Jet pT minimum"};
  Configurable<float> cfgJetR{"cfgJetR", 0.4, "Set Jet radius parameter"};
  Configurable<bool> cfgSingleJet{"cfgSingleJet", false, "Enforces strict phi-jet correspondance"};
  Configurable<bool> cfgReqJets{"cfgReqJets", false, "False: MB, True: Inside Jets"};
  Configurable<std::string> cfgRealTriggerMasks{"cfgRealTriggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<std::string> cfgTriggerMasksTest1{"cfgTriggerMasksTest1", "", "possible JE Trigger masks Test1"};
  Configurable<std::string> cfgTriggerMasksTest2{"cfgTriggerMasksTest2", "", "possible JE Trigger masks Test2"};
  Configurable<std::string> cfgTriggerMasksTest3{"cfgTriggerMasksTest3", "", "possible JE Trigger masks Test3"};
  Configurable<bool> cfgForceTrueINELgt{"cfgForceTrueINELgt", true, "Check generated event with True INELgt0"};

  Configurable<bool> cfgIsKstar{"cfgIsKstar", false, "Swaps Phi for Kstar analysis"};
  Configurable<bool> cfgBR{"cfgBR", false, "Forces Gen. Charged BR Only"};

  std::vector<int> eventSelectionBits;
  std::vector<int> RealTriggerMaskBits;
  std::vector<int> triggerMaskBitsTest1;
  std::vector<int> triggerMaskBitsTest2;
  std::vector<int> triggerMaskBitsTest3;

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPhi{200, -1, +7, "#phi"};
    const AxisSpec ptAxis = {200, 0, 20.0};
    const AxisSpec pidAxis = {120, -6, 6};
    const AxisSpec minvAxis = {cfgMinvNBins, cfgMinvMin, cfgMinvMax};
    const AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    const AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};
    const AxisSpec dRAxis = {100, 0, 100};

    if (cfgEventCutQA) {
      histos.add("hEvent_Cut", "Number of event after cuts", kTH1D, {{13, -0.5, 12.5}});
      histos.add("hPosZ_BC", "hPosZ_BC", kTH1F, {{300, -15.0, 15.0}});
      histos.add("hPosZ_AC", "hPosZ_AC", kTH1F, {{300, -15.0, 15.0}});
      histos.add("hcentFT0C_BC", "centFT0C_BC", kTH1F, {{110, 0.0, 110.0}});
      histos.add("hcentFT0C_AC", "centFT0C_AC", kTH1F, {{110, 0.0, 110.0}});

      std::shared_ptr<TH1> hCutFlow = histos.get<TH1>(HIST("hEvent_Cut"));
      std::vector<std::string> eventCutLabels = {
        "All Events",
        "sel8",
        Form("|Vz| < %.1f", cfgEventVtxCut.value),
        "kIsGoodZvtxFT0vsPV",
        "kNoSameBunchPileup",
        "kNoTimeFrameBorder",
        "kNoITSROFrameBorder",
        "kNoCollInTimeRangeStandard",
        "kIsGoodITSLayersAll",
        Form("Occupancy < %.0f", cfgOccupancyMax.value),
        "All passed events"};
      for (size_t i = 0; i < eventCutLabels.size(); ++i) {
        hCutFlow->GetXaxis()->SetBinLabel(i + 1, eventCutLabels[i].c_str());
      }
    }
    if (cfgTrackCutQA) {
      histos.add("hDCArToPv_BC", "DCArToPv_BC", kTH1F, {axisDCAxy});
      histos.add("hDCAzToPv_BC", "DCAzToPv_BC", kTH1F, {axisDCAz});
      histos.add("hIsPrim_BC", "hIsPrim_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsGood_BC", "hIsGood_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsPrimCont_BC", "hIsPrimCont_BC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hFindableTPCClusters_BC", "hFindableTPCClusters_BC", kTH1F, {{200, 0, 200}});
      histos.add("hFindableTPCRows_BC", "hFindableTPCRows_BC", kTH1F, {{200, 0, 200}});
      histos.add("hClustersVsRows_BC", "hClustersVsRows_BC", kTH1F, {{200, 0, 2}});
      histos.add("hTPCChi2_BC", "hTPCChi2_BC", kTH1F, {{200, 0, 100}});
      histos.add("hITSChi2_BC", "hITSChi2_BC", kTH1F, {{200, 0, 100}});
      histos.add("QA_nSigma_pion_TPC_BC", "QA_nSigma_pion_TPC_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_pion_TOF_BC", "QA_nSigma_pion_TOF_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_pion_TPC_TOF_BC", "QA_pion_TPC_TOF_BC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TPC_BC", "QA_nSigma_kaon_TPC_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TOF_BC", "QA_nSigma_kaon_TOF_BC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_kaon_TPC_TOF_BC", "QA_kaon_TPC_TOF_BC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_track_pT_BC", "QA_track_pT_BC", kTH1F, {{13, 0.0, 13.0}});

      histos.add("hDCArToPv_AC", "DCArToPv_AC", kTH1F, {axisDCAxy});
      histos.add("hDCAzToPv_AC", "DCAzToPv_AC", kTH1F, {axisDCAz});
      histos.add("hIsPrim_AC", "hIsPrim_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsGood_AC", "hIsGood_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hIsPrimCont_AC", "hIsPrimCont_AC", kTH1F, {{2, -0.5, 1.5}});
      histos.add("hFindableTPCClusters_AC", "hFindableTPCClusters_AC", kTH1F, {{200, 0, 200}});
      histos.add("hFindableTPCRows_AC", "hFindableTPCRows_AC", kTH1F, {{200, 0, 200}});
      histos.add("hClustersVsRows_AC", "hClustersVsRows_AC", kTH1F, {{200, 0, 2}});
      histos.add("hTPCChi2_AC", "hTPCChi2_AC", kTH1F, {{200, 0, 100}});
      histos.add("hITSChi2_AC", "hITSChi2_AC", kTH1F, {{200, 0, 100}});
      histos.add("QA_nSigma_pion_TPC_AC", "QA_nSigma_pion_TPC_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_pion_TOF_AC", "QA_nSigma_pion_TOF_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_pion_TPC_TOF_AC", "QA_pion_TPC_TOF_AC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TPC_AC", "QA_nSigma_kaon_TPC_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_nSigma_kaon_TOF_AC", "QA_nSigma_kaon_TOF_AC", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA_kaon_TPC_TOF_AC", "QA_kaon_TPC_TOF_AC", {HistType::kTH2F, {pidAxis, pidAxis}});
      histos.add("QA_track_pT_AC", "QA_track_pT_AC", kTH1F, {{13, 0.0, 13.0}});
    }
    if (cfgJetQAHistos) {
      histos.add("nTriggerQA", "nTriggerQA", kTH1F, {{8, 0.0, 8.0}});

      histos.add("JetpT", "Jet pT (GeV/c)", kTH1F, {{4000, 0., 200.}});
      histos.add("JetEta", "Jet Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("JetPhi", "Jet Phi", kTH1F, {{80, -1.0, 7.0}});

      histos.add("nGoodJets", "The number of good jets", kTH1F, {{6, -0.5, 5.5}});
      histos.add("nJetsPerEvent", "The number of jet per event", kTH1F, {{6, -0.5, 5.5}});

      histos.add("rawDimpT", "rawDimpT", kTH2F, {{1000, 0.0, 10.0}, {100, -0.5, 0.5}});

      histos.add("jetTrackEta", "Jet Track Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("jetTrackPhi", "Jet Track Phi", kTH1F, {{80, -1.0, 7.0}});
    }
    if (cfgJetdRHistos) {
      histos.add("mcpjet_pt", "MC particle-level Jet pT (GeV/c)", kTH1F, {{2000, 0., 100.}});
      histos.add("mcpjet_eta", "MC particle-level Jet Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("mcpjet_phi", "MC particle-level Jet Phi", kTH1F, {{80, -1.0, 7.0}});

      histos.add("Gen_particle_All_BR", "Generated particle All pT (GeV/c)", kTH1F, {{2000, 0., 100.}});
      histos.add("Gen_particle_BR", "Gen_particle_BR", kTH1F, {minvAxis});
      histos.add("Gen_particle_pT", "Generated particle pT (GeV/c)", kTH1F, {{2000, 0., 100.}});

      histos.add("dR_taggedjet_kaon", "dR between tagged jet and kaon wo pT cut", {HistType::kTH2F, {{60, 0., 1.5}, {100, 0., 20.}}});
      histos.add("dR_taggedjet_pion", "dR between tagged jet and pion wo pT cut", {HistType::kTH2F, {{60, 0., 1.5}, {100, 0., 20.}}});
      histos.add("dR_taggedjet_all", "dR between tagged jet and kpi", {HistType::kTH2F, {{60, 0., 1.5}, {100, 0., 20.}}});
      histos.add("dR_taggedjet_all_6_8", "dR between tagged jet and kpi 6 < jet pT < 8", {HistType::kTH2F, {{60, 0., 1.5}, {100, 0., 20.}}});

      histos.add("missed_kpi_INJets_6_8", "missed kpi In Jets with 6 < jetPt < 8", {HistType::kTH2F, {{120, 0.0, 1.2}, {100, 0., 20.}}});
      histos.add("missed_kpi_INJets_8_10", "missed kpi In Jets with 8 < jetPt < 10", {HistType::kTH2F, {{120, 0.0, 1.2}, {100, 0., 20.}}});
      histos.add("missed_kpi_INJets_10_12", "missed kpi In Jets with 10 < jetPt < 12", {HistType::kTH2F, {{120, 0.0, 1.2}, {100, 0., 20.}}});
      histos.add("missed_kpi_INJets_12_15", "missed kpi In Jets with 12 < jetPt < 15", {HistType::kTH2F, {{120, 0.0, 1.2}, {100, 0., 20.}}});
      histos.add("missed_kpi_INJets_15_25", "missed kpi In Jets with 15 < jetPt < 25", {HistType::kTH2F, {{120, 0.0, 1.2}, {100, 0., 20.}}});
      histos.add("missed_kpi_INJets_25_infinite", "missed kpi In Jets with 25 < jetPt < infinite", {HistType::kTH2F, {{120, 0.0, 1.2}, {100, 0., 20.}}});
      histos.add("missed_kpi_INJets_8_infinite", "missed kpi In Jets with 8 < jetPt < infinite", {HistType::kTH2F, {{120, 0.0, 1.2}, {100, 0., 20.}}});

      histos.add("recoveredJetpT_6_8to8_10", "recovered Jet pT", kTH1F, {{2000, 0., 100.}});

      histos.add("JetMigration", "bin to bin migration", {HistType::kTH2F, {{100, 0.0, 50.0, "True jet pT (GeV/c)"}, {100, 0., 50., "Recovered jet pT (GeV/c)"}}});
    }

    ////////////////////////////////////
    histos.add("nEvents", "nEvents", kTH1F, {{7, 0.0, 7.0}});
    histos.add("hUSS_KPi", "hUSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hUSS_PiK", "hUSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hLSS_KPi", "hLSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hLSS_PiK", "hLSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

    if (cfgMixedHistos) {
      histos.add("hUSS_KPi_Mix", "hUSS_KPi_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_PiK_Mix", "hUSS_PiK_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    }
    if (cfgMCHistos) {
      histos.add("nEvents_Gen", "nEvents_Gen", kTH1F, {{7, 0.0, 7.0}});
      histos.add("hUSS_TrueRec", "hUSS_TrueRec", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hGen_pT_Kstar", "Gen_pT_Kstar (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hRec_pT_Kstar", "Rec_pT_Kstar", kTHnSparseF, {cfgCentAxis, ptAxis});
    }
    if (cfgJetHistos) {
      histos.add("hUSS_KPi_INSIDE", "hUSS_KPi_INSIDE", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hUSS_PiK_INSIDE", "hUSS_PiK_INSIDE", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hLSS_KPi_INSIDE", "hLSS_KPi_INSIDE", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hLSS_PiK_INSIDE", "hLSS_PiK_INSIDE", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
    }
    if (cfgJetMCHistos) {
      histos.add("nEvents_Gen", "nEvents_Gen", kTH1F, {{7, 0.0, 7.0}});
      histos.add("hUSS_TrueRec", "hUSS_TrueRec", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_TrueRec_INSIDE", "hUSS_TrueRec_INSIDE", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});

      histos.add("hGen_pT_Kstar", "Gen_pT_Kstar (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hRec_pT_Kstar", "Rec_pT_Kstar (GeV/c)", kTH1F, {{800, 0., 40.}});

      histos.add("hEffRec_pT", "EffRec_pT (GeV/c)", kTH1F, {{1600, 0., 80.}});
      histos.add("hEffRecTest0_pT", "EffRecTest0_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest1_pT", "EffRecTest1_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest2_pT", "EffRecTest2_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest3_pT", "EffRecTest3_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest4_pT", "EffRecTest4_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest5_pT", "EffRecTest5_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest6_pT", "EffRecTest6_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest7_pT", "EffRecTest7_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest8_pT", "EffRecTest8_pT (GeV/c)", kTH1F, {{800, 0., 40.}});

      histos.add("hEffRecLowPtKstarDaughter", "EffRecLow Kstar daughter_pT (GeV/c)", kTH1F, {{800, 0., 40.}});

      histos.add("hEffGen_pT", "EffGen_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hMotherPdg1", "hMotherPdg1", kTH1F, {{5000, 0., 5000.}});
      histos.add("hMotherPdg2", "hMotherPdg2", kTH1F, {{5000, 0., 5000.}});
    }

    // Jet
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(cfgEventSelections));
    RealTriggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(cfgRealTriggerMasks);
    triggerMaskBitsTest1 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest1);
    triggerMaskBitsTest2 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest2);
    triggerMaskBitsTest3 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest3);
  } // end of init

  //===============
  //| LF Analysis
  //===============
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Cs>; //, aod::CentFT0Ms, aod::CentFT0As
  using EventCandidatesTrue = aod::McCollisions;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection,
                                      aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;

  //=================
  //|| Jet Analysis
  //=================
  Filter JEPosZFilter = nabs(aod::jcollision::posZ) < cfgEventVtxCut;
  Filter JEMCPosZFilter = nabs(aod::jmccollision::posZ) < cfgEventVtxCut;
  Filter JetCuts = aod::jet::pt > cfgJetpT&& aod::jet::r == nround(cfgJetR.node() * 100.0f);

  using JetTrackCandidatesMC = soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>;

  // For Mixed Event
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  Partition<TrackCandidates> kaon = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfgTrackTPCPIDnSig);
  Partition<TrackCandidates> pion = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaPi) <= cfgTrackTPCPIDnSig);
  Partition<TrackCandidatesMC> kaonMC = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaKa) <= cfgTrackTPCPIDnSig);
  Partition<TrackCandidatesMC> pionMC = !cfgTrackTPCPID || (nabs(aod::pidtpc::tpcNSigmaPi) <= cfgTrackTPCPIDnSig);

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  //==================================
  //||
  //||       Helper Templates
  //||
  //==================================
  template <typename objType>
  void fillQA(const bool pass, const objType& obj, const int objecttype = 0)
  {
    if (objecttype == 1) {
      if constexpr (requires { obj.posZ(); }) {
        if (!pass) {
          histos.fill(HIST("hPosZ_BC"), obj.posZ());
          histos.fill(HIST("hcentFT0C_BC"), obj.centFT0C());
        } else {
          histos.fill(HIST("hPosZ_AC"), obj.posZ());
          histos.fill(HIST("hcentFT0C_AC"), obj.centFT0C());
        }
      }
    } // eventSelection histogram
    if (objecttype == 2) {
      if constexpr (requires { obj.posZ(); }) {
        if (!pass) {
          histos.fill(HIST("hPosZ_BC"), obj.posZ());
          histos.fill(HIST("hcentFT0C_BC"), obj.centFT0C());
        } else {
          histos.fill(HIST("hPosZ_AC"), obj.posZ());
          histos.fill(HIST("hcentFT0C_AC"), obj.centFT0C());
        }
      }
    } // Jet eventSelection histogram
    if constexpr (requires { obj.tpcCrossedRowsOverFindableCls(); }) {
      if (objecttype == 3) {
        if (!pass) {
          histos.fill(HIST("hDCArToPv_BC"), obj.dcaXY());
          histos.fill(HIST("hDCAzToPv_BC"), obj.dcaZ());
          histos.fill(HIST("hIsPrim_BC"), obj.isPrimaryTrack());
          histos.fill(HIST("hIsGood_BC"), obj.isGlobalTrackWoDCA());
          histos.fill(HIST("hIsPrimCont_BC"), obj.isPVContributor());
          histos.fill(HIST("hFindableTPCClusters_BC"), obj.tpcNClsFindable());
          histos.fill(HIST("hFindableTPCRows_BC"), obj.tpcNClsCrossedRows());
          histos.fill(HIST("hClustersVsRows_BC"), obj.tpcCrossedRowsOverFindableCls());
          histos.fill(HIST("hTPCChi2_BC"), obj.tpcChi2NCl());
          histos.fill(HIST("hITSChi2_BC"), obj.itsChi2NCl());
          histos.fill(HIST("QA_track_pT_BC"), obj.pt());
        } else {
          histos.fill(HIST("hDCArToPv_AC"), obj.dcaXY());
          histos.fill(HIST("hDCAzToPv_AC"), obj.dcaZ());
          histos.fill(HIST("hIsPrim_AC"), obj.isPrimaryTrack());
          histos.fill(HIST("hIsGood_AC"), obj.isGlobalTrackWoDCA());
          histos.fill(HIST("hIsPrimCont_AC"), obj.isPVContributor());
          histos.fill(HIST("hFindableTPCClusters_AC"), obj.tpcNClsFindable());
          histos.fill(HIST("hFindableTPCRows_AC"), obj.tpcNClsCrossedRows());
          histos.fill(HIST("hClustersVsRows_AC"), obj.tpcCrossedRowsOverFindableCls());
          histos.fill(HIST("hTPCChi2_AC"), obj.tpcChi2NCl());
          histos.fill(HIST("hITSChi2_AC"), obj.itsChi2NCl());
          histos.fill(HIST("QA_track_pT_AC"), obj.pt());
        }
      }
    } // trackSelection
    if (objecttype == 4) {
      if constexpr (requires { obj.pt(); }) {
        if (!pass) {
          histos.fill(HIST("QA_nSigma_kaon_TPC_BC"), obj.pt(), obj.tpcNSigmaKa());
          histos.fill(HIST("QA_nSigma_kaon_TOF_BC"), obj.pt(), obj.tofNSigmaKa());
          histos.fill(HIST("QA_kaon_TPC_TOF_BC"), obj.tpcNSigmaKa(), obj.tofNSigmaKa());
        } else {
          histos.fill(HIST("QA_nSigma_kaon_TPC_AC"), obj.pt(), obj.tpcNSigmaKa());
          histos.fill(HIST("QA_nSigma_kaon_TOF_AC"), obj.pt(), obj.tofNSigmaKa());
          histos.fill(HIST("QA_kaon_TPC_TOF_AC"), obj.tpcNSigmaKa(), obj.tofNSigmaKa());
        }
      }
    } // kaon pid Selection
    if (objecttype == 5) {
      if constexpr (requires { obj.pt(); }) {
        if (!pass) {
          histos.fill(HIST("QA_nSigma_pion_TPC_BC"), obj.pt(), obj.tpcNSigmaPi());
          histos.fill(HIST("QA_nSigma_pion_TOF_BC"), obj.pt(), obj.tofNSigmaPi());
          histos.fill(HIST("QA_pion_TPC_TOF_BC"), obj.tpcNSigmaPi(), obj.tofNSigmaPi());
        } else {
          histos.fill(HIST("QA_nSigma_pion_TPC_AC"), obj.pt(), obj.tpcNSigmaPi());
          histos.fill(HIST("QA_nSigma_pion_TOF_AC"), obj.pt(), obj.tofNSigmaPi());
          histos.fill(HIST("QA_pion_TPC_TOF_AC"), obj.tpcNSigmaPi(), obj.tofNSigmaPi());
        }
      }
    } // pion pid Selection
  } // fill QA

  enum class objectType { MB,
                          Jets,
                          MBRecParticle,
                          JetsRecParticle };
  template <typename TrackType>
  void fillMinv(objectType type, const TrackType& trk1, const TrackType& trk2, const ROOT::Math::PxPyPzMVector& lReso, double centrality, double jetpt, bool IsMix, bool flip)
  {
    double conjugate = trk1.sign() * trk2.sign();
    switch (type) {
      case objectType::MB:
        if (IsMix && cfgMixedHistos) {
          if (conjugate < 0) {
            if (!flip)
              histos.fill(HIST("hUSS_KPi_Mix"), centrality, lReso.Pt(), lReso.M());
            else
              histos.fill(HIST("hUSS_PiK_Mix"), centrality, lReso.Pt(), lReso.M());
          }
        } else {
          if (conjugate < 0) {
            if (!flip)
              histos.fill(HIST("hUSS_KPi"), centrality, lReso.Pt(), lReso.M());
            else
              histos.fill(HIST("hUSS_PiK"), centrality, lReso.Pt(), lReso.M());
          } else if (conjugate > 0) {
            if (!flip)
              histos.fill(HIST("hLSS_KPi"), centrality, lReso.Pt(), lReso.M());
            else
              histos.fill(HIST("hLSS_PiK"), centrality, lReso.Pt(), lReso.M());
          }
        }
        break;

      case objectType::Jets:
        if (!IsMix && cfgJetHistos) {
          if (conjugate < 0) {
            if (!flip)
              histos.fill(HIST("hUSS_KPi_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
            else
              histos.fill(HIST("hUSS_PiK_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
          } else {
            if (conjugate > 0) {
              if (!flip)
                histos.fill(HIST("hLSS_KPi_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
              else
                histos.fill(HIST("hLSS_PiK_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
            }
          }
        }
        break;

      case objectType::MBRecParticle:
        if (cfgMCHistos) {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_TrueRec"), centrality, lReso.Pt(), lReso.M());
          }
        }
        break;

      case objectType::JetsRecParticle:
        if (cfgJetMCHistos) {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_TrueRec_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
          }
        }
        break;
    } // switch
  } // fillMinv
  //======================================================================

  template <typename EventType>
  std::pair<bool, int> eventSelection(const EventType event, const bool QA)
  {
    if (cfgEventCutQA && QA) {
      fillQA(false, event, 1);
      histos.fill(HIST("hEvent_Cut"), 0);
    }

    if (!event.sel8())
      return {false, 1};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 1);
    }
    if (std::abs(event.posZ()) > cfgEventVtxCut)
      return {false, 2};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 2);
    }
    if (!event.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return {false, 3};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 3);
    }
    if (!event.selection_bit(aod::evsel::kNoSameBunchPileup))
      return {false, 4};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 4);
    }
    if (!event.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return {false, 5};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 5);
    }
    if (!event.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return {false, 6};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 6);
    }
    if (!event.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))
      return {false, 7};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 7);
    }
    if (!event.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return {false, 8};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 8);
    }
    if (cfgOccupancySel && (event.trackOccupancyInTimeRange() > cfgOccupancyMax || event.trackOccupancyInTimeRange() < cfgOccupancyMin))
      return {false, 9};
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 9);
    }

    if (cfgEventCutQA && QA) {
      fillQA(true, event, 1);
      histos.fill(HIST("hEvent_Cut"), 10);
    }
    return {true, 11};
  };

  template <typename EventType>
  std::pair<bool, int> JeteventSelection(const EventType event, const bool QA)
  {
    if (cfgEventCutQA && QA) {
      fillQA(false, event, 2);
    }

    if (!jetderiveddatautilities::selectCollision(event, eventSelectionBits)) { // sel8
      return {false, 1};
    }
    if (std::abs(event.posZ()) > cfgEventVtxCut)
      return {false, 2};

    if (cfgEventCutQA && QA) {
      fillQA(true, event, 2);
    }
    return {true, 8};
  };

  template <typename TracksType>
  bool trackSelection(const TracksType track, const bool QA)
  {
    if (cfgTrackCutQA && QA) {
      fillQA(false, track, 3);
    }

    if (cfgTrackGlobalSel && !track.isGlobalTrack())
      return false;
    if (track.pt() < cfgTrackMinPt)
      return false;
    if (std::abs(track.eta()) > cfgTrackMaxEta)
      return false;
    if (std::abs(track.dcaXY()) > cfgTrackMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cfgTrackMaxDCAzToPVcut)
      return false;
    if (cfgTrackPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgTrackGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgTrackFindableTPCClusters > 0 && track.tpcNClsFindable() < cfgTrackFindableTPCClusters)
      return false;
    if (track.tpcNClsCrossedRows() < cfgTrackTPCCrossedRows)
      return false;
    if (cfgTrackRowsOverFindable > 0 && track.tpcCrossedRowsOverFindableCls() > cfgTrackRowsOverFindable)
      return false;
    if (track.tpcChi2NCl() > cfgTrackTPCChi2)
      return false;
    if (track.itsChi2NCl() > cfgTrackITSChi2)
      return false;
    if (cfgTrackConnectedToPV && !track.isPVContributor())
      return false;

    if (cfgTrackCutQA && QA) {
      fillQA(true, track, 3);
    }
    return true;
  };

  template <typename TrackPID>
  bool trackPIDKaon(const TrackPID& candidate, const bool QA)
  {
    double tpcpid = 0;
    double tofpid = 0;
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (cfgTrackCutQA && QA) {
      fillQA(false, candidate, 4);
    }

    // TPC
    if (cfgTrackSquarePIDCut) {
      if (std::abs(candidate.tpcNSigmaKa()) < cfgTrackTPCPIDnSig)
        tpcPIDPassed = true;
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaKa()) < cfgTrackTOFPIDnSig) {
          tofPIDPassed = true;
        }
      } else {
        if (!cfgTrackTOFHard) {
          tofPIDPassed = true;
        } else {
          tofPIDPassed = false;
        }
      }
    } // end of square cut
    if (cfgTrackCirclePIDCut) {
      if (std::abs(candidate.tpcNSigmaKa()) < cfgTrackTPCPIDnSig)
        tpcpid = std::abs(candidate.tpcNSigmaKa());
      tofpid = 0;

      if (candidate.hasTOF()) {
        tofpid = std::abs(candidate.tofNSigmaKa());
      } else {
        if (cfgTrackTOFHard) {
          tofpid = 999;
        }
      }
      if (tpcpid * tpcpid + tofpid * tofpid < cfgTrackCircleValue) {
        tpcPIDPassed = true;
        tofPIDPassed = true;
      }
    } // circular cut

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed) {
      if (cfgTrackCutQA && QA) {
        fillQA(true, candidate, 4);
      }
      return true;
    }
    return false;
  }

  template <typename TrackPID>
  bool trackPIDPion(const TrackPID& candidate, const bool QA)
  {
    double tpcpid = 0;
    double tofpid = 0;
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (cfgTrackCutQA && QA) {
      fillQA(false, candidate, 5);
    }

    // TPC
    if (cfgTrackSquarePIDCut) {
      if (std::abs(candidate.tpcNSigmaPi()) < cfgTrackTPCPIDnSig)
        tpcPIDPassed = true;
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaPi()) < cfgTrackTOFPIDnSig) {
          tofPIDPassed = true;
        }
      } else {
        if (!cfgTrackTOFHard) {
          tofPIDPassed = true;
        } else {
          tofPIDPassed = false;
        }
      }
    } // end of square cut
    if (cfgTrackCirclePIDCut) {
      if (std::abs(candidate.tpcNSigmaPi()) < cfgTrackTPCPIDnSig)
        tpcpid = std::abs(candidate.tpcNSigmaPi());
      tofpid = 0;

      if (candidate.hasTOF()) {
        tofpid = std::abs(candidate.tofNSigmaPi());
      } else {
        if (cfgTrackTOFHard) {
          tofpid = 999;
        }
      }
      if (tpcpid * tpcpid + tofpid * tofpid < cfgTrackCircleValue) {
        tpcPIDPassed = true;
        tofPIDPassed = true;
      }
    } // circular cut

    // TPC & TOF
    if (tpcPIDPassed && tofPIDPassed) {
      if (cfgTrackCutQA && QA) {
        fillQA(true, candidate, 5);
      }
      return true;
    }
    return false;
  }

  template <typename TracksType>
  ROOT::Math::PxPyPzMVector minvReconstruction(const TracksType& trk1, const TracksType& trk2, const bool QA, const bool flip)
  {
    if (!trackSelection(trk1, false) || !trackSelection(trk2, false))
      return {};
    if (!trackPIDKaon(trk1, QA) || !trackPIDPion(trk2, QA))
      return {};
    if (trk1.globalIndex() >= trk2.globalIndex())
      return {};

    ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonance;
    if (!flip) {
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
    } else {
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massKa);
    }
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Eta()) > cfgTrackMaxEta)
      return {};

    return {lResonance};
  } // minvReconstruction

  template <typename TracksType>
  ROOT::Math::PxPyPzMVector TrueReconstruction(const TracksType& trk1, const TracksType& trk2)
  {
    double conjugate = trk1.sign() * trk2.sign();
    if (conjugate > 0)
      return {};

    auto particle1 = trk1.mcParticle();
    auto particle2 = trk2.mcParticle();

    if (!particle1.has_mothers() || !particle2.has_mothers()) {
      return {};
    }
    if (std::abs(particle1.eta()) > cfgTrackMaxEta || std::abs(particle2.eta()) > cfgTrackMaxEta) {
      return {};
    }

    std::vector<int> mothers1{};
    std::vector<int> mothers1PDG{};
    for (auto& particle1_mom : particle1.template mothers_as<aod::McParticles>()) {
      mothers1.push_back(particle1_mom.globalIndex());
      mothers1PDG.push_back(particle1_mom.pdgCode());
    }

    std::vector<int> mothers2{};
    std::vector<int> mothers2PDG{};
    for (auto& particle2_mom : particle2.template mothers_as<aod::McParticles>()) {
      mothers2.push_back(particle2_mom.globalIndex());
      mothers2PDG.push_back(particle2_mom.pdgCode());
    }

    if (mothers1PDG[0] != 313)
      return {}; // mother not K*0
    if (mothers2PDG[0] != 313)
      return {}; // mothers not K*0
    if (mothers1[0] != mothers2[0])
      return {}; // Kaon and pion not from the same K*0

    if (std::abs(particle1.pdgCode()) != 211 && std::abs(particle1.pdgCode()) != 321)
      return {};
    if (std::abs(particle2.pdgCode()) != 211 && std::abs(particle2.pdgCode()) != 321)
      return {};

    double track1_mass, track2_mass;
    if (std::abs(particle1.pdgCode()) == 211) {
      track1_mass = massPi;
    } else {
      track1_mass = massKa;
    }

    if (std::abs(particle2.pdgCode()) == 211) {
      track2_mass = massPi;
    } else {
      track2_mass = massKa;
    }

    if (track1_mass == track2_mass) {
      return {};
    }

    ROOT::Math::PxPyPzMVector lTrueDaughter1, lTrueDaughter2, lTrueReso;
    lTrueDaughter1 = ROOT::Math::PxPyPzMVector(particle1.px(), particle1.py(), particle1.pz(), track1_mass);
    lTrueDaughter2 = ROOT::Math::PxPyPzMVector(particle2.px(), particle2.py(), particle2.pz(), track2_mass);
    lTrueReso = lTrueDaughter1 + lTrueDaughter2;

    if (lTrueReso.M() < 0)
      return {};
    if (std::abs(lTrueReso.Eta()) > cfgTrackMaxEta)
      return {};

    return {lTrueReso};
  } // TrueReconstruction

  template <typename CollisionType, typename TracksType>
  void TrackSlicing(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool IsMix, const bool QA)
  {
    auto tracks1 = kaon->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = pion->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      for (bool flip : {false, true}) {
        auto lReso = minvReconstruction(trk1, trk2, QA, flip);
        if (lReso.M() < 0)
          continue;

        fillMinv(objectType::MB, trk1, trk2, lReso, centrality, -1.0, IsMix, flip);
      } // flip
    } // for
  } // TrackSlicing

  template <typename CollisionType, typename TracksType>
  void TrackSlicingMC(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool IsMix, const bool QA)
  {
    auto tracks1 = kaonMC->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = pionMC->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
      if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
        continue;

      for (bool flip : {false, true}) {
        auto lReso = minvReconstruction(trk1, trk2, QA, flip);
        if (lReso.M() < 0)
          continue;

        fillMinv(objectType::MB, trk1, trk2, lReso, centrality, -1.0, IsMix, flip);
      } // flip

      //============================
      //|  True Reconstruction
      //============================
      auto lTrueReso = TrueReconstruction(trk1, trk2);
      fillMinv(objectType::MBRecParticle, trk1, trk2, lTrueReso, centrality, -1.0, IsMix, false);
    } // tracks
  } // TrackSlicingMC

  template <typename JetType>
  double DistinguishJets(const JetType& jets, ROOT::Math::PxPyPzMVector lResonance)
  {
    if (cDebugLevel > 0)
      std::cout << "Finded multiple jets to the same phi." << std::endl;

    double bestR = 0;
    double bestJetpT = 0;
    for (auto const& jet : jets) {
      double phidiff = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
      double etadiff = jet.eta() - lResonance.Eta();
      double R = TMath::Sqrt((phidiff * phidiff) + (etadiff * etadiff));
      if (R < cfgJetR && bestR == 0) {
        bestR = R;
        bestJetpT = jet.pt();
      } else if (R < bestR) {
        bestR = R;
        bestJetpT = jet.pt();
      }
    } // jet
    return bestJetpT;
  }

  template <typename TracksType, typename JetType>
  void JetTrackSlicing(aod::JetCollision const& collision, TracksType const& jetTracks, const JetType& chargedjets, const bool IsMix, const bool QA)
  {
    //=============
    //| Inclusive
    //=============
    auto centrality = collision.centFT0C();
    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsUpperIndexPolicy(jetTracks, jetTracks))) {
      auto trk1 = track1.template track_as<TrackCandidates>();
      auto trk2 = track2.template track_as<TrackCandidates>();

      for (bool flip : {false, true}) {
        auto lResonance = minvReconstruction(trk1, trk2, QA, flip);
        if (lResonance.M() < 0)
          continue;

        fillMinv(objectType::MB, trk1, trk2, lResonance, centrality, -1.0, IsMix, flip);

        //========
        //| Jets
        //========
        bool jetFlag = false;
        int goodjets = 0;
        double jetpt = 0;
        double R = 0;
        for (auto const& jet : chargedjets) {
          double phidiff = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
          double etadiff = jet.eta() - lResonance.Eta();
          R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));

          if (R < cfgJetR) {
            jetFlag = true;
            jetpt = jet.pt();
            goodjets++;
          }
        } // The loop of chargedjets
        if (cfgJetQAHistos) {
          histos.fill(HIST("nGoodJets"), goodjets);
        }

        if (!cfgSingleJet) {
          if (goodjets > 1) {
            jetpt = DistinguishJets<JetType>(chargedjets, lResonance);
          }
        }

        if (jetFlag) {
          fillMinv(objectType::Jets, trk1, trk2, lResonance, centrality, jetpt, IsMix, flip);
        } // jetFlag
      } // flip
    } // jetTracks
  } // JetTrackSlicing

  template <typename TracksType, typename JetType>
  void JetTrackSlicingMC(aod::JetCollision const& collision, TracksType const& jetTracks, const JetType& mcdjets, const bool IsMix, const bool QA)
  {
    //============================
    //| MB: Track Reconstruction
    //============================
    auto centrality = collision.centFT0C();
    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsUpperIndexPolicy(jetTracks, jetTracks))) {
      auto trk1 = track1.template track_as<TrackCandidatesMC>();
      auto trk2 = track2.template track_as<TrackCandidatesMC>();

      if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
        continue;

      for (bool flip : {false, true}) {
        auto lResonance = minvReconstruction(trk1, trk2, QA, flip);
        if (lResonance.M() < 0)
          continue;

        fillMinv(objectType::MB, trk1, trk2, lResonance, centrality, -1.0, IsMix, flip);
        //==============================
        //| Jets: Track Reconstruction
        //==============================
        bool jetFlag = false;
        int goodjets = 0;
        double jetpt = 0;
        for (auto const& jet : mcdjets) {
          double phidiff = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
          double etadiff = jet.eta() - lResonance.Eta();
          double R = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));
          if (R < cfgJetR) {
            jetFlag = true;
            jetpt = jet.pt();
            goodjets++;
          }
        }
        if (cfgJetQAHistos) {
          histos.fill(HIST("nGoodJets"), goodjets);
        }
        if (!cfgSingleJet) {
          if (goodjets > 1) {
            jetpt = DistinguishJets<JetType>(mcdjets, lResonance);
          }
        }

        if (jetFlag) {
          fillMinv(objectType::Jets, trk1, trk2, lResonance, centrality, jetpt, IsMix, flip);
        } // jetFlag
      } // filp
    } // Tracks loop
  } // JetTrackSlicingMC

  //=======================================================
  //|
  //|                 JET DATA STUFF
  //|
  //=======================================================
  using JetTrackCandidates = soa::Join<aod::JTracks, aod::JTrackPIs>;
  using JetFilteredJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;

  int nJetEvents = 0;
  void processDataJets(o2::aod::JetCollision const& collision, JetFilteredJets const& chargedjets, JetTrackCandidates const& jetTracks, TrackCandidates const&)
  {
    if (cDebugLevel > 0) {
      nJetEvents++;
      if ((nJetEvents + 1) % 10000 == 0) {
        std::cout << "Processed Jet Data Events: " << nJetEvents << std::endl;
      }
    }
    histos.fill(HIST("nEvents"), 0.5); // Raw event

    bool INELgt0 = false;
    for (auto& jetTrack : jetTracks) {
      if (std::abs(jetTrack.eta()) < cfgEventMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;
    histos.fill(HIST("nEvents"), 1.5); // INEL>0 event

    if (std::abs(collision.posZ()) > cfgEventVtxCut)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    histos.fill(HIST("nEvents"), 2.5); // After selection

    // Trigger before we start jet finding
    //=====================
    //|    Trigger
    //=====================
    if (cfgCutonTrig) {
      bool RT = false;
      bool VTtest1 = false;
      bool VTtest2 = false;
      bool VTtest3 = false;
      if (jetderiveddatautilities::selectTrigger(collision, RealTriggerMaskBits)) {
        RT = true;
        if (cfgJetQAHistos) {
          histos.fill(HIST("nTriggerQA"), 0.5);
        }
      }
      if (jetderiveddatautilities::selectTrigger(collision, triggerMaskBitsTest1)) {
        VTtest1 = true;
        if (cfgJetQAHistos) {
          histos.fill(HIST("nTriggerQA"), 1.5);
        }
      }
      if (jetderiveddatautilities::selectTrigger(collision, triggerMaskBitsTest2)) {
        VTtest2 = true;
        if (cfgJetQAHistos) {
          histos.fill(HIST("nTriggerQA"), 2.5);
        }
      }
      if (jetderiveddatautilities::selectTrigger(collision, triggerMaskBitsTest3)) {
        VTtest3 = true;
        if (cfgJetQAHistos) {
          histos.fill(HIST("nTriggerQA"), 3.5);
        }
      }
      if (!RT) {
        return;
      } else if (RT && (VTtest1 || VTtest2 || VTtest3)) {
        return;
      }
    } // Trigger cut
    histos.fill(HIST("nEvents"), 3.5); // After trigger cut

    std::vector<double> jetpT{};
    std::vector<double> jetEta{};
    std::vector<double> jetPhi{};
    bool HasJets = false;
    int nJets = 0;
    for (auto chargedjet : chargedjets) {
      jetpT.push_back(chargedjet.pt());
      jetEta.push_back(chargedjet.eta());
      jetPhi.push_back(chargedjet.phi());
      nJets++;
      if (cfgJetQAHistos) {
        histos.fill(HIST("JetpT"), chargedjet.pt());
        histos.fill(HIST("JetEta"), chargedjet.eta());
        histos.fill(HIST("JetPhi"), chargedjet.phi());
      }
      if (chargedjet.pt() > cfgJetpT)
        HasJets = true;
    }
    if (cfgJetQAHistos) {
      histos.fill(HIST("nJetsPerEvent"), nJets);
    }

    //====================
    //|    Has Jets
    //====================
    if (cfgReqJets) {
      if (!HasJets)
        return;
    }
    histos.fill(HIST("nEvents"), 4.5); // Has jets

    for (auto& jetTrack : jetTracks) {
      auto originTrack = jetTrack.track_as<TrackCandidates>();
      if (!trackSelection(originTrack, true))
        continue;

      if (cfgJetQAHistos) {
        histos.fill(HIST("rawDimpT"), jetTrack.pt(), jetTrack.pt() - originTrack.pt());
        histos.fill(HIST("jetTrackEta"), jetTrack.eta());
        histos.fill(HIST("jetTrackPhi"), jetTrack.phi());
      }
    } // jetTrack loop

    JetTrackSlicing(collision, jetTracks, chargedjets, false, true);
  } // ProcessDataJets
  PROCESS_SWITCH(kstarInOO, processDataJets, "process Data Jets", false);

  //=======================================================
  //|
  //|                 JET MC STUFF
  //|
  //=======================================================
  int nJetMCEvents = 0;
  void processMCJets(o2::aod::JetCollision const& collision, JetTrackCandidatesMC const& jetTracks, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, TrackCandidatesMC const&, aod::McParticles const&, aod::JetParticles const&)
  {
    if (cDebugLevel > 0) {
      nJetMCEvents++;
      if ((nJetMCEvents + 1) % 10000 == 0) {
        std::cout << "Processed Jet MC Events: " << nJetMCEvents << std::endl;
      }
    }
    histos.fill(HIST("nEvents"), 0.5); // Gen event

    bool INELgt0 = false;
    for (auto& jetTrack : jetTracks) {
      if (std::abs(jetTrack.eta()) < cfgEventMaxEta) {
        INELgt0 = true;
        break;
      }
    } // jetTrack loop
    if (!INELgt0)
      return;
    histos.fill(HIST("nEvents"), 1.5); // INEL>0 event

    if (std::abs(collision.posZ()) > cfgEventVtxCut)
      return;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits))
      return;
    histos.fill(HIST("nEvents"), 2.5); // Inclusive event

    // The trigger option doesn't work
    // if (!jetderiveddatautilities::selectTrigger(collision, RealTriggerMaskBits))
    //   return;
    // histos.fill(HIST("nEvents"), 3.5); // Events for Inclusive

    std::vector<double> mcdjetpT{};
    std::vector<double> mcdjetEta{};
    std::vector<double> mcdjetPhi{};

    bool HasJets = false;
    int nJets = 0;
    for (auto mcdjet : mcdjets) {
      mcdjetpT.push_back(mcdjet.pt());
      mcdjetEta.push_back(mcdjet.eta());
      mcdjetPhi.push_back(mcdjet.phi());
      nJets++;
      if (cfgJetQAHistos) {
        histos.fill(HIST("JetpT"), mcdjet.pt());
        histos.fill(HIST("JetEta"), mcdjet.eta());
        histos.fill(HIST("JetPhi"), mcdjet.phi());
      }
      if (mcdjet.pt() > cfgJetpT) {
        HasJets = true;
      }
    }
    if (cfgJetQAHistos) {
      histos.fill(HIST("nJetsPerEvent"), nJets);
    }

    //====================
    //|    Has Jets
    //====================
    if (cfgReqJets) {
      if (!HasJets) {
        return;
      }
    }
    histos.fill(HIST("nEvents"), 3.5); // Jet event

    //  JetTrackSlicingMC(collision, jetTracks, mcdjets, false, true);

    // Instead of using "JetTrackSlicingMC", Once we have to test showing the distribution each condition
    for (auto& [track1, track2] : combinations(o2::soa::CombinationsUpperIndexPolicy(jetTracks, jetTracks))) {
      auto trk1 = track1.track_as<TrackCandidatesMC>();
      auto trk2 = track2.track_as<TrackCandidatesMC>();

      // Each section, test pT. which section we lost the entries
      ROOT::Math::PxPyPzMVector lDecayDaughterTest1, lDecayDaughterTest2, lResonanceTest1;
      lDecayDaughterTest1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughterTest2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
      lResonanceTest1 = lDecayDaughterTest1 + lDecayDaughterTest2;

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest0_pT"), lResonanceTest1.Pt());
      }

      if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
        continue;
      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest1_pT"), lResonanceTest1.Pt());
      }

      //////////////////////////////
      //////////////////////////////
      if (cfgReqMcEffTrackQA) { // false, true
        if (!trackSelection(trk1, true) || !trackSelection(trk2, false))
          continue;
      }
      //////////////////////////////
      //////////////////////////////

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest2_pT"), lResonanceTest1.Pt());
      }

      if (cfgReqMcEffPID) { // false, true
        if (!trackPIDKaon(trk1, true) || !trackPIDPion(trk2, true))
          continue;
      }
      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest3_pT"), lResonanceTest1.Pt());
      }

      auto particle1 = trk1.mcParticle();
      auto particle2 = trk2.mcParticle();
      ////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////
      if (!particle1.has_mothers() || !particle2.has_mothers()) {
        continue;
      }
      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest4_pT"), lResonanceTest1.Pt());
      }
      ////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////

      std::vector<int> mothers1{};
      std::vector<int> mothers1PDG{};
      for (auto& particle1_mom : particle1.template mothers_as<aod::McParticles>()) {
        mothers1.push_back(particle1_mom.globalIndex());
        mothers1PDG.push_back(particle1_mom.pdgCode());
      }

      if (cfgJetMCHistos) {
        histos.fill(HIST("hMotherPdg1"), std::abs(mothers1PDG[0]));
      }
      std::vector<int> mothers2{};
      std::vector<int> mothers2PDG{};
      for (auto& particle2_mom : particle2.template mothers_as<aod::McParticles>()) {
        mothers2.push_back(particle2_mom.globalIndex());
        mothers2PDG.push_back(particle2_mom.pdgCode());
      }
      if (cfgJetMCHistos) {
        histos.fill(HIST("hMotherPdg2"), std::abs(mothers2PDG[0]));
      }

      if (mothers1[0] != mothers2[0])
        continue; // Kaon and pion not from the same K*0

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest5_pT"), lResonanceTest1.Pt());
      }

      if (std::abs(particle1.pdgCode()) != 321) // kaon
        continue;

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest6_pT"), lResonanceTest1.Pt());
      }

      if (std::abs(particle2.pdgCode()) != 211) // pion
        continue;

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest7_pT"), lResonanceTest1.Pt());
      }

      if (std::abs(mothers1PDG[0]) != 313)
        continue; // mother not K*0
      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest8_pT"), lResonanceTest1.Pt());
        if (lResonanceTest1.Pt() < 0.4) {
          histos.fill(HIST("hEffRecLowPtKstarDaughter"), lDecayDaughterTest1.Pt());
          histos.fill(HIST("hEffRecLowPtKstarDaughter"), lDecayDaughterTest2.Pt());
        }
      }

      if (std::abs(mothers2PDG[0]) != 313)
        continue; // mothers not K*0

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRec_pT"), lResonanceTest1.Pt());
      }

    } // track loop
  } // process loop
  PROCESS_SWITCH(kstarInOO, processMCJets, "process MC Jets", false);

  //=======================================================
  //|
  //|                  DATA STUFF (SE)
  //|
  //=======================================================
  int nEvents = 0;
  void processDataSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks)
  {
    if (cDebugLevel > 0) {
      nEvents++;
      if ((nEvents + 1) % 10000 == 0) {
        std::cout << "Processed Data Events: " << nEvents << std::endl;
      }
    }
    histos.fill(HIST("nEvents"), 0.5);

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (std::abs(track.eta()) < cfgEventMaxEta) { // cfgTrackMaxEta --> cfgEventMaxEta
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;
    histos.fill(HIST("nEvents"), 1.5);

    auto [goodEv, code] = eventSelection(collision, true);
    if (!goodEv)
      return;
    histos.fill(HIST("nEvents"), 2.5);

    // for trackQA plot
    for (const auto& track : tracks) {
      if (!trackSelection(track, true))
        continue;
    }

    TrackSlicing(collision, tracks, collision, tracks, false, true);

  } // processSameEvents
  PROCESS_SWITCH(kstarInOO, processDataSameEvent, "process Data Same Event", false);

  //=======================================================
  //|
  //|                  DATA STUFF (ME)
  //|
  //=======================================================
  int nEventsMix = 0;
  void processDataMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType colBinning{{cfgBinsMixVtx, cfgBinsMixMult}, true}; // true is for 'ignore overflows' (true by default)
    SameKindPair<EventCandidates, TrackCandidates, BinningType> pairs{colBinning, cfgMixNMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        nEventsMix++;
        if ((nEventsMix + 1) % 10000 == 0) {
          std::cout << "Processed DATA Mixed Events : " << nEventsMix << std::endl;
        }
      }
      auto [goodEv1, code1] = eventSelection(collision1, false);
      auto [goodEv2, code2] = eventSelection(collision2, false);
      bool VtxMixFlag = false;
      bool CentMixFlag = false;
      // bool OccupanacyMixFlag = false;
      if (std::abs(collision1.posZ() - collision2.posZ()) <= cfgVtxMixCut) // set default to maybe 10
        VtxMixFlag = true;
      if (std::abs(collision1.centFT0C() - collision2.centFT0C()) <= cfgVtxMixCut) // set default to maybe 10
        CentMixFlag = true;

      if (!goodEv1 || !goodEv2)
        continue;
      if (!CentMixFlag)
        continue;
      if (!VtxMixFlag)
        continue;

      TrackSlicing(collision1, tracks1, collision2, tracks2, true, false);
    }
  }
  PROCESS_SWITCH(kstarInOO, processDataMixedEvent, "process DATA Mixed Event", false);

  //=======================================================
  //|
  //|                  MC STUFF (SE)
  //|
  //=========================================================
  int nEventsMC = 0;
  void processSameEventMC(EventCandidates::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    if (cDebugLevel > 0) {
      nEventsMC++;
      if ((nEventsMC + 1) % 10000 == 0) {
        double histmem = histos.getSize();
        std::cout << histmem << std::endl;
        std::cout << "process_SameEvent_MC: " << nEventsMC << std::endl;
      }
    }
    histos.fill(HIST("nEvents"), 0.5);

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (std::abs(track.eta()) < cfgEventMaxEta) {
        INELgt0 = true;
        break;
      }
    }
    if (!INELgt0)
      return;
    histos.fill(HIST("nEvents"), 1.5);

    auto [goodEv, code] = eventSelection(collision, true); // sel8 & vtx
    if (!goodEv)
      return;
    histos.fill(HIST("nEvents"), 2.5);

    // for trackQA plot
    for (const auto& track : tracks) {
      if (!trackSelection(track, true))
        continue;
    }

    TrackSlicingMC(collision, tracks, collision, tracks, false, true);
  } // processSameEvents_MC
  PROCESS_SWITCH(kstarInOO, processSameEventMC, "process Same Event MC", false);

  //=======================================================
  //|
  //|                  MC STUFF (ME)
  //|
  //=======================================================
  int nEventsMCMix = 0;
  void processMixedEventMC(EventCandidates const& collisions, TrackCandidatesMC const& tracks, aod::McParticles const&)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType colBinning{{cfgBinsMixVtx, cfgBinsMixMult}, true}; // true is for 'ignore overflows' (true by default)
    SameKindPair<EventCandidates, TrackCandidatesMC, BinningType> pairs{colBinning, cfgMixNMixedEvents, -1, collisions, tracksTuple, &cache};
    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (cDebugLevel > 0) {
        nEventsMCMix++;
        if ((nEventsMCMix + 1) % 10000 == 0) {
          std::cout << "Processed Mixed Events: " << nEventsMCMix << std::endl;
        }
      }
      auto [goodEv1, code1] = eventSelection(collision1, false);
      auto [goodEv2, code2] = eventSelection(collision2, false);
      if (!goodEv1 || !goodEv2) {
        continue;
      }

      TrackSlicingMC(collision1, tracks1, collision2, tracks2, true, false);
    } // mixing
  } // processMixedEvent_MC
  PROCESS_SWITCH(kstarInOO, processMixedEventMC, "process Mixed Event MC", false);

  //======================================
  //|
  //|           GENERATED STUFF
  //|
  //======================================
  int nEventsGen = 0;
  void processGen(EventCandidatesTrue::iterator const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, EventCandidates>> const& recocolls, aod::McParticles const& mcParticles)
  {
    if (cDebugLevel > 0) {
      ++nEventsGen;
      if (nEventsGen % 10000 == 0) {
        std::cout << "Processed MC (GEN) Events: " << nEventsGen << std::endl;
      }
    }
    if (cfgMCHistos) {
      histos.fill(HIST("nEvents_Gen"), 0.5); // Gen events
    }

    //=======================
    //| Event & Signal loss
    //=======================
    bool INELgt0 = false;
    for (auto& particle : mcParticles) {
      if (std::abs(particle.eta()) > cfgEventMaxEta)
        continue;
      if (particle.pt() <= 0.0)
        continue;
      INELgt0 = true;
      break;
    }
    if (cfgForceTrueINELgt) {
      if (!INELgt0) {
        return;
      }
    }
    if (cfgMCHistos) {
      histos.fill(HIST("nEvents_Gen"), 1.5); // INEL>0 Gen events & EL: denominator
    }

    if (std::abs(collision.posZ()) > cfgEventVtxCut)
      return;

    for (auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != 313)
        continue;
      if (std::abs(particle.eta()) > cfgTrackMaxEta)
        continue;

      if (cfgMCHistos) {
        histos.fill(HIST("hGen_pT_Kstar"), particle.pt()); // SL: denominator
      }
    } // Gen colls

    if (recocolls.size() <= 0) { // not reconstructed
      if (cfgForceGenReco) {
        return;
      }
    }

    //=================
    //| Efficiency
    //=================
    double centrality = -1;
    bool hasGoodEv = false;
    for (auto& recocoll : recocolls) { // poorly reconstructed
      auto [goodEv, code] = eventSelection(recocoll, true);
      if (recocoll.posZ() > cfgEventVtxCut)
        goodEv = false;
      if (!goodEv)
        continue;

      hasGoodEv = true;
      centrality = recocoll.centFT0C();
      break;
    } // recocolls
    if (!hasGoodEv)
      return;
    if (cfgMCHistos) {
      histos.fill(HIST("nEvents_Gen"), 2.5); // EL: numerator
    }

    for (auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != 313)
        continue; // Not K*0
      if (std::abs(particle.eta()) > cfgTrackMaxEta)
        continue;
      if (particle.pt() < cfgTrackMinPt)
        continue;

      if (cfgMCHistos) {
        histos.fill(HIST("hRec_pT_Kstar"), centrality, particle.pt()); // SL: numerator // eff: denominator
      }

    } // loop over particles
  } // processGen
  PROCESS_SWITCH(kstarInOO, processGen, "process Generated Particles", false);

  //==============================================
  //|
  //|    GENERATED STUFF (INCLUSIVE & JETS)
  //|
  //==============================================
  int nprocessGenEvents = 0;
  void processGenJets(o2::aod::JetMcCollision const& collision, soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JetCollisions>> const& recocolls, aod::JetParticles const& mcParticles)
  {
    if (cDebugLevel > 0) {
      ++nprocessGenEvents;
      if (nprocessGenEvents % 10000 == 0) {
        std::cout << "Processed MC (GEN) Events: " << nprocessGenEvents << std::endl;
      }
    }
    if (cfgJetMCHistos) {
      histos.fill(HIST("nEvents_Gen"), 0.5);
    }
    //=======================
    //| Event & Signal loss
    //=======================
    bool INELgt0 = false;
    for (auto& particle : mcParticles) {
      if (std::abs(particle.eta()) > cfgEventMaxEta)
        continue;
      if (particle.pt() <= 0.0)
        continue;
      INELgt0 = true;
      break;
    }
    if (cfgForceTrueINELgt) {
      if (!INELgt0) {
        return;
      }
    }
    if (cfgJetMCHistos) {
      histos.fill(HIST("nEvents_Gen"), 1.5); // EL: denominator
    }

    if (std::abs(collision.posZ()) > cfgEventVtxCut)
      return;

    for (auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != 313)
        continue;
      if (std::abs(particle.eta()) > cfgTrackMaxEta)
        continue;

      if (cfgJetMCHistos) {
        histos.fill(HIST("hGen_pT_Kstar"), particle.pt()); // SL: denominator
      }
    } // Unreco.

    if (recocolls.size() <= 0) { // not reconstructed
      return;
    }
    //=================
    //| Efficiency
    //=================
    for (auto& recocoll : recocolls) { // poorly reconstructed
      auto goodEv = jetderiveddatautilities::selectCollision(recocoll, eventSelectionBits);
      if (goodEv) {
        goodEv = jetderiveddatautilities::selectTrigger(recocoll, RealTriggerMaskBits);
      }
      if (std::abs(recocoll.posZ()) > cfgEventVtxCut)
        goodEv = false;

      if (!goodEv)
        return;
    } // reco.coll
    if (cfgJetMCHistos) {
      histos.fill(HIST("nEvents_Gen"), 2.5); // EL: numerator
    }

    for (auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != 313)
        continue;
      if (std::abs(particle.eta()) > cfgTrackMaxEta)
        continue;
      if (particle.pt() < cfgTrackMinPt)
        continue;

      /* // Not Yet
   if (cfg_Force_BR) {
   bool baddecay = false;
   for (auto& phidaughter : particle.daughters_as<aod::McParticles>()) {
   if (std::abs(phidaughter.pdgCode()) != 321) {
   baddecay = true;
   break;
   }
   if (cfg_Force_Kaon_Acceptence) {
   if (std::abs(phidaughter.eta()) > cfg_Track_MaxEta) {
   baddecay = true;
   break;
   }
   }
   } // loop over daughters

   if (baddecay)
   continue;
   } // enforce BR restriction
      */

      if (cfgJetMCHistos) {
        histos.fill(HIST("hRec_pT_Kstar"), particle.pt()); // SL: numerator and eff: denominator
      }

    } // loop over particles
  } // end of process
  PROCESS_SWITCH(kstarInOO, processGenJets, "Process Generated Particles Inclusive&Jets", false);

  int ndRtest = 0;
  void processJetQA(o2::aod::JetMcCollision const& collision, soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets, aod::JetParticles const& mcParticles, aod::McParticles const&)
  {
    if (cDebugLevel > 0) {
      ++ndRtest;
      if (ndRtest % 10000 == 0) {
        std::cout << "Processed dR test: " << ndRtest << std::endl;
      }
    }

    bool INELgt0 = false;
    for (auto& mcpjet : mcpjets) {
      if (std::abs(mcpjet.eta()) > cfgEventMaxEta)
        continue;
      if (mcpjet.pt() <= 0.0)
        continue;
      INELgt0 = true;
      break;
    } // Selection event through mcpjet loop
    if (!INELgt0)
      return;

    if (std::abs(collision.posZ()) > cfgEventVtxCut)
      return;
    histos.fill(HIST("nEvents"), 0.5);

    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.eta()) > cfgTrackMaxEta)
        continue;
      if (!mcParticle.has_daughters())
        continue;

      ROOT::Math::PxPyPzEVector lResonance;
      lResonance = ROOT::Math::PxPyPzEVector(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());

      int GenPID = 0;
      if (!cfgIsKstar)
        GenPID = 333;
      else
        GenPID = 313;

      if (std::abs(mcParticle.pdgCode()) != GenPID)
        continue;

      if (cfgJetdRHistos) {
        histos.fill(HIST("Gen_particle_All_BR"), mcParticle.pt());
      }

      bool skip = false;
      int daughter_kaon = 0;
      int daughter_pion = 0;
      if (!cfgIsKstar) {
        for (auto& daughter : mcParticle.daughters_as<aod::JetParticles>()) {
          if (std::abs(daughter.pdgCode()) != 321)
            skip = true;
        }
      } // phi(1020)
      else {
        for (auto& daughter : mcParticle.daughters_as<aod::JetParticles>()) {
          if (std::abs(daughter.pdgCode()) == 321)
            ++daughter_kaon;
          else if (std::abs(daughter.pdgCode()) == 211)
            ++daughter_pion;
        }
        if (daughter_kaon != 1 || daughter_pion != 1)
          skip = true;
      } // K*(892)

      if (skip && cfgBR)
        continue;

      if (cfgJetdRHistos) {
        histos.fill(HIST("Gen_particle_BR"), lResonance.M());
        histos.fill(HIST("Gen_particle_pT"), mcParticle.pt());
      }

      //==================
      // Distinguish Jets
      double bestR = 999;
      double bestJetpT = 0;
      double bestJetPhi = 0;
      double bestJetEta = 0;
      for (auto& mcpjet : mcpjets) {
        if (mcpjet.pt() < cfgJetpT)
          continue;

        if (cfgJetdRHistos) {
          histos.fill(HIST("mcpjet_eta"), mcpjet.eta());
          histos.fill(HIST("mcpjet_phi"), mcpjet.phi());
          histos.fill(HIST("mcpjet_pt"), mcpjet.pt());
        }

        double dphi = TVector2::Phi_mpi_pi(mcpjet.phi() - lResonance.Phi());
        double deta = mcpjet.eta() - lResonance.Eta();
        double R = TMath::Sqrt((dphi * dphi) + (deta * deta));
        if (R < bestR) {
          bestR = R;
          bestJetpT = mcpjet.pt();
          bestJetPhi = mcpjet.phi();
          bestJetEta = mcpjet.eta();
        }
      } // mcpJets
      if (bestR > cfgJetR)
        continue;

      //==================
      // daughters
      double missing_pt = 0;
      double dR_kaon, dR_pion;
      bool kaon_out = false;
      bool pion_out = false;
      for (auto& daughter : mcParticle.daughters_as<aod::JetParticles>()) {
        if (cfgIsKstar) {
          if (std::abs(daughter.pdgCode()) == 321) {

            double dphi_kaon = TVector2::Phi_mpi_pi(bestJetPhi - daughter.phi());
            double deta_kaon = bestJetEta - daughter.eta();
            dR_kaon = TMath::Sqrt((dphi_kaon * dphi_kaon) + (deta_kaon * deta_kaon));

            if (bestR < cfgJetR) {
              if (cfgJetdRHistos) {
                histos.fill(HIST("dR_taggedjet_kaon"), dR_kaon, lResonance.Pt());
                histos.fill(HIST("dR_taggedjet_all"), dR_kaon, lResonance.Pt());
              }
              if (dR_kaon > cfgJetR) {
                kaon_out = true;
                missing_pt += daughter.pt();
              }
            } // INSIDE Jets
          } // kaon daughter
          if (std::abs(daughter.pdgCode()) == 211) {

            double dphi_pion = TVector2::Phi_mpi_pi(bestJetPhi - daughter.phi());
            double deta_pion = bestJetEta - daughter.eta();
            dR_pion = TMath::Sqrt((dphi_pion * dphi_pion) + (deta_pion * deta_pion));

            if (bestR < cfgJetR) {
              if (cfgJetdRHistos) {
                histos.fill(HIST("dR_taggedjet_pion"), dR_pion, lResonance.Pt());
                histos.fill(HIST("dR_taggedjet_all"), dR_pion, lResonance.Pt());

                if (bestJetpT > 6.0 && bestJetpT < 8.0)
                  histos.fill(HIST("dR_taggedjet_all_6_8"), dR_pion, lResonance.Pt());
              }
              if (dR_pion > cfgJetR) {
                pion_out = true;
                missing_pt += daughter.pt();
              }
            } // INSIDE Jets
          } // pion daughter
        } // K*(892)0
        else {
          if (std::abs(daughter.pdgCode()) == 321) {
            double dphi_kaon = TVector2::Phi_mpi_pi(bestJetPhi - daughter.phi());
            double deta_kaon = bestJetEta - daughter.eta();
            dR_kaon = TMath::Sqrt((dphi_kaon * dphi_kaon) + (deta_kaon * deta_kaon));

            if (bestR < cfgJetR) {
              if (cfgJetdRHistos) {
                histos.fill(HIST("dR_taggedjet_kaon"), dR_kaon, lResonance.Pt());
                histos.fill(HIST("dR_taggedjet_all"), dR_kaon, lResonance.Pt());

                if (bestJetpT > 6.0 && bestJetpT < 8.0)
                  histos.fill(HIST("dR_taggedjet_all_6_8"), dR_kaon, lResonance.Pt());
              }

              if (dR_kaon > cfgJetR) {
                kaon_out = true;
                missing_pt = daughter.pt();
              }
            }
          } // kaon daughter
        } // phi(1020)
      } // daughter

      if (kaon_out || pion_out) {
        double recoveredJetpT = bestJetpT + missing_pt;
        if (cfgJetdRHistos) {
          if (bestJetpT > 6.0 && bestJetpT < 8.0) {
            histos.fill(HIST("missed_kpi_INJets_6_8"), (bestJetpT - missing_pt) / bestJetpT, lResonance.Pt());
            if (recoveredJetpT > 8.0)
              histos.fill(HIST("recoveredJetpT_6_8to8_10"), recoveredJetpT);
          }
          if (bestJetpT > 8.0 && bestJetpT < 10.0)
            histos.fill(HIST("missed_kpi_INJets_8_10"), (bestJetpT - missing_pt) / bestJetpT, lResonance.Pt());
          if (bestJetpT > 10.0 && bestJetpT < 12.0)
            histos.fill(HIST("missed_kpi_INJets_10_12"), (bestJetpT - missing_pt) / bestJetpT, lResonance.Pt());
          if (bestJetpT > 12.0 && bestJetpT < 15.0)
            histos.fill(HIST("missed_kpi_INJets_12_15"), (bestJetpT - missing_pt) / bestJetpT, lResonance.Pt());
          if (bestJetpT > 15.0 && bestJetpT < 25.0)
            histos.fill(HIST("missed_kpi_INJets_15_25"), (bestJetpT - missing_pt) / bestJetpT, lResonance.Pt());
          if (bestJetpT > 25.0)
            histos.fill(HIST("missed_kpi_INJets_25_infinite"), (bestJetpT - missing_pt) / bestJetpT, lResonance.Pt());

          if (bestJetpT > 8.0)
            histos.fill(HIST("missed_kpi_INJets_8_infinite"), (bestJetpT - missing_pt) / bestJetpT, lResonance.Pt());

          histos.fill(HIST("JetMigration"), bestJetpT, recoveredJetpT);

        } // cfgJetdRHistos
      } // kaon_out || pion_out
    } // mcParticles
  };

  PROCESS_SWITCH(kstarInOO, processJetQA, "Process dR of K*0 Inclusive and Inside jet", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarInOO>(cfgc)};
}
