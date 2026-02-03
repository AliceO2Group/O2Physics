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
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include <Framework/ASoAHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

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
  //||
  //||         Selection
  //||
  //==================================

  // Event Selection
  Configurable<std::string> cfgEventSelections{"cfgEventSelections", "sel8", "Set event selection"};
  Configurable<float> cfgEventVtxCut{"cfgEventVtxCut", 10.0, "V_z cut selection"};
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
  Configurable<bool> cfgTrackPrimaryTrack{"cfgTrackPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgTrackConnectedToPV{"cfgTrackConnectedToPV", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgTrackGlobalWoDCATrack{"cfgTrackGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
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

  Configurable<bool> cfgDataHistos{"cfgDataHistos", false, "Enable Data Hists"};
  Configurable<bool> cfgMcHistos{"cfgMcHistos", false, "Enable MC Hists"};

  Configurable<bool> cfgJetDataHistos{"cfgJetDataHistos", false, "Enable Jet Data Histos"};
  Configurable<bool> cfgJetMCHistos{"cfgJetMCHistos", false, "Enable Jet MC Histos"};
  Configurable<bool> cfgCutonTrig{"cfgCutonTrig", false, "Enable Jet Cut on Trig"};

  //======================
  //||
  //||        JET
  //||
  //======================
  Configurable<float> cfgJetpT{"cfgJetpT", 8.0, "Set Jet pT minimum"};
  Configurable<float> cfgJetR{"cfgJetR", 0.4, "Set Jet radius parameter"};
  Configurable<bool> cfgSingleJet{"cfgSingleJet", false, "Enforces strict phi-jet correspondance"};
  Configurable<bool> cfgReqJets{"cfgReqJets", false, "False: MB, True: Inside Jets"};
  Configurable<std::string> cfgRealTriggerMasks{"cfgRealTriggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<std::string> cfgTriggerMasksTest1{"cfgTriggerMasksTest1", "", "possible JE Trigger masks Test1"};
  Configurable<std::string> cfgTriggerMasksTest2{"cfgTriggerMasksTest2", "", "possible JE Trigger masks Test2"};
  Configurable<std::string> cfgTriggerMasksTest3{"cfgTriggerMasksTest3", "", "possible JE Trigger masks Test3"};

  std::vector<int> eventSelectionBits;
  std::vector<int> RealTriggerMaskBits;
  std::vector<int> triggerMaskBitsTest1;
  std::vector<int> triggerMaskBitsTest2;
  std::vector<int> triggerMaskBitsTest3;
  // Main
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
      histos.add("hEvent_Cut", "Number of event after cuts", kTH1D, {{12, 0, 12}});
      histos.add("hPosZ_BC", "hPosZ_BC", kTH1F, {{300, -15.0, 15.0}});
      histos.add("hPosZ_AC", "hPosZ_AC", kTH1F, {{300, -15.0, 15.0}});
      histos.add("hcentFT0C_BC", "centFT0C_BC", kTH1F, {{110, 0.0, 110.0}});
      histos.add("hcentFT0C_AC", "centFT0C_AC", kTH1F, {{110, 0.0, 110.0}});
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

    if (cfgDataHistos) {
      histos.add("nEvents", "nEvents", kTH1F, {{4, 0.0, 4.0}});
      histos.add("hUSS_KPi", "hUSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_PiK", "hUSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hLSS_KPi", "hLSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hLSS_PiK", "hLSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_KPi_Mix", "hUSS_KPi_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_PiK_Mix", "hUSS_PiK_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    }

    if (cfgMcHistos) {
      histos.add("nEvents_MC", "nEvents_MC", kTH1F, {{4, 0.0, 4.0}});
      histos.add("nEvents_MC_True", "nEvents_MC_True", kTH1F, {{4, 0.0, 4.0}});
      histos.add("hMC_kstar_True", "hMC_kstar_True", kTHnSparseF, {cfgCentAxis, ptAxis});

      histos.add("hMC_USS_True", "hMC_USS_True", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_USS_KPi", "hMC_USS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_USS_PiK", "hMC_USS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_LSS_KPi", "hMC_LSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_LSS_PiK", "hMC_LSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

      histos.add("hMC_USS_KPi_Mix", "hMC_USS_KPi_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_USS_PiK_Mix", "hMC_USS_PiK_Mix", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

      histos.add("hMC_USS_KPi_True", "hMC_USS_KPi_True", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMC_USS_PiK_True", "hMC_USS_PiK_True", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    }

    if (cfgJetQAHistos) {
      histos.add("nTriggerQA", "nTriggerQA", kTH1F, {{7, 0.0, 7.0}});
      histos.add("JetpT", "Jet pT (GeV/c)", kTH1F, {{4000, 0., 200.}});
      histos.add("JetEta", "Jet Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("JetPhi", "Jet Phi", kTH1F, {{80, -1.0, 7.0}});

      histos.add("rawDimpT", "rawDimpT", kTH2F, {{1000, 0.0, 10.0}, {100, -0.5, 0.5}});
      histos.add("jetTrackEta", "Jet Track Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("jetTrackPhi", "Jet Track Phi", kTH1F, {{80, -1.0, 7.0}});

      histos.add("nJetsPerEvent", "nJetsPerEvent", kTH1I, {{4, -0.5, 3.5}});
      histos.add("nGoodJets", "nGoodJets", kTH1I, {{4, -0.5, 3.5}});
    }
    if (cfgJetDataHistos) {
      histos.add("nEvents", "nEvents", kTH1F, {{7, 0.0, 7.0}});

      histos.add("hMB_USS_KPi", "hMB_USS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMB_USS_PiK", "hMB_USS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMB_LSS_KPi", "hMB_LSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMB_LSS_PiK", "hMB_LSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

      histos.add("hUSS_INSIDE_KPi", "hUSS_INSIDE_KPi", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hUSS_INSIDE_PiK", "hUSS_INSIDE_PiK", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hLSS_INSIDE_KPi", "hLSS_INSIDE_KPi", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hLSS_INSIDE_PiK", "hLSS_INSIDE_PiK", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
    }
    if (cfgJetMCHistos) {
      histos.add("nEvents_MC", "nEvents_MC", kTH1F, {{7, -.0, 7.0}});
      histos.add("nEvents_MC_True", "nEvents_MC_True", kTH1F, {{7, -.0, 7.0}});

      histos.add("hMB_USS_KPi_MC", "hMB_USS_KPi_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMB_USS_PiK_MC", "hMB_USS_PiK_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMB_LSS_KPi_MC", "hMB_LSS_KPi_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hMB_LSS_PiK_MC", "hMB_LSS_PiK_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

      histos.add("hUSS_INSIDE_KPi_MC", "hUSS_INSIDE_KPi_MC", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hUSS_INSIDE_PiK_MC", "hUSS_INSIDE_PiK_MC", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hLSS_INSIDE_KPi_MC", "hLSS_INSIDE_KPi_MC", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});
      histos.add("hLSS_INSIDE_PiK_MC", "hLSS_INSIDE_PiK_MC", kTHnSparseF, {cfgCentAxis, dRAxis, ptAxis, minvAxis});

      histos.add("hUSS_True_MC", "hUSS_True_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_KPi_True_MC", "hUSS_KPi_True_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
      histos.add("hUSS_PiK_True_MC", "hUSS_PiK_True_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

      histos.add("hUSS_TrueRec_INSIDE_MC", "hUSS_TrueRec_INSIDE_MC", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

      histos.add("hEffRec_pT", "EffRec_pT (GeV/c)", kTH1F, {{1600, 0., 80.}});
      histos.add("hEffRecTest1_pT", "EffRecTest1_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest2_pT", "EffRecTest2_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest3_pT", "EffRecTest3_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest4_pT", "EffRecTest4_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest5_pT", "EffRecTest5_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest6_pT", "EffRecTest6_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest7_pT", "EffRecTest7_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffRecTest8_pT", "EffRecTest8_pT (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hEffGen_pT", "EffGen_pT (GeV/c)", kTH1F, {{800, 0., 40.}});

      histos.add("hMotherPdg1", "hMotherPdg1", kTH1F, {{5000, 0., 5000.}});
      histos.add("hMotherPdg2", "hMotherPdg2", kTH1F, {{5000, 0., 5000.}});
    }

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

    // Jet
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(cfgEventSelections));
    RealTriggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(cfgRealTriggerMasks);
    triggerMaskBitsTest1 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest1);
    triggerMaskBitsTest2 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest2);
    triggerMaskBitsTest3 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest3);
  } // end of init

  //======================
  //|| For LF Analysis
  //======================
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Cs>; //, aod::CentFT0Ms, aod::CentFT0As
  using EventCandidatesTrue = aod::McCollisions;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection,
                                      aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;

  //==============
  //|| For jets
  //==============
  Filter JEPosZFilter = nabs(aod::jcollision::posZ) < cfgEventVtxCut;
  Filter JEMCPosZFilter = nabs(aod::jmccollision::posZ) < cfgEventVtxCut;
  Filter jetCuts = aod::jet::pt > cfgJetpT&& aod::jet::r == nround(cfgJetR.node() * 100.0f);

  using JetTrackCandidates = soa::Join<aod::JTracks, aod::JTrackPIs>;
  using JetTrackCandidatesMC = soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>;

  using JetFilteredJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;

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
  }

  //======================================================================

  template <typename EventType>
  std::pair<bool, int> eventSelection(const EventType event, const bool QA)
  {
    if (cfgEventCutQA && QA) {
      fillQA(false, event, 1);
    }

    if (!event.sel8())
      return {false, 1};
    if (std::abs(event.posZ()) > cfgEventVtxCut)
      return {false, 2};
    if (!event.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return {false, 3};
    if (!event.selection_bit(aod::evsel::kNoSameBunchPileup))
      return {false, 4};
    if (!event.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return {false, 5};
    if (!event.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return {false, 6};
    if (!event.selection_bit(aod::evsel::kNoCollInTimeRangeStandard))
      return {false, 7};
    if (!event.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return {false, 8};
    if (cfgOccupancySel && (event.trackOccupancyInTimeRange() > cfgOccupancyMax || event.trackOccupancyInTimeRange() < cfgOccupancyMin))
      return {false, 9};

    if (cfgEventCutQA && QA) {
      fillQA(true, event, 1);
    }
    return {true, 8};
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
  }

  template <typename CollisionType, typename TracksType>
  void TrackSlicing(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&, const bool IsMix, const bool QA)
  {
    auto tracks1 = kaon->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
    auto tracks2 = pion->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
    auto centrality = collision1.centFT0C();

    for (const auto& [trk1, trk2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

      auto lResoKpi = minvReconstruction(trk1, trk2, QA, false);
      auto lResopiK = minvReconstruction(trk1, trk2, QA, true);

      if (lResoKpi.M() < 0)
        continue;

      double conjugate = trk1.sign() * trk2.sign();
      if (cfgDataHistos) {
        if (!IsMix) {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_KPi"), centrality, lResoKpi.Pt(), lResoKpi.M());
            histos.fill(HIST("hUSS_PiK"), centrality, lResopiK.Pt(), lResopiK.M());
          } else if (conjugate > 0) {
            histos.fill(HIST("hLSS_KPi"), centrality, lResoKpi.Pt(), lResoKpi.M());
            histos.fill(HIST("hLSS_PiK"), centrality, lResopiK.Pt(), lResopiK.M());
          }
        } else {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_KPi_Mix"), centrality, lResoKpi.Pt(), lResoKpi.M());
            histos.fill(HIST("hUSS_PiK_Mix"), centrality, lResopiK.Pt(), lResopiK.M());
          }
        }
      }
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

      auto lResoKpi = minvReconstruction(trk1, trk2, QA, false);
      auto lResopiK = minvReconstruction(trk1, trk2, QA, true);

      if (lResoKpi.M() < 0)
        continue;

      double conjugate = trk1.sign() * trk2.sign();
      if (cfgMcHistos) {
        if (!IsMix) {
          if (conjugate < 0) {
            histos.fill(HIST("hMC_USS_KPi"), centrality, lResoKpi.Pt(), lResoKpi.M());
            histos.fill(HIST("hMC_USS_PiK"), centrality, lResopiK.Pt(), lResopiK.M());
          } else if (conjugate > 0) {
            histos.fill(HIST("hMC_LSS_KPi"), centrality, lResoKpi.Pt(), lResoKpi.M());
            histos.fill(HIST("hMC_LSS_PiK"), centrality, lResopiK.Pt(), lResopiK.M());
          }
        } else {
          if (conjugate < 0) {
            histos.fill(HIST("hMC_USS_KPi_Mix"), centrality, lResoKpi.Pt(), lResoKpi.M());
            histos.fill(HIST("hMC_USS_PiK_Mix"), centrality, lResopiK.Pt(), lResopiK.M());
          }
        }
      }
      //======================
      // Gen MC
      std::vector<int> mcMemory;

      auto particle1 = trk1.mcParticle();
      auto particle2 = trk2.mcParticle();

      if (!particle1.has_mothers() || !particle2.has_mothers()) {
        continue;
      }
      int mcindex1 = trk1.globalIndex();
      int mcindex2 = trk2.globalIndex();

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
        continue; // mother not K*0
      if (mothers2PDG[0] != 313)
        continue; // mothers not K*0

      if (mothers1[0] != mothers2[0])
        continue; // Kaon and pion not from the same K*0

      if (std::fabs(particle1.pdgCode()) != 211 && std::fabs(particle1.pdgCode()) != 321)
        continue;
      if (std::fabs(particle2.pdgCode()) != 211 && std::fabs(particle2.pdgCode()) != 321)
        continue;

      double track1_mass, track2_mass;
      bool track1f{false}; // true means pion

      if (std::fabs(particle1.pdgCode()) == 211) {
        track1f = true;
        track1_mass = massPi;
      } else {
        track1_mass = massKa;
      }

      if (std::fabs(particle2.pdgCode()) == 211) {
        track2_mass = massPi;
      } else {
        track2_mass = massKa;
      }

      if (track1_mass == track2_mass) {
        return;
      }

      bool exists1 = std::find(mcMemory.begin(), mcMemory.end(), mcindex1) != mcMemory.end();
      bool exists2 = std::find(mcMemory.begin(), mcMemory.end(), mcindex2) != mcMemory.end();
      if (exists1 || exists2) {
        continue;
      } else {
        mcMemory.push_back(trk1.globalIndex());
        mcMemory.push_back(trk2.globalIndex());
      }

      ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonance;
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), track1_mass);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), track2_mass);
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (cfgMcHistos) {
        histos.fill(HIST("hMC_USS_True"), centrality, lResonance.Pt(), lResonance.M());
        if (track1f) {
          histos.fill(HIST("hMC_USS_PiK_True"), centrality, lResonance.Pt(), lResonance.M());
        } else {
          histos.fill(HIST("hMC_USS_KPi_True"), centrality, lResonance.Pt(), lResonance.M());
        }
      }
      //======================
    } // for
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
  ROOT::Math::PxPyPzMVector JetminvReconstruction(aod::JetCollision const& collision, const TracksType& trk1, const TracksType& trk2, const JetType& jets, const bool QA, const bool flip)
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

    // auto lResoKpi = minvReconstruction(trk1, trk2, QA, false);
    // auto lResopiK = minvReconstruction(trk1, trk2, QA, true);

    //======================
    //| MinBias Event M_inv
    //======================
    auto mult = collision.centFT0C();
    double conjugate = trk1.sign() * trk2.sign();
    if (cfgJetDataHistos) {
      if (!flip) {
        if (conjugate < 0) {
          histos.fill(HIST("hMB_USS_KPi"), mult, lResonance.Pt(), lResonance.M());
        } else if (conjugate > 0) {
          histos.fill(HIST("hMB_LSS_KPi"), mult, lResonance.Pt(), lResonance.M());
        }
      } else {
        if (conjugate < 0) {
          histos.fill(HIST("hMB_USS_PiK"), mult, lResonance.Pt(), lResonance.M());
        } else if (conjugate > 0) {
          histos.fill(HIST("hMB_LSS_PiK"), mult, lResonance.Pt(), lResonance.M());
        }
      }
    } // cfgJetDataHistos

    //======================
    //| Inside jets M_inv
    //======================
    bool jetFlag = false;
    int goodjets = 0;
    double jetpt = 0;
    for (auto const& jet : jets) {
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
        jetpt = DistinguishJets<JetType>(jets, lResonance);
      }
    }

    if (jetFlag) {
      if (cfgJetDataHistos) {
        if (!flip) {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_INSIDE_KPi"), mult, jetpt, lResonance.Pt(), lResonance.M());
          } else if (conjugate > 0) {
            histos.fill(HIST("hLSS_INSIDE_KPi"), mult, jetpt, lResonance.Pt(), lResonance.M());
          }
        } else {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_INSIDE_PiK"), mult, jetpt, lResonance.Pt(), lResonance.M());
          } else if (conjugate > 0) {
            histos.fill(HIST("hLSS_INSIDE_PiK"), mult, jetpt, lResonance.Pt(), lResonance.M());
          }
        }
      } // cfgJetDataHistos
    } // jetFlag

    if (lResonance.M() > 0.85 && lResonance.M() < 0.95) {
      if (jetFlag)
        return {};
      if (goodjets > 0)
        return {};
      return {};
    } else {
      return {};
    }
  } // JetminvReconstruction

  template <typename TracksType, typename JetType>
  std::pair<double, double> JetminvReconstructionMC(o2::aod::JetCollision const& collision, const TracksType& trk1, const TracksType& trk2, const JetType& jets, const bool QA, const bool flip)
  {
    if (!trackSelection(trk1, false) || !trackSelection(trk2, false))
      return {-1.0, -1.0};

    if (!trackPIDKaon(trk1, QA) || !trackPIDPion(trk2, QA))
      return {-1.0, -1.0};

    if (trk1.globalIndex() >= trk2.globalIndex())
      return {-1.0, -1.0};

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    if (!flip) {
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);
    } else {
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
    }
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Eta()) > cfgTrackMaxEta)
      return {-1.0, -1.0};

    //======================
    //| MinBias Event M_inv
    //======================
    auto mult = collision.centFT0C();
    double conjugate = trk1.sign() * trk2.sign();
    if (cfgJetMCHistos) {
      if (!flip) {
        if (conjugate < 0) {
          histos.fill(HIST("hMB_USS_KPi_MC"), mult, lResonance.Pt(), lResonance.M());
        } else if (conjugate > 0) {
          histos.fill(HIST("hMB_LSS_KPi_MC"), mult, lResonance.Pt(), lResonance.M());
        }
      } else {
        if (conjugate < 0) {
          histos.fill(HIST("hMB_USS_PiK_MC"), mult, lResonance.Pt(), lResonance.M());
        } else if (conjugate > 0) {
          histos.fill(HIST("hMB_LSS_PiK_MC"), mult, lResonance.Pt(), lResonance.M());
        }
      }
    } // cfgJetMCHistos

    //======================
    //| Inside jets M_inv
    //======================
    bool jetFlag = false;
    int goodjets = 0;
    double jetpt = 0;
    for (auto const& jet : jets) {
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
        jetpt = DistinguishJets<JetType>(jets, lResonance);
      }
    }

    if (jetFlag) {
      if (cfgJetMCHistos) {
        if (!flip) {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_INSIDE_KPi_MC"), mult, jetpt, lResonance.Pt(), lResonance.M());
          } else if (conjugate > 0) {
            histos.fill(HIST("hLSS_INSIDE_KPi_MC"), mult, jetpt, lResonance.Pt(), lResonance.M());
          }
        } else {
          if (conjugate < 0) {
            histos.fill(HIST("hUSS_INSIDE_PiK_MC"), mult, jetpt, lResonance.Pt(), lResonance.M());
          } else if (conjugate > 0) {
            histos.fill(HIST("hLSS_INSIDE_PiK_MC"), mult, jetpt, lResonance.Pt(), lResonance.M());
          }
        }
      } // cfgJetDataHistos
    } // jetFlag

    //======================
    //| MinBias True M_inv
    //======================
    if (trk1.has_mcParticle() && trk2.has_mcParticle()) {
      auto particle1 = trk1.mcParticle();
      auto particle2 = trk2.mcParticle();

      if (!particle1.has_mothers() || !particle2.has_mothers()) {
        return {-1.0, -1.0};
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
        return {-1.0, -1.0}; // mother not K*0
      if (mothers2PDG[0] != 313)
        return {-1.0, -1.0}; // mothers not K*0

      if (mothers1[0] != mothers2[0])
        return {-1.0, -1.0}; // Kaon and pion not from the same K*0

      if (std::fabs(particle1.pdgCode()) != 211 && std::fabs(particle1.pdgCode()) != 321)
        return {-1.0, -1.0};
      if (std::fabs(particle2.pdgCode()) != 211 && std::fabs(particle2.pdgCode()) != 321)
        return {-1.0, -1.0};

      double track1_mass, track2_mass;
      bool track1f{false}; // true means pion

      if (std::fabs(particle1.pdgCode()) == 211) {
        track1f = true;
        track1_mass = massPi;
      } else {
        track1_mass = massKa;
      }

      if (std::fabs(particle2.pdgCode()) == 211) {
        track2_mass = massPi;
      } else {
        track2_mass = massKa;
      }

      if (track1_mass == track2_mass) {
        return {-1.0, -1.0};
      }
      int mcindex1 = trk1.globalIndex();
      int mcindex2 = trk2.globalIndex();
      std::vector<int> mcMemory;

      bool exists1 = std::find(mcMemory.begin(), mcMemory.end(), mcindex1) != mcMemory.end();
      bool exists2 = std::find(mcMemory.begin(), mcMemory.end(), mcindex2) != mcMemory.end();
      if (exists1 || exists2) {
        return {-1.0, -1.0};
      } else {
        mcMemory.push_back(trk1.globalIndex());
        mcMemory.push_back(trk2.globalIndex());
      }

      TLorentzVector lTrueDaughter1, lTrueDaughter2, lTrueReso;
      lTrueDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), track1_mass);
      lTrueDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), track2_mass);
      lTrueReso = lTrueDaughter1 + lTrueDaughter2;

      if (cfgJetMCHistos) {
        histos.fill(HIST("hUSS_True_MC"), mult, lTrueReso.Pt(), lTrueReso.M());
        if (track1f) {
          histos.fill(HIST("hUSS_PiK_True_MC"), mult, lTrueReso.Pt(), lTrueReso.M());
        } else {
          histos.fill(HIST("hUSS_KPi_True_MC"), mult, lTrueReso.Pt(), lTrueReso.M());
        }
      }

      //===========================
      // INSIDE REC True Closure
      //===========================
      if (jetFlag) {
        if (conjugate < 0) {
          histos.fill(HIST("hUSS_TrueRec_INSIDE_MC"), mult, lTrueReso.Pt(), lTrueReso.M());
        }
      }

    } // has_mcParticle

    if (lResonance.M() > 0.85 && lResonance.M() < 0.95) {
      if (jetFlag)
        return {3.0, 3.0};
      if (goodjets > 0)
        return {2.0, 2.0};
      return {1.0, 1.0};
    } else {
      return {-1.0, -1.0};
    }

  } // JetMCminvReconstruction

  //=======================================================
  //|
  //|                 JET DATA STUFF
  //|
  //=======================================================
  int nJetEvents = 0;
  void processDataJets(o2::aod::JetCollision const& collision, JetFilteredJets const& chargedjets, JetTrackCandidates const& jetTracks, TrackCandidates const&)
  {
    if (cDebugLevel > 0) {
      nJetEvents++;
      if ((nJetEvents + 1) % 10000 == 0) {
        std::cout << "Processed Jet Data Events: " << nJetEvents << std::endl;
      }
    }
    if (cfgJetDataHistos) {
      histos.fill(HIST("nEvents"), 0.5); // Raw event
    }
    auto [goodEv, code] = JeteventSelection(collision, true);
    if (!goodEv)
      return;

    // Trigger before we start jet finding
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

    if (cfgJetDataHistos) {
      histos.fill(HIST("nEvents"), 1.5); // Before passing the condition
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (cfgJetDataHistos) {
      histos.fill(HIST("nEvents"), 2.5); // Events after event quality selection for Inclusive
    }

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
    //||    Has Jets
    //====================
    if (cfgReqJets) {
      if (!HasJets)
        return;
    }
    if (cfgJetDataHistos) {
      histos.fill(HIST("nEvents"), 3.5); // Has jets
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(jetTracks, jetTracks))) {
      auto trk1 = track1.track_as<TrackCandidates>();
      auto trk2 = track2.track_as<TrackCandidates>();

      JetminvReconstruction(collision, trk1, trk2, chargedjets, true, false);
      JetminvReconstruction(collision, trk1, trk2, chargedjets, true, true);
    }

    bool INELgt0 = false;
    for (auto& jetTrack : jetTracks) {
      auto originTrack = jetTrack.track_as<TrackCandidates>();
      if (!trackSelection(originTrack, true))
        continue;
      INELgt0 = true;

      if (cfgJetQAHistos) {
        histos.fill(HIST("rawDimpT"), jetTrack.pt(), jetTrack.pt() - originTrack.pt());
        histos.fill(HIST("jetTrackEta"), jetTrack.eta());
        histos.fill(HIST("jetTrackPhi"), jetTrack.phi());
      }
    } // jetTrack loop
    if (!INELgt0)
      return;
  } // ProcessDataJets
  PROCESS_SWITCH(kstarInOO, processDataJets, "process Data Jets", false);

  //=======================================================
  //|
  //|                 JET MC STUFF
  //|
  //=======================================================
  int nJetMCEvents = 0;
  void processMCJets(o2::aod::JetCollision const& collision, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets, JetTrackCandidatesMC const& jetTracks, TrackCandidatesMC const&, aod::McParticles const&)
  {
    if (cDebugLevel > 0) {
      nJetMCEvents++;
      if ((nJetMCEvents + 1) % 10000 == 0) {
        std::cout << "Processed Jet MC Events: " << nJetMCEvents << std::endl;
      }
    }

    if (cfgJetMCHistos) {
      histos.fill(HIST("nEvents_MC"), 0.5); // Raw event
    }

    if (std::abs(collision.posZ()) > cfgEventVtxCut)
      return;

    if (!jetderiveddatautilities::selectTrigger(collision, RealTriggerMaskBits))
      return;

    if (cfgJetMCHistos) {
      histos.fill(HIST("nEvents_MC"), 1.5); // Before passing the condition
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    bool INELgt0 = false;
    for (auto& jetTrack : jetTracks) {
      if (std::fabs(jetTrack.eta()) < cfgTrackMaxEta) {
        INELgt0 = true;
        break;
      }
    } // jetTrack loop
    if (!INELgt0)
      return;

    if (cfgJetMCHistos) {
      histos.fill(HIST("nEvents_MC"), 2.5); // Events after event quality selection for Inclusive
    }

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
    //||    Has Jets
    //====================
    if (cfgReqJets) {
      if (!HasJets) {
        return;
      }
    }

    if (cfgJetMCHistos) {
      histos.fill(HIST("nEvents_MC"), 3.5); // Has jets
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsUpperIndexPolicy(jetTracks, jetTracks))) {
      auto trk1 = track1.track_as<TrackCandidatesMC>();
      auto trk2 = track2.track_as<TrackCandidatesMC>();

      // Each section, test pT. which section we lost the entries
      ROOT::Math::PxPyPzMVector lDecayDaughterTest1, lDecayDaughterTest2, lResonanceTest1;
      lDecayDaughterTest1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughterTest2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
      lResonanceTest1 = lDecayDaughterTest1 + lDecayDaughterTest2;

      if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
        continue;
      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest1_pT"), lResonanceTest1.Pt());
      }

      if (!trackSelection(trk1, true) || !trackSelection(trk2, false))
        continue;
      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest2_pT"), lResonanceTest1.Pt());
      }

      if (cfgReqMcEffPID) {
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

      // std::cout<<mothers1PDG[0]<<std::endl;
      if (cfgJetMCHistos) {
        histos.fill(HIST("hMotherPdg1"), std::fabs(mothers1PDG[0]));
      }
      std::vector<int> mothers2{};
      std::vector<int> mothers2PDG{};
      for (auto& particle2_mom : particle2.template mothers_as<aod::McParticles>()) {
        mothers2.push_back(particle2_mom.globalIndex());
        mothers2PDG.push_back(particle2_mom.pdgCode());
      }
      if (cfgJetMCHistos) {
        histos.fill(HIST("hMotherPdg2"), std::fabs(mothers2PDG[0]));
      }

      if (mothers1[0] != mothers2[0])
        continue; // Kaon and pion not from the same K*0

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest5_pT"), lResonanceTest1.Pt());
      }

      if (std::fabs(particle1.pdgCode()) != 321) // kaon
        continue;

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest6_pT"), lResonanceTest1.Pt());
      }

      if (std::fabs(particle2.pdgCode()) != 211) // pion
        continue;

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest7_pT"), lResonanceTest1.Pt());
      }

      if (std::fabs(mothers1PDG[0]) != 313)
        continue; // mother not K*0
      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRecTest8_pT"), lResonanceTest1.Pt());
      }

      if (std::fabs(mothers2PDG[0]) != 313)
        continue; // mothers not K*0

      if (cfgJetMCHistos) {
        histos.fill(HIST("hEffRec_pT"), lResonanceTest1.Pt());
      }
    } // track loop
  } // process loop
  PROCESS_SWITCH(kstarInOO, processMCJets, "process MC Jets", false);

  //======================================================
  //|
  //|             Efficiency JET MC STUFF
  //|
  //======================================================
  int nprocessEffiEvents = 0;
  void processEff(o2::aod::JetMcCollision const& collision, soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JetCollisions>> const& recocolls, aod::JetParticles const& mcParticles)
  {
    if (cDebugLevel > 0) {
      ++nprocessEffiEvents;
      if (nprocessEffiEvents % 10000 == 0) {
        std::cout << "Processed MC (GEN) Events: " << nprocessEffiEvents << std::endl;
      }
    }

    if (fabs(collision.posZ()) > cfgEventVtxCut)
      return;
    if (recocolls.size() <= 0) { // not reconstructed
      return;
    }

    for (auto& recocoll : recocolls) { // poorly reconstructed
      auto goodEv = jetderiveddatautilities::selectCollision(recocoll, eventSelectionBits);
      if (goodEv) {
        goodEv = jetderiveddatautilities::selectTrigger(recocoll, RealTriggerMaskBits);
      }
      if (cfgJetMCHistos) {
        histos.fill(HIST("nEvents_MC_True"), 0.5);
      }
      if (!goodEv)
        return;
    }

    for (auto& particle : mcParticles) {
      if (particle.pdgCode() != 313)
        continue;
      if (std::fabs(particle.eta()) > cfgTrackMaxEta)
        continue;
      if (particle.pt() < cfgTrackMinPt)
        continue;

      /* // Not Yet
      if (cfg_Force_BR) {
        bool baddecay = false;
        for (auto& phidaughter : particle.daughters_as<aod::McParticles>()) {
          if (std::fabs(phidaughter.pdgCode()) != 321) {
            baddecay = true;
            break;
          }
          if (cfg_Force_Kaon_Acceptence) {
            if (std::fabs(phidaughter.eta()) > cfg_Track_MaxEta) {
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
        histos.fill(HIST("hEffGen_pT"), particle.pt());
      }
    } // loop over particles

  } // end of process
  PROCESS_SWITCH(kstarInOO, processEff, "Process Particles", false);

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

    auto [goodEv, code] = eventSelection(collision, true);
    if (cfgDataHistos) {
      histos.fill(HIST("nEvents"), 0.5);
    }
    if (!goodEv)
      return;

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (!trackSelection(track, true))
        continue;
      if (std::fabs(track.eta()) < cfgTrackMaxEta) {
        INELgt0 = true;
      }
    }
    if (!INELgt0)
      return;

    if (cfgDataHistos) {
      histos.fill(HIST("nEvents"), 1.5);
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
      if (std::fabs(collision1.posZ() - collision2.posZ()) <= cfgVtxMixCut) // set default to maybe 10
        VtxMixFlag = true;
      if (std::fabs(collision1.centFT0C() - collision2.centFT0C()) <= cfgVtxMixCut) // set default to maybe 10
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
    auto [goodEv, code] = eventSelection(collision, true);
    if (cfgMcHistos) {
      histos.fill(HIST("nEvents_MC"), 0.5);
    }
    if (!goodEv)
      return;

    bool INELgt0 = false;
    for (const auto& track : tracks) {
      if (!trackSelection(track, true))
        continue;
      if (std::fabs(track.eta()) < cfgTrackMaxEta) {
        INELgt0 = true;
      }
    }
    if (!INELgt0)
      return;

    if (cfgMcHistos) {
      histos.fill(HIST("nEvents_MC"), 1.5);
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

  //=======================================================
  //|
  //|             GENERATED MC STUFF (TRUE)
  //|
  //=======================================================
  int nEventsTrue = 0;
  void processMCTrue(EventCandidatesTrue::iterator const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, EventCandidates>> const& recocolls, aod::McParticles const& particles)
  {
    if (cDebugLevel > 0) {
      ++nEventsTrue;
    }

    if (fabs(collision.posZ()) > cfgEventVtxCut)
      return;
    if (recocolls.size() <= 0) { // not reconstructed
      if (cfgForceGenReco) {
        return;
      }
    }

    double centrality = -1;
    for (auto& recocoll : recocolls) {
      centrality = recocoll.centFT0C();
      auto [goodEv, code] = eventSelection(recocoll, true);

      if (cfgMcHistos) {
        histos.fill(HIST("nEvents_MC_True"), 0.5);
      }
      if (!goodEv)
        continue;
    } // for

    for (auto& particle : particles) {
      if (particle.pdgCode() != 313)
        continue; // Not K*0
      if (std::fabs(particle.eta()) > cfgTrackMaxEta)
        continue;

      if (cfgMcHistos) {
        histos.fill(HIST("hMC_kstar_True"), centrality, particle.pt());
      }
      if (cfgMcHistos) {
        histos.fill(HIST("nEvents_MC_True"), 1.5);
      }

    } // loop over particles
  } // processMCTrue
  PROCESS_SWITCH(kstarInOO, processMCTrue, "process MC True", false);

  void processEventsDummy(EventCandidates::iterator const&, TrackCandidates const&)
  {
    return;
  }
  PROCESS_SWITCH(kstarInOO, processEventsDummy, "dummy", false);
}; // kstarInOO
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarInOO>(cfgc)};
};
