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
/// \file kstarInJets.cxx
/// \brief the pT spectra of k*0(892) resonance analysis in Jets collisions
/// \author Jimun Lee <jimun.lee@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TMath.h>
#include <TVector2.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include <stdlib.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace ROOT::Math;

struct kstarInJets {
  SliceCache cache;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //==================================
  //|
  //|         Selection
  //|
  //==================================
  Configurable<int> cDebugLevel{"cDebugLevel", 0, "Resolution of Debug"};

  // Event Selection
  Configurable<std::string> cfgEventSelections{"cfgEventSelections", "sel8", "Set event selection"};
  Configurable<float> cfgEventVtxCut{"cfgEventVtxCut", 10.0, "V_z cut selection"};
  Configurable<double> cfgEventMaxEta{"cfgEventMaxEta", 1.0, "set INEL event eta cut"};
  ConfigurableAxis cfgCentAxis{"cfgCentAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0}, "Binning of the centrality axis"};

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
  // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNTCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  // TPC
  Configurable<double> cfgTrackFindableTPCClusters{"cfgTrackFindableTPCClusters", 0, "nFindable TPC Clusters"}; // origin = 50
  Configurable<double> cfgTrackTPCCrossedRows{"cfgTrackTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgTrackRowsOverFindable{"cfgTrackRowsOverFindable", 0, "nRowsOverFindable TPC CLusters"}; // origin = 1.2
  Configurable<double> cfgTrackTPCChi2{"cfgTrackTPCChi2", 4.0, "nTPC Chi2 per Cluster"};

  // ITS
  Configurable<double> cfgTrackITSChi2{"cfgTrackITSChi2", 36.0, "nITS Chi2 per Cluster"};

  // PID
  Configurable<bool> cfgTrackTPCPID{"cfgTrackTPCPID", true, "Enables TPC PID"};
  Configurable<bool> cfgTrackTOFPID{"cfgTrackTOFPID", true, "Enables TOF PID"};
  Configurable<bool> cfgTrackSquarePIDCut{"cfgTrackSqurePIDCut", true, "Enables PID cut shape square switch"};
  Configurable<bool> cfgTrackCirclePIDCut{"cfgTrackCirclePIDCut", false, "Enables PID cut shape circle switch"};
  Configurable<int> cfgTrackCircleValue{"cfgTrackCircleValue", 2, "Enables TOF TPC PID circle cut value"};
  Configurable<bool> cfgTrackTOFHard{"cfgTrackTOFHard", false, "Enables TOF Hard"};
  Configurable<float> cfgTrackTPCPIDnSig{"cfgTrackTPCPIDnSig", 4.0, "nTPC PID sigma"};
  Configurable<float> cfgTrackTOFPIDnSig{"cfgTrackTOFPIDnSig", 4.0, "nTOF PID sigma"};

  // MCGen
  Configurable<bool> cfgForceGenReco{"cfgForceGenReco", false, "Only consider events which are reconstructed (neglect event-loss)"};
  Configurable<bool> cfgReqMcEffPID{"cfgReqMcEffPID", true, "Request McEfficiency PID"};

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
  Configurable<bool> cfgJetInvHistos{"cfgJetInvHistos", false, "Enable Jet Inv Histos"};
  Configurable<bool> cfgJetMCHistos{"cfgJetMCHistos", false, "Enable Jet MC Histos"};
  Configurable<bool> cfgJetdRHistos{"cfgJetdRHistos", false, "Enable Jet dR QA Histos"};
  Configurable<bool> cfgMatchedGenHistos{"cfgMatchedGenHistos", false, "Enable Matched Gen Histos"};
  Configurable<bool> cfgPriorHistos{"cfgPriorHistos", false, "Enable Prior Gen and Rec Histos"};
  Configurable<bool> cfgFullClosureTest{"cfgFullClosureTest", false, "Enable Prior Gen for full closure"};

  Configurable<bool> cfgCutonTrig{"cfgCutonTrig", false, "Enable Jet Cut on Trig"};

  //======================
  //|
  //|        JET
  //|
  //======================
  Configurable<float> cfgJetpT{"cfgJetpT", 8.0, "Set Jet pT minimum"};
  Configurable<float> cfgJetR{"cfgJetR", 0.4, "Anti-kT Radius"};
  Configurable<float> cfgJetdR{"cfgJetdR", 0.4, "Set Jet radius parameter"};
  Configurable<float> cfgJetMaxEta{"cfgJetMaxEta", 0.9, "Set Jet Max Eta"};

  Configurable<bool> cfgReqJets{"cfgReqJets", false, "False: MB, True: In Jets"};
  Configurable<std::string> cfgRealTriggerMasks{"cfgRealTriggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<std::string> cfgTriggerMasksTest1{"cfgTriggerMasksTest1", "", "possible JE Trigger masks Test1"};
  Configurable<std::string> cfgTriggerMasksTest2{"cfgTriggerMasksTest2", "", "possible JE Trigger masks Test2"};
  Configurable<std::string> cfgTriggerMasksTest3{"cfgTriggerMasksTest3", "", "possible JE Trigger masks Test3"};

  std::vector<int> eventSelectionBits;
  std::vector<int> RealTriggerMaskBits;
  std::vector<int> triggerMaskBitsTest1;
  std::vector<int> triggerMaskBitsTest2;
  std::vector<int> triggerMaskBitsTest3;

  void init(o2::framework::InitContext&)
  {
    // HISTOGRAMS
    const AxisSpec axisPt{100, 0, 100, "#it{p}_{T}"}; // Jet pT
    const AxisSpec ptAxis = {1000, 0, 100};           // Particle pT
    const AxisSpec pidAxis = {120, -6, 6};
    const AxisSpec minvAxis = {cfgMinvNBins, cfgMinvMin, cfgMinvMax};
    const AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    const AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};

    if (cfgEventCutQA) {
      histos.add("hEvent_Cut", "Number of event after cuts", kTH1D, {{13, -0.5, 12.5}});
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

    if (cfgJetQAHistos) {
      histos.add("nTriggerQA", "nTriggerQA", kTH1F, {{8, 0.0, 8.0}});

      histos.add("JetpT", "Jet pT (GeV/c)", kTH1F, {{4000, 0., 200.}});
      histos.add("JetEta", "Jet Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("JetPhi", "Jet Phi", kTH1F, {{80, -1.0, 7.0}});

      histos.add("nGoodJets", "The number of good jets", kTH1F, {{6, -0.5, 5.5}});
      histos.add("nJetsPerEvent", "The number of jet per event", kTH1F, {{6, -0.5, 5.5}});
      histos.add("jetTrackEta", "Jet Track Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("jetTrackPhi", "Jet Track Phi", kTH1F, {{80, -1.0, 7.0}});
    }
    if (cfgJetdRHistos) {
      histos.add("mcpjet_pt", "MC particle-level Jet pT (GeV/c)", kTH1F, {{2000, 0., 100.}});
      histos.add("mcpjet_eta", "MC particle-level Jet Eta", kTH1F, {{100, -1.0, 1.0}});
      histos.add("mcpjet_phi", "MC particle-level Jet Phi", kTH1F, {{80, -1.0, 7.0}});

      histos.add("Gen_Kstar_BR", "Generated K*0 BR", kTH1F, {minvAxis});
      histos.add("Gen_KstarPt_BR", "Generated K*0 pT (GeV/c)", kTH1F, {{2000, 0., 100.}});

      histos.add("dR_kstar_pion", "dR between kstar and pion", {HistType::kTH2F, {{100, 0., 3.0}, {150, 0., 30.}}});
      histos.add("dR_kstar_kaon", "dR between kstar and kaon", {HistType::kTH2F, {{100, 0., 3.0}, {150, 0., 30.}}});
      histos.add("dR_taggedjet_kaon", "dR between tagged jet and kaon wo pT cut", {HistType::kTH2F, {{100, 0., 3.0}, {150, 0., 30.}}});
      histos.add("dR_taggedjet_pion", "dR between tagged jet and pion wo pT cut", {HistType::kTH2F, {{100, 0., 3.0}, {150, 0., 30.}}});
      histos.add("dR_taggedjet_all", "dR between tagged jet and kpi", {HistType::kTH2F, {{100, 0., 3.0}, {150, 0., 30.}}});
      histos.add("dR_taggedjet_all_6_8", "dR between tagged jet and kpi 6 < jet pT < 8", {HistType::kTH2F, {{100, 0., 3.0}, {150, 0., 30.}}});

      histos.add("Kstar_pT_INJet_6_8", "K*0 pT Inside 6-8 Jet", kTH1F, {{500, 0., 50.}});
      histos.add("Kstar_pT_INJet_8_10", "K*0 pT Inside 8-10 Jet", kTH1F, {{500, 0., 50.}});
      histos.add("Kstar_pT_INJet_10_12", "K*0 pT Inside 10-12 Jet", kTH1F, {{500, 0., 50.}});

      histos.add("recoveredJetpT_6_8to8", "Recovered Jet pT 6-8 --> 8", kTH1F, {{2000, 0., 100.}});
      histos.add("Kstar_pT_6_8to8_IN_recoveredJet", "K*0 pT In recovered Jet pT 6-8 --> 8", kTH1F, {{1000, 0., 50.}});
      histos.add("Kstar_pT_6_8to8_10_IN_recoveredJet", "K*0 pT In recovered Jet pT 6-8 --> 8-10", kTH1F, {{1000, 0., 50.}});

      histos.add("JetMigration", "bin to bin migration", {HistType::kTH2F, {{150, 0.0, 15.0, "True jet pT (GeV/c)"}, {150, 0., 15., "Recovered jet pT (GeV/c)"}}});
    }

    ///////////////////////////////////
    histos.add("nEvents", "nEvents", kTH1F, {{7, 0.0, 7.0}});
    histos.add("hUSS_KPi", "hUSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hLSS_KPi", "hLSS_KPi", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hUSS_PiK", "hUSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});
    histos.add("hLSS_PiK", "hLSS_PiK", kTHnSparseF, {cfgCentAxis, ptAxis, minvAxis});

    if (cfgJetInvHistos) {
      histos.add("hUSS_KPi_INSIDE", "hUSS_KPi_INSIDE", kTHnSparseF, {cfgCentAxis, axisPt, ptAxis, minvAxis});
      histos.add("hUSS_PiK_INSIDE", "hUSS_PiK_INSIDE", kTHnSparseF, {cfgCentAxis, axisPt, ptAxis, minvAxis});
      histos.add("hLSS_KPi_INSIDE", "hLSS_KPi_INSIDE", kTHnSparseF, {cfgCentAxis, axisPt, ptAxis, minvAxis});
      histos.add("hLSS_PiK_INSIDE", "hLSS_PiK_INSIDE", kTHnSparseF, {cfgCentAxis, axisPt, ptAxis, minvAxis});
    }
    if (cfgJetMCHistos) {
      histos.add("nEvents_Gen", "nEvents_Gen", kTH1F, {{7, 0.0, 7.0}});
      histos.add("hKstar_VtxEvt", "Kstar pT AC Vtx only (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hKstar_Gen", "Kstar pT AC good Event (GeV/c)", kTH1F, {{800, 0., 40.}});
      histos.add("hKstar_EffGen", "Gen. Kstar pT (GeV/c)", kTH1F, {{800, 0., 40.}});

      histos.add("Gen_daughter", "Gen_daughter", kTH1F, {{1000, -500., 500}});
      histos.add("Gen_mother", "Gen_mother", kTH1F, {{1000, -500., 500}});

      histos.add("nEvents_Rec", "nEvents_Rec", kTH1F, {{7, 0.0, 7.0}});
      histos.add("hEffRec_pT", "EffRec_pT (GeV/c)", kTH1F, {{1600, 0., 80.}});

      histos.add("lResonance_Rec_Inv", "lResonance_Rec_Inv", kTH1F, {minvAxis});
    }
    if (cfgMatchedGenHistos) {
      histos.add("nEvents_MCGen_MATCHED", "nEvents_MCGen_MATCHED", kTH1F, {{7, 0.0, 7.0}});
      histos.add("matchedJets_pT", "#it{p}_{T,jet part}; #it{p}_{T,jet det}", {HistType::kTH2F, {{100, 0., 20.}, {400, -20., 20.0}}});
      histos.add("matchedJets_eta", "#eta_{jet part}; #eta_{jet det}", {HistType::kTH2F, {{100, -1.0, 1.0}, {400, -20.0, 20.0}}});
      histos.add("matchedJets_phi", "#phi_{jet part}; #phi_{jet det}", {HistType::kTH2F, {{80, -1.0, 7.}, {400, -20.0, 20.}}});

      histos.add("GeneratedJets_pT", "GeneratedJets_pT", kTH1F, {axisPt});
      histos.add("NoMatchedJets_pT", "Gen. No Matched detector-level Jet", kTH1F, {axisPt});
      histos.add("MatchedJets_pT", "Gen. Matched detector-level Jet", kTH1F, {axisPt});

      histos.add("RespGen_Matched_rand0", "Gen. Response Matrix Matched 0", kTHnSparseD, {ptAxis, axisPt, ptAxis, axisPt});
      histos.add("RespGen_Matched_rand1", "Gen. Response Matrix Matched 1", kTHnSparseD, {ptAxis, axisPt, ptAxis, axisPt});

      histos.add("NoMatched", "NoMatched", kTHnSparseD, {ptAxis, axisPt});
      histos.add("Matched", "Matched", kTHnSparseD, {ptAxis, axisPt});
    }

    if (cfgPriorHistos) {
      histos.add("nEvents_Prior", "nEvents_Prior", kTH1F, {{7, 0.0, 7.0}});
      histos.add("Prior", "Prior", kTHnSparseD, {ptAxis, axisPt});

      histos.add("nEvents_PriorRec", "nEvents_PriorRec", kTH1F, {{7, 0.0, 7.0}});
      histos.add("Prior_Rec", "Prior_Rec", kTHnSparseD, {ptAxis, axisPt});
    }

    // Jet
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(cfgEventSelections));
    RealTriggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(cfgRealTriggerMasks);
    triggerMaskBitsTest1 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest1);
    triggerMaskBitsTest2 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest2);
    triggerMaskBitsTest3 = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasksTest3);
  } // end of init

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection,
                                      aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>;
  using JetTrackCandidatesMC = soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>;

  Preslice<aod::JTracks> tracksPerRecoCollision = aod::jtrack::collisionId;
  Filter JEPosZFilter = nabs(aod::jcollision::posZ) < cfgEventVtxCut;
  Filter JEMCPosZFilter = nabs(aod::jmccollision::posZ) < cfgEventVtxCut;
  Filter JetCuts = aod::jet::pt > cfgJetpT&& aod::jet::r == nround(cfgJetR.node() * 100.0f);

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  static constexpr int Kstar0PDG = 313;
  static constexpr int KaonPDG = 321;
  static constexpr int PionPDG = 211;

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
    if constexpr (requires { obj.tpcCrossedRowsOverFindableCls(); }) {
      if (objecttype == 2) {
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
    if (objecttype == 3) {
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
    if (objecttype == 4) {
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

  enum class objectType { Inclusive,
                          Jets };
  template <typename TrackType>
  void fillMinv(objectType type, const TrackType& trk1, const TrackType& trk2, const ROOT::Math::PxPyPzMVector& lReso, double centrality, double jetpt, bool flip)
  {
    double conjugate = trk1.sign() * trk2.sign();
    switch (type) {
      case objectType::Inclusive:
        if (conjugate < 0) {
          if (!flip) {
            histos.fill(HIST("hUSS_KPi"), centrality, lReso.Pt(), lReso.M());
          } else {
            histos.fill(HIST("hUSS_PiK"), centrality, lReso.Pt(), lReso.M());
          }
        } else if (conjugate > 0) {
          if (!flip) {
            histos.fill(HIST("hLSS_KPi"), centrality, lReso.Pt(), lReso.M());
          } else {
            histos.fill(HIST("hLSS_PiK"), centrality, lReso.Pt(), lReso.M());
          }
        }
        break;

      case objectType::Jets:
        if (cfgJetInvHistos) {
          if (conjugate < 0) {
            if (!flip) {
              histos.fill(HIST("hUSS_KPi_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
            } else {
              histos.fill(HIST("hUSS_PiK_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
            }
          } else if (conjugate > 0) {
            if (!flip) {
              histos.fill(HIST("hLSS_KPi_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
            } else {
              histos.fill(HIST("hLSS_PiK_INSIDE"), centrality, jetpt, lReso.Pt(), lReso.M());
            }
          }
        }
        break;

    } // switch
  } // fillMinv
  //======================================================================

  template <typename Table>
  bool hasINELgt0(Table const& objects)
  {
    for (const auto& obj : objects) {
      if (std::abs(obj.eta()) < cfgEventMaxEta && obj.pt() > 0.0) {
        return true;
      }
    }
    return false;
  };

  template <typename EventType>
  std::pair<bool, int> JetEventSelection(const EventType& event, const bool QA)
  {
    if (cfgEventCutQA && QA) {
      fillQA(false, event, 1);
      histos.fill(HIST("hEvent_Cut"), 0);
    }

    if (!jetderiveddatautilities::selectCollision(event, eventSelectionBits)) { // sel8
      return {false, 1};
    }
    if (cfgEventCutQA && QA) {
      histos.fill(HIST("hEvent_Cut"), 1);
    }
    if (std::abs(event.posZ()) > cfgEventVtxCut) {
      return {false, 2};
    }
    if (cfgEventCutQA && QA) {
      fillQA(true, event, 1);
      histos.fill(HIST("hEvent_Cut"), 2);
    }
    return {true, 8};
  };

  template <typename TracksType>
  bool trackSelection(const TracksType& track, const bool QA)
  {
    if (cfgTrackCutQA && QA) {
      fillQA(false, track, 2);
    }

    if (cfgTrackGlobalSel && !track.isGlobalTrack()) {
      return false;
    }
    if (track.pt() < cfgTrackMinPt) {
      return false;
    }
    if (std::abs(track.eta()) > cfgTrackMaxEta) {
      return false;
    }
    if (std::abs(track.dcaXY()) > cfgTrackMaxDCArToPVcut) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cfgTrackMaxDCAzToPVcut) {
      return false;
    }
    if (cfgTrackPrimaryTrack && !track.isPrimaryTrack()) {
      return false;
    }
    if (cfgTrackGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (cfgTrackFindableTPCClusters > 0 && track.tpcNClsFindable() < cfgTrackFindableTPCClusters) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < cfgTrackTPCCrossedRows) {
      return false;
    }
    if (cfgTrackRowsOverFindable > 0 && track.tpcCrossedRowsOverFindableCls() > cfgTrackRowsOverFindable) {
      return false;
    }
    if (track.tpcChi2NCl() > cfgTrackTPCChi2) {
      return false;
    }
    if (track.itsChi2NCl() > cfgTrackITSChi2) {
      return false;
    }
    if (cfgTrackConnectedToPV && !track.isPVContributor()) {
      return false;
    }

    if (cfgTrackCutQA && QA) {
      fillQA(true, track, 2);
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
      fillQA(false, candidate, 3);
    }

    // TPC
    if (cfgTrackSquarePIDCut) {
      if (std::abs(candidate.tpcNSigmaKa()) < cfgTrackTPCPIDnSig) {
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
      } // TPCPID cut
    } // end of square cut

    if (cfgTrackCirclePIDCut) {
      if (std::abs(candidate.tpcNSigmaKa()) < cfgTrackTPCPIDnSig) {
        tpcpid = std::abs(candidate.tpcNSigmaKa());
      }
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
      if (cfgTrackCutQA && QA)
        fillQA(true, candidate, 3);
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
      fillQA(false, candidate, 4);
    }

    // TPC
    if (cfgTrackSquarePIDCut) {
      if (std::abs(candidate.tpcNSigmaPi()) < cfgTrackTPCPIDnSig) {
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
      } // TPCPID cut
    } // end of square cut

    if (cfgTrackCirclePIDCut) {
      if (std::abs(candidate.tpcNSigmaPi()) < cfgTrackTPCPIDnSig) {
        tpcpid = std::abs(candidate.tpcNSigmaPi());
      }
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
      if (cfgTrackCutQA && QA)
        fillQA(true, candidate, 4);
      return true;
    }
    return false;
  }

  double getdR(double Eta1, double Phi1, double Eta2, double Phi2)
  {
    const double dphi = TVector2::Phi_mpi_pi(Phi1 - Phi2);
    const double deta = Eta1 - Eta2;
    return TMath::Sqrt(dphi * dphi + deta * deta);
  }

  template <typename JetType>
  double DistinguishJets(const JetType& jets, ROOT::Math::PxPyPzMVector& lResonance)
  {
    if (cDebugLevel > 0) {
      LOG(info) << "Found multiple jets to the same K*0.";
    }

    double bestR = 0;
    double bestJetpT = 0;
    for (auto const& jet : jets) {

      const double dR = getdR(jet.eta(), jet.phi(), lResonance.Eta(), lResonance.Phi());
      if (dR < cfgJetR && (bestR == 0 || dR < bestR)) {
        bestR = dR;
        bestJetpT = jet.pt();
      }
    } // jet
    return bestJetpT;
  }

  template <typename TracksType>
  ROOT::Math::PxPyPzMVector minvReconstruction(const TracksType& trk1, const TracksType& trk2, const bool QA, const bool flip)
  {
    if (!trackSelection(trk1, false) || !trackSelection(trk2, false)) {
      return {};
    }
    if (!trackPIDKaon(trk1, QA) || !trackPIDPion(trk2, QA)) {
      return {};
    }
    if (trk1.globalIndex() == trk2.globalIndex()) {
      return {};
    }

    ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonance;
    if (!flip) {
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
    } else {
      lDecayDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massKa);
    }
    lResonance = lDecayDaughter1 + lDecayDaughter2;

    if (std::abs(lResonance.Eta()) > cfgTrackMaxEta) {
      return {};
    }

    return {lResonance};
  } // minvReconstruction

  template <typename ParticleType>
  bool KPidaughters(ParticleType const& mcParticle)
  {
    if (std::abs(mcParticle.eta()) > cfgTrackMaxEta) {
      return {};
    }
    if (std::abs(mcParticle.pdgCode()) != Kstar0PDG) {
      return {};
    }
    if (!mcParticle.has_daughters()) {
      return {};
    }

    int nKaon = 0;
    int nPion = 0;

    for (const auto& daughter : mcParticle.template daughters_as<aod::JetParticles>()) {
      if (std::abs(daughter.eta()) > cfgTrackMaxEta) {
        continue;
      }
      if (std::abs(daughter.pdgCode()) == KaonPDG) {
        ++nKaon;
      } else if (std::abs(daughter.pdgCode()) == PionPDG) {
        ++nPion;
      }
    } // daughter loop
    return nKaon == 1 && nPion == 1;
  };

  using RecoCandidate = std::pair<int, ROOT::Math::PxPyPzMVector>;
  template <typename TracksType>
  std::vector<RecoCandidate> recoResonance(TracksType const& tracks)
  {
    std::vector<RecoCandidate> results;
    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      auto trk1 = track1.template track_as<TrackCandidatesMC>();
      auto trk2 = track2.template track_as<TrackCandidatesMC>();

      if (!trk1.has_mcParticle() || !trk2.has_mcParticle()) {
        continue;
      }
      if (!trackSelection(trk1, true) || !trackSelection(trk2, false)) {
        continue;
      }
      if (!trackPIDKaon(trk1, true) || !trackPIDPion(trk2, true)) {
        continue;
      }
      if (trk1.globalIndex() == trk2.globalIndex()) {
        continue;
      }

      auto particle1 = trk1.mcParticle();
      auto particle2 = trk2.mcParticle();

      if (!particle1.has_mothers() || !particle2.has_mothers()) {
        continue;
      }
      if (std::abs(particle1.pdgCode()) != KaonPDG) {
        continue;
      }
      if (std::abs(particle2.pdgCode()) != PionPDG) {
        continue;
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

      if (mothers1[0] != mothers2[0]) {
        continue; // Kaon and pion not from the same K*0
      }
      if (std::abs(mothers1PDG[0]) != Kstar0PDG) {
        continue;
      }
      if (std::abs(mothers2PDG[0]) != Kstar0PDG) {
        continue;
      }

      ROOT::Math::PxPyPzMVector lDaughter1, lDaughter2, lResonance;
      lDaughter1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDaughter2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
      lResonance = lDaughter1 + lDaughter2;

      if (std::abs(lResonance.Eta()) > cfgTrackMaxEta) {
        continue;
      }
      results.emplace_back(mothers1[0], lResonance);

    } // tracks loop
    return results;
  };

  template <typename TracksType, typename JetType>
  void JetTrackSlicing(aod::JetCollision const& collision, TracksType const& jetTracks, const JetType& chargedjets, const bool QA)
  {
    //=============
    //| Inclusive
    //=============
    auto centrality = collision.centFT0C();
    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(jetTracks, jetTracks))) {
      auto trk1 = track1.template track_as<TrackCandidates>();
      auto trk2 = track2.template track_as<TrackCandidates>();

      for (bool flip : {false, true}) {
        auto lResonance = minvReconstruction(trk1, trk2, QA, flip);
        if (lResonance.M() <= 0) {
          continue;
        }
        fillMinv(objectType::Inclusive, trk1, trk2, lResonance, centrality, -1.0, flip);

        //========
        //| Jets
        //========
        bool jetFlag = false;
        int goodjets = 0;
        double jetpt = 0;
        double dR = 0;
        for (auto const& jet : chargedjets) {
          double dphi = TVector2::Phi_mpi_pi(jet.phi() - lResonance.Phi());
          double deta = jet.eta() - lResonance.Eta();
          dR = TMath::Sqrt((deta * deta) + (dphi * dphi));

          if (dR < cfgJetR) {
            jetFlag = true;
            jetpt = jet.pt();
            goodjets++;
          }
        } // The loop of chargedjets
        if (cfgJetQAHistos) {
          histos.fill(HIST("nGoodJets"), goodjets);
        }

        if (goodjets > 1) {
          jetpt = DistinguishJets<JetType>(chargedjets, lResonance);
        }

        if (jetFlag) {
          fillMinv(objectType::Jets, trk1, trk2, lResonance, centrality, jetpt, flip);
        } // jetFlag
      } // flip
    } // jetTracks
  } // JetTrackSlicing

  template <typename TracksType, typename JetType>
  void JetTrackSlicingMC(aod::JetCollision const& collision, TracksType const& jetTracks, const JetType& mcdjets, const bool QA)
  {
    //============================
    //| MB: Track Reconstruction
    //============================
    auto centrality = collision.centFT0C();
    for (const auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(jetTracks, jetTracks))) {
      auto trk1 = track1.template track_as<TrackCandidatesMC>();
      auto trk2 = track2.template track_as<TrackCandidatesMC>();

      if (!trk1.has_mcParticle() || !trk2.has_mcParticle()) {
        continue;
      }

      for (bool flip : {false, true}) {
        auto lResonance = minvReconstruction(trk1, trk2, QA, flip);
        if (lResonance.M() <= 0) {
          continue;
        }
        fillMinv(objectType::Inclusive, trk1, trk2, lResonance, centrality, -1.0, flip);

        //==============================
        //| Jets: Track Reconstruction
        //==============================
        bool jetFlag = false;
        int goodjets = 0;
        double jetpt = 0;
        for (auto const& jet : mcdjets) {
          const double dR = getdR(jet.eta(), jet.phi(), lResonance.Eta(), lResonance.Phi());

          if (dR < cfgJetR) {
            jetFlag = true;
            jetpt = jet.pt();
            goodjets++;
          }
        } // mcdjets
        if (cfgJetQAHistos) {
          histos.fill(HIST("nGoodJets"), goodjets);
        }

        if (goodjets > 1) {
          jetpt = DistinguishJets<JetType>(mcdjets, lResonance);
        }

        if (jetFlag) {
          fillMinv(objectType::Jets, trk1, trk2, lResonance, centrality, jetpt, flip);
        }
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

  void processDataJets(o2::aod::JetCollision const& collision, JetFilteredJets const& chargedjets, JetTrackCandidates const& jetTracks, TrackCandidates const&)
  {
    static int nJetEvents = 0;
    if (cDebugLevel > 0) {
      nJetEvents++;
      if ((nJetEvents + 1) % 10000 == 0) {
        LOG(info) << "Processed Jet Data Events: " << nJetEvents;
      }
    }
    histos.fill(HIST("nEvents"), 0.5); // Raw event

    if (!hasINELgt0(jetTracks)) {
      return;
    }
    histos.fill(HIST("nEvents"), 1.5); // INEL>0 event

    auto [goodEv, code] = JetEventSelection(collision, true);
    if (!goodEv) {
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
      }
      if (VTtest1 || VTtest2 || VTtest3) {
        return;
      }
    } // Trigger cut
    histos.fill(HIST("nEvents"), 3.5); // After trigger cut

    bool HasJets = false;
    int nJets = 0;
    for (auto& chargedjet : chargedjets) {
      if (std::abs(chargedjet.eta()) > (cfgJetMaxEta - cfgJetdR)) {
        continue;
      }

      nJets++;
      if (cfgJetQAHistos) {
        histos.fill(HIST("JetpT"), chargedjet.pt());
        histos.fill(HIST("JetEta"), chargedjet.eta());
        histos.fill(HIST("JetPhi"), chargedjet.phi());
      }
      if (chargedjet.pt() > cfgJetpT) {
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
    histos.fill(HIST("nEvents"), 4.5); // Has jets

    // Just for checking track QA plot
    for (auto& jetTrack : jetTracks) {
      auto originTrack = jetTrack.track_as<TrackCandidates>();
      if (!trackSelection(originTrack, true)) {
        continue;
      }

      if (cfgJetQAHistos) {
        histos.fill(HIST("jetTrackEta"), jetTrack.eta());
        histos.fill(HIST("jetTrackPhi"), jetTrack.phi());
      }
    } // jetTrack loop

    JetTrackSlicing(collision, jetTracks, chargedjets, true);

  } // ProcessDataJets
  PROCESS_SWITCH(kstarInJets, processDataJets, "process Data Jets", false);

  //=======================================================
  //|
  //|                 JET MC STUFF
  //|
  //=======================================================

  void processMCJets(o2::aod::JetCollision const& collision,
                     JetTrackCandidatesMC const& jetTracks, TrackCandidatesMC const&,
                     soa::Filtered<aod::ChargedMCDetectorLevelJets> const& mcdjets,
                     aod::McParticles const&, aod::JetParticles const&)
  {
    static int nJetMCEvents = 0;
    if (cDebugLevel > 0) {
      nJetMCEvents++;
      if ((nJetMCEvents + 1) % 10000 == 0) {
        LOG(info) << "Processed Jet MC Events: " << nJetMCEvents;
      }
    }
    if (!cfgJetMCHistos) {
      return;
    }
    histos.fill(HIST("nEvents"), 0.5); // Gen event

    if (!hasINELgt0(jetTracks)) {
      return;
    }
    histos.fill(HIST("nEvents"), 1.5); // INEL>0 event

    auto [goodEv, code] = JetEventSelection(collision, true);
    if (!goodEv) {
      return;
    }
    histos.fill(HIST("nEvents"), 2.5); // Inclusive event

    // The trigger option doesn't work
    // if (!jetderiveddatautilities::selectTrigger(collision, RealTriggerMaskBits))
    //   return;
    // histos.fill(HIST("nEvents"), 3.5); // Events for Inclusive

    bool HasJets = false;
    int nJets = 0;
    for (auto& mcdjet : mcdjets) {
      if (std::abs(mcdjet.eta()) > cfgJetMaxEta - cfgJetdR) {
        continue;
      }
      nJets++;
      if (cfgJetQAHistos) {
        histos.fill(HIST("JetpT"), mcdjet.pt());
        histos.fill(HIST("JetEta"), mcdjet.eta());
        histos.fill(HIST("JetPhi"), mcdjet.phi());
      }
      if (mcdjet.pt() > cfgJetpT) {
        HasJets = true;
      }
    } // mcdjets
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
    JetTrackSlicingMC(collision, jetTracks, mcdjets, true);
  } // process loop
  PROCESS_SWITCH(kstarInJets, processMCJets, "process MC Jets", false);

  //========================
  //|
  //|    GENERATED STUFF
  //|
  //========================
  void processGen(o2::aod::JetMcCollision const& collision,
                  soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JetCollisions>> const& recocolls,
                  aod::JetParticles const& mcParticles, aod::McParticles const&)
  {
    static int nprocessGenEvents = 0;
    if (cDebugLevel > 0) {
      ++nprocessGenEvents;
      if (nprocessGenEvents % 10000 == 0) {
        LOG(info) << "Processed MC (GEN) Events: " << nprocessGenEvents;
      }
    }
    if (!cfgJetMCHistos) {
      return;
    }
    histos.fill(HIST("nEvents_Gen"), 0.5);

    //=======================
    //| Event & Signal loss
    //=======================
    if (!hasINELgt0(mcParticles)) {
      return;
    }
    histos.fill(HIST("nEvents_Gen"), 1.5); // EL: denominator

    if (std::abs(collision.posZ()) > cfgEventVtxCut) {
      return;
    }

    //  Find Gen K*0
    for (const auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.pdgCode()) != Kstar0PDG) {
        continue;
      }
      if (std::abs(mcParticle.eta()) > cfgTrackMaxEta) {
        continue;
      }
      histos.fill(HIST("hKstar_VtxEvt"), mcParticle.pt()); // SL: denominator
    } // mcParticles in GenEvents

    if (recocolls.size() <= 0) { // not reconstructed
      return;
    }

    for (auto& recocoll : recocolls) { // poorly reconstructed
      auto goodEv = jetderiveddatautilities::selectCollision(recocoll, eventSelectionBits);
      if (std::abs(recocoll.posZ()) > cfgEventVtxCut) {
        goodEv = false;
      }
      if (!goodEv) {
        return;
      }
    } // reco.coll
    histos.fill(HIST("nEvents_Gen"), 2.5); // EL: numerator

    //  Find Gen K*0 for SL
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.pdgCode()) != Kstar0PDG) {
        continue;
      }
      if (std::abs(mcParticle.eta()) > cfgTrackMaxEta) {
        continue;
      }
      if (mcParticle.pt() < cfgTrackMinPt) {
        continue;
      }
      histos.fill(HIST("hKstar_Gen"), mcParticle.pt()); // SL: numorator
      histos.fill(HIST("Gen_mother"), mcParticle.pdgCode());
    }

    //==============
    //| Efficiency
    //==============
    //  Find Gen K*0 for eff
    for (const auto& mcParticle : mcParticles) {
      if (!KPidaughters(mcParticle)) {
        continue;
      }
      histos.fill(HIST("hKstar_EffGen"), mcParticle.pt()); // eff: denominator
    }

  } // end of process
  PROCESS_SWITCH(kstarInJets, processGen, "Process Generated Particles", false);

  //========================
  //|
  //|    Rec. STUFF
  //|
  //========================
  void processRec(o2::aod::JetCollision const& collision,
                  JetTrackCandidatesMC const& tracks, TrackCandidatesMC const&,
                  aod::JetParticles const&, aod::McParticles const&)
  {
    static int nprocessRecEvents = 0;
    if (cDebugLevel > 0) {
      ++nprocessRecEvents;
      if (nprocessRecEvents % 10000 == 0) {
        LOG(info) << "Processed MC (REC) Events: " << nprocessRecEvents;
      }
    }
    if (!cfgJetMCHistos) {
      return;
    }
    histos.fill(HIST("nEvents_Rec"), 0.5);

    if (!hasINELgt0(tracks)) {
      return;
    }
    histos.fill(HIST("nEvents_Rec"), 1.5);

    auto [goodEv, code] = JetEventSelection(collision, true);
    if (!goodEv)
      return;
    histos.fill(HIST("nEvents_Rec"), 2.5);

    for (const auto& candidate : recoResonance(tracks)) {
      const auto& lResonance = candidate.second;

      histos.fill(HIST("hEffRec_pT"), lResonance.Pt());
      histos.fill(HIST("lResonance_Rec_Inv"), lResonance.M());

    } // recoResonances
  } // processRec
  PROCESS_SWITCH(kstarInJets, processRec, "Process Reconstructed Particles", false);

  void processJetBias(o2::aod::JetMcCollision const& collision,
                      soa::Filtered<aod::ChargedMCParticleLevelJets> const& mcpjets,
                      aod::JetParticles const& mcParticles, aod::McParticles const&)
  {
    static int ndRtest = 0;
    if (cDebugLevel > 0) {
      ++ndRtest;
      if (ndRtest % 10000 == 0) {
        LOG(info) << "Processed dR test: " << ndRtest;
      }
    }
    if (!cfgJetdRHistos) {
      return;
    }

    if (!hasINELgt0(mcpjets)) {
      return;
    }
    if (std::abs(collision.posZ()) > cfgEventVtxCut) {
      return;
    }
    histos.fill(HIST("nEvents"), 0.5);

    for (auto& mcParticle : mcParticles) {

      if (!KPidaughters(mcParticle)) {
        continue;
      }
      histos.fill(HIST("Gen_KstarPt_BR"), mcParticle.pt());

      //==================
      // Distinguish Jets
      double bestR = 999;
      ROOT::Math::PxPyPzEVector bestJet;
      ROOT::Math::PxPyPzEVector lResonance = ROOT::Math::PxPyPzEVector(mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e());

      for (auto& mcpjet : mcpjets) {
        if (mcpjet.pt() < cfgJetpT) {
          continue;
        }
        histos.fill(HIST("mcpjet_eta"), mcpjet.eta());
        histos.fill(HIST("mcpjet_phi"), mcpjet.phi());
        histos.fill(HIST("mcpjet_pt"), mcpjet.pt());

        const double dR = getdR(mcpjet.eta(), mcpjet.phi(), lResonance.Eta(), lResonance.Phi());
        if (dR < bestR) {
          bestR = dR;
          bestJet = ROOT::Math::PxPyPzEVector(mcpjet.px(), mcpjet.py(), mcpjet.pz(), mcpjet.energy());
        } // bestR
      } // mcpJets
      if (bestR > cfgJetdR) {
        continue;
      }

      //==================
      // daughters
      double dR_kaon = 0;
      double dR_pion = 0;
      double dR_k_kstar = 0;
      double dR_pi_kstar = 0;

      bool isKaonInJet = false;
      bool isPionInJet = false;
      ROOT::Math::PxPyPzEVector Daughter, missing_daughter;

      for (auto& daughter : mcParticle.daughters_as<aod::JetParticles>()) {
        Daughter = ROOT::Math::PxPyPzEVector(daughter.px(), daughter.py(), daughter.pz(), daughter.energy());

        if (std::abs(daughter.pdgCode()) == KaonPDG) {

          // K*0 and kaon
          dR_k_kstar = getdR(lResonance.Eta(), lResonance.Phi(), daughter.eta(), daughter.phi());
          dR_kaon = getdR(bestJet.Eta(), bestJet.Phi(), daughter.eta(), daughter.phi());

          if (dR_kaon <= cfgJetdR) {
            isKaonInJet = true;
          } else {
            missing_daughter += Daughter;
          }
        } else if (std::abs(daughter.pdgCode()) == PionPDG) {

          // K*0 and pion
          dR_pi_kstar = getdR(lResonance.Eta(), lResonance.Phi(), daughter.eta(), daughter.phi());
          dR_pion = getdR(bestJet.Eta(), bestJet.Phi(), daughter.eta(), daughter.phi());

          if (dR_pion <= cfgJetdR) {
            isPionInJet = true;
          } else {
            missing_daughter += Daughter;
          }
        }
      } // daughter

      histos.fill(HIST("dR_kstar_kaon"), dR_k_kstar, lResonance.Pt());
      histos.fill(HIST("dR_kstar_pion"), dR_pi_kstar, lResonance.Pt());

      if (bestR < cfgJetdR) {
        histos.fill(HIST("dR_taggedjet_all"), dR_kaon, lResonance.Pt());
        histos.fill(HIST("dR_taggedjet_all"), dR_pion, lResonance.Pt());
        histos.fill(HIST("dR_taggedjet_kaon"), dR_kaon, lResonance.Pt());
        histos.fill(HIST("dR_taggedjet_pion"), dR_pion, lResonance.Pt());

        if (bestJet.Pt() > 6.0 && bestJet.Pt() < 8.0) {
          histos.fill(HIST("dR_taggedjet_all_6_8"), dR_kaon, lResonance.Pt());
          histos.fill(HIST("dR_taggedjet_all_6_8"), dR_pion, lResonance.Pt());
          histos.fill(HIST("Kstar_pT_INJet_6_8"), lResonance.Pt());
        }

        if (bestJet.Pt() > 8.0 && bestJet.Pt() < 10.0) {
          histos.fill(HIST("Kstar_pT_INJet_8_10"), lResonance.Pt());
        }
        if (bestJet.Pt() > 10.0 && bestJet.Pt() < 12.0) {
          histos.fill(HIST("Kstar_pT_INJet_10_12"), lResonance.Pt());
        }

        // Daughters out of jets
        if (!isKaonInJet || !isPionInJet) {
          ROOT::Math::PxPyPzEVector recoveredJet = bestJet + missing_daughter;
          if (bestJet.Pt() > 6.0 && bestJet.Pt() < 8.0) {

            if (recoveredJet.Pt() > 8.0) {
              histos.fill(HIST("recoveredJetpT_6_8to8"), recoveredJet.Pt());
              histos.fill(HIST("Kstar_pT_6_8to8_IN_recoveredJet"), lResonance.Pt());
            }
            if (recoveredJet.Pt() > 8.0 && recoveredJet.Pt() < 10.0) {
              histos.fill(HIST("Kstar_pT_6_8to8_10_IN_recoveredJet"), lResonance.Pt());
            }
          }
          histos.fill(HIST("JetMigration"), bestJet.Pt(), recoveredJet.Pt());
        } // Daughters out of jet
      } // K*0 InJet
    } // mcParticles
  };
  PROCESS_SWITCH(kstarInJets, processJetBias, "Process dR of K*0 Inclusive and Inside jet", false);

  //========================
  //
  //   Response matrix
  //
  //========================

  using JetMCDTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>>;
  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;

  void processMatched(aod::JetMcCollision const& collision,
                      soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JetCollisions>> const& recocolls,
                      JetMCPTable const& mcpjets, JetMCDTable const&,
                      JetTrackCandidatesMC const& tracks, TrackCandidatesMC const&,
                      aod::JetParticles const& mcParticles, aod::McParticles const&)
  {
    static int nProcessMatched = 0;
    if (cDebugLevel > 0) {
      ++nProcessMatched;
      if (nProcessMatched % 10000 == 0) {
        LOG(info) << "Processed Matched: " << nProcessMatched;
      }
    }
    if (!cfgMatchedGenHistos) {
      return;
    }
    histos.fill(HIST("nEvents_MCGen_MATCHED"), 0.5);

    if (std::fabs(collision.posZ()) > cfgEventVtxCut) {
      return;
    }
    if (recocolls.size() <= 0) {
      return;
    }
    for (auto& recocoll : recocolls) { // poorly reconstructed
      auto goodEv = jetderiveddatautilities::selectCollision(recocoll, eventSelectionBits);
      if (std::abs(recocoll.posZ()) > cfgEventVtxCut) {
        goodEv = false;
      }
      if (!goodEv) {
        return;
      }
    } // reco.coll

    if (!hasINELgt0(mcpjets)) {
      return;
    }
    histos.fill(HIST("nEvents_MCGen_MATCHED"), 1.5);

    // 1. Gen jet --> Reco jet
    std::vector<double> mcd_pt{};
    std::vector<double> mcd_phi{};
    std::vector<double> mcd_eta{};
    std::vector<double> mcp_pt{};
    std::vector<double> mcp_phi{};
    std::vector<double> mcp_eta{};
    std::vector<double> NOMATCHED_mcp_pt{};
    std::vector<double> NOMATCHED_mcp_phi{};
    std::vector<double> NOMATCHED_mcp_eta{};

    bool hasJets = false;
    for (const auto& mcpjet : mcpjets) {
      histos.fill(HIST("GeneratedJets_pT"), mcpjet.pt());

      int nMatches = 0;
      if (mcpjet.has_matchedJetGeo()) {
        for (const auto& mcdjet : mcpjet.template matchedJetGeo_as<JetMCDTable>()) {
          if (!mcdjet.has_matchedJetGeo()) {
            continue;
          }
          ++nMatches;
          hasJets = true;

          histos.fill(HIST("matchedJets_pT"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
          histos.fill(HIST("matchedJets_phi"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
          histos.fill(HIST("matchedJets_eta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

          mcd_pt.push_back(mcdjet.pt());
          mcd_phi.push_back(mcdjet.phi());
          mcd_eta.push_back(mcdjet.eta());

          mcp_pt.push_back(mcpjet.pt());
          mcp_phi.push_back(mcpjet.phi());
          mcp_eta.push_back(mcpjet.eta());
        } // mcd matching
      } // mcp matching

      if (nMatches > 0) {
        histos.fill(HIST("MatchedJets_pT"), mcpjet.pt());
      } else {
        histos.fill(HIST("NoMatchedJets_pT"), mcpjet.pt());

        NOMATCHED_mcp_pt.push_back(mcpjet.pt());
        NOMATCHED_mcp_phi.push_back(mcpjet.phi());
        NOMATCHED_mcp_eta.push_back(mcpjet.eta());
      }
    } // mcpjets loop

    if (hasJets) {
      histos.fill(HIST("nEvents_MCGen_MATCHED"), 2.5);
    }

    // Filtered the meaningful tracks from recocoll

    std::vector<RecoCandidate> recoCandidates;
    for (const auto& recocoll : recocolls) {
      const auto collisionTracks = tracks.sliceBy(tracksPerRecoCollision, recocoll.globalIndex());
      const auto candidates = recoResonance(collisionTracks);

      recoCandidates.insert(recoCandidates.end(), candidates.begin(), candidates.end());
    }

    // Gen K*0 loop
    for (const auto& mcParticle : mcParticles) {
      if (!KPidaughters(mcParticle)) {
        continue;
      }

      double best_dR = cfgJetdR;
      //==============
      // Matched Jet
      //==============
      auto bestMatched_index = mcp_pt.size();
      for (std::vector<double>::size_type i = 0; i < mcp_pt.size(); ++i) {
        if (mcp_pt[i] < cfgJetpT) {
          continue;
        }
        const double dR = getdR(mcp_eta[i], mcp_phi[i], mcParticle.eta(), mcParticle.phi());
        if (dR < best_dR) {
          best_dR = dR;
          bestMatched_index = i;
        }
      } // mcpjet loop

      //===============
      // NoMatched Jet
      //===============
      auto bestNoMatched_index = NOMATCHED_mcp_pt.size();
      for (std::vector<double>::size_type i = 0; i < NOMATCHED_mcp_pt.size(); ++i) {
        if (NOMATCHED_mcp_pt[i] < cfgJetpT) {
          continue;
        }
        const double dR = getdR(NOMATCHED_mcp_eta[i], NOMATCHED_mcp_phi[i], mcParticle.eta(), mcParticle.phi());
        if (dR < best_dR) {
          best_dR = dR;
          bestNoMatched_index = i;
        }
      } // NOMATCHED mcpjet loop

      if (bestNoMatched_index != NOMATCHED_mcp_pt.size()) {
        histos.fill(HIST("NoMatched"), mcParticle.pt(), NOMATCHED_mcp_pt[bestNoMatched_index]);
        continue;
      }

      if (bestMatched_index == mcp_pt.size()) {
        continue;
      }

      const auto i_jet = bestMatched_index;
      const double jetpt_mcp = mcp_pt[i_jet];

      histos.fill(HIST("Matched"), mcParticle.pt(), jetpt_mcp);

      const auto gen_Index = mcParticle.globalIndex();
      for (const auto& [motherIndex, lResonance] : recoCandidates) {
        if (motherIndex != gen_Index) {
          continue;
        }

        // Matched Rec. K*0 --> Matching detector-level jet

        const double dR = getdR(mcd_eta[i_jet], mcd_phi[i_jet], lResonance.Eta(), lResonance.Phi());
        const bool recoJetFlag = dR < cfgJetdR;

        if (!recoJetFlag) {
          continue;
        }
        const double jetpt_mcd = mcd_pt[i_jet];

        unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
        int dice = rand_r(&seed) % 2;

        if (dice > 0) {
          // Train
          histos.fill(HIST("RespGen_Matched_rand0"), lResonance.Pt(), jetpt_mcd, mcParticle.pt(), jetpt_mcp);
        } else {
          // Test
          histos.fill(HIST("RespGen_Matched_rand1"), lResonance.Pt(), jetpt_mcd, mcParticle.pt(), jetpt_mcp);
        }
        break;
      } // recoCandidates loop
    } // mcParticles
  }; // processmatchedJet
  PROCESS_SWITCH(kstarInJets, processMatched, "Process Matched K*0 jets", false);

  void processPrior(aod::JetMcCollision const& collision,
                    soa::SmallGroups<soa::Join<aod::JMcCollisionLbs, aod::JetCollisions>> const& recocolls,
                    JetMCPTable const& mcpjets,
                    aod::JetParticles const& mcParticles, aod::McParticles const&)
  {
    static int nPrior = 0;
    if (cDebugLevel > 0) {
      ++nPrior;
      if (nPrior % 10000 == 0) {
        LOG(info) << "Processed Prior: " << nPrior;
      }
    }
    if (!cfgPriorHistos) {
      return;
    }
    histos.fill(HIST("nEvents_Prior"), 0.5);

    if (std::fabs(collision.posZ()) > cfgEventVtxCut) {
      return;
    }
    if (!hasINELgt0(mcParticles)) {
      return;
    }

    if (!cfgFullClosureTest) {
      if (recocolls.size() <= 0) {
        return;
      }
      for (auto& recocoll : recocolls) {
        if (std::abs(recocoll.posZ()) > cfgEventVtxCut) {
          return;
        }
        if (!jetderiveddatautilities::selectCollision(recocoll, eventSelectionBits)) {
          return;
        }
      }
    } // Full Closure Test
    histos.fill(HIST("nEvents_Prior"), 1.5);

    for (const auto& mcParticle : mcParticles) {
      if (!KPidaughters(mcParticle)) {
        continue;
      }

      double best_dR = cfgJetdR;
      ROOT::Math::PxPyPzEVector bestJet;

      for (const auto& mcpjet : mcpjets) {
        if (mcpjet.pt() < cfgJetpT) {
          continue;
        }

        double dR = getdR(mcpjet.eta(), mcpjet.phi(), mcParticle.eta(), mcParticle.phi());
        if (dR < best_dR) {
          best_dR = dR;
          bestJet = ROOT::Math::PxPyPzEVector(mcpjet.px(), mcpjet.py(), mcpjet.pz(), mcpjet.energy());
        } // bestJet
      } // mcpjets loop

      if (best_dR < cfgJetdR) {
        histos.fill(HIST("Prior"), mcParticle.pt(), bestJet.Pt());
      }
    } // mcParticles loop
  };
  PROCESS_SWITCH(kstarInJets, processPrior, "Process Prior K*0 jets", false);

  void processRecInJet(aod::JetCollision const& collision,
                       JetMCDTable const& mcdjets,
                       JetTrackCandidatesMC const& tracks, TrackCandidatesMC const&,
                       aod::JetParticles const& mcParticles, aod::McParticles const&)
  {
    static int nRec = 0;
    if (cDebugLevel > 0) {
      ++nRec;
      if (nRec % 10000 == 0) {
        LOG(info) << "Processed Rec In Jet: " << nRec;
      }
    }
    if (!cfgPriorHistos) {
      return;
    }
    histos.fill(HIST("nEvents_PriorRec"), 0.5);

    if (std::fabs(collision.posZ()) > cfgEventVtxCut) {
      return;
    }
    if (!hasINELgt0(mcParticles)) {
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    histos.fill(HIST("nEvents_PriorRec"), 1.5);

    for (const auto& candidate : recoResonance(tracks)) {
      const auto& lResonance = candidate.second;

      double best_dR = cfgJetdR;
      ROOT::Math::PxPyPzEVector bestJet;

      for (const auto& mcdjet : mcdjets) {
        if (mcdjet.pt() < cfgJetpT) {
          continue;
        }

        double dR = getdR(mcdjet.eta(), mcdjet.phi(), lResonance.eta(), lResonance.phi());
        if (dR < best_dR) {
          best_dR = dR;
          bestJet = ROOT::Math::PxPyPzEVector(mcdjet.px(), mcdjet.py(), mcdjet.pz(), mcdjet.energy());
        }
      } // mcdjets loop

      if (best_dR < cfgJetdR) {
        histos.fill(HIST("Prior_Rec"), lResonance.pt(), bestJet.Pt());
      }

    } // tracks loop
  };
  PROCESS_SWITCH(kstarInJets, processRecInJet, "Process Rec K*0 In Jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarInJets>(cfgc)};
}
