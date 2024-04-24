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

/// \file lambda1520analysis.cxx
/// \brief This task reconstructs track-track decay lambda(1520) resonance candidate
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

#include "TLorentzVector.h"
#include "TF1.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::constants::physics;

struct lambda1520analysis {
  // Define slice per Resocollision
  SliceCache cache;
  Preslice<aod::ResoTracks> perResoCollision = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  // switches
  Configurable<bool> cEtaAssym{"cEtaAssym", false, "Turn on/off EtaAssym calculation"};
  Configurable<bool> isFilladditionalQA{"isFilladditionalQA", false, "Turn on/off additional QA plots"};
  Configurable<bool> cAddlTrackcut{"cAddlTrackcut", false, "Switch to turn on/off Additional track cut"};
  Configurable<bool> cOldPIDcut{"cOldPIDcut", false, "Switch to turn on/off old PID cut to apply pt dependent cut"};
  Configurable<bool> cDCAr7SigCut{"cDCAr7SigCut", false, "Track DCAr 7 Sigma cut to PV Maximum"};
  Configurable<bool> cKinCuts{"cKinCuts", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<bool> cTPCNClsFound{"cTPCNClsFound", false, "Switch to turn on/off TPCNClsFound cut"};
  Configurable<bool> additionalEvsel{"additionalEvsel", true, "Additional event selcection"};

  // Pre-selection Track cuts
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15f, "Minimal pT for tracks"};
  Configurable<float> cMinRtpccut{"cMinRtpccut", 0.8f, "minimum ratio of number of Xrows to findable clusters in TPC"};
  Configurable<float> cMaxChi2ITScut{"cMaxChi2ITScut", 36.0f, "Maximal pT for Chi2/cluster for ITS"};
  Configurable<float> cMaxChi2TPCcut{"cMaxChi2TPCcut", 4.0f, "Maximal pT for Chi2/cluster for TPC"};
  Configurable<float> cMinTPCNClsFound{"cMinTPCNClsFound", 120, "minimum TPCNClsFound value for good track"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cMinTPCncr{"cMinTPCncr", 70, "Minimum number of TPC X rows"};

  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1f, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};

  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor

  /// PID Selections
  Configurable<float> cRejNsigmaTpc{"cRejNsigmaTpc", 3.0, "Reject tracks to improve purity of TPC PID"}; // Reject missidentified particles when tpc bands merge
  Configurable<float> cRejNsigmaTof{"cRejNsigmaTof", 3.0, "Reject tracks to improve purity of TOF PID"}; // Reject missidentified particles when tpc bands merge
  Configurable<bool> cUseRejNsigma{"cUseRejNsigma", false, "Switch on/off track rejection method to improve purity"};

  // Kaon
  // Old PID use case
  Configurable<std::vector<double>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
  Configurable<double> cMaxTPCnSigmaKaonVETO{"cMaxTPCnSigmaKaonVETO", 3.0, "TPC nSigma VETO cut for Kaon"}; // TPC

  // New PID use case
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};                // TPC
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"};   // Combined
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Pion"};                // TOF
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection

  // Proton
  // Old PID use case
  Configurable<std::vector<double>> protonTPCPIDpTintv{"protonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> protonTPCPIDcuts{"protonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> protonTOFPIDpTintv{"protonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> protonTOFPIDcuts{"protonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
  Configurable<double> cMaxTPCnSigmaProtonVETO{"cMaxTPCnSigmaProtonVETO", 3.0, "TPC nSigma VETO cut for Proton"}; // TPC

  // New PID use case
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 3.0, "TPC nSigma cut for Proton"};              // TPC
  Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", 3.0, "Combined nSigma cut for Proton"}; // Combined
  Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 3.0, "TOF nSigma cut for Pion"};                // TOF
  Configurable<bool> cUseOnlyTOFTrackPr{"cUseOnlyTOFTrackPr", false, "Use only TOF track for PID selection"};     // Use only TOF track for Pion PID selection

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - multiplicity"};

  /// Figures
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsEta{"binsEta", {100, -1, 1}, ""};
  ConfigurableAxis binsMass{"binsMass", {1700, 1.3, 3.0}, "Invariant Mass (GeV/#it{c}^2)"};
  ConfigurableAxis binsMult{"binsMult", {110, 0.0, 110.0}, "mult_{FT0M}"};
  ConfigurableAxis binsDCAz{"binsDCAz", {600, -3, 3}, ""};
  ConfigurableAxis binsDCAxy{"binsDCAxy", {300, -1.5, 1.5}, ""};
  ConfigurableAxis binsTPCXrows{"binsTPCXrows", {200, 0, 200}, ""};
  ConfigurableAxis binsnSigma{"binsnSigma", {130, -6.5, 6.5}, ""};
  ConfigurableAxis binsnTPCSignal{"binsnTPCSignal", {1000, 0, 1000}, ""};

  // Event selection cuts - (Temporary, need to fix!)
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(o2::framework::InitContext&)
  {
    // axes
    AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{binsEta, "#eta"};
    AxisSpec axisMassLambda1520{binsMass, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec axisMult{binsMult, "mult_{V0M}"};
    AxisSpec axisDCAz{binsDCAz, "DCA_{z}"};
    AxisSpec axisDCAxy{binsDCAxy, "DCA_{XY}"};
    AxisSpec axisTPCXrow{binsTPCXrows, "#Xrows_{TPC}"};
    AxisSpec pidQAAxis = {binsnSigma, "#sigma"};
    AxisSpec axisTPCSignal = {binsnTPCSignal, ""};

    if (doprocessData) {
      // Track QA before cuts
      //  --- Track
      histos.add("QA/QAbefore/Track/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPC_Nsigma_ka_only", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/TPC_Nsigma_pr_only", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAbefore/Track/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAbefore/Track/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAbefore/Track/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAbefore/Track/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAbefore/Track/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});

      if (isFilladditionalQA) {
        // TPC ncluster distirbutions
        histos.add("TPCncluster/TPCnclusterpr", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
        histos.add("TPCncluster/TPCnclusterka", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
        histos.add("TPCncluster/TPCnclusterPhipr", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});
        histos.add("TPCncluster/TPCnclusterPhika", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});

        // Multiplicity correlation calibrations
        histos.add("MultCalib/centglopr_before", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
        histos.add("MultCalib/centgloka_before", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
        histos.add("MultCalib/GloPVpr_before", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
        histos.add("MultCalib/GloPVka_before", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});

        // Multiplicity correlation calibrations
        histos.add("MultCalib/centglopr_after", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
        histos.add("MultCalib/centgloka_after", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
        histos.add("MultCalib/GloPVpr_after", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
        histos.add("MultCalib/GloPVka_after", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
      }

      // PID QA after cuts
      //  --- Kaon
      histos.add("QA/QAafter/Kaon/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Kaon/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Kaon/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Kaon/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Kaon/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Kaon/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Kaon/TPC_Signal_ka_all", "TPC Signal for Kaon;#it{p}_{T} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});

      //  --- Proton
      histos.add("QA/QAafter/Proton/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QA/QAafter/Proton/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH2F, {axisPt, pidQAAxis}});
      histos.add("QA/QAafter/Proton/dcaZ", "DCA_{Z} distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm);", HistType::kTH2F, {axisPt, axisDCAz});
      histos.add("QA/QAafter/Proton/dcaXY", "DCA_{XY} momentum distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCAxy});
      histos.add("QA/QAafter/Proton/TPC_CR", "# TPC Xrows distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
      histos.add("QA/QAafter/Proton/pT", "pT distribution of Protons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("QA/QAafter/Proton/eta", "#eta distribution of Protons; #eta; Counts;", {HistType::kTH1F, {axisEta}});
      histos.add("QA/QAafter/Proton/TPC_Signal_pr_all", "TPC Signal for Proton;#it{p}_{T} (GeV/#it{c});TPC Signal (A.U.)", {HistType::kTH2F, {axisPt, axisTPCSignal}});

      //  Mass QA 1D for quick check
      histos.add("Result/Data/lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/antilambda1520invmass", "Invariant mass of #Lambda(1520) K^{#mp}p^{#pm}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});   // K+ + Pr
      histos.add("Result/Data/lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}}); // K- + anti-Pr

      // 3d histogram
      histos.add("Result/Data/h3lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3antilambda1520invmass", "Invariant mass of #Lambda(1520) K^{#mp}p^{#pm}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});   // K+ + Pr
      histos.add("Result/Data/h3lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520}); // K- + anti-Pr
    }
    if (doprocessME) {
      histos.add("Result/Data/lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/h3lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
    }

    if (cEtaAssym) {
      histos.add("Result/Data/hlambda1520invmassUnlikeSignAside", "Invariant mass of #Lambda(1520) Unlike Sign A side", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/hlambda1520invmassLikeSignAside", "Invariant mass of #Lambda(1520) Like Sign A side", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/hlambda1520invmassUnlikeSignCside", "Invariant mass of #Lambda(1520) Unlike Sign C side", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/hlambda1520invmassLikeSignCside", "Invariant mass of #Lambda(1520) Like Sign C side", {HistType::kTH1F, {axisMassLambda1520}});

      histos.add("Result/Data/h3lambda1520invmassUnlikeSignAside", "Invariant mass of #Lambda(1520) Unlike Sign A side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLikeSignAside", "Invariant mass of #Lambda(1520) Like Sign A side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassUnlikeSignCside", "Invariant mass of #Lambda(1520) Unlike Sign C side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLikeSignCside", "Invariant mass of #Lambda(1520) Like Sign C side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      if (doprocessME) {
        histos.add("Result/Data/hlambda1520invmassMixedAside", "Invariant mass of #Lambda(1520) Mixed A side", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/hlambda1520invmassMixedCside", "Invariant mass of #Lambda(1520) Mixed C side", {HistType::kTH1F, {axisMassLambda1520}});

        histos.add("Result/Data/h3lambda1520invmassMixedAside", "Invariant mass of #Lambda(1520) Mixed A side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassMixedCside", "Invariant mass of #Lambda(1520) Mixed C side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      }
    }
    //}
    // MC QA
    if (doprocessMCTrue) {
      histos.add("Result/MC/Genlambda1520pt", "pT distribution of True MC #Lambda(1520)0", kTH2F, {axisPt, axisMult});
      histos.add("Result/MC/Genantilambda1520pt", "pT distribution of True MC Anti-#Lambda(1520)0", kTH2F, {axisPt, axisMult});
    }
    if (doprocessMC) {
      histos.add("QA/MC/trkDCAxy_pr", "DCAxy distribution of proton track candidates", HistType::kTH1F, {axisDCAxy});
      histos.add("QA/MC/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {axisDCAxy});
      histos.add("QA/MC/trkDCAz_pr", "DCAz distribution of proton track candidates", HistType::kTH1F, {axisDCAz});
      histos.add("QA/MC/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {axisDCAz});
      histos.add("QA/MC/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/MC/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/MC/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});
      histos.add("QA/MC/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {axisMult, axisPt, pidQAAxis}});

      histos.add("Result/MC/h3lambda1520Recoinvmass", "Invariant mass of Reconstructed MC #Lambda(1520)0", kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/MC/h3antilambda1520Recoinvmass", "Invariant mass of Reconstructed MC Anti-#Lambda(1520)0", kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/MC/lambda1520Reco", "pT distribution of Reconstructed MC #Lambda(1520)0", kTH2F, {axisPt, axisMult});
      histos.add("Result/MC/antilambda1520Reco", "pT distribution of Reconstructed MC Anti-#Lambda(1520)0", kTH2F, {axisPt, axisMult});
      histos.add("Result/MC/hlambda1520Recoinvmass", "Inv mass distribution of Reconstructed MC #Lambda(1520)", kTH1F, {axisMassLambda1520});
      histos.add("Result/MC/hantilambda1520Recoinvmass", "Inv mass distribution of Reconstructed MC Anti-#Lambda(1520)", kTH1F, {axisMassLambda1520});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    // if (collision.alias_bit(kTVXinTRD)) {
    //   // TRD triggered
    //   // return 0;
    // }
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    // if (multTrk < fMultCutLow->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultCutHigh->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
    //  return 0;

    return 1;
  }

  double massKa = MassKaonCharged;
  double massPr = MassProton;

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (cDCAr7SigCut) {
      if (std::abs(track.dcaXY()) > (0.004f + 0.0130f / (track.pt()))) // 7 - Sigma cut
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
        return false;
    }
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (track.tpcNClsCrossedRows() < cMinTPCncr)
      return false;
    if (cAddlTrackcut) {
      if (!track.passedITSRefit() || !track.passedTPCRefit())
        return false;
      if (track.tpcCrossedRowsOverFindableCls() < cMinRtpccut)
        return false;
      if (track.itsChi2NCl() > cMaxChi2ITScut)
        return false;
      if (track.tpcChi2NCl() > cMaxChi2TPCcut)
        return false;
      if (track.itsNCls() < cfgITScluster)
        return false;
    }
    if (cTPCNClsFound && (track.tpcNClsFound() < cMinTPCNClsFound))
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;

    return true;
  }

  // PID selection new PID method
  template <typename T>
  bool selectionnewPIDProton(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPr()) < cMaxTPCnSigmaProton) {
      if (cUseRejNsigma) {
        if (candidate.tpcNSigmaPi() > cRejNsigmaTpc && candidate.tpcNSigmaKa() > cRejNsigmaTpc) {
          tpcPIDPassed = true;
        }
      } else {
        tpcPIDPassed = true;
      }
    }
    if (candidate.hasTOF()) {
      if ((nsigmaCutCombinedProton > 0) && ((candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (nsigmaCutCombinedProton * nsigmaCutCombinedProton))) {
        if (cUseRejNsigma) {
          if (candidate.tofNSigmaPi() > cRejNsigmaTof && candidate.tofNSigmaKa() > cRejNsigmaTof) {
            tofPIDPassed = true;
          }
        } else {
          tofPIDPassed = true;
        }
      } else if (std::abs(candidate.tofNSigmaPr()) < cMaxTOFnSigmaProton) {
        if (cUseRejNsigma) {
          if (candidate.tofNSigmaPi() > cRejNsigmaTof && candidate.tofNSigmaKa() > cRejNsigmaTof) {
            tofPIDPassed = true;
          }
        } else {
          tofPIDPassed = true;
        }
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionnewPIDKaon(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) {
      if (cUseRejNsigma) {
        if (candidate.tpcNSigmaPi() > cRejNsigmaTpc && candidate.tpcNSigmaPr() > cRejNsigmaTpc) {
          tpcPIDPassed = true;
        }
      } else {
        tpcPIDPassed = true;
      }
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon) {
        if (cUseRejNsigma) {
          if (candidate.tofNSigmaPi() > cRejNsigmaTof && candidate.tofNSigmaPr() > cRejNsigmaTof) {
            tofPIDPassed = true;
          }
        } else {
          tofPIDPassed = true;
        }
      }
      if ((nsigmaCutCombinedKaon > 0) && ((candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < (nsigmaCutCombinedKaon * nsigmaCutCombinedKaon))) {
        if (cUseRejNsigma) {
          if (candidate.tofNSigmaPi() > cRejNsigmaTof && candidate.tofNSigmaPr() > cRejNsigmaTof) {
            tofPIDPassed = true;
          }
        } else {
          tofPIDPassed = true;
        }
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  // PID selection old PID method
  template <typename T>
  bool selectionoldPIDProton(const T& candidate)
  {
    auto vProtonTPCPIDpTintv = static_cast<std::vector<double>>(protonTPCPIDpTintv);
    vProtonTPCPIDpTintv.insert(vProtonTPCPIDpTintv.begin(), cMinPtcut);
    auto vProtonTPCPIDcuts = static_cast<std::vector<double>>(protonTPCPIDcuts);
    auto vProtonTOFPIDpTintv = static_cast<std::vector<double>>(protonTOFPIDpTintv);
    auto vProtonTOFPIDcuts = static_cast<std::vector<double>>(protonTOFPIDcuts);
    auto lengthOfprotonTPCPIDpTintv = static_cast<int>(vProtonTPCPIDpTintv.size());
    auto lengthOfprotonTOFPIDpTintv = static_cast<int>(vProtonTOFPIDpTintv.size());

    bool isTrk1Selected{true};

    // For Proton candidate:
    if (candidate.hasTOF()) {
      if (lengthOfprotonTOFPIDpTintv > 0) {
        if (candidate.pt() > vProtonTOFPIDpTintv[lengthOfprotonTOFPIDpTintv - 1]) {
          isTrk1Selected = false;
        } else {
          for (int i = 0; i < lengthOfprotonTOFPIDpTintv; i++) {
            if (candidate.pt() < vProtonTOFPIDpTintv[i]) {
              if (std::abs(candidate.tofNSigmaPr()) > vProtonTOFPIDcuts[i])
                isTrk1Selected = false;
              if (std::abs(candidate.tpcNSigmaPr()) > cMaxTPCnSigmaProtonVETO)
                isTrk1Selected = false;
            }
          }
        }
      }
    } else {
      if (lengthOfprotonTPCPIDpTintv > 0) {
        if (candidate.pt() > vProtonTPCPIDpTintv[lengthOfprotonTPCPIDpTintv - 1]) {
          isTrk1Selected = false;
        } else {
          for (int i = 0; i < lengthOfprotonTPCPIDpTintv; i++) {
            if (candidate.pt() > vProtonTPCPIDpTintv[i] && candidate.pt() < vProtonTPCPIDpTintv[i + 1]) {
              if (std::abs(candidate.tpcNSigmaPr()) > vProtonTPCPIDcuts[i])
                isTrk1Selected = false;
            }
          }
        }
      }
    }
    return isTrk1Selected;
  }

  template <typename T>
  bool selectionoldPIDKaon(const T& candidate)
  {
    auto vKaonTPCPIDpTintv = static_cast<std::vector<double>>(kaonTPCPIDpTintv);
    vKaonTPCPIDpTintv.insert(vKaonTPCPIDpTintv.begin(), cMinPtcut);
    auto vKaonTPCPIDcuts = static_cast<std::vector<double>>(kaonTPCPIDcuts);
    auto vKaonTOFPIDpTintv = static_cast<std::vector<double>>(kaonTOFPIDpTintv);
    auto vKaonTOFPIDcuts = static_cast<std::vector<double>>(kaonTOFPIDcuts);
    auto lengthOfkaonTPCPIDpTintv = static_cast<int>(vKaonTPCPIDpTintv.size());
    auto lengthOfkaonTOFPIDpTintv = static_cast<int>(vKaonTOFPIDpTintv.size());

    bool isTrk2Selected{true};

    // For Kaon candidate:
    if (candidate.hasTOF()) {
      if (lengthOfkaonTOFPIDpTintv > 0) {
        if (candidate.pt() > vKaonTOFPIDpTintv[lengthOfkaonTOFPIDpTintv - 1]) {
          isTrk2Selected = false;
        } else {
          for (int i = 0; i < lengthOfkaonTOFPIDpTintv; i++) {
            if (candidate.pt() < vKaonTOFPIDpTintv[i]) {
              if (std::abs(candidate.tofNSigmaKa()) > vKaonTOFPIDcuts[i])
                isTrk2Selected = false;
              if (std::abs(candidate.tpcNSigmaKa()) > cMaxTPCnSigmaKaonVETO)
                isTrk2Selected = false;
            }
          }
        }
      }
    } else {
      if (lengthOfkaonTPCPIDpTintv > 0) {
        if (candidate.pt() > vKaonTPCPIDpTintv[lengthOfkaonTPCPIDpTintv - 1]) {
          isTrk2Selected = false;
        } else {
          for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
            if (candidate.pt() > vKaonTPCPIDpTintv[i] && candidate.pt() < vKaonTPCPIDpTintv[i + 1]) {
              if (std::abs(candidate.tpcNSigmaKa()) > vKaonTPCPIDcuts[i])
                isTrk2Selected = false;
            }
          }
        }
      }
    }
    return isTrk2Selected;
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.cent();

    // Multiplicity correlation calibration plots
    if (isFilladditionalQA) {
      if constexpr (!IsMix) {
        histos.fill(HIST("MultCalib/centglopr_before"), multiplicity, dTracks1.size());
        histos.fill(HIST("MultCalib/centgloka_before"), multiplicity, dTracks2.size());
        histos.fill(HIST("MultCalib/GloPVpr_before"), dTracks1.size(), collision.multNTracksPV());
        histos.fill(HIST("MultCalib/GloPVka_before"), dTracks2.size(), collision.multNTracksPV());
      }
    }

    if (additionalEvsel && !eventSelected(collision, multiplicity)) {
      return;
    }
    if (isFilladditionalQA) {
      if constexpr (!IsMix) {
        histos.fill(HIST("MultCalib/centglopr_after"), multiplicity, dTracks1.size());
        histos.fill(HIST("MultCalib/centgloka_after"), multiplicity, dTracks2.size());
        histos.fill(HIST("MultCalib/GloPVpr_after"), dTracks1.size(), collision.multNTracksPV());
        histos.fill(HIST("MultCalib/GloPVka_after"), dTracks2.size(), collision.multNTracksPV());
      }
    }

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      //// Initialize variables
      // Trk1: Proton, Trk2: Kaon
      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();

      auto trk1ptPr = trk1.pt();
      auto trk1NSigmaPrTPC = trk1.tpcNSigmaPr();
      auto trk1NSigmaPrTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPr() : -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

      //// QA plots before the selection
      //  --- Track QA all
      if constexpr (!IsMix) {
        histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_pr_only"), trk1ptPr, trk1NSigmaPrTPC);
        }
        histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_ka_only"), trk2ptKa, trk2NSigmaKaTPC);
        }

        histos.fill(HIST("QA/QAbefore/Track/dcaZ"), trk1ptPr, trk1.dcaZ());
        histos.fill(HIST("QA/QAbefore/Track/dcaXY"), trk1ptPr, trk1.dcaXY());
        histos.fill(HIST("QA/QAbefore/Track/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAbefore/Track/pT"), trk1ptPr);
        histos.fill(HIST("QA/QAbefore/Track/eta"), trk1.eta());
      }

      //// Apply the pid selection
      if (cUseOnlyTOFTrackPr && !isTrk1hasTOF)
        continue;
      if (cUseOnlyTOFTrackKa && !isTrk2hasTOF)
        continue;
      if (cOldPIDcut) {
        if (!selectionoldPIDProton(trk1) || !selectionoldPIDKaon(trk2))
          continue;
      } else {
        if (!selectionnewPIDProton(trk1) || !selectionnewPIDKaon(trk2))
          continue;
      }

      //// QA plots after the selection
      if constexpr (!IsMix) { //  --- PID QA Proton
        histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
        histos.fill(HIST("QA/QAafter/Proton/TPC_Signal_pr_all"), trk1ptPr, trk1.tpcSignal());
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          histos.fill(HIST("QA/QAafter/Proton/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
        }
        if (!isTrk1hasTOF) {
          histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly"), trk1ptPr, trk1NSigmaPrTPC);
        }
        histos.fill(HIST("QA/QAafter/Proton/dcaZ"), trk1ptPr, trk1.dcaZ());
        histos.fill(HIST("QA/QAafter/Proton/dcaXY"), trk1ptPr, trk1.dcaXY());
        histos.fill(HIST("QA/QAafter/Proton/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Proton/pT"), trk1ptPr);
        histos.fill(HIST("QA/QAafter/Proton/eta"), trk1.eta());

        //  --- PID QA Kaon
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Signal_ka_all"), trk2ptKa, trk2.tpcSignal());
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/QAafter/Kaon/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        if (!isTrk2hasTOF) {
          histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly"), trk2ptKa, trk2NSigmaKaTPC);
        }
        histos.fill(HIST("QA/QAafter/Kaon/dcaZ"), trk2ptKa, trk2.dcaZ());
        histos.fill(HIST("QA/QAafter/Kaon/dcaXY"), trk2ptKa, trk2.dcaXY());
        histos.fill(HIST("QA/QAafter/Kaon/TPC_CR"), trk2ptKa, trk2.tpcNClsCrossedRows());
        histos.fill(HIST("QA/QAafter/Kaon/pT"), trk2ptKa);
        histos.fill(HIST("QA/QAafter/Kaon/eta"), trk2.eta());

        if (isFilladditionalQA) {
          // TPCncluster distributions
          histos.fill(HIST("TPCncluster/TPCnclusterpr"), trk1.tpcNClsFound());
          histos.fill(HIST("TPCncluster/TPCnclusterka"), trk2.tpcNClsFound());
          histos.fill(HIST("TPCncluster/TPCnclusterPhipr"), trk1.tpcNClsFound(), trk1.phi());
          histos.fill(HIST("TPCncluster/TPCnclusterPhika"), trk2.tpcNClsFound(), trk2.phi());
        }
      }

      // Apply kinematic cuts.
      if (cKinCuts) {
        TVector3 v1(trk1.px(), trk1.py(), trk1.pz());
        TVector3 v2(trk2.px(), trk2.py(), trk2.pz());
        float alpha = v1.Angle(v2);
        if (alpha > 1.4 && alpha < 2.4)
          continue;
      }

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPr);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (abs(lResonance.Rapidity()) > 0.5)
        continue;
      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (!IsMix) {
          if (trk1.sign() < 0) {
            histos.fill(HIST("Result/Data/lambda1520invmass"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmass"), collision.cent(), lResonance.Pt(), lResonance.M());
          } else {
            histos.fill(HIST("Result/Data/antilambda1520invmass"), lResonance.M());
            histos.fill(HIST("Result/Data/h3antilambda1520invmass"), collision.cent(), lResonance.Pt(), lResonance.M());
          }
          if (cEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassUnlikeSignAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassUnlikeSignAside"), collision.cent(), lResonance.Pt(), lResonance.M());
          } else if (cEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassUnlikeSignCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassUnlikeSignCside"), collision.cent(), lResonance.Pt(), lResonance.M());
          }
        } else {
          histos.fill(HIST("Result/Data/lambda1520invmassME"), lResonance.M());
          histos.fill(HIST("Result/Data/h3lambda1520invmassME"), collision.cent(), lResonance.Pt(), lResonance.M());
          if (cEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassMixedAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassMixedAside"), collision.cent(), lResonance.Pt(), lResonance.M());
          } else if (cEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassMixedCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassMixedCside"), collision.cent(), lResonance.Pt(), lResonance.M());
          }
        }

        // MC
        if constexpr (IsMC) {
          // LOG(info) << "trk1 pdgcode: " << trk1.pdgCode() << "trk2 pdgcode: " << trk2.pdgCode() << std::endl;

          if (abs(trk1.pdgCode()) != 2212 || abs(trk2.pdgCode()) != 321)
            continue;
          if (trk1.motherId() != trk2.motherId()) // Same mother
            continue;
          if (abs(trk1.motherPDG()) != 102134)
            continue;

          // Track selection check.
          histos.fill(HIST("QA/MC/trkDCAxy_pr"), trk1.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAxy_ka"), trk2.dcaXY());
          histos.fill(HIST("QA/MC/trkDCAz_pr"), trk1.dcaZ());
          histos.fill(HIST("QA/MC/trkDCAz_ka"), trk2.dcaZ());

          histos.fill(HIST("QA/MC/TPC_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_pr_all"), multiplicity, trk1ptPr, trk1NSigmaPrTOF);
          }
          histos.fill(HIST("QA/MC/TPC_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QA/MC/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          }

          // MC histograms
          if (trk1.motherPDG() > 0) {
            histos.fill(HIST("Result/MC/lambda1520Reco"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("Result/MC/hlambda1520Recoinvmass"), lResonance.M());
            histos.fill(HIST("Result/MC/h3lambda1520Recoinvmass"), collision.cent(), lResonance.Pt(), lResonance.M());
          } else {
            histos.fill(HIST("Result/MC/antilambda1520Reco"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("Result/MC/hantilambda1520Recoinvmass"), lResonance.M());
            histos.fill(HIST("Result/MC/h3antilambda1520Recoinvmass"), collision.cent(), lResonance.Pt(), lResonance.M());
          }
        }
      } else {
        if constexpr (!IsMix) {
          if (cEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassLikeSignAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLikeSignAside"), collision.cent(), lResonance.Pt(), lResonance.M());
          } else if (cEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassLikeSignCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLikeSignCside"), collision.cent(), lResonance.Pt(), lResonance.M());
          }
          // Like sign pair ++
          if (trk1.sign() > 0) {
            histos.fill(HIST("Result/Data/lambda1520invmassLSPP"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSPP"), collision.cent(), lResonance.Pt(), lResonance.M());
          } else { // Like sign pair --
            histos.fill(HIST("Result/Data/lambda1520invmassLSMM"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSMM"), collision.cent(), lResonance.Pt(), lResonance.M());
          }
        }
      }
    }
  }

  void processData(aod::ResoCollision& collision,
                   aod::ResoTracks const& resotracks)
  {
    fillHistograms<false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(lambda1520analysis, processData, "Process Event for data without partition", false);

  void processMC(aod::ResoCollision& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks)
  {
    fillHistograms<true, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(lambda1520analysis, processMC, "Process Event for MC Light without partition", false);

  void processMCTrue(aod::ResoCollision& collision, aod::ResoMCParents& resoParents)
  {
    auto multiplicity = collision.cent();
    // Not related to the real collisions
    for (auto& part : resoParents) {     // loop over all MC particles
      if (abs(part.pdgCode()) != 102134) // Lambda1520(0)
        continue;
      if (abs(part.y()) > 0.5) // rapidity cut
        continue;
      bool pass1 = false;
      bool pass2 = false;
      if (abs(part.daughterPDG1()) == 321 || abs(part.daughterPDG2()) == 321) { // At least one decay to Kaon
        pass2 = true;
      }
      if (abs(part.daughterPDG1()) == 2212 || abs(part.daughterPDG2()) == 2212) { // At least one decay to Proton
        pass1 = true;
      }
      if (!pass1 || !pass2) // If we have both decay products
        continue;
      // histos.fill(HIST("Result/MC/GenTruelambda1520pt"), part.pt());
      if (part.pdgCode() > 0)
        histos.fill(HIST("Result/MC/Genlambda1520pt"), part.pt(), multiplicity);
      else
        histos.fill(HIST("Result/MC/Genantilambda1520pt"), part.pt(), multiplicity);
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processMCTrue, "Process Event for MC only", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      fillHistograms<false, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(lambda1520analysis, processME, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambda1520analysis>(cfgc, TaskName{"lf-lambda1520analysis"})};
}
