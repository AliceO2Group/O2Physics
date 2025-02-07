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
//
// This is a task that reads sigma0 tables (from sigma0builder) to perform analysis.
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  Sigma0 analysis task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFSigmaTables.h"
#include "CCDB/BasicCCDBManager.h"
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using std::array;
using V0MCSigmas = soa::Join<aod::Sigma0Cores, aod::SigmaPhotonExtras, aod::SigmaLambdaExtras, aod::SigmaMCCores>;
using V0Sigmas = soa::Join<aod::Sigma0Cores, aod::SigmaPhotonExtras, aod::SigmaLambdaExtras>;

struct sigmaanalysis {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Interaction rate selection:
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  Configurable<float> minIR{"minIR", -1, "Min Interaction Rate (kHz). Leave -1 if no selection desired."};
  Configurable<float> maxIR{"maxIR", -1, "Max Interaction Rate (kHz). Leave -1 if no selection desired."};

  // Analysis strategy:
  Configurable<bool> fUseMLSel{"fUseMLSel", false, "Flag to use ML selection. If False, the standard selection is applied."};
  Configurable<bool> fProcessMonteCarlo{"fProcessMonteCarlo", false, "Flag to process MC data."};
  Configurable<bool> fselLambdaTPCPID{"fselLambdaTPCPID", true, "Flag to select lambda-like candidates using TPC NSigma."};
  Configurable<bool> fselLambdaTOFPID{"fselLambdaTOFPID", false, "Flag to select lambda-like candidates using TOF NSigma."};
  Configurable<bool> doMCAssociation{"doMCAssociation", false, "Flag to process only signal candidates. Use only with processMonteCarlo!"};

  // For ML Selection
  Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
  Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
  Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};

  // For Standard Selection:
  //// Lambda standard criteria::
  Configurable<float> LambdaMinDCANegToPv{"LambdaMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> LambdaMinDCAPosToPv{"LambdaMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> ALambdaMinDCANegToPv{"ALambdaMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> ALambdaMinDCAPosToPv{"ALambdaMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> LambdaMaxDCAV0Dau{"LambdaMaxDCAV0Dau", 2.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> LambdaMinv0radius{"LambdaMinv0radius", 0.0, "Min V0 radius (cm)"};
  Configurable<float> LambdaMaxv0radius{"LambdaMaxv0radius", 40, "Max V0 radius (cm)"};
  Configurable<float> LambdaMinQt{"LambdaMinQt", 0.01, "Min lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMaxQt{"LambdaMaxQt", 0.17, "Max lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMinAlpha{"LambdaMinAlpha", 0.25, "Min lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMaxAlpha{"LambdaMaxAlpha", 1.0, "Max lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMinv0cospa{"LambdaMinv0cospa", 0.95, "Min V0 CosPA"};
  Configurable<float> LambdaMaxLifeTime{"LambdaMaxLifeTime", 30, "Max lifetime"};
  Configurable<float> LambdaWindow{"LambdaWindow", 0.015, "Mass window around expected (in GeV/c2)"};
  Configurable<float> LambdaMaxRap{"LambdaMaxRap", 0.8, "Max lambda rapidity"};
  Configurable<float> LambdaMaxDauEta{"LambdaMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
  Configurable<float> LambdaMaxTPCNSigmas{"LambdaMaxTPCNSigmas", 1e+9, "Max TPC NSigmas for daughters"};
  Configurable<float> LambdaPrMaxTOFNSigmas{"LambdaPrMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
  Configurable<float> LambdaPiMaxTOFNSigmas{"LambdaPiMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
  Configurable<int> LambdaMinTPCCrossedRows{"LambdaMinTPCCrossedRows", 50, "Min daughter TPC Crossed Rows"};
  Configurable<int> LambdaMinITSclusters{"LambdaMinITSclusters", 1, "minimum ITS clusters"};
  Configurable<bool> LambdaRejectPosITSafterburner{"LambdaRejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
  Configurable<bool> LambdaRejectNegITSafterburner{"LambdaRejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};

  //// Photon standard criteria:
  Configurable<int> Photonv0TypeSel{"Photonv0TypeSel", 7, "select on a certain V0 type (leave negative if no selection desired)"};
  Configurable<float> PhotonDauMinPt{"PhotonDauMinPt", 0.0, "Min daughter pT (GeV/c)"};
  Configurable<float> PhotonMinDCADauToPv{"PhotonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
  Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<int> PhotonMinTPCCrossedRows{"PhotonMinTPCCrossedRows", 30, "Min daughter TPC Crossed Rows"};
  Configurable<float> PhotonMinTPCNSigmas{"PhotonMinTPCNSigmas", -7, "Min TPC NSigmas for daughters"};
  Configurable<float> PhotonMaxTPCNSigmas{"PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
  Configurable<float> PiMaxTPCNSigmas{"PiMaxTPCNSigmas", 1, "Max TPC NSigmas for pi rejection"};
  Configurable<float> piMaxpT{"piMaxpT", 3.5, "Max pT for pi rejection"};
  Configurable<float> PhotonMinPt{"PhotonMinPt", 0.0, "Min photon pT (GeV/c)"};
  Configurable<float> PhotonMaxPt{"PhotonMaxPt", 50.0, "Max photon pT (GeV/c)"};
  Configurable<float> PhotonMaxRap{"PhotonMaxRap", 0.5, "Max photon rapidity"};
  Configurable<float> PhotonMinRadius{"PhotonMinRadius", 3.0, "Min photon conversion radius (cm)"};
  Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 115, "Max photon conversion radius (cm)"};
  Configurable<float> PhotonMaxZ{"PhotonMaxZ", 240, "Max photon conversion point z value (cm)"};
  Configurable<float> PhotonMaxQt{"PhotonMaxQt", 0.05, "Max photon qt value (AP plot) (GeV/c)"};
  Configurable<float> PhotonMaxAlpha{"PhotonMaxAlpha", 0.95, "Max photon alpha absolute value (AP plot)"};
  Configurable<float> PhotonMinV0cospa{"PhotonMinV0cospa", 0.80, "Min V0 CosPA"};
  Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};
  Configurable<float> PhotonPsiPairMax{"PhotonPsiPairMax", 1e+9, "maximum psi angle of the track pair"};
  Configurable<float> PhotonMaxDauEta{"PhotonMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
  Configurable<float> PhotonLineCutZ0{"PhotonLineCutZ0", 7.0, "The offset for the linecute used in the Z vs R plot"};

  Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.5, "Max sigma0 rapidity"};

  // Axis
  // base properties
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisDeltaPt{"axisDeltaPt", {100, 0.0, +1.0}, "#Delta(p_{T})"};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  ConfigurableAxis axisIRBinning{"axisIRBinning", {5000, 0, 1500}, "Binning for the interaction rate (kHz)"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {1000, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.05f, 1.151f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {600, -0.1f, 0.5f}, "M_{#Gamma}"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality, PID and other axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisNCls{"axisNCls", {8, -0.5, 7.5}, "NCls"};
  ConfigurableAxis axisChi2PerNcl{"axisChi2PerNcl", {80, -40, 40}, "Chi2 Per Ncl"};
  ConfigurableAxis axisTPCNSigma{"axisTPCNSigma", {120, -30, 30}, "TPC NSigma"};
  ConfigurableAxis axisTOFNSigma{"axisTOFNSigma", {120, -30, 30}, "TOF NSigma"};
  ConfigurableAxis axisLifetime{"axisLifetime", {200, 0, 200}, "Chi2 Per Ncl"};

  // topological variable QA axes
  ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisCosPA{"axisCosPA", {200, 0.5f, 1.0f}, "Cosine of pointing angle"};
  ConfigurableAxis axisPsiPair{"axisPsiPair", {500, -5.0f, 5.0f}, "Psipair for photons"};
  ConfigurableAxis axisCandSel{"axisCandSel", {31, 0.5f, +31.5f}, "Candidate Selection"};

  // ML
  ConfigurableAxis MLProb{"MLOutput", {100, 0.0f, 1.0f}, ""};
  int nSigmaCandidates = 0;
  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    // All candidates received
    histos.add("GeneralQA/hInteractionRate", "hInteractionRate", kTH1F, {axisIRBinning});
    histos.add("GeneralQA/hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2F, {axisCentrality, axisIRBinning});
    histos.add("GeneralQA/h2dArmenterosBeforeSel", "h2dArmenterosBeforeSel", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/h2dArmenterosAfterSel", "h2dArmenterosAfterSel", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hMassSigma0BeforeSel", "hMassSigma0BeforeSel", kTH1F, {axisSigmaMass});

    // Candidates Counters
    histos.add("GeneralQA/hCandidateAnalysisSelection", "hCandidateAnalysisSelection", kTH1F, {axisCandSel});
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(2, "Photon V0Type");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(3, "Photon Mass");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(4, "Photon DauPt");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(5, "Photon DCAToPV");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(6, "Photon DCADau");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(7, "Photon TPCCrossedRows");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(8, "Photon TPCNSigmaEl");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(9, "Photon TPCNSigmaPi");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(10, "Photon Pt");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(11, "Photon Y/Eta");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(12, "Photon Radius");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(13, "Photon RZ line");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(14, "Photon QT");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(15, "Photon Alpha");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(16, "Photon CosPA");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(17, "Photon PsiPair");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(18, "Lambda Radius");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(19, "Lambda DCADau");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(20, "Lambda QT");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(21, "Lambda Alpha");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(22, "Lambda CosPA");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(23, "Lambda Y/Eta");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(24, "Lambda TPCCrossedRows");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(25, "Lambda ITSNCls");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(26, "Lambda Lifetime");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(27, "Lambda/ALambda PID");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(28, "Lambda/ALambda DCAToPV");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(29, "Lambda/ALambda Mass");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(30, "Sigma Y");

    // Photon Selection QA histos
    histos.add("GeneralQA/hPhotonV0Type", "hPhotonV0Type", kTH1F, {{8, 0.5f, 8.5f}});
    histos.add("GeneralQA/hPhotonMass", "hPhotonMass", kTH1F, {axisPhotonMass});
    histos.add("GeneralQA/hPhotonNegpT", "hPhotonNegpT", kTH1F, {axisPt});
    histos.add("GeneralQA/hPhotonPospT", "hPhotonPospT", kTH1F, {axisPt});
    histos.add("GeneralQA/hPhotonDCANegToPV", "hPhotonDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hPhotonDCAPosToPV", "hPhotonDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hPhotonDCADau", "hPhotonDCADau", kTH1F, {axisDCAdau});
    histos.add("GeneralQA/hPhotonPosTPCCR", "hPhotonPosTPCCR", kTH1F, {axisTPCrows});
    histos.add("GeneralQA/hPhotonNegTPCCR", "hPhotonNegTPCCR", kTH1F, {axisTPCrows});
    histos.add("GeneralQA/h2dPhotonPosTPCNSigmaEl", "h2dPhotonPosTPCNSigmaEl", {HistType::kTH2F, {axisPt, {30, -15.0f, 15.0f}}});
    histos.add("GeneralQA/h2dPhotonNegTPCNSigmaEl", "h2dPhotonNegTPCNSigmaEl", {HistType::kTH2F, {axisPt, {30, -15.0f, 15.0f}}});
    histos.add("GeneralQA/h2dPhotonPosTPCNSigmaPi", "h2dPhotonPosTPCNSigmaPi", {HistType::kTH2F, {axisPt, {30, -15.0f, 15.0f}}});
    histos.add("GeneralQA/h2dPhotonNegTPCNSigmaPi", "h2dPhotonNegTPCNSigmaPi", {HistType::kTH2F, {axisPt, {30, -15.0f, 15.0f}}});
    histos.add("GeneralQA/hPhotonpT", "hPhotonpT", kTH1F, {axisPt});
    histos.add("GeneralQA/hPhotonY", "hPhotonY", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hPhotonPosEta", "hPhotonPosEta", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hPhotonNegEta", "hPhotonNegEta", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hPhotonRadius", "hPhotonRadius", kTH1F, {axisRadius});
    histos.add("GeneralQA/hPhotonZ", "hPhotonZ", kTH1F, {{240, 0.0f, 120.0f}});
    histos.add("GeneralQA/h2dRZCut", "h2dRZCut", {HistType::kTH2F, {{240, -120.0f, 120.0f}, axisRadius}});
    histos.add("GeneralQA/h2dRZPlane", "h2dRZPlane", {HistType::kTH2F, {{240, -120.0f, 120.0f}, axisRadius}});

    histos.add("GeneralQA/h2dPhotonArmenteros", "h2dPhotonArmenteros", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hPhotonCosPA", "hPhotonCosPA", kTH1F, {axisCosPA});
    histos.add("GeneralQA/hPhotonPsiPair", "hPhotonPsiPair", kTH1F, {axisPsiPair});

    // Lambda Selection QA histos
    histos.add("GeneralQA/hLambdaRadius", "hLambdaRadius", kTH1F, {axisRadius});
    histos.add("GeneralQA/hLambdaDCADau", "hLambdaDCADau", kTH1F, {axisDCAdau});
    histos.add("GeneralQA/h2dLambdaArmenteros", "h2dLambdaArmenteros", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hLambdaCosPA", "hLambdaCosPA", kTH1F, {axisCosPA});
    histos.add("GeneralQA/hLambdaY", "hLambdaY", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hLambdaPosEta", "hLambdaPosEta", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hLambdaNegEta", "hLambdaNegEta", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hLambdaPosTPCCR", "hLambdaPosTPCCR", kTH1F, {axisTPCrows});
    histos.add("GeneralQA/hLambdaNegTPCCR", "hLambdaNegTPCCR", kTH1F, {axisTPCrows});
    histos.add("GeneralQA/hLambdaPosITSCls", "hLambdaPosITSCls", kTH1F, {axisNCls});
    histos.add("GeneralQA/hLambdaNegITSCls", "hLambdaNegITSCls", kTH1F, {axisNCls});
    histos.add("GeneralQA/hLambdaPosChi2PerNc", "hLambdaPosChi2PerNc", kTH1F, {axisChi2PerNcl});
    histos.add("GeneralQA/hLambdaNegChi2PerNc", "hLambdaNegChi2PerNc", kTH1F, {axisChi2PerNcl});
    histos.add("GeneralQA/hLambdaLifeTime", "hLambdaLifeTime", kTH1F, {axisLifetime});

    histos.add("GeneralQA/h2dTPCvsTOFNSigma_LambdaPr", "h2dTPCvsTOFNSigma_LambdaPr", {HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});
    histos.add("GeneralQA/h2dTPCvsTOFNSigma_LambdaPi", "h2dTPCvsTOFNSigma_LambdaPi", {HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});
    histos.add("GeneralQA/hLambdaDCANegToPV", "hLambdaDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hLambdaDCAPosToPV", "hLambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hLambdaMass", "hLambdaMass", kTH1F, {axisLambdaMass});
    histos.add("GeneralQA/h2dTPCvsTOFNSigma_ALambdaPr", "h2dTPCvsTOFNSigma_ALambdaPr", {HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});
    histos.add("GeneralQA/h2dTPCvsTOFNSigma_ALambdaPi", "h2dTPCvsTOFNSigma_ALambdaPi", {HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});
    histos.add("GeneralQA/hALambdaDCANegToPV", "hALambdaDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hALambdaDCAPosToPV", "hALambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hAntiLambdaMass", "hAntiLambdaMass", kTH1F, {axisLambdaMass});

    histos.add("GeneralQA/hPhotonMassSelected", "hPhotonMassSelected", kTH1F, {axisPhotonMass});
    histos.add("GeneralQA/hLambdaMassSelected", "hLambdaMassSelected", kTH1F, {axisLambdaMass});
    histos.add("GeneralQA/hAntiLambdaMassSelected", "hAntiLambdaMassSelected", kTH1F, {axisLambdaMass});

    histos.add("GeneralQA/hSigmaY", "hSigmaY", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hSigmaOPAngle", "hSigmaOPAngle", kTH1F, {{140, 0.0f, +7.0f}});

    // Specific sigma0 QA
    histos.add("SigmaMassQA/h2dPhotonV0Type", "h2dPhotonV0Type", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonMass", "h2dPhotonMass", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonDaupT", "h2dPhotonDaupT", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonDCADauToPV", "h2dPhotonDCADauToPV", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonDCADau", "h2dPhotonDCADau", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonDauTPCCR", "h2dPhotonDauTPCCR", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonTPCNSigmaEl", "h2dPhotonTPCNSigmaEl", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonpT", "h2dPhotonpT", {HistType::kTH2F, {axisPt, axisSigmaMass}}); //
    histos.add("SigmaMassQA/h2dPhotonY", "h2dPhotonY", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonRadius", "h2dPhotonRadius", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dRZCut", "h2dRZCut", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonArmenteros", "h2dPhotonArmenteros", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonCosPA", "h2dPhotonCosPA", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dPhotonPsiPair", "h2dPhotonPsiPair", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaRadius", "h2dLambdaRadius", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaDCADau", "h2dLambdaDCADau", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaArmenteros", "h2dLambdaArmenteros", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaCosPA", "h2dLambdaCosPA", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaY", "h2dLambdaY", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaDauTPCCR", "h2dLambdaDauTPCCR", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaDauITSCls", "h2dLambdaDauITSCls", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaLifeTime", "h2dLambdaLifeTime", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dTPCvsTOFNSigma_Lambda", "h2dTPCvsTOFNSigma_Lambda", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaDCADauToPV", "h2dLambdaDCADauToPV", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dLambdaMass", "h2dLambdaMass", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dTPCvsTOFNSigma_ALambda", "h2dTPCvsTOFNSigma_ALambda", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dALambdaDCADauToPV", "h2dALambdaDCADauToPV", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dAntiLambdaMass", "h2dAntiLambdaMass", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaMassQA/h2dSigmaY", "h2dSigmaY", {HistType::kTH2F, {axisPt, axisSigmaMass}});

    // For Signal Extraction

    // Sigma0
    histos.add("Sigma0/h3dMassSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("Sigma0/hMassSigma0", "hMassSigma0", kTH1F, {axisSigmaMass});
    histos.add("Sigma0/hPtSigma0", "hPtSigma0", kTH1F, {axisPt});
    histos.add("Sigma0/hRapiditySigma0", "hRapiditySigma0", kTH1F, {axisRapidity});

    // AntiSigma0
    histos.add("AntiSigma0/h3dMassAntiSigma0", "h3dMassAntiSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("AntiSigma0/hMassAntiSigma0", "hMassAntiSigma0", kTH1F, {axisSigmaMass});
    histos.add("AntiSigma0/hPtAntiSigma0", "hPtAntiSigma0", kTH1F, {axisPt});
    histos.add("AntiSigma0/hRapidityAntiSigma0", "hRapidityAntiSigma0", kTH1F, {axisRapidity});

    if (fProcessMonteCarlo) {

      // Kinematic
      histos.add("MC/h3dMassSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
      histos.add("MC/h3dMassAntiSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});

      histos.add("MC/h2dArmenterosBeforeSel", "h2dArmenterosBeforeSel", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
      histos.add("MC/h2dArmenterosAfterSel", "h2dArmenterosAfterSel", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});

      // Sigma0 QA
      histos.add("MC/hMassSigma0BeforeSel", "hMassSigma0BeforeSel", kTH1F, {axisSigmaMass});
      histos.add("MC/hPtSigma0BeforeSel", "hPtSigma0BeforeSel", kTH1F, {axisPt});
      histos.add("MC/hMassSigma0", "hMassSigma0", kTH1F, {axisSigmaMass});
      histos.add("MC/hPtSigma0", "hPtSigma0", kTH1F, {axisPt});
      histos.add("MC/hMassAntiSigma0", "hMassAntiSigma0", kTH1F, {axisSigmaMass});
      histos.add("MC/hPtAntiSigma0", "hPtAntiSigma0", kTH1F, {axisPt});

      // For background decomposition
      histos.add("MC/h2dPtVsMassSigma_SignalBkg", "h2dPtVsMassSigma_SignalBkg", kTH2D, {axisPt, axisSigmaMass});
      histos.add("MC/h2dPtVsMassSigma_SignalOnly", "h2dPtVsMassSigma_SignalOnly", kTH2D, {axisPt, axisSigmaMass});
      histos.add("MC/h2dPtVsMassSigma_TrueDaughters", "h2dPtVsMassSigma_TrueDaughters", kTH2D, {axisPt, axisSigmaMass});
      histos.add("MC/h2dPtVsMassSigma_TrueGammaFakeLambda", "h2dPtVsMassSigma_TrueGammaFakeLambda", kTH2D, {axisPt, axisSigmaMass});
      histos.add("MC/h2dPtVsMassSigma_FakeGammaTrueLambda", "h2dPtVsMassSigma_FakeGammaTrueLambda", kTH2D, {axisPt, axisSigmaMass});
      histos.add("MC/h2dPtVsMassSigma_FakeDaughters", "h2dPtVsMassSigma_FakeDaughters", kTH2D, {axisPt, axisSigmaMass});
      histos.add("MC/h2dTrueDaughtersMatrix", "h2dTrueDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});

      // For new selection studies:
      //// Opening angle between daughters
      histos.add("MC/h2dPtVsOPAngle_SignalOnly", "h2dPtVsOPAngle_SignalOnly", kTH2D, {axisPt, {140, 0.0f, +7.0f}});
      histos.add("MC/h2dPtVsOPAngle_TrueDaughters", "h2dPtVsOPAngle_TrueDaughters", kTH2D, {axisPt, {140, 0.0f, +7.0f}});
      histos.add("MC/h2dPtVsMassSigma_AfterOPAngleSel", "h2dPtVsMassSigma_AfterOPAngleSel", kTH2D, {axisPt, axisSigmaMass});

      // For efficiency/Purity studies
      // Before any selection
      histos.add("MC/hPtTrueLambda_BeforeSel", "hPtTrueLambda_BeforeSel", kTH1F, {axisPt}); // Signal only
      histos.add("MC/hPtTrueAntiLambda_BeforeSel", "hPtTrueAntiLambda_BeforeSel", kTH1F, {axisPt}); // Signal only
      histos.add("MC/hPtTrueGamma_BeforeSel", "hPtTrueGamma_BeforeSel", kTH1F, {axisPt});   // Signal only
      histos.add("MC/hPtTrueSigma_BeforeSel", "hPtTrueSigma_BeforeSel", kTH1F, {axisPt});   // Signal only
      histos.add("MC/hPtTrueAntiSigma_BeforeSel", "hPtTrueAntiSigma_BeforeSel", kTH1F, {axisPt}); // Signal only
      histos.add("MC/hPtLambdaCand_BeforeSel", "hPtLambdaCand_BeforeSel", kTH1F, {axisPt}); // Bkg + Signal
      histos.add("MC/hPtGammaCand_BeforeSel", "hPtGammaCand_BeforeSel", kTH1F, {axisPt});   // Bkg + Signal
      histos.add("MC/hPtSigmaCand_BeforeSel", "hPtGammaCand_BeforeSel", kTH1F, {axisPt});   // Bkg + Signal

      // After analysis selections
      histos.add("MC/hPtTrueLambda_AfterSel", "hPtTrueLambda_AfterSel", kTH1F, {axisPt}); // Signal only
      histos.add("MC/hPtTrueAntiLambda_AfterSel", "hPtTrueAntiLambda_AfterSel", kTH1F, {axisPt}); // Signal only
      histos.add("MC/hPtTrueGamma_AfterSel", "hPtTrueGamma_AfterSel", kTH1F, {axisPt});   // Signal only
      histos.add("MC/hPtTrueSigma_AfterSel", "hPtTrueSigma_AfterSel", kTH1F, {axisPt});   // Signal only

      histos.add("MC/hPtLambdaCand_AfterSel", "hPtLambdaCand_AfterSel", kTH1F, {axisPt});
      histos.add("MC/hPtGammaCand_AfterSel", "hPtGammaCand_AfterSel", kTH1F, {axisPt});
      histos.add("MC/hPtSigmaCand_AfterSel", "hPtSigmaCand_AfterSel", kTH1F, {axisPt});

      // TPC vs TOF N Sigmas distributions
      histos.add("MC/h3dTPCvsTOFNSigma_LambdaPr", "h3dTPCvsTOFNSigma_LambdaPr", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
      histos.add("MC/h3dTPCvsTOFNSigma_LambdaPi", "h3dTPCvsTOFNSigma_LambdaPi", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
      histos.add("MC/h3dTPCvsTOFNSigma_TrueLambdaPr", "h3dTPCvsTOFNSigma_TrueLambdaPr", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
      histos.add("MC/h3dTPCvsTOFNSigma_TrueLambdaPi", "h3dTPCvsTOFNSigma_TrueLambdaPi", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
      histos.add("MC/h3dTPCvsTOFNSigma_TrueALambdaPr", "h3dTPCvsTOFNSigma_TrueALambdaPr", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
      histos.add("MC/h3dTPCvsTOFNSigma_TrueALambdaPi", "h3dTPCvsTOFNSigma_TrueALambdaPi", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});

      // QA of PID selections:
      //// TPC PID
      histos.add("MC/hPtTrueLambda_passedTPCPID", "hPtTrueLambda_passedTPCPID", kTH1F, {axisPt});
      histos.add("MC/hPtLambdaCandidates_passedTPCPID", "hPtLambdaCandidates_passedTPCPID", kTH1F, {axisPt});

      //// TOF PID
      histos.add("MC/hPtTrueLambda_passedTOFPID", "hPtTrueLambda_passedTOFPID", kTH1F, {axisPt});
      histos.add("MC/hPtLambdaCandidates_passedTOFPID", "hPtLambdaCandidates_passedTOFPID", kTH1F, {axisPt});

      //// TPC+TOF PID
      histos.add("MC/hPtTrueLambda_passedTPCTOFPID", "hPtTrueLambda_passedTPCTOFPID", kTH1F, {axisPt});
      histos.add("MC/hPtLambdaCandidates_passedTPCTOFPID", "hPtLambdaCandidates_passedTPCTOFPID", kTH1F, {axisPt});

      // pT Resolution:
      histos.add("MC/h3dLambdaPtResolution", "h3dLambdaPtResolution", kTH3F, {axisPt, axisDeltaPt, axisLambdaMass});
      histos.add("MC/h3dAntiLambdaPtResolution", "h3dAntiLambdaPtResolution", kTH3F, {axisPt, axisDeltaPt, axisLambdaMass});
      histos.add("MC/h3dGammaPtResolution", "h3dGammaPtResolution", kTH3F, {axisPt, axisDeltaPt, axisPhotonMass});
      histos.add("MC/h3dSigma0PtResolution", "h3dSigma0PtResolution", kTH3F, {axisPt, axisDeltaPt, axisSigmaMass});
    }
  }

  // Apply selections in sigma candidates
  template <typename TV0Object>
  bool processSigmaCandidate(TV0Object const& cand, bool isLambdalike)
  {
    if (fUseMLSel) {
      if ((cand.gammaBDTScore() == -1) || (cand.lambdaBDTScore() == -1) || (cand.antilambdaBDTScore() == -1)) {
        LOGF(fatal, "ML Score is not available! Please, enable gamma and lambda selection with ML in sigmabuilder!");
      }
      // Gamma selection:
      if (cand.gammaBDTScore() <= Gamma_MLThreshold)
        return false;

      // Lambda selection:
      if (cand.lambdaBDTScore() <= Lambda_MLThreshold)
        return false;

      // AntiLambda selection:
      if (cand.antilambdaBDTScore() <= AntiLambda_MLThreshold)
        return false;
    } else {

      // Photon Selections
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 1.);
      histos.fill(HIST("GeneralQA/hPhotonV0Type"), cand.photonV0Type());
      if (cand.photonV0Type() != Photonv0TypeSel && Photonv0TypeSel > -1)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonV0Type"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hPhotonMass"), cand.photonMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 2.);
      if (TMath::Abs(cand.photonMass()) > PhotonMaxMass)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonMass"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hPhotonNegpT"), cand.photonNegPt());
      histos.fill(HIST("GeneralQA/hPhotonPospT"), cand.photonPosPt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 3.);
      if ((cand.photonPosPt() < PhotonDauMinPt) || (cand.photonNegPt() < PhotonDauMinPt))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonDaupT"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hPhotonDCANegToPV"), cand.photonDCANegPV());
      histos.fill(HIST("GeneralQA/hPhotonDCAPosToPV"), cand.photonDCAPosPV());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 4.);
      if ((TMath::Abs(cand.photonDCAPosPV()) < PhotonMinDCADauToPv) || (TMath::Abs(cand.photonDCANegPV()) < PhotonMinDCADauToPv))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonDCADauToPV"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 5.);
      histos.fill(HIST("GeneralQA/hPhotonDCADau"), cand.photonDCADau());
      if (TMath::Abs(cand.photonDCADau()) > PhotonMaxDCAV0Dau)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonDCADau"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hPhotonPosTPCCR"), cand.photonPosTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hPhotonNegTPCCR"), cand.photonNegTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 6.);
      if ((cand.photonPosTPCCrossedRows() < PhotonMinTPCCrossedRows) || (cand.photonNegTPCCrossedRows() < PhotonMinTPCCrossedRows))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonDauTPCCR"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/h2dPhotonPosTPCNSigmaEl"), cand.photonPosPt(), cand.photonPosTPCNSigmaEl());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 7.);
      if (((cand.photonPosTPCNSigmaEl() < PhotonMinTPCNSigmas) || (cand.photonPosTPCNSigmaEl() > PhotonMaxTPCNSigmas)))
        return false;
      histos.fill(HIST("GeneralQA/h2dPhotonNegTPCNSigmaEl"), cand.photonNegPt(), cand.photonNegTPCNSigmaEl());
      if (((cand.photonNegTPCNSigmaEl() < PhotonMinTPCNSigmas) || (cand.photonNegTPCNSigmaEl() > PhotonMaxTPCNSigmas)))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonTPCNSigmaEl"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/h2dPhotonPosTPCNSigmaPi"), cand.photonPosPt(), cand.photonPosTPCNSigmaPi());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 8.);
      if (((TMath::Abs(cand.photonPosTPCNSigmaPi()) < PiMaxTPCNSigmas) && cand.photonPosPt() <= piMaxpT))
        return false;
      histos.fill(HIST("GeneralQA/h2dPhotonNegTPCNSigmaPi"), cand.photonNegPt(), cand.photonNegTPCNSigmaPi());
      if (((TMath::Abs(cand.photonNegTPCNSigmaPi()) < PiMaxTPCNSigmas) && cand.photonNegPt() <= piMaxpT))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonpT"), cand.photonPt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 9.);
      if ((cand.photonPt() < PhotonMinPt) || (cand.photonPt() > PhotonMaxPt))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonpT"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hPhotonY"), cand.photonY());
      histos.fill(HIST("GeneralQA/hPhotonPosEta"), cand.photonPosEta());
      histos.fill(HIST("GeneralQA/hPhotonNegEta"), cand.photonNegEta());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 10.);
      if ((TMath::Abs(cand.photonY()) > PhotonMaxRap) || (TMath::Abs(cand.photonPosEta()) > PhotonMaxDauEta) || (TMath::Abs(cand.photonNegEta()) > PhotonMaxDauEta))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonY"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hPhotonRadius"), cand.photonRadius());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 11.);
      if ((cand.photonRadius() < PhotonMinRadius) || (cand.photonRadius() > PhotonMaxRadius))
        return false;
      float photonRZLineCut = TMath::Abs(cand.photonZconv()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-PhotonMaxDauEta))) - PhotonLineCutZ0;
      histos.fill(HIST("SigmaMassQA/h2dPhotonRadius"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hPhotonZ"), cand.photonZconv());
      histos.fill(HIST("GeneralQA/h2dRZCut"), cand.photonRadius(), photonRZLineCut);
      histos.fill(HIST("GeneralQA/h2dRZPlane"), cand.photonZconv(), cand.photonRadius());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 12.);
      if ((TMath::Abs(cand.photonRadius()) < photonRZLineCut) || (TMath::Abs(cand.photonZconv()) > PhotonMaxZ))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dRZCut"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/h2dPhotonArmenteros"), cand.photonAlpha(), cand.photonQt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 13.);
      if (cand.photonQt() > PhotonMaxQt)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 14.);
      if (TMath::Abs(cand.photonAlpha()) > PhotonMaxAlpha)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonArmenteros"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 15.);
      histos.fill(HIST("GeneralQA/hPhotonCosPA"), cand.photonCosPA());
      if (cand.photonCosPA() < PhotonMinV0cospa)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonCosPA"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 16.);
      histos.fill(HIST("GeneralQA/hPhotonPsiPair"), cand.photonPsiPair());
      if (TMath::Abs(cand.photonPsiPair()) > PhotonPsiPairMax)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dPhotonPsiPair"), cand.sigmapT(), cand.sigmaMass());

      // Lambda selections
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 17.);
      histos.fill(HIST("GeneralQA/hLambdaRadius"), cand.lambdaRadius());
      if ((cand.lambdaRadius() < LambdaMinv0radius) || (cand.lambdaRadius() > LambdaMaxv0radius))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaRadius"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hLambdaDCADau"), cand.lambdaDCADau());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 18.);
      if (TMath::Abs(cand.lambdaDCADau()) > LambdaMaxDCAV0Dau)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaDCADau"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/h2dLambdaArmenteros"), cand.lambdaAlpha(), cand.lambdaQt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 19.);
      if ((cand.lambdaQt() < LambdaMinQt) || (cand.lambdaQt() > LambdaMaxQt))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 20.);
      if ((TMath::Abs(cand.lambdaAlpha()) < LambdaMinAlpha) || (TMath::Abs(cand.lambdaAlpha()) > LambdaMaxAlpha))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaArmenteros"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hLambdaCosPA"), cand.lambdaCosPA());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 21.);
      if (cand.lambdaCosPA() < LambdaMinv0cospa)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaCosPA"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hLambdaY"), cand.lambdaY());
      histos.fill(HIST("GeneralQA/hLambdaPosEta"), cand.lambdaPosEta());
      histos.fill(HIST("GeneralQA/hLambdaNegEta"), cand.lambdaNegEta());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 22.);
      if ((TMath::Abs(cand.lambdaY()) > LambdaMaxRap) || (TMath::Abs(cand.lambdaPosEta()) > LambdaMaxDauEta) || (TMath::Abs(cand.lambdaNegEta()) > LambdaMaxDauEta))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaY"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hLambdaPosTPCCR"), cand.lambdaPosTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hLambdaNegTPCCR"), cand.lambdaNegTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 23.);
      if ((cand.lambdaPosTPCCrossedRows() < LambdaMinTPCCrossedRows) || (cand.lambdaNegTPCCrossedRows() < LambdaMinTPCCrossedRows))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaDauTPCCR"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hLambdaPosITSCls"), cand.lambdaPosITSCls());
      histos.fill(HIST("GeneralQA/hLambdaNegITSCls"), cand.lambdaNegITSCls());
      histos.fill(HIST("GeneralQA/hLambdaPosChi2PerNc"), cand.lambdaPosChi2PerNcl());
      histos.fill(HIST("GeneralQA/hLambdaNegChi2PerNc"), cand.lambdaNegChi2PerNcl());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 24.);
      // check minimum number of ITS clusters + reject ITS afterburner tracks if requested
      bool posIsFromAfterburner = cand.lambdaPosChi2PerNcl() < 0;
      bool negIsFromAfterburner = cand.lambdaNegChi2PerNcl() < 0;
      if (cand.lambdaPosITSCls() < LambdaMinITSclusters && (!LambdaRejectPosITSafterburner || posIsFromAfterburner))
        return false;
      if (cand.lambdaNegITSCls() < LambdaMinITSclusters && (!LambdaRejectNegITSafterburner || negIsFromAfterburner))
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaDauITSCls"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hLambdaLifeTime"), cand.lambdaLifeTime());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 25.);
      if (cand.lambdaLifeTime() > LambdaMaxLifeTime)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dLambdaLifeTime"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 26.);
      // Separating lambda and antilambda selections:
      if (isLambdalike) { // Lambda selection
        histos.fill(HIST("GeneralQA/h2dTPCvsTOFNSigma_LambdaPr"), cand.lambdaPosPrTPCNSigma(), cand.lambdaPrTOFNSigma());
        histos.fill(HIST("GeneralQA/h2dTPCvsTOFNSigma_LambdaPi"), cand.lambdaNegPiTPCNSigma(), cand.lambdaPiTOFNSigma());

        // TPC Selection
        if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaPosPrTPCNSigma()) > LambdaMaxTPCNSigmas))
          return false;
        if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaNegPiTPCNSigma()) > LambdaMaxTPCNSigmas))
          return false;

        // TOF Selection
        if (fselLambdaTOFPID && (TMath::Abs(cand.lambdaPrTOFNSigma()) > LambdaPrMaxTOFNSigmas))
          return false;
        if (fselLambdaTOFPID && (TMath::Abs(cand.lambdaPiTOFNSigma()) > LambdaPiMaxTOFNSigmas))
          return false;

        histos.fill(HIST("SigmaMassQA/h2dTPCvsTOFNSigma_Lambda"), cand.sigmapT(), cand.sigmaMass());
        // DCA Selection
        histos.fill(HIST("GeneralQA/hLambdaDCANegToPV"), cand.lambdaDCANegPV());
        histos.fill(HIST("GeneralQA/hLambdaDCAPosToPV"), cand.lambdaDCAPosPV());
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 27.);
        if ((TMath::Abs(cand.lambdaDCAPosPV()) < LambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < LambdaMinDCANegToPv))
          return false;
        histos.fill(HIST("SigmaMassQA/h2dLambdaDCADauToPV"), cand.sigmapT(), cand.sigmaMass());

        // Mass Selection
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 28.);
        histos.fill(HIST("GeneralQA/hLambdaMass"), cand.lambdaMass());
        if (TMath::Abs(cand.lambdaMass() - 1.115683) > LambdaWindow)
          return false;
        histos.fill(HIST("SigmaMassQA/h2dLambdaMass"), cand.sigmapT(), cand.sigmaMass());

        // if constexpr (requires { cand.lambdaCandPDGCode(); }) {
        //   if (doPIDQA && passedTPC) {
        //     histos.fill(HIST("MC/hPtLambdaCandidates_passedTPCPID"), cand.lambdaPt());
        //     if (cand.lambdaCandPDGCode() == 3122)
        //       histos.fill(HIST("MC/hPtTrueLambda_passedTPCPID"), cand.lambdaPt());
        //   }
        //   if (doPIDQA && passedTOF) {
        //     histos.fill(HIST("MC/hPtLambdaCandidates_passedTOFPID"), cand.lambdaPt());
        //     if (cand.lambdaCandPDGCode() == 3122)
        //       histos.fill(HIST("MC/hPtTrueLambda_passedTOFPID"), cand.lambdaPt());
        //   }
        //   if (doPIDQA && passedTPC && passedTOF) {
        //     histos.fill(HIST("MC/hPtLambdaCandidates_passedTPCTOFPID"), cand.lambdaPt());
        //     if (cand.lambdaCandPDGCode() == 3122)
        //       histos.fill(HIST("MC/hPtTrueLambda_passedTPCTOFPID"), cand.lambdaPt());
        //   }
        // }
      } else { // AntiLambda selection
        histos.fill(HIST("GeneralQA/h2dTPCvsTOFNSigma_ALambdaPr"), cand.lambdaNegPrTPCNSigma(), cand.aLambdaPrTOFNSigma());
        histos.fill(HIST("GeneralQA/h2dTPCvsTOFNSigma_ALambdaPi"), cand.lambdaPosPiTPCNSigma(), cand.aLambdaPiTOFNSigma());

        // TPC Selection
        if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaPosPiTPCNSigma()) > LambdaMaxTPCNSigmas))
          return false;
        if (fselLambdaTPCPID && (TMath::Abs(cand.lambdaNegPrTPCNSigma()) > LambdaMaxTPCNSigmas))
          return false;

        // TOF Selection
        if (fselLambdaTOFPID && (TMath::Abs(cand.aLambdaPrTOFNSigma()) > LambdaPrMaxTOFNSigmas))
          return false;
        if (fselLambdaTOFPID && (TMath::Abs(cand.aLambdaPiTOFNSigma()) > LambdaPiMaxTOFNSigmas))
          return false;

        histos.fill(HIST("SigmaMassQA/h2dTPCvsTOFNSigma_ALambda"), cand.sigmapT(), cand.sigmaMass());
        // DCA Selection
        histos.fill(HIST("GeneralQA/hALambdaDCANegToPV"), cand.lambdaDCANegPV());
        histos.fill(HIST("GeneralQA/hALambdaDCAPosToPV"), cand.lambdaDCAPosPV());
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 27.);
        if ((TMath::Abs(cand.lambdaDCAPosPV()) < ALambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < ALambdaMinDCANegToPv))
          return false;

        histos.fill(HIST("SigmaMassQA/h2dALambdaDCADauToPV"), cand.sigmapT(), cand.sigmaMass());
        // Mass Selection
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 28.);
        histos.fill(HIST("GeneralQA/hAntiLambdaMass"), cand.antilambdaMass());
        if (TMath::Abs(cand.antilambdaMass() - 1.115683) > LambdaWindow)
          return false;
        histos.fill(HIST("SigmaMassQA/h2dAntiLambdaMass"), cand.sigmapT(), cand.sigmaMass());
      }

      // Sigma0 selection
      histos.fill(HIST("GeneralQA/hSigmaY"), cand.sigmaRapidity());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 29.);
      if (TMath::Abs(cand.sigmaRapidity()) > SigmaMaxRap)
        return false;
      histos.fill(HIST("SigmaMassQA/h2dSigmaY"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hSigmaOPAngle"), cand.sigmaOPAngle());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 30.);
    }
    return true;
  }

  void processMonteCarlo(V0MCSigmas const& sigmas)
  {
    for (auto& sigma : sigmas) { // selecting Sigma0-like candidates
      if (doMCAssociation && !(sigma.isSigma() || sigma.isAntiSigma())) {
        continue;
      }
      if (fGetIR) {
        double interactionRate = rateFetcher.fetch(ccdb.service, sigma.sigmaTimestamp(), sigma.sigmaRunNumber(), irSource) * 1.e-3;
        histos.fill(HIST("GeneralQA/hInteractionRate"), interactionRate);
        histos.fill(HIST("GeneralQA/hCentralityVsInteractionRate"), sigma.sigmaCentrality(), interactionRate);
        if ((maxIR != -1) && (minIR != -1) && ((interactionRate <= minIR) || (interactionRate >= maxIR))) {
          continue;
        }
      }

      // Filling histos before analysis selection
      histos.fill(HIST("MC/h2dArmenterosBeforeSel"), sigma.photonAlpha(), sigma.photonQt());
      histos.fill(HIST("MC/h2dArmenterosBeforeSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
      histos.fill(HIST("MC/hMassSigma0BeforeSel"), sigma.sigmaMass());
      histos.fill(HIST("MC/hPtSigma0BeforeSel"), sigma.sigmapT());
      histos.fill(HIST("MC/hPtGammaCand_BeforeSel"), sigma.photonPt());
      histos.fill(HIST("MC/hPtLambdaCand_BeforeSel"), sigma.lambdaPt());
      histos.fill(HIST("MC/hPtSigmaCand_BeforeSel"), sigma.sigmapT());

      if (sigma.photonCandPDGCode() == 22)
        histos.fill(HIST("MC/hPtTrueGamma_BeforeSel"), sigma.photonPt());
      if (sigma.lambdaCandPDGCode() == 3122)
        histos.fill(HIST("MC/hPtTrueLambda_BeforeSel"), sigma.lambdaPt());
      if (sigma.lambdaCandPDGCode() == -3122)
        histos.fill(HIST("MC/hPtTrueAntiLambda_BeforeSel"), sigma.lambdaPt());
      if (sigma.isSigma())
        histos.fill(HIST("MC/hPtTrueSigma_BeforeSel"), sigma.sigmapT());
      if (sigma.isAntiSigma())
        histos.fill(HIST("MC/hPtTrueAntiSigma_BeforeSel"), sigma.sigmapT());

      if (sigma.lambdaAlpha() > 0) { // Lambda Analysis
        if (!processSigmaCandidate(sigma, true))
          continue;

        // For Lambda PID Studies
        histos.fill(HIST("MC/hPtLambdaCand_AfterSel"), sigma.lambdaPt());
        histos.fill(HIST("MC/h3dTPCvsTOFNSigma_LambdaPr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma(), sigma.lambdaPt());
        histos.fill(HIST("MC/h3dTPCvsTOFNSigma_LambdaPi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma(), sigma.lambdaPt());

        if (sigma.lambdaCandPDGCode() == 3122) {
          histos.fill(HIST("MC/hPtTrueLambda_AfterSel"), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dLambdaPtResolution"), sigma.lambdaPt(), TMath::Abs((sigma.lambdaMCPt() - sigma.lambdaPt()) / sigma.lambdaMCPt()), sigma.lambdaMass()); // pT resolution
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueLambdaPr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma(), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueLambdaPi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma(), sigma.lambdaPt());
        }
        histos.fill(HIST("GeneralQA/hLambdaMassSelected"), sigma.lambdaMass());
        if (sigma.isSigma()) { // Signal study
          histos.fill(HIST("MC/h2dArmenterosAfterSel"), sigma.photonAlpha(), sigma.photonQt());
          histos.fill(HIST("MC/h2dArmenterosAfterSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
          histos.fill(HIST("MC/hMassSigma0"), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtSigma0"), sigma.sigmapT());
          histos.fill(HIST("MC/h3dMassSigma0"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
          histos.fill(HIST("MC/h3dSigma0PtResolution"), sigma.sigmapT(), TMath::Abs((sigma.sigmaMCPt() - sigma.sigmapT()) / sigma.sigmaMCPt()), sigma.sigmaMass()); // pT resolution
          histos.fill(HIST("MC/h2dPtVsMassSigma_SignalOnly"), sigma.sigmapT(), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtTrueSigma_AfterSel"), sigma.sigmapT());
        }
      } else { // AntiLambda Analysis
        if (!processSigmaCandidate(sigma, false))
          continue;

        if (sigma.lambdaCandPDGCode() == -3122) {
          histos.fill(HIST("MC/hPtTrueAntiLambda_AfterSel"), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dAntiLambdaPtResolution"), sigma.lambdaPt(), TMath::Abs((sigma.lambdaMCPt() - sigma.lambdaPt()) / sigma.lambdaMCPt()), sigma.antilambdaMass()); // pT resolution
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueALambdaPr"), sigma.lambdaNegPrTPCNSigma(), sigma.aLambdaPrTOFNSigma(), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueALambdaPi"), sigma.lambdaPosPiTPCNSigma(), sigma.aLambdaPiTOFNSigma(), sigma.lambdaPt());
        }

        histos.fill(HIST("GeneralQA/hAntiLambdaMassSelected"), sigma.antilambdaMass());
        if (sigma.isAntiSigma()) { // Signal study
          histos.fill(HIST("MC/h2dArmenterosAfterSel"), sigma.photonAlpha(), sigma.photonQt());
          histos.fill(HIST("MC/h2dArmenterosAfterSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
          histos.fill(HIST("MC/hMassAntiSigma0"), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtAntiSigma0"), sigma.sigmapT());
          histos.fill(HIST("MC/h3dMassAntiSigma0"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
          histos.fill(HIST("MC/h2dPtVsMassSigma_SignalOnly"), sigma.sigmapT(), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtTrueSigma_AfterSel"), sigma.sigmapT());
        }
      }

      // Fill histos after selection, please
      histos.fill(HIST("MC/hPtGammaCand_AfterSel"), sigma.photonPt());
      histos.fill(HIST("GeneralQA/hPhotonMassSelected"), sigma.photonMass());
      histos.fill(HIST("MC/hPtSigmaCand_AfterSel"), sigma.sigmapT());

      if (sigma.photonCandPDGCode() == 22) {
        histos.fill(HIST("MC/hPtTrueGamma_AfterSel"), sigma.photonPt());
        histos.fill(HIST("MC/h3dGammaPtResolution"), sigma.photonPt(), TMath::Abs((sigma.photonMCPt() - sigma.photonPt()) / sigma.photonMCPt()), sigma.photonMass()); // pT resolution
      }

      // For background studies:
      histos.fill(HIST("MC/h2dPtVsMassSigma_SignalBkg"), sigma.sigmapT(), sigma.sigmaMass());
      // Real Gamma x Real Lambda - but not from the same sigma0/antisigma0!
      if ((sigma.photonCandPDGCode() == 22) && ((sigma.lambdaCandPDGCode() == 3122) || (sigma.lambdaCandPDGCode() == -3122)) && !(sigma.isSigma()) && !(sigma.isAntiSigma())) {
        histos.fill(HIST("MC/h2dPtVsMassSigma_TrueDaughters"), sigma.sigmapT(), sigma.sigmaMass());
        histos.fill(HIST("MC/h2dTrueDaughtersMatrix"), sigma.lambdaCandPDGCodeMother(), sigma.photonCandPDGCodeMother());
        histos.fill(HIST("MC/h2dPtVsOPAngle_TrueDaughters"), sigma.sigmapT(), sigma.sigmaOPAngle());
      }
      // Real Gamma x fake Lambda
      if ((sigma.photonCandPDGCode() == 22) && (sigma.lambdaCandPDGCode() != 3122) && (sigma.lambdaCandPDGCode() != -3122))
        histos.fill(HIST("MC/h2dPtVsMassSigma_TrueGammaFakeLambda"), sigma.sigmapT(), sigma.sigmaMass());

      // Fake Gamma x Real Lambda
      if ((sigma.photonCandPDGCode() != 22) && ((sigma.lambdaCandPDGCode() == 3122) || (sigma.lambdaCandPDGCode() == -3122)))
        histos.fill(HIST("MC/h2dPtVsMassSigma_FakeGammaTrueLambda"), sigma.sigmapT(), sigma.sigmaMass());

      // Fake Gamma x Fake Lambda
      if ((sigma.photonCandPDGCode() != 22) && (sigma.lambdaCandPDGCode() != 3122) && (sigma.lambdaCandPDGCode() != -3122))
        histos.fill(HIST("MC/h2dPtVsMassSigma_FakeDaughters"), sigma.sigmapT(), sigma.sigmaMass());
    }
  }

  void processRealData(V0Sigmas const& sigmas)
  {
    for (auto& sigma : sigmas) { // selecting Sigma0-like candidates
      if (fGetIR) {
        double interactionRate = rateFetcher.fetch(ccdb.service, sigma.sigmaTimestamp(), sigma.sigmaRunNumber(), irSource) * 1.e-3;
        histos.fill(HIST("GeneralQA/hInteractionRate"), interactionRate);
        histos.fill(HIST("GeneralQA/hCentralityVsInteractionRate"), sigma.sigmaCentrality(), interactionRate);
        if ((maxIR != -1) && (minIR != -1) && ((interactionRate <= minIR) || (interactionRate >= maxIR))) {
          continue;
        }
      }
      histos.fill(HIST("GeneralQA/h2dArmenterosBeforeSel"), sigma.photonAlpha(), sigma.photonQt());
      histos.fill(HIST("GeneralQA/h2dArmenterosBeforeSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
      histos.fill(HIST("GeneralQA/hMassSigma0BeforeSel"), sigma.sigmaMass());

      nSigmaCandidates++;
      if (nSigmaCandidates % 100000 == 0) {
        LOG(info) << "Sigma0-like Candidates processed: " << nSigmaCandidates;
      }

      if (sigma.lambdaAlpha() > 0) {
        // Perform analysis selection for sigma0
        if (!processSigmaCandidate(sigma, true))
          continue;

        histos.fill(HIST("GeneralQA/h2dArmenterosAfterSel"), sigma.photonAlpha(), sigma.photonQt());
        histos.fill(HIST("GeneralQA/h2dArmenterosAfterSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
        histos.fill(HIST("GeneralQA/hLambdaMassSelected"), sigma.lambdaMass());
        histos.fill(HIST("Sigma0/hMassSigma0"), sigma.sigmaMass());
        histos.fill(HIST("Sigma0/hPtSigma0"), sigma.sigmapT());
        histos.fill(HIST("Sigma0/hRapiditySigma0"), sigma.sigmaRapidity());
        histos.fill(HIST("Sigma0/h3dMassSigma0"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());

        // Low-IR
        // Normal IR
        // High-IR
      } else {

        // Perform analysis selection for antisigma0
        if (!processSigmaCandidate(sigma, false))
          continue;

        histos.fill(HIST("GeneralQA/h2dArmenterosAfterSel"), sigma.photonAlpha(), sigma.photonQt());
        histos.fill(HIST("GeneralQA/h2dArmenterosAfterSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
        histos.fill(HIST("GeneralQA/hAntiLambdaMassSelected"), sigma.antilambdaMass());
        histos.fill(HIST("AntiSigma0/hMassAntiSigma0"), sigma.sigmaMass());
        histos.fill(HIST("AntiSigma0/hPtAntiSigma0"), sigma.sigmapT());
        histos.fill(HIST("AntiSigma0/hRapidityAntiSigma0"), sigma.sigmaRapidity());
        histos.fill(HIST("AntiSigma0/h3dMassAntiSigma0"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
      }
      histos.fill(HIST("GeneralQA/hPhotonMassSelected"), sigma.photonMass());
    }
  }

  // PROCESS_SWITCH(sigmaanalysis, processCounterQA, "Check standard counter correctness", true);
  PROCESS_SWITCH(sigmaanalysis, processMonteCarlo, "Do Monte-Carlo-based analysis", false);
  PROCESS_SWITCH(sigmaanalysis, processRealData, "Do real data analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigmaanalysis>(cfgc)};
}
