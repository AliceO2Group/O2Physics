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
// This is a task that employs the standard V0 tables and attempts to combine
// two V0s into a Sigma0 -> Lambda + gamma candidate.
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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFSigmaTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
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
using std::cout;
using std::endl;

using V0MCSigmas = soa::Join<aod::V0SigmaCandidates, aod::V0Sigma0CollRefs, aod::V0SigmaPhotonExtras, aod::V0SigmaLambdaExtras, aod::V0SigmaMCCandidates>;
using V0Sigmas = soa::Join<aod::V0SigmaCandidates, aod::V0Sigma0CollRefs, aod::V0SigmaPhotonExtras, aod::V0SigmaLambdaExtras>;

struct sigmaanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  // Configurable<bool> analyseSigma{"analyseSigma", false, "process Sigma-like candidates"};
  // Configurable<bool> analyseAntiSigma{"analyseAntiSigma", false, "process AntiSigma-like candidates"};

  // Analysis strategy:
  Configurable<bool> fUseMLSel{"fUseMLSel", true, "Flag to use ML selection. If False, the standard selection is applied."};
  Configurable<bool> fProcessMonteCarlo{"fProcessMonteCarlo", false, "Flag to process MC data."};

  // For ML Selection
  Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
  Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
  Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};

  // For Standard Selection:
  //// Lambda standard criteria::
  Configurable<float> LambdaMinDCANegToPv{"LambdaMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> LambdaMinDCAPosToPv{"LambdaMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> LambdaMaxDCAV0Dau{"LambdaMaxDCAV0Dau", 2.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> LambdaMinv0radius{"LambdaMinv0radius", 0.5, "Min V0 radius (cm)"};
  Configurable<float> LambdaMaxv0radius{"LambdaMaxv0radius", 180, "Max V0 radius (cm)"};
  Configurable<float> LambdaMinQt{"LambdaMinQt", 0.01, "Min lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMaxQt{"LambdaMaxQt", 0.17, "Max lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMinAlpha{"LambdaMinAlpha", 0.2, "Min lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMaxAlpha{"LambdaMaxAlpha", 0.9, "Max lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMinv0cospa{"LambdaMinv0cospa", 0.95, "Min V0 CosPA"};
  Configurable<float> LambdaWindow{"LambdaWindow", 0.01, "Mass window around expected (in GeV/c2)"};

  //// Photon standard criteria:
  Configurable<float> PhotonMaxDauPseudoRap{"PhotonMaxDauPseudoRap", 0.9, "Max pseudorapidity of daughter tracks"};
  Configurable<float> PhotonDauMinPt{"PhotonDauMinPt", 0.05, "Min daughter pT (GeV/c)"};
  Configurable<float> PhotonMinDCADauToPv{"PhotonMinDCADauToPv", 0.05, "Min DCA daughter To PV (cm)"};
  Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 2.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> PhotonMinTPCCrossedRows{"PhotonMinTPCCrossedRows", 55, "Min daughter TPC Crossed Rows"};
  Configurable<float> PhotonMinTPCNSigmas{"PhotonMinTPCNSigmas", -6, "Min TPC NSigmas for daughters"};
  Configurable<float> PhotonMaxTPCNSigmas{"PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
  Configurable<float> PhotonMinPt{"PhotonMinPt", 0.02, "Min photon pT (GeV/c)"};
  Configurable<float> PhotonMaxPt{"PhotonMaxPt", 50.0, "Max photon pT (GeV/c)"};
  Configurable<float> PhotonMaxPseudoRap{"PhotonMaxPseudoRap", 0.9, "Max photon pseudorapidity"};
  Configurable<float> PhotonMinRadius{"PhotonMinRadius", 0.5, "Min photon conversion radius (cm)"};
  Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 180, "Max photon conversion radius (cm)"};
  Configurable<float> PhotonMaxZ{"PhotonMaxZ", 240, "Max photon conversion point z value (cm)"};
  Configurable<float> PhotonMaxQt{"PhotonMaxQt", 0.06, "Max photon qt value (AP plot) (GeV/c)"};
  Configurable<float> PhotonMaxAlpha{"PhotonMaxAlpha", 0.95, "Max photon alpha absolute value (AP plot)"};
  Configurable<float> PhotonMinV0cospa{"PhotonMinV0cospa", 0.90, "Min V0 CosPA"};
  Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.15, "Max photon mass (GeV/c^{2})"};
  // TODO: Include PsiPair selection

  // Axis
  // base properties
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisRadius{"axisRadius", {200, 0.0f, 100.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -1.0f, 1.0f}, "Rapidity"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {200, 1.16f, 1.23f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.1f}, "M_{#Gamma}"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisITSclus{"axisITSclus", {7, 0.0f, 7.0f}, "N ITS Clusters"};

  // topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {20, 0.0f, 1.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {20, 0.0f, 2.0f}, "DCA (cm)"};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {20, 0.0f, 2.0f}, "pointing angle (rad)"};
  ConfigurableAxis axisV0Radius{"axisV0Radius", {20, 0.0f, 60.0f}, "V0 2D radius (cm)"};

  // ML
  ConfigurableAxis MLProb{"MLOutput", {100, 0.0f, 1.0f}, ""};
  int nSigmaCandidates = 0;
  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {axisCentrality});

    // For Signal Extraction
    histos.add("h3dMassSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("h3dMassAntiSigma0", "h3dMassAntiSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});

    // All candidates received
    // histos.add("GeneralQA/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
    // histos.add("GeneralQA/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
    // histos.add("GeneralQA/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
    // histos.add("GeneralQA/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
    // histos.add("GeneralQA/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    // histos.add("GeneralQA/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    // histos.add("GeneralQA/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("GeneralQA/h2dArmenterosAll", "h2dArmenterosAll", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/h2dArmenterosSelected", "h2dArmenterosSelected", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hMassSigma0All", "hMassSigma0All", kTH1F, {axisSigmaMass});

    // Candidates Counters
    histos.add("GeneralQA/hCandidateAnalysisSelection", "hCandidateAnalysisSelection", kTH1F, {{22, -0.5f, +21.5f}});
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(1, "Photon Mass Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(2, "Photon DauEta Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(3, "Photon DauPt Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(4, "Photon DCAToPV Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(5, "Photon DCADau Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(6, "Photon TPCCrossedRows Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(7, "Photon PosTPCSigma Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(8, "Photon NegTPCSigma Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(9, "Photon Pt Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(10, "Photon Eta Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(11, "Photon Radius Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(12, "Photon Zconv Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(13, "Photon QT Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(14, "Photon Alpha Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(15, "Photon CosPA Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(16, "Lambda Mass Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(17, "Lambda DCAToPV Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(18, "Lambda Radius Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(19, "Lambda DCADau Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(20, "Lambda QT Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(21, "Lambda Alpha Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(22, "Lambda CosPA Cut");

    // Sigma0 QA
    histos.add("Sigma0/hMassSigma0", "hMassSigma0", kTH1F, {axisSigmaMass});
    histos.add("Sigma0/hPtSigma0", "hPtSigma0", kTH1F, {axisPt});
    histos.add("Sigma0/hRapiditySigma0", "hRapiditySigma0", kTH1F, {axisRapidity});

    // Sigma0 QA
    histos.add("AntiSigma0/hMassSigma0", "hMassSigma0", kTH1F, {axisSigmaMass});
    histos.add("AntiSigma0/hPtSigma0", "hPtSigma0", kTH1F, {axisPt});
    histos.add("AntiSigma0/hRapiditySigma0", "hRapiditySigma0", kTH1F, {axisRapidity});

    // Gamma QA
    histos.add("Gamma/hMassGamma", "hMassGamma", kTH1F, {axisSigmaMass});
    histos.add("Gamma/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
    histos.add("Gamma/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
    histos.add("Gamma/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
    histos.add("Gamma/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
    histos.add("Gamma/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("Gamma/hMLOutputGamma", "hMLOutputGamma", kTH1F, {MLProb});
    histos.add("Gamma/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("Gamma/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("Gamma/h2dPhotonRadiusVsPt", "hPhotonRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});
    histos.add("Gamma/h2dPhotonMassVsMLScore", "h2dPhotonMassVsMLScore", {HistType::kTH2F, {MLProb, axisSigmaMass}});

    // Lambda QA
    histos.add("Lambda/hMassLambda", "hMassLambda", kTH1F, {axisLambdaMass});
    histos.add("Lambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
    histos.add("Lambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
    histos.add("Lambda/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
    histos.add("Lambda/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
    histos.add("Lambda/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("Lambda/hMLOutputLambda", "hMLOutputLambda", kTH1F, {MLProb});
    histos.add("Lambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("Lambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("Lambda/h2dLambdaRadiusVsPt", "hLambdaRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});
    histos.add("Lambda/h2dLambdaMassVsMLScore", "h2dLambdaMassVsMLScore", {HistType::kTH2F, {MLProb, axisLambdaMass}});

    // AntiLambda QA
    histos.add("AntiLambda/hMassAntiLambda", "hMassAntiLambda", kTH1F, {axisLambdaMass});
    histos.add("AntiLambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
    histos.add("AntiLambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
    histos.add("AntiLambda/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
    histos.add("AntiLambda/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
    histos.add("AntiLambda/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("AntiLambda/hMLOutputAntiLambda", "hMLOutputAntiLambda", kTH1F, {MLProb});
    histos.add("AntiLambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("AntiLambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    histos.add("AntiLambda/h2dAntiLambdaRadiusVsPt", "h2dAntiLambdaRadiusVsPt", {HistType::kTH2F, {axisPt, axisRadius}});
    histos.add("AntiLambda/h2dAntiLambdaMassVsMLScore", "h2dAntiLambdaMassVsMLScore", {HistType::kTH2F, {MLProb, axisLambdaMass}});

    if (fProcessMonteCarlo) {
      // Event counter
      histos.add("hMCEventCentrality", "hMCEventCentrality", kTH1F, {axisCentrality});

      // For Signal Extraction
      histos.add("h3dMCMassSigma0", "h3dMCMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});

      histos.add("GeneralQA/h2dMCArmenterosAll", "h2dMCArmenterosAll", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
      histos.add("GeneralQA/h2dMCArmenterosSelected", "h2dMCArmenterosSelected", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});

      // Sigma0 QA
      histos.add("GeneralQA/hMCMassSigma0All", "hMCMassSigma0All", kTH1F, {axisSigmaMass});
      histos.add("GeneralQA/hMCPtSigma0All", "hMCPtSigma0All", kTH1F, {axisPt});
      histos.add("Sigma0/hMCMassSigma0", "hMCMassSigma0", kTH1F, {axisSigmaMass});
      histos.add("Sigma0/hMCPtSigma0", "hMCPtSigma0", kTH1F, {axisPt});
    }
  }

  // Apply selections in sigma candidates
  template <typename TV0Object>
  bool processSigmaCandidate(TV0Object const& cand)
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
      if (TMath::Abs(cand.photonMass()) > PhotonMaxMass)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 0.);
      if ((TMath::Abs(cand.photonPosEta()) > PhotonMaxDauPseudoRap) || (TMath::Abs(cand.photonNegEta()) > PhotonMaxDauPseudoRap))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 1.);
      if ((cand.photonPosPt() < PhotonDauMinPt) || (cand.photonNegPt() < PhotonDauMinPt))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 2.);
      if ((TMath::Abs(cand.photonDCAPosPV()) < PhotonMinDCADauToPv) || (TMath::Abs(cand.photonDCANegPV()) < PhotonMinDCADauToPv))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 3.);
      if (cand.photonDCADau() > PhotonMaxDCAV0Dau)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 4.);
      if ((cand.photonPosTPCCrossedRows() < PhotonMinTPCCrossedRows) || (cand.photonNegTPCCrossedRows() < PhotonMinTPCCrossedRows))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 5.);
      if ((cand.photonPosTPCNSigma() < PhotonMinTPCNSigmas) || (cand.photonPosTPCNSigma() > PhotonMaxTPCNSigmas))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 6.);
      if ((cand.photonNegTPCNSigma() < PhotonMinTPCNSigmas) || (cand.photonNegTPCNSigma() > PhotonMaxTPCNSigmas))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 7.);
      if ((cand.photonPt() < PhotonMinPt) || (cand.photonPt() > PhotonMaxPt))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 8.);
      if (TMath::Abs(cand.photonEta()) > PhotonMaxPseudoRap)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 9.);
      if ((cand.photonRadius() < PhotonMinRadius) || (cand.photonRadius() > PhotonMaxRadius))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 10.);
      if (TMath::Abs(cand.photonZconv()) > PhotonMaxZ)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 11.);
      if (cand.photonQt() > PhotonMaxQt)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 12.);
      if (TMath::Abs(cand.photonAlpha()) > PhotonMaxAlpha)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 13.);
      if (cand.photonCosPA() < PhotonMinV0cospa)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 14.);

      // Lambda selection
      if (TMath::Abs(cand.lambdaMass() - 1.115683) > LambdaWindow)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 15.);
      if ((TMath::Abs(cand.lambdaDCAPosPV()) < LambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < LambdaMinDCANegToPv))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 16.);
      if ((cand.lambdaRadius() < LambdaMinv0radius) || (cand.lambdaRadius() > LambdaMaxv0radius))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 17.);
      if (TMath::Abs(cand.lambdaDCADau()) > LambdaMaxDCAV0Dau)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 18.);
      if ((cand.lambdaQt() < LambdaMinQt) || (cand.lambdaQt() > LambdaMaxQt))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 19.);
      if ((TMath::Abs(cand.lambdaAlpha()) < LambdaMinAlpha) || (TMath::Abs(cand.lambdaAlpha()) > LambdaMaxAlpha))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 20.);
      if (cand.lambdaCosPA() < LambdaMinv0cospa)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 21.);
    }

    return true;
  }

  // This process function cross-checks index correctness
  // void processCounterQA(V0Sigmas const& v0s)
  // {
  //   for (auto& gamma : v0s) {
  //     histos.fill(HIST("hGammaIndices"), gamma.globalIndex());
  //     histos.fill(HIST("hCollIndices"), gamma.straCollisionId());
  //     histos.fill(HIST("h2dIndices"), gamma.straCollisionId(), gamma.globalIndex());
  //   }
  // }

  void processMonteCarlo(aod::Sigma0Collision const& coll, V0MCSigmas const& v0s)
  {
    histos.fill(HIST("hMCEventCentrality"), coll.centFT0C());
    for (auto& sigma : v0s) { // selecting Sigma0-like candidates
      if (sigma.isSigma()) {
        histos.fill(HIST("GeneralQA/h2dMCArmenterosAll"), sigma.photonAlpha(), sigma.photonQt());
        histos.fill(HIST("GeneralQA/h2dMCArmenterosAll"), sigma.lambdaAlpha(), sigma.lambdaQt());
        histos.fill(HIST("GeneralQA/hMCMassSigma0All"), sigma.sigmaMass());
        histos.fill(HIST("GeneralQA/hMCPtSigma0All"), sigma.sigmapT());

        if (!processSigmaCandidate(sigma))
          continue;

        histos.fill(HIST("GeneralQA/h2dMCArmenterosSelected"), sigma.photonAlpha(), sigma.photonQt());
        histos.fill(HIST("GeneralQA/h2dMCArmenterosSelected"), sigma.lambdaAlpha(), sigma.lambdaQt());

        histos.fill(HIST("Sigma0/hMCMassSigma0"), sigma.sigmaMass());
        histos.fill(HIST("Sigma0/hMCPtSigma0"), sigma.sigmapT());

        histos.fill(HIST("h3dMCMassSigma0"), coll.centFT0C(), sigma.sigmapT(), sigma.sigmaMass());
      }
    }
  }

  void processRealData(aod::Sigma0Collision const& coll, V0Sigmas const& v0s)
  {
    histos.fill(HIST("hEventCentrality"), coll.centFT0C());
    for (auto& sigma : v0s) { // selecting Sigma0-like candidates
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), sigma.photonAlpha(), sigma.photonQt());
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), sigma.lambdaAlpha(), sigma.lambdaQt());
      histos.fill(HIST("GeneralQA/hMassSigma0All"), sigma.sigmaMass());

      nSigmaCandidates++;
      if (nSigmaCandidates % 50000 == 0) {
        LOG(info) << "Sigma0 Candidates processed: " << nSigmaCandidates;
      }
      if (!processSigmaCandidate(sigma))
        continue;

      histos.fill(HIST("GeneralQA/h2dArmenterosSelected"), sigma.photonAlpha(), sigma.photonQt());
      histos.fill(HIST("GeneralQA/h2dArmenterosSelected"), sigma.lambdaAlpha(), sigma.lambdaQt());

      histos.fill(HIST("Sigma0/hMassSigma0"), sigma.sigmaMass());
      histos.fill(HIST("Sigma0/hPtSigma0"), sigma.sigmapT());
      histos.fill(HIST("Sigma0/hRapiditySigma0"), sigma.sigmaRapidity());  
      histos.fill(HIST("h3dMassSigma0"), coll.centFT0C(), sigma.sigmapT(), sigma.sigmaMass());

      histos.fill(HIST("Gamma/hMassGamma"), sigma.lambdaMass());
      histos.fill(HIST("Lambda/hMassLambda"), sigma.photonMass());
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
