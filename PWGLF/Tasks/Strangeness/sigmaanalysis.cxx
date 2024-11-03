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

using V0MCSigmas = soa::Join<aod::Sigma0Cores, aod::Sigma0CollRefs, aod::SigmaPhotonExtras, aod::SigmaLambdaExtras, aod::SigmaMCCores>;
using V0Sigmas = soa::Join<aod::Sigma0Cores, aod::Sigma0CollRefs, aod::SigmaPhotonExtras, aod::SigmaLambdaExtras>;

struct sigmaanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  // Configurable<bool> analyseSigma{"analyseSigma", false, "process Sigma-like candidates"};
  // Configurable<bool> analyseAntiSigma{"analyseAntiSigma", false, "process AntiSigma-like candidates"};

  // Analysis strategy:
  Configurable<bool> fUseMLSel{"fUseMLSel", false, "Flag to use ML selection. If False, the standard selection is applied."};
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
  Configurable<float> LambdaMinv0radius{"LambdaMinv0radius", 0.0, "Min V0 radius (cm)"};
  Configurable<float> LambdaMaxv0radius{"LambdaMaxv0radius", 40, "Max V0 radius (cm)"};
  Configurable<float> LambdaMinQt{"LambdaMinQt", 0.01, "Min lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMaxQt{"LambdaMaxQt", 0.17, "Max lambda qt value (AP plot) (GeV/c)"};
  Configurable<float> LambdaMinAlpha{"LambdaMinAlpha", 0.25, "Min lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMaxAlpha{"LambdaMaxAlpha", 1.0, "Max lambda alpha absolute value (AP plot)"};
  Configurable<float> LambdaMinv0cospa{"LambdaMinv0cospa", 0.95, "Min V0 CosPA"};
  Configurable<float> LambdaWindow{"LambdaWindow", 0.015, "Mass window around expected (in GeV/c2)"};
  Configurable<float> LambdaMaxRap{"LambdaMaxRap", 0.8, "Max lambda rapidity"};

  //// Photon standard criteria:
  // Configurable<float> PhotonMaxDauPseudoRap{"PhotonMaxDauPseudoRap", 0.9, "Max pseudorapidity of daughter tracks"};
  Configurable<float> PhotonDauMinPt{"PhotonDauMinPt", 0.0, "Min daughter pT (GeV/c)"};
  Configurable<float> PhotonMinDCADauToPv{"PhotonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
  Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> PhotonMinTPCCrossedRows{"PhotonMinTPCCrossedRows", 0, "Min daughter TPC Crossed Rows"};
  Configurable<float> PhotonMinTPCNSigmas{"PhotonMinTPCNSigmas", -7, "Min TPC NSigmas for daughters"};
  Configurable<float> PhotonMaxTPCNSigmas{"PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
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
  // TODO: Include PsiPair selection

  Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.5, "Max sigma0 rapidity"};

  // Axis
  // base properties
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {1000, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.05f, 1.151f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {600, -0.1f, 0.5f}, "M_{#Gamma}"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};

  // topological variable QA axes
  ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisCosPA{"axisCosPA", {200, 0.5f, 1.0f}, "Cosine of pointing angle"};
  ConfigurableAxis axisCandSel{"axisCandSel", {25, 0.5f, +25.5f}, "Candidate Selection"};

  // ML
  ConfigurableAxis MLProb{"MLOutput", {100, 0.0f, 1.0f}, ""};
  int nSigmaCandidates = 0;
  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {axisCentrality});

    // All candidates received
    histos.add("GeneralQA/h2dArmenterosAll", "h2dArmenterosAll", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/h2dArmenterosSelected", "h2dArmenterosSelected", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hMassSigma0All", "hMassSigma0All", kTH1F, {axisSigmaMass});

    // Candidates Counters
    histos.add("GeneralQA/hCandidateAnalysisSelection", "hCandidateAnalysisSelection", kTH1F, {axisCandSel});
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(2, "Photon Mass Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(3, "Photon DauPt Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(4, "Photon DCAToPV Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(5, "Photon DCADau Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(6, "Photon TPCCrossedRows Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(7, "Photon PosTPCSigma Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(8, "Photon NegTPCSigma Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(9, "Photon Pt Cut");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(10, "Photon Y Cut");
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
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(23, "Lambda Y Cut");

    // Photon Selection QA histos
    histos.add("GeneralQA/hPhotonMass", "hPhotonMass", kTH1F, {axisPhotonMass});
    histos.add("GeneralQA/hPhotonNegpT", "hPhotonNegpT", kTH1F, {axisPt});
    histos.add("GeneralQA/hPhotonPospT", "hPhotonPospT", kTH1F, {axisPt});
    histos.add("GeneralQA/hPhotonDCANegToPV", "hPhotonDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hPhotonDCAPosToPV", "hPhotonDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hPhotonDCADau", "hPhotonDCADau", kTH1F, {axisDCAdau});
    histos.add("GeneralQA/hPhotonPosTPCCR", "hPhotonPosTPCCR", kTH1F, {axisTPCrows});
    histos.add("GeneralQA/hPhotonNegTPCCR", "hPhotonNegTPCCR", kTH1F, {axisTPCrows});
    histos.add("GeneralQA/hPhotonPosTPCNSigma", "hPhotonPosTPCNSigma", kTH1F, {{30, -15.0f, 15.0f}});
    histos.add("GeneralQA/hPhotonNegTPCNSigma", "hPhotonNegTPCNSigma", kTH1F, {{30, -15.0f, 15.0f}});
    histos.add("GeneralQA/hPhotonpT", "hPhotonpT", kTH1F, {axisPt});
    histos.add("GeneralQA/hPhotonY", "hPhotonY", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hPhotonRadius", "hPhotonRadius", kTH1F, {axisRadius});
    histos.add("GeneralQA/hPhotonZ", "hPhotonZ", kTH1F, {{240, 0.0f, 120.0f}});
    histos.add("GeneralQA/h2dPhotonArmenteros", "h2dPhotonArmenteros", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hPhotonCosPA", "hPhotonCosPA", kTH1F, {axisCosPA});

    // Lambda Selection QA histos
    histos.add("GeneralQA/hLambdaMass", "hLambdaMass", kTH1F, {axisLambdaMass});
    histos.add("GeneralQA/hAntiLambdaMass", "hAntiLambdaMass", kTH1F, {axisLambdaMass});
    histos.add("GeneralQA/hLambdaDCANegToPV", "hLambdaDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hLambdaDCAPosToPV", "hLambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hLambdaRadius", "hLambdaRadius", kTH1F, {axisRadius});
    histos.add("GeneralQA/hLambdaDCADau", "hLambdaDCADau", kTH1F, {axisDCAdau});
    histos.add("GeneralQA/h2dLambdaArmenteros", "h2dLambdaArmenteros", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hLambdaCosPA", "hLambdaCosPA", kTH1F, {axisCosPA});
    histos.add("GeneralQA/hLambdaY", "hLambdaY", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hSigmaY", "hSigmaY", kTH1F, {axisRapidity});

    histos.add("GeneralQA/hPhotonMassSelected", "hPhotonMassSelected", kTH1F, {axisPhotonMass});
    histos.add("GeneralQA/hLambdaMassSelected", "hLambdaMassSelected", kTH1F, {axisLambdaMass});
    histos.add("GeneralQA/hAntiLambdaMassSelected", "hAntiLambdaMassSelected", kTH1F, {axisLambdaMass});
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
      // Event counter
      histos.add("hMCEventCentrality", "hMCEventCentrality", kTH1F, {axisCentrality});

      // Kinematic
      histos.add("MC/h3dMassSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
      histos.add("MC/h3dMassAntiSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});

      histos.add("MC/h2dArmenterosAll", "h2dArmenterosAll", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
      histos.add("MC/h2dArmenterosSelected", "h2dArmenterosSelected", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});

      // Sigma0 QA
      histos.add("MC/hMassSigma0All", "hMassSigma0All", kTH1F, {axisSigmaMass});
      histos.add("MC/hPtSigma0All", "hPtSigma0All", kTH1F, {axisPt});
      histos.add("MC/hMassSigma0", "hMassSigma0", kTH1F, {axisSigmaMass});
      histos.add("MC/hPtSigma0", "hPtSigma0", kTH1F, {axisPt});
      histos.add("MC/hMassAntiSigma0", "hMassAntiSigma0", kTH1F, {axisSigmaMass});
      histos.add("MC/hPtAntiSigma0", "hPtAntiSigma0", kTH1F, {axisPt});
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

      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 1.);
      histos.fill(HIST("GeneralQA/hPhotonMass"), cand.photonMass());
      if (TMath::Abs(cand.photonMass()) > PhotonMaxMass)
        return false;
      histos.fill(HIST("GeneralQA/hPhotonNegpT"), cand.photonNegPt());
      histos.fill(HIST("GeneralQA/hPhotonPospT"), cand.photonPosPt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 2.);
      if ((cand.photonPosPt() < PhotonDauMinPt) || (cand.photonNegPt() < PhotonDauMinPt))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonDCANegToPV"), cand.photonDCANegPV());
      histos.fill(HIST("GeneralQA/hPhotonDCAPosToPV"), cand.photonDCAPosPV());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 3.);
      if ((TMath::Abs(cand.photonDCAPosPV()) < PhotonMinDCADauToPv) || (TMath::Abs(cand.photonDCANegPV()) < PhotonMinDCADauToPv))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 4.);
      histos.fill(HIST("GeneralQA/hPhotonDCADau"), cand.photonDCADau());
      if (TMath::Abs(cand.photonDCADau()) > PhotonMaxDCAV0Dau)
        return false;
      histos.fill(HIST("GeneralQA/hPhotonPosTPCCR"), cand.photonPosTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hPhotonNegTPCCR"), cand.photonNegTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 5.);
      if ((cand.photonPosTPCCrossedRows() < PhotonMinTPCCrossedRows) || (cand.photonNegTPCCrossedRows() < PhotonMinTPCCrossedRows))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonPosTPCNSigma"), cand.photonPosTPCNSigma());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 6.);
      if ((cand.photonPosTPCNSigma() < PhotonMinTPCNSigmas) || (cand.photonPosTPCNSigma() > PhotonMaxTPCNSigmas))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonNegTPCNSigma"), cand.photonNegTPCNSigma());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 7.);
      if ((cand.photonNegTPCNSigma() < PhotonMinTPCNSigmas) || (cand.photonNegTPCNSigma() > PhotonMaxTPCNSigmas))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonpT"), cand.photonPt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 8.);
      if ((cand.photonPt() < PhotonMinPt) || (cand.photonPt() > PhotonMaxPt))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonY"), cand.photonY());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 9.);
      if ((TMath::Abs(cand.photonY()) > PhotonMaxRap))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonRadius"), cand.photonRadius());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 10.);
      if ((cand.photonRadius() < PhotonMinRadius) || (cand.photonRadius() > PhotonMaxRadius))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonZ"), cand.photonZconv());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 11.);
      if (TMath::Abs(cand.photonZconv()) > PhotonMaxZ)
        return false;
      histos.fill(HIST("GeneralQA/h2dPhotonArmenteros"), cand.photonAlpha(), cand.photonQt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 12.);
      if (cand.photonQt() > PhotonMaxQt)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 13.);
      if (TMath::Abs(cand.photonAlpha()) > PhotonMaxAlpha)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 14.);
      histos.fill(HIST("GeneralQA/hPhotonCosPA"), cand.photonCosPA());
      if (cand.photonCosPA() < PhotonMinV0cospa)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 15.);

      // Lambda selection
      histos.fill(HIST("GeneralQA/hLambdaMass"), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hAntiLambdaMass"), cand.antilambdaMass());
      if ((TMath::Abs(cand.lambdaMass() - 1.115683) > LambdaWindow) && (TMath::Abs(cand.antilambdaMass() - 1.115683) > LambdaWindow))
        return false;
      histos.fill(HIST("GeneralQA/hLambdaDCANegToPV"), cand.lambdaDCANegPV());
      histos.fill(HIST("GeneralQA/hLambdaDCAPosToPV"), cand.lambdaDCAPosPV());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 16.);
      if ((TMath::Abs(cand.lambdaDCAPosPV()) < LambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < LambdaMinDCANegToPv))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 17.);
      histos.fill(HIST("GeneralQA/hLambdaRadius"), cand.lambdaRadius());
      if ((cand.lambdaRadius() < LambdaMinv0radius) || (cand.lambdaRadius() > LambdaMaxv0radius))
        return false;
      histos.fill(HIST("GeneralQA/hLambdaDCADau"), cand.lambdaDCADau());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 18.);
      if (TMath::Abs(cand.lambdaDCADau()) > LambdaMaxDCAV0Dau)
        return false;
      histos.fill(HIST("GeneralQA/h2dLambdaArmenteros"), cand.lambdaAlpha(), cand.lambdaQt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 19.);
      if ((cand.lambdaQt() < LambdaMinQt) || (cand.lambdaQt() > LambdaMaxQt))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 20.);
      if ((TMath::Abs(cand.lambdaAlpha()) < LambdaMinAlpha) || (TMath::Abs(cand.lambdaAlpha()) > LambdaMaxAlpha))
        return false;
      histos.fill(HIST("GeneralQA/hLambdaCosPA"), cand.lambdaCosPA());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 21.);
      if (cand.lambdaCosPA() < LambdaMinv0cospa)
        return false;
      histos.fill(HIST("GeneralQA/hLambdaY"), cand.lambdaY());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 22.);
      if (TMath::Abs(cand.lambdaY()) > LambdaMaxRap)
        return false;
      histos.fill(HIST("GeneralQA/hSigmaY"), cand.sigmaRapidity());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 23.);
      if (TMath::Abs(cand.sigmaRapidity()) > SigmaMaxRap)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 24.);
      histos.fill(HIST("GeneralQA/hPhotonMassSelected"), cand.photonMass());
      histos.fill(HIST("GeneralQA/hLambdaMassSelected"), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hAntiLambdaMassSelected"), cand.antilambdaMass());
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
      if (sigma.isSigma() || sigma.isAntiSigma()) {
        histos.fill(HIST("MC/h2dArmenterosAll"), sigma.photonAlpha(), sigma.photonQt());
        histos.fill(HIST("MC/h2dArmenterosAll"), sigma.lambdaAlpha(), sigma.lambdaQt());
        histos.fill(HIST("MC/hMassSigma0All"), sigma.sigmaMass());
        histos.fill(HIST("MC/hPtSigma0All"), sigma.sigmapT());

        if (!processSigmaCandidate(sigma))
          continue;

        histos.fill(HIST("MC/h2dArmenterosSelected"), sigma.photonAlpha(), sigma.photonQt());
        histos.fill(HIST("MC/h2dArmenterosSelected"), sigma.lambdaAlpha(), sigma.lambdaQt());

        if (sigma.isSigma()) {
          histos.fill(HIST("MC/hMassSigma0"), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtSigma0"), sigma.sigmapT());
          histos.fill(HIST("MC/h3dMassSigma0"), coll.centFT0C(), sigma.sigmapT(), sigma.sigmaMass());
        } else {
          histos.fill(HIST("MC/hMassAntiSigma0"), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtAntiSigma0"), sigma.sigmapT());
          histos.fill(HIST("MC/h3dMassAntiSigma0"), coll.centFT0C(), sigma.sigmapT(), sigma.sigmaMass());
        }
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

      if (sigma.lambdaAlpha() > 0) {
        histos.fill(HIST("Sigma0/hMassSigma0"), sigma.sigmaMass());
        histos.fill(HIST("Sigma0/hPtSigma0"), sigma.sigmapT());
        histos.fill(HIST("Sigma0/hRapiditySigma0"), sigma.sigmaRapidity());
        histos.fill(HIST("Sigma0/h3dMassSigma0"), coll.centFT0C(), sigma.sigmapT(), sigma.sigmaMass());
      } else {
        histos.fill(HIST("AntiSigma0/hMassAntiSigma0"), sigma.sigmaMass());
        histos.fill(HIST("AntiSigma0/hPtAntiSigma0"), sigma.sigmapT());
        histos.fill(HIST("AntiSigma0/hRapidityAntiSigma0"), sigma.sigmaRapidity());
        histos.fill(HIST("AntiSigma0/h3dMassAntiSigma0"), coll.centFT0C(), sigma.sigmapT(), sigma.sigmaMass());
      }
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
