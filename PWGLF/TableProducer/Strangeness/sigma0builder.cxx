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
//  Sigma0 builder task
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
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFSigmaTables.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/CCDB/ctpRateFetcher.h"
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
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0StandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0TOFStandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0TOFDerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;

struct sigma0builder {    
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;
  
  //__________________________________________________
  // Sigma0 specific
  Produces<aod::Sigma0Cores> sigma0cores;                         // sigma0 candidates info for analysis  
  Produces<aod::Sigma0PhotonExtras> sigmaPhotonExtras;            // photons from sigma0 candidates info
  Produces<aod::Sigma0LambdaExtras> sigmaLambdaExtras;            // lambdas from sigma0 candidates info
  Produces<aod::SigmaCollRef> sigma0CollRefs;                     // references collisions from Sigma0Cores
  Produces<aod::Sigma0MCCores> sigma0mccores;                     // Reco sigma0 MC properties
  Produces<aod::Sigma0Gens> sigma0Gens;                           // Generated sigma0s
  Produces<aod::SigmaGenCollRef> sigma0GenCollRefs;               // references collisions from sigma0Gens

  //__________________________________________________
  // Pi0 specific
  Produces<aod::Pi0Cores> pi0cores;                    // pi0 candidates info for analysis  
  Produces<aod::Pi0CollRef> pi0coresRefs;          // references collisions from photonpair
  Produces<aod::Pi0CoresMC> pi0coresmc;                 // Reco pi0 MC properties
  Produces<aod::Pi0Gens> pi0Gens;                            // Generated pi0s
  Produces<aod::Pi0GenCollRef> pi0GenCollRefs;               // references collisions from pi0Gens

  //__________________________________________________
  // pack track quality but separte also afterburner
  // dynamic range: 0-31
  enum selection : int { hasTPC = 0,
                         hasITSTracker,
                         hasITSAfterburner,
                         hasTRD,
                         hasTOF };

  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> fillQAhistos{"fillQAhistos", false, "if true, fill QA histograms"};  
  Configurable<bool> doAssocStudy{"doAssocStudy", false, "Do v0 to collision association study."};
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<bool> fIRCrashOnNull{"fIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  // Tables to fill
  Configurable<bool> fillPi0Tables{"fillPi0Tables", false, "fill pi0 tables for QA"};
  Configurable<bool> fillSigma0Tables{"fillSigma0Tables", true, "fill sigma0 tables for analysis"};  
  
  // For ML Selection
  Configurable<bool> useMLScores{"useMLScores", false, "use ML scores to select candidates"};
  Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
  Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
  Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};

  // For standard approach:
  //// Lambda criteria:
  Configurable<float> V0Rapidity{"V0Rapidity", 0.5, "v0 rapidity"};

  Configurable<float> LambdaDauPseudoRap{"LambdaDauPseudoRap", 1.5, "Max pseudorapidity of daughter tracks"};
  Configurable<float> LambdaMinDCANegToPv{"LambdaMinDCANegToPv", 0.0, "min DCA Neg To PV (cm)"};
  Configurable<float> LambdaMinDCAPosToPv{"LambdaMinDCAPosToPv", 0.0, "min DCA Pos To PV (cm)"};
  Configurable<float> LambdaMaxDCAV0Dau{"LambdaMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> LambdaMinv0radius{"LambdaMinv0radius", 0.0, "Min V0 radius (cm)"};
  Configurable<float> LambdaMaxv0radius{"LambdaMaxv0radius", 60, "Max V0 radius (cm)"};
  Configurable<float> LambdaWindow{"LambdaWindow", 0.05, "Mass window around expected (in GeV/c2)"};

  //// Photon criteria:
  Configurable<float> PhotonMaxDauPseudoRap{"PhotonMaxDauPseudoRap", 1.5, "Max pseudorapidity of daughter tracks"};
  Configurable<float> PhotonMinDCAToPv{"PhotonMinDCAToPv", 0.0, "Min DCA daughter To PV (cm)"};
  Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> PhotonMinRadius{"PhotonMinRadius", 0.0, "Min photon conversion radius (cm)"};
  Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 240, "Max photon conversion radius (cm)"};
  Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.3, "Max photon mass (GeV/c^{2})"};

  //// Sigma0 criteria:
  Configurable<float> Sigma0Window{"Sigma0Window", 0.1, "Mass window around expected (in GeV/c2)"};
  Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.8, "Max sigma0 rapidity"};
  
  //// Extras:
  Configurable<float> Pi0MaxRap{"Pi0MaxRap", 0.8, "Max Pi0 Rapidity"};
  Configurable<float> Pi0MassWindow{"Pi0MassWindow", 0.115, "Mass window around expected (in GeV/c2)"};


  Configurable<float> GenMaxRap{"GenMaxRap", 0.5, "Max generated particle rapidity"};

  // Axis
  // base properties
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {500, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.05f, 1.151f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.5f}, "M_{#Gamma}"};
  ConfigurableAxis axisPi0Mass{"axisPi0Mass", {200, 0.08f, 0.18f}, "M_{#Pi^{0}}"};
  ConfigurableAxis axisK0SMass{"axisK0SMass", {200, 0.4f, 0.6f}, "M_{K^{0}}"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};

  // topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisXY{"axisXY", {120, -120.0f, 120.0f}, "XY axis"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisPA{"axisPA", {100, 0.0f, 1}, "Pointing angle"};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  ConfigurableAxis axisCandSel{"axisCandSel", {7, 0.5f, +7.5f}, "Candidate Selection"};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "N_{ch}"};
  ConfigurableAxis axisIRBinning{"axisIRBinning", {151, -10, 1500}, "Binning for the interaction rate (kHz)"};
  
  void init(InitContext const&)
  {        
    LOGF(info, "Initializing now: cross-checking correctness...");
    if (doprocessRealData +
        doprocessRealDataWithTOF +
        doprocessMonteCarlo +
        doprocessMonteCarloWithTOF > 1) {
      LOGF(fatal, "You have enabled more than one process function. Please check your configuration! Aborting now.");
    }

    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    histos.add("hEventCentrality", "hEventCentrality", kTH1D, {axisCentrality});

    histos.add("PhotonSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "Photon Mass Cut");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "Photon Eta/Y Cut");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(4, "Photon DCAToPV Cut");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(5, "Photon DCADau Cut");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(6, "Photon Radius Cut");

    histos.add("PhotonSel/hPhotonMass", "hPhotonMass", kTH1F, {axisPhotonMass});
    histos.add("PhotonSel/hPhotonNegEta", "hPhotonNegEta", kTH1F, {axisRapidity});
    histos.add("PhotonSel/hPhotonPosEta", "hPhotonPosEta", kTH1F, {axisRapidity});
    histos.add("PhotonSel/hPhotonY", "hPhotonY", kTH1F, {axisRapidity});
    histos.add("PhotonSel/hPhotonDCANegToPV", "hPhotonDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("PhotonSel/hPhotonDCAPosToPV", "hPhotonDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("PhotonSel/hPhotonDCADau", "hPhotonDCADau", kTH1F, {axisDCAdau});
    histos.add("PhotonSel/hPhotonRadius", "hPhotonRadius", kTH1F, {axisRadius});
    histos.add("PhotonSel/h3dPhotonMass", "h3dPhotonMass", kTH3D, {axisCentrality, axisPt, axisPhotonMass});

    histos.add("LambdaSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "Lambda Mass Cut");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "Lambda Eta/Y Cut");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(4, "Lambda DCAToPV Cut");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(5, "Lambda Radius Cut");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(6, "Lambda DCADau Cut");

    histos.add("LambdaSel/hLambdaMass", "hLambdaMass", kTH1F, {axisLambdaMass});
    histos.add("LambdaSel/hAntiLambdaMass", "hAntiLambdaMass", kTH1F, {axisLambdaMass});
    histos.add("LambdaSel/hLambdaNegEta", "hLambdaNegEta", kTH1F, {axisRapidity});
    histos.add("LambdaSel/hLambdaPosEta", "hLambdaPosEta", kTH1F, {axisRapidity});
    histos.add("LambdaSel/hLambdaY", "hLambdaY", kTH1F, {axisRapidity});
    histos.add("LambdaSel/hLambdaDCANegToPV", "hLambdaDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("LambdaSel/hLambdaDCAPosToPV", "hLambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("LambdaSel/hLambdaDCADau", "hLambdaDCADau", kTH1F, {axisDCAdau});
    histos.add("LambdaSel/hLambdaRadius", "hLambdaRadius", kTH1F, {axisRadius});
    histos.add("LambdaSel/h3dLambdaMass", "h3dLambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});
    histos.add("LambdaSel/h3dALambdaMass", "h3dALambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});

    histos.add("SigmaSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
    histos.get<TH1>(HIST("SigmaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("SigmaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "Sigma Mass Window");
    histos.get<TH1>(HIST("SigmaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "Sigma Y Window");

    // For selection:
    histos.add("SigmaSel/h3dMassSigma0BeforeSel", "h3dMassSigma0BeforeSel", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("SigmaSel/hSigmaMass", "hSigmaMass", kTH1F, {axisSigmaMass});
    histos.add("SigmaSel/hSigmaMassWindow", "hSigmaMassWindow", kTH1F, {{200, -0.09f, 0.11f}});
    histos.add("SigmaSel/hSigmaY", "hSigmaY", kTH1F, {axisRapidity});
    histos.add("SigmaSel/hSigmaMassSelected", "hSigmaMassSelected", kTH1F, {axisSigmaMass});
    histos.add("SigmaSel/h3dMassSigma0AfterSel", "h3dMassSigma0AfterSel", kTH3D, {axisCentrality, axisPt, axisSigmaMass});

    if (fillQAhistos) {
      histos.add("GeneralQA/h2dMassGammaVsK0S", "h2dMassGammaVsK0S", kTH2D, {axisPhotonMass, axisK0SMass});
      histos.add("GeneralQA/h2dMassLambdaVsK0S", "h2dMassLambdaVsK0S", kTH2D, {axisLambdaMass, axisK0SMass});
      histos.add("GeneralQA/h2dMassGammaVsLambda", "h2dMassGammaVsLambda", kTH2D, {axisPhotonMass, axisLambdaMass});
      histos.add("GeneralQA/h2dMassLambdaVsGamma", "h2dMassLambdaVsGamma", kTH2D, {axisLambdaMass, axisPhotonMass});
      histos.add("GeneralQA/h3dMassSigma0VsDaupTs", "h3dMassSigma0VsDaupTs", kTH3F, {axisPt, axisPt, axisSigmaMass});
      histos.add("GeneralQA/h2dMassGammaVsK0SAfterMassSel", "h2dMassGammaVsK0SAfterMassSel", kTH2D, {axisPhotonMass, axisK0SMass});
      histos.add("GeneralQA/h2dMassLambdaVsK0SAfterMassSel", "h2dMassLambdaVsK0SAfterMassSel", kTH2D, {axisLambdaMass, axisK0SMass});
      histos.add("GeneralQA/h2dMassGammaVsLambdaAfterMassSel", "h2dMassGammaVsLambdaAfterMassSel", kTH2D, {axisPhotonMass, axisLambdaMass});
      histos.add("GeneralQA/h2dV0XY", "h2dV0XY", kTH2F, {axisXY, axisXY});
    }

    if (doAssocStudy && (doprocessMonteCarlo || doprocessMonteCarloWithTOF)) {
      histos.add("V0AssoQA/h2dIRVsPt_TrueGamma", "h2dIRVsPt_TrueGamma", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueGamma", "h3dPAVsIRVsPt_TrueGamma", kTH3F, {axisPA, axisIRBinning, axisPt});
      histos.add("V0AssoQA/h2dIRVsPt_TrueGamma_BadCollAssig", "h2dIRVsPt_TrueGamma_BadCollAssig", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueGamma_BadCollAssig", "h3dPAVsIRVsPt_TrueGamma_BadCollAssig", kTH3F, {axisPA, axisIRBinning, axisPt});

      histos.add("V0AssoQA/h2dIRVsPt_TrueLambda", "h2dIRVsPt_TrueLambda", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueLambda", "h3dPAVsIRVsPt_TrueLambda", kTH3F, {axisPA, axisIRBinning, axisPt});
      histos.add("V0AssoQA/h2dIRVsPt_TrueLambda_BadCollAssig", "h2dIRVsPt_TrueLambda_BadCollAssig", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueLambda_BadCollAssig", "h3dPAVsIRVsPt_TrueLambda_BadCollAssig", kTH3F, {axisPA, axisIRBinning, axisPt});
    }

    // MC
    if (doprocessMonteCarlo || doprocessMonteCarloWithTOF) {
      histos.add("MC/h2dPtVsCentralityBeforeSel_MCAssocGamma", "h2dPtVsCentralityBeforeSel_MCAssocGamma", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dPtVsCentralityBeforeSel_MCAssocLambda", "h2dPtVsCentralityBeforeSel_MCAssocLambda", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dPtVsCentralityBeforeSel_MCAssocALambda", "h2dPtVsCentralityBeforeSel_MCAssocALambda", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dPtVsCentralityBeforeSel_MCAssocSigma0", "h2dPtVsCentralityBeforeSel_MCAssocSigma0", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dPtVsCentralityBeforeSel_MCAssocASigma0", "h2dPtVsCentralityBeforeSel_MCAssocASigma0", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dSigma0PtVsLambdaPtBeforeSel_MCAssoc", "h2dSigma0PtVsLambdaPtBeforeSel_MCAssoc", kTH2D, {axisPt, axisPt});
      histos.add("MC/h2dSigma0PtVsGammaPtBeforeSel_MCAssoc", "h2dSigma0PtVsGammaPtBeforeSel_MCAssoc", kTH2D, {axisPt, axisPt});
      histos.add("MC/h2dPtVsCentralityAfterSel_MCAssocSigma0", "h2dPtVsCentralityAfterSel_MCAssocSigma0", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dPtVsCentralityAfterSel_MCAssocASigma0", "h2dPtVsCentralityAfterSel_MCAssocASigma0", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dGammaXYConversion", "h2dGammaXYConversion", kTH2F, {axisXY, axisXY});
    }

    if (doprocessGeneratedRun3) {
      
      histos.add("Gen/hGenGamma", "hGenGamma", kTH1D, {axisPt});
      histos.add("Gen/hGenLambda", "hGenLambda", kTH1D, {axisPt});
      histos.add("Gen/hGenAntiLambda", "hGenAntiLambda", kTH1D, {axisPt});
      histos.add("Gen/hGenSigma0", "hGenSigma0", kTH1D, {axisPt});
      histos.add("Gen/hGenAntiSigma0", "hGenAntiSigma0", kTH1D, {axisPt});
      histos.add("Gen/hGenPi0", "hGenPi0", kTH1D, {axisPt});

      auto hGenSpecies = histos.add<TH1>("Gen/hGenSpecies", "hGenSpecies", kTH1D, {{4, -0.5f, 3.5f}});
      hGenSpecies->GetXaxis()->SetBinLabel(1, "Primary Lambda");
      hGenSpecies->GetXaxis()->SetBinLabel(2, "Primary ALambda");
      hGenSpecies->GetXaxis()->SetBinLabel(3, "Sigma0");
      hGenSpecies->GetXaxis()->SetBinLabel(4, "ASigma0");
      
      auto hPrimaryPi0s = histos.add<TH1>("Gen/hPrimaryPi0s", "hPrimaryPi0s", kTH1D, {{2, -0.5f, 1.5f}});
      hPrimaryPi0s->GetXaxis()->SetBinLabel(1, "All Pi0s");
      hPrimaryPi0s->GetXaxis()->SetBinLabel(2, "Primary Pi0s");

      auto h2DGenSigma0TypeVsProcess = histos.add<TH2>("Gen/h2DGenSigma0TypeVsProcess", "h2DGenSigma0TypeVsProcess", kTH2D, {{4, -0.5f, 3.5f}, {50, -0.5f, 49.5f}});
      h2DGenSigma0TypeVsProcess->GetXaxis()->SetBinLabel(1, "All Sigma0s");
      h2DGenSigma0TypeVsProcess->GetXaxis()->SetBinLabel(2, "Sterile");
      h2DGenSigma0TypeVsProcess->GetXaxis()->SetBinLabel(3, "Lambda+Gamma");
      h2DGenSigma0TypeVsProcess->GetXaxis()->SetBinLabel(4, "Others");
    }
  }

  // ______________________________________________________
  // Helper struct to store sigma0 reco and MC properties
  struct {
    float Mass;
    float pT;
    float Y;
    float SigmaOPAngle;
    bool fIsSigma;
    bool fIsAntiSigma;
    bool fIsPhotonCorrectlyAssign;
    bool fIsLambdaCorrectlyAssign;
    bool fIsPhotonPrimary;
    bool fIsLambdaPrimary;
    int PhotonCandPDGCode;
    int PhotonCandPDGCodeMother;
    int LambdaCandPDGCode;
    int LambdaCandPDGCodeMother;         
    float SigmaMCpT;
    float LambdaMCpT;
    float PhotonMCpT;

  } sigmaCandidate;
  
  // ______________________________________________________
  // MC-specific
  // Analyze v0-to-collision association
  template <typename TCollision, typename TV0Object>
  void analyzeV0CollAssoc(TCollision const& collision, TV0Object const& fullv0s, std::vector<int> selV0Indices, bool isPhotonAnalysis)
  {
    auto v0MCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
    float IR = (fGetIR) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource, fIRCrashOnNull) * 1.e-3 : -1;

    for (size_t i = 0; i < selV0Indices.size(); ++i) {
      auto v0 = fullv0s.rawIteratorAt(selV0Indices[i]);
      auto v0MC = v0.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

      float V0MCpT = RecoDecay::pt(array<float, 2>{v0MC.pxMC(), v0MC.pyMC()});
      float V0PA = TMath::ACos(v0.v0cosPA());
      bool fIsV0CorrectlyAssigned = (v0MC.straMCCollisionId() == v0MCCollision.globalIndex());
      bool isPrimary = v0MC.isPhysicalPrimary();

      if ((v0MC.pdgCode() == 22) && isPhotonAnalysis && isPrimary) { // True Gamma
        histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueGamma"), IR, V0MCpT);
        histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueGamma"), V0PA, IR, V0MCpT);

        if (!fIsV0CorrectlyAssigned) {
          histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueGamma_BadCollAssig"), IR, V0MCpT);
          histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueGamma_BadCollAssig"), V0PA, IR, V0MCpT);
        }
      }
      if ((v0MC.pdgCode() == 3122) && !isPhotonAnalysis && isPrimary) { // True Lambda
        histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueLambda"), IR, V0MCpT);
        histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueLambda"), V0PA, IR, V0MCpT);

        if (!fIsV0CorrectlyAssigned) {
          histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueLambda_BadCollAssig"), IR, V0MCpT);
          histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueLambda_BadCollAssig"), V0PA, IR, V0MCpT);
        }
      }
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)  
  template <typename TMCParticles>
  void fillGeneratedTable(TMCParticles const& mcParticles)
  {    
    for (auto& mcParticle : mcParticles) {
      if (TMath::Abs(mcParticle.y()) > GenMaxRap)
        continue;

      // Calculating properties               
      float ptmc = mcParticle.pt();                
      bool isV0Photon = mcParticle.pdgCode() == 22;
      bool isV0Lambda = mcParticle.pdgCode() == 3122;
      bool isV0AntiLambda = mcParticle.pdgCode() == -3122;
      bool isSigma0 = mcParticle.pdgCode() == 3212;
      bool isAntiSigma0 = mcParticle.pdgCode() == -3212;
      bool isPi0 = mcParticle.pdgCode() == 111;        
      bool isPrimary = mcParticle.isPhysicalPrimary();
      int sigma0Type = 0;
      int mccollisionid = -1;

      if (mcParticle.has_mcCollision()) 
        mccollisionid = mcParticle.mcCollisionId(); // save this reference, please
      
      if (isV0Photon && isPrimary) {                  
        histos.fill(HIST("Gen/hGenGamma"), ptmc);
      }
      if (isV0Lambda && isPrimary) {        
        histos.fill(HIST("Gen/hGenLambda"), ptmc);
      }
      if (isV0AntiLambda && isPrimary) {        
        histos.fill(HIST("Gen/hGenAntiLambda"), ptmc);
      }                                              
      if (isPi0){        
        histos.fill(HIST("Gen/hGenPi0"), ptmc);
        histos.fill(HIST("Gen/hPrimaryPi0s"), 0);
        if (isPrimary) histos.fill(HIST("Gen/hPrimaryPi0s"), 1);
        
        pi0Gens(ptmc); // optional table to store generated pi0 candidates. Be careful, this is a large table!
        pi0GenCollRefs(mccollisionid); // link to stramccollision table
      }  
      
      // Sigma0-specific
      if (isSigma0 || isAntiSigma0){     
        
        // Checking decay mode
        auto daughtersIDs = mcParticle.daughtersIds();
        auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
                
        if (daughters.size() == 2) {
          bool hasPhoton = false;
          bool hasLambdaOrAntiLambda = false;

          for (auto& daughter : daughters) {
            int daupdg = daughter.pdgCode();
            if (daupdg == 22)
              hasPhoton = true;
            if (TMath::Abs(daupdg) == 3122)
              hasLambdaOrAntiLambda = true;
          }

          if (hasPhoton && hasLambdaOrAntiLambda)
            sigma0Type = 1;
          else
            sigma0Type = 2;

        } else if (daughters.size() > 0) {
          sigma0Type = 2;
        }

        // Fill QA histograms and tables
        if (isSigma0)          
          histos.fill(HIST("Gen/hGenSigma0"), ptmc);        
        if (isAntiSigma0)
          histos.fill(HIST("Gen/hGenAntiSigma0"), ptmc);                    
              
        histos.fill(HIST("Gen/h2DGenSigma0TypeVsProcess"), 0, mcParticle.getProcess()); 
        if (sigma0Type==0)          
          histos.fill(HIST("Gen/h2DGenSigma0TypeVsProcess"), 1, mcParticle.getProcess());         
        if (sigma0Type==1)
          histos.fill(HIST("Gen/h2DGenSigma0TypeVsProcess"), 2, mcParticle.getProcess());         
        if (sigma0Type==2)          
          histos.fill(HIST("Gen/h2DGenSigma0TypeVsProcess"), 3, mcParticle.getProcess());         

        if (fillSigma0Tables){
          sigma0Gens(isSigma0, ptmc, sigma0Type); 
          sigma0GenCollRefs(mccollisionid); // link to stramccollision table
        }                        
      }
      
      if (isV0Lambda && isPrimary) histos.fill(HIST("Gen/hGenSpecies"), 0);
      if (isV0AntiLambda && isPrimary) histos.fill(HIST("Gen/hGenSpecies"), 1);
      if (isSigma0) histos.fill(HIST("Gen/hGenSpecies"), 2);
      if (isAntiSigma0) histos.fill(HIST("Gen/hGenSpecies"), 3);
    }      
  }

  //_______________________________________________
  // Process photon candidate
  template <typename TV0Object, typename TCollision>
  bool processPhotonCandidate(TV0Object const& gamma, TCollision const& collision)
  {
    // V0 type selection
    if (gamma.v0Type() == 0)
      return false;

    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    //_______________________________________________
    // MC Processing
    if constexpr (requires { gamma.motherMCPartId();}) {
      if (!gamma.has_v0MCCore())
        return false;

      auto gammaMC = gamma.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();      

      if (gammaMC.pdgCode() == 22) {
        histos.fill(HIST("MC/h2dGammaXYConversion"), gamma.x(), gamma.y());
        float GammaY = TMath::Abs(RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma));
        if (GammaY <= 0.5) {                                                                                                                // rapidity selection
          histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocGamma"), centrality, gamma.pt());                                            // isgamma
        }
      }
    }
        
    if (useMLScores) {      
      if (gamma.gammaBDTScore() <= Gamma_MLThreshold)
        return false;

    } else {
      // Standard selection
      // Gamma basic selection criteria:
      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 1.);
      histos.fill(HIST("PhotonSel/hPhotonMass"), gamma.mGamma());
      if ((gamma.mGamma() < 0) || (gamma.mGamma() > PhotonMaxMass))
        return false;
      float PhotonY = RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma);
      histos.fill(HIST("PhotonSel/hPhotonNegEta"), gamma.negativeeta());
      histos.fill(HIST("PhotonSel/hPhotonPosEta"), gamma.positiveeta());
      histos.fill(HIST("PhotonSel/hPhotonY"), PhotonY);
      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 2.);
      if ((TMath::Abs(PhotonY) > V0Rapidity) || (TMath::Abs(gamma.negativeeta()) > PhotonMaxDauPseudoRap) || (TMath::Abs(gamma.positiveeta()) > PhotonMaxDauPseudoRap))
        return false;
      histos.fill(HIST("PhotonSel/hPhotonDCANegToPV"), TMath::Abs(gamma.dcanegtopv()));
      histos.fill(HIST("PhotonSel/hPhotonDCAPosToPV"), TMath::Abs(gamma.dcapostopv()));
      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 3.);
      if ((TMath::Abs(gamma.dcapostopv()) < PhotonMinDCAToPv) || (TMath::Abs(gamma.dcanegtopv()) < PhotonMinDCAToPv))
        return false;
      histos.fill(HIST("PhotonSel/hPhotonDCADau"), TMath::Abs(gamma.dcaV0daughters()));
      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 4.);
      if (TMath::Abs(gamma.dcaV0daughters()) > PhotonMaxDCAV0Dau)
        return false;
      histos.fill(HIST("PhotonSel/hPhotonRadius"), gamma.v0radius());
      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 5.);
      if ((gamma.v0radius() < PhotonMinRadius) || (gamma.v0radius() > PhotonMaxRadius))
        return false;
      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 6.);
    }
    
    histos.fill(HIST("PhotonSel/h3dPhotonMass"), centrality, gamma.pt(), gamma.mGamma());
    return true;
  }

  //_______________________________________________
  // Process lambda candidate
  template <typename TV0Object, typename TCollision>
  bool processLambdaCandidate(TV0Object const& lambda, TCollision const& collision)
  {    
    // V0 type selection
    if (lambda.v0Type() != 1)
      return false;

    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    
    //_______________________________________________
    // MC Processing
    if constexpr (requires { lambda.motherMCPartId();}) {
      if (!lambda.has_v0MCCore())
        return false;
      
      auto lambdaMC = lambda.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
      
      if (TMath::Abs(lambda.yLambda()) <= 0.5) {
        if (lambdaMC.pdgCode() == 3122) // Is Lambda
          histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocLambda"), centrality, lambda.pt());
        if (lambdaMC.pdgCode() == -3122) // Is AntiLambda
          histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocALambda"), centrality, lambda.pt());
      }
    }
        
    if (useMLScores) {
      if ((lambda.lambdaBDTScore() <= Lambda_MLThreshold) && (lambda.antiLambdaBDTScore() <= AntiLambda_MLThreshold))
        return false;

    } else {
      // Lambda basic selection criteria:
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 1.);
      histos.fill(HIST("LambdaSel/hLambdaMass"), lambda.mLambda());
      histos.fill(HIST("LambdaSel/hAntiLambdaMass"), lambda.mAntiLambda());
      if ((TMath::Abs(lambda.mLambda() - o2::constants::physics::MassLambda0) > LambdaWindow) && (TMath::Abs(lambda.mAntiLambda() - o2::constants::physics::MassLambda0) > LambdaWindow))
        return false;
      histos.fill(HIST("LambdaSel/hLambdaNegEta"), lambda.negativeeta());
      histos.fill(HIST("LambdaSel/hLambdaPosEta"), lambda.positiveeta());
      histos.fill(HIST("LambdaSel/hLambdaY"), lambda.yLambda());
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 2.);
      if ((TMath::Abs(lambda.yLambda()) > V0Rapidity) || (TMath::Abs(lambda.negativeeta()) > LambdaDauPseudoRap) || (TMath::Abs(lambda.positiveeta()) > LambdaDauPseudoRap))
        return false;
      histos.fill(HIST("LambdaSel/hLambdaDCANegToPV"), lambda.dcanegtopv());
      histos.fill(HIST("LambdaSel/hLambdaDCAPosToPV"), lambda.dcapostopv());
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 3.);
      if ((TMath::Abs(lambda.dcapostopv()) < LambdaMinDCAPosToPv) || (TMath::Abs(lambda.dcanegtopv()) < LambdaMinDCANegToPv))
        return false;
      histos.fill(HIST("LambdaSel/hLambdaRadius"), lambda.v0radius());
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 4.);
      if ((lambda.v0radius() < LambdaMinv0radius) || (lambda.v0radius() > LambdaMaxv0radius))
        return false;
      histos.fill(HIST("LambdaSel/hLambdaDCADau"), lambda.dcaV0daughters());
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 5.);
      if (TMath::Abs(lambda.dcaV0daughters()) > LambdaMaxDCAV0Dau)
        return false;
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 6.);
    }
    
    histos.fill(HIST("LambdaSel/h3dLambdaMass"), centrality, lambda.pt(), lambda.mLambda());
    histos.fill(HIST("LambdaSel/h3dALambdaMass"), centrality, lambda.pt(), lambda.mAntiLambda());

    return true;
  }

  //_______________________________________________
  // Build pi0 candidate for QA
  template <typename TV0Object, typename TCollision>
  bool buildPi0(TV0Object const& gamma1, TV0Object const& gamma2, TCollision collision)
  {
    //_______________________________________________
    // Check if both V0s are made of the same tracks
    if (gamma1.posTrackExtraId() == gamma2.posTrackExtraId() ||
        gamma1.negTrackExtraId() == gamma2.negTrackExtraId()) {
      return false;
    }

    //_______________________________________________
    // Calculate pi0 properties
    std::array<float, 3> pVecGamma1{gamma1.px(), gamma1.py(), gamma1.pz()};
    std::array<float, 3> pVecGamma2{gamma2.px(), gamma2.py(), gamma2.pz()};
    std::array arrpi0{pVecGamma1, pVecGamma2};
    float pi0Mass = RecoDecay::m(arrpi0, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassPhoton});
    float pi0Pt = RecoDecay::pt(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py()});
    float pi0Y = RecoDecay::y(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py(), gamma1.pz() + gamma2.pz()}, o2::constants::physics::MassPi0);    

    //_______________________________________________
    // Pi0-specific selections:
    if (TMath::Abs(pi0Y) > Pi0MaxRap) 
      return false;

    if (TMath::Abs(pi0Mass - o2::constants::physics::MassPi0) > Pi0MassWindow)
      return false;
    
    // Fill optional tables for QA
    // Define the table!    
    auto posTrackGamma1 = gamma1.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma1 = gamma1.template negTrackExtra_as<dauTracks>();
    auto posTrackGamma2 = gamma2.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma2 = gamma2.template negTrackExtra_as<dauTracks>();

    float fPhoton1Y = RecoDecay::y(std::array{gamma1.px(), gamma1.py(), gamma1.pz()}, o2::constants::physics::MassGamma);
    float fPhoton2Y = RecoDecay::y(std::array{gamma2.px(), gamma2.py(), gamma2.pz()}, o2::constants::physics::MassGamma);

    // Check if MC data and populate corresponding table
    if constexpr (requires { gamma1.motherMCPartId(); gamma2.motherMCPartId(); }) {
      bool fIsPi0 = false;
      bool fhasMCColl = false;
      bool fIsPhoton1CorrectlyAssign = false;
      bool fIsPhoton2CorrectlyAssign = false;
      
      if (collision.has_straMCCollision())
        fhasMCColl = true;

      if (!gamma1.has_v0MCCore() || !gamma2.has_v0MCCore())
        return false;
      
      auto gamma1MC = gamma1.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
      auto gamma2MC = gamma2.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

      if (gamma1MC.pdgCode() == 22 && gamma2MC.pdgCode() == 22 &&
          gamma1MC.pdgCodeMother() == 111 && gamma2MC.pdgCodeMother() == 111 &&
          gamma1.motherMCPartId() == gamma2.motherMCPartId()) {
        fIsPi0 = true;        
      }

      if (fhasMCColl) {
        auto MCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();        
        fIsPhoton1CorrectlyAssign = (gamma1MC.straMCCollisionId() == MCCollision.globalIndex());      
        fIsPhoton2CorrectlyAssign = (gamma2MC.straMCCollisionId() == MCCollision.globalIndex());
      }

      bool fIsPhoton1Primary = gamma1MC.isPhysicalPrimary();
      bool fIsPhoton2Primary = gamma2MC.isPhysicalPrimary();

      int Photon1PDGCode = gamma1MC.pdgCode();
      int Photon1PDGCodeMother = gamma1MC.pdgCodeMother();
      int Photon2PDGCode = gamma2MC.pdgCode();
      int Photon2PDGCodeMother = gamma2MC.pdgCodeMother();
      
      float pi0MCpT = RecoDecay::pt(array{gamma1MC.pxMC() + gamma2MC.pxMC(), gamma1MC.pyMC() + gamma2MC.pyMC()});
      float Photon1MCpT = RecoDecay::pt(array{gamma1MC.pxMC(), gamma1MC.pyMC()});
      float Photon2MCpT = RecoDecay::pt(array{gamma2MC.pxMC(), gamma2MC.pyMC()});

      // Fill table with MC info
      pi0coresmc(fIsPi0, pi0MCpT, 
                    Photon1MCpT, Photon1PDGCode, Photon1PDGCodeMother, fIsPhoton1Primary, fIsPhoton1CorrectlyAssign,
                    Photon2MCpT, Photon2PDGCode, Photon2PDGCodeMother, fIsPhoton2Primary, fIsPhoton2CorrectlyAssign);             
    }
    
    pi0cores(pi0Pt, pi0Mass, pi0Y,
                gamma1.mGamma(), gamma1.pt(), gamma1.qtarm(), gamma1.alpha(), gamma1.dcapostopv(), gamma1.dcanegtopv(), gamma1.dcaV0daughters(), gamma1.negativeeta(), gamma1.positiveeta(), 
                gamma1.v0cosPA(), gamma1.v0radius(), gamma1.z(), posTrackGamma1.tpcCrossedRows(), negTrackGamma1.tpcCrossedRows(), posTrackGamma1.tpcNSigmaEl(), negTrackGamma1.tpcNSigmaEl(), fPhoton1Y, gamma1.v0Type(), 
                gamma2.mGamma(), gamma2.pt(), gamma2.qtarm(), gamma2.alpha(), gamma2.dcapostopv(), gamma2.dcanegtopv(), gamma2.dcaV0daughters(), gamma2.negativeeta(), gamma2.positiveeta(), 
                gamma2.v0cosPA(), gamma2.v0radius(), gamma2.z(), posTrackGamma2.tpcCrossedRows(), negTrackGamma2.tpcCrossedRows(), posTrackGamma2.tpcNSigmaEl(), negTrackGamma2.tpcNSigmaEl(), fPhoton2Y, gamma2.v0Type()); 
  
    pi0coresRefs(collision.globalIndex());   

    return true;                     
  }

  //_______________________________________________
  // Build sigma0 candidate 
  template <typename TV0Object, typename TCollision>
  bool buildSigma0(TV0Object const& lambda, TV0Object const& gamma, TCollision const& collision)
  {
    //_______________________________________________
    // Initial setup
    // Checking if both V0s are made of the very same tracks
    if (gamma.posTrackExtraId() == lambda.posTrackExtraId() ||
        gamma.negTrackExtraId() == lambda.negTrackExtraId()) {
      return false;
    }

    // Sigma0 candidate properties
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};

    TVector3 v1(gamma.px(), gamma.py(), gamma.pz());
    TVector3 v2(lambda.px(), lambda.py(), lambda.pz());

    auto arrMom = std::array{pVecPhotons, pVecLambda};
    float sigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
    float sigmaY = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);
    float fSigmaOPAngle = v1.Angle(v2);
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    //_______________________________________________
    // MC-specific
    bool fIsSigma = false;
    bool fIsAntiSigma = false;
    bool fhasMCColl = false;
    bool fIsPhotonCorrectlyAssign = false;
    bool fIsLambdaCorrectlyAssign = false;
    bool fIsPhotonPrimary = false;
    bool fIsLambdaPrimary = false;
    int PhotonCandPDGCode = -1;
    int PhotonCandPDGCodeMother = -1;
    int LambdaCandPDGCode = -1;
    int LambdaCandPDGCodeMother = -1;
    float SigmaMCpT = -1;
    float LambdaMCpT = -1;
    float PhotonMCpT = -1;

    if constexpr (requires { gamma.motherMCPartId(); lambda.motherMCPartId(); }) {

      if (collision.has_straMCCollision())
        fhasMCColl = true;

      if (!gamma.has_v0MCCore() || !lambda.has_v0MCCore())
        return false;
      
      auto gammaMC = gamma.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
      auto lambdaMC = lambda.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

      if (fhasMCColl) {
        auto gammaMCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
        auto lambdaMCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
        fIsPhotonCorrectlyAssign = (gammaMC.straMCCollisionId() == gammaMCCollision.globalIndex());      
        fIsLambdaCorrectlyAssign = (lambdaMC.straMCCollisionId() == lambdaMCCollision.globalIndex());
      }

      fIsPhotonPrimary = gammaMC.isPhysicalPrimary();
      fIsLambdaPrimary = lambdaMC.isPhysicalPrimary();

      PhotonCandPDGCode = gammaMC.pdgCode();
      PhotonCandPDGCodeMother = gammaMC.pdgCodeMother();
      LambdaCandPDGCode = lambdaMC.pdgCode();
      LambdaCandPDGCodeMother = lambdaMC.pdgCodeMother();
      SigmaMCpT = RecoDecay::pt(array{gammaMC.pxMC() + lambdaMC.pxMC(), gammaMC.pyMC() + lambdaMC.pyMC()});
      LambdaMCpT = RecoDecay::pt(array{lambdaMC.pxMC(), lambdaMC.pyMC()});
      PhotonMCpT = RecoDecay::pt(array{gammaMC.pxMC(), gammaMC.pyMC()});

      if ((PhotonCandPDGCode == 22) && (PhotonCandPDGCodeMother == 3212) && (LambdaCandPDGCode == 3122) && (LambdaCandPDGCodeMother == 3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()))
        fIsSigma = true;
      if ((PhotonCandPDGCode == 22) && (PhotonCandPDGCodeMother == -3212) && (LambdaCandPDGCode == -3122) && (LambdaCandPDGCodeMother == -3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()))
        fIsAntiSigma = true;

      if (TMath::Abs(sigmaY) <= 0.5) {
        if (fIsSigma) {
          histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocSigma0"), centrality, SigmaMCpT);
          histos.fill(HIST("MC/h2dSigma0PtVsLambdaPtBeforeSel_MCAssoc"), SigmaMCpT, LambdaMCpT);
          histos.fill(HIST("MC/h2dSigma0PtVsGammaPtBeforeSel_MCAssoc"), SigmaMCpT, PhotonMCpT);
        }
        if (fIsAntiSigma)
          histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocASigma0"), centrality, SigmaMCpT);
      }
    }
      
    //_______________________________________________
    // Before any selection
    histos.fill(HIST("SigmaSel/h3dMassSigma0BeforeSel"), centrality, SigmapT, sigmaMass);
    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 1.);
    histos.fill(HIST("SigmaSel/hSigmaMass"), sigmaMass);
    histos.fill(HIST("SigmaSel/hSigmaMassWindow"), sigmaMass - o2::constants::physics::MassSigma0);

    if (fillQAhistos) {
      histos.fill(HIST("GeneralQA/h2dMassGammaVsK0S"), gamma.mGamma(), gamma.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassLambdaVsK0S"), lambda.mLambda(), lambda.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassGammaVsLambda"), gamma.mGamma(), gamma.mLambda());
      histos.fill(HIST("GeneralQA/h2dMassLambdaVsGamma"), lambda.mLambda(), lambda.mGamma());
      histos.fill(HIST("GeneralQA/h3dMassSigma0VsDaupTs"), gamma.pt(), lambda.pt(), sigmaMass);
    }

    //_______________________________________________
    // Sigma0 selections
    if (TMath::Abs(sigmaMass - o2::constants::physics::MassSigma0) > Sigma0Window)
      return false;

    histos.fill(HIST("SigmaSel/hSigmaY"), sigmaY);
    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 2.);

    if (TMath::Abs(sigmaY) > SigmaMaxRap)
      return false;

    histos.fill(HIST("SigmaSel/hSigmaMassSelected"), sigmaMass);
    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 3.);

    //_______________________________________________
    // After selections
    if (fillQAhistos) {
      histos.fill(HIST("GeneralQA/h2dMassGammaVsK0SAfterMassSel"), gamma.mGamma(), gamma.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassLambdaVsK0SAfterMassSel"), lambda.mLambda(), lambda.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassGammaVsLambdaAfterMassSel"), gamma.mGamma(), lambda.mLambda());
      histos.fill(HIST("GeneralQA/h2dV0XY"), gamma.x(), gamma.y());
    }
    
    // MC-specific
    if constexpr (requires { gamma.motherMCPartId();}) {
      if (TMath::Abs(sigmaY) <= 0.5) {
        if (fIsSigma)
          histos.fill(HIST("MC/h2dPtVsCentralityAfterSel_MCAssocSigma0"), centrality, SigmaMCpT);
        if (fIsAntiSigma)
          histos.fill(HIST("MC/h2dPtVsCentralityAfterSel_MCAssocASigma0"), centrality, SigmaMCpT);              
      }      
    }
    
    histos.fill(HIST("SigmaSel/h3dMassSigma0AfterSel"), centrality, SigmapT, sigmaMass);

    // Store Reco properties in struct
    sigmaCandidate.Mass = sigmaMass;
    sigmaCandidate.pT = SigmapT;
    sigmaCandidate.Y = sigmaY;
    sigmaCandidate.SigmaOPAngle = fSigmaOPAngle;

    // Store MC properties if available
    sigmaCandidate.fIsSigma = fIsSigma;
    sigmaCandidate.fIsAntiSigma = fIsAntiSigma;
    sigmaCandidate.fIsPhotonCorrectlyAssign = fIsPhotonCorrectlyAssign;
    sigmaCandidate.fIsLambdaCorrectlyAssign = fIsLambdaCorrectlyAssign;
    sigmaCandidate.fIsPhotonPrimary = fIsPhotonPrimary;
    sigmaCandidate.fIsLambdaPrimary = fIsLambdaPrimary;
    sigmaCandidate.PhotonCandPDGCode = PhotonCandPDGCode;
    sigmaCandidate.PhotonCandPDGCodeMother = PhotonCandPDGCodeMother;
    sigmaCandidate.LambdaCandPDGCode = LambdaCandPDGCode;
    sigmaCandidate.LambdaCandPDGCodeMother = LambdaCandPDGCodeMother;    
    sigmaCandidate.SigmaMCpT = SigmaMCpT;
    sigmaCandidate.LambdaMCpT = LambdaMCpT;
    sigmaCandidate.PhotonMCpT = PhotonMCpT;

    return true;
  }

  //_______________________________________________
  // Fill tables with reconstructed sigma0 candidate
  template <typename TV0Object, typename TCollision>
  void buildSigma0Tables(TV0Object const& lambda, TV0Object const& gamma, TCollision const& coll)
  {
    //_______________________________________________    
    // Photon properties
    auto posTrackGamma = gamma.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma = gamma.template negTrackExtra_as<dauTracks>();

    float fPhotonPt = gamma.pt();
    float fPhotonMass = gamma.mGamma();
    float fPhotonQt = gamma.qtarm();
    float fPhotonAlpha = gamma.alpha();
    float fPhotonRadius = gamma.v0radius();
    float fPhotonCosPA = gamma.v0cosPA();
    float fPhotonDCADau = gamma.dcaV0daughters();
    float fPhotonDCANegPV = gamma.dcanegtopv();
    float fPhotonDCAPosPV = gamma.dcapostopv();
    float fPhotonZconv = gamma.z();
    float fPhotonEta = gamma.eta();
    float fPhotonY = RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma);
    float fPhotonPhi = RecoDecay::phi(gamma.px(), gamma.py());
    float fPhotonPosTPCNSigmaEl = posTrackGamma.tpcNSigmaEl();
    float fPhotonNegTPCNSigmaEl = negTrackGamma.tpcNSigmaEl();
    float fPhotonPosTPCNSigmaPi = posTrackGamma.tpcNSigmaPi();
    float fPhotonNegTPCNSigmaPi = negTrackGamma.tpcNSigmaPi();
    uint8_t fPhotonPosTPCCrossedRows = posTrackGamma.tpcCrossedRows();
    uint8_t fPhotonNegTPCCrossedRows = negTrackGamma.tpcCrossedRows();
    float fPhotonPosPt = gamma.positivept();
    float fPhotonNegPt = gamma.negativept();
    float fPhotonPosEta = gamma.positiveeta();
    float fPhotonNegEta = gamma.negativeeta();
    float fPhotonPosY = RecoDecay::y(std::array{gamma.pxpos(), gamma.pypos(), gamma.pzpos()}, o2::constants::physics::MassElectron);
    float fPhotonNegY = RecoDecay::y(std::array{gamma.pxneg(), gamma.pyneg(), gamma.pzneg()}, o2::constants::physics::MassElectron);
    float fPhotonPsiPair = gamma.psipair();
    int fPhotonPosITSCls = posTrackGamma.itsNCls();
    int fPhotonNegITSCls = negTrackGamma.itsNCls();
    float fPhotonPosITSChi2PerNcl = posTrackGamma.itsChi2PerNcl();
    float fPhotonNegITSChi2PerNcl = negTrackGamma.itsChi2PerNcl();
    uint8_t fPhotonV0Type = gamma.v0Type();

    uint8_t fPhotonPosTrackCode = ((uint8_t(posTrackGamma.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrackGamma.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrackGamma.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrackGamma.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrackGamma.hasTOF()) << hasTOF));

    uint8_t fPhotonNegTrackCode = ((uint8_t(negTrackGamma.hasTPC()) << hasTPC) |
                                   (uint8_t(negTrackGamma.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(negTrackGamma.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(negTrackGamma.hasTRD()) << hasTRD) |
                                   (uint8_t(negTrackGamma.hasTOF()) << hasTOF));

    float GammaBDTScore = gamma.gammaBDTScore();

    //_______________________________________________    
    // Lambda properties
    auto posTrackLambda = lambda.template posTrackExtra_as<dauTracks>();
    auto negTrackLambda = lambda.template negTrackExtra_as<dauTracks>();

    float fLambdaPt = lambda.pt();
    float fLambdaMass = lambda.mLambda();
    float fAntiLambdaMass = lambda.mAntiLambda();
    float fLambdaQt = lambda.qtarm();
    float fLambdaAlpha = lambda.alpha();
    float fLambdaLifeTime = lambda.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * o2::constants::physics::MassLambda0;
    float fLambdaRadius = lambda.v0radius();
    float fLambdaCosPA = lambda.v0cosPA();
    float fLambdaDCADau = lambda.dcaV0daughters();
    float fLambdaDCANegPV = lambda.dcanegtopv();
    float fLambdaDCAPosPV = lambda.dcapostopv();
    float fLambdaEta = lambda.eta();
    float fLambdaY = lambda.yLambda();
    float fLambdaPhi = RecoDecay::phi(lambda.px(), lambda.py());
    float fLambdaPosPrTPCNSigma = posTrackLambda.tpcNSigmaPr();
    float fLambdaPosPiTPCNSigma = posTrackLambda.tpcNSigmaPi();
    float fLambdaNegPrTPCNSigma = negTrackLambda.tpcNSigmaPr();
    float fLambdaNegPiTPCNSigma = negTrackLambda.tpcNSigmaPi();

    float fLambdaPrTOFNSigma = -999.f;
    float fLambdaPiTOFNSigma = -999.f;
    float fALambdaPrTOFNSigma = -999.f;
    float fALambdaPiTOFNSigma = -999.f;

    if constexpr (requires { lambda.tofNSigmaLaPr();}) { // If TOF info avaiable
      fLambdaPrTOFNSigma = lambda.tofNSigmaLaPr();
      fLambdaPiTOFNSigma = lambda.tofNSigmaLaPi();
      fALambdaPrTOFNSigma = lambda.tofNSigmaALaPr();
      fALambdaPiTOFNSigma = lambda.tofNSigmaALaPi();
    }    

    uint8_t fLambdaPosTPCCrossedRows = posTrackLambda.tpcCrossedRows();
    uint8_t fLambdaNegTPCCrossedRows = negTrackLambda.tpcCrossedRows();
    float fLambdaPosPt = lambda.positivept();
    float fLambdaNegPt = lambda.negativept();
    float fLambdaPosEta = lambda.positiveeta();
    float fLambdaNegEta = lambda.negativeeta();
    float fLambdaPosPrY = RecoDecay::y(std::array{lambda.pxpos(), lambda.pypos(), lambda.pzpos()}, o2::constants::physics::MassProton);
    float fLambdaPosPiY = RecoDecay::y(std::array{lambda.pxpos(), lambda.pypos(), lambda.pzpos()}, o2::constants::physics::MassPionCharged);
    float fLambdaNegPrY = RecoDecay::y(std::array{lambda.pxneg(), lambda.pyneg(), lambda.pzneg()}, o2::constants::physics::MassProton);
    float fLambdaNegPiY = RecoDecay::y(std::array{lambda.pxneg(), lambda.pyneg(), lambda.pzneg()}, o2::constants::physics::MassPionCharged);
    int fLambdaPosITSCls = posTrackLambda.itsNCls();
    int fLambdaNegITSCls = negTrackLambda.itsNCls();
    float fLambdaPosITSChi2PerNcl = posTrackLambda.itsChi2PerNcl();
    float fLambdaNegITSChi2PerNcl = negTrackLambda.itsChi2PerNcl();
    uint8_t fLambdaV0Type = lambda.v0Type();

    uint8_t fLambdaPosTrackCode = ((uint8_t(posTrackLambda.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrackLambda.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrackLambda.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrackLambda.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrackLambda.hasTOF()) << hasTOF));

    uint8_t fLambdaNegTrackCode = ((uint8_t(negTrackLambda.hasTPC()) << hasTPC) |
                                   (uint8_t(negTrackLambda.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(negTrackLambda.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(negTrackLambda.hasTRD()) << hasTRD) |
                                   (uint8_t(negTrackLambda.hasTOF()) << hasTOF));

    float LambdaBDTScore = lambda.lambdaBDTScore();
    float AntiLambdaBDTScore = lambda.antiLambdaBDTScore();                                   

    //_______________________________________________    
    // Fill the tables, please
    sigma0cores(sigmaCandidate.pT, sigmaCandidate.Mass, sigmaCandidate.Y, sigmaCandidate.SigmaOPAngle);

    sigmaPhotonExtras(fPhotonPt, fPhotonMass, fPhotonQt, fPhotonAlpha, fPhotonRadius,
                      fPhotonCosPA, fPhotonDCADau, fPhotonDCANegPV, fPhotonDCAPosPV, fPhotonZconv,
                      fPhotonEta, fPhotonY, fPhotonPhi, fPhotonPosTPCNSigmaEl, fPhotonNegTPCNSigmaEl, fPhotonPosTPCNSigmaPi, fPhotonNegTPCNSigmaPi, fPhotonPosTPCCrossedRows,
                      fPhotonNegTPCCrossedRows, fPhotonPosPt, fPhotonNegPt, fPhotonPosEta,
                      fPhotonNegEta, fPhotonPosY, fPhotonNegY, fPhotonPsiPair,
                      fPhotonPosITSCls, fPhotonNegITSCls, fPhotonPosITSChi2PerNcl, fPhotonNegITSChi2PerNcl, fPhotonPosTrackCode, fPhotonNegTrackCode,
                      fPhotonV0Type, GammaBDTScore);

    sigmaLambdaExtras(fLambdaPt, fLambdaMass, fAntiLambdaMass, fLambdaQt, fLambdaAlpha, fLambdaLifeTime,
                      fLambdaRadius, fLambdaCosPA, fLambdaDCADau, fLambdaDCANegPV,
                      fLambdaDCAPosPV, fLambdaEta, fLambdaY, fLambdaPhi, fLambdaPosPrTPCNSigma,
                      fLambdaPosPiTPCNSigma, fLambdaNegPrTPCNSigma, fLambdaNegPiTPCNSigma,
                      fLambdaPrTOFNSigma, fLambdaPiTOFNSigma, fALambdaPrTOFNSigma, fALambdaPiTOFNSigma,
                      fLambdaPosTPCCrossedRows, fLambdaNegTPCCrossedRows, fLambdaPosPt, fLambdaNegPt, fLambdaPosEta,
                      fLambdaNegEta, fLambdaPosPrY, fLambdaPosPiY, fLambdaNegPrY, fLambdaNegPiY,
                      fLambdaPosITSCls, fLambdaNegITSCls, fLambdaPosITSChi2PerNcl, fLambdaNegITSChi2PerNcl, fLambdaPosTrackCode, fLambdaNegTrackCode,
                      fLambdaV0Type, LambdaBDTScore, AntiLambdaBDTScore);

    if (doprocessMonteCarlo || doprocessMonteCarloWithTOF){

      sigma0mccores(sigmaCandidate.fIsSigma, sigmaCandidate.fIsAntiSigma, sigmaCandidate.SigmaMCpT,
                    sigmaCandidate.PhotonCandPDGCode, sigmaCandidate.PhotonCandPDGCodeMother, sigmaCandidate.fIsPhotonPrimary, sigmaCandidate.PhotonMCpT, sigmaCandidate.fIsPhotonCorrectlyAssign,
                    sigmaCandidate.LambdaCandPDGCode, sigmaCandidate.LambdaCandPDGCodeMother, sigmaCandidate.fIsLambdaPrimary, sigmaCandidate.LambdaMCpT, sigmaCandidate.fIsLambdaCorrectlyAssign);
      
    }

    sigma0CollRefs(coll.globalIndex());
       
  }

  // Process photon and lambda candidates to build sigma0 candidates
  template <typename TCollision, typename TV0s>
  void dataProcess(TCollision const& collisions, TV0s const& fullV0s)
  {
    //_______________________________________________
    // Initial setup
    // Auxiliary vectors to store best candidates
    std::vector<int> bestGammasArray;
    std::vector<int> bestLambdasArray;

    // Custom grouping
    std::vector<std::vector<int>> v0grouped(collisions.size());

    for (const auto& v0 : fullV0s) {
      v0grouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }

    //_______________________________________________
    // Collisions loop
    for (const auto& coll : collisions) {
      // Clear vectors
      bestGammasArray.clear();
      bestLambdasArray.clear();

      float centrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();
      histos.fill(HIST("hEventCentrality"), centrality);

      //_______________________________________________
      // V0s loop
      for (size_t i = 0; i < v0grouped[coll.globalIndex()].size(); i++) {
        auto v0 = fullV0s.rawIteratorAt(v0grouped[coll.globalIndex()][i]);

        if (processPhotonCandidate(v0, coll))          // selecting photons
          bestGammasArray.push_back(v0.globalIndex()); // Save indices of best gamma candidates

        if (processLambdaCandidate(v0, coll))           // selecting lambdas
          bestLambdasArray.push_back(v0.globalIndex()); // Save indices of best lambda candidates
      }

      //_______________________________________________      
      // Wrongly collision association study (MC-specific)
      if constexpr (requires { coll.StraMCCollisionId(); }) {
        if (doAssocStudy) {          
          analyzeV0CollAssoc(coll, fullV0s, bestGammasArray, true);   // Photon-analysis
          analyzeV0CollAssoc(coll, fullV0s, bestLambdasArray, false); // Lambda-analysis
        }
      }

      //_______________________________________________
      // V0 nested loop
      for (size_t i = 0; i < bestGammasArray.size(); ++i) {
        auto gamma1 = fullV0s.rawIteratorAt(bestGammasArray[i]);
                
        //_______________________________________________
        // Sigma0 loop
        if (fillSigma0Tables){
          for (size_t j = 0; j < bestLambdasArray.size(); ++j) {
            auto lambda = fullV0s.rawIteratorAt(bestLambdasArray[j]);
                            
            // Building sigma0 candidate
            if (!buildSigma0(lambda, gamma1, coll))
              continue;
          
            // Fill tables with accepted candidates           
            buildSigma0Tables(lambda, gamma1, coll);          
          }
        }
        
        //_______________________________________________
        // pi0 loop
        if (fillPi0Tables){
          for (size_t j = i + 1; j < bestGammasArray.size(); ++j) {
            auto gamma2 = fullV0s.rawIteratorAt(bestGammasArray[j]);
            if (!buildPi0(gamma1, gamma2, coll))
              continue;
          }
        }
      }
    }
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0StandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    dataProcess(collisions, fullV0s);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0DerivedMCDatas const& fullV0s, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    dataProcess(collisions, fullV0s);
  }

  void processRealDataWithTOF(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0TOFStandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    dataProcess(collisions, fullV0s);
  }

  void processMonteCarloWithTOF(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0TOFDerivedMCDatas const& fullV0s, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    dataProcess(collisions, fullV0s);
  }

  // Simulated processing in Run 3 - run this over original AO2Ds
  void processGeneratedRun3(aod::McParticles const& mcParticles)
  {
    fillGeneratedTable(mcParticles);
  }
  
  PROCESS_SWITCH(sigma0builder, processRealData, "process as if real data", true);
  PROCESS_SWITCH(sigma0builder, processRealDataWithTOF, "process as if real data", false);
  PROCESS_SWITCH(sigma0builder, processMonteCarlo, "process as if MC data", false);
  PROCESS_SWITCH(sigma0builder, processMonteCarloWithTOF, "process as if MC data, uses TOF PID info", false);
  PROCESS_SWITCH(sigma0builder, processGeneratedRun3, "process generated MC info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigma0builder>(cfgc)};
}
