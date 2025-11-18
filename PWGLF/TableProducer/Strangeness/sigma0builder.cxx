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

#include "PWGLF/DataModel/LFSigmaTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector3D.h"
#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

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
  Produces<aod::Sigma0Cores> sigma0cores;              // sigma0 candidates info for analysis
  Produces<aod::Sigma0PhotonExtras> sigmaPhotonExtras; // photons from sigma0 candidates info
  Produces<aod::Sigma0LambdaExtras> sigmaLambdaExtras; // lambdas from sigma0 candidates info
  Produces<aod::SigmaCollRef> sigma0CollRefs;          // references collisions from Sigma0Cores
  Produces<aod::Sigma0MCCores> sigma0mccores;          // Reco sigma0 MC properties
  Produces<aod::Sigma0Gens> sigma0Gens;                // Generated sigma0s
  Produces<aod::SigmaGenCollRef> sigma0GenCollRefs;    // references collisions from sigma0Gens

  //__________________________________________________
  // Pi0 specific
  Produces<aod::Pi0Cores> pi0cores;            // pi0 candidates info for analysis
  Produces<aod::Pi0CollRef> pi0coresRefs;      // references collisions from photonpair
  Produces<aod::Pi0CoresMC> pi0coresmc;        // Reco pi0 MC properties
  Produces<aod::Pi0Gens> pi0Gens;              // Generated pi0s
  Produces<aod::Pi0GenCollRef> pi0GenCollRefs; // references collisions from pi0Gens

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

  //// Photon criteria (for sigma0s and pi0s):
  Configurable<float> PhotonMaxDauPseudoRap{"PhotonMaxDauPseudoRap", 1.5, "Max pseudorapidity of daughter tracks"};
  Configurable<float> PhotonMinDCAToPv{"PhotonMinDCAToPv", 0.0, "Min DCA daughter To PV (cm)"};
  Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> PhotonMinRadius{"PhotonMinRadius", 0.0, "Min photon conversion radius (cm)"};
  Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 240, "Max photon conversion radius (cm)"};
  Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.3, "Max photon mass (GeV/c^{2})"};

  //// Sigma0 criteria:
  Configurable<float> Sigma0Window{"Sigma0Window", 0.1, "Mass window around expected (in GeV/c2)"};
  Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.8, "Max sigma0 rapidity"};

  //// Pi0 criteria::
  Configurable<float> Pi0MaxRap{"Pi0MaxRap", 0.8, "Max Pi0 Rapidity"};
  Configurable<float> Pi0MassWindow{"Pi0MassWindow", 0.115, "Mass window around expected (in GeV/c2)"};

  //// Generated particles criteria:
  struct : ConfigurableGroup {
    Configurable<bool> doQA{"doQA", true, "If True, fill QA histos"};
    Configurable<bool> mc_keepOnlyFromGenerator{"mc_keepOnlyFromGenerator", false, "Keep only mcparticles from the generator"};
    Configurable<bool> mc_keepOnlyFromTransport{"mc_keepOnlyFromTransport", false, "Keep only mcparticles from the transport code"};
    Configurable<int> mc_selectMCProcess{"mc_selectMCProcess", -1, "Keep only mcparticles produced in the selected MC process"};
    Configurable<float> mc_rapidityWindow{"mc_rapidityWindow", 0.5, "Max generated particle rapidity"};
  } genSelections;

  // Axis
  // base properties
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {500, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.05f, 1.151f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.5f}, "M_{#Gamma}"};
  ConfigurableAxis axisK0SMass{"axisK0SMass", {200, 0.4f, 0.6f}, "M_{K^{0}}"};

  // topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisXY{"axisXY", {120, -120.0f, 120.0f}, "XY axis"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisPA{"axisPA", {100, 0.0f, 1}, "Pointing angle"};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  ConfigurableAxis axisCandSel{"axisCandSel", {7, 0.5f, +7.5f}, "Candidate Selection"};
  ConfigurableAxis axisIRBinning{"axisIRBinning", {151, -10, 1500}, "Binning for the interaction rate (kHz)"};

  void init(InitContext const&)
  {
    LOGF(info, "Initializing now: cross-checking correctness...");
    if (doprocessRealData +
          doprocessRealDataWithTOF +
          doprocessMonteCarlo +
          doprocessMonteCarloWithTOF >
        1) {
      LOGF(fatal, "You have enabled more than one process function. Please check your configuration! Aborting now.");
    }

    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    histos.add("hEventCentrality", "hEventCentrality", kTH1D, {axisCentrality});

    histos.add("PhotonSel/h2dMassGammaVsK0S", "h2dMassGammaVsK0S", kTH2D, {axisPhotonMass, axisK0SMass});
    histos.add("PhotonSel/h2dMassGammaVsLambda", "h2dMassGammaVsLambda", kTH2D, {axisPhotonMass, axisLambdaMass});
    histos.add("PhotonSel/h2dV0XY", "h2dV0XY", kTH2F, {axisXY, axisXY});

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

    histos.add("LambdaSel/h2dMassLambdaVsK0S", "h2dMassLambdaVsK0S", kTH2D, {axisLambdaMass, axisK0SMass});
    histos.add("LambdaSel/h2dMassLambdaVsGamma", "h2dMassLambdaVsGamma", kTH2D, {axisLambdaMass, axisPhotonMass});
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
    histos.add("SigmaSel/hSigmaMassSelected", "hSigmaMassSelected", kTH1F, {axisSigmaMass});

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
      histos.add("MC/h2dGammaXYConversion", "h2dGammaXYConversion", kTH2F, {axisXY, axisXY});
      histos.add("MC/h2dPtVsCentrality_MCAssocGamma", "h2dPtVsCentrality_MCAssocGamma", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dPtVsCentrality_MCAssocLambda", "h2dPtVsCentrality_MCAssocLambda", kTH2D, {axisCentrality, axisPt});
      histos.add("MC/h2dPtVsCentrality_MCAssocALambda", "h2dPtVsCentrality_MCAssocALambda", kTH2D, {axisCentrality, axisPt});
    }

    if (doprocessGeneratedRun3 && genSelections.doQA) {

      // Pi0s
      histos.add("GenQA/hGenPi0", "hGenPi0", kTH1D, {axisPt});

      auto hPrimaryPi0s = histos.add<TH1>("GenQA/hPrimaryPi0s", "hPrimaryPi0s", kTH1D, {{2, -0.5f, 1.5f}});
      hPrimaryPi0s->GetXaxis()->SetBinLabel(1, "All Pi0s");
      hPrimaryPi0s->GetXaxis()->SetBinLabel(2, "Primary Pi0s");

      histos.add("GenQA/h2dPi0MCSourceVsPDGMother", "h2dPi0MCSourceVsPDGMother", kTHnSparseD, {{2, -0.5f, 1.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("GenQA/h2dPi0NDaughtersVsPDG", "h2dPi0NDaughtersVsPDG", kTHnSparseD, {{10, -0.5f, +9.5f}, {10001, -5000.5f, +5000.5f}});

      auto h2DGenPi0TypeVsProducedByGen = histos.add<TH2>("GenQA/h2DGenPi0TypeVsProducedByGen", "h2DGenPi0TypeVsProducedByGen", kTH2D, {{2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
      h2DGenPi0TypeVsProducedByGen->GetXaxis()->SetBinLabel(1, "Sterile");
      h2DGenPi0TypeVsProducedByGen->GetXaxis()->SetBinLabel(2, "Non-Sterile");
      h2DGenPi0TypeVsProducedByGen->GetYaxis()->SetBinLabel(1, "Generator");
      h2DGenPi0TypeVsProducedByGen->GetYaxis()->SetBinLabel(2, "Transport");

      // ______________________________________________________
      // Sigma0s
      histos.add("GenQA/hGenSigma0", "hGenSigma0", kTH1D, {axisPt});
      histos.add("GenQA/hGenAntiSigma0", "hGenAntiSigma0", kTH1D, {axisPt});

      histos.add("GenQA/h2dGenSigma0xy_Generator", "hGenSigma0xy_Generator", kTH2D, {axisXY, axisXY});
      histos.add("GenQA/h2dGenSigma0xy_Transport", "hGenSigma0xy_Transport", kTH2D, {axisXY, axisXY});
      histos.add("GenQA/hGenSigma0Radius_Generator", "hGenSigma0Radius_Generator", kTH1D, {axisRadius});
      histos.add("GenQA/hGenSigma0Radius_Transport", "hGenSigma0Radius_Transport", kTH1D, {axisRadius});

      histos.add("GenQA/h2dSigma0MCSourceVsPDGMother", "h2dSigma0MCSourceVsPDGMother", kTHnSparseD, {{2, -0.5f, 1.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("GenQA/h2dSigma0NDaughtersVsPDG", "h2dSigma0NDaughtersVsPDG", kTHnSparseD, {{10, -0.5f, +9.5f}, {10001, -5000.5f, +5000.5f}});

      auto hPrimarySigma0s = histos.add<TH1>("GenQA/hPrimarySigma0s", "hPrimarySigma0s", kTH1D, {{2, -0.5f, 1.5f}});
      hPrimarySigma0s->GetXaxis()->SetBinLabel(1, "All Sigma0s");
      hPrimarySigma0s->GetXaxis()->SetBinLabel(2, "Primary Sigma0s");

      auto hGenSpecies = histos.add<TH1>("GenQA/hGenSpecies", "hGenSpecies", kTH1D, {{4, -0.5f, 3.5f}});
      hGenSpecies->GetXaxis()->SetBinLabel(1, "All Prim. Lambda");
      hGenSpecies->GetXaxis()->SetBinLabel(2, "All Prim. ALambda");
      hGenSpecies->GetXaxis()->SetBinLabel(5, "All Sigma0s");
      hGenSpecies->GetXaxis()->SetBinLabel(6, "All ASigma0s");

      histos.add("GenQA/hSigma0NDau", "hSigma0NDau", kTH1D, {{10, -0.5f, +9.5f}});
      histos.add("GenQA/h2dSigma0NDauVsProcess", "h2dSigma0NDauVsProcess", kTH2D, {{10, -0.5f, +9.5f}, {50, -0.5f, 49.5f}});

      auto h2DGenSigma0TypeVsProducedByGen = histos.add<TH2>("GenQA/h2DGenSigma0TypeVsProducedByGen", "h2DGenSigma0TypeVsProducedByGen", kTH2D, {{2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
      h2DGenSigma0TypeVsProducedByGen->GetXaxis()->SetBinLabel(1, "Sterile");
      h2DGenSigma0TypeVsProducedByGen->GetXaxis()->SetBinLabel(2, "Non-Sterile");
      h2DGenSigma0TypeVsProducedByGen->GetYaxis()->SetBinLabel(1, "Generator");
      h2DGenSigma0TypeVsProducedByGen->GetYaxis()->SetBinLabel(2, "Transport");
    }
  }

  // ______________________________________________________
  // Struct to store V0Pair properties
  struct V0PairTopoInfo {
    float X = -999.f;
    float Y = -999.f;
    float Z = -999.f;
    float DCADau = -999.f;
    float CosPA = -1.f;
  };

  // ______________________________________________________
  // Struct to store V0Pair MC properties
  struct V0PairMCInfo {
    bool fIsV01CorrectlyAssign = false;
    bool fIsV02CorrectlyAssign = false;
    bool fIsV01Primary = false;
    bool fIsV02Primary = false;
    bool fV0PairProducedByGenerator = false;
    int V01PDGCode = 0;
    int V02PDGCode = 0;
    int V01PDGCodeMother = 0;
    int V02PDGCodeMother = 0;
    int V0PairPDGCode = 0;
    int V0PairPDGCodeMother = 0;
    int V0PairMCProcess = -1;
    int V0PairMCParticleID = -1;
    float V01MCpx = -999.f;
    float V01MCpy = -999.f;
    float V01MCpz = -999.f;
    float V02MCpx = -999.f;
    float V02MCpy = -999.f;
    float V02MCpz = -999.f;
    float V0PairMCRadius = -999.f;
  };

  // ______________________________________________________
  // Struct to store V0Pair Generated properties
  struct V0PairGenInfo {
    bool IsPrimary = false;
    bool IsV0Lambda = false;
    bool IsV0AntiLambda = false;
    bool IsPi0 = false;
    bool IsSigma0 = false;
    bool IsAntiSigma0 = false;
    bool IsProducedByGenerator = false;
    bool IsSterile = false;
    int MCProcess = -1;
    int MCCollId = -1;
    int PDGCodeMother = 0;
    int NDaughters = -1;
    float MCPt = -999.f;
    float MCvx = 999.f;
    float MCvy = 999.f;
  };

  template <typename TV01, typename TV02>
  V0PairTopoInfo propagateV0PairToDCA(TV01 const& v01, TV02 const& v02)
  {
    V0PairTopoInfo info;

    // Positions
    ROOT::Math::XYZVector v01position(v01.x(), v01.y(), v01.z());
    ROOT::Math::XYZVector v02position(v02.x(), v02.y(), v02.z());

    // Momenta
    ROOT::Math::XYZVector v01momentum(v01.px(), v01.py(), v01.pz());
    ROOT::Math::XYZVector v02momentum(v02.px(), v02.py(), v02.pz());

    // Momenta (normalized)
    ROOT::Math::XYZVector v01momentumNorm(v01.px() / v01.p(), v01.py() / v01.p(), v01.pz() / v01.p());
    ROOT::Math::XYZVector v02momentumNorm(v02.px() / v02.p(), v02.py() / v02.p(), v02.pz() / v02.p());

    // DCADau calculation (using full momenta for precision)
    ROOT::Math::XYZVector posdiff = v02position - v01position;
    ROOT::Math::XYZVector cross = v01momentum.Cross(v02momentum);

    float d = 1.0f - TMath::Power(v01momentumNorm.Dot(v02momentumNorm), 2);
    float t = posdiff.Dot(v01momentumNorm - v01momentumNorm.Dot(v02momentumNorm) * v02momentumNorm) / d;
    float s = -posdiff.Dot(v02momentumNorm - v01momentumNorm.Dot(v02momentumNorm) * v01momentumNorm) / d;

    ROOT::Math::XYZVector pointOn1 = v01position + t * v01momentumNorm;
    ROOT::Math::XYZVector pointOn2 = v02position + s * v02momentumNorm;
    ROOT::Math::XYZVector PCA = 0.5 * (pointOn1 + pointOn2);

    // Calculate properties and fill struct
    info.DCADau = (cross.Mag2() > 0) ? std::abs(posdiff.Dot(cross)) / cross.R() : 999.f;
    info.CosPA = v01momentumNorm.Dot(v02momentumNorm);

    if (d < 1e-5f) {                  // Parallel or nearly parallel lines
      info.X = info.Y = info.Z = 0.f; // should we use another dummy value? Perhaps 999.f?
      return info;
    }

    info.X = PCA.X();
    info.Y = PCA.Y();
    info.Z = PCA.Z();

    return info;
  }

  template <typename TV01, typename TV02, typename TCollision, typename TMCParticles>
  V0PairMCInfo getV0PairMCInfo(TV01 const& v01, TV02 const& v02, TCollision const& collision, TMCParticles const& mcparticles)
  {
    V0PairMCInfo MCinfo;

    if (!v01.has_v0MCCore() || !v02.has_v0MCCore())
      return MCinfo;

    auto v01MC = v01.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
    auto v02MC = v02.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

    if (collision.has_straMCCollision()) {
      auto MCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
      MCinfo.fIsV01CorrectlyAssign = (v01MC.straMCCollisionId() == MCCollision.globalIndex());
      MCinfo.fIsV02CorrectlyAssign = (v02MC.straMCCollisionId() == MCCollision.globalIndex());
    }

    MCinfo.V01MCpx = v01MC.pxMC();
    MCinfo.V01MCpy = v01MC.pyMC();
    MCinfo.V01MCpz = v01MC.pzMC();
    MCinfo.V02MCpx = v02MC.pxMC();
    MCinfo.V02MCpy = v02MC.pyMC();
    MCinfo.V02MCpz = v02MC.pzMC();

    // Get corresponding entries in MCParticles table
    auto MCParticle_v01 = mcparticles.rawIteratorAt(v01MC.particleIdMC());
    auto MCParticle_v02 = mcparticles.rawIteratorAt(v02MC.particleIdMC());

    // Get MC Mothers
    auto const& MCMothersList_v01 = MCParticle_v01.template mothers_as<aod::McParticles>();
    auto const& MCMothersList_v02 = MCParticle_v02.template mothers_as<aod::McParticles>();

    if (!MCMothersList_v01.empty() && !MCMothersList_v02.empty()) { // Are there mothers?

      auto const& MCMother_v01 = MCMothersList_v01.front(); // First mother
      auto const& MCMother_v02 = MCMothersList_v02.front(); // First mother

      if (MCMother_v01.globalIndex() == MCMother_v02.globalIndex()) { // Is it the same mother?

        MCinfo.fV0PairProducedByGenerator = MCMother_v01.producedByGenerator();
        MCinfo.V0PairPDGCode = MCMother_v01.pdgCode();
        MCinfo.V0PairMCProcess = MCMother_v01.getProcess();
        MCinfo.V0PairMCParticleID = MCMother_v01.globalIndex();
        MCinfo.V0PairMCRadius = std::hypot(MCMother_v01.vx(), MCMother_v01.vy()); // production position radius

        auto const& v0pairmothers = MCMother_v01.template mothers_as<aod::McParticles>(); // Get mothers
        if (!v0pairmothers.empty()) {
          auto& v0PairMother = v0pairmothers.front(); // V0Pair mother, V0s grandmother
          MCinfo.V0PairPDGCodeMother = v0PairMother.pdgCode();
        }
      }
    }

    MCinfo.fIsV01Primary = v01MC.isPhysicalPrimary();
    MCinfo.fIsV02Primary = v02MC.isPhysicalPrimary();
    MCinfo.V01PDGCode = v01MC.pdgCode();
    MCinfo.V02PDGCode = v02MC.pdgCode();
    MCinfo.V01PDGCodeMother = v01MC.pdgCodeMother();
    MCinfo.V02PDGCodeMother = v02MC.pdgCodeMother();

    return MCinfo;
  }

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

  template <typename TMCParticle>
  V0PairGenInfo getV0PairGenInfo(TMCParticle const& mcParticle)
  {
    V0PairGenInfo GenInfo; // auxiliary struct to store info

    // Fill with properties
    GenInfo.IsPrimary = mcParticle.isPhysicalPrimary();
    GenInfo.IsV0Lambda = mcParticle.pdgCode() == 3122;
    GenInfo.IsV0AntiLambda = mcParticle.pdgCode() == -3122;
    GenInfo.IsPi0 = mcParticle.pdgCode() == 111;
    GenInfo.IsSigma0 = mcParticle.pdgCode() == 3212;
    GenInfo.IsAntiSigma0 = mcParticle.pdgCode() == -3212;
    GenInfo.IsProducedByGenerator = mcParticle.producedByGenerator();
    GenInfo.MCProcess = mcParticle.getProcess();
    GenInfo.MCPt = mcParticle.pt();
    GenInfo.MCvx = mcParticle.vx(); // production position X
    GenInfo.MCvy = mcParticle.vy(); // production position Y

    if (mcParticle.has_mcCollision())
      GenInfo.MCCollId = mcParticle.mcCollisionId(); // save this reference, please

    // Checking decay mode if sigma0
    if (GenInfo.IsSigma0 || GenInfo.IsAntiSigma0 || GenInfo.IsPi0) {

      // This is a costly operation, so we do it only for pi0s and sigma0s
      auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
      GenInfo.NDaughters = daughters.size();
      GenInfo.IsSterile = daughters.size() == 0;

      auto const& GenMothersList = mcParticle.template mothers_as<aod::McParticles>();
      GenInfo.PDGCodeMother = (!GenMothersList.empty()) ? GenMothersList.front().pdgCode() : 0;

      if ((GenInfo.IsSigma0 || GenInfo.IsAntiSigma0) && genSelections.doQA) {
        histos.fill(HIST("GenQA/h2dSigma0MCSourceVsPDGMother"), GenInfo.IsProducedByGenerator, GenInfo.PDGCodeMother);
        for (auto& daughter : daughters) // checking decay modes
          histos.fill(HIST("GenQA/h2dSigma0NDaughtersVsPDG"), daughters.size(), daughter.pdgCode());
      }

      if (GenInfo.IsPi0 && genSelections.doQA) {
        histos.fill(HIST("GenQA/h2dPi0MCSourceVsPDGMother"), GenInfo.IsProducedByGenerator, GenInfo.PDGCodeMother);
        for (auto& daughter : daughters) // checking decay modes
          histos.fill(HIST("GenQA/h2dPi0NDaughtersVsPDG"), daughters.size(), daughter.pdgCode());
      }
    }
    return GenInfo;
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void fillGenQAHistos(V0PairGenInfo const& GenInfo)
  {
    if (GenInfo.IsPi0) {
      histos.fill(HIST("GenQA/hGenPi0"), GenInfo.MCPt);
      histos.fill(HIST("GenQA/hPrimaryPi0s"), 0);
      if (GenInfo.IsPrimary)
        histos.fill(HIST("GenQA/hPrimaryPi0s"), 1);

      if (GenInfo.IsSterile) {
        if (GenInfo.IsProducedByGenerator)
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 0, 0);
        else
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 0, 1);
      } else {
        if (GenInfo.IsProducedByGenerator)
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 1, 0);
        else
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 1, 1);
      }
    }

    if (GenInfo.IsV0Lambda && GenInfo.IsPrimary)
      histos.fill(HIST("GenQA/hGenSpecies"), 0);
    if (GenInfo.IsV0AntiLambda && GenInfo.IsPrimary)
      histos.fill(HIST("GenQA/hGenSpecies"), 1);

    // Checking decay mode
    if (GenInfo.IsSigma0 || GenInfo.IsAntiSigma0) {
      histos.fill(HIST("GenQA/hSigma0NDau"), GenInfo.NDaughters);
      histos.fill(HIST("GenQA/h2dSigma0NDauVsProcess"), GenInfo.NDaughters, GenInfo.MCProcess);

      const auto radius = std::hypot(GenInfo.MCvx, GenInfo.MCvy);
      // Sigma0 XY and radius (separate histos for Gen/Transport)
      if (GenInfo.IsProducedByGenerator) {
        histos.fill(HIST("GenQA/h2dGenSigma0xy_Generator"), GenInfo.MCvx, GenInfo.MCvy);
        histos.fill(HIST("GenQA/hGenSigma0Radius_Generator"), radius);
      } else {
        histos.fill(HIST("GenQA/h2dGenSigma0xy_Transport"), GenInfo.MCvx, GenInfo.MCvy);
        histos.fill(HIST("GenQA/hGenSigma0Radius_Transport"), radius);
      }

      // Sigma0 type vs origin (single 2D histo)
      const int genIndex = GenInfo.IsProducedByGenerator ? 0 : 1; // 0 = Generator, 1 = Transport
      const int typeIndex = GenInfo.IsSterile ? 0 : 1;            // 0 = Sterile,   1 = Normal
      histos.fill(HIST("GenQA/h2DGenSigma0TypeVsProducedByGen"), typeIndex, genIndex);

      // Fill histograms
      if (GenInfo.IsSigma0) {
        histos.fill(HIST("GenQA/hGenSpecies"), 2);
        histos.fill(HIST("GenQA/hGenSigma0"), GenInfo.MCPt);

        histos.fill(HIST("GenQA/hPrimarySigma0s"), 0);
        if (GenInfo.IsPrimary)
          histos.fill(HIST("GenQA/hPrimarySigma0s"), 1);
      }
      if (GenInfo.IsAntiSigma0) {
        histos.fill(HIST("GenQA/hGenSpecies"), 3);
        histos.fill(HIST("GenQA/hGenAntiSigma0"), GenInfo.MCPt);
      }
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  template <typename TMCParticles>
  void genProcess(TMCParticles const& mcParticles)
  {
    for (auto& mcParticle : mcParticles) {
      // Rapidity selection
      if (TMath::Abs(mcParticle.y()) > genSelections.mc_rapidityWindow)
        continue;

      // Selection on the source (generator/transport)
      if (genSelections.mc_keepOnlyFromGenerator && !genSelections.mc_keepOnlyFromTransport) {
        if (!mcParticle.producedByGenerator())
          continue;
      }

      if (genSelections.mc_keepOnlyFromTransport && !genSelections.mc_keepOnlyFromGenerator) {
        if (mcParticle.producedByGenerator())
          continue;
      }

      // MC Process selection
      if ((genSelections.mc_selectMCProcess >= 0) && (genSelections.mc_selectMCProcess != mcParticle.getProcess()))
        continue;

      // Get generated particle info
      auto MCGenInfo = getV0PairGenInfo(mcParticle);

      // Fill QA histos
      if (genSelections.doQA)
        fillGenQAHistos(MCGenInfo);

      // Fill tables
      // Pi0
      if (fillPi0Tables && MCGenInfo.IsPi0) {
        pi0Gens(MCGenInfo.IsProducedByGenerator, MCGenInfo.MCPt); // optional table to store generated pi0 candidates. Be careful, this is a large table!
        pi0GenCollRefs(MCGenInfo.MCCollId);                       // link to stramccollision table
      }

      // Sigma0/ASigma0
      if (fillSigma0Tables && (MCGenInfo.IsSigma0 || MCGenInfo.IsAntiSigma0)) {
        sigma0Gens(MCGenInfo.IsSigma0, MCGenInfo.IsProducedByGenerator, MCGenInfo.MCPt);
        sigma0GenCollRefs(MCGenInfo.MCCollId); // link to stramccollision table
      }
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
    float PhotonY = RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma);
    histos.fill(HIST("PhotonSel/h2dMassGammaVsK0S"), gamma.mGamma(), gamma.mK0Short());
    histos.fill(HIST("PhotonSel/h2dMassGammaVsLambda"), gamma.mGamma(), gamma.mLambda());

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

    histos.fill(HIST("PhotonSel/h2dV0XY"), gamma.x(), gamma.y());
    histos.fill(HIST("PhotonSel/h3dPhotonMass"), centrality, gamma.pt(), gamma.mGamma());

    //_______________________________________________
    // MC Processing
    if constexpr (requires { gamma.motherMCPartId(); }) {
      if (gamma.has_v0MCCore()) {
        auto gammaMC = gamma.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
        if (gammaMC.pdgCode() == 22) {
          histos.fill(HIST("MC/h2dGammaXYConversion"), gamma.x(), gamma.y());
          histos.fill(HIST("MC/h2dPtVsCentrality_MCAssocGamma"), centrality, gamma.pt());
        }
      }
    }

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
    histos.fill(HIST("LambdaSel/h2dMassLambdaVsK0S"), lambda.mLambda(), lambda.mK0Short());
    histos.fill(HIST("LambdaSel/h2dMassLambdaVsGamma"), lambda.mLambda(), lambda.mGamma());

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

    //_______________________________________________
    // MC Processing (if available)
    if constexpr (requires { lambda.motherMCPartId(); }) {
      if (lambda.has_v0MCCore()) {
        auto lambdaMC = lambda.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
        if (lambdaMC.pdgCode() == 3122) // Is Lambda
          histos.fill(HIST("MC/h2dPtVsCentrality_MCAssocLambda"), centrality, lambda.pt());
        if (lambdaMC.pdgCode() == -3122) // Is AntiLambda
          histos.fill(HIST("MC/h2dPtVsCentrality_MCAssocALambda"), centrality, lambda.pt());
      }
    }

    return true;
  }

  //_______________________________________________
  // Build pi0 candidate for QA
  template <typename TV0Object, typename TCollision, typename TMCParticles>
  bool buildPi0(TV0Object const& gamma1, TV0Object const& gamma2, TCollision const& collision, TMCParticles const& mcparticles)
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

    // Calculate Pi0 topological info
    auto pi0TopoInfo = propagateV0PairToDCA(gamma1, gamma2);

    // Check if MC data and populate corresponding table
    if constexpr (requires { gamma1.motherMCPartId(); gamma2.motherMCPartId(); }) {
      auto pi0MCInfo = getV0PairMCInfo(gamma1, gamma2, collision, mcparticles);

      pi0coresmc(pi0MCInfo.V0PairMCRadius, pi0MCInfo.V0PairPDGCode, pi0MCInfo.V0PairPDGCodeMother, pi0MCInfo.V0PairMCProcess, pi0MCInfo.fV0PairProducedByGenerator,
                 pi0MCInfo.V01MCpx, pi0MCInfo.V01MCpy, pi0MCInfo.V01MCpz,
                 pi0MCInfo.fIsV01Primary, pi0MCInfo.V01PDGCode, pi0MCInfo.V01PDGCodeMother, pi0MCInfo.fIsV01CorrectlyAssign,
                 pi0MCInfo.V02MCpx, pi0MCInfo.V02MCpy, pi0MCInfo.V02MCpz,
                 pi0MCInfo.fIsV02Primary, pi0MCInfo.V02PDGCode, pi0MCInfo.V02PDGCodeMother, pi0MCInfo.fIsV02CorrectlyAssign);
    }

    pi0cores(pi0TopoInfo.X, pi0TopoInfo.Y, pi0TopoInfo.Z, pi0TopoInfo.DCADau, pi0TopoInfo.CosPA,
             gamma1.px(), gamma1.py(), gamma1.pz(),
             gamma1.mGamma(), gamma1.qtarm(), gamma1.alpha(), gamma1.dcapostopv(), gamma1.dcanegtopv(), gamma1.dcaV0daughters(),
             gamma1.negativeeta(), gamma1.positiveeta(), gamma1.v0cosPA(), gamma1.v0radius(), gamma1.z(),
             posTrackGamma1.tpcCrossedRows(), negTrackGamma1.tpcCrossedRows(), posTrackGamma1.tpcNSigmaEl(), negTrackGamma1.tpcNSigmaEl(), gamma1.v0Type(),
             gamma2.px(), gamma2.py(), gamma2.pz(),
             gamma2.mGamma(), gamma2.qtarm(), gamma2.alpha(), gamma2.dcapostopv(), gamma2.dcanegtopv(), gamma2.dcaV0daughters(),
             gamma2.negativeeta(), gamma2.positiveeta(), gamma2.v0cosPA(), gamma2.v0radius(), gamma2.z(),
             posTrackGamma2.tpcCrossedRows(), negTrackGamma2.tpcCrossedRows(), posTrackGamma2.tpcNSigmaEl(), negTrackGamma2.tpcNSigmaEl(), gamma2.v0Type());

    pi0coresRefs(collision.globalIndex());

    return true;
  }

  //_______________________________________________
  // Build sigma0 candidate
  template <typename TV0Object, typename TCollision, typename TMCParticles>
  bool buildSigma0(TV0Object const& lambda, TV0Object const& gamma, TCollision const& collision, TMCParticles const& mcparticles)
  {
    //_______________________________________________
    // Checking if both V0s are made of the very same tracks
    if (gamma.posTrackExtraId() == lambda.posTrackExtraId() ||
        gamma.negTrackExtraId() == lambda.negTrackExtraId()) {
      return false;
    }

    //_______________________________________________
    // Sigma0 pre-selections
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};

    auto arrMom = std::array{pVecPhotons, pVecLambda};
    float sigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    float sigmaY = -999.f;

    if constexpr (requires { gamma.pxMC(); lambda.pxMC(); }) // If MC
      sigmaY = RecoDecay::y(std::array{gamma.pxMC() + lambda.pxMC(), gamma.pyMC() + lambda.pyMC(), gamma.pzMC() + lambda.pzMC()}, o2::constants::physics::MassSigma0);
    else // If DATA
      sigmaY = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);

    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 1.);
    if (TMath::Abs(sigmaMass - o2::constants::physics::MassSigma0) > Sigma0Window)
      return false;

    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 2.);
    if (TMath::Abs(sigmaY) > SigmaMaxRap)
      return false;

    histos.fill(HIST("SigmaSel/hSigmaMassSelected"), sigmaMass);
    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 3.);
    //_______________________________________________
    // Calculate properties & Fill tables

    // Sigma0 topological info
    auto sigma0TopoInfo = propagateV0PairToDCA(gamma, lambda);

    sigma0cores(sigma0TopoInfo.X, sigma0TopoInfo.Y, sigma0TopoInfo.Z, sigma0TopoInfo.DCADau,
                gamma.px(), gamma.py(), gamma.pz(), gamma.mGamma(), lambda.px(), lambda.py(), lambda.pz(), lambda.mLambda(), lambda.mAntiLambda());

    // MC properties
    if constexpr (requires { gamma.motherMCPartId(); lambda.motherMCPartId(); }) {
      auto sigma0MCInfo = getV0PairMCInfo(gamma, lambda, collision, mcparticles);

      sigma0mccores(sigma0MCInfo.V0PairMCRadius, sigma0MCInfo.V0PairPDGCode, sigma0MCInfo.V0PairPDGCodeMother, sigma0MCInfo.V0PairMCProcess, sigma0MCInfo.fV0PairProducedByGenerator,
                    sigma0MCInfo.V01MCpx, sigma0MCInfo.V01MCpy, sigma0MCInfo.V01MCpz,
                    sigma0MCInfo.fIsV01Primary, sigma0MCInfo.V01PDGCode, sigma0MCInfo.V01PDGCodeMother, sigma0MCInfo.fIsV01CorrectlyAssign,
                    sigma0MCInfo.V02MCpx, sigma0MCInfo.V02MCpy, sigma0MCInfo.V02MCpz,
                    sigma0MCInfo.fIsV02Primary, sigma0MCInfo.V02PDGCode, sigma0MCInfo.V02PDGCodeMother, sigma0MCInfo.fIsV02CorrectlyAssign);
    }

    // Sigma0s -> stracollisions link
    sigma0CollRefs(collision.globalIndex());

    //_______________________________________________
    // Photon extra properties
    auto posTrackGamma = gamma.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma = gamma.template negTrackExtra_as<dauTracks>();

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

    sigmaPhotonExtras(gamma.qtarm(), gamma.alpha(), gamma.v0cosPA(), gamma.dcaV0daughters(), gamma.dcanegtopv(), gamma.dcapostopv(), gamma.v0radius(), gamma.z(),
                      posTrackGamma.tpcNSigmaEl(), negTrackGamma.tpcNSigmaEl(), posTrackGamma.tpcCrossedRows(), negTrackGamma.tpcCrossedRows(),
                      gamma.positiveeta(), gamma.negativeeta(), gamma.psipair(), posTrackGamma.itsNCls(), negTrackGamma.itsNCls(), posTrackGamma.itsChi2PerNcl(), negTrackGamma.itsChi2PerNcl(),
                      fPhotonPosTrackCode, fPhotonNegTrackCode, gamma.v0Type());

    //_______________________________________________
    // Lambda extra properties
    auto posTrackLambda = lambda.template posTrackExtra_as<dauTracks>();
    auto negTrackLambda = lambda.template negTrackExtra_as<dauTracks>();

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

    float fLambdaLifeTime = lambda.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    float fLambdaPrTOFNSigma = -999.f;
    float fLambdaPiTOFNSigma = -999.f;
    float fALambdaPrTOFNSigma = -999.f;
    float fALambdaPiTOFNSigma = -999.f;

    if constexpr (requires { lambda.tofNSigmaLaPr(); }) { // If TOF info avaiable
      fLambdaPrTOFNSigma = lambda.tofNSigmaLaPr();
      fLambdaPiTOFNSigma = lambda.tofNSigmaLaPi();
      fALambdaPrTOFNSigma = lambda.tofNSigmaALaPr();
      fALambdaPiTOFNSigma = lambda.tofNSigmaALaPi();
    }

    sigmaLambdaExtras(lambda.qtarm(), lambda.alpha(), fLambdaLifeTime, lambda.v0radius(), lambda.v0cosPA(), lambda.dcaV0daughters(), lambda.dcanegtopv(), lambda.dcapostopv(),
                      posTrackLambda.tpcNSigmaPr(), posTrackLambda.tpcNSigmaPi(), negTrackLambda.tpcNSigmaPr(), negTrackLambda.tpcNSigmaPi(),
                      fLambdaPrTOFNSigma, fLambdaPiTOFNSigma, fALambdaPrTOFNSigma, fALambdaPiTOFNSigma,
                      posTrackLambda.tpcCrossedRows(), negTrackLambda.tpcCrossedRows(),
                      lambda.positiveeta(), lambda.negativeeta(),
                      posTrackLambda.itsNCls(), negTrackLambda.itsNCls(), posTrackLambda.itsChi2PerNcl(), negTrackLambda.itsChi2PerNcl(),
                      fLambdaPosTrackCode, fLambdaNegTrackCode, lambda.v0Type());

    return true;
  }

  // Process photon and lambda candidates to build sigma0 candidates
  template <typename TCollision, typename TV0s, typename TMCParticles>
  void dataProcess(TCollision const& collisions, TV0s const& fullV0s, TMCParticles const& mcparticles)
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
        if (fillSigma0Tables) {
          for (size_t j = 0; j < bestLambdasArray.size(); ++j) {
            auto lambda = fullV0s.rawIteratorAt(bestLambdasArray[j]);

            // Building sigma0 candidate & filling tables
            if (!buildSigma0(lambda, gamma1, coll, mcparticles))
              continue;
          }
        }

        //_______________________________________________
        // pi0 loop
        if (fillPi0Tables) {
          for (size_t j = i + 1; j < bestGammasArray.size(); ++j) {
            auto gamma2 = fullV0s.rawIteratorAt(bestGammasArray[j]);

            // Building pi0 candidate & filling tables
            if (!buildPi0(gamma1, gamma2, coll, mcparticles))
              continue;
          }
        }
      }
    }
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0StandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    dataProcess(collisions, fullV0s, nullptr);
  }

  void processRealDataWithTOF(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0TOFStandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    dataProcess(collisions, fullV0s, nullptr);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0DerivedMCDatas const& fullV0s, aod::McParticles const& mcParticles, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    dataProcess(collisions, fullV0s, mcParticles);
  }

  void processMonteCarloWithTOF(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0TOFDerivedMCDatas const& fullV0s, aod::McParticles const& mcParticles, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    dataProcess(collisions, fullV0s, mcParticles);
  }

  // Simulated processing in Run 3 - run this over original AO2Ds
  void processGeneratedRun3(aod::McParticles const& mcParticles)
  {
    genProcess(mcParticles);
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
