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
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCDatas>;
using V0MLDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0GammaMLScores, aod::V0AntiLambdaMLScores>;
using V0StandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;

struct sigma0builder {
  SliceCache cache;

  Produces<aod::Sigma0Collision> sigma0Coll;          // characterises collisions
  Produces<aod::Sigma0CollRefs> sigma0CollRefs;       // characterises collisions
  Produces<aod::Sigma0Cores> sigma0cores;             // save sigma0 candidates for analysis
  Produces<aod::SigmaPhotonExtras> sigmaPhotonExtras; // save sigma0 candidates for analysis
  Produces<aod::SigmaLambdaExtras> sigmaLambdaExtras; // save sigma0 candidates for analysis
  Produces<aod::SigmaMCCores> sigma0mccores;

  // For manual sliceBy
  Preslice<V0DerivedMCDatas> perCollisionMCDerived = o2::aod::v0data::straCollisionId;
  Preslice<V0StandardDerivedDatas> perCollisionSTDDerived = o2::aod::v0data::straCollisionId;
  Preslice<V0MLDerivedDatas> perCollisionMLDerived = o2::aod::v0data::straCollisionId;

  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // For ML Selection
  Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
  Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
  Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};

  // For standard approach:
  //// Lambda criteria:
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

  // Axis
  // base properties
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
  ConfigurableAxis axisDeltaPt{"axisDeltaPt", {100, -1.0, +1.0}, "#Delta(p_{T})"};

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
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  ConfigurableAxis axisCandSel{"axisCandSel", {13, 0.5f, +13.5f}, "Candidate Selection"};

  int nSigmaCandidates = 0;
  void init(InitContext const&)
  {
    // Event counter
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {axisCentrality});
    histos.add("hCandidateBuilderSelection", "hCandidateBuilderSelection", kTH1F, {axisCandSel});
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(2, "Photon Mass Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(3, "Photon DauEta Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(4, "Photon DCAToPV Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(5, "Photon DCADau Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(6, "Photon Radius Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(7, "Lambda Mass Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(8, "Lambda DauEta Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(9, "Lambda DCAToPV Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(10, "Lambda Radius Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(11, "Lambda DCADau Cut");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(12, "Sigma Mass Window");
    histos.get<TH1>(HIST("hCandidateBuilderSelection"))->GetXaxis()->SetBinLabel(13, "Sigma Y Window");

    // For QA:
    histos.add("Selection/hPhotonMass", "hPhotonMass", kTH1F, {axisPhotonMass});
    histos.add("Selection/hPhotonNegEta", "hPhotonNegEta", kTH1F, {axisRapidity});
    histos.add("Selection/hPhotonPosEta", "hPhotonPosEta", kTH1F, {axisRapidity});
    histos.add("Selection/hPhotonDCANegToPV", "hPhotonDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("Selection/hPhotonDCAPosToPV", "hPhotonDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("Selection/hPhotonDCADau", "hPhotonDCADau", kTH1F, {axisDCAdau});
    histos.add("Selection/hPhotonRadius", "hPhotonRadius", kTH1F, {axisRadius});
    histos.add("Selection/hLambdaMass", "hLambdaMass", kTH1F, {axisLambdaMass});
    histos.add("Selection/hAntiLambdaMass", "hAntiLambdaMass", kTH1F, {axisLambdaMass});
    histos.add("Selection/hLambdaNegEta", "hLambdaNegEta", kTH1F, {axisRapidity});
    histos.add("Selection/hLambdaPosEta", "hLambdaPosEta", kTH1F, {axisRapidity});
    histos.add("Selection/hLambdaDCANegToPV", "hLambdaDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("Selection/hLambdaDCAPosToPV", "hLambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("Selection/hLambdaDCADau", "hLambdaDCADau", kTH1F, {axisDCAdau});
    histos.add("Selection/hLambdaRadius", "hLambdaRadius", kTH1F, {axisRadius});
    histos.add("Selection/hSigmaMass", "hSigmaMass", kTH1F, {axisSigmaMass});
    histos.add("Selection/hSigmaMassWindow", "hSigmaMassWindow", kTH1F, {{1000, -0.09f, 0.11f}});
    histos.add("Selection/hSigmaY", "hSigmaY", kTH1F, {axisRapidity});

    histos.add("GeneralQA/h2dMassGammaVsK0S", "h2dMassGammaVsK0S", kTH2D, {axisPhotonMass, {200, 0.4f, 0.6f}});
    histos.add("GeneralQA/h2dMassLambdaVsK0S", "h2dMassLambdaVsK0S", kTH2D, {axisLambdaMass, {200, 0.4f, 0.6f}});
    histos.add("GeneralQA/h2dMassGammaVsLambda", "h2dMassGammaVsLambda", kTH2D, {axisPhotonMass, axisLambdaMass});
    histos.add("GeneralQA/h3dMassSigma0VsDaupTs", "h3dMassSigma0VsDaupTs", kTH3F, {axisPt, axisPt, axisSigmaMass});
    histos.add("GeneralQA/h2dMassGammaVsK0SAfterMassSel", "h2dMassGammaVsK0SAfterMassSel", kTH2D, {axisPhotonMass, {200, 0.4f, 0.6f}});
    histos.add("GeneralQA/h2dMassLambdaVsK0SAfterMassSel", "h2dMassLambdaVsK0SAfterMassSel", kTH2D, {axisLambdaMass, {200, 0.4f, 0.6f}});
    histos.add("GeneralQA/h2dMassGammaVsLambdaAfterMassSel", "h2dMassGammaVsLambdaAfterMassSel", kTH2D, {axisPhotonMass, axisLambdaMass});

    // MC
    histos.add("MC/h2dPtVsCentrality_GammaBeforeSel", "h2dPtVsCentrality_GammaBeforeSel", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_LambdaBeforeSel", "h2dPtVsCentrality_LambdaBeforeSel", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_AntiLambdaBeforeSel", "h2dPtVsCentrality_AntiLambdaBeforeSel", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_GammaSigma0", "h2dPtVsCentrality_GammaSigma0", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_LambdaSigma0", "h2dPtVsCentrality_LambdaSigma0", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_Sigma0BeforeSel", "h2dPtVsCentrality_Sigma0BeforeSel", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_Sigma0AfterSel", "h2dPtVsCentrality_Sigma0AfterSel", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_AntiSigma0BeforeSel", "h2dPtVsCentrality_AntiSigma0BeforeSel", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_GammaAntiSigma0", "h2dPtVsCentrality_GammaAntiSigma0", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_LambdaAntiSigma0", "h2dPtVsCentrality_LambdaAntiSigma0", kTH2D, {axisCentrality, axisPt});
    histos.add("MC/h2dPtVsCentrality_AntiSigma0AfterSel", "h2dPtVsCentrality_AntiSigma0AfterSel", kTH2D, {axisCentrality, axisPt});

    // Sigma vs Daughters pT
    histos.add("MC/h2dSigmaPtVsLambdaPt", "h2dSigmaPtVsLambdaPt", kTH2D, {axisPt, axisPt});
    histos.add("MC/h2dSigmaPtVsGammaPt", "h2dSigmaPtVsGammaPt", kTH2D, {axisPt, axisPt});

    // pT Resolution:
    histos.add("MC/h2dLambdaPtResolution", "h2dLambdaPtResolution", kTH2D, {axisPt, axisDeltaPt});
    histos.add("MC/h2dGammaPtResolution", "h2dGammaPtResolution", kTH2D, {axisPt, axisDeltaPt});

    // For background decomposition
    histos.add("MC/h2dPtVsMassSigma_All", "h2dPtVsMassSigma_All", kTH2D, {axisPt, axisSigmaMass});
    histos.add("MC/h2dPtVsMassSigma_SignalOnly", "h2dPtVsMassSigma_SignalOnly", kTH2D, {axisPt, axisSigmaMass});
    histos.add("MC/h2dPtVsMassSigma_TrueDaughters", "h2dPtVsMassSigma_TrueDaughters", kTH2D, {axisPt, axisSigmaMass});
    histos.add("MC/h2dPtVsMassSigma_TrueGammaFakeLambda", "h2dPtVsMassSigma_TrueGammaFakeLambda", kTH2D, {axisPt, axisSigmaMass});
    histos.add("MC/h2dPtVsMassSigma_FakeGammaTrueLambda", "h2dPtVsMassSigma_FakeGammaTrueLambda", kTH2D, {axisPt, axisSigmaMass});
    histos.add("MC/h2dPtVsMassSigma_FakeDaughters", "h2dPtVsMassSigma_FakeDaughters", kTH2D, {axisPt, axisSigmaMass});
    histos.add("MC/h2dTrueDaughtersMatrix", "h2dTrueDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
    histos.add("MC/h2dTrueGammaFakeLambdaMatrix", "h2dTrueGammaFakeLambdaMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
    histos.add("MC/h2dFakeGammaTrueLambdaMatrix", "h2dFakeGammaTrueLambdaMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
    histos.add("MC/h2dFakeDaughtersMatrix", "h2dFakeDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});

    histos.add("h3dMassSigmasBeforeSel", "h3dMassSigmasBeforeSel", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("h3dMassSigmasAfterSel", "h3dMassSigmasAfterSel", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
  }

  // Process sigma candidate and store properties in object
  template <typename TV0Object>
  bool processSigmaCandidate(TV0Object const& lambda, TV0Object const& gamma)
  {
    if ((lambda.v0Type() == 0) || (gamma.v0Type() == 0))
      return false;

    // Checking if both V0s are made of the very same tracks
    if ((gamma.posTrackExtraId() == lambda.posTrackExtraId()) || (gamma.negTrackExtraId() == lambda.negTrackExtraId()) || (gamma.posTrackExtraId() == lambda.negTrackExtraId()) || (gamma.negTrackExtraId() == lambda.posTrackExtraId()) || (gamma.posTrackExtraId() == lambda.negTrackExtraId()))
      return false;

    if constexpr (
      requires { gamma.gammaBDTScore(); } &&
      requires { lambda.lambdaBDTScore(); } &&
      requires { lambda.antiLambdaBDTScore(); }) {

      LOGF(info, "X-check: ML Selection is on!");
      // Gamma selection:
      if (gamma.gammaBDTScore() <= Gamma_MLThreshold)
        return false;

      // Lambda and AntiLambda selection
      if ((lambda.lambdaBDTScore() <= Lambda_MLThreshold) && (lambda.antiLambdaBDTScore() <= AntiLambda_MLThreshold))
        return false;

    } else {
      // Standard selection
      // Gamma basic selection criteria:
      histos.fill(HIST("hCandidateBuilderSelection"), 1.);
      histos.fill(HIST("Selection/hPhotonMass"), gamma.mGamma());
      if ((gamma.mGamma() < 0) || (gamma.mGamma() > PhotonMaxMass))
        return false;
      histos.fill(HIST("Selection/hPhotonNegEta"), gamma.negativeeta());
      histos.fill(HIST("Selection/hPhotonPosEta"), gamma.positiveeta());
      histos.fill(HIST("hCandidateBuilderSelection"), 2.);
      if ((TMath::Abs(gamma.negativeeta()) > PhotonMaxDauPseudoRap) || (TMath::Abs(gamma.positiveeta()) > PhotonMaxDauPseudoRap))
        return false;
      histos.fill(HIST("Selection/hPhotonDCANegToPV"), TMath::Abs(gamma.dcanegtopv()));
      histos.fill(HIST("Selection/hPhotonDCAPosToPV"), TMath::Abs(gamma.dcapostopv()));
      histos.fill(HIST("hCandidateBuilderSelection"), 3.);
      if ((TMath::Abs(gamma.dcapostopv()) < PhotonMinDCAToPv) || (TMath::Abs(gamma.dcanegtopv()) < PhotonMinDCAToPv))
        return false;
      histos.fill(HIST("Selection/hPhotonDCADau"), TMath::Abs(gamma.dcaV0daughters()));
      histos.fill(HIST("hCandidateBuilderSelection"), 4.);
      if (TMath::Abs(gamma.dcaV0daughters()) > PhotonMaxDCAV0Dau)
        return false;
      histos.fill(HIST("Selection/hPhotonRadius"), gamma.v0radius());
      histos.fill(HIST("hCandidateBuilderSelection"), 5.);
      if ((gamma.v0radius() < PhotonMinRadius) || (gamma.v0radius() > PhotonMaxRadius))
        return false;

      histos.fill(HIST("hCandidateBuilderSelection"), 6.);
      histos.fill(HIST("Selection/hLambdaMass"), lambda.mLambda());
      histos.fill(HIST("Selection/hAntiLambdaMass"), lambda.mAntiLambda());
      // Lambda basic selection criteria:
      if ((TMath::Abs(lambda.mLambda() - 1.115683) > LambdaWindow) && (TMath::Abs(lambda.mAntiLambda() - 1.115683) > LambdaWindow))
        return false;
      histos.fill(HIST("Selection/hLambdaNegEta"), lambda.negativeeta());
      histos.fill(HIST("Selection/hLambdaPosEta"), lambda.positiveeta());
      histos.fill(HIST("hCandidateBuilderSelection"), 7.);
      if ((TMath::Abs(lambda.negativeeta()) > LambdaDauPseudoRap) || (TMath::Abs(lambda.positiveeta()) > LambdaDauPseudoRap))
        return false;
      histos.fill(HIST("Selection/hLambdaDCANegToPV"), lambda.dcanegtopv());
      histos.fill(HIST("Selection/hLambdaDCAPosToPV"), lambda.dcapostopv());
      histos.fill(HIST("hCandidateBuilderSelection"), 8.);
      if ((TMath::Abs(lambda.dcapostopv()) < LambdaMinDCAPosToPv) || (TMath::Abs(lambda.dcanegtopv()) < LambdaMinDCANegToPv))
        return false;
      histos.fill(HIST("Selection/hLambdaRadius"), lambda.v0radius());
      histos.fill(HIST("hCandidateBuilderSelection"), 9.);
      if ((lambda.v0radius() < LambdaMinv0radius) || (lambda.v0radius() > LambdaMaxv0radius))
        return false;
      histos.fill(HIST("Selection/hLambdaDCADau"), lambda.dcaV0daughters());
      histos.fill(HIST("hCandidateBuilderSelection"), 10.);
      if (TMath::Abs(lambda.dcaV0daughters()) > LambdaMaxDCAV0Dau)
        return false;
      histos.fill(HIST("hCandidateBuilderSelection"), 11.);
    }
    // Sigma0 candidate properties
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
    auto arrMom = std::array{pVecPhotons, pVecLambda};
    float sigmamass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    float sigmarap = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);
    float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});

    histos.fill(HIST("Selection/hSigmaMass"), sigmamass);
    histos.fill(HIST("Selection/hSigmaMassWindow"), sigmamass - 1.192642);
    histos.fill(HIST("GeneralQA/h2dMassGammaVsK0S"), gamma.mGamma(), gamma.mK0Short());
    histos.fill(HIST("GeneralQA/h2dMassLambdaVsK0S"), lambda.mLambda(), lambda.mK0Short());
    histos.fill(HIST("GeneralQA/h2dMassGammaVsLambda"), gamma.mGamma(), lambda.mLambda());
    histos.fill(HIST("GeneralQA/h3dMassSigma0VsDaupTs"), gamma.pt(), lambda.pt(), sigmamass);

    if constexpr (requires { gamma.pdgCode(); } && requires { lambda.pdgCode(); }) {

      histos.fill(HIST("MC/h2dPtVsMassSigma_All"), SigmapT, sigmamass);

      // Real Gamma x Real Lambda - but not from the same sigma0/antisigma0!
      if ((gamma.pdgCode() == 22) && ((lambda.pdgCode() == 3122) || (lambda.pdgCode() == -3122)) && (gamma.motherMCPartId() != lambda.motherMCPartId())) {
        histos.fill(HIST("MC/h2dPtVsMassSigma_TrueDaughters"), SigmapT, sigmamass);
        histos.fill(HIST("MC/h2dTrueDaughtersMatrix"), lambda.pdgCodeMother(), gamma.pdgCodeMother());
      }

      // Real Gamma x fake Lambda
      if ((gamma.pdgCode() == 22) && (lambda.pdgCode() != 3122) && (lambda.pdgCode() != -3122)) {
        histos.fill(HIST("MC/h2dPtVsMassSigma_TrueGammaFakeLambda"), SigmapT, sigmamass);
        histos.fill(HIST("MC/h2dTrueGammaFakeLambdaMatrix"), lambda.pdgCodeMother(), gamma.pdgCodeMother());
      }

      // Fake Gamma x Real Lambda
      if ((gamma.pdgCode() != 22) && ((lambda.pdgCode() == 3122) || (lambda.pdgCode() == -3122))) {
        histos.fill(HIST("MC/h2dPtVsMassSigma_FakeGammaTrueLambda"), SigmapT, sigmamass);
        histos.fill(HIST("MC/h2dFakeGammaTrueLambdaMatrix"), lambda.pdgCodeMother(), gamma.pdgCodeMother());
      }

      // Fake Gamma x Fake Lambda
      if ((gamma.pdgCode() != 22) && (lambda.pdgCode() != 3122) && (lambda.pdgCode() != -3122)) {
        histos.fill(HIST("MC/h2dPtVsMassSigma_FakeDaughters"), SigmapT, sigmamass);
        histos.fill(HIST("MC/h2dFakeDaughtersMatrix"), lambda.pdgCodeMother(), gamma.pdgCodeMother());
      }
    }

    if (TMath::Abs(sigmamass - 1.192642) > Sigma0Window)
      return false;

    histos.fill(HIST("GeneralQA/h2dMassGammaVsK0SAfterMassSel"), gamma.mGamma(), gamma.mK0Short());
    histos.fill(HIST("GeneralQA/h2dMassLambdaVsK0SAfterMassSel"), lambda.mLambda(), lambda.mK0Short());
    histos.fill(HIST("GeneralQA/h2dMassGammaVsLambdaAfterMassSel"), gamma.mGamma(), lambda.mLambda());
    histos.fill(HIST("Selection/hSigmaY"), sigmarap);
    histos.fill(HIST("hCandidateBuilderSelection"), 12.);

    if (TMath::Abs(sigmarap) > SigmaMaxRap)
      return false;

    histos.fill(HIST("hCandidateBuilderSelection"), 13.);
    return true;
  }
  // Helper struct to pass v0 information
  struct {
    float mass;
    float pT;
    float Rapidity;
    float OPAngle;
    float DeltaEta;
    float DeltaPhi;
  } sigmaCandidate;

  // Fill tables with reconstructed sigma0 candidate
  template <typename TV0Object>
  void fillTables(TV0Object const& lambda, TV0Object const& gamma)
  {

    float GammaBDTScore = -1;
    float LambdaBDTScore = -1;
    float AntiLambdaBDTScore = -1;

    if constexpr (
      requires { gamma.gammaBDTScore(); } &&
      requires { lambda.lambdaBDTScore(); } &&
      requires { lambda.antiLambdaBDTScore(); }) {

      GammaBDTScore = gamma.gammaBDTScore();
      LambdaBDTScore = lambda.lambdaBDTScore();
      AntiLambdaBDTScore = lambda.antiLambdaBDTScore();
    }

    // Daughters related
    /// Photon
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
    float fPhotonPosTPCNSigma = posTrackGamma.tpcNSigmaEl();
    float fPhotonNegTPCNSigma = negTrackGamma.tpcNSigmaEl();
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
    uint32_t fPhotonPosITSClSize = posTrackGamma.itsClusterSizes();
    uint32_t fPhotonNegITSClSize = negTrackGamma.itsClusterSizes();
    uint8_t fPhotonV0Type = gamma.v0Type();

    // Lambda
    auto posTrackLambda = lambda.template posTrackExtra_as<dauTracks>();
    auto negTrackLambda = lambda.template negTrackExtra_as<dauTracks>();

    float fLambdaPt = lambda.pt();
    float fLambdaMass = lambda.mLambda();
    float fAntiLambdaMass = lambda.mAntiLambda();
    float fLambdaQt = lambda.qtarm();
    float fLambdaAlpha = lambda.alpha();
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

    float fLambdaPrTOFNSigma = lambda.tofNSigmaLaPr();
    float fLambdaPiTOFNSigma = lambda.tofNSigmaLaPi();
    float fALambdaPrTOFNSigma = lambda.tofNSigmaALaPr();
    float fALambdaPiTOFNSigma = lambda.tofNSigmaALaPi();

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
    uint32_t fLambdaPosITSClSize = posTrackLambda.itsClusterSizes();
    uint32_t fLambdaNegITSClSize = negTrackLambda.itsClusterSizes();
    uint8_t fLambdaV0Type = lambda.v0Type();

    // Sigma0 candidate properties
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
    auto arrMom = std::array{pVecPhotons, pVecLambda};
    TVector3 v1(gamma.px(), gamma.py(), gamma.pz());
    TVector3 v2(lambda.px(), lambda.py(), lambda.pz());

    sigmaCandidate.mass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    sigmaCandidate.pT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
    sigmaCandidate.Rapidity = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);
    sigmaCandidate.OPAngle = v1.Angle(v2);
    sigmaCandidate.DeltaEta = fLambdaEta - fPhotonEta;
    sigmaCandidate.DeltaPhi = fLambdaPhi - fPhotonPhi;

    // Sigma related
    float fSigmapT = sigmaCandidate.pT;
    float fSigmaMass = sigmaCandidate.mass;
    float fSigmaRap = sigmaCandidate.Rapidity;
    float fSigmaOPAngle = sigmaCandidate.OPAngle;
    float fSigmaDeltaEta = sigmaCandidate.DeltaEta;
    float fSigmaDeltaPhi = sigmaCandidate.DeltaPhi;

    // Filling TTree for ML analysis
    sigma0cores(fSigmapT, fSigmaMass, fSigmaRap, fSigmaOPAngle, fSigmaDeltaEta, fSigmaDeltaPhi);

    sigmaPhotonExtras(fPhotonPt, fPhotonMass, fPhotonQt, fPhotonAlpha, fPhotonRadius,
                      fPhotonCosPA, fPhotonDCADau, fPhotonDCANegPV, fPhotonDCAPosPV, fPhotonZconv,
                      fPhotonEta, fPhotonY, fPhotonPhi, fPhotonPosTPCNSigma, fPhotonNegTPCNSigma, fPhotonPosTPCCrossedRows,
                      fPhotonNegTPCCrossedRows, fPhotonPosPt, fPhotonNegPt, fPhotonPosEta,
                      fPhotonNegEta, fPhotonPosY, fPhotonNegY, fPhotonPsiPair,
                      fPhotonPosITSCls, fPhotonNegITSCls, fPhotonPosITSClSize, fPhotonNegITSClSize,
                      fPhotonV0Type, GammaBDTScore);

    sigmaLambdaExtras(fLambdaPt, fLambdaMass, fAntiLambdaMass, fLambdaQt, fLambdaAlpha,
                      fLambdaRadius, fLambdaCosPA, fLambdaDCADau, fLambdaDCANegPV,
                      fLambdaDCAPosPV, fLambdaEta, fLambdaY, fLambdaPhi, fLambdaPosPrTPCNSigma,
                      fLambdaPosPiTPCNSigma, fLambdaNegPrTPCNSigma, fLambdaNegPiTPCNSigma,
                      fLambdaPrTOFNSigma, fLambdaPiTOFNSigma, fALambdaPrTOFNSigma, fALambdaPiTOFNSigma,
                      fLambdaPosTPCCrossedRows, fLambdaNegTPCCrossedRows, fLambdaPosPt, fLambdaNegPt, fLambdaPosEta,
                      fLambdaNegEta, fLambdaPosPrY, fLambdaPosPiY, fLambdaNegPrY, fLambdaNegPiY,
                      fLambdaPosITSCls, fLambdaNegITSCls, fLambdaPosITSClSize, fLambdaNegITSClSize,
                      fLambdaV0Type, LambdaBDTScore, AntiLambdaBDTScore);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents> const& collisions, V0DerivedMCDatas const& V0s)
  {
    for (const auto& coll : collisions) {
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = coll.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionMCDerived, collIdx);

      // V0 table sliced
      for (auto& gamma : V0Table_thisCollision) { // selecting photons from Sigma0
        float centrality = coll.centFT0C();

        // Auxiliary histograms:
        if (gamma.pdgCode() == 22) {
          float GammaY = TMath::Abs(RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma));

          if (GammaY < 0.5) {                                                                                                                // rapidity selection
            histos.fill(HIST("MC/h2dPtVsCentrality_GammaBeforeSel"), centrality, gamma.pt());                                                // isgamma
            histos.fill(HIST("MC/h2dGammaPtResolution"), gamma.pt(), gamma.pt() - RecoDecay::pt(array{gamma.pxMC(), gamma.pyMC()}));         // pT resolution

            if (gamma.pdgCodeMother() == 3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_GammaSigma0"), centrality, gamma.pt()); // isgamma from sigma
            }
            if (gamma.pdgCodeMother() == -3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_GammaAntiSigma0"), centrality, gamma.pt()); // isgamma from sigma
            }
          }
        }
        if (gamma.pdgCode() == 3122) { // Is Lambda
          float LambdaY = TMath::Abs(RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassLambda));
          if (LambdaY < 0.5) { // rapidity selection
            histos.fill(HIST("MC/h2dPtVsCentrality_LambdaBeforeSel"), centrality, gamma.pt());
            histos.fill(HIST("MC/h2dLambdaPtResolution"), gamma.pt(), gamma.pt() - RecoDecay::pt(array{gamma.pxMC(), gamma.pyMC()})); // pT resolution
            if (gamma.pdgCodeMother() == 3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_LambdaSigma0"), centrality, gamma.pt());
            }
          }
        }
        if (gamma.pdgCode() == -3122) { // Is AntiLambda
          float AntiLambdaY = TMath::Abs(RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassLambda));
          if (AntiLambdaY < 0.5) { // rapidity selection
            histos.fill(HIST("MC/h2dPtVsCentrality_AntiLambdaBeforeSel"), centrality, gamma.pt());
            if (gamma.pdgCodeMother() == -3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_LambdaAntiSigma0"), centrality, gamma.pt()); // isantilambda from antisigma
            }
          }
        }

        for (auto& lambda : V0Table_thisCollision) { // selecting lambdas from Sigma0
          // Sigma0 candidate properties
          std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
          std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
          auto arrMom = std::array{pVecPhotons, pVecLambda};
          float SigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
          float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
          float SigmaY = TMath::Abs(RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0));

          if ((gamma.pdgCode() == 22) && (gamma.pdgCodeMother() == 3212) && (lambda.pdgCode() == 3122) && (lambda.pdgCodeMother() == 3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()) && (SigmaY < 0.5)) {
            histos.fill(HIST("MC/h2dPtVsCentrality_Sigma0BeforeSel"), centrality, RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()}));
            histos.fill(HIST("MC/h2dSigmaPtVsLambdaPt"), SigmapT, lambda.pt());
            histos.fill(HIST("MC/h2dSigmaPtVsGammaPt"), SigmapT, gamma.pt());
          }
          if ((gamma.pdgCode() == 22) && (gamma.pdgCodeMother() == -3212) && (lambda.pdgCode() == -3122) && (lambda.pdgCodeMother() == -3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()) && (SigmaY < 0.5))
            histos.fill(HIST("MC/h2dPtVsCentrality_AntiSigma0BeforeSel"), centrality, SigmapT);

          if (!processSigmaCandidate(lambda, gamma)) // basic selection
            continue;

          bool fIsSigma = false;
          bool fIsAntiSigma = false;
          bool fIsPhotonPrimary = gamma.isPhysicalPrimary();
          int PhotonCandPDGCode = gamma.pdgCode();
          int PhotonCandPDGCodeMother = gamma.pdgCodeMother();
          bool fIsLambdaPrimary = lambda.isPhysicalPrimary();
          int LambdaCandPDGCode = lambda.pdgCode();
          int LambdaCandPDGCodeMother = lambda.pdgCodeMother();

          if ((gamma.pdgCode() == 22) && (gamma.pdgCodeMother() == 3212) && (lambda.pdgCode() == 3122) && (lambda.pdgCodeMother() == 3212) && (gamma.motherMCPartId() == lambda.motherMCPartId())) {
            fIsSigma = true;
            histos.fill(HIST("MC/h2dPtVsCentrality_Sigma0AfterSel"), centrality, RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()}));
          }
          if ((gamma.pdgCode() == 22) && (gamma.pdgCodeMother() == -3212) && (lambda.pdgCode() == -3122) && (lambda.pdgCodeMother() == -3212) && (gamma.motherMCPartId() == lambda.motherMCPartId())) {
            fIsAntiSigma = true;
            histos.fill(HIST("MC/h2dPtVsCentrality_AntiSigma0AfterSel"), centrality, RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()}));
            // TH3D Mass histogram
          }
          sigma0mccores(fIsSigma, fIsAntiSigma,
                        PhotonCandPDGCode, PhotonCandPDGCodeMother, fIsPhotonPrimary,
                        LambdaCandPDGCode, LambdaCandPDGCodeMother, fIsLambdaPrimary);

          // QA histograms
          // Signal only (sigma0+antisigma0)
          if (fIsSigma || fIsAntiSigma)
            histos.fill(HIST("MC/h2dPtVsMassSigma_SignalOnly"), SigmapT, SigmaMass);
        }
      }
    }
  }

  void processSTDSelection(soa::Join<aod::StraCollisions, aod::StraCents> const& collisions, V0StandardDerivedDatas const& V0s, dauTracks const&)
  {
    for (const auto& coll : collisions) {
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = coll.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionSTDDerived, collIdx);

      histos.fill(HIST("hEventCentrality"), coll.centFT0C());
      sigma0Coll(coll.posX(), coll.posY(), coll.posZ(), coll.centFT0M(), coll.centFT0A(), coll.centFT0C(), coll.centFV0A());

      // V0 table sliced
      for (auto& gamma : V0Table_thisCollision) {    // selecting photons from Sigma0
        for (auto& lambda : V0Table_thisCollision) { // selecting lambdas from Sigma0
          std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
          std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
          auto arrMom = std::array{pVecPhotons, pVecLambda};
          float SigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
          float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
          float SigmaY = TMath::Abs(RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0));
          histos.fill(HIST("h3dMassSigmasBeforeSel"), coll.centFT0C(), SigmapT, SigmaMass);

          if (!processSigmaCandidate(lambda, gamma)) // applying selection for reconstruction
            continue;

          histos.fill(HIST("h3dMassSigmasAfterSel"), coll.centFT0C(), SigmapT, SigmaMass);

          sigma0CollRefs(collIdx);
          fillTables(lambda, gamma); // filling tables with accepted candidates

          nSigmaCandidates++;
          if (nSigmaCandidates % 5000 == 0) {
            LOG(info) << "Sigma0 Candidates built: " << nSigmaCandidates;
          }
        }
      }
    }
  }

  void processMLSelection(soa::Join<aod::StraCollisions, aod::StraCents> const& collisions, V0MLDerivedDatas const& V0s, dauTracks const&)
  {
    for (const auto& coll : collisions) {
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = coll.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionMLDerived, collIdx);

      histos.fill(HIST("hEventCentrality"), coll.centFT0C());
      sigma0Coll(coll.posX(), coll.posY(), coll.posZ(), coll.centFT0M(), coll.centFT0A(), coll.centFT0C(), coll.centFV0A());

      // V0 table sliced
      for (auto& gamma : V0Table_thisCollision) {    // selecting photons from Sigma0
        for (auto& lambda : V0Table_thisCollision) { // selecting lambdas from Sigma0
          if (!processSigmaCandidate(lambda, gamma))
            continue;

          nSigmaCandidates++;
          if (nSigmaCandidates % 5000 == 0) {
            LOG(info) << "Sigma0 Candidates built: " << nSigmaCandidates;
          }
          sigma0CollRefs(collIdx);
          fillTables(lambda, gamma); // filling tables with accepted candidates
        }
      }
    }
  }
  PROCESS_SWITCH(sigma0builder, processMonteCarlo, "Fill sigma0 MC table", false);
  PROCESS_SWITCH(sigma0builder, processSTDSelection, "Select gammas and lambdas with standard cuts", true);
  PROCESS_SWITCH(sigma0builder, processMLSelection, "Select gammas and lambdas with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigma0builder>(cfgc)};
}
