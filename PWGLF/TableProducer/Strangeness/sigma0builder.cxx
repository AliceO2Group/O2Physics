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
// using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCDatas>;
using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0StandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;

struct sigma0builder {
  SliceCache cache;

  Produces<aod::Sigma0Cores> sigma0cores;             // save sigma0 candidates for analysis
  Produces<aod::SigmaPhotonExtras> sigmaPhotonExtras; // save sigma0 candidates for analysis
  Produces<aod::SigmaLambdaExtras> sigmaLambdaExtras; // save sigma0 candidates for analysis
  Produces<aod::SigmaMCCores> sigma0mccores;

  // For manual sliceBy
  Preslice<V0DerivedMCDatas> perCollisionMCDerived = o2::aod::v0data::straCollisionId;
  Preslice<V0StandardDerivedDatas> perCollisionSTDDerived = o2::aod::v0data::straCollisionId;

  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event selection
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};

  struct : ConfigurableGroup {
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", true, "require events with at least one ITS-TPC track"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", false, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeVzDep{"requireNoCollInTimeRangeVzDep", false, "reject collisions corrupted by the cannibalism, with other collisions with pvZ of drifting TPC tracks from past/future collisions within 2.5 cm the current pvZ"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF"};
    Configurable<bool> requireINEL0{"requireINEL0", false, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};
    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};
    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
  } eventSelections;

  // For ML Selection
  Configurable<bool> useMLScores{"useMLScores", false, "use ML scores to select candidates"};
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

  //// Extras:
  Configurable<bool> doPi0QA{"doPi0QA", true, "Flag to fill QA histos for pi0 rejection study."};
  Configurable<float> Pi0PhotonMinDCADauToPv{"Pi0PhotonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
  Configurable<float> Pi0PhotonMaxDCAV0Dau{"Pi0PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<int> Pi0PhotonMinTPCCrossedRows{"Pi0PhotonMinTPCCrossedRows", 0, "Min daughter TPC Crossed Rows"};
  Configurable<int> Pi0PhotonMaxTPCNSigmas{"Pi0PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
  Configurable<float> Pi0PhotonMaxEta{"Pi0PhotonMaxEta", 0.8, "Max photon rapidity"};
  Configurable<float> Pi0PhotonMinRadius{"Pi0PhotonMinRadius", 3.0, "Min photon conversion radius (cm)"};
  Configurable<float> Pi0PhotonMaxRadius{"Pi0PhotonMaxRadius", 115, "Max photon conversion radius (cm)"};
  Configurable<float> Pi0PhotonMaxQt{"Pi0PhotonMaxQt", 0.05, "Max photon qt value (AP plot) (GeV/c)"};
  Configurable<float> Pi0PhotonMaxAlpha{"Pi0PhotonMaxAlpha", 0.95, "Max photon alpha absolute value (AP plot)"};
  Configurable<float> Pi0PhotonMinV0cospa{"Pi0PhotonMinV0cospa", 0.80, "Min V0 CosPA"};
  Configurable<float> Pi0PhotonMaxMass{"Pi0PhotonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};

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
    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{20, -0.5f, +18.5f}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "posZ cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(7, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTOFmatched");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(10, "kIsVertexTRDmatched");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(11, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeStd");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeStrict");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(14, "kNoCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(15, "kNoCollInRofStd");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(16, "kNoCollInRofStrict");
    if (doPPAnalysis) {
      histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(17, "INEL>0");
      histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(18, "INEL>1");
    } else {
      histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(17, "Below min occup.");
      histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(18, "Above max occup.");
    }

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
    histos.add("GeneralQA/h2dPtVsMassPi0BeforeSel_Candidates", "h2dPtVsMassPi0BeforeSel_Candidates", kTH2D, {axisPt, {500, 0.08f, 0.18f}});
    histos.add("GeneralQA/h2dPtVsMassPi0AfterSel_Candidates", "h2dPtVsMassPi0AfterSel_Candidates", kTH2D, {axisPt, {500, 0.08f, 0.18f}});

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

    // For Pi0 QA
    histos.add("MC/h2dPtVsMassPi0BeforeSel_SignalOnly", "h2dPtVsMassPi0BeforeSel_SignalOnly", kTH2D, {axisPt, {500, 0.08f, 0.18f}});
    histos.add("MC/h2dPtVsMassPi0AfterSel_SignalOnly", "h2dPtVsMassPi0AfterSel_SignalOnly", kTH2D, {axisPt, {500, 0.08f, 0.18f}});

    histos.add("h3dMassSigmasBeforeSel", "h3dMassSigmasBeforeSel", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("h3dMassSigmasAfterSel", "h3dMassSigmasAfterSel", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
  }

  template <typename TCollision>
  bool IsEventAccepted(TCollision collision, bool fillHists)
  // check whether the collision passes our collision selections
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);
    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
    if (doPPAnalysis) { // we are in pp
      if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 16 /* INEL > 0 */);
      if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 17 /* INEL > 1 */);
    } else { // we are in Pb-Pb
      float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 16 /* Below min occupancy */);
      if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 17 /* Above max occupancy */);
    }
    return true;
  }

  template <typename TV0Object>
  void runPi0QA(TV0Object const& gamma1, TV0Object const& gamma2)
  {

    // Check if both V0s are made of the same tracks
    if (gamma1.posTrackExtraId() == gamma2.posTrackExtraId() ||
        gamma1.negTrackExtraId() == gamma2.negTrackExtraId() ||
        gamma1.posTrackExtraId() == gamma2.negTrackExtraId() ||
        gamma1.negTrackExtraId() == gamma2.posTrackExtraId()) {
      return;
    }

    // Calculate pi0 properties
    std::array<float, 3> pVecGamma1{gamma1.px(), gamma1.py(), gamma1.pz()};
    std::array<float, 3> pVecGamma2{gamma2.px(), gamma2.py(), gamma2.pz()};
    std::array arrpi0{pVecGamma1, pVecGamma2};
    float pi0Mass = RecoDecay::m(arrpi0, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassPhoton});
    float pi0Pt = RecoDecay::pt(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py()});
    float pi0Y = RecoDecay::y(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py(), gamma1.pz() + gamma2.pz()}, o2::constants::physics::MassPi0);

    // MC-specific variables
    bool fIsPi0 = false, fIsMC = false;

    // Check if MC data and populate fIsMC, fIsPi0
    if constexpr (requires { gamma1.motherMCPartId(); gamma2.motherMCPartId(); }) {
      if (gamma1.has_v0MCCore() && gamma2.has_v0MCCore()) {
        fIsMC = true;
        auto gamma1MC = gamma1.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
        auto gamma2MC = gamma2.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

        if (gamma1MC.pdgCode() == 22 && gamma2MC.pdgCode() == 22 &&
            gamma1MC.pdgCodeMother() == 111 && gamma2MC.pdgCodeMother() == 111 &&
            gamma1.motherMCPartId() == gamma2.motherMCPartId()) {
          fIsPi0 = true;
          histos.fill(HIST("MC/h2dPtVsMassPi0BeforeSel_SignalOnly"), pi0Pt, pi0Mass);
        }
      }
    } else {
      histos.fill(HIST("GeneralQA/h2dPtVsMassPi0BeforeSel_Candidates"), pi0Pt, pi0Mass);
    }

    // Photon-specific selections
    auto posTrackGamma1 = gamma1.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma1 = gamma1.template negTrackExtra_as<dauTracks>();
    auto posTrackGamma2 = gamma2.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma2 = gamma2.template negTrackExtra_as<dauTracks>();

    // Gamma1 Selection
    bool passedTPCGamma1 = (posTrackGamma1.tpcNSigmaEl() == -999.f || TMath::Abs(posTrackGamma1.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas) &&
                           (negTrackGamma1.tpcNSigmaEl() == -999.f || TMath::Abs(negTrackGamma1.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas);

    if (TMath::Abs(gamma1.mGamma()) > Pi0PhotonMaxMass ||
        gamma1.qtarm() >= Pi0PhotonMaxQt ||
        TMath::Abs(gamma1.alpha()) >= Pi0PhotonMaxAlpha ||
        TMath::Abs(gamma1.dcapostopv()) < Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma1.dcanegtopv()) < Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma1.dcaV0daughters()) > Pi0PhotonMaxDCAV0Dau ||
        TMath::Abs(gamma1.negativeeta()) >= Pi0PhotonMaxEta ||
        TMath::Abs(gamma1.positiveeta()) >= Pi0PhotonMaxEta ||
        gamma1.v0cosPA() <= Pi0PhotonMinV0cospa ||
        gamma1.v0radius() <= Pi0PhotonMinRadius ||
        gamma1.v0radius() >= Pi0PhotonMaxRadius ||
        posTrackGamma1.tpcCrossedRows() < Pi0PhotonMinTPCCrossedRows ||
        negTrackGamma1.tpcCrossedRows() < Pi0PhotonMinTPCCrossedRows ||
        !passedTPCGamma1) {
      return;
    }

    // Gamma2 Selection
    bool passedTPCGamma2 = (posTrackGamma2.tpcNSigmaEl() == -999.f || TMath::Abs(posTrackGamma2.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas) &&
                           (negTrackGamma2.tpcNSigmaEl() == -999.f || TMath::Abs(negTrackGamma2.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas);

    if (TMath::Abs(gamma2.mGamma()) > Pi0PhotonMaxMass ||
        gamma2.qtarm() >= Pi0PhotonMaxQt ||
        TMath::Abs(gamma2.alpha()) >= Pi0PhotonMaxAlpha ||
        TMath::Abs(gamma2.dcapostopv()) < Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma2.dcanegtopv()) < Pi0PhotonMinDCADauToPv ||
        TMath::Abs(gamma2.dcaV0daughters()) > Pi0PhotonMaxDCAV0Dau ||
        TMath::Abs(gamma2.negativeeta()) >= Pi0PhotonMaxEta ||
        TMath::Abs(gamma2.positiveeta()) >= Pi0PhotonMaxEta ||
        gamma2.v0cosPA() <= Pi0PhotonMinV0cospa ||
        gamma2.v0radius() <= Pi0PhotonMinRadius ||
        gamma2.v0radius() >= Pi0PhotonMaxRadius ||
        posTrackGamma2.tpcCrossedRows() < Pi0PhotonMinTPCCrossedRows ||
        negTrackGamma2.tpcCrossedRows() < Pi0PhotonMinTPCCrossedRows ||
        !passedTPCGamma2) {
      return;
    }

    // Pi0-specific selections:
    if (TMath::Abs(pi0Y) > 0.5) {
      return;
    }

    // Fill histograms
    if (fIsMC) {
      if (fIsPi0)
        histos.fill(HIST("MC/h2dPtVsMassPi0AfterSel_SignalOnly"), pi0Pt, pi0Mass);
    } else {
      histos.fill(HIST("GeneralQA/h2dPtVsMassPi0AfterSel_Candidates"), pi0Pt, pi0Mass);
    }
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

    if (useMLScores) {
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

      // Lambda basic selection criteria:
      histos.fill(HIST("hCandidateBuilderSelection"), 6.);
      histos.fill(HIST("Selection/hLambdaMass"), lambda.mLambda());
      histos.fill(HIST("Selection/hAntiLambdaMass"), lambda.mAntiLambda());
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

  // Fill tables with reconstructed sigma0 candidate
  template <typename TV0Object, typename TCollision>
  void fillTables(TV0Object const& lambda, TV0Object const& gamma, TCollision const& coll)
  {

    float GammaBDTScore = gamma.gammaBDTScore();
    float LambdaBDTScore = lambda.lambdaBDTScore();
    float AntiLambdaBDTScore = lambda.antiLambdaBDTScore();

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

    // Lambda
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
    float fLambdaPosITSChi2PerNcl = posTrackLambda.itsChi2PerNcl();
    float fLambdaNegITSChi2PerNcl = negTrackLambda.itsChi2PerNcl();
    uint8_t fLambdaV0Type = lambda.v0Type();

    // Sigma0 candidate properties
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
    auto arrMom = std::array{pVecPhotons, pVecLambda};
    TVector3 v1(gamma.px(), gamma.py(), gamma.pz());
    TVector3 v2(lambda.px(), lambda.py(), lambda.pz());

    // Sigma related
    float fSigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
    float fSigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    float fSigmaRap = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);
    float fSigmaOPAngle = v1.Angle(v2);
    float fSigmaCentrality = coll.centFT0C();
    float fSigmaTimeStamp = coll.timestamp();
    float fSigmaRunNumber = coll.runNumber();

    // Filling TTree for ML analysis
    sigma0cores(fSigmapT, fSigmaMass, fSigmaRap, fSigmaOPAngle, fSigmaCentrality, fSigmaRunNumber, fSigmaTimeStamp);

    sigmaPhotonExtras(fPhotonPt, fPhotonMass, fPhotonQt, fPhotonAlpha, fPhotonRadius,
                      fPhotonCosPA, fPhotonDCADau, fPhotonDCANegPV, fPhotonDCAPosPV, fPhotonZconv,
                      fPhotonEta, fPhotonY, fPhotonPhi, fPhotonPosTPCNSigmaEl, fPhotonNegTPCNSigmaEl, fPhotonPosTPCNSigmaPi, fPhotonNegTPCNSigmaPi, fPhotonPosTPCCrossedRows,
                      fPhotonNegTPCCrossedRows, fPhotonPosPt, fPhotonNegPt, fPhotonPosEta,
                      fPhotonNegEta, fPhotonPosY, fPhotonNegY, fPhotonPsiPair,
                      fPhotonPosITSCls, fPhotonNegITSCls, fPhotonPosITSChi2PerNcl, fPhotonNegITSChi2PerNcl,
                      fPhotonV0Type, GammaBDTScore);

    sigmaLambdaExtras(fLambdaPt, fLambdaMass, fAntiLambdaMass, fLambdaQt, fLambdaAlpha, fLambdaLifeTime,
                      fLambdaRadius, fLambdaCosPA, fLambdaDCADau, fLambdaDCANegPV,
                      fLambdaDCAPosPV, fLambdaEta, fLambdaY, fLambdaPhi, fLambdaPosPrTPCNSigma,
                      fLambdaPosPiTPCNSigma, fLambdaNegPrTPCNSigma, fLambdaNegPiTPCNSigma,
                      fLambdaPrTOFNSigma, fLambdaPiTOFNSigma, fALambdaPrTOFNSigma, fALambdaPiTOFNSigma,
                      fLambdaPosTPCCrossedRows, fLambdaNegTPCCrossedRows, fLambdaPosPt, fLambdaNegPt, fLambdaPosEta,
                      fLambdaNegEta, fLambdaPosPrY, fLambdaPosPiY, fLambdaNegPrY, fLambdaNegPiY,
                      fLambdaPosITSCls, fLambdaNegITSCls, fLambdaPosITSChi2PerNcl, fLambdaNegITSChi2PerNcl,
                      fLambdaV0Type, LambdaBDTScore, AntiLambdaBDTScore);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0DerivedMCDatas const& V0s, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    for (const auto& coll : collisions) {
      if (!IsEventAccepted(coll, true)) {
        continue;
      }
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = coll.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionMCDerived, collIdx);

      histos.fill(HIST("hEventCentrality"), coll.centFT0C());
      // V0 table sliced
      for (auto& gamma : V0Table_thisCollision) { // selecting photons from Sigma0
        float centrality = coll.centFT0C();

        if (!gamma.has_v0MCCore())
          continue;

        auto gammaMC = gamma.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

        // Auxiliary histograms:
        if (gammaMC.pdgCode() == 22) {
          float GammaY = TMath::Abs(RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma));

          if (GammaY < 0.5) {                                                                                                                // rapidity selection
            histos.fill(HIST("MC/h2dPtVsCentrality_GammaBeforeSel"), centrality, gamma.pt());                                                // isgamma
            histos.fill(HIST("MC/h2dGammaPtResolution"), gamma.pt(), gamma.pt() - RecoDecay::pt(array{gammaMC.pxMC(), gammaMC.pyMC()}));     // pT resolution

            if (gammaMC.pdgCodeMother() == 3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_GammaSigma0"), centrality, gamma.pt()); // isgamma from sigma
            }
            if (gammaMC.pdgCodeMother() == -3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_GammaAntiSigma0"), centrality, gamma.pt()); // isgamma from sigma
            }
          }
        }
        if (gammaMC.pdgCode() == 3122) { // Is Lambda
          float LambdaY = TMath::Abs(RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassLambda));
          if (LambdaY < 0.5) { // rapidity selection
            histos.fill(HIST("MC/h2dPtVsCentrality_LambdaBeforeSel"), centrality, gamma.pt());
            histos.fill(HIST("MC/h2dLambdaPtResolution"), gamma.pt(), gamma.pt() - RecoDecay::pt(array{gammaMC.pxMC(), gammaMC.pyMC()})); // pT resolution
            if (gammaMC.pdgCodeMother() == 3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_LambdaSigma0"), centrality, gamma.pt());
            }
          }
        }
        if (gammaMC.pdgCode() == -3122) { // Is AntiLambda
          float AntiLambdaY = TMath::Abs(RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassLambda));
          if (AntiLambdaY < 0.5) { // rapidity selection
            histos.fill(HIST("MC/h2dPtVsCentrality_AntiLambdaBeforeSel"), centrality, gamma.pt());
            if (gammaMC.pdgCodeMother() == -3212) {
              histos.fill(HIST("MC/h2dPtVsCentrality_LambdaAntiSigma0"), centrality, gamma.pt()); // isantilambda from antisigma
            }
          }
        }

        for (auto& lambda : V0Table_thisCollision) { // selecting lambdas from Sigma0
          if (!lambda.has_v0MCCore())
            continue;

          if (lambda.v0Type() != 1) { // safeguard to avoid TPC-only photons
            continue;
          }

          auto lambdaMC = lambda.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

          if (doPi0QA)                               // Pi0 QA study
            runPi0QA(gamma, lambda);

          // Sigma0 candidate properties
          std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
          std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
          auto arrMom = std::array{pVecPhotons, pVecLambda};
          float SigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
          float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
          float SigmaY = TMath::Abs(RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0));

          if ((gammaMC.pdgCode() == 22) && (gammaMC.pdgCodeMother() == 3212) && (lambdaMC.pdgCode() == 3122) && (lambdaMC.pdgCodeMother() == 3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()) && (SigmaY < 0.5)) {
            histos.fill(HIST("MC/h2dPtVsCentrality_Sigma0BeforeSel"), centrality, RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()}));
            histos.fill(HIST("MC/h2dSigmaPtVsLambdaPt"), SigmapT, lambda.pt());
            histos.fill(HIST("MC/h2dSigmaPtVsGammaPt"), SigmapT, gamma.pt());
          }
          if ((gammaMC.pdgCode() == 22) && (gammaMC.pdgCodeMother() == -3212) && (lambdaMC.pdgCode() == -3122) && (lambdaMC.pdgCodeMother() == -3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()) && (SigmaY < 0.5))
            histos.fill(HIST("MC/h2dPtVsCentrality_AntiSigma0BeforeSel"), centrality, SigmapT);

          if (!processSigmaCandidate(lambda, gamma)) // basic selection
            continue;

          bool fIsSigma = false;
          bool fIsAntiSigma = false;
          float SigmaMCpT = RecoDecay::pt(array{gammaMC.pxMC() + lambdaMC.pxMC(), gammaMC.pyMC() + lambdaMC.pyMC()});
          bool fIsPhotonPrimary = gammaMC.isPhysicalPrimary();
          int PhotonCandPDGCode = gammaMC.pdgCode();
          int PhotonCandPDGCodeMother = gammaMC.pdgCodeMother();
          float PhotonMCpT = RecoDecay::pt(array{gammaMC.pxMC(), gammaMC.pyMC()});
          bool fIsLambdaPrimary = lambdaMC.isPhysicalPrimary();
          int LambdaCandPDGCode = lambdaMC.pdgCode();
          int LambdaCandPDGCodeMother = lambdaMC.pdgCodeMother();
          float LambdaMCpT = RecoDecay::pt(array{lambdaMC.pxMC(), lambdaMC.pyMC()});

          if ((gammaMC.pdgCode() == 22) && (gammaMC.pdgCodeMother() == 3212) && (lambdaMC.pdgCode() == 3122) && (lambdaMC.pdgCodeMother() == 3212) && (gamma.motherMCPartId() == lambda.motherMCPartId())) {
            fIsSigma = true;
            histos.fill(HIST("MC/h2dPtVsCentrality_Sigma0AfterSel"), centrality, RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()}));
          }
          if ((gammaMC.pdgCode() == 22) && (gammaMC.pdgCodeMother() == -3212) && (lambdaMC.pdgCode() == -3122) && (lambdaMC.pdgCodeMother() == -3212) && (gamma.motherMCPartId() == lambda.motherMCPartId())) {
            fIsAntiSigma = true;
            histos.fill(HIST("MC/h2dPtVsCentrality_AntiSigma0AfterSel"), centrality, RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()}));
            // TH3D Mass histogram
          }
          sigma0mccores(fIsSigma, fIsAntiSigma, SigmaMCpT,
                        PhotonCandPDGCode, PhotonCandPDGCodeMother, fIsPhotonPrimary, PhotonMCpT,
                        LambdaCandPDGCode, LambdaCandPDGCodeMother, fIsLambdaPrimary, LambdaMCpT);

          fillTables(lambda, gamma, coll); // filling tables with accepted candidates

          nSigmaCandidates++;
          if (nSigmaCandidates % 5000 == 0) {
            LOG(info) << "Sigma0 Candidates built: " << nSigmaCandidates;
          }

          // QA histograms
          // Signal only (sigma0+antisigma0)
          if (fIsSigma || fIsAntiSigma)
            histos.fill(HIST("MC/h2dPtVsMassSigma_SignalOnly"), SigmapT, SigmaMass);
        }
      }
    }
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0StandardDerivedDatas const& V0s, dauTracks const&)
  {
    for (const auto& coll : collisions) {
      if (!IsEventAccepted(coll, true)) {
        continue;
      }
      // Do analysis with collision-grouped V0s, retain full collision information
      const uint64_t collIdx = coll.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollisionSTDDerived, collIdx);

      histos.fill(HIST("hEventCentrality"), coll.centFT0C());
      // V0 table sliced
      for (auto& gamma : V0Table_thisCollision) {    // selecting photons from Sigma0
        for (auto& lambda : V0Table_thisCollision) { // selecting lambdas from Sigma0
          if (doPi0QA)                               // Pi0 QA study
            runPi0QA(gamma, lambda);

          if (lambda.v0Type() != 1) { // safeguard to avoid TPC-only photons
            continue;
          }

          // Sigma0 candidate properties
          std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
          std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
          auto arrMom = std::array{pVecPhotons, pVecLambda};
          float SigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
          float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
          // float SigmaY = TMath::Abs(RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0));
          histos.fill(HIST("h3dMassSigmasBeforeSel"), coll.centFT0C(), SigmapT, SigmaMass);

          if (!processSigmaCandidate(lambda, gamma)) // applying selection for reconstruction
            continue;

          histos.fill(HIST("h3dMassSigmasAfterSel"), coll.centFT0C(), SigmapT, SigmaMass);

          fillTables(lambda, gamma, coll); // filling tables with accepted candidates

          nSigmaCandidates++;
          if (nSigmaCandidates % 5000 == 0) {
            LOG(info) << "Sigma0 Candidates built: " << nSigmaCandidates;
          }
        }
      }
    }
  }

  PROCESS_SWITCH(sigma0builder, processMonteCarlo, "process as if MC data", false);
  PROCESS_SWITCH(sigma0builder, processRealData, "process as if real data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigma0builder>(cfgc)};
}
