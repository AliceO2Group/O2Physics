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
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0StandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;

struct sigma0builder {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  // SliceCache cache;

  Produces<aod::Sigma0Cores> sigma0cores;             // save sigma0 candidates for analysis
  Produces<aod::SigmaPhotonExtras> sigmaPhotonExtras; // save sigma0 candidates for analysis
  Produces<aod::SigmaLambdaExtras> sigmaLambdaExtras; // save sigma0 candidates for analysis
  Produces<aod::SigmaMCCores> sigma0mccores;

  // For manual sliceBy
  // PresliceUnsorted<V0DerivedMCDatas> perCollisionMCDerived = o2::aod::v0data::straCollisionId;
  // PresliceUnsorted<V0StandardDerivedDatas> perCollisionSTDDerived = o2::aod::v0data::straCollisionId;
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

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
  Configurable<bool> fillBkgQAhistos{"fillBkgQAhistos", false, "if true, fill MC QA histograms for Bkg study"};
  Configurable<bool> doPi0QA{"doPi0QA", true, "Flag to fill QA histos for pi0 rejection study."};
  Configurable<bool> doAssocStudy{"doAssocStudy", false, "Do v0 to collision association study."};

  // Event level
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<bool> fIRCrashOnNull{"fIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

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
    Configurable<bool> useEvtSelInDenomEff{"useEvtSelInDenomEff", false, "Consider event selections in the recoed <-> gen collision association for the denominator (or numerator) of the acc. x eff. (or signal loss)?"};
    Configurable<bool> applyZVtxSelOnMCPV{"applyZVtxSelOnMCPV", false, "Apply Z-vtx cut on the PV of the generated collision?"};
    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};

    // fast check on interaction rate
    Configurable<float> minIR{"minIR", -1, "minimum IR collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR collisions"};

  } eventSelections;

  // For ML Selection
  Configurable<bool> useMLScores{"useMLScores", false, "use ML scores to select candidates"};
  Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
  Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
  Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};

  // For standard approach:
  //// Lambda criteria:
  Configurable<float> V0Rapidity{"V0Rapidity", 0.8, "v0 rapidity"};

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
  ConfigurableAxis axisIRBinning{"axisIRBinning", {150, 0, 1500}, "Binning for the interaction rate (kHz)"};

  int nSigmaCandidates = 0;
  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1D, {{21, -0.5f, +20.5f}});
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
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(19, "Below min IR");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(20, "Above max IR");

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

    if (fGetIR) {
      histos.add("GeneralQA/hRunNumberNegativeIR", "", kTH1D, {{1, 0., 1.}});
      histos.add("GeneralQA/hInteractionRate", "hInteractionRate", kTH1F, {axisIRBinning});
      histos.add("GeneralQA/hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2F, {axisCentrality, axisIRBinning});
    }

    if (doAssocStudy && doprocessMonteCarlo) {
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
    if (doprocessMonteCarlo) {
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

    // For background decomposition
    if (fillBkgQAhistos && doprocessMonteCarlo) {
      histos.add("BkgStudy/h2dPtVsMassSigma_All", "h2dPtVsMassSigma_All", kTH2D, {axisPt, axisSigmaMass});
      histos.add("BkgStudy/h2dPtVsMassSigma_TrueDaughters", "h2dPtVsMassSigma_TrueDaughters", kTH2D, {axisPt, axisSigmaMass});
      histos.add("BkgStudy/h2dPtVsMassSigma_TrueGammaFakeLambda", "h2dPtVsMassSigma_TrueGammaFakeLambda", kTH2D, {axisPt, axisSigmaMass});
      histos.add("BkgStudy/h2dPtVsMassSigma_FakeGammaTrueLambda", "h2dPtVsMassSigma_FakeGammaTrueLambda", kTH2D, {axisPt, axisSigmaMass});
      histos.add("BkgStudy/h2dPtVsMassSigma_FakeDaughters", "h2dPtVsMassSigma_FakeDaughters", kTH2D, {axisPt, axisSigmaMass});
      histos.add("BkgStudy/h2dTrueDaughtersMatrix", "h2dTrueDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("BkgStudy/h2dTrueGammaFakeLambdaMatrix", "h2dTrueGammaFakeLambdaMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("BkgStudy/h2dFakeGammaTrueLambdaMatrix", "h2dFakeGammaTrueLambdaMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("BkgStudy/h2dFakeDaughtersMatrix", "h2dFakeDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
    }

    // For Pi0 QA
    if (doPi0QA) {
      histos.add("Pi0QA/h3dMassPi0BeforeSel_MCAssoc", "h3dMassPi0BeforeSel_MCAssoc", kTH3D, {axisCentrality, axisPt, axisPi0Mass});
      histos.add("Pi0QA/h3dMassPi0AfterSel_MCAssoc", "h3dMassPi0AfterSel_MCAssoc", kTH3D, {axisCentrality, axisPt, axisPi0Mass});
      histos.add("Pi0QA/h3dMassPi0BeforeSel_Candidates", "h3dMassPi0BeforeSel_Candidates", kTH3D, {axisCentrality, axisPt, axisPi0Mass});
      histos.add("Pi0QA/h3dMassPi0AfterSel_Candidates", "h3dMassPi0AfterSel_Candidates", kTH3D, {axisCentrality, axisPt, axisPi0Mass});
    }

    if (doprocessGeneratedRun3) {

      histos.add("Gen/hGenEvents", "hGenEvents", kTH2F, {{axisNch}, {2, -0.5f, +1.5f}});
      histos.get<TH2>(HIST("Gen/hGenEvents"))->GetYaxis()->SetBinLabel(1, "All gen. events");
      histos.get<TH2>(HIST("Gen/hGenEvents"))->GetYaxis()->SetBinLabel(2, "Gen. with at least 1 rec. events");

      histos.add("Gen/hGenEventCentrality", "hGenEventCentrality", kTH1F, {{101, 0.0f, 101.0f}});
      histos.add("Gen/hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", kTH2F, {axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("Gen/hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", kTH2F, {axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("Gen/hCentralityVsMultMC", "hCentralityVsMultMC", kTH2F, {{101, 0.0f, 101.0f}, axisNch});
      histos.add("Gen/h2dGenGamma", "h2dGenGamma", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/h2dGenLambda", "h2dGenLambda", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/h2dGenAntiLambda", "h2dGenAntiLambda", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/h2dGenGammaVsMultMC_RecoedEvt", "h2dGenGammaVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenLambdaVsMultMC_RecoedEvt", "h2dGenLambdaVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenAntiLambdaVsMultMC_RecoedEvt", "h2dGenAntiLambdaVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenGammaVsMultMC", "h2dGenGammaVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenLambdaVsMultMC", "h2dGenLambdaVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenAntiLambdaVsMultMC", "h2dGenAntiLambdaVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("Gen/hEventPVzMC", "hEventPVzMC", kTH1F, {{100, -20.0f, +20.0f}});
      histos.add("Gen/hCentralityVsPVzMC", "hCentralityVsPVzMC", kTH2F, {{101, 0.0f, 101.0f}, {100, -20.0f, +20.0f}});

      auto hPrimaryV0s = histos.add<TH1>("Gen/hPrimaryV0s", "hPrimaryV0s", kTH1D, {{2, -0.5f, 1.5f}});
      hPrimaryV0s->GetXaxis()->SetBinLabel(1, "All V0s");
      hPrimaryV0s->GetXaxis()->SetBinLabel(2, "Primary V0s");
    }
  }

  template <typename TCollision>
  bool IsEventAccepted(TCollision const& collision, bool fillHists)
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
    // Fetch interaction rate only if required (in order to limit ccdb calls)
    double interactionRate = (eventSelections.minIR >= 0 || eventSelections.maxIR >= 0) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource, fIRCrashOnNull) * 1.e-3 : -1;
    if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 18 /* Below min IR */);

    if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 19 /* Above max IR */);

    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    histos.fill(HIST("hEventCentrality"), centrality);
    return true;
  }

  void runBkgAnalysis(bool fIsSigma, bool fIsAntiSigma, int PhotonPDGCode, int PhotonPDGCodeMother, int LambdaPDGCode, int LambdaPDGCodeMother, float sigmapT, float sigmaMass)
  {
    histos.fill(HIST("BkgStudy/h2dPtVsMassSigma_All"), sigmapT, sigmaMass);

    // Real Gamma x Real Lambda - but not from the same sigma0/antisigma0!
    if ((PhotonPDGCode == 22) && ((LambdaPDGCode == 3122) || (LambdaPDGCode == -3122)) && (!fIsSigma && !fIsAntiSigma)) {
      histos.fill(HIST("BkgStudy/h2dPtVsMassSigma_TrueDaughters"), sigmapT, sigmaMass);
      histos.fill(HIST("BkgStudy/h2dTrueDaughtersMatrix"), LambdaPDGCodeMother, PhotonPDGCodeMother);
    }

    // Real Gamma x fake Lambda
    if ((PhotonPDGCode == 22) && (LambdaPDGCode != 3122) && (LambdaPDGCode != -3122)) {
      histos.fill(HIST("BkgStudy/h2dPtVsMassSigma_TrueGammaFakeLambda"), sigmapT, sigmaMass);
      histos.fill(HIST("BkgStudy/h2dTrueGammaFakeLambdaMatrix"), LambdaPDGCodeMother, PhotonPDGCodeMother);
    }

    // Fake Gamma x Real Lambda
    if ((PhotonPDGCode != 22) && ((LambdaPDGCode == 3122) || (LambdaPDGCode == -3122))) {
      histos.fill(HIST("BkgStudy/h2dPtVsMassSigma_FakeGammaTrueLambda"), sigmapT, sigmaMass);
      histos.fill(HIST("BkgStudy/h2dFakeGammaTrueLambdaMatrix"), LambdaPDGCodeMother, PhotonPDGCodeMother);
    }

    // Fake Gamma x Fake Lambda
    if ((PhotonPDGCode != 22) && (LambdaPDGCode != 3122) && (LambdaPDGCode != -3122)) {
      histos.fill(HIST("BkgStudy/h2dPtVsMassSigma_FakeDaughters"), sigmapT, sigmaMass);
      histos.fill(HIST("BkgStudy/h2dFakeDaughtersMatrix"), LambdaPDGCodeMother, PhotonPDGCodeMother);
    }
  }

  template <typename TCollision, typename TV0Object>
  void analyzeV0CollAssoc(TCollision const& collision, TV0Object const& fullv0s, std::vector<int> selV0Indices, float IR, bool isPhotonAnalysis)
  {
    auto v0MCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();

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
  // Simulated processing
  // Return the list of indices to the recoed collision associated to a given MC collision.
  template <typename TMCollisions, typename TCollisions>
  std::vector<int> getListOfRecoCollIndices(TMCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      for (auto const& collision : groupedCollisions) {
        // consider event selections in the recoed <-> gen collision association, for the denominator (or numerator) of the efficiency (or signal loss)?
        if (eventSelections.useEvtSelInDenomEff) {
          if (!IsEventAccepted(collision, false)) {
            continue;
          }
        }
        // Find the collision with the biggest nbr of PV contributors
        // Follows what was done here: https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/mcCollsExtra.cxx#L93
        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          bestCollisionIndex = collision.globalIndex();
        }
      }
      listBestCollisionIdx[mcCollision.globalIndex()] = bestCollisionIndex;
    }
    return listBestCollisionIdx;
  }

  // ______________________________________________________
  // Simulated processing
  // Fill generated event information (for event loss/splitting estimation)
  template <typename TMCCollisions, typename TCollisions>
  void fillGeneratedEventProperties(TMCCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      // Apply selections on MC collisions
      if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
        continue;
      }
      if (doPPAnalysis) { // we are in pp
        if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
          continue;
        }

        if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
          continue;
        }
      }

      histos.fill(HIST("Gen/hGenEvents"), mcCollision.multMCNParticlesEta05(), 0 /* all gen. events*/);

      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      // Check if there is at least one of the reconstructed collisions associated to this MC collision
      // If so, we consider it
      bool atLeastOne = false;
      int biggestNContribs = -1;
      float centrality = 100.5f;
      int nCollisions = 0;
      for (auto const& collision : groupedCollisions) {

        if (!IsEventAccepted(collision, false)) {
          continue;
        }

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
        }

        nCollisions++;
        atLeastOne = true;
      }

      histos.fill(HIST("Gen/hCentralityVsNcoll_beforeEvSel"), centrality, groupedCollisions.size());
      histos.fill(HIST("Gen/hCentralityVsNcoll_afterEvSel"), centrality, nCollisions);
      histos.fill(HIST("Gen/hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta05());
      histos.fill(HIST("Gen/hCentralityVsPVzMC"), centrality, mcCollision.posZ());
      histos.fill(HIST("Gen/hEventPVzMC"), mcCollision.posZ());

      if (atLeastOne) {
        histos.fill(HIST("Gen/hGenEvents"), mcCollision.multMCNParticlesEta05(), 1 /* at least 1 rec. event*/);
        histos.fill(HIST("Gen/hGenEventCentrality"), centrality);
      }
    }
    return;
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  template <typename TMCCollisions, typename TV0MCs, typename TCollisions>
  void analyzeGeneratedV0s(TMCCollisions const& mcCollisions, TV0MCs const& V0MCCores, TCollisions const& collisions)
  {
    fillGeneratedEventProperties(mcCollisions, collisions);
    std::vector<int> listBestCollisionIdx = getListOfRecoCollIndices(mcCollisions, collisions);
    for (auto const& v0MC : V0MCCores) {
      if (!v0MC.has_straMCCollision())
        continue;

      histos.fill(HIST("Gen/hPrimaryV0s"), 0);
      if (!v0MC.isPhysicalPrimary())
        continue;

      histos.fill(HIST("Gen/hPrimaryV0s"), 1);

      // TODO: get generated sigma0s

      float ptmc = v0MC.ptMC();
      float ymc = 1e3;
      if (v0MC.pdgCode() == 22)
        ymc = RecoDecay::y(std::array{v0MC.pxMC(), v0MC.pyMC(), v0MC.pzMC()}, o2::constants::physics::MassGamma);

      else if (std::abs(v0MC.pdgCode()) == 3122)
        ymc = v0MC.rapidityMC(1);

      if (std::abs(ymc) > V0Rapidity)
        continue;

      auto mcCollision = v0MC.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
      if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
        continue;
      }
      if (doPPAnalysis) { // we are in pp
        if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
          continue;
        }

        if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
          continue;
        }
      }

      float centrality = 100.5f;
      if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
        auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
        centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

        if (v0MC.pdgCode() == 22) {
          histos.fill(HIST("Gen/h2dGenGammaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (v0MC.pdgCode() == 3122) {
          histos.fill(HIST("Gen/h2dGenLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (v0MC.pdgCode() == -3122) {
          histos.fill(HIST("Gen/h2dGenAntiLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
      }
      if (v0MC.pdgCode() == 22) {
        histos.fill(HIST("Gen/h2dGenGamma"), centrality, ptmc);
        histos.fill(HIST("Gen/h2dGenGammaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (v0MC.pdgCode() == 3122) {
        histos.fill(HIST("Gen/h2dGenLambda"), centrality, ptmc);
        histos.fill(HIST("Gen/h2dGenLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (v0MC.pdgCode() == -3122) {
        histos.fill(HIST("Gen/h2dGenAntiLambda"), centrality, ptmc);
        histos.fill(HIST("Gen/h2dGenAntiLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    }
  }

  template <typename TV0Object, typename TCollision>
  void runPi0QA(TV0Object const& gamma1, TV0Object const& gamma2, TCollision collision)
  {
    // Check if both V0s are made of the same tracks
    if (gamma1.posTrackExtraId() == gamma2.posTrackExtraId() ||
        gamma1.negTrackExtraId() == gamma2.negTrackExtraId()) {
      return;
    }

    // Calculate pi0 properties
    std::array<float, 3> pVecGamma1{gamma1.px(), gamma1.py(), gamma1.pz()};
    std::array<float, 3> pVecGamma2{gamma2.px(), gamma2.py(), gamma2.pz()};
    std::array arrpi0{pVecGamma1, pVecGamma2};
    float pi0Mass = RecoDecay::m(arrpi0, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassPhoton});
    float pi0Pt = RecoDecay::pt(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py()});
    float pi0Y = RecoDecay::y(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py(), gamma1.pz() + gamma2.pz()}, o2::constants::physics::MassPi0);
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

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
          histos.fill(HIST("Pi0QA/h3dMassPi0BeforeSel_MCAssoc"), centrality, pi0Pt, pi0Mass);
        }
      }
    }

    histos.fill(HIST("Pi0QA/h3dMassPi0BeforeSel_Candidates"), centrality, pi0Pt, pi0Mass);

    // Photon-specific selections
    auto posTrackGamma1 = gamma1.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma1 = gamma1.template negTrackExtra_as<dauTracks>();
    auto posTrackGamma2 = gamma2.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma2 = gamma2.template negTrackExtra_as<dauTracks>();

    // Gamma1 Selection
    bool passedTPCGamma1 = (TMath::Abs(posTrackGamma1.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas) ||
                           (TMath::Abs(negTrackGamma1.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas);

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
    bool passedTPCGamma2 = (TMath::Abs(posTrackGamma2.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas) ||
                           (TMath::Abs(negTrackGamma2.tpcNSigmaEl()) < Pi0PhotonMaxTPCNSigmas);

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
    histos.fill(HIST("Pi0QA/h3dMassPi0AfterSel_Candidates"), centrality, pi0Pt, pi0Mass);
    if (fIsMC && fIsPi0)
      histos.fill(HIST("Pi0QA/h3dMassPi0AfterSel_MCAssoc"), centrality, pi0Pt, pi0Mass);
  }

  // Process photon candidate
  template <typename TV0Object, typename TCollision>
  bool processPhotonCandidate(TV0Object const& gamma, TCollision collision)
  {
    if (gamma.v0Type() == 0)
      return false;

    if (useMLScores) {
      // Gamma selection:
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
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    histos.fill(HIST("PhotonSel/h3dPhotonMass"), centrality, gamma.pt(), gamma.mGamma());
    return true;
  }

  // Process photon candidate
  template <typename TV0Object, typename TCollision>
  bool processLambdaCandidate(TV0Object const& lambda, TCollision collision)
  {
    if (lambda.v0Type() != 1)
      return false;

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

    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    histos.fill(HIST("LambdaSel/h3dLambdaMass"), centrality, lambda.pt(), lambda.mLambda());
    histos.fill(HIST("LambdaSel/h3dALambdaMass"), centrality, lambda.pt(), lambda.mAntiLambda());

    return true;
  }
  ///////////
  // Process sigma candidate and store properties in object
  template <typename TV0Object, typename TCollision>
  bool buildSigma0(TV0Object const& lambda, TV0Object const& gamma, TCollision collision)
  {
    // Checking if both V0s are made of the very same tracks
    if (gamma.posTrackExtraId() == lambda.posTrackExtraId() ||
        gamma.negTrackExtraId() == lambda.negTrackExtraId()) {
      return false;
    }

    // Sigma0 candidate properties
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
    auto arrMom = std::array{pVecPhotons, pVecLambda};
    float sigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    float sigmaY = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);
    float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    // Before any selection
    histos.fill(HIST("SigmaSel/h3dMassSigma0BeforeSel"), centrality, SigmapT, sigmaMass);

    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 1.);
    histos.fill(HIST("SigmaSel/hSigmaMass"), sigmaMass);
    histos.fill(HIST("SigmaSel/hSigmaMassWindow"), sigmaMass - 1.192642);

    if (fillQAhistos) {
      histos.fill(HIST("GeneralQA/h2dMassGammaVsK0S"), gamma.mGamma(), gamma.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassLambdaVsK0S"), lambda.mLambda(), lambda.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassGammaVsLambda"), gamma.mGamma(), gamma.mLambda());
      histos.fill(HIST("GeneralQA/h2dMassLambdaVsGamma"), lambda.mLambda(), lambda.mGamma());
      histos.fill(HIST("GeneralQA/h3dMassSigma0VsDaupTs"), gamma.pt(), lambda.pt(), sigmaMass);
    }

    if (TMath::Abs(sigmaMass - 1.192642) > Sigma0Window)
      return false;

    histos.fill(HIST("SigmaSel/hSigmaY"), sigmaY);
    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 2.);

    if (TMath::Abs(sigmaY) > SigmaMaxRap)
      return false;

    histos.fill(HIST("SigmaSel/hSigmaMassSelected"), sigmaMass);
    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 3.);

    if (fillQAhistos) {
      histos.fill(HIST("GeneralQA/h2dMassGammaVsK0SAfterMassSel"), gamma.mGamma(), gamma.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassLambdaVsK0SAfterMassSel"), lambda.mLambda(), lambda.mK0Short());
      histos.fill(HIST("GeneralQA/h2dMassGammaVsLambdaAfterMassSel"), gamma.mGamma(), lambda.mLambda());
      histos.fill(HIST("GeneralQA/h2dV0XY"), gamma.x(), gamma.y());
    }

    histos.fill(HIST("SigmaSel/h3dMassSigma0AfterSel"), centrality, SigmapT, sigmaMass);

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
    float fSigmaCentrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();
    uint64_t fSigmaTimeStamp = coll.timestamp();
    int fSigmaRunNumber = coll.runNumber();

    // Filling TTree for ML analysis
    sigma0cores(fSigmapT, fSigmaMass, fSigmaRap, fSigmaOPAngle, fSigmaCentrality, fSigmaRunNumber, fSigmaTimeStamp);

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
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0DerivedMCDatas const& fullV0s, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    // Initialize auxiliary vectors
    std::vector<int> bestGammasArray;
    std::vector<int> bestLambdasArray;

    // brute force grouped index construction
    std::vector<std::vector<int>> v0grouped(collisions.size());

    for (const auto& v0 : fullV0s) {
      v0grouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }

    for (const auto& coll : collisions) {
      // Clear vectors
      bestGammasArray.clear();
      bestLambdasArray.clear();

      if (!IsEventAccepted(coll, true))
        continue;

      float centrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();

      bool fhasMCColl = false;
      if (coll.has_straMCCollision())
        fhasMCColl = true;

      //_______________________________________________
      // Retrieving IR info
      float interactionRate = -1;
      if (fGetIR) {
        interactionRate = rateFetcher.fetch(ccdb.service, coll.timestamp(), coll.runNumber(), irSource, fIRCrashOnNull) * 1.e-3;
        if (interactionRate < 0)
          histos.get<TH1>(HIST("GeneralQA/hRunNumberNegativeIR"))->Fill(Form("%d", coll.runNumber()), 1);

        histos.fill(HIST("GeneralQA/hInteractionRate"), interactionRate);
        histos.fill(HIST("GeneralQA/hCentralityVsInteractionRate"), centrality, interactionRate);
      }

      //_______________________________________________
      // V0s loop
      for (size_t i = 0; i < v0grouped[coll.globalIndex()].size(); i++) {
        auto v0 = fullV0s.rawIteratorAt(v0grouped[coll.globalIndex()][i]);

        if (!v0.has_v0MCCore())
          continue;

        auto v0MC = v0.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

        if (v0MC.pdgCode() == 22) {
          histos.fill(HIST("MC/h2dGammaXYConversion"), v0.x(), v0.y());
          float GammaY = TMath::Abs(RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma));
          if (GammaY < 0.5) {                                                                                                                // rapidity selection
            histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocGamma"), centrality, v0.pt());                                            // isgamma
          }
        }

        float lambdaY = TMath::Abs(RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassLambda));
        if (lambdaY < 0.5) {
          if (v0MC.pdgCode() == 3122) // Is Lambda
            histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocLambda"), centrality, v0.pt());
          if (v0MC.pdgCode() == -3122) // Is AntiLambda
            histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocALambda"), centrality, v0.pt());
        }

        if (processPhotonCandidate(v0, coll))          // selecting photons
          bestGammasArray.push_back(v0.globalIndex()); // Save indices of best gamma candidates

        if (processLambdaCandidate(v0, coll))           // selecting lambdas
          bestLambdasArray.push_back(v0.globalIndex()); // Save indices of best lambda candidates
      }

      //_______________________________________________
      // Pi0 optional loop
      if (doPi0QA) {
        for (size_t i = 0; i < bestGammasArray.size(); ++i) {
          auto gamma1 = fullV0s.rawIteratorAt(bestGammasArray[i]);
          for (size_t j = i + 1; j < bestGammasArray.size(); ++j) {
            auto gamma2 = fullV0s.rawIteratorAt(bestGammasArray[j]);
            runPi0QA(gamma1, gamma2, coll);
          }
        }
      }

      //_______________________________________________
      // Wrongly collision association study
      if (doAssocStudy && fhasMCColl) {
        analyzeV0CollAssoc(coll, fullV0s, bestGammasArray, interactionRate, true);   // Gamma
        analyzeV0CollAssoc(coll, fullV0s, bestLambdasArray, interactionRate, false); // Lambda
      }

      //_______________________________________________
      // Sigma0 loop
      for (size_t i = 0; i < bestGammasArray.size(); ++i) {
        auto gamma = fullV0s.rawIteratorAt(bestGammasArray[i]);

        if (!gamma.has_v0MCCore())
          continue;

        auto gammaMC = gamma.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

        bool fIsPhotonCorrectlyAssign = false;
        if (fhasMCColl) {
          auto gammaMCCollision = coll.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
          fIsPhotonCorrectlyAssign = (gammaMC.straMCCollisionId() == gammaMCCollision.globalIndex());
        }

        for (size_t j = 0; j < bestLambdasArray.size(); ++j) {
          auto lambda = fullV0s.rawIteratorAt(bestLambdasArray[j]);

          if (!lambda.has_v0MCCore())
            continue;

          auto lambdaMC = lambda.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

          // Sigma0 candidate properties
          std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
          std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};
          auto arrMom = std::array{pVecPhotons, pVecLambda};
          float SigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
          float SigmapT = RecoDecay::pt(array{gamma.px() + lambda.px(), gamma.py() + lambda.py()});
          float SigmaY = TMath::Abs(RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0));

          // MC properties
          bool fIsSigma = false;
          bool fIsAntiSigma = false;
          bool fIsPhotonPrimary = gammaMC.isPhysicalPrimary();
          bool fIsLambdaPrimary = lambdaMC.isPhysicalPrimary();
          bool fIsLambdaCorrectlyAssign = false;

          int PhotonCandPDGCode = gammaMC.pdgCode();
          int PhotonCandPDGCodeMother = gammaMC.pdgCodeMother();
          int LambdaCandPDGCode = lambdaMC.pdgCode();
          int LambdaCandPDGCodeMother = lambdaMC.pdgCodeMother();

          float SigmaMCpT = RecoDecay::pt(array{gammaMC.pxMC() + lambdaMC.pxMC(), gammaMC.pyMC() + lambdaMC.pyMC()});
          float PhotonMCpT = RecoDecay::pt(array{gammaMC.pxMC(), gammaMC.pyMC()});
          float LambdaMCpT = RecoDecay::pt(array{lambdaMC.pxMC(), lambdaMC.pyMC()});

          if (fhasMCColl) {
            auto lambdaMCCollision = coll.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
            fIsLambdaCorrectlyAssign = (lambdaMC.straMCCollisionId() == lambdaMCCollision.globalIndex());
          }

          if ((PhotonCandPDGCode == 22) && (PhotonCandPDGCodeMother == 3212) && (LambdaCandPDGCode == 3122) && (LambdaCandPDGCodeMother == 3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()))
            fIsSigma = true;
          if ((PhotonCandPDGCode == 22) && (PhotonCandPDGCodeMother == -3212) && (LambdaCandPDGCode == -3122) && (LambdaCandPDGCodeMother == -3212) && (gamma.motherMCPartId() == lambda.motherMCPartId()))
            fIsAntiSigma = true;

          if (SigmaY < 0.5) {
            if (fIsSigma) {
              histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocSigma0"), centrality, SigmaMCpT);
              histos.fill(HIST("MC/h2dSigma0PtVsLambdaPtBeforeSel_MCAssoc"), SigmaMCpT, LambdaMCpT);
              histos.fill(HIST("MC/h2dSigma0PtVsGammaPtBeforeSel_MCAssoc"), SigmaMCpT, PhotonMCpT);
            }
            if (fIsAntiSigma)
              histos.fill(HIST("MC/h2dPtVsCentralityBeforeSel_MCAssocASigma0"), centrality, SigmaMCpT);
          }

          // Build sigma0 candidate, please
          if (!buildSigma0(lambda, gamma, coll))
            continue;

          if (SigmaY < 0.5) {
            if (fIsSigma)
              histos.fill(HIST("MC/h2dPtVsCentralityAfterSel_MCAssocSigma0"), centrality, SigmaMCpT);
            if (fIsAntiSigma)
              histos.fill(HIST("MC/h2dPtVsCentralityAfterSel_MCAssocASigma0"), centrality, SigmaMCpT);
          }

          if (fillBkgQAhistos)
            runBkgAnalysis(fIsSigma, fIsAntiSigma, PhotonCandPDGCode, PhotonCandPDGCodeMother, LambdaCandPDGCode, LambdaCandPDGCodeMother, SigmapT, SigmaMass);

          // Fill Tables please
          sigma0mccores(fIsSigma, fIsAntiSigma, SigmaMCpT,
                        PhotonCandPDGCode, PhotonCandPDGCodeMother, fIsPhotonPrimary, PhotonMCpT, fIsPhotonCorrectlyAssign,
                        LambdaCandPDGCode, LambdaCandPDGCodeMother, fIsLambdaPrimary, LambdaMCpT, fIsLambdaCorrectlyAssign);

          // Filling tables with accepted candidates
          fillTables(lambda, gamma, coll);

          nSigmaCandidates++;
          if (nSigmaCandidates % 10000 == 0)
            LOG(info) << "Sigma0 Candidates built: " << nSigmaCandidates;
        }
      }
    }
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0StandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    // Initialize auxiliary vectors
    std::vector<int> bestGammasArray;
    std::vector<int> bestLambdasArray;

    // brute force grouped index construction
    std::vector<std::vector<int>> v0grouped(collisions.size());

    for (const auto& v0 : fullV0s) {
      v0grouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }

    for (const auto& coll : collisions) {
      // Clear vectors
      bestGammasArray.clear();
      bestLambdasArray.clear();

      if (!IsEventAccepted(coll, true))
        continue;

      float centrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();

      //_______________________________________________
      // Retrieving IR info
      float interactionRate = -1;
      if (fGetIR) {
        interactionRate = rateFetcher.fetch(ccdb.service, coll.timestamp(), coll.runNumber(), irSource, fIRCrashOnNull) * 1.e-3;
        if (interactionRate < 0)
          histos.get<TH1>(HIST("GeneralQA/hRunNumberNegativeIR"))->Fill(Form("%d", coll.runNumber()), 1);

        histos.fill(HIST("GeneralQA/hInteractionRate"), interactionRate);
        histos.fill(HIST("GeneralQA/hCentralityVsInteractionRate"), centrality, interactionRate);
      }

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
      // Pi0 optional loop
      if (doPi0QA) {
        for (size_t i = 0; i < bestGammasArray.size(); ++i) {
          auto gamma1 = fullV0s.rawIteratorAt(bestGammasArray[i]);
          for (size_t j = i + 1; j < bestGammasArray.size(); ++j) {
            auto gamma2 = fullV0s.rawIteratorAt(bestGammasArray[j]);
            runPi0QA(gamma1, gamma2, coll);
          }
        }
      }

      //_______________________________________________
      // Sigma0 nested loop
      for (size_t i = 0; i < bestGammasArray.size(); ++i) {
        auto gamma = fullV0s.rawIteratorAt(bestGammasArray[i]);

        for (size_t j = 0; j < bestLambdasArray.size(); ++j) {
          auto lambda = fullV0s.rawIteratorAt(bestLambdasArray[j]);

          // Building sigma0 candidate
          if (!buildSigma0(lambda, gamma, coll))
            continue;

          // Filling tables with accepted candidates
          fillTables(lambda, gamma, coll);

          nSigmaCandidates++;
          if (nSigmaCandidates % 10000 == 0)
            LOG(info) << "Sigma0 Candidates built: " << nSigmaCandidates;
        }
      }
    }
  }

  // Simulated processing in Run 3 (subscribes to MC information too)
  void processGeneratedRun3(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const& V0MCCores, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions)
  {
    analyzeGeneratedV0s(mcCollisions, V0MCCores, collisions);
  }

  PROCESS_SWITCH(sigma0builder, processMonteCarlo, "process as if MC data", false);
  PROCESS_SWITCH(sigma0builder, processRealData, "process as if real data", true);
  PROCESS_SWITCH(sigma0builder, processGeneratedRun3, "process generated MC collisions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigma0builder>(cfgc)};
}
