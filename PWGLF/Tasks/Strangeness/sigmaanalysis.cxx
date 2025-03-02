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
  // Event selection
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  Configurable<float> minIR{"minIR", -1, "Min Interaction Rate (kHz). Leave -1 if no selection desired."};
  Configurable<float> maxIR{"maxIR", -1, "Max Interaction Rate (kHz). Leave -1 if no selection desired."};

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

  // Analysis strategy:
  Configurable<bool> fUseMLSel{"fUseMLSel", false, "Flag to use ML selection. If False, the standard selection is applied."};
  Configurable<bool> fProcessMonteCarlo{"fProcessMonteCarlo", false, "Flag to process MC data."};
  Configurable<bool> fselLambdaTPCPID{"fselLambdaTPCPID", true, "Flag to select lambda-like candidates using TPC NSigma."};
  Configurable<bool> fselLambdaTOFPID{"fselLambdaTOFPID", false, "Flag to select lambda-like candidates using TOF NSigma."};
  Configurable<bool> doMCAssociation{"doMCAssociation", false, "Flag to process only signal candidates. Use only with processMonteCarlo!"};
  Configurable<bool> doPhotonLambdaSelQA{"doPhotonLambdaSelQA", false, "Flag to fill photon and lambda QA histos!"};

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
  Configurable<float> PhotonPhiMin1{"PhotonPhiMin1", -1, "Phi min value for photons, region 1 (leave negative if no selection desired)"};
  Configurable<float> PhotonPhiMax1{"PhotonPhiMax1", -1, "Phi max value for photons, region 1 (leave negative if no selection desired)"};
  Configurable<float> PhotonPhiMin2{"PhotonPhiMin2", -1, "Phi max value for photons, region 2 (leave negative if no selection desired)"};
  Configurable<float> PhotonPhiMax2{"PhotonPhiMax2", -1, "Phi min value for photons, region 2 (leave negative if no selection desired)"};

  Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.5, "Max sigma0 rapidity"};

  // Axis
  // base properties
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisInvPt{"axisInvPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis axisDeltaPt{"axisDeltaPt", {400, -50.0, 50.0}, ""};
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
  ConfigurableAxis axisCandSel{"axisCandSel", {32, 0.5f, +32.5f}, "Candidate Selection"};

  // ML
  ConfigurableAxis MLProb{"MLOutput", {100, 0.0f, 1.0f}, ""};
  int nSigmaCandidates = 0;
  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

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

    // All candidates received
    histos.add("GeneralQA/hInteractionRate", "hInteractionRate", kTH1F, {axisIRBinning});
    histos.add("GeneralQA/hInteractionRatePerColl", "hInteractionRatePerColl", kTH1F, {axisIRBinning});
    histos.add("GeneralQA/hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2F, {axisCentrality, axisIRBinning});
    histos.add("GeneralQA/hCentralityVsInteractionRatePerColl", "hCentralityVsInteractionRatePerColl", kTH2F, {axisCentrality, axisIRBinning});
    histos.add("GeneralQA/h2dArmenterosBeforeSel", "h2dArmenterosBeforeSel", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/h2dArmenterosAfterSel", "h2dArmenterosAfterSel", {HistType::kTH2F, {axisAPAlpha, axisAPQt}});
    histos.add("GeneralQA/hMassSigma0BeforeSel", "hMassSigma0BeforeSel", kTH1F, {axisSigmaMass});

    // Candidates Counters
    histos.add("GeneralQA/hCandidateAnalysisSelection", "hCandidateAnalysisSelection", kTH1F, {axisCandSel});
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(2, "Photon V0Type");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(3, "Photon DauPt");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(4, "Photon DCAToPV");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(5, "Photon DCADau");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(6, "Photon TPCCrossedRows");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(7, "Photon TPCNSigmaEl");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(8, "Photon TPCNSigmaPi");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(9, "Photon Pt");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(10, "Photon Y/Eta");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(11, "Photon Radius");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(12, "Photon RZ line");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(13, "Photon QT");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(14, "Photon Alpha");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(15, "Photon CosPA");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(16, "Photon PsiPair");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(17, "Photon Phi");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(18, "Photon Mass");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(19, "Lambda Radius");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(20, "Lambda DCADau");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(21, "Lambda QT");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(22, "Lambda Alpha");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(23, "Lambda CosPA");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(24, "Lambda Y/Eta");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(25, "Lambda TPCCrossedRows");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(26, "Lambda ITSNCls");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(27, "Lambda Lifetime");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(28, "Lambda/ALambda PID");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(29, "Lambda/ALambda DCAToPV");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(30, "Lambda/ALambda Mass");
    histos.get<TH1>(HIST("GeneralQA/hCandidateAnalysisSelection"))->GetXaxis()->SetBinLabel(31, "Sigma Y");

    // Photon Selection QA histos
    histos.add("GeneralQA/hPhotonV0Type", "hPhotonV0Type", kTH1F, {{8, 0.5f, 8.5f}});
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
    histos.add("GeneralQA/hPhotonPhi", "hPhotonPhi", kTH1F, {{200, 0, 2 * o2::constants::math::PI}});
    histos.add("GeneralQA/h3dPhotonMass", "h3dPhotonMass", kTH3F, {axisCentrality, axisPt, axisPhotonMass});

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
    histos.add("GeneralQA/h3dLambdaMass", "h3dLambdaMass", kTH3F, {axisCentrality, axisPt, axisLambdaMass});
    histos.add("GeneralQA/h2dTPCvsTOFNSigma_ALambdaPr", "h2dTPCvsTOFNSigma_ALambdaPr", {HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});
    histos.add("GeneralQA/h2dTPCvsTOFNSigma_ALambdaPi", "h2dTPCvsTOFNSigma_ALambdaPi", {HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});
    histos.add("GeneralQA/hALambdaDCANegToPV", "hALambdaDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/hALambdaDCAPosToPV", "hALambdaDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("GeneralQA/h3dAntiLambdaMass", "h3dAntiLambdaMass", kTH3F, {axisCentrality, axisPt, axisLambdaMass});

    histos.add("GeneralQA/hPhotonMassSelected", "hPhotonMassSelected", kTH1F, {axisPhotonMass});
    histos.add("GeneralQA/hLambdaMassSelected", "hLambdaMassSelected", kTH1F, {axisLambdaMass});
    histos.add("GeneralQA/hAntiLambdaMassSelected", "hAntiLambdaMassSelected", kTH1F, {axisLambdaMass});

    histos.add("GeneralQA/hSigmaY", "hSigmaY", kTH1F, {axisRapidity});
    histos.add("GeneralQA/hSigmaOPAngle", "hSigmaOPAngle", kTH1F, {{140, 0.0f, +7.0f}});

    // Specific sigma0 QA
    histos.add("SigmaSelQA/h2dPhotonV0Type", "h2dPhotonV0Type", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonMass", "h2dPhotonMass", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonDaupT", "h2dPhotonDaupT", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonDCADauToPV", "h2dPhotonDCADauToPV", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonDCADau", "h2dPhotonDCADau", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonDauTPCCR", "h2dPhotonDauTPCCR", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonTPCNSigmaEl", "h2dPhotonTPCNSigmaEl", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonpT", "h2dPhotonpT", {HistType::kTH2F, {axisPt, axisSigmaMass}}); //
    histos.add("SigmaSelQA/h2dPhotonY", "h2dPhotonY", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonRadius", "h2dPhotonRadius", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dRZCut", "h2dRZCut", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonArmenteros", "h2dPhotonArmenteros", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonCosPA", "h2dPhotonCosPA", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dPhotonPsiPair", "h2dPhotonPsiPair", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaRadius", "h2dLambdaRadius", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaDCADau", "h2dLambdaDCADau", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaArmenteros", "h2dLambdaArmenteros", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaCosPA", "h2dLambdaCosPA", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaY", "h2dLambdaY", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaDauTPCCR", "h2dLambdaDauTPCCR", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaDauITSCls", "h2dLambdaDauITSCls", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaLifeTime", "h2dLambdaLifeTime", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dTPCvsTOFNSigma_Lambda", "h2dTPCvsTOFNSigma_Lambda", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaDCADauToPV", "h2dLambdaDCADauToPV", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dLambdaMass", "h2dLambdaMass", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dTPCvsTOFNSigma_ALambda", "h2dTPCvsTOFNSigma_ALambda", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dALambdaDCADauToPV", "h2dALambdaDCADauToPV", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dAntiLambdaMass", "h2dAntiLambdaMass", {HistType::kTH2F, {axisPt, axisSigmaMass}});
    histos.add("SigmaSelQA/h2dSigmaY", "h2dSigmaY", {HistType::kTH2F, {axisPt, axisSigmaMass}});

    // Specific photon QA
    histos.add("PhotonSelQA/h2dPhotonBaseline", "h2dPhotonBaseline", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonV0Type", "h2dPhotonV0Type", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonMass", "h2dPhotonMass", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonDaupT", "h2dPhotonDaupT", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonDCADauToPV", "h2dPhotonDCADauToPV", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonDCADau", "h2dPhotonDCADau", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonDauTPCCR", "h2dPhotonDauTPCCR", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonTPCNSigmaEl", "h2dPhotonTPCNSigmaEl", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonpT", "h2dPhotonpT", {HistType::kTH2F, {axisPt, axisPhotonMass}}); //
    histos.add("PhotonSelQA/h2dPhotonY", "h2dPhotonY", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonRadius", "h2dPhotonRadius", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dRZCut", "h2dRZCut", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonArmenteros", "h2dPhotonArmenteros", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonCosPA", "h2dPhotonCosPA", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonPsiPair", "h2dPhotonPsiPair", {HistType::kTH2F, {axisPt, axisPhotonMass}});
    histos.add("PhotonSelQA/h2dPhotonPhi", "h2dPhotonPhi", {HistType::kTH2F, {axisPt, axisPhotonMass}});

    // Specific Lambda/ALambda QA
    histos.add("LambdaSelQA/h2dLambdaBaseline", "h2dLambdaBaseline", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaRadius", "h2dLambdaRadius", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaDCADau", "h2dLambdaDCADau", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaArmenteros", "h2dLambdaArmenteros", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaCosPA", "h2dLambdaCosPA", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaY", "h2dLambdaY", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaDauTPCCR", "h2dLambdaDauTPCCR", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaDauITSCls", "h2dLambdaDauITSCls", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaLifeTime", "h2dLambdaLifeTime", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dTPCvsTOFNSigma_Lambda", "h2dTPCvsTOFNSigma_Lambda", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaDCADauToPV", "h2dLambdaDCADauToPV", {HistType::kTH2F, {axisPt, axisLambdaMass}});
    histos.add("LambdaSelQA/h2dLambdaMass", "h2dLambdaMass", {HistType::kTH2F, {axisPt, axisLambdaMass}});

    // For Signal Extraction

    // Sigma0
    histos.add("Sigma0/h3dMassSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("Sigma0/h3dPhotonRadiusVsMassSigma0", "h3dPhotonRadiusVsMassSigma0", kTH3F, {axisCentrality, axisRadius, axisSigmaMass});
    histos.add("Sigma0/hMassSigma0", "hMassSigma0", kTH1F, {axisSigmaMass});
    histos.add("Sigma0/hPtSigma0", "hPtSigma0", kTH1F, {axisPt});
    histos.add("Sigma0/hRapiditySigma0", "hRapiditySigma0", kTH1F, {axisRapidity});

    // AntiSigma0
    histos.add("AntiSigma0/h3dMassAntiSigma0", "h3dMassAntiSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
    histos.add("AntiSigma0/h3dPhotonRadiusVsMassAntiSigma0", "h3dPhotonRadiusVsMassAntiSigma0", kTH3F, {axisCentrality, axisRadius, axisSigmaMass});
    histos.add("AntiSigma0/hMassAntiSigma0", "hMassAntiSigma0", kTH1F, {axisSigmaMass});
    histos.add("AntiSigma0/hPtAntiSigma0", "hPtAntiSigma0", kTH1F, {axisPt});
    histos.add("AntiSigma0/hRapidityAntiSigma0", "hRapidityAntiSigma0", kTH1F, {axisRapidity});

    if (fProcessMonteCarlo) {

      // Kinematic
      histos.add("MC/h3dMassSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
      histos.add("MC/h3dMassAntiSigma0", "h3dMassSigma0", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
      histos.add("MC/h3dMassAllSigma0sBeforesel", "h3dMassAllSigma0sBeforesel", kTH3F, {axisCentrality, axisPt, axisSigmaMass});
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

      // 1/pT Resolution:
      histos.add("MC/h2dLambdaPtResolution", "h2dLambdaPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
      histos.add("MC/h2dAntiLambdaPtResolution", "h2dAntiLambdaPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
      histos.add("MC/h2dGammaPtResolution", "h2dGammaPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
      histos.add("MC/h2dSigma0PtResolution", "h2dSigma0PtResolution", kTH2D, {axisInvPt, axisDeltaPt});
      histos.add("MC/h2dAntiSigma0PtResolution", "h2dAntiSigma0PtResolution", kTH2D, {axisInvPt, axisDeltaPt});
      histos.add("MC/h3dLambdaPtResoVsTPCCR", "h3dLambdaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, {320, -160.0f, 160.0f}});
      histos.add("MC/h3dAntiLambdaPtResoVsTPCCR", "h3dAntiLambdaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, {320, -160.0f, 160.0f}});
      histos.add("MC/h3dGammaPtResoVsTPCCR", "h3dGammaPtResoVsTPCCR", kTH3F, {axisInvPt, axisDeltaPt, {320, -160.0f, 160.0f}});

      // pTMC info:
      histos.add("MC/h2dSigmaMCPtVsLambdaMCPt", "h2dSigmaMCPtVsLambdaMCPt", kTH2D, {axisPt, axisPt});
      histos.add("MC/h2dSigmaMCPtVsGammaMCPt", "h2dSigmaMCPtVsGammaMCPt", kTH2D, {axisPt, axisPt});
    }
  }

  template <typename TCollision>
  bool IsEventAccepted(TCollision collision)
  // check whether the collision passes our collision selections
  {
    histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);
    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
    if (doPPAnalysis) { // we are in pp
      if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
        return false;
      }
      histos.fill(HIST("hEventSelection"), 16 /* INEL > 0 */);
      if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
        return false;
      }
      histos.fill(HIST("hEventSelection"), 17 /* INEL > 1 */);
    } else { // we are in Pb-Pb
      float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
        return false;
      }
      histos.fill(HIST("hEventSelection"), 16 /* Below min occupancy */);
      if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
        return false;
      }
      histos.fill(HIST("hEventSelection"), 17 /* Above max occupancy */);
    }
    return true;
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

      bool isMCTrueLambda = false;
      bool isMCTruePhoton = false;
      if constexpr (requires { cand.lambdaCandPDGCode(); cand.photonCandPDGCode(); }) {
        if (cand.photonCandPDGCode() == 22)
          isMCTruePhoton = true;
        if (cand.lambdaCandPDGCode() == 3122)
          isMCTrueLambda = true;
      }

      // Photon Selections
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 1.);
      histos.fill(HIST("GeneralQA/hPhotonV0Type"), cand.photonV0Type());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonBaseline"), cand.photonPt(), cand.photonMass());
      if (cand.photonV0Type() != Photonv0TypeSel && Photonv0TypeSel > -1)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonV0Type"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonV0Type"), cand.photonPt(), cand.photonMass());

      // histos.fill(HIST("GeneralQA/h3dPhotonMass"), cand.sigmaCentrality(), cand.photonPt(), cand.photonMass());
      // histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 2.);
      // if (TMath::Abs(cand.photonMass()) > PhotonMaxMass)
      //   return false;
      // histos.fill(HIST("SigmaSelQA/h2dPhotonMass"), cand.sigmapT(), cand.sigmaMass());
      // if (isMCTruePhoton || doPhotonLambdaSelQA) histos.fill(HIST("PhotonSelQA/h2dPhotonMass"), cand.photonPt(), cand.photonMass());

      histos.fill(HIST("GeneralQA/hPhotonNegpT"), cand.photonNegPt());
      histos.fill(HIST("GeneralQA/hPhotonPospT"), cand.photonPosPt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 2.);
      if ((cand.photonPosPt() < PhotonDauMinPt) || (cand.photonNegPt() < PhotonDauMinPt))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonDaupT"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonDaupT"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hPhotonDCANegToPV"), cand.photonDCANegPV());
      histos.fill(HIST("GeneralQA/hPhotonDCAPosToPV"), cand.photonDCAPosPV());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 3.);
      if ((TMath::Abs(cand.photonDCAPosPV()) < PhotonMinDCADauToPv) || (TMath::Abs(cand.photonDCANegPV()) < PhotonMinDCADauToPv))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonDCADauToPV"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonDCADauToPV"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 4.);
      histos.fill(HIST("GeneralQA/hPhotonDCADau"), cand.photonDCADau());
      if (TMath::Abs(cand.photonDCADau()) > PhotonMaxDCAV0Dau)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonDCADau"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonDCADau"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hPhotonPosTPCCR"), cand.photonPosTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hPhotonNegTPCCR"), cand.photonNegTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 5.);
      if ((cand.photonPosTPCCrossedRows() < PhotonMinTPCCrossedRows) || (cand.photonNegTPCCrossedRows() < PhotonMinTPCCrossedRows))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonDauTPCCR"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonDauTPCCR"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/h2dPhotonPosTPCNSigmaEl"), cand.photonPosPt(), cand.photonPosTPCNSigmaEl());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 6.);
      if (((cand.photonPosTPCNSigmaEl() < PhotonMinTPCNSigmas) || (cand.photonPosTPCNSigmaEl() > PhotonMaxTPCNSigmas)))
        return false;
      histos.fill(HIST("GeneralQA/h2dPhotonNegTPCNSigmaEl"), cand.photonNegPt(), cand.photonNegTPCNSigmaEl());
      if (((cand.photonNegTPCNSigmaEl() < PhotonMinTPCNSigmas) || (cand.photonNegTPCNSigmaEl() > PhotonMaxTPCNSigmas)))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonTPCNSigmaEl"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonTPCNSigmaEl"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/h2dPhotonPosTPCNSigmaPi"), cand.photonPosPt(), cand.photonPosTPCNSigmaPi());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 7.);
      if (((TMath::Abs(cand.photonPosTPCNSigmaPi()) < PiMaxTPCNSigmas) && cand.photonPosPt() <= piMaxpT))
        return false;
      histos.fill(HIST("GeneralQA/h2dPhotonNegTPCNSigmaPi"), cand.photonNegPt(), cand.photonNegTPCNSigmaPi());
      if (((TMath::Abs(cand.photonNegTPCNSigmaPi()) < PiMaxTPCNSigmas) && cand.photonNegPt() <= piMaxpT))
        return false;
      histos.fill(HIST("GeneralQA/hPhotonpT"), cand.photonPt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 8.);
      if ((cand.photonPt() < PhotonMinPt) || (cand.photonPt() > PhotonMaxPt))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonpT"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonpT"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hPhotonY"), cand.photonY());
      histos.fill(HIST("GeneralQA/hPhotonPosEta"), cand.photonPosEta());
      histos.fill(HIST("GeneralQA/hPhotonNegEta"), cand.photonNegEta());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 9.);
      if ((TMath::Abs(cand.photonY()) > PhotonMaxRap) || (TMath::Abs(cand.photonPosEta()) > PhotonMaxDauEta) || (TMath::Abs(cand.photonNegEta()) > PhotonMaxDauEta))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonY"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonY"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hPhotonRadius"), cand.photonRadius());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 10.);
      if ((cand.photonRadius() < PhotonMinRadius) || (cand.photonRadius() > PhotonMaxRadius))
        return false;
      float photonRZLineCut = TMath::Abs(cand.photonZconv()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-PhotonMaxDauEta))) - PhotonLineCutZ0;
      histos.fill(HIST("SigmaSelQA/h2dPhotonRadius"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonRadius"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hPhotonZ"), cand.photonZconv());
      histos.fill(HIST("GeneralQA/h2dRZCut"), cand.photonRadius(), photonRZLineCut);
      histos.fill(HIST("GeneralQA/h2dRZPlane"), cand.photonZconv(), cand.photonRadius());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 11.);
      if ((TMath::Abs(cand.photonRadius()) < photonRZLineCut) || (TMath::Abs(cand.photonZconv()) > PhotonMaxZ))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dRZCut"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dRZCut"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/h2dPhotonArmenteros"), cand.photonAlpha(), cand.photonQt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 12.);
      if (cand.photonQt() > PhotonMaxQt)
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 13.);
      if (TMath::Abs(cand.photonAlpha()) > PhotonMaxAlpha)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonArmenteros"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonArmenteros"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 14.);
      histos.fill(HIST("GeneralQA/hPhotonCosPA"), cand.photonCosPA());
      if (cand.photonCosPA() < PhotonMinV0cospa)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonCosPA"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonCosPA"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 15.);
      histos.fill(HIST("GeneralQA/hPhotonPsiPair"), cand.photonPsiPair());
      if (TMath::Abs(cand.photonPsiPair()) > PhotonPsiPairMax)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonPsiPair"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonPsiPair"), cand.photonPt(), cand.photonMass());

      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 16.);
      histos.fill(HIST("GeneralQA/hPhotonPhi"), cand.photonPhi());
      if (((cand.photonPhi() < PhotonPhiMin1) || ((cand.photonPhi() > PhotonPhiMax1) && (cand.photonPhi() < PhotonPhiMin2)) || ((cand.photonPhi() > PhotonPhiMax2) && (PhotonPhiMax2 != -1))))
        return false;
      if ((isMCTruePhoton || doPhotonLambdaSelQA) && (TMath::Abs(cand.photonMass()) <= PhotonMaxMass))
        histos.fill(HIST("PhotonSelQA/h2dPhotonPhi"), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/h3dPhotonMass"), cand.sigmaCentrality(), cand.photonPt(), cand.photonMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 17.);
      if (TMath::Abs(cand.photonMass()) > PhotonMaxMass)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dPhotonMass"), cand.sigmapT(), cand.sigmaMass());
      if (isMCTruePhoton || doPhotonLambdaSelQA)
        histos.fill(HIST("PhotonSelQA/h2dPhotonMass"), cand.photonPt(), cand.photonMass());

      // ####
      //  Lambda selections
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 18.);
      histos.fill(HIST("GeneralQA/hLambdaRadius"), cand.lambdaRadius());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaBaseline"), cand.lambdaPt(), cand.lambdaMass());
      if ((cand.lambdaRadius() < LambdaMinv0radius) || (cand.lambdaRadius() > LambdaMaxv0radius))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaRadius"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaRadius"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hLambdaDCADau"), cand.lambdaDCADau());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 19.);
      if (TMath::Abs(cand.lambdaDCADau()) > LambdaMaxDCAV0Dau)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaDCADau"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaDCADau"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/h2dLambdaArmenteros"), cand.lambdaAlpha(), cand.lambdaQt());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 20.);
      if ((cand.lambdaQt() < LambdaMinQt) || (cand.lambdaQt() > LambdaMaxQt))
        return false;
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 21.);
      if ((TMath::Abs(cand.lambdaAlpha()) < LambdaMinAlpha) || (TMath::Abs(cand.lambdaAlpha()) > LambdaMaxAlpha))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaArmenteros"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaArmenteros"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hLambdaCosPA"), cand.lambdaCosPA());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 22.);
      if (cand.lambdaCosPA() < LambdaMinv0cospa)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaCosPA"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaCosPA"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hLambdaY"), cand.lambdaY());
      histos.fill(HIST("GeneralQA/hLambdaPosEta"), cand.lambdaPosEta());
      histos.fill(HIST("GeneralQA/hLambdaNegEta"), cand.lambdaNegEta());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 23.);
      if ((TMath::Abs(cand.lambdaY()) > LambdaMaxRap) || (TMath::Abs(cand.lambdaPosEta()) > LambdaMaxDauEta) || (TMath::Abs(cand.lambdaNegEta()) > LambdaMaxDauEta))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaY"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaY"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hLambdaPosTPCCR"), cand.lambdaPosTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hLambdaNegTPCCR"), cand.lambdaNegTPCCrossedRows());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 24.);
      if ((cand.lambdaPosTPCCrossedRows() < LambdaMinTPCCrossedRows) || (cand.lambdaNegTPCCrossedRows() < LambdaMinTPCCrossedRows))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaDauTPCCR"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaDauTPCCR"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hLambdaPosITSCls"), cand.lambdaPosITSCls());
      histos.fill(HIST("GeneralQA/hLambdaNegITSCls"), cand.lambdaNegITSCls());
      histos.fill(HIST("GeneralQA/hLambdaPosChi2PerNc"), cand.lambdaPosChi2PerNcl());
      histos.fill(HIST("GeneralQA/hLambdaNegChi2PerNc"), cand.lambdaNegChi2PerNcl());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 25.);
      // check minimum number of ITS clusters + reject ITS afterburner tracks if requested
      bool posIsFromAfterburner = cand.lambdaPosChi2PerNcl() < 0;
      bool negIsFromAfterburner = cand.lambdaNegChi2PerNcl() < 0;
      if (cand.lambdaPosITSCls() < LambdaMinITSclusters && (!LambdaRejectPosITSafterburner || posIsFromAfterburner))
        return false;
      if (cand.lambdaNegITSCls() < LambdaMinITSclusters && (!LambdaRejectNegITSafterburner || negIsFromAfterburner))
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaDauITSCls"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaDauITSCls"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hLambdaLifeTime"), cand.lambdaLifeTime());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 26.);
      if (cand.lambdaLifeTime() > LambdaMaxLifeTime)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dLambdaLifeTime"), cand.sigmapT(), cand.sigmaMass());
      if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
        histos.fill(HIST("LambdaSelQA/h2dLambdaLifeTime"), cand.lambdaPt(), cand.lambdaMass());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 27.);
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

        histos.fill(HIST("SigmaSelQA/h2dTPCvsTOFNSigma_Lambda"), cand.sigmapT(), cand.sigmaMass());
        if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
          histos.fill(HIST("LambdaSelQA/h2dTPCvsTOFNSigma_Lambda"), cand.lambdaPt(), cand.lambdaMass());
        // DCA Selection
        histos.fill(HIST("GeneralQA/hLambdaDCANegToPV"), cand.lambdaDCANegPV());
        histos.fill(HIST("GeneralQA/hLambdaDCAPosToPV"), cand.lambdaDCAPosPV());
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 28.);
        if ((TMath::Abs(cand.lambdaDCAPosPV()) < LambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < LambdaMinDCANegToPv))
          return false;
        histos.fill(HIST("SigmaSelQA/h2dLambdaDCADauToPV"), cand.sigmapT(), cand.sigmaMass());
        if ((isMCTrueLambda || doPhotonLambdaSelQA) && (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) < LambdaWindow))
          histos.fill(HIST("LambdaSelQA/h2dLambdaDCADauToPV"), cand.lambdaPt(), cand.lambdaMass());

        // Mass Selection
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 29.);
        histos.fill(HIST("GeneralQA/h3dLambdaMass"), cand.sigmaCentrality(), cand.lambdaPt(), cand.lambdaMass());
        if (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) > LambdaWindow)
          return false;
        histos.fill(HIST("SigmaSelQA/h2dLambdaMass"), cand.sigmapT(), cand.sigmaMass());
        if (isMCTrueLambda || doPhotonLambdaSelQA)
          histos.fill(HIST("LambdaSelQA/h2dLambdaMass"), cand.lambdaPt(), cand.lambdaMass());

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

        histos.fill(HIST("SigmaSelQA/h2dTPCvsTOFNSigma_ALambda"), cand.sigmapT(), cand.sigmaMass());
        // DCA Selection
        histos.fill(HIST("GeneralQA/hALambdaDCANegToPV"), cand.lambdaDCANegPV());
        histos.fill(HIST("GeneralQA/hALambdaDCAPosToPV"), cand.lambdaDCAPosPV());
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 28.);
        if ((TMath::Abs(cand.lambdaDCAPosPV()) < ALambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < ALambdaMinDCANegToPv))
          return false;

        histos.fill(HIST("SigmaSelQA/h2dALambdaDCADauToPV"), cand.sigmapT(), cand.sigmaMass());
        // Mass Selection
        histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 29.);
        histos.fill(HIST("GeneralQA/h3dAntiLambdaMass"), cand.sigmaCentrality(), cand.lambdaPt(), cand.antilambdaMass());
        if (TMath::Abs(cand.antilambdaMass() - o2::constants::physics::MassLambda0) > LambdaWindow)
          return false;
        histos.fill(HIST("SigmaSelQA/h2dAntiLambdaMass"), cand.sigmapT(), cand.sigmaMass());
      }

      // Sigma0 selection
      histos.fill(HIST("GeneralQA/hSigmaY"), cand.sigmaRapidity());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 30.);
      if (TMath::Abs(cand.sigmaRapidity()) > SigmaMaxRap)
        return false;
      histos.fill(HIST("SigmaSelQA/h2dSigmaY"), cand.sigmapT(), cand.sigmaMass());
      histos.fill(HIST("GeneralQA/hSigmaOPAngle"), cand.sigmaOPAngle());
      histos.fill(HIST("GeneralQA/hCandidateAnalysisSelection"), 31.);
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
      histos.fill(HIST("MC/h2dSigmaMCPtVsLambdaMCPt"), sigma.sigmaMCPt(), sigma.lambdaMCPt());
      histos.fill(HIST("MC/h2dSigmaMCPtVsGammaMCPt"), sigma.sigmaMCPt(), sigma.photonMCPt());
      histos.fill(HIST("MC/h3dMassAllSigma0sBeforesel"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
      histos.fill(HIST("MC/h2dArmenterosBeforeSel"), sigma.photonAlpha(), sigma.photonQt());
      histos.fill(HIST("MC/h2dArmenterosBeforeSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
      histos.fill(HIST("MC/hMassSigma0BeforeSel"), sigma.sigmaMass());
      histos.fill(HIST("MC/hPtSigma0BeforeSel"), sigma.sigmapT());
      histos.fill(HIST("MC/hPtGammaCand_BeforeSel"), sigma.photonPt());
      histos.fill(HIST("MC/hPtLambdaCand_BeforeSel"), sigma.lambdaPt());
      histos.fill(HIST("MC/hPtSigmaCand_BeforeSel"), sigma.sigmapT());

      if (sigma.photonCandPDGCode() == 22) {
        histos.fill(HIST("MC/hPtTrueGamma_BeforeSel"), sigma.photonPt());
        if (sigma.photonMCPt() > 0) {
          histos.fill(HIST("MC/h3dGammaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), -1 * sigma.photonNegTPCCrossedRows()); // 1/pT resolution
          histos.fill(HIST("MC/h3dGammaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), sigma.photonPosTPCCrossedRows());      // 1/pT resolution
          histos.fill(HIST("MC/h2dGammaPtResolution"), 1.f / sigma.photonMCPt(), 1.f / sigma.photonPt() - 1.f / sigma.photonMCPt());                                        // pT resolution
        }
      }
      if (sigma.lambdaCandPDGCode() == 3122) {
        histos.fill(HIST("MC/hPtTrueLambda_BeforeSel"), sigma.lambdaPt());
        if (sigma.lambdaMCPt() > 0) {
          histos.fill(HIST("MC/h2dLambdaPtResolution"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt());                                        // 1/pT resolution
          histos.fill(HIST("MC/h3dLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), -1 * sigma.lambdaNegTPCCrossedRows()); // 1/pT resolution
          histos.fill(HIST("MC/h3dLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), sigma.lambdaPosTPCCrossedRows());      // 1/pT resolution
        }
      }
      if (sigma.lambdaCandPDGCode() == -3122) {
        if (sigma.lambdaMCPt() > 0) {
          histos.fill(HIST("MC/h2dAntiLambdaPtResolution"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt());                                        // pT resolution
          histos.fill(HIST("MC/h3dAntiLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), -1 * sigma.lambdaNegTPCCrossedRows()); // 1/pT resolution
          histos.fill(HIST("MC/h3dAntiLambdaPtResoVsTPCCR"), 1.f / sigma.lambdaMCPt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdaMCPt(), sigma.lambdaPosTPCCrossedRows());      // 1/pT resolution
        }
        histos.fill(HIST("MC/hPtTrueAntiLambda_BeforeSel"), sigma.lambdaPt());
      }
      if (sigma.isSigma()) {
        histos.fill(HIST("MC/hPtTrueSigma_BeforeSel"), sigma.sigmapT());
        if (sigma.sigmaMCPt() > 0)
          histos.fill(HIST("MC/h2dSigma0PtResolution"), 1.f / sigma.sigmaMCPt(), 1.f / sigma.sigmapT() - 1.f / sigma.sigmaMCPt()); // pT resolution
      }
      if (sigma.isAntiSigma()) {
        histos.fill(HIST("MC/hPtTrueAntiSigma_BeforeSel"), sigma.sigmapT());
        if (sigma.sigmaMCPt() > 0)
          histos.fill(HIST("MC/h2dAntiSigma0PtResolution"), 1.f / sigma.sigmaMCPt(), 1.f / sigma.sigmapT() - 1.f / sigma.sigmaMCPt()); // pT resolution
      }

      if (sigma.lambdaAlpha() > 0) { // Lambda Analysis
        if (!processSigmaCandidate(sigma, true))
          continue;

        // For Lambda PID Studies
        histos.fill(HIST("MC/hPtLambdaCand_AfterSel"), sigma.lambdaPt());
        histos.fill(HIST("MC/h3dTPCvsTOFNSigma_LambdaPr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma(), sigma.lambdaPt());
        histos.fill(HIST("MC/h3dTPCvsTOFNSigma_LambdaPi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma(), sigma.lambdaPt());

        if (sigma.lambdaCandPDGCode() == 3122) {
          histos.fill(HIST("MC/hPtTrueLambda_AfterSel"), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueLambdaPr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma(), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueLambdaPi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma(), sigma.lambdaPt());
        }
        histos.fill(HIST("GeneralQA/hLambdaMassSelected"), sigma.lambdaMass());
        histos.fill(HIST("Sigma0/hMassSigma0"), sigma.sigmaMass());
        histos.fill(HIST("Sigma0/hPtSigma0"), sigma.sigmapT());
        histos.fill(HIST("Sigma0/hRapiditySigma0"), sigma.sigmaRapidity());
        histos.fill(HIST("Sigma0/h3dMassSigma0"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
        histos.fill(HIST("Sigma0/h3dPhotonRadiusVsMassSigma0"), sigma.sigmaCentrality(), sigma.photonRadius(), sigma.sigmaMass());

        if (sigma.isSigma()) { // Signal study
          histos.fill(HIST("MC/h2dArmenterosAfterSel"), sigma.photonAlpha(), sigma.photonQt());
          histos.fill(HIST("MC/h2dArmenterosAfterSel"), sigma.lambdaAlpha(), sigma.lambdaQt());
          histos.fill(HIST("MC/hMassSigma0"), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtSigma0"), sigma.sigmapT());
          histos.fill(HIST("MC/h3dMassSigma0"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
          histos.fill(HIST("MC/h2dPtVsMassSigma_SignalOnly"), sigma.sigmapT(), sigma.sigmaMass());
          histos.fill(HIST("MC/hPtTrueSigma_AfterSel"), sigma.sigmapT());
        }
      } else { // AntiLambda Analysis
        if (!processSigmaCandidate(sigma, false))
          continue;

        if (sigma.lambdaCandPDGCode() == -3122) {
          histos.fill(HIST("MC/hPtTrueAntiLambda_AfterSel"), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueALambdaPr"), sigma.lambdaNegPrTPCNSigma(), sigma.aLambdaPrTOFNSigma(), sigma.lambdaPt());
          histos.fill(HIST("MC/h3dTPCvsTOFNSigma_TrueALambdaPi"), sigma.lambdaPosPiTPCNSigma(), sigma.aLambdaPiTOFNSigma(), sigma.lambdaPt());
        }
        histos.fill(HIST("GeneralQA/hAntiLambdaMassSelected"), sigma.antilambdaMass());
        histos.fill(HIST("AntiSigma0/hMassAntiSigma0"), sigma.sigmaMass());
        histos.fill(HIST("AntiSigma0/hPtAntiSigma0"), sigma.sigmapT());
        histos.fill(HIST("AntiSigma0/hRapidityAntiSigma0"), sigma.sigmaRapidity());
        histos.fill(HIST("AntiSigma0/h3dMassAntiSigma0"), sigma.sigmaCentrality(), sigma.sigmapT(), sigma.sigmaMass());
        histos.fill(HIST("AntiSigma0/h3dPhotonRadiusVsMassAntiSigma0"), sigma.sigmaCentrality(), sigma.photonRadius(), sigma.sigmaMass());
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
        histos.fill(HIST("Sigma0/h3dPhotonRadiusVsMassSigma0"), sigma.sigmaCentrality(), sigma.photonRadius(), sigma.sigmaMass());

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
        histos.fill(HIST("AntiSigma0/h3dPhotonRadiusVsMassAntiSigma0"), sigma.sigmaCentrality(), sigma.photonRadius(), sigma.sigmaMass());
      }
      histos.fill(HIST("GeneralQA/hPhotonMassSelected"), sigma.photonMass());
    }
  }

  void processCollisions(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions)
  {
    for (const auto& coll : collisions) {
      if (!IsEventAccepted(coll)) {
        continue;
      }
      double interactionRate = rateFetcher.fetch(ccdb.service, coll.timestamp(), coll.runNumber(), irSource) * 1.e-3;
      histos.fill(HIST("GeneralQA/hInteractionRatePerColl"), interactionRate);
      histos.fill(HIST("GeneralQA/hCentralityVsInteractionRatePerColl"), coll.centFT0C(), interactionRate);
    }
  }

  // PROCESS_SWITCH(sigmaanalysis, processCounterQA, "Check standard counter correctness", true);
  PROCESS_SWITCH(sigmaanalysis, processMonteCarlo, "Do Monte-Carlo-based analysis", false);
  PROCESS_SWITCH(sigmaanalysis, processRealData, "Do real data analysis", true);
  PROCESS_SWITCH(sigmaanalysis, processCollisions, "Process collisions to retrieve IR info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigmaanalysis>(cfgc)};
}
