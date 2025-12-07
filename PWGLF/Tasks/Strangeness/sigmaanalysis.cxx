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
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

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
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using std::array;
using MCSigma0s = soa::Join<aod::Sigma0Cores, aod::Sigma0PhotonExtras, aod::Sigma0LambdaExtras, aod::Sigma0MCCores, aod::SigmaCollRef>;
using Sigma0s = soa::Join<aod::Sigma0Cores, aod::Sigma0PhotonExtras, aod::Sigma0LambdaExtras, aod::SigmaCollRef>;

static const std::vector<std::string> PhotonSels = {"NoSel", "V0Type", "DCADauToPV",
                                                    "DCADau", "DauTPCCR", "TPCNSigmaEl", "V0pT",
                                                    "Y", "V0Radius", "RZCut", "Armenteros", "CosPA",
                                                    "PsiPair", "Phi", "Mass"};

static const std::vector<std::string> LambdaSels = {"NoSel", "V0Radius", "DCADau", "Armenteros",
                                                    "CosPA", "Y", "TPCCR", "DauITSCls", "Lifetime",
                                                    "TPCTOFPID", "DCADauToPV", "Mass"};

static const std::vector<std::string> DirList = {"BeforeSel", "AfterSel"};

struct sigmaanalysis {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  //__________________________________________________
  // For manual sliceBy
  // SliceCache cache;
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Event level
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<bool> fIRCrashOnNull{"fIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
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

  // Generated Sigma0s
  struct : ConfigurableGroup {
    std::string prefix = "genSelections"; // JSON group name
    Configurable<bool> mc_keepOnlyFromGenerator{"mc_keepOnlyFromGenerator", true, "if true, consider only particles from generator to calculate efficiency."};
    Configurable<float> mc_rapidityMin{"mc_rapidityMin", -0.5, "Min generated particle rapidity"};
    Configurable<float> mc_rapidityMax{"mc_rapidityMax", 0.5, "Max generated particle rapidity"};
  } genSelections;

  // QA
  Configurable<bool> fdoSigma0QA{"doSigma0QA", false, "if true, perform Sigma0 QA analysis. Only works with MC."};
  Configurable<bool> fillBkgQAhistos{"fillBkgQAhistos", false, "if true, fill MC QA histograms for Bkg study. Only works with MC."};
  Configurable<bool> fillResoQAhistos{"fillResoQAhistos", false, "if true, fill MC QA histograms for pT resolution study. Only works with MC."};
  Configurable<bool> fillSelhistos{"fillSelhistos", true, "if true, fill QA histos for selections."};

  // Analysis strategy:
  Configurable<bool> doMCAssociation{"doMCAssociation", false, "Flag to process only signal candidates. Use only with processMonteCarlo!"};
  Configurable<bool> selRecoFromGenerator{"selRecoFromGenerator", false, "Flag to process only signal candidates from generator"};

  // For Selection:
  //// Lambda criteria::
  struct : ConfigurableGroup {
    std::string prefix = "lambdaSelections"; // JSON group name
    Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
    Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};
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
    Configurable<float> LambdaMinRapidity{"LambdaMinRapidity", -0.5, "v0 min rapidity"};
    Configurable<float> LambdaMaxRapidity{"LambdaMaxRapidity", 0.5, "v0 max rapidity"};
    Configurable<float> LambdaMinDauEta{"LambdaMinDauEta", -0.8, "Min pseudorapidity of daughter tracks"};
    Configurable<float> LambdaMaxDauEta{"LambdaMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<bool> fselLambdaTPCPID{"fselLambdaTPCPID", true, "Flag to select lambda-like candidates using TPC NSigma."};
    Configurable<bool> fselLambdaTOFPID{"fselLambdaTOFPID", false, "Flag to select lambda-like candidates using TOF NSigma."};
    Configurable<float> LambdaMaxTPCNSigmas{"LambdaMaxTPCNSigmas", 1e+9, "Max TPC NSigmas for daughters"};
    Configurable<float> LambdaPrMaxTOFNSigmas{"LambdaPrMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
    Configurable<float> LambdaPiMaxTOFNSigmas{"LambdaPiMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
    Configurable<int> LambdaMinTPCCrossedRows{"LambdaMinTPCCrossedRows", 50, "Min daughter TPC Crossed Rows"};
    Configurable<int> LambdaMinITSclusters{"LambdaMinITSclusters", 1, "minimum ITS clusters"};
    Configurable<bool> LambdaRejectPosITSafterburner{"LambdaRejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
    Configurable<bool> LambdaRejectNegITSafterburner{"LambdaRejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};

  } lambdaSelections;

  //// Photon criteria:
  struct : ConfigurableGroup {
    std::string prefix = "photonSelections"; // JSON group name
    Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
    Configurable<int> Photonv0TypeSel{"Photonv0TypeSel", 7, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> PhotonMinDCADauToPv{"PhotonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
    Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
    Configurable<int> PhotonMinTPCCrossedRows{"PhotonMinTPCCrossedRows", 30, "Min daughter TPC Crossed Rows"};
    Configurable<float> PhotonMinTPCNSigmas{"PhotonMinTPCNSigmas", -7, "Min TPC NSigmas for daughters"};
    Configurable<float> PhotonMaxTPCNSigmas{"PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
    Configurable<float> PhotonMinPt{"PhotonMinPt", 0.0, "Min photon pT (GeV/c)"};
    Configurable<float> PhotonMaxPt{"PhotonMaxPt", 50.0, "Max photon pT (GeV/c)"};
    Configurable<float> PhotonMinRapidity{"PhotonMinRapidity", -0.5, "v0 min rapidity"};
    Configurable<float> PhotonMaxRapidity{"PhotonMaxRapidity", 0.5, "v0 max rapidity"};
    Configurable<float> PhotonMinDauEta{"PhotonMinDauEta", -0.8, "Min pseudorapidity of daughter tracks"};
    Configurable<float> PhotonMaxDauEta{"PhotonMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<float> PhotonMinRadius{"PhotonMinRadius", 3.0, "Min photon conversion radius (cm)"};
    Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 115, "Max photon conversion radius (cm)"};
    Configurable<float> PhotonMinZ{"PhotonMinZ", -240, "Min photon conversion point z value (cm)"};
    Configurable<float> PhotonMaxZ{"PhotonMaxZ", 240, "Max photon conversion point z value (cm)"};
    Configurable<float> PhotonMaxQt{"PhotonMaxQt", 0.05, "Max photon qt value (AP plot) (GeV/c)"};
    Configurable<float> PhotonMaxAlpha{"PhotonMaxAlpha", 0.95, "Max photon alpha absolute value (AP plot)"};
    Configurable<float> PhotonMinV0cospa{"PhotonMinV0cospa", 0.80, "Min V0 CosPA"};
    Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};
    Configurable<float> PhotonPsiPairMax{"PhotonPsiPairMax", 1e+9, "maximum psi angle of the track pair"};
    Configurable<float> PhotonLineCutZ0{"PhotonLineCutZ0", 7.0, "The offset for the linecute used in the Z vs R plot"};
    Configurable<float> PhotonPhiMin1{"PhotonPhiMin1", -1, "Phi min value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> PhotonPhiMax1{"PhotonPhiMax1", -1, "Phi max value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> PhotonPhiMin2{"PhotonPhiMin2", -1, "Phi max value to reject photons, region 2 (leave negative if no selection desired)"};
    Configurable<float> PhotonPhiMax2{"PhotonPhiMax2", -1, "Phi min value to reject photons, region 2 (leave negative if no selection desired)"};
  } photonSelections;

  struct : ConfigurableGroup {
    std::string prefix = "sigma0Selections"; // JSON group name
    Configurable<float> Sigma0MinRapidity{"Sigma0MinRapidity", -0.5, "sigma0 min rapidity"};
    Configurable<float> Sigma0MaxRapidity{"Sigma0MaxRapidity", 0.5, "sigma0 max rapidity"};
    Configurable<float> Sigma0MaxRadius{"Sigma0MaxRadius", 200, "Max sigma0 decay radius"};
    Configurable<float> Sigma0MaxDCADau{"Sigma0MaxDCADau", 50, "Max sigma0 DCA between daughters"};
    Configurable<float> Sigma0MaxOPAngle{"Sigma0MaxOPAngle", 7, "Max sigma0 OP Angle between daughters"};
  } sigma0Selections;

  struct : ConfigurableGroup {
    std::string prefix = "pi0Selections"; // JSON group name
    Configurable<float> Pi0MinRapidity{"Pi0MinRapidity", -0.5, "pi0 min rapidity"};
    Configurable<float> Pi0MaxRapidity{"Pi0MaxRapidity", 0.5, "pi0 max rapidity"};
    Configurable<float> Pi0MaxRadius{"Pi0MaxRadius", 200, "Max sigma0 decay radius"};
    Configurable<float> Pi0MaxDCADau{"Pi0MaxDCADau", 50, "Max sigma0 DCA between daughters"};
  } pi0Selections;

  // Axis
  // base properties
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};
  ConfigurableAxis axisInvPt{"axisInvPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis axisDeltaPt{"axisDeltaPt", {400, -50.0, 50.0}, ""};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  ConfigurableAxis axisIRBinning{"axisIRBinning", {150, 0, 1500}, "Binning for the interaction rate (kHz)"};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "N_{ch}"};
  ConfigurableAxis axisGeneratorIds{"axisGeneratorIds", {256, -0.5f, 255.5f}, "axis for generatorIds"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {500, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.05f, 1.151f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.5f}, "M_{#Gamma}"};
  ConfigurableAxis axisPi0Mass{"axisPi0Mass", {200, 0.08f, 0.18f}, "M_{#Pi^{0}}"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality, PID and other axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisNCls{"axisNCls", {8, -0.5, 7.5}, "NCls"};
  ConfigurableAxis axisChi2PerNcl{"axisChi2PerNcl", {80, -40, 40}, "Chi2 Per Ncl"};
  ConfigurableAxis axisTPCNSigma{"axisTPCNSigma", {120, -30, 30}, "TPC NSigma"};
  ConfigurableAxis axisTOFNSigma{"axisTOFNSigma", {120, -30, 30}, "TOF NSigma"};
  ConfigurableAxis axisLifetime{"axisLifetime", {100, 0, 100}, "Chi2 Per Ncl"};

  // topological variable QA axes
  ConfigurableAxis axisV0Radius{"axisV0Radius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisV0PairRadius{"axisV0PairRadius", {200, 0.0f, 20.0f}, "V0Pair radius (cm)"};
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisCosPA{"axisCosPA", {200, 0.5f, 1.0f}, "Cosine of pointing angle"};
  ConfigurableAxis axisPA{"axisPA", {100, 0.0f, 1}, "Pointing angle"};
  ConfigurableAxis axisPsiPair{"axisPsiPair", {250, -5.0f, 5.0f}, "Psipair for photons"};
  ConfigurableAxis axisPhi{"axisPhi", {200, 0, 2 * o2::constants::math::PI}, "Phi for photons"};
  ConfigurableAxis axisZ{"axisZ", {120, -120.0f, 120.0f}, "V0 Z position (cm)"};

  ConfigurableAxis axisCandSel{"axisCandSel", {20, 0.5f, +20.5f}, "Candidate Selection"};

  // ML
  ConfigurableAxis MLProb{"MLOutput", {100, 0.0f, 1.0f}, ""};

  int NSigma0Cand = 0;
  void init(InitContext const&)
  {
    LOGF(info, "Initializing now: cross-checking correctness...");
    if ((doprocessRealData + doprocessMonteCarlo + doprocessPi0RealData + doprocessPi0MonteCarlo > 1) ||
        (doprocessGeneratedRun3 + doprocessPi0GeneratedRun3 > 1)) {
      LOGF(fatal, "You have enabled more than one process function. Please check your configuration! Aborting now.");
    }

    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    // Event Counters
    histos.add("hEventCentrality", "hEventCentrality", kTH1D, {axisCentrality});

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

    if (fGetIR) {
      histos.add("GeneralQA/hRunNumberNegativeIR", "", kTH1D, {{1, 0., 1.}});
      histos.add("GeneralQA/hInteractionRate", "hInteractionRate", kTH1D, {axisIRBinning});
      histos.add("GeneralQA/hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2D, {axisCentrality, axisIRBinning});
    }

    if (doprocessRealData || doprocessMonteCarlo) {
      for (const auto& histodir : DirList) {

        if (fillSelhistos) {
          histos.add(histodir + "/Photon/hTrackCode", "hTrackCode", kTH1D, {{11, 0.5f, 11.5f}});
          histos.add(histodir + "/Photon/hV0Type", "hV0Type", kTH1D, {{8, 0.5f, 8.5f}});
          histos.add(histodir + "/Photon/hDCANegToPV", "hDCANegToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/Photon/hDCAPosToPV", "hDCAPosToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/Photon/hDCADau", "hDCADau", kTH1D, {axisDCAdau});
          histos.add(histodir + "/Photon/hPosTPCCR", "hPosTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/Photon/hNegTPCCR", "hNegTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/Photon/hPosTPCNSigmaEl", "hPosTPCNSigmaEl", kTH1D, {axisTPCNSigma});
          histos.add(histodir + "/Photon/hNegTPCNSigmaEl", "hNegTPCNSigmaEl", kTH1D, {axisTPCNSigma});

          histos.add(histodir + "/Photon/hpT", "hpT", kTH1D, {axisPt});
          histos.add(histodir + "/Photon/hY", "hY", kTH1D, {axisRapidity});
          histos.add(histodir + "/Photon/hPosEta", "hPosEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/Photon/hNegEta", "hNegEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/Photon/hRadius", "hRadius", kTH1D, {axisV0Radius});
          histos.add(histodir + "/Photon/hZ", "hZ", kTH1D, {axisZ});
          histos.add(histodir + "/Photon/h2dRZCut", "h2dRZCut", kTH2D, {axisZ, axisV0Radius});
          histos.add(histodir + "/Photon/h2dRZPlane", "h2dRZPlane", kTH2D, {axisZ, axisV0Radius});
          histos.add(histodir + "/Photon/hCosPA", "hCosPA", kTH1D, {axisCosPA});
          histos.add(histodir + "/Photon/hPsiPair", "hPsiPair", kTH1D, {axisPsiPair});
          histos.add(histodir + "/Photon/hPhi", "hPhi", kTH1D, {axisPhi});
          histos.add(histodir + "/Photon/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisPhotonMass});
          histos.add(histodir + "/Photon/hMass", "hMass", kTH1D, {axisPhotonMass});

          histos.add(histodir + "/Lambda/hTrackCode", "hTrackCode", kTH1D, {{11, 0.5f, 11.5f}});
          histos.add(histodir + "/Lambda/hRadius", "hRadius", kTH1D, {axisV0Radius});
          histos.add(histodir + "/Lambda/hDCADau", "hDCADau", kTH1D, {axisDCAdau});
          histos.add(histodir + "/Lambda/hCosPA", "hCosPA", kTH1D, {axisCosPA});
          histos.add(histodir + "/Lambda/hY", "hY", kTH1D, {axisRapidity});
          histos.add(histodir + "/Lambda/hPosEta", "hPosEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/Lambda/hNegEta", "hNegEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/Lambda/hPosTPCCR", "hPosTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/Lambda/hNegTPCCR", "hNegTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/Lambda/hPosITSCls", "hPosITSCls", kTH1D, {axisNCls});
          histos.add(histodir + "/Lambda/hNegITSCls", "hNegITSCls", kTH1D, {axisNCls});
          histos.add(histodir + "/Lambda/hPosChi2PerNc", "hPosChi2PerNc", kTH1D, {axisChi2PerNcl});
          histos.add(histodir + "/Lambda/hNegChi2PerNc", "hNegChi2PerNc", kTH1D, {axisChi2PerNcl});
          histos.add(histodir + "/Lambda/hLifeTime", "hLifeTime", kTH1D, {axisLifetime});
          histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_LambdaPr", "h2dTPCvsTOFNSigma_LambdaPr", kTH2D, {axisTPCNSigma, axisTOFNSigma});
          histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_LambdaPi", "h2dTPCvsTOFNSigma_LambdaPi", kTH2D, {axisTPCNSigma, axisTOFNSigma});
          histos.add(histodir + "/Lambda/hLambdaDCANegToPV", "hLambdaDCANegToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/Lambda/hLambdaDCAPosToPV", "hLambdaDCAPosToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/Lambda/hLambdapT", "hLambdapT", kTH1D, {axisPt});
          histos.add(histodir + "/Lambda/hLambdaMass", "hLambdaMass", kTH1D, {axisLambdaMass});
          histos.add(histodir + "/Lambda/h3dLambdaMass", "h3dLambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});
          histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_ALambdaPr", "h2dTPCvsTOFNSigma_ALambdaPr", kTH2D, {axisTPCNSigma, axisTOFNSigma});
          histos.add(histodir + "/Lambda/h2dTPCvsTOFNSigma_ALambdaPi", "h2dTPCvsTOFNSigma_ALambdaPi", kTH2D, {axisTPCNSigma, axisTOFNSigma});
          histos.add(histodir + "/Lambda/hALambdaDCANegToPV", "hALambdaDCANegToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/Lambda/hALambdaDCAPosToPV", "hALambdaDCAPosToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/Lambda/hALambdapT", "hALambdapT", kTH1D, {axisPt});
          histos.add(histodir + "/Lambda/hAntiLambdaMass", "hAntiLambdaMass", kTH1D, {axisLambdaMass});
          histos.add(histodir + "/Lambda/h3dAntiLambdaMass", "h3dAntiLambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});
        }

        histos.add(histodir + "/h2dArmenteros", "h2dArmenteros", kTH2D, {axisAPAlpha, axisAPQt});

        histos.add(histodir + "/Sigma0/hMass", "hMass", kTH1D, {axisSigmaMass});
        histos.add(histodir + "/Sigma0/hPt", "hPt", kTH1D, {axisPt});
        histos.add(histodir + "/Sigma0/hY", "hY", kTH1D, {axisRapidity});
        histos.add(histodir + "/Sigma0/hRadius", "hRadius", kTH1D, {axisV0PairRadius});
        histos.add(histodir + "/Sigma0/h2dRadiusVspT", "h2dRadiusVspT", kTH2D, {axisV0PairRadius, axisPt});
        histos.add(histodir + "/Sigma0/hDCAPairDau", "hDCAPairDau", kTH1D, {axisDCAdau});
        histos.add(histodir + "/Sigma0/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisSigmaMass});
        histos.add(histodir + "/Sigma0/h3dPhotonYMass", "h3dPhotonYMass", kTH3D, {axisRapidity, axisPt, axisSigmaMass});
        histos.add(histodir + "/Sigma0/h3dPhotonRadiusMass", "h3dPhotonRadiusMass", kTH3D, {axisV0Radius, axisPt, axisSigmaMass});
        histos.add(histodir + "/Sigma0/h3dOPAngleVsMass", "h3dOPAngleVsMass", kTH3D, {{140, 0.0f, +7.0f}, axisPt, axisSigmaMass});

        histos.add(histodir + "/ASigma0/hMass", "hMass", kTH1D, {axisSigmaMass});
        histos.add(histodir + "/ASigma0/hPt", "hPt", kTH1D, {axisPt});
        histos.add(histodir + "/ASigma0/hY", "hY", kTH1D, {axisRapidity});
        histos.add(histodir + "/ASigma0/hRadius", "hRadius", kTH1D, {axisV0PairRadius});
        histos.add(histodir + "/ASigma0/h2dRadiusVspT", "h2dRadiusVspT", kTH2D, {axisV0PairRadius, axisPt});
        histos.add(histodir + "/ASigma0/hDCAPairDau", "hDCAPairDau", kTH1D, {axisDCAdau});
        histos.add(histodir + "/ASigma0/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisSigmaMass});
        histos.add(histodir + "/ASigma0/h3dPhotonYMass", "h3dPhotonYMass", kTH3D, {axisRapidity, axisPt, axisSigmaMass});
        histos.add(histodir + "/ASigma0/h3dPhotonRadiusMass", "h3dPhotonRadiusMass", kTH3D, {axisV0Radius, axisPt, axisSigmaMass});
        histos.add(histodir + "/ASigma0/h3dOPAngleVsMass", "h3dOPAngleVsMass", kTH3D, {{140, 0.0f, +7.0f}, axisPt, axisSigmaMass});

        // Process MC
        if (doprocessMonteCarlo) {

          if (fillSelhistos) {
            histos.add(histodir + "/MC/Photon/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1D, {{2, 0.0f, 2.0f}});
            histos.add(histodir + "/MC/Photon/hPt", "hPt", kTH1D, {axisPt});
            histos.add(histodir + "/MC/Photon/hMCPt", "hMCPt", kTH1D, {axisPt});
            histos.add(histodir + "/MC/Photon/hPosTPCNSigmaEl", "hPosTPCNSigmaEl", kTH1D, {axisTPCNSigma});
            histos.add(histodir + "/MC/Photon/hNegTPCNSigmaEl", "hNegTPCNSigmaEl", kTH1D, {axisTPCNSigma});
            histos.add(histodir + "/MC/Photon/h2dPAVsPt", "h2dPAVsPt", kTH2D, {axisPA, axisPt});
            histos.add(histodir + "/MC/Photon/hPt_BadCollAssig", "hPt_BadCollAssig", kTH1D, {axisPt});
            histos.add(histodir + "/MC/Photon/h2dPAVsPt_BadCollAssig", "h2dPAVsPt_BadCollAssig", kTH2D, {axisPA, axisPt});

            histos.add(histodir + "/MC/Lambda/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1D, {{2, 0.0f, 2.0f}});
            histos.add(histodir + "/MC/Lambda/hPt", "hPt", kTH1D, {axisPt});
            histos.add(histodir + "/MC/Lambda/hMCPt", "hMCPt", kTH1D, {axisPt});
            histos.add(histodir + "/MC/Lambda/h3dTPCvsTOFNSigma_Pr", "h3dTPCvsTOFNSigma_Pr", kTH3D, {axisTPCNSigma, axisTOFNSigma, axisPt});
            histos.add(histodir + "/MC/Lambda/h3dTPCvsTOFNSigma_Pi", "h3dTPCvsTOFNSigma_Pi", kTH3D, {axisTPCNSigma, axisTOFNSigma, axisPt});

            histos.add(histodir + "/MC/ALambda/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1D, {{2, 0.0f, 2.0f}});
            histos.add(histodir + "/MC/ALambda/hPt", "hPt", kTH1D, {axisPt});
            histos.add(histodir + "/MC/ALambda/hMCPt", "hMCPt", kTH1D, {axisPt});
            histos.add(histodir + "/MC/ALambda/h3dTPCvsTOFNSigma_Pr", "h3dTPCvsTOFNSigma_Pr", kTH3D, {axisTPCNSigma, axisTOFNSigma, axisPt});
            histos.add(histodir + "/MC/ALambda/h3dTPCvsTOFNSigma_Pi", "h3dTPCvsTOFNSigma_Pi", kTH3D, {axisTPCNSigma, axisTOFNSigma, axisPt});
          }

          histos.add(histodir + "/MC/h2dArmenteros", "h2dArmenteros", kTH2D, {axisAPAlpha, axisAPQt});

          histos.add(histodir + "/MC/Sigma0/hPt", "hPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/Sigma0/hMCPt", "hMCPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/Sigma0/hMass", "hMass", kTH1D, {axisSigmaMass});
          histos.add(histodir + "/MC/Sigma0/hMCProcess", "hMCProcess", kTH1D, {{50, -0.5f, 49.5f}});
          histos.add(histodir + "/MC/Sigma0/hGenRadius", "hGenRadius", kTH1D, {axisV0PairRadius});
          histos.add(histodir + "/MC/Sigma0/h2dMCPtVsLambdaMCPt", "h2dMCPtVsLambdaMCPt", kTH2D, {axisPt, axisPt});
          histos.add(histodir + "/MC/Sigma0/h2dMCPtVsPhotonMCPt", "h2dMCPtVsPhotonMCPt", kTH2D, {axisPt, axisPt});
          histos.add(histodir + "/MC/Sigma0/h2dMCProcessVsGenRadius", "h2dMCProcessVsGenRadius", kTH2D, {{50, -0.5f, 49.5f}, axisV0PairRadius});
          histos.add(histodir + "/MC/Sigma0/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisSigmaMass});
          histos.add(histodir + "/MC/Sigma0/h3dMCProcess", "h3dMCProcess", kTH3D, {{50, -0.5f, 49.5f}, axisPt, axisSigmaMass});

          histos.add(histodir + "/MC/ASigma0/hPt", "hPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/ASigma0/hMCPt", "hMCPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/ASigma0/hMass", "hMass", kTH1D, {axisSigmaMass});
          histos.add(histodir + "/MC/ASigma0/hMCProcess", "hMCProcess", kTH1D, {{50, -0.5f, 49.5f}});
          histos.add(histodir + "/MC/ASigma0/hGenRadius", "hGenRadius", kTH1D, {axisV0PairRadius});
          histos.add(histodir + "/MC/ASigma0/h2dMCPtVsLambdaMCPt", "h2dMCPtVsLambdaMCPt", kTH2D, {axisPt, axisPt});
          histos.add(histodir + "/MC/ASigma0/h2dMCPtVsPhotonMCPt", "h2dMCPtVsPhotonMCPt", kTH2D, {axisPt, axisPt});
          histos.add(histodir + "/MC/ASigma0/h2dMCProcessVsGenRadius", "h2dMCProcessVsGenRadius", kTH2D, {{50, -0.5f, 49.5f}, axisV0PairRadius});
          histos.add(histodir + "/MC/ASigma0/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisSigmaMass});
          histos.add(histodir + "/MC/ASigma0/h3dMCProcess", "h3dMCProcess", kTH3D, {{50, -0.5f, 49.5f}, axisPt, axisSigmaMass});

          // 1/pT Resolution:
          if (fillResoQAhistos && histodir == "BeforeSel") {
            histos.add(histodir + "/MC/Reso/h3dGammaPtResoVsTPCCR", "h3dGammaPtResoVsTPCCR", kTH3D, {axisInvPt, axisDeltaPt, axisTPCrows});
            histos.add(histodir + "/MC/Reso/h3dGammaPtResoVsTPCCR", "h3dGammaPtResoVsTPCCR", kTH3D, {axisInvPt, axisDeltaPt, axisTPCrows});
            histos.add(histodir + "/MC/Reso/h2dGammaPtResolution", "h2dGammaPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h2dLambdaPtResolution", "h2dLambdaPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h3dLambdaPtResoVsTPCCR", "h3dLambdaPtResoVsTPCCR", kTH3D, {axisInvPt, axisDeltaPt, axisTPCrows});
            histos.add(histodir + "/MC/Reso/h3dLambdaPtResoVsTPCCR", "h3dLambdaPtResoVsTPCCR", kTH3D, {axisInvPt, axisDeltaPt, axisTPCrows});
            histos.add(histodir + "/MC/Reso/h2dAntiLambdaPtResolution", "h2dAntiLambdaPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h3dAntiLambdaPtResoVsTPCCR", "h3dAntiLambdaPtResoVsTPCCR", kTH3D, {axisInvPt, axisDeltaPt, axisTPCrows});
            histos.add(histodir + "/MC/Reso/h3dAntiLambdaPtResoVsTPCCR", "h3dAntiLambdaPtResoVsTPCCR", kTH3D, {axisInvPt, axisDeltaPt, axisTPCrows});
            histos.add(histodir + "/MC/Reso/h2dSigma0PtResolution", "h2dSigma0PtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h2dAntiSigma0PtResolution", "h2dAntiSigma0PtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h2dSigma0RadiusResolution", "h2dSigma0RadiusResolution", kTH2D, {axisPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h2dASigma0RadiusResolution", "h2dASigma0RadiusResolution", kTH2D, {axisPt, axisDeltaPt});
          }

          // For background decomposition study
          if (fillBkgQAhistos) {
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_All", "h2dPtVsMassSigma_All", kTH2D, {axisPt, axisSigmaMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_TrueDaughters", "h2dPtVsMassSigma_TrueDaughters", kTH2D, {axisPt, axisSigmaMass});
            histos.add(histodir + "/MC/BkgStudy/h2dTrueDaughtersMatrix", "h2dTrueDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_TrueGammaFakeLambda", "h2dPtVsMassSigma_TrueGammaFakeLambda", kTH2D, {axisPt, axisSigmaMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_FakeGammaTrueLambda", "h2dPtVsMassSigma_FakeGammaTrueLambda", kTH2D, {axisPt, axisSigmaMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassSigma_FakeDaughters", "h2dPtVsMassSigma_FakeDaughters", kTH2D, {axisPt, axisSigmaMass});
          }
        }
      }

      // Selections
      histos.add("Selection/Photon/hCandidateSel", "hCandidateSel", kTH1D, {axisCandSel});
      histos.add("Selection/Lambda/hCandidateSel", "hCandidateSel", kTH1D, {axisCandSel});

      // For background decomposition study
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

      if (fillSelhistos) {
        for (size_t i = 0; i < PhotonSels.size(); ++i) {
          const auto& sel = PhotonSels[i];

          histos.add(Form("Selection/Photon/h2d%s", sel.c_str()), ("h2d" + sel).c_str(), kTH2D, {axisPt, axisPhotonMass});
          histos.get<TH1>(HIST("Selection/Photon/hCandidateSel"))->GetXaxis()->SetBinLabel(i + 1, sel.c_str());
          histos.add(Form("Selection/Sigma0/h2dPhoton%s", sel.c_str()), ("h2dPhoton" + sel).c_str(), kTH2D, {axisPt, axisSigmaMass});
        }

        for (size_t i = 0; i < LambdaSels.size(); ++i) {
          const auto& sel = LambdaSels[i];

          histos.add(Form("Selection/Lambda/h2d%s", sel.c_str()), ("h2d" + sel).c_str(), kTH2D, {axisPt, axisLambdaMass});
          histos.get<TH1>(HIST("Selection/Lambda/hCandidateSel"))->GetXaxis()->SetBinLabel(i + 1, sel.c_str());
          histos.add(Form("Selection/Sigma0/h2dLambda%s", sel.c_str()), ("h2dLambda" + sel).c_str(), kTH2D, {axisPt, axisSigmaMass});
        }
      }
    }

    if (doprocessPi0RealData || doprocessPi0MonteCarlo) {
      histos.add("Pi0/hMass", "hMass", kTH1D, {axisPi0Mass});
      histos.add("Pi0/hPt", "hPt", kTH1D, {axisPt});
      histos.add("Pi0/hY", "hY", kTH1D, {axisRapidity});
      histos.add("Pi0/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisPi0Mass});
      histos.add("Pi0/h3dMass_MCAssociated", "h3dMass_MCAssociated", kTH3D, {axisCentrality, axisPt, axisPi0Mass});
    }

    if (doprocessGeneratedRun3 || doprocessPi0GeneratedRun3) {

      histos.add("Gen/hGenEvents", "hGenEvents", kTH2D, {{axisNch}, {2, -0.5f, +1.5f}});
      histos.get<TH2>(HIST("Gen/hGenEvents"))->GetYaxis()->SetBinLabel(1, "All gen. events");
      histos.get<TH2>(HIST("Gen/hGenEvents"))->GetYaxis()->SetBinLabel(2, "Gen. with at least 1 rec. events");

      histos.add("Gen/hGenEventCentrality", "hGenEventCentrality", kTH1D, {axisCentrality});
      histos.add("Gen/hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", kTH2D, {axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("Gen/hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", kTH2D, {axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("Gen/hCentralityVsMultMC", "hCentralityVsMultMC", kTH2D, {axisCentrality, axisNch});

      histos.add("Gen/hEventPVzMC", "hEventPVzMC", kTH1D, {{100, -20.0f, +20.0f}});
      histos.add("Gen/hCentralityVsPVzMC", "hCentralityVsPVzMC", kTH2D, {{101, 0.0f, 101.0f}, {100, -20.0f, +20.0f}});

      // Sigma0 specific
      if (doprocessGeneratedRun3) {
        histos.add("Gen/hGenParticlesNoMCColl", "hGenParticlesNoMCColl", kTH1D, {{3, 0.5f, +3.5f}});
        histos.add("Gen/h2dGenSigma0", "h2dGenSigma0", kTH2D, {axisCentrality, axisPt});
        histos.add("Gen/h2dGenAntiSigma0", "h2dGenAntiSigma0", kTH2D, {axisCentrality, axisPt});
        histos.add("Gen/h2dGenSigma0VsMultMC_RecoedEvt", "h2dGenSigma0VsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
        histos.add("Gen/h2dGenAntiSigma0VsMultMC_RecoedEvt", "h2dGenAntiSigma0VsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
        histos.add("Gen/h2dGenSigma0VsMultMC", "h2dGenSigma0VsMultMC", kTH2D, {axisNch, axisPt});
        histos.add("Gen/h2dGenAntiSigma0VsMultMC", "h2dGenAntiSigma0VsMultMC", kTH2D, {axisNch, axisPt});
      }
      if (doprocessPi0GeneratedRun3) { // Pi0 specific
        histos.add("Gen/h2dGenPi0VsMultMC_RecoedEvt", "h2dGenPi0VsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
        histos.add("Gen/h2dGenPi0", "h2dGenPi0", kTH2D, {axisCentrality, axisPt});
        histos.add("Gen/h2dGenPi0VsMultMC", "h2dGenPi0VsMultMC", kTH2D, {axisNch, axisPt});
      }
    }

    // Sigma0 QA
    if (fdoSigma0QA) {
      histos.add("Sigma0QA/h2dAllSigma0CandMCRapVsRecoRap", "h2dAllSigma0CandMCRapVsRecoRap", kTH2D, {axisRapidity, axisRapidity});
      histos.add("Sigma0QA/h2dSigma0MCRapVsRecoRap", "h2dSigma0MCRapVsRecoRap", kTH2D, {axisRapidity, axisRapidity});
      histos.add("Sigma0QA/h2dASigma0MCRapVsRecoRap", "h2dASigma0MCRapVsRecoRap", kTH2D, {axisRapidity, axisRapidity});

      histos.add("Sigma0QA/hDuplicates", "hDuplicates", kTH1D, {{10, -0.5f, +9.5f}});
      histos.add("Sigma0QA/hSigma0Duplicates", "hSigma0Duplicates", kTH1D, {{10, -0.5f, +9.5f}});
      histos.add("Sigma0QA/hASigma0Duplicates", "hASigma0Duplicates", kTH1D, {{10, -0.5f, +9.5f}});
    }

    // inspect histogram sizes, please
    histos.print();
  }

  // ______________________________________________________
  // Check whether the collision passes our collision selections
  // Should work with collisions, mccollisions, stracollisions and stramccollisions tables!
  template <typename TCollision>
  bool IsEventAccepted(TCollision const& collision, bool fillHists)
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
    float interactionRate = (fGetIR) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource, fIRCrashOnNull) * 1.e-3 : -1;
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    if (fGetIR) {
      if (interactionRate < 0)
        histos.get<TH1>(HIST("GeneralQA/hRunNumberNegativeIR"))->Fill(Form("%d", collision.runNumber()), 1); // This lists all run numbers without IR info!

      histos.fill(HIST("GeneralQA/hInteractionRate"), interactionRate);
      histos.fill(HIST("GeneralQA/hCentralityVsInteractionRate"), centrality, interactionRate);
    }

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

    // Fill centrality histogram after event selection
    if (fillHists)
      histos.fill(HIST("hEventCentrality"), centrality);

    return true;
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
  template <typename TMCCollisions, typename TCollisions, typename TGenParticles>
  void analyzeGenerated(TMCCollisions const& mcCollisions, TCollisions const& collisions, TGenParticles const& genParticles)
  {
    fillGeneratedEventProperties(mcCollisions, collisions);
    std::vector<int> listBestCollisionIdx = getListOfRecoCollIndices(mcCollisions, collisions);

    for (auto& genParticle : genParticles) {
      float centrality = 100.5f;

      // Has MC collision
      if (!genParticle.has_straMCCollision()) {
        histos.fill(HIST("Gen/hGenParticlesNoMCColl"), 1);
        continue;
      }

      // Rapidity selection
      if ((genParticle.mcy() < genSelections.mc_rapidityMin) || (genParticle.mcy() > genSelections.mc_rapidityMax))
        continue;

      // Selection on the source (generator/transport)
      if (!genParticle.producedByGenerator() && genSelections.mc_keepOnlyFromGenerator)
        continue;

      // Select corresponding mc collision && Basic event selection
      auto mcCollision = genParticle.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
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

      //______________________________________________________________________________
      // Generated Sigma0 processing
      if constexpr (requires { genParticle.isSigma0(); }) {
        if (doprocessGeneratedRun3) {

          float ptmc = genParticle.mcpt();

          if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
            auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
            centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

            if (genParticle.isSigma0())
              histos.fill(HIST("Gen/h2dGenSigma0VsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);

            else
              histos.fill(HIST("Gen/h2dGenAntiSigma0VsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
          }

          if (genParticle.isSigma0()) {
            histos.fill(HIST("Gen/h2dGenSigma0"), centrality, ptmc);
            histos.fill(HIST("Gen/h2dGenSigma0VsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
          } else {
            histos.fill(HIST("Gen/h2dGenAntiSigma0"), centrality, ptmc);
            histos.fill(HIST("Gen/h2dGenAntiSigma0VsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
          }
        }
      }
      //______________________________________________________________________________
      // Generated Pi0 processing
      if (doprocessPi0GeneratedRun3) {
        float ptmc = genParticle.mcpt();

        if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
          auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
          centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
          histos.fill(HIST("Gen/h2dGenPi0VsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }

        histos.fill(HIST("Gen/h2dGenPi0"), centrality, ptmc);
        histos.fill(HIST("Gen/h2dGenPi0VsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    }
  }

  //__________________________________________
  template <bool isGamma, typename TSigma0Object>
  int retrieveV0TrackCode(TSigma0Object const& sigma)
  {

    int TrkCode = 10; // 1: TPC-only, 2: TPC+Something, 3: ITS-Only, 4: ITS+TPC + Something, 10: anything else

    if (isGamma) {
      if (sigma.photonPosTrackCode() == 1 && sigma.photonNegTrackCode() == 1)
        TrkCode = 1;
      if ((sigma.photonPosTrackCode() != 1 && sigma.photonNegTrackCode() == 1) || (sigma.photonPosTrackCode() == 1 && sigma.photonNegTrackCode() != 1))
        TrkCode = 2;
      if (sigma.photonPosTrackCode() == 3 && sigma.photonNegTrackCode() == 3)
        TrkCode = 3;
      if (sigma.photonPosTrackCode() == 2 || sigma.photonNegTrackCode() == 2)
        TrkCode = 4;
    } else {
      if (sigma.lambdaPosTrackCode() == 1 && sigma.lambdaNegTrackCode() == 1)
        TrkCode = 1;
      if ((sigma.lambdaPosTrackCode() != 1 && sigma.lambdaNegTrackCode() == 1) || (sigma.lambdaPosTrackCode() == 1 && sigma.lambdaNegTrackCode() != 1))
        TrkCode = 2;
      if (sigma.lambdaPosTrackCode() == 3 && sigma.lambdaNegTrackCode() == 3)
        TrkCode = 3;
      if (sigma.lambdaPosTrackCode() == 2 || sigma.lambdaNegTrackCode() == 2)
        TrkCode = 4;
    }

    return TrkCode;
  }

  template <typename TSigma0Object>
  void getResolution(TSigma0Object const& sigma)
  {

    //_______________________________________
    // Gamma MC association
    if (sigma.photonPDGCode() == 22) {
      if (sigma.photonmcpt() > 0) {
        histos.fill(HIST("BeforeSel/MC/Reso/h3dGammaPtResoVsTPCCR"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt(), -1 * sigma.photonNegTPCCrossedRows()); // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/Reso/h3dGammaPtResoVsTPCCR"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt(), sigma.photonPosTPCCrossedRows());      // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/Reso/h2dGammaPtResolution"), 1.f / sigma.photonmcpt(), 1.f / sigma.photonPt() - 1.f / sigma.photonmcpt());                                        // pT resolution
      }
    }

    //_______________________________________
    // Lambda MC association
    if (sigma.lambdaPDGCode() == 3122) {
      if (sigma.lambdamcpt() > 0) {
        histos.fill(HIST("BeforeSel/MC/Reso/h2dLambdaPtResolution"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt());                                        // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/Reso/h3dLambdaPtResoVsTPCCR"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt(), -1 * sigma.lambdaNegTPCCrossedRows()); // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/Reso/h3dLambdaPtResoVsTPCCR"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt(), sigma.lambdaPosTPCCrossedRows());      // 1/pT resolution
      }
    }

    //_______________________________________
    // AntiLambda MC association
    if (sigma.lambdaPDGCode() == -3122) {
      if (sigma.lambdamcpt() > 0) {
        histos.fill(HIST("BeforeSel/MC/Reso/h2dAntiLambdaPtResolution"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt());                                        // pT resolution
        histos.fill(HIST("BeforeSel/MC/Reso/h3dAntiLambdaPtResoVsTPCCR"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt(), -1 * sigma.lambdaNegTPCCrossedRows()); // 1/pT resolution
        histos.fill(HIST("BeforeSel/MC/Reso/h3dAntiLambdaPtResoVsTPCCR"), 1.f / sigma.lambdamcpt(), 1.f / sigma.lambdaPt() - 1.f / sigma.lambdamcpt(), sigma.lambdaPosTPCCrossedRows());      // 1/pT resolution
      }
    }

    //_______________________________________
    // Sigma and AntiSigma MC association
    if (sigma.isSigma0()) {
      histos.fill(HIST("BeforeSel/MC/Reso/h2dSigma0RadiusResolution"), sigma.mcpt(), sigma.radius() - sigma.mcradius()); // pT resolution
      if (sigma.mcpt() > 0)
        histos.fill(HIST("BeforeSel/MC/Reso/h2dSigma0PtResolution"), 1.f / sigma.mcpt(), 1.f / sigma.pt() - 1.f / sigma.mcpt()); // pT resolution
    }
    if (sigma.isAntiSigma0()) {
      histos.fill(HIST("BeforeSel/MC/Reso/h2dASigma0RadiusResolution"), sigma.mcpt(), sigma.radius() - sigma.mcradius()); // pT resolution
      if (sigma.mcpt() > 0)
        histos.fill(HIST("BeforeSel/MC/Reso/h2dAntiSigma0PtResolution"), 1.f / sigma.mcpt(), 1.f / sigma.pt() - 1.f / sigma.mcpt()); // pT resolution
    }
  }

  // To save histograms for background analysis
  template <int mode, typename TSigma0Object>
  void runBkgAnalysis(TSigma0Object const& sigma)
  {
    // Check whether it is before or after selections
    static constexpr std::string_view MainDir[] = {"BeforeSel", "AfterSel"};

    bool fIsSigma = sigma.isSigma0();
    bool fIsAntiSigma = sigma.isAntiSigma0();
    int PhotonPDGCode = sigma.photonPDGCode();
    int PhotonPDGCodeMother = sigma.photonPDGCodeMother();
    int LambdaPDGCode = sigma.lambdaPDGCode();
    int LambdaPDGCodeMother = sigma.lambdaPDGCodeMother();
    float sigmapT = sigma.pt();
    float sigmaMass = sigma.sigma0Mass();

    histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_All"), sigmapT, sigmaMass);

    //_______________________________________
    // Real Gamma x Real Lambda - but not from the same sigma0/antisigma0!
    if ((PhotonPDGCode == 22) && ((LambdaPDGCode == 3122) || (LambdaPDGCode == -3122)) && (!fIsSigma && !fIsAntiSigma)) {
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_TrueDaughters"), sigmapT, sigmaMass);
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dTrueDaughtersMatrix"), LambdaPDGCodeMother, PhotonPDGCodeMother);
    }

    //_______________________________________
    // Real Gamma x fake Lambda
    if ((PhotonPDGCode == 22) && (LambdaPDGCode != 3122) && (LambdaPDGCode != -3122))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_TrueGammaFakeLambda"), sigmapT, sigmaMass);

    //_______________________________________
    // Fake Gamma x Real Lambda
    if ((PhotonPDGCode != 22) && ((LambdaPDGCode == 3122) || (LambdaPDGCode == -3122)))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_FakeGammaTrueLambda"), sigmapT, sigmaMass);

    //_______________________________________
    // Fake Gamma x Fake Lambda
    if ((PhotonPDGCode != 22) && (LambdaPDGCode != 3122) && (LambdaPDGCode != -3122))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassSigma_FakeDaughters"), sigmapT, sigmaMass);
  }

  template <int mode, typename TSigma0Object, typename TCollision>
  void fillHistos(TSigma0Object const& sigma, TCollision const& collision)
  {

    // Check whether it is before or after selections
    static constexpr std::string_view MainDir[] = {"BeforeSel", "AfterSel"};

    // Get V0trackCode
    int GammaTrkCode = retrieveV0TrackCode<true>(sigma);
    int LambdaTrkCode = retrieveV0TrackCode<false>(sigma);

    float photonRZLineCut = TMath::Abs(sigma.photonZconv()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-photonSelections.PhotonMaxDauEta))) - photonSelections.PhotonLineCutZ0;
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    if (fillSelhistos) {
      //_______________________________________
      // Photon
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hTrackCode"), GammaTrkCode);
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hV0Type"), sigma.photonV0Type());

      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCANegToPV"), sigma.photonDCANegPV());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCAPosToPV"), sigma.photonDCAPosPV());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCADau"), sigma.photonDCADau());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosTPCCR"), sigma.photonPosTPCCrossedRows());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegTPCCR"), sigma.photonNegTPCCrossedRows());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosTPCNSigmaEl"), sigma.photonPosTPCNSigmaEl());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegTPCNSigmaEl"), sigma.photonNegTPCNSigmaEl());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hpT"), sigma.photonPt());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hY"), sigma.photonY());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosEta"), sigma.photonPosEta());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegEta"), sigma.photonNegEta());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hRadius"), sigma.photonRadius());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hZ"), sigma.photonZconv());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dRZCut"), sigma.photonRadius(), photonRZLineCut);
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dRZPlane"), sigma.photonZconv(), sigma.photonRadius());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hCosPA"), sigma.photonCosPA());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPsiPair"), sigma.photonPsiPair());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPhi"), sigma.photonPhi());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h3dMass"), centrality, sigma.photonPt(), sigma.photonMass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hMass"), sigma.photonMass());

      //_______________________________________
      // Lambdas
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hTrackCode"), LambdaTrkCode);
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hRadius"), sigma.lambdaRadius());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hDCADau"), sigma.lambdaDCADau());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hCosPA"), sigma.lambdaCosPA());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hY"), sigma.lambdaY());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosEta"), sigma.lambdaPosEta());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegEta"), sigma.lambdaNegEta());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosTPCCR"), sigma.lambdaPosTPCCrossedRows());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegTPCCR"), sigma.lambdaNegTPCCrossedRows());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosITSCls"), sigma.lambdaPosITSCls());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegITSCls"), sigma.lambdaNegITSCls());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hPosChi2PerNc"), sigma.lambdaPosChi2PerNcl());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hNegChi2PerNc"), sigma.lambdaNegChi2PerNcl());
      histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLifeTime"), sigma.lambdaLifeTime());
    }

    //_______________________________________
    // Sigmas and Lambdas
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenteros"), sigma.photonAlpha(), sigma.photonQt());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenteros"), sigma.lambdaAlpha(), sigma.lambdaQt());

    if (sigma.lambdaAlpha() > 0) {
      if (fillSelhistos) {
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_LambdaPr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_LambdaPi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdaDCANegToPV"), sigma.lambdaDCANegPV());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdaDCAPosToPV"), sigma.lambdaDCAPosPV());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdapT"), sigma.lambdaPt());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hLambdaMass"), sigma.lambdaMass());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h3dLambdaMass"), centrality, sigma.lambdaPt(), sigma.lambdaMass());
      }

      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hMass"), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hPt"), sigma.pt());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hY"), sigma.sigma0Y());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hRadius"), sigma.radius());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h2dRadiusVspT"), sigma.radius(), sigma.pt());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/hDCAPairDau"), sigma.dcadaughters());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h3dMass"), centrality, sigma.pt(), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h3dPhotonYMass"), sigma.photonY(), sigma.pt(), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h3dPhotonRadiusMass"), sigma.photonRadius(), sigma.pt(), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/Sigma0/h3dOPAngleVsMass"), sigma.opAngle(), sigma.pt(), sigma.sigma0Mass());
    } else {
      if (fillSelhistos) {
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_ALambdaPr"), sigma.lambdaNegPrTPCNSigma(), sigma.aLambdaPrTOFNSigma());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h2dTPCvsTOFNSigma_ALambdaPi"), sigma.lambdaPosPiTPCNSigma(), sigma.aLambdaPiTOFNSigma());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hALambdaDCANegToPV"), sigma.lambdaDCANegPV());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hALambdaDCAPosToPV"), sigma.lambdaDCAPosPV());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hALambdapT"), sigma.lambdaPt());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/hAntiLambdaMass"), sigma.antilambdaMass());
        histos.fill(HIST(MainDir[mode]) + HIST("/Lambda/h3dAntiLambdaMass"), centrality, sigma.lambdaPt(), sigma.antilambdaMass());
      }

      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hMass"), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hPt"), sigma.pt());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hY"), sigma.sigma0Y());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hRadius"), sigma.radius());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h2dRadiusVspT"), sigma.radius(), sigma.pt());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/hDCAPairDau"), sigma.dcadaughters());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h3dMass"), centrality, sigma.pt(), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h3dPhotonYMass"), sigma.photonY(), sigma.pt(), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h3dPhotonRadiusMass"), sigma.photonRadius(), sigma.pt(), sigma.sigma0Mass());
      histos.fill(HIST(MainDir[mode]) + HIST("/ASigma0/h3dOPAngleVsMass"), sigma.opAngle(), sigma.pt(), sigma.sigma0Mass());
    }

    //_______________________________________
    // MC specific
    if (doprocessMonteCarlo) {
      if constexpr (requires { sigma.lambdaPDGCode(); sigma.photonPDGCode(); }) {

        if (fillSelhistos) {
          //_______________________________________
          // Gamma MC association
          if (sigma.photonPDGCode() == 22) {
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hV0ToCollAssoc"), sigma.photonIsCorrectlyAssoc());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hPt"), sigma.photonPt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hMCPt"), sigma.photonmcpt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hPosTPCNSigmaEl"), sigma.photonPosTPCNSigmaEl());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hNegTPCNSigmaEl"), sigma.photonNegTPCNSigmaEl());

            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dPAVsPt"), TMath::ACos(sigma.photonCosPA()), sigma.photonmcpt());

            if (!sigma.photonIsCorrectlyAssoc()) {
              histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hPt_BadCollAssig"), sigma.photonmcpt());
              histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dPAVsPt_BadCollAssig"), TMath::ACos(sigma.photonCosPA()), sigma.photonmcpt());
            }
          }

          //_______________________________________
          // Lambda MC association
          if (sigma.lambdaPDGCode() == 3122) {
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/hV0ToCollAssoc"), sigma.lambdaIsCorrectlyAssoc());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/hPt"), sigma.lambdaPt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/hMCPt"), sigma.lambdamcpt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/h3dTPCvsTOFNSigma_Pr"), sigma.lambdaPosPrTPCNSigma(), sigma.lambdaPrTOFNSigma(), sigma.lambdaPt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Lambda/h3dTPCvsTOFNSigma_Pi"), sigma.lambdaNegPiTPCNSigma(), sigma.lambdaPiTOFNSigma(), sigma.lambdaPt());
          }

          //_______________________________________
          // AntiLambda MC association
          if (sigma.lambdaPDGCode() == -3122) {
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/hV0ToCollAssoc"), sigma.lambdaIsCorrectlyAssoc());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/hPt"), sigma.lambdaPt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/hMCPt"), sigma.lambdamcpt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/h3dTPCvsTOFNSigma_Pr"), sigma.lambdaNegPrTPCNSigma(), sigma.aLambdaPrTOFNSigma(), sigma.lambdaPt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/ALambda/h3dTPCvsTOFNSigma_Pi"), sigma.lambdaPosPiTPCNSigma(), sigma.aLambdaPiTOFNSigma(), sigma.lambdaPt());
          }
        }
        //_______________________________________
        // Sigma0 MC association
        if (sigma.isSigma0()) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.photonAlpha(), sigma.photonQt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.lambdaAlpha(), sigma.lambdaQt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hPt"), sigma.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hMCPt"), sigma.mcpt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h2dMCPtVsLambdaMCPt"), sigma.mcpt(), sigma.lambdamcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h2dMCPtVsPhotonMCPt"), sigma.mcpt(), sigma.photonmcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hMass"), sigma.sigma0Mass());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h3dMass"), centrality, sigma.mcpt(), sigma.sigma0Mass());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hMCProcess"), sigma.mcprocess());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/hGenRadius"), sigma.mcradius());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h2dMCProcessVsGenRadius"), sigma.mcprocess(), sigma.mcradius());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Sigma0/h3dMCProcess"), sigma.mcprocess(), sigma.mcpt(), sigma.sigma0Mass());
        }

        //_______________________________________
        // AntiSigma0 MC association
        if (sigma.isAntiSigma0()) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.photonAlpha(), sigma.photonQt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), sigma.lambdaAlpha(), sigma.lambdaQt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hPt"), sigma.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hMCPt"), sigma.mcpt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h2dMCPtVsLambdaMCPt"), sigma.mcpt(), sigma.lambdamcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h2dMCPtVsPhotonMCPt"), sigma.mcpt(), sigma.photonmcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hMass"), sigma.sigma0Mass());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h3dMass"), centrality, sigma.mcpt(), sigma.sigma0Mass());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hMCProcess"), sigma.mcprocess());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/hGenRadius"), sigma.mcradius());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h2dMCProcessVsGenRadius"), sigma.mcprocess(), sigma.mcradius());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/ASigma0/h3dMCProcess"), sigma.mcprocess(), sigma.mcpt(), sigma.sigma0Mass());
        }

        // For background studies:
        if (fillBkgQAhistos)
          runBkgAnalysis<mode>(sigma);

        //_______________________________________
        // pT resolution histos
        if ((mode == 0) && fillResoQAhistos)
          getResolution(sigma);
      }
    }
  }

  template <int selection_index, typename TSigma0Object>
  void fillSelHistos(TSigma0Object const& sigma, int PDGRequired)
  {

    static constexpr std::string_view PhotonSelsLocal[] = {"NoSel", "V0Type", "DCADauToPV",
                                                           "DCADau", "DauTPCCR", "TPCNSigmaEl", "V0pT",
                                                           "Y", "V0Radius", "RZCut", "Armenteros", "CosPA",
                                                           "PsiPair", "Phi", "Mass"};

    static constexpr std::string_view LambdaSelsLocal[] = {"NoSel", "V0Radius", "DCADau", "Armenteros",
                                                           "CosPA", "Y", "TPCCR", "DauITSCls", "Lifetime",
                                                           "TPCTOFPID", "DCADauToPV", "Mass"};

    if (PDGRequired == 22) {
      if constexpr (selection_index >= 0 && selection_index < static_cast<int>(std::size(PhotonSelsLocal))) {
        histos.fill(HIST("Selection/Photon/hCandidateSel"), selection_index);
        histos.fill(HIST("Selection/Photon/h2d") + HIST(PhotonSelsLocal[selection_index]), sigma.photonPt(), sigma.photonMass());
        histos.fill(HIST("Selection/Sigma0/h2dPhoton") + HIST(PhotonSelsLocal[selection_index]), sigma.pt(), sigma.sigma0Mass());
      }
    }

    if (PDGRequired == 3122) {
      if constexpr (selection_index >= 0 && selection_index < static_cast<int>(std::size(LambdaSelsLocal))) {
        histos.fill(HIST("Selection/Lambda/hCandidateSel"), selection_index);
        histos.fill(HIST("Selection/Lambda/h2d") + HIST(LambdaSelsLocal[selection_index]), sigma.lambdaPt(), sigma.lambdaMass());
        histos.fill(HIST("Selection/Sigma0/h2dLambda") + HIST(LambdaSelsLocal[selection_index]), sigma.pt(), sigma.sigma0Mass());
      }
    }
  }

  // Sigma0 QA analysis function
  template <typename TSigma0s>
  void doSigma0QA(TSigma0s const& fullSigma0s)
  {
    std::unordered_map<int, int> sigma0Counts;
    std::unordered_map<int, int> sigma0Index;

    for (const auto& sigma0 : fullSigma0s) {
      // Crosscheck reco rapidity and MC rapidity
      histos.fill(HIST("Sigma0QA/h2dAllSigma0CandMCRapVsRecoRap"), sigma0.sigma0MCY(), sigma0.sigma0Y());

      if (sigma0.isSigma0())
        histos.fill(HIST("Sigma0QA/h2dSigma0MCRapVsRecoRap"), sigma0.sigma0MCY(), sigma0.sigma0Y());
      if (sigma0.isAntiSigma0())
        histos.fill(HIST("Sigma0QA/h2dASigma0MCRapVsRecoRap"), sigma0.sigma0MCY(), sigma0.sigma0Y());

      // grouping per mcparticle
      if (sigma0.has_mcParticle()) {
        sigma0Counts[sigma0.mcParticleId()] += 1;                  // duplicate count
        sigma0Index[sigma0.mcParticleId()] = sigma0.globalIndex(); // saving index
      }
    }
    for (const auto& [mcid, NDuplicates] : sigma0Counts) {
      auto sigma0mc = fullSigma0s.rawIteratorAt(sigma0Index[mcid]);
      histos.fill(HIST("Sigma0QA/hDuplicates"), NDuplicates); // how many times a mc sigma0 was reconstructed

      if (sigma0mc.isSigma0())
        histos.fill(HIST("Sigma0QA/hSigma0Duplicates"), NDuplicates); // how many times a mc sigma0 was reconstructed

      if (sigma0mc.isAntiSigma0())
        histos.fill(HIST("Sigma0QA/hASigma0Duplicates"), NDuplicates); // how many times a mc sigma0 was reconstructed
    }
  }

  // Apply specific selections for photons
  template <typename TV0Object>
  bool selectPhoton(TV0Object const& cand)
  {
    fillSelHistos<0>(cand, 22);
    if (cand.photonV0Type() != photonSelections.Photonv0TypeSel && photonSelections.Photonv0TypeSel > -1)
      return false;

    fillSelHistos<1>(cand, 22);
    if ((TMath::Abs(cand.photonDCAPosPV()) < photonSelections.PhotonMinDCADauToPv) || (TMath::Abs(cand.photonDCANegPV()) < photonSelections.PhotonMinDCADauToPv))
      return false;

    fillSelHistos<2>(cand, 22);
    if (TMath::Abs(cand.photonDCADau()) > photonSelections.PhotonMaxDCAV0Dau)
      return false;

    fillSelHistos<3>(cand, 22);
    if ((cand.photonPosTPCCrossedRows() < photonSelections.PhotonMinTPCCrossedRows) || (cand.photonNegTPCCrossedRows() < photonSelections.PhotonMinTPCCrossedRows))
      return false;

    fillSelHistos<4>(cand, 22);
    if (((cand.photonPosTPCNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) || (cand.photonPosTPCNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)))
      return false;

    if (((cand.photonNegTPCNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) || (cand.photonNegTPCNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)))
      return false;

    fillSelHistos<5>(cand, 22);
    if ((cand.photonPt() < photonSelections.PhotonMinPt) || (cand.photonPt() > photonSelections.PhotonMaxPt))
      return false;

    fillSelHistos<6>(cand, 22);
    if ((cand.photonY() < photonSelections.PhotonMinRapidity) || (cand.photonY() > photonSelections.PhotonMaxRapidity))
      return false;
    if ((cand.photonPosEta() < photonSelections.PhotonMinDauEta) || (cand.photonPosEta() > photonSelections.PhotonMaxDauEta))
      return false;
    if ((cand.photonNegEta() < photonSelections.PhotonMinDauEta) || (cand.photonNegEta() > photonSelections.PhotonMaxDauEta))
      return false;

    fillSelHistos<7>(cand, 22);
    if ((cand.photonRadius() < photonSelections.PhotonMinRadius) || (cand.photonRadius() > photonSelections.PhotonMaxRadius))
      return false;

    fillSelHistos<8>(cand, 22);
    if ((cand.photonZconv() < photonSelections.PhotonMinZ) || (cand.photonZconv() > photonSelections.PhotonMaxZ))
      return false;

    float photonRZLineCut = TMath::Abs(cand.photonZconv()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-photonSelections.PhotonMaxDauEta))) - photonSelections.PhotonLineCutZ0;
    if ((TMath::Abs(cand.photonRadius()) < photonRZLineCut))
      return false;

    fillSelHistos<9>(cand, 22);
    if (cand.photonQt() > photonSelections.PhotonMaxQt)
      return false;

    if (TMath::Abs(cand.photonAlpha()) > photonSelections.PhotonMaxAlpha)
      return false;

    fillSelHistos<10>(cand, 22);
    if (cand.photonCosPA() < photonSelections.PhotonMinV0cospa)
      return false;

    fillSelHistos<11>(cand, 22);
    if (TMath::Abs(cand.photonPsiPair()) > photonSelections.PhotonPsiPairMax)
      return false;

    fillSelHistos<12>(cand, 22);
    if ((((cand.photonPhi() > photonSelections.PhotonPhiMin1) && (cand.photonPhi() < photonSelections.PhotonPhiMax1)) || ((cand.photonPhi() > photonSelections.PhotonPhiMin2) && (cand.photonPhi() < photonSelections.PhotonPhiMax2))) && ((photonSelections.PhotonPhiMin1 != -1) && (photonSelections.PhotonPhiMax1 != -1) && (photonSelections.PhotonPhiMin2 != -1) && (photonSelections.PhotonPhiMax2 != -1)))
      return false;

    fillSelHistos<13>(cand, 22);
    if (TMath::Abs(cand.photonMass()) > photonSelections.PhotonMaxMass)
      return false;

    fillSelHistos<14>(cand, 22);
    return true;
  }

  // Apply specific selections for lambdas
  template <typename TV0Object>
  bool selectLambda(TV0Object const& cand)
  {
    fillSelHistos<0>(cand, 3122);
    if ((cand.lambdaRadius() < lambdaSelections.LambdaMinv0radius) || (cand.lambdaRadius() > lambdaSelections.LambdaMaxv0radius))
      return false;

    fillSelHistos<1>(cand, 3122);
    if (TMath::Abs(cand.lambdaDCADau()) > lambdaSelections.LambdaMaxDCAV0Dau)
      return false;

    fillSelHistos<2>(cand, 3122);
    if ((cand.lambdaQt() < lambdaSelections.LambdaMinQt) || (cand.lambdaQt() > lambdaSelections.LambdaMaxQt))
      return false;

    if ((TMath::Abs(cand.lambdaAlpha()) < lambdaSelections.LambdaMinAlpha) || (TMath::Abs(cand.lambdaAlpha()) > lambdaSelections.LambdaMaxAlpha))
      return false;

    fillSelHistos<3>(cand, 3122);
    if (cand.lambdaCosPA() < lambdaSelections.LambdaMinv0cospa)
      return false;

    fillSelHistos<4>(cand, 3122);
    if ((cand.lambdaY() < lambdaSelections.LambdaMinRapidity) || (cand.lambdaY() > lambdaSelections.LambdaMaxRapidity))
      return false;
    if ((cand.lambdaPosEta() < lambdaSelections.LambdaMinDauEta) || (cand.lambdaPosEta() > lambdaSelections.LambdaMaxDauEta))
      return false;
    if ((cand.lambdaNegEta() < lambdaSelections.LambdaMinDauEta) || (cand.lambdaNegEta() > lambdaSelections.LambdaMaxDauEta))
      return false;

    fillSelHistos<5>(cand, 3122);
    if ((cand.lambdaPosTPCCrossedRows() < lambdaSelections.LambdaMinTPCCrossedRows) || (cand.lambdaNegTPCCrossedRows() < lambdaSelections.LambdaMinTPCCrossedRows))
      return false;

    fillSelHistos<6>(cand, 3122);
    // check minimum number of ITS clusters + reject ITS afterburner tracks if requested
    bool posIsFromAfterburner = cand.lambdaPosChi2PerNcl() < 0;
    bool negIsFromAfterburner = cand.lambdaNegChi2PerNcl() < 0;
    if (cand.lambdaPosITSCls() < lambdaSelections.LambdaMinITSclusters && (!lambdaSelections.LambdaRejectPosITSafterburner || posIsFromAfterburner))
      return false;
    if (cand.lambdaNegITSCls() < lambdaSelections.LambdaMinITSclusters && (!lambdaSelections.LambdaRejectNegITSafterburner || negIsFromAfterburner))
      return false;

    fillSelHistos<7>(cand, 3122);
    if (cand.lambdaLifeTime() > lambdaSelections.LambdaMaxLifeTime)
      return false;

    // Separating lambda and antilambda selections:
    fillSelHistos<8>(cand, 3122);
    if (cand.lambdaAlpha() > 0) { // Lambda selection
      // TPC Selection
      if (lambdaSelections.fselLambdaTPCPID && (TMath::Abs(cand.lambdaPosPrTPCNSigma()) > lambdaSelections.LambdaMaxTPCNSigmas))
        return false;
      if (lambdaSelections.fselLambdaTPCPID && (TMath::Abs(cand.lambdaNegPiTPCNSigma()) > lambdaSelections.LambdaMaxTPCNSigmas))
        return false;

      // TOF Selection
      if (lambdaSelections.fselLambdaTOFPID && (TMath::Abs(cand.lambdaPrTOFNSigma()) > lambdaSelections.LambdaPrMaxTOFNSigmas))
        return false;
      if (lambdaSelections.fselLambdaTOFPID && (TMath::Abs(cand.lambdaPiTOFNSigma()) > lambdaSelections.LambdaPiMaxTOFNSigmas))
        return false;

      // DCA Selection
      fillSelHistos<9>(cand, 3122);
      if ((TMath::Abs(cand.lambdaDCAPosPV()) < lambdaSelections.LambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < lambdaSelections.LambdaMinDCANegToPv))
        return false;

      // Mass Selection
      fillSelHistos<10>(cand, 3122);
      if (TMath::Abs(cand.lambdaMass() - o2::constants::physics::MassLambda0) > lambdaSelections.LambdaWindow)
        return false;

      fillSelHistos<11>(cand, 3122);

    } else { // AntiLambda selection

      // TPC Selection
      if (lambdaSelections.fselLambdaTPCPID && (TMath::Abs(cand.lambdaPosPiTPCNSigma()) > lambdaSelections.LambdaMaxTPCNSigmas))
        return false;
      if (lambdaSelections.fselLambdaTPCPID && (TMath::Abs(cand.lambdaNegPrTPCNSigma()) > lambdaSelections.LambdaMaxTPCNSigmas))
        return false;

      // TOF Selection
      if (lambdaSelections.fselLambdaTOFPID && (TMath::Abs(cand.aLambdaPrTOFNSigma()) > lambdaSelections.LambdaPrMaxTOFNSigmas))
        return false;
      if (lambdaSelections.fselLambdaTOFPID && (TMath::Abs(cand.aLambdaPiTOFNSigma()) > lambdaSelections.LambdaPiMaxTOFNSigmas))
        return false;

      // DCA Selection
      fillSelHistos<9>(cand, 3122);
      if ((TMath::Abs(cand.lambdaDCAPosPV()) < lambdaSelections.ALambdaMinDCAPosToPv) || (TMath::Abs(cand.lambdaDCANegPV()) < lambdaSelections.ALambdaMinDCANegToPv))
        return false;

      // Mass Selection
      fillSelHistos<10>(cand, 3122);
      if (TMath::Abs(cand.antilambdaMass() - o2::constants::physics::MassLambda0) > lambdaSelections.LambdaWindow)
        return false;

      fillSelHistos<11>(cand, 3122);
    }

    return true;
  }

  // Apply selections in sigma0 candidates
  template <typename TSigma0Object>
  bool processSigma0Candidate(TSigma0Object const& cand)
  {
    // Photon specific selections
    if (!selectPhoton(cand))
      return false;

    // Lambda specific selections
    if (!selectLambda(cand))
      return false;

    // Sigma0 specific selections
    float sigma0Y = cand.sigma0Y();
    if constexpr (requires { cand.sigma0MCY(); }) { // If MC
      sigma0Y = cand.sigma0MCY();
    }

    // Rapidity
    if ((sigma0Y < sigma0Selections.Sigma0MinRapidity) || (sigma0Y > sigma0Selections.Sigma0MaxRapidity))
      return false;

    // V0Pair Radius
    if (cand.radius() > sigma0Selections.Sigma0MaxRadius)
      return false;

    // DCA V0Pair Daughters
    if (cand.dcadaughters() > sigma0Selections.Sigma0MaxDCADau)
      return false;

    // Opening Angle
    if (cand.opAngle() > sigma0Selections.Sigma0MaxOPAngle)
      return false;

    return true;
  }

  // Main analysis function
  template <typename TCollisions, typename TSigma0s>
  void analyzeRecoeSigma0s(TCollisions const& collisions, TSigma0s const& fullSigma0s)
  {
    // Custom grouping
    std::vector<std::vector<int>> sigma0grouped(collisions.size());

    for (const auto& sigma0 : fullSigma0s) {
      sigma0grouped[sigma0.straCollisionId()].push_back(sigma0.globalIndex());
    }

    // Collisions loop
    for (const auto& coll : collisions) {

      // Event selection
      if (!IsEventAccepted(coll, true))
        continue;

      // Sigma0s loop
      for (size_t i = 0; i < sigma0grouped[coll.globalIndex()].size(); i++) {
        auto sigma0 = fullSigma0s.rawIteratorAt(sigma0grouped[coll.globalIndex()][i]);

        // if MC
        if constexpr (requires { sigma0.isSigma0(); sigma0.isAntiSigma0(); }) {
          if (doMCAssociation && !(sigma0.isSigma0() || sigma0.isAntiSigma0()))
            continue;

          if (selRecoFromGenerator && !sigma0.isProducedByGenerator())
            continue;
        }

        // Fill histos before any selection
        fillHistos<0>(sigma0, coll);

        // Select sigma0 candidates
        if (!processSigma0Candidate(sigma0))
          continue;

        // Fill histos after all selections
        fillHistos<1>(sigma0, coll);
      }
    }

    // Optionally run QA analysis for reco sigma0
    // if (fdoSigma0QA) doSigma0QA(fullSigma0s); // TODO: improve this to run sigma0 QA
  }

  // Apply selections in sigma0 candidates
  template <typename TPi0Object>
  bool processPi0Candidate(TPi0Object const& cand)
  {
    if ((cand.photon1V0Type() != photonSelections.Photonv0TypeSel || cand.photon2V0Type() != photonSelections.Photonv0TypeSel) && photonSelections.Photonv0TypeSel > -1)
      return false;

    if ((TMath::Abs(cand.photon1DCAPosPV()) < photonSelections.PhotonMinDCADauToPv) ||
        (TMath::Abs(cand.photon2DCAPosPV()) < photonSelections.PhotonMinDCADauToPv) ||
        (TMath::Abs(cand.photon1DCANegPV()) < photonSelections.PhotonMinDCADauToPv) ||
        (TMath::Abs(cand.photon2DCANegPV()) < photonSelections.PhotonMinDCADauToPv))
      return false;

    if ((TMath::Abs(cand.photon1DCADau()) > photonSelections.PhotonMaxDCAV0Dau) || (TMath::Abs(cand.photon2DCADau()) > photonSelections.PhotonMaxDCAV0Dau))
      return false;

    if ((cand.photon1PosTPCCrossedRows() < photonSelections.PhotonMinTPCCrossedRows) ||
        (cand.photon2PosTPCCrossedRows() < photonSelections.PhotonMinTPCCrossedRows) ||
        (cand.photon1NegTPCCrossedRows() < photonSelections.PhotonMinTPCCrossedRows) ||
        (cand.photon2NegTPCCrossedRows() < photonSelections.PhotonMinTPCCrossedRows))
      return false;

    if (((cand.photon1PosTPCNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) ||
         (cand.photon1PosTPCNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)) ||
        ((cand.photon2PosTPCNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) ||
         (cand.photon2PosTPCNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)))
      return false;

    if (((cand.photon1NegTPCNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) ||
         (cand.photon1NegTPCNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)) ||
        ((cand.photon2NegTPCNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) ||
         (cand.photon2NegTPCNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)))
      return false;

    if (((cand.photon1Pt() < photonSelections.PhotonMinPt) ||
         (cand.photon1Pt() > photonSelections.PhotonMaxPt)) ||
        ((cand.photon2Pt() < photonSelections.PhotonMinPt) ||
         (cand.photon2Pt() > photonSelections.PhotonMaxPt)))
      return false;

    // Rapidity selection
    if ((cand.photon1Y() < photonSelections.PhotonMinRapidity) || (cand.photon1Y() > photonSelections.PhotonMaxRapidity))
      return false;
    if ((cand.photon1PosEta() < photonSelections.PhotonMinDauEta) || (cand.photon1PosEta() > photonSelections.PhotonMaxDauEta))
      return false;
    if ((cand.photon1NegEta() < photonSelections.PhotonMinDauEta) || (cand.photon1NegEta() > photonSelections.PhotonMaxDauEta))
      return false;

    if ((cand.photon2Y() < photonSelections.PhotonMinRapidity) || (cand.photon2Y() > photonSelections.PhotonMaxRapidity))
      return false;
    if ((cand.photon2PosEta() < photonSelections.PhotonMinDauEta) || (cand.photon2PosEta() > photonSelections.PhotonMaxDauEta))
      return false;
    if ((cand.photon2NegEta() < photonSelections.PhotonMinDauEta) || (cand.photon2NegEta() > photonSelections.PhotonMaxDauEta))
      return false;

    // Z conversion
    if ((cand.photon1Zconv() < photonSelections.PhotonMinZ) || (cand.photon1Zconv() > photonSelections.PhotonMaxZ))
      return false;

    if ((cand.photon2Zconv() < photonSelections.PhotonMinZ) || (cand.photon2Zconv() > photonSelections.PhotonMaxZ))
      return false;

    if (((cand.photon1Radius() < photonSelections.PhotonMinRadius) || (cand.photon1Radius() > photonSelections.PhotonMaxRadius)) ||
        ((cand.photon2Radius() < photonSelections.PhotonMinRadius) || (cand.photon2Radius() > photonSelections.PhotonMaxRadius)))
      return false;

    if ((cand.photon1Qt() > photonSelections.PhotonMaxQt) || cand.photon2Qt() > photonSelections.PhotonMaxQt)
      return false;

    if ((TMath::Abs(cand.photon1Alpha()) > photonSelections.PhotonMaxAlpha) || (TMath::Abs(cand.photon2Alpha()) > photonSelections.PhotonMaxAlpha))
      return false;

    if ((cand.photon1CosPA() < photonSelections.PhotonMinV0cospa) || (cand.photon2CosPA() < photonSelections.PhotonMinV0cospa))
      return false;

    if ((TMath::Abs(cand.photon1Mass()) > photonSelections.PhotonMaxMass) || (TMath::Abs(cand.photon2Mass()) > photonSelections.PhotonMaxMass))
      return false;

    // Pi0 specific selections
    float pi0Y = cand.pi0Y();
    if constexpr (requires { cand.pi0MCY(); }) { // If MC
      pi0Y = cand.pi0MCY();
    }

    // Rapidity
    if ((pi0Y < pi0Selections.Pi0MinRapidity) || (pi0Y > pi0Selections.Pi0MaxRapidity))
      return false;

    // V0Pair Radius
    if (cand.radius() > pi0Selections.Pi0MaxRadius)
      return false;

    // DCA V0Pair Daughters
    if (cand.dcadaughters() > pi0Selections.Pi0MaxDCADau)
      return false;

    return true;
  }

  // Main Pi0 QA analysis function
  template <typename TCollisions, typename TPi0s>
  void analyzeRecoePi0s(TCollisions const& collisions, TPi0s const& fullPi0s)
  {
    // Custom grouping
    std::vector<std::vector<int>> pi0grouped(collisions.size());

    for (const auto& pi0 : fullPi0s) {
      pi0grouped[pi0.straCollisionId()].push_back(pi0.globalIndex());
    }

    // Collisions loop
    for (const auto& coll : collisions) {

      // Event selection
      if (!IsEventAccepted(coll, true))
        continue;

      // Pi0s loop
      float centrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();

      for (size_t i = 0; i < pi0grouped[coll.globalIndex()].size(); i++) {
        auto pi0 = fullPi0s.rawIteratorAt(pi0grouped[coll.globalIndex()][i]);

        // Select sigma0 candidates
        if (!processPi0Candidate(pi0))
          continue;

        // If MC
        if constexpr (requires { pi0.isPi0(); }) {
          if (selRecoFromGenerator && !pi0.isProducedByGenerator())
            continue;

          if (pi0.isPi0())
            histos.fill(HIST("Pi0/h3dMass_MCAssociated"), centrality, pi0.mcpt(), pi0.pi0Mass());
        }

        // Fill histos after all selections
        histos.fill(HIST("Pi0/hMass"), pi0.pi0Mass());
        histos.fill(HIST("Pi0/hPt"), pi0.pt());
        histos.fill(HIST("Pi0/hY"), pi0.pi0Y());
        histos.fill(HIST("Pi0/h3dMass"), centrality, pi0.pt(), pi0.pi0Mass());
      }
    }
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, Sigma0s const& fullSigma0s)
  {
    analyzeRecoeSigma0s(collisions, fullSigma0s);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, MCSigma0s const& fullSigma0s)
  {
    analyzeRecoeSigma0s(collisions, fullSigma0s);
  }

  // Simulated processing in Run 3
  void processGeneratedRun3(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, soa::Join<aod::Sigma0Gens, aod::SigmaGenCollRef> const& Sigma0Gens)
  {
    analyzeGenerated(mcCollisions, collisions, Sigma0Gens);
  }

  // _____________________________________________________
  // Pi0 QA
  void processPi0RealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, soa::Join<aod::Pi0Cores, aod::Pi0CollRef> const& fullPi0s)
  {
    analyzeRecoePi0s(collisions, fullPi0s);
  }

  void processPi0MonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, soa::Join<aod::Pi0Cores, aod::Pi0CoresMC, aod::Pi0CollRef> const& fullPi0s)
  {
    analyzeRecoePi0s(collisions, fullPi0s);
  }

  void processPi0GeneratedRun3(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, soa::Join<aod::Pi0Gens, aod::Pi0GenCollRef> const& Pi0Gens)
  {
    analyzeGenerated(mcCollisions, collisions, Pi0Gens);
  }

  // _____________________________________________________
  PROCESS_SWITCH(sigmaanalysis, processRealData, "Do real data analysis", true);
  PROCESS_SWITCH(sigmaanalysis, processMonteCarlo, "Do Monte-Carlo-based analysis", false);
  PROCESS_SWITCH(sigmaanalysis, processGeneratedRun3, "process MC generated Run 3", false);
  PROCESS_SWITCH(sigmaanalysis, processPi0RealData, "Do real data analysis for pi0 QA", false);
  PROCESS_SWITCH(sigmaanalysis, processPi0MonteCarlo, "Do Monte-Carlo-based analysis for pi0 QA", false);
  PROCESS_SWITCH(sigmaanalysis, processPi0GeneratedRun3, "process MC generated Run 3 for pi0 QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigmaanalysis>(cfgc)};
}
