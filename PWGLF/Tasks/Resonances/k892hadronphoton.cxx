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

/// \file k892hadronphoton.cxx
/// \brief This is a task that reads kstar tables (from sigma0builder) to perform analysis.
/// \author Gianni Shigeru Setoue Liveraro
/// \author Oussama Benchikhi

#include "PWGLF/DataModel/LFSigmaTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using KStars = soa::Join<aod::KStarCores, aod::KStarPhotonExtras, aod::KShortExtras, aod::KStarCollRef>;
using MCKStars = soa::Join<aod::KStarCores, aod::KStarPhotonExtras, aod::KShortExtras, aod::KStarMCCores, aod::KStarCollRef>;

static const std::vector<std::string> photonSels = {"NoSel", "V0Type", "DCADauToPV",
                                                    "DCADau", "DauTPCCR", "TPCNSigmaEl", "V0pT",
                                                    "Y", "V0Radius", "RZCut", "Armenteros", "CosPA",
                                                    "PsiPair", "Phi", "Mass"};

static const std::vector<std::string> kshortSels = {"NoSel", "V0Radius", "DCADau", "Armenteros",
                                                    "CosPA", "Y", "TPCCR", "DauITSCls", "Lifetime",
                                                    "TPCTOFPID", "DCADauToPV", "Mass"};

static const std::vector<std::string> DirList = {"BeforeSel", "AfterSel"};

struct k892hadronphoton {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  //__________________________________________________
  // For manual sliceBy
  // SliceCache cache;
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraEvSelExtras, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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

  // generated
  Configurable<bool> mcKeepOnlyFromGenerator{"mcKeepOnlyFromGenerator", true, "if true, consider only particles from generator to calculate efficiency."};

  // QA
  Configurable<bool> fillBkgQAhistos{"fillBkgQAhistos", false, "if true, fill MC QA histograms for Bkg study. Only works with MC."};
  Configurable<bool> fillResoQAhistos{"fillResoQAhistos", false, "if true, fill MC QA histograms for pT resolution study. Only works with MC."};

  // Analysis strategy:
  Configurable<bool> doMCAssociation{"doMCAssociation", false, "Flag to process only signal candidates. Use only with processMonteCarlo!"};
  Configurable<bool> selRecoFromGenerator{"selRecoFromGenerator", false, "Flag to process only signal candidates from generator"};
  Configurable<bool> doPhotonKShortSelQA{"doPhotonKShortSelQA", false, "Flag to fill photon and kshort QA histos!"};

  //// K0Short criteria:
  struct : ConfigurableGroup {
    Configurable<float> kshortMLThreshold{"kshortMLThreshold", 0.1, "Decision Threshold value to select kshorts"};
    Configurable<float> kshortMinDCANegToPv{"kshortMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> kshortMinDCAPosToPv{"kshortMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> kshortMaxDCAV0Dau{"kshortMaxDCAV0Dau", 2.5, "Max DCA V0 Daughters (cm)"};
    Configurable<float> kshortMinv0radius{"kshortMinv0radius", 0.0, "Min V0 radius (cm)"};
    Configurable<float> kshortMaxv0radius{"kshortMaxv0radius", 40, "Max V0 radius (cm)"};
    Configurable<float> kshortMinQt{"kshortMinQt", 0.1, "Min kshort qt value (AP plot) (GeV/c)"};
    Configurable<float> kshortMaxQt{"kshortMaxQt", 2.5, "Max kshort qt value (AP plot) (GeV/c)"};
    Configurable<float> kshortMinAlpha{"kshortMinAlpha", 0.0, "Min kshort alpha absolute value (AP plot)"};
    Configurable<float> kshortMaxAlpha{"kshortMaxAlpha", 1.0, "Max kshort alpha absolute value (AP plot)"};
    Configurable<float> kshortMinv0cospa{"kshortMinv0cospa", 0.95, "Min V0 CosPA"};
    Configurable<float> kshortMaxLifeTime{"kshortMaxLifeTime", 30, "Max lifetime"};
    Configurable<float> kshortWindow{"kshortWindow", 0.015, "Mass window around expected (in GeV/c2)"};
    Configurable<float> kshortMaxRap{"kshortMaxRap", 0.8, "Max kshort rapidity"};
    Configurable<float> kshortMaxDauEta{"kshortMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<bool> fselKShortTPCPID{"fselKShortTPCPID", true, "Flag to select kshort-like candidates using TPC NSigma."};
    Configurable<bool> fselKShortTOFPID{"fselKShortTOFPID", false, "Flag to select kshort-like candidates using TOF NSigma."};
    Configurable<float> kshortMaxTPCNSigmas{"kshortMaxTPCNSigmas", 1e+9, "Max TPC NSigmas for daughters"};
    // Configurable<float> kshortMaxTOFNSigmas{"kshortMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
    Configurable<float> kshortPiMaxTOFNSigmas{"kshortPiMaxTOFNSigmas", 1e+9, "Max TOF NSigmas for daughters"};
    Configurable<int> kshortMinTPCCrossedRows{"kshortMinTPCCrossedRows", 50, "Min daughter TPC Crossed Rows"};
    Configurable<int> kshortMinITSclusters{"kshortMinITSclusters", 1, "minimum ITS clusters"};
    Configurable<bool> kshortRejectPosITSafterburner{"kshortRejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
    Configurable<bool> kshortRejectNegITSafterburner{"kshortRejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};
  } kshortSelections;

  //// Photon criteria:
  struct : ConfigurableGroup {
    Configurable<float> gammaMLThreshold{"gammaMLThreshold", 0.1, "Decision Threshold value to select gammas"};
    Configurable<int> photonv0TypeSel{"photonv0TypeSel", 7, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> photonMinDCADauToPv{"photonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
    Configurable<float> photonMaxDCAV0Dau{"photonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
    Configurable<int> photonMinTPCCrossedRows{"photonMinTPCCrossedRows", 30, "Min daughter TPC Crossed Rows"};
    Configurable<float> photonMinTPCNSigmas{"photonMinTPCNSigmas", -7, "Min TPC NSigmas for daughters"};
    Configurable<float> photonMaxTPCNSigmas{"photonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
    Configurable<float> photonMinPt{"photonMinPt", 0.0, "Min photon pT (GeV/c)"};
    Configurable<float> photonMaxPt{"photonMaxPt", 50.0, "Max photon pT (GeV/c)"};
    Configurable<float> photonMaxRap{"photonMaxRap", 0.5, "Max photon rapidity"};
    Configurable<float> photonMinRadius{"photonMinRadius", 3.0, "Min photon conversion radius (cm)"};
    Configurable<float> photonMaxRadius{"photonMaxRadius", 115, "Max photon conversion radius (cm)"};
    Configurable<float> photonMaxZ{"photonMaxZ", 240, "Max photon conversion point z value (cm)"};
    Configurable<float> photonMaxQt{"photonMaxQt", 0.05, "Max photon qt value (AP plot) (GeV/c)"};
    Configurable<float> photonMaxAlpha{"photonMaxAlpha", 0.95, "Max photon alpha absolute value (AP plot)"};
    Configurable<float> photonMinV0cospa{"photonMinV0cospa", 0.80, "Min V0 CosPA"};
    Configurable<float> photonMaxMass{"photonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};
    Configurable<float> photonPsiPairMax{"photonPsiPairMax", 1e+9, "maximum psi angle of the track pair"};
    Configurable<float> photonMaxDauEta{"photonMaxDauEta", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<float> photonLineCutZ0{"photonLineCutZ0", 7.0, "The offset for the linecute used in the Z vs R plot"};
    Configurable<float> photonPhiMin1{"photonPhiMin1", -1, "Phi min value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> photonPhiMax1{"photonPhiMax1", -1, "Phi max value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> photonPhiMin2{"photonPhiMin2", -1, "Phi max value to reject photons, region 2 (leave negative if no selection desired)"};
    Configurable<float> photonPhiMax2{"photonPhiMax2", -1, "Phi min value to reject photons, region 2 (leave negative if no selection desired)"};
  } photonSelections;

  struct : ConfigurableGroup {
    Configurable<float> kstarMaxRap{"kstarMaxRap", 0.5, "Max kstar rapidity"};
    Configurable<float> kstarMaxRadius{"kstarMaxRadius", 200, "Max kstar decay radius"};
    Configurable<float> kstarMaxDCADau{"kstarMaxDCADau", 50, "Max kstar DCA between daughters"};
    Configurable<float> kstarMaxOPAngle{"kstarMaxOPAngle", 7, "Max kstar OP Angle between daughters"};
  } kstarSelections;

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

  // invariant mass
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, -0.1f, 0.5f}, "M_{#Gamma}"};
  ConfigurableAxis axisKShortMass{"axisKShortMass", {200, 0.3f, 0.6f}, "M_{#K_{s}^{0}}"};
  ConfigurableAxis axisKStarMass{"axisKStarMass", {200, 0.08f, 1.5f}, "M_{#K^{*}}"};

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
  ConfigurableAxis mlProb{"mlOutput", {100, 0.0f, 1.0f}, ""};

  void init(InitContext const&)
  {
    LOGF(info, "Initializing now: cross-checking correctness...");
    if (doprocessRealData + doprocessMonteCarlo > 1) {
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

        histos.add(histodir + "/KShort/hTrackCode", "hTrackCode", kTH1D, {{11, 0.5f, 11.5f}});
        histos.add(histodir + "/KShort/hRadius", "hRadius", kTH1D, {axisV0Radius});
        histos.add(histodir + "/KShort/hDCADau", "hDCADau", kTH1D, {axisDCAdau});
        histos.add(histodir + "/KShort/hCosPA", "hCosPA", kTH1D, {axisCosPA});
        histos.add(histodir + "/KShort/hY", "hY", kTH1D, {axisRapidity});
        histos.add(histodir + "/KShort/hPosEta", "hPosEta", kTH1D, {axisRapidity});
        histos.add(histodir + "/KShort/hNegEta", "hNegEta", kTH1D, {axisRapidity});
        histos.add(histodir + "/KShort/hPosTPCCR", "hPosTPCCR", kTH1D, {axisTPCrows});
        histos.add(histodir + "/KShort/hNegTPCCR", "hNegTPCCR", kTH1D, {axisTPCrows});
        histos.add(histodir + "/KShort/hPosITSCls", "hPosITSCls", kTH1D, {axisNCls});
        histos.add(histodir + "/KShort/hNegITSCls", "hNegITSCls", kTH1D, {axisNCls});
        histos.add(histodir + "/KShort/hPosChi2PerNc", "hPosChi2PerNc", kTH1D, {axisChi2PerNcl});
        histos.add(histodir + "/KShort/hNegChi2PerNc", "hNegChi2PerNc", kTH1D, {axisChi2PerNcl});
        histos.add(histodir + "/KShort/hLifeTime", "hLifeTime", kTH1D, {axisLifetime});
        histos.add(histodir + "/KShort/h2dTPCvsTOFNSigma_KShortPi", "h2dTPCvsTOFNSigma_KShortPi", kTH2D, {axisTPCNSigma, axisTOFNSigma});
        histos.add(histodir + "/KShort/hKShortDCANegToPV", "hKShortDCANegToPV", kTH1D, {axisDCAtoPV});
        histos.add(histodir + "/KShort/hKShortDCAPosToPV", "hKShortDCAPosToPV", kTH1D, {axisDCAtoPV});
        histos.add(histodir + "/KShort/hKShortpT", "hKShortpT", kTH1D, {axisPt});
        histos.add(histodir + "/KShort/hKShortMass", "hKShortMass", kTH1D, {axisKShortMass});
        histos.add(histodir + "/KShort/h3dKShortMass", "h3dKShortMass", kTH3D, {axisCentrality, axisPt, axisKShortMass});

        histos.add(histodir + "/h2dArmenteros", "h2dArmenteros", kTH2D, {axisAPAlpha, axisAPQt});

        histos.add(histodir + "/KStar/hMass", "hMass", kTH1D, {axisKStarMass});
        histos.add(histodir + "/KStar/hPt", "hPt", kTH1D, {axisPt});
        histos.add(histodir + "/KStar/hY", "hY", kTH1D, {axisRapidity});
        histos.add(histodir + "/KStar/hRadius", "hRadius", kTH1D, {axisV0PairRadius});
        histos.add(histodir + "/KStar/h2dRadiusVspT", "h2dRadiusVspT", kTH2D, {axisV0PairRadius, axisPt});
        histos.add(histodir + "/KStar/hDCAPairDau", "hDCAPairDau", kTH1D, {axisDCAdau});
        histos.add(histodir + "/KStar/hDCAPairDauVsPt", "hDCAPairDauVsPt", kTH2D, {axisDCAdau, axisPt});

        histos.add(histodir + "/KStar/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisKStarMass});
        histos.add(histodir + "/KStar/h3dOPAngleVsMass", "h3dOPAngleVsMass", kTH3D, {{140, 0.0f, +7.0f}, axisPt, axisKStarMass});
        histos.add(histodir + "/KStar/h2dOPAngleVsPt", "h2dOPAngleVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
        histos.add(histodir + "/KStar/h2dOPAnglePhotonVsPt", "h2dOPAnglePhotonVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
        histos.add(histodir + "/KStar/h2dOPAngleKShortVsPt", "h2dOPAngleKShortVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});

        if (doprocessMonteCarlo) {

          histos.add(histodir + "/MC/Photon/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1D, {{2, 0.0f, 2.0f}});
          histos.add(histodir + "/MC/Photon/hPt", "hPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/Photon/hMCPt", "hMCPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/Photon/hPosTPCNSigmaEl", "hPosTPCNSigmaEl", kTH1D, {axisTPCNSigma});
          histos.add(histodir + "/MC/Photon/hNegTPCNSigmaEl", "hNegTPCNSigmaEl", kTH1D, {axisTPCNSigma});
          histos.add(histodir + "/MC/Photon/h2dPAVsPt", "h2dPAVsPt", kTH2D, {axisPA, axisPt});
          histos.add(histodir + "/MC/Photon/hPt_BadCollAssig", "hPt_BadCollAssig", kTH1D, {axisPt});
          histos.add(histodir + "/MC/Photon/h2dPAVsPt_BadCollAssig", "h2dPAVsPt_BadCollAssig", kTH2D, {axisPA, axisPt});

          histos.add(histodir + "/MC/KShort/hV0ToCollAssoc", "hV0ToCollAssoc", kTH1D, {{2, 0.0f, 2.0f}});
          histos.add(histodir + "/MC/KShort/hPt", "hPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/KShort/hMass", "hMass", kTH1D, {axisKShortMass});
          histos.add(histodir + "/MC/KShort/hMCPt", "hMCPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/KShort/h3dTPCvsTOFNSigma_Pi", "h3dTPCvsTOFNSigma_Pi", kTH3D, {axisTPCNSigma, axisTOFNSigma, axisPt});

          histos.add(histodir + "/MC/h2dArmenteros", "h2dArmenteros", kTH2D, {axisAPAlpha, axisAPQt});

          histos.add(histodir + "/MC/KStar/hPt", "hPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/KStar/hMCPt", "hMCPt", kTH1D, {axisPt});
          histos.add(histodir + "/MC/KStar/hMass", "hMass", kTH1D, {axisKStarMass});
          histos.add(histodir + "/MC/KStar/hMCProcess", "hMCProcess", kTH1D, {{50, -0.5f, 49.5f}});
          histos.add(histodir + "/MC/KStar/hGenRadius", "hGenRadius", kTH1D, {axisV0PairRadius});
          histos.add(histodir + "/MC/KStar/h2dMCPtVsKShortMCPt", "h2dMCPtVsKShortMCPt", kTH2D, {axisPt, axisPt});
          histos.add(histodir + "/MC/KStar/h2dMCPtVsPhotonMCPt", "h2dMCPtVsPhotonMCPt", kTH2D, {axisPt, axisPt});
          histos.add(histodir + "/MC/KStar/h2dMCProcessVsGenRadius", "h2dMCProcessVsGenRadius", kTH2D, {{50, -0.5f, 49.5f}, axisV0PairRadius});
          histos.add(histodir + "/MC/KStar/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisKStarMass});
          histos.add(histodir + "/MC/KStar/h3dMCProcess", "h3dMCProcess", kTH3D, {{50, -0.5f, 49.5f}, axisPt, axisKStarMass});
          histos.add(histodir + "/MC/KStar/h2dOPAngleVsPt", "h2dOPAngleVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
          histos.add(histodir + "/MC/KStar/h2dOPAngleKShortVsPt", "h2dOPAngleKShortVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
          histos.add(histodir + "/MC/KStar/h2dOPAnglePhotonVsPt", "h2dOPAnglePhotonVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
          histos.add(histodir + "/MC/KStar/h2dRadiusVspT", "h2dRadiusVspT", kTH2D, {axisV0PairRadius, axisPt});
          histos.add(histodir + "/MC/KStar/hDCAPairDauVsPt", "hDCAPairDauVsPt", kTH2D, {axisDCAdau, axisPt});

          // K*(892)0 -> pi K candidates
          histos.add(histodir + "/PionKaon/Photon/hTrackCode", "hTrackCode", kTH1D, {{11, 0.5f, 11.5f}});
          histos.add(histodir + "/PionKaon/Photon/hV0Type", "hV0Type", kTH1D, {{8, 0.5f, 8.5f}});
          histos.add(histodir + "/PionKaon/Photon/hDCANegToPV", "hDCANegToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/PionKaon/Photon/hDCAPosToPV", "hDCAPosToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/PionKaon/Photon/hDCADau", "hDCADau", kTH1D, {axisDCAdau});
          histos.add(histodir + "/PionKaon/Photon/hPosTPCCR", "hPosTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/PionKaon/Photon/hNegTPCCR", "hNegTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/PionKaon/Photon/hPosTPCNSigmaEl", "hPosTPCNSigmaEl", kTH1D, {axisTPCNSigma});
          histos.add(histodir + "/PionKaon/Photon/hNegTPCNSigmaEl", "hNegTPCNSigmaEl", kTH1D, {axisTPCNSigma});
          histos.add(histodir + "/PionKaon/Photon/hpT", "hpT", kTH1D, {axisPt});
          histos.add(histodir + "/PionKaon/Photon/hY", "hY", kTH1D, {axisRapidity});
          histos.add(histodir + "/PionKaon/Photon/hPosEta", "hPosEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/PionKaon/Photon/hNegEta", "hNegEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/PionKaon/Photon/hRadius", "hRadius", kTH1D, {axisV0Radius});
          histos.add(histodir + "/PionKaon/Photon/hZ", "hZ", kTH1D, {axisZ});
          histos.add(histodir + "/PionKaon/Photon/h2dRZCut", "h2dRZCut", kTH2D, {axisZ, axisV0Radius});
          histos.add(histodir + "/PionKaon/Photon/h2dRZPlane", "h2dRZPlane", kTH2D, {axisZ, axisV0Radius});
          histos.add(histodir + "/PionKaon/Photon/hCosPA", "hCosPA", kTH1D, {axisCosPA});
          histos.add(histodir + "/PionKaon/Photon/hPsiPair", "hPsiPair", kTH1D, {axisPsiPair});
          histos.add(histodir + "/PionKaon/Photon/hPhi", "hPhi", kTH1D, {axisPhi});
          histos.add(histodir + "/PionKaon/Photon/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisPhotonMass});
          histos.add(histodir + "/PionKaon/Photon/hMass", "hMass", kTH1D, {axisPhotonMass});

          histos.add(histodir + "/PionKaon/h2dArmenteros", "h2dArmenteros", kTH2D, {axisAPAlpha, axisAPQt});

          histos.add(histodir + "/PionKaon/KShort/hTrackCode", "hTrackCode", kTH1D, {{11, 0.5f, 11.5f}});
          histos.add(histodir + "/PionKaon/KShort/hRadius", "hRadius", kTH1D, {axisV0Radius});
          histos.add(histodir + "/PionKaon/KShort/hDCADau", "hDCADau", kTH1D, {axisDCAdau});
          histos.add(histodir + "/PionKaon/KShort/hCosPA", "hCosPA", kTH1D, {axisCosPA});
          histos.add(histodir + "/PionKaon/KShort/hY", "hY", kTH1D, {axisRapidity});
          histos.add(histodir + "/PionKaon/KShort/hPosEta", "hPosEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/PionKaon/KShort/hNegEta", "hNegEta", kTH1D, {axisRapidity});
          histos.add(histodir + "/PionKaon/KShort/hPosTPCCR", "hPosTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/PionKaon/KShort/hNegTPCCR", "hNegTPCCR", kTH1D, {axisTPCrows});
          histos.add(histodir + "/PionKaon/KShort/hPosITSCls", "hPosITSCls", kTH1D, {axisNCls});
          histos.add(histodir + "/PionKaon/KShort/hNegITSCls", "hNegITSCls", kTH1D, {axisNCls});
          histos.add(histodir + "/PionKaon/KShort/hPosChi2PerNc", "hPosChi2PerNc", kTH1D, {axisChi2PerNcl});
          histos.add(histodir + "/PionKaon/KShort/hNegChi2PerNc", "hNegChi2PerNc", kTH1D, {axisChi2PerNcl});
          histos.add(histodir + "/PionKaon/KShort/hLifeTime", "hLifeTime", kTH1D, {axisLifetime});
          histos.add(histodir + "/PionKaon/KShort/h2dTPCvsTOFNSigma_KShortPi", "h2dTPCvsTOFNSigma_KShortPi", kTH2D, {axisTPCNSigma, axisTOFNSigma});
          histos.add(histodir + "/PionKaon/KShort/hKShortDCANegToPV", "hKShortDCANegToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/PionKaon/KShort/hKShortDCAPosToPV", "hKShortDCAPosToPV", kTH1D, {axisDCAtoPV});
          histos.add(histodir + "/PionKaon/KShort/hKShortpT", "hKShortpT", kTH1D, {axisPt});
          histos.add(histodir + "/PionKaon/KShort/hKShortMass", "hKShortMass", kTH1D, {axisKShortMass});
          histos.add(histodir + "/PionKaon/KShort/h3dKShortMass", "h3dKShortMass", kTH3D, {axisCentrality, axisPt, axisKShortMass});

          histos.add(histodir + "/PionKaon/KStar/hMass", "hMass", kTH1D, {axisKStarMass});
          histos.add(histodir + "/PionKaon/KStar/hPt", "hPt", kTH1D, {axisPt});
          histos.add(histodir + "/PionKaon/KStar/hY", "hY", kTH1D, {axisRapidity});
          histos.add(histodir + "/PionKaon/KStar/hRadius", "hRadius", kTH1D, {axisV0PairRadius});
          histos.add(histodir + "/PionKaon/KStar/h2dRadiusVspT", "h2dRadiusVspT", kTH2D, {axisV0PairRadius, axisPt});
          histos.add(histodir + "/PionKaon/KStar/hDCAPairDau", "hDCAPairDau", kTH1D, {axisDCAdau});
          histos.add(histodir + "/PionKaon/KStar/hDCAPairDauVsPt", "hDCAPairDauVsPt", kTH2D, {axisDCAdau, axisPt});
          histos.add(histodir + "/PionKaon/KStar/h3dMass", "h3dMass", kTH3D, {axisCentrality, axisPt, axisKStarMass});
          histos.add(histodir + "/PionKaon/KStar/h3dOPAngleVsMass", "h3dOPAngleVsMass", kTH3D, {{140, 0.0f, +7.0f}, axisPt, axisKStarMass});
          histos.add(histodir + "/PionKaon/KStar/h2dOPAngleVsPt", "h2dOPAngleVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
          histos.add(histodir + "/PionKaon/KStar/h2dOPAngleKShortVsPt", "h2dOPAngleKShortVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
          histos.add(histodir + "/PionKaon/KStar/h2dOPAnglePhotonVsPt", "h2dOPAnglePhotonVsPt", kTH2D, {{140, 0.0f, +7.0f}, axisPt});
          // 1/pT Resolution:

          if (fillResoQAhistos) {
            histos.add(histodir + "/MC/Reso/h2dKShortPtResolution", "h2dKShortPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h2dGammaPtResolution", "h2dGammaPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h2dKStarPtResolution", "h2dKStarPtResolution", kTH2D, {axisInvPt, axisDeltaPt});
            histos.add(histodir + "/MC/Reso/h2dKStarRadiusResolution", "h2dKStarRadiusResolution", kTH2D, {axisPt, axisDeltaPt});
          }

          // For background decomposition study
          if (fillBkgQAhistos) {
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_All", "h2dPtVsMassKStar_All", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_TrueDaughters", "h2dPtVsMassKStar_TrueDaughters", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_PhotonOmega", "h2dPtVsMassKStar_PhotonOmega", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_PhotonRho", "h2dPtVsMassKStar_PhotonRho", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_PhotonEta", "h2dPtVsMassKStar_PhotonEta", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_KShortKCharged", "h2dPtVsMassKStar_KShortKCharged", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_KStarPionKaon", "h2dPtVsMassKStar_KStarPionKaon", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_PionKaon", "h2dPtVsMassKStar_PionKaon", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_TrueGammaFakeKShort", "h2dPtVsMassKStar_TrueGammaFakeKShort", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_FakeGammaTrueKShort", "h2dPtVsMassKStar_FakeGammaTrueKShort", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dPtVsMassKStar_FakeDaughters", "h2dPtVsMassKStar_FakeDaughters", kTH2D, {axisPt, axisKStarMass});
            histos.add(histodir + "/MC/BkgStudy/h2dTrueDaughtersMatrix", "h2dTrueDaughtersMatrix", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
            histos.add(histodir + "/MC/BkgStudy/h2dPionKaonMother", "h2dPionKaonMother", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
            histos.add(histodir + "/MC/BkgStudy/h2dTrueDaughtersMatrixGrandMother", "h2dTrueDaughtersMatrixGrandMother", kTHnSparseD, {{10001, -5000.5f, +5000.5f}, {10001, -5000.5f, +5000.5f}});
          }
        }
      }

      // Selections
      histos.add("Selection/Photon/hCandidateSel", "hCandidateSel", kTH1D, {axisCandSel});
      histos.add("Selection/KShort/hCandidateSel", "hCandidateSel", kTH1D, {axisCandSel});

      for (size_t i = 0; i < photonSels.size(); ++i) {
        const auto& sel = photonSels[i];

        histos.add(Form("Selection/Photon/h2d%s", sel.c_str()), ("h2d" + sel).c_str(), kTH2D, {axisPt, axisPhotonMass});
        histos.get<TH1>(HIST("Selection/Photon/hCandidateSel"))->GetXaxis()->SetBinLabel(i + 1, sel.c_str());
        histos.add(Form("Selection/KStar/h2dPhoton%s", sel.c_str()), ("h2dPhoton" + sel).c_str(), kTH2D, {axisPt, axisKStarMass});
      }

      for (size_t i = 0; i < kshortSels.size(); ++i) {
        const auto& sel = kshortSels[i];

        histos.add(Form("Selection/KShort/h2d%s", sel.c_str()), ("h2d" + sel).c_str(), kTH2D, {axisPt, axisKShortMass});
        histos.get<TH1>(HIST("Selection/KShort/hCandidateSel"))->GetXaxis()->SetBinLabel(i + 1, sel.c_str());
        histos.add(Form("Selection/KStar/h2dKShort%s", sel.c_str()), ("h2dKShort" + sel).c_str(), kTH2D, {axisPt, axisKStarMass});
      }
    }

    if (doprocessGeneratedRun3) {

      histos.add("Gen/hGenEvents", "hGenEvents", kTH2D, {{axisNch}, {2, -0.5f, +1.5f}});
      histos.get<TH2>(HIST("Gen/hGenEvents"))->GetYaxis()->SetBinLabel(1, "All gen. events");
      histos.get<TH2>(HIST("Gen/hGenEvents"))->GetYaxis()->SetBinLabel(2, "Gen. with at least 1 rec. events");

      histos.add("Gen/hGenEventCentrality", "hGenEventCentrality", kTH1D, {axisCentrality});
      histos.add("Gen/hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", kTH2D, {axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("Gen/hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", kTH2D, {axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("Gen/hCentralityVsMultMC", "hCentralityVsMultMC", kTH2D, {axisCentrality, axisNch});

      histos.add("Gen/hEventPVzMC", "hEventPVzMC", kTH1D, {{100, -20.0f, +20.0f}});
      histos.add("Gen/hCentralityVsPVzMC", "hCentralityVsPVzMC", kTH2D, {{101, 0.0f, 101.0f}, {100, -20.0f, +20.0f}});

      // KStar specific
      histos.add("Gen/h2dGenKStar", "h2dGenKStar", kTH2D, {axisCentrality, axisPt});
      histos.add("Gen/h2dGenKStarVsMultMC_RecoedEvt", "h2dGenKStarVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
      histos.add("Gen/h2dGenKStarVsMultMC", "h2dGenKStarVsMultMC", kTH2D, {axisNch, axisPt});
    }
  }

  // ______________________________________________________
  // Check whether the collision passes our collision selections
  // Should work with collisions, mccollisions, stracollisions and stramccollisions tables!
  template <typename TCollision>
  bool isEventAccepted(TCollision const& collision, bool fillHists)
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

    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 19 /* Above max IR */);
      // Fill centrality histogram after event selection
      histos.fill(HIST("hEventCentrality"), centrality);
    }

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
          if (!isEventAccepted(collision, false)) {
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

        if (!isEventAccepted(collision, false)) {
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

    for (const auto& genParticle : genParticles) {
      float centrality = 100.5f;

      // Has MC collision
      if (!genParticle.has_straMCCollision())
        continue;

      // Selection on the source (generator/transport)
      if (!genParticle.producedByGenerator() && mcKeepOnlyFromGenerator)
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
      // Generated KStar processing
      if constexpr (requires { genParticle.kstarMCPt(); }) {

        float ptmc = genParticle.kstarMCPt();

        if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
          auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
          centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

          if (genParticle.isKStar())
            histos.fill(HIST("Gen/h2dGenKStarVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }

        if (genParticle.isKStar()) {
          histos.fill(HIST("Gen/h2dGenKStar"), centrality, ptmc);
          histos.fill(HIST("Gen/h2dGenKStarVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
      }
    }
  }

  //__________________________________________
  template <bool isGamma, typename TKStarObject>
  int retrieveV0TrackCode(TKStarObject const& kstar)
  {

    int trkCode = 10; // 1: TPC-only, 2: TPC+Something, 3: ITS-Only, 4: ITS+TPC + Something, 10: anything else

    if (isGamma) {
      if (kstar.photonPosTrackCode() == 1 && kstar.photonNegTrackCode() == 1)
        trkCode = 1;
      if ((kstar.photonPosTrackCode() != 1 && kstar.photonNegTrackCode() == 1) || (kstar.photonPosTrackCode() == 1 && kstar.photonNegTrackCode() != 1))
        trkCode = 2;
      if (kstar.photonPosTrackCode() == 3 && kstar.photonNegTrackCode() == 3)
        trkCode = 3;
      if (kstar.photonPosTrackCode() == 2 || kstar.photonNegTrackCode() == 2)
        trkCode = 4;
    } else {
      if (kstar.kshortPosTrackCode() == 1 && kstar.kshortNegTrackCode() == 1)
        trkCode = 1;
      if ((kstar.kshortPosTrackCode() != 1 && kstar.kshortNegTrackCode() == 1) || (kstar.kshortPosTrackCode() == 1 && kstar.kshortNegTrackCode() != 1))
        trkCode = 2;
      if (kstar.kshortPosTrackCode() == 3 && kstar.kshortNegTrackCode() == 3)
        trkCode = 3;
      if (kstar.kshortPosTrackCode() == 2 || kstar.kshortNegTrackCode() == 2)
        trkCode = 4;
    }

    return trkCode;
  }

  template <int mode, typename TKStarObject>
  void getResolution(TKStarObject const& kstar)
  {
    // Check whether it is before or after selections
    static constexpr std::string_view MainDir[] = {"BeforeSel", "AfterSel"};

    //_______________________________________
    // Gamma MC association
    if (std::abs(kstar.photonPDGCode()) == PDG_t::kGamma) {
      if (kstar.photonmcpt() > 0) {
        // histos.fill(HIST(MainDir[mode]) + HIST("/MC/Reso/h2dGammaPtResolution"), 1.f / kstar.photonmcpt(), kstar.photonPt() - kstar.photonmcpt());
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/Reso/h2dGammaPtResolution"), 1.f / kstar.photonmcpt(), (kstar.photonPt() - kstar.photonmcpt()) / kstar.photonmcpt());
      }
    }

    //_______________________________________
    // KShort MC association
    if (std::abs(kstar.kshortPDGCode()) == PDG_t::kK0Short) {
      if (kstar.kshortmcpt() > 0) {
        // histos.fill(HIST(MainDir[mode]) + HIST("/MC/Reso/h2dKShortPtResolution"), 1.f / kstar.kshortmcpt(), kstar.kshortPt() - kstar.kshortmcpt());
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/Reso/h2dKShortPtResolution"), 1.f / kstar.kshortmcpt(), (kstar.kshortPt() - kstar.kshortmcpt()) / kstar.kshortmcpt());
      }
    }

    //_______________________________________
    // KStar MC association
    if (kstar.isKStar()) {
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/Reso/h2dKStarRadiusResolution"), kstar.mcpt(), kstar.radius() - kstar.mcradius()); // pT resolution
      if (kstar.mcpt() > 0)
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/Reso/h2dKStarPtResolution"), 1.f / kstar.mcpt(), (1.f / kstar.pt() - 1.f / kstar.mcpt()) / (1.f / kstar.mcpt())); // pT resolution
    }
  }

  // To save histograms for background analysis
  template <int mode, typename TKStarObject>
  void runBkgAnalysis(TKStarObject const& kstar)
  {
    // Check whether it is before or after selections
    static constexpr std::string_view MainDir[] = {"BeforeSel", "AfterSel"};

    bool fIsKStar = kstar.isKStar();
    int photonPDGCode = kstar.photonPDGCode();
    int photonPDGCodeMother = kstar.photonPDGCodeMother();
    int photonPDGCodeGrandMother = kstar.photonPDGCodeGrandMother();
    int photonGlobalIndexGrandMother = kstar.photonGlobalIndexGrandMother();
    int kshortPDGCode = kstar.kshortPDGCode();
    int kshortPDGCodeMother = kstar.kshortPDGCodeMother();
    int kshortPDGCodeGrandMother = kstar.kshortPDGCodeGrandMother();
    int kshortGlobalIndexGrandMother = kstar.kshortGlobalIndexGrandMother();
    float kstarpT = kstar.pt();
    float kstarMass = kstar.kstarMass();

    histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_All"), kstarpT, kstarMass);
    // coming from the same kstar (coupled), but via a different decay channel (K* -> K0 pi0 -> K0s gamma)

    //_______________________________________
    // Real Gamma x Real KShort - but not direct decays from the same kstar!
    if ((!fIsKStar) && (std::abs(photonPDGCode) == PDG_t::kGamma) && (std::abs(kshortPDGCode) == PDG_t::kK0Short)) {
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_TrueDaughters"), kstarpT, kstarMass);
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dTrueDaughtersMatrix"), kshortPDGCodeMother, photonPDGCodeMother);
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dTrueDaughtersMatrixGrandMother"), kshortPDGCodeGrandMother, photonPDGCodeGrandMother);

      // coming from the same kstar (coupled), but via a different decay channel (K* -> K0 pi0 -> K0s gamma)
      if (std::abs(kshortPDGCodeGrandMother) == o2::constants::physics::Pdg::kK0Star892 && std::abs(photonPDGCodeGrandMother) == o2::constants::physics::Pdg::kK0Star892 && (photonGlobalIndexGrandMother == kshortGlobalIndexGrandMother)) // K*(892)0
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_KStarPionKaon"), kstarpT, kstarMass);

      if (photonGlobalIndexGrandMother == kshortGlobalIndexGrandMother) {
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_PionKaon"), kstarpT, kstarMass);
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPionKaonMother"), kshortPDGCodeGrandMother, photonPDGCodeGrandMother);
      }

      // Break down of different photon / kshort sources
      if (std::abs(photonPDGCodeGrandMother) == o2::constants::physics::Pdg::kOmega) // omega(782)
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_PhotonOmega"), kstarpT, kstarMass);

      if (std::abs(photonPDGCodeGrandMother) == PDG_t::kRho770Plus) // rho+-(770)
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_PhotonRho"), kstarpT, kstarMass);

      if (std::abs(photonPDGCodeMother) == o2::constants::physics::Pdg::kEta) // eta
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_PhotonEta"), kstarpT, kstarMass);

      if (std::abs(kshortPDGCodeGrandMother) == o2::constants::physics::Pdg::kKPlusStar892) // K*(892)+-
        histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_KShortKCharged"), kstarpT, kstarMass);
    }

    //_______________________________________
    // Real Gamma x fake KShort
    if ((std::abs(photonPDGCode) == PDG_t::kGamma) && (std::abs(kshortPDGCode) != PDG_t::kK0Short))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_TrueGammaFakeKShort"), kstarpT, kstarMass);

    //_______________________________________
    // Fake Gamma x Real KShort
    if ((std::abs(photonPDGCode) != PDG_t::kGamma) && ((std::abs(kshortPDGCode) == PDG_t::kK0Short)))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_FakeGammaTrueKShort"), kstarpT, kstarMass);

    //_______________________________________
    // Fake Gamma x Fake KShort
    if ((std::abs(photonPDGCode) != PDG_t::kGamma) && (std::abs(kshortPDGCode) != PDG_t::kK0Short))
      histos.fill(HIST(MainDir[mode]) + HIST("/MC/BkgStudy/h2dPtVsMassKStar_FakeDaughters"), kstarpT, kstarMass);
  }

  template <int mode, typename TKStarObject, typename TCollision>
  void fillHistos(TKStarObject const& kstar, TCollision const& collision)
  {

    // Check whether it is before or after selections
    static constexpr std::string_view MainDir[] = {"BeforeSel", "AfterSel"};

    // Get V0trackCode
    int gammaTrkCode = retrieveV0TrackCode<true>(kstar);
    int kshortTrkCode = retrieveV0TrackCode<false>(kstar);

    float photonRZLineCut = std::abs(kstar.photonZconv()) * std::tan(2 * std::atan(std::exp(-photonSelections.photonMaxDauEta))) - photonSelections.photonLineCutZ0;
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
    //_______________________________________
    // Photon
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hTrackCode"), gammaTrkCode);
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hV0Type"), kstar.photonV0Type());

    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCANegToPV"), kstar.photonDCANegPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCAPosToPV"), kstar.photonDCAPosPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hDCADau"), kstar.photonDCADau());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosTPCCR"), kstar.photonPosTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegTPCCR"), kstar.photonNegTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosTPCNSigmaEl"), kstar.photonPosTPCNSigmaEl());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegTPCNSigmaEl"), kstar.photonNegTPCNSigmaEl());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hpT"), kstar.photonPt());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hY"), kstar.photonY());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPosEta"), kstar.photonPosEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hNegEta"), kstar.photonNegEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hRadius"), kstar.photonRadius());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hZ"), kstar.photonZconv());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dRZCut"), kstar.photonRadius(), photonRZLineCut);
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h2dRZPlane"), kstar.photonZconv(), kstar.photonRadius());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hCosPA"), kstar.photonCosPA());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPsiPair"), kstar.photonPsiPair());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hPhi"), kstar.photonPhi());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/h3dMass"), centrality, kstar.photonPt(), kstar.photonMass());
    histos.fill(HIST(MainDir[mode]) + HIST("/Photon/hMass"), kstar.photonMass());

    //_______________________________________
    // KShorts
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hTrackCode"), kshortTrkCode);
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hRadius"), kstar.kshortRadius());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hDCADau"), kstar.kshortDCADau());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hCosPA"), kstar.kshortCosPA());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hY"), kstar.kshortY());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hPosEta"), kstar.kshortPosEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hNegEta"), kstar.kshortNegEta());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hPosTPCCR"), kstar.kshortPosTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hNegTPCCR"), kstar.kshortNegTPCCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hPosITSCls"), kstar.kshortPosITSCls());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hNegITSCls"), kstar.kshortNegITSCls());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hPosChi2PerNc"), kstar.kshortPosChi2PerNcl());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hNegChi2PerNc"), kstar.kshortNegChi2PerNcl());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hLifeTime"), kstar.kshortLifeTime());

    //_______________________________________
    // Sigmas and KShorts
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenteros"), kstar.photonAlpha(), kstar.photonQt());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenteros"), kstar.kshortAlpha(), kstar.kshortQt());

    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/h2dTPCvsTOFNSigma_KShortPi"), kstar.kshortPosPiTPCNSigma(), kstar.kshortPiTOFNSigma());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/h2dTPCvsTOFNSigma_KShortPi"), kstar.kshortNegPiTPCNSigma(), kstar.kshortPiTOFNSigma());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hKShortDCANegToPV"), kstar.kshortDCANegPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hKShortDCAPosToPV"), kstar.kshortDCAPosPV());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hKShortpT"), kstar.kshortPt());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/hKShortMass"), kstar.kshortMass());
    histos.fill(HIST(MainDir[mode]) + HIST("/KShort/h3dKShortMass"), centrality, kstar.kshortPt(), kstar.kshortMass());

    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/hMass"), kstar.kstarMass());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/hPt"), kstar.pt());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/hY"), kstar.kstarY());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/hRadius"), kstar.radius());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/h2dRadiusVspT"), kstar.radius(), kstar.pt());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/hDCAPairDau"), kstar.dcadaughters());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/hDCAPairDauVsPt"), kstar.dcadaughters(), kstar.pt());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/h3dMass"), centrality, kstar.pt(), kstar.kstarMass());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/h3dOPAngleVsMass"), kstar.opAngle(), kstar.pt(), kstar.kstarMass());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/h2dOPAngleVsPt"), kstar.opAngle(), kstar.pt());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/h2dOPAngleKShortVsPt"), kstar.opAngleKShortKStar(), kstar.pt());
    histos.fill(HIST(MainDir[mode]) + HIST("/KStar/h2dOPAnglePhotonVsPt"), kstar.opAnglePhotonKStar(), kstar.pt());

    //_______________________________________
    // MC specific
    if (doprocessMonteCarlo) {
      if constexpr (requires { kstar.kshortPDGCode(); kstar.photonPDGCode(); }) {

        //_______________________________________
        // Gamma MC association
        if (kstar.photonPDGCode() == PDG_t::kGamma) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hV0ToCollAssoc"), kstar.photonIsCorrectlyAssoc());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hPt"), kstar.photonPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hMCPt"), kstar.photonmcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hPosTPCNSigmaEl"), kstar.photonPosTPCNSigmaEl());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hNegTPCNSigmaEl"), kstar.photonNegTPCNSigmaEl());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dPAVsPt"), std::acos(kstar.photonCosPA()), kstar.photonmcpt());

          if (!kstar.photonIsCorrectlyAssoc()) {
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/hPt_BadCollAssig"), kstar.photonmcpt());
            histos.fill(HIST(MainDir[mode]) + HIST("/MC/Photon/h2dPAVsPt_BadCollAssig"), std::acos(kstar.photonCosPA()), kstar.photonmcpt());
          }
        }

        //_______________________________________
        // KShort MC association
        if (kstar.kshortPDGCode() == PDG_t::kK0Short) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KShort/hV0ToCollAssoc"), kstar.kshortIsCorrectlyAssoc());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KShort/hPt"), kstar.kshortPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KShort/hMCPt"), kstar.kshortmcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KShort/hMass"), kstar.kshortMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KShort/h3dTPCvsTOFNSigma_Pi"), kstar.kshortPosPiTPCNSigma(), kstar.kshortPiTOFNSigma(), kstar.kshortPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KShort/h3dTPCvsTOFNSigma_Pi"), kstar.kshortNegPiTPCNSigma(), kstar.kshortPiTOFNSigma(), kstar.kshortPt());
        }

        //_______________________________________
        // KStar MC association
        if (kstar.isKStar()) {
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), kstar.photonAlpha(), kstar.photonQt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/h2dArmenteros"), kstar.kshortAlpha(), kstar.kshortQt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/hPt"), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/hMCPt"), kstar.mcpt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h2dMCPtVsKShortMCPt"), kstar.mcpt(), kstar.kshortmcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h2dMCPtVsPhotonMCPt"), kstar.mcpt(), kstar.photonmcpt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/hMass"), kstar.kstarMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h3dMass"), centrality, kstar.mcpt(), kstar.kstarMass());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/hMCProcess"), kstar.mcprocess());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/hGenRadius"), kstar.mcradius());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h2dMCProcessVsGenRadius"), kstar.mcprocess(), kstar.mcradius());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h3dMCProcess"), kstar.mcprocess(), kstar.mcpt(), kstar.kstarMass());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h2dOPAngleVsPt"), kstar.opAngle(), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h2dRadiusVspT"), kstar.radius(), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/hDCAPairDauVsPt"), kstar.dcadaughters(), kstar.pt());

          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h2dOPAngleKShortVsPt"), kstar.opAngleKShortKStar(), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/MC/KStar/h2dOPAnglePhotonVsPt"), kstar.opAnglePhotonKStar(), kstar.pt());
        }

        if ((kstar.kshortPDGCode() == PDG_t::kK0Short) && (kstar.photonPDGCode() == PDG_t::kGamma) && (std::abs(kstar.kshortPDGCodeGrandMother()) == o2::constants::physics::Pdg::kK0Star892 && std::abs(kstar.photonPDGCodeGrandMother()) == o2::constants::physics::Pdg::kK0Star892 && (kstar.photonGlobalIndexGrandMother() == kstar.kshortGlobalIndexGrandMother()))) { // K*(892)0
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hTrackCode"), gammaTrkCode);
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hV0Type"), kstar.photonV0Type());

          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hDCANegToPV"), kstar.photonDCANegPV());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hDCAPosToPV"), kstar.photonDCAPosPV());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hDCADau"), kstar.photonDCADau());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hPosTPCCR"), kstar.photonPosTPCCrossedRows());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hNegTPCCR"), kstar.photonNegTPCCrossedRows());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hPosTPCNSigmaEl"), kstar.photonPosTPCNSigmaEl());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hNegTPCNSigmaEl"), kstar.photonNegTPCNSigmaEl());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hpT"), kstar.photonPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hY"), kstar.photonY());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hPosEta"), kstar.photonPosEta());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hNegEta"), kstar.photonNegEta());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hRadius"), kstar.photonRadius());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hZ"), kstar.photonZconv());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/h2dRZCut"), kstar.photonRadius(), photonRZLineCut);
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/h2dRZPlane"), kstar.photonZconv(), kstar.photonRadius());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hCosPA"), kstar.photonCosPA());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hPsiPair"), kstar.photonPsiPair());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hPhi"), kstar.photonPhi());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/h3dMass"), centrality, kstar.photonPt(), kstar.photonMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/Photon/hMass"), kstar.photonMass());

          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hTrackCode"), kshortTrkCode);
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hRadius"), kstar.kshortRadius());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hDCADau"), kstar.kshortDCADau());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hCosPA"), kstar.kshortCosPA());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hY"), kstar.kshortY());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hPosEta"), kstar.kshortPosEta());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hNegEta"), kstar.kshortNegEta());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hPosTPCCR"), kstar.kshortPosTPCCrossedRows());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hNegTPCCR"), kstar.kshortNegTPCCrossedRows());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hPosITSCls"), kstar.kshortPosITSCls());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hNegITSCls"), kstar.kshortNegITSCls());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hPosChi2PerNc"), kstar.kshortPosChi2PerNcl());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hNegChi2PerNc"), kstar.kshortNegChi2PerNcl());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hLifeTime"), kstar.kshortLifeTime());

          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/h2dArmenteros"), kstar.photonAlpha(), kstar.photonQt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/h2dArmenteros"), kstar.kshortAlpha(), kstar.kshortQt());

          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/h2dTPCvsTOFNSigma_KShortPi"), kstar.kshortPosPiTPCNSigma(), kstar.kshortPiTOFNSigma());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/h2dTPCvsTOFNSigma_KShortPi"), kstar.kshortNegPiTPCNSigma(), kstar.kshortPiTOFNSigma());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hKShortDCANegToPV"), kstar.kshortDCANegPV());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hKShortDCAPosToPV"), kstar.kshortDCAPosPV());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hKShortpT"), kstar.kshortPt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/hKShortMass"), kstar.kshortMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KShort/h3dKShortMass"), centrality, kstar.kshortPt(), kstar.kshortMass());

          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/hMass"), kstar.kstarMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/hPt"), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/hY"), kstar.kstarY());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/hRadius"), kstar.radius());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/h2dRadiusVspT"), kstar.radius(), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/hDCAPairDau"), kstar.dcadaughters());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/hDCAPairDauVsPt"), kstar.dcadaughters(), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/h3dMass"), centrality, kstar.pt(), kstar.kstarMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/h3dOPAngleVsMass"), kstar.opAngle(), kstar.pt(), kstar.kstarMass());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/h2dOPAngleVsPt"), kstar.opAngle(), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/h2dOPAngleKShortVsPt"), kstar.opAngleKShortKStar(), kstar.pt());
          histos.fill(HIST(MainDir[mode]) + HIST("/PionKaon/KStar/h2dOPAnglePhotonVsPt"), kstar.opAnglePhotonKStar(), kstar.pt());
        }

        // For background studies:
        if (fillBkgQAhistos)
          runBkgAnalysis<mode>(kstar);

        //_______________________________________
        // pT resolution histos
        if (fillResoQAhistos)
          getResolution<mode>(kstar);
      }
    }
  }

  template <int selection_index, typename TKStarObject>
  void fillSelHistos(TKStarObject const& kstar, int PDGRequired)
  {

    static constexpr std::string_view photonSelsLocal[] = {"NoSel", "V0Type", "DCADauToPV",
                                                           "DCADau", "DauTPCCR", "TPCNSigmaEl", "V0pT",
                                                           "Y", "V0Radius", "RZCut", "Armenteros", "CosPA",
                                                           "PsiPair", "Phi", "Mass"};

    static constexpr std::string_view kshortSelsLocal[] = {"NoSel", "V0Radius", "DCADau", "Armenteros",
                                                           "CosPA", "Y", "TPCCR", "DauITSCls", "Lifetime",
                                                           "TPCTOFPID", "DCADauToPV", "Mass"};

    if (std::abs(PDGRequired) == PDG_t::kGamma) {
      if constexpr (selection_index >= 0 && selection_index < static_cast<int>(std::size(photonSelsLocal))) {
        histos.fill(HIST("Selection/Photon/hCandidateSel"), selection_index);
        histos.fill(HIST("Selection/Photon/h2d") + HIST(photonSelsLocal[selection_index]), kstar.photonPt(), kstar.photonMass());
        histos.fill(HIST("Selection/KStar/h2dPhoton") + HIST(photonSelsLocal[selection_index]), kstar.pt(), kstar.kstarMass());
      }
    }

    if (std::abs(PDGRequired) == PDG_t::kK0Short) {
      if constexpr (selection_index >= 0 && selection_index < static_cast<int>(std::size(kshortSelsLocal))) {
        histos.fill(HIST("Selection/KShort/hCandidateSel"), selection_index);
        histos.fill(HIST("Selection/KShort/h2d") + HIST(kshortSelsLocal[selection_index]), kstar.kshortPt(), kstar.kshortMass());
        histos.fill(HIST("Selection/KStar/h2dKShort") + HIST(kshortSelsLocal[selection_index]), kstar.pt(), kstar.kstarMass());
      }
    }
  }

  // Apply specific selections for photons
  template <typename TV0Object>
  bool selectPhoton(TV0Object const& cand)
  {
    fillSelHistos<0>(cand, PDG_t::kGamma);
    if (cand.photonV0Type() != photonSelections.photonv0TypeSel && photonSelections.photonv0TypeSel > -1)
      return false;

    fillSelHistos<1>(cand, PDG_t::kGamma);
    if ((std::abs(cand.photonDCAPosPV()) < photonSelections.photonMinDCADauToPv) || (std::abs(cand.photonDCANegPV()) < photonSelections.photonMinDCADauToPv))
      return false;

    fillSelHistos<2>(cand, PDG_t::kGamma);
    if (std::abs(cand.photonDCADau()) > photonSelections.photonMaxDCAV0Dau)
      return false;

    fillSelHistos<3>(cand, PDG_t::kGamma);
    if ((cand.photonPosTPCCrossedRows() < photonSelections.photonMinTPCCrossedRows) || (cand.photonNegTPCCrossedRows() < photonSelections.photonMinTPCCrossedRows))
      return false;

    fillSelHistos<4>(cand, PDG_t::kGamma);
    if (((cand.photonPosTPCNSigmaEl() < photonSelections.photonMinTPCNSigmas) || (cand.photonPosTPCNSigmaEl() > photonSelections.photonMaxTPCNSigmas)))
      return false;

    if (((cand.photonNegTPCNSigmaEl() < photonSelections.photonMinTPCNSigmas) || (cand.photonNegTPCNSigmaEl() > photonSelections.photonMaxTPCNSigmas)))
      return false;

    fillSelHistos<5>(cand, PDG_t::kGamma);
    if ((cand.photonPt() < photonSelections.photonMinPt) || (cand.photonPt() > photonSelections.photonMaxPt))
      return false;

    fillSelHistos<6>(cand, PDG_t::kGamma);
    if ((std::abs(cand.photonY()) > photonSelections.photonMaxRap) || (std::abs(cand.photonPosEta()) > photonSelections.photonMaxDauEta) || (std::abs(cand.photonNegEta()) > photonSelections.photonMaxDauEta))
      return false;

    fillSelHistos<7>(cand, PDG_t::kGamma);
    if ((cand.photonRadius() < photonSelections.photonMinRadius) || (cand.photonRadius() > photonSelections.photonMaxRadius))
      return false;

    fillSelHistos<8>(cand, PDG_t::kGamma);
    float photonRZLineCut = std::abs(cand.photonZconv()) * std::tan(2 * std::atan(std::exp(-photonSelections.photonMaxDauEta))) - photonSelections.photonLineCutZ0;
    if ((std::abs(cand.photonRadius()) < photonRZLineCut) || (std::abs(cand.photonZconv()) > photonSelections.photonMaxZ))
      return false;

    fillSelHistos<9>(cand, PDG_t::kGamma);
    if (cand.photonQt() > photonSelections.photonMaxQt)
      return false;

    if (std::abs(cand.photonAlpha()) > photonSelections.photonMaxAlpha)
      return false;

    fillSelHistos<10>(cand, PDG_t::kGamma);
    if (cand.photonCosPA() < photonSelections.photonMinV0cospa)
      return false;

    fillSelHistos<11>(cand, PDG_t::kGamma);
    if (std::abs(cand.photonPsiPair()) > photonSelections.photonPsiPairMax)
      return false;

    fillSelHistos<12>(cand, PDG_t::kGamma);
    if ((((cand.photonPhi() > photonSelections.photonPhiMin1) && (cand.photonPhi() < photonSelections.photonPhiMax1)) || ((cand.photonPhi() > photonSelections.photonPhiMin2) && (cand.photonPhi() < photonSelections.photonPhiMax2))) && ((photonSelections.photonPhiMin1 != -1) && (photonSelections.photonPhiMax1 != -1) && (photonSelections.photonPhiMin2 != -1) && (photonSelections.photonPhiMax2 != -1)))
      return false;

    fillSelHistos<13>(cand, PDG_t::kGamma);
    if (std::abs(cand.photonMass()) > photonSelections.photonMaxMass)
      return false;

    fillSelHistos<14>(cand, PDG_t::kGamma);
    return true;
  }

  // Apply specific selections for kshortrs
  template <typename TV0Object>
  bool selectKShort(TV0Object const& cand)
  {
    fillSelHistos<0>(cand, PDG_t::kK0Short);
    if ((cand.kshortRadius() < kshortSelections.kshortMinv0radius) || (cand.kshortRadius() > kshortSelections.kshortMaxv0radius))
      return false;

    fillSelHistos<1>(cand, PDG_t::kK0Short);
    if (std::abs(cand.kshortDCADau()) > kshortSelections.kshortMaxDCAV0Dau)
      return false;

    fillSelHistos<2>(cand, PDG_t::kK0Short);
    if ((cand.kshortQt() < kshortSelections.kshortMinQt) || (cand.kshortQt() > kshortSelections.kshortMaxQt))
      return false;

    if ((std::abs(cand.kshortAlpha()) < kshortSelections.kshortMinAlpha) || (std::abs(cand.kshortAlpha()) > kshortSelections.kshortMaxAlpha))
      return false;

    fillSelHistos<3>(cand, PDG_t::kK0Short);
    if (cand.kshortCosPA() < kshortSelections.kshortMinv0cospa)
      return false;

    fillSelHistos<4>(cand, PDG_t::kK0Short);
    if ((std::abs(cand.kshortY()) > kshortSelections.kshortMaxRap) || (std::abs(cand.kshortPosEta()) > kshortSelections.kshortMaxDauEta) || (std::abs(cand.kshortNegEta()) > kshortSelections.kshortMaxDauEta))
      return false;

    fillSelHistos<5>(cand, PDG_t::kK0Short);
    if ((cand.kshortPosTPCCrossedRows() < kshortSelections.kshortMinTPCCrossedRows) || (cand.kshortNegTPCCrossedRows() < kshortSelections.kshortMinTPCCrossedRows))
      return false;

    fillSelHistos<6>(cand, PDG_t::kK0Short);

    // check minimum number of ITS clusters + reject ITS afterburner tracks if requested
    bool posIsFromAfterburner = cand.kshortPosChi2PerNcl() < 0;
    bool negIsFromAfterburner = cand.kshortNegChi2PerNcl() < 0;
    if (cand.kshortPosITSCls() < kshortSelections.kshortMinITSclusters && (!kshortSelections.kshortRejectPosITSafterburner || posIsFromAfterburner))
      return false;
    if (cand.kshortNegITSCls() < kshortSelections.kshortMinITSclusters && (!kshortSelections.kshortRejectNegITSafterburner || negIsFromAfterburner))
      return false;

    fillSelHistos<7>(cand, PDG_t::kK0Short);
    if (cand.kshortLifeTime() > kshortSelections.kshortMaxLifeTime)
      return false;

    // Separating kshort selections:
    fillSelHistos<8>(cand, PDG_t::kK0Short);

    // TPC Selection
    if (kshortSelections.fselKShortTPCPID && (std::abs(cand.kshortPosPiTPCNSigma()) > kshortSelections.kshortMaxTPCNSigmas))
      return false;
    if (kshortSelections.fselKShortTPCPID && (std::abs(cand.kshortNegPiTPCNSigma()) > kshortSelections.kshortMaxTPCNSigmas))
      return false;

    // // TOF Selection
    // if (kshortSelections.fselKShortTOFPID && (std::abs(cand.kshortPiTOFNSigma()) > kshortSelections.kshortPiMaxTOFNSigmas))
    //   return false;
    // if (kshortSelections.fselKShortTOFPID && (std::abs(cand.lambdaPiTOFNSigma()) > kshortSelections.kshortPiMaxTOFNSigmas))
    //   return false;
    // DCA Selection
    fillSelHistos<9>(cand, PDG_t::kK0Short);
    if ((std::abs(cand.kshortDCAPosPV()) < kshortSelections.kshortMinDCAPosToPv) || (std::abs(cand.kshortDCANegPV()) < kshortSelections.kshortMinDCANegToPv))
      return false;

    // Mass Selection
    fillSelHistos<10>(cand, PDG_t::kK0Short);
    if (std::abs(cand.kshortMass() - o2::constants::physics::MassK0Short) > kshortSelections.kshortWindow)
      return false;

    fillSelHistos<11>(cand, PDG_t::kK0Short);

    return true;
  }

  // Apply selections in kdyst candidates
  template <typename TKStarObject>
  bool processKStarCandidate(TKStarObject const& cand)
  {
    // Photon specific selections
    if (!selectPhoton(cand))
      return false;

    // KShort specific selections
    if (!selectKShort(cand))
      return false;

    // KStar specific selections
    // Rapidity
    if constexpr (requires { cand.kstarMCY(); }) { // MC
      if (std::abs(cand.kstarMCY()) > kstarSelections.kstarMaxRap)
        return false;
    } else { // Real data
      if (std::abs(cand.kstarY()) > kstarSelections.kstarMaxRap)
        return false;
    }

    // V0Pair Radius
    if (cand.radius() > kstarSelections.kstarMaxRadius)
      return false;

    // DCA V0Pair Daughters
    if (cand.dcadaughters() > kstarSelections.kstarMaxDCADau)
      return false;

    // Opening Angle
    if (cand.opAngle() > kstarSelections.kstarMaxOPAngle)
      return false;

    return true;
  }

  // Main analysis function
  template <typename TCollisions, typename TKStars>
  void analyzeRecoeKStars(TCollisions const& collisions, TKStars const& fullKStars)
  {
    // Custom grouping
    std::vector<std::vector<int>> kstargrouped(collisions.size());

    for (const auto& kstar : fullKStars) {
      kstargrouped[kstar.straCollisionId()].push_back(kstar.globalIndex());
    }

    // Collisions loop
    for (const auto& coll : collisions) {

      // Event selection
      if (!isEventAccepted(coll, true))
        continue;

      // KStars loop
      for (size_t i = 0; i < kstargrouped[coll.globalIndex()].size(); i++) {
        auto kstar = fullKStars.rawIteratorAt(kstargrouped[coll.globalIndex()][i]);

        // if MC
        if constexpr (requires { kstar.isKStar(); }) {
          if (doMCAssociation && !(kstar.isKStar()))
            continue;

          if (selRecoFromGenerator && !kstar.isProducedByGenerator())
            continue;
        }

        // Fill histos before any selection
        fillHistos<0>(kstar, coll);

        // Select kstar candidates
        if (!processKStarCandidate(kstar))
          continue;

        // Fill histos after all selections
        fillHistos<1>(kstar, coll);
      }
    }
  }

  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraEvSelExtras, aod::StraStamps> const& collisions, KStars const& fullKStars)
  {
    analyzeRecoeKStars(collisions, fullKStars);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraEvSelExtras, aod::StraStamps, aod::StraCollLabels> const& collisions, MCKStars const& fullKStars)
  {
    analyzeRecoeKStars(collisions, fullKStars);
  }

  // Simulated processing in Run 3
  void processGeneratedRun3(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraEvSelExtras, aod::StraStamps, aod::StraCollLabels> const& collisions, soa::Join<aod::KStarGens, aod::KStarGenCollRef> const& KStarGens)
  {
    analyzeGenerated(mcCollisions, collisions, KStarGens);
  }

  // _____________________________________________________
  PROCESS_SWITCH(k892hadronphoton, processRealData, "Do real data analysis", true);
  PROCESS_SWITCH(k892hadronphoton, processMonteCarlo, "Do Monte-Carlo-based analysis", false);
  PROCESS_SWITCH(k892hadronphoton, processGeneratedRun3, "process MC generated Run 3", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892hadronphoton>(cfgc)};
}
