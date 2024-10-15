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
/// \ post processing for Cascade analysis runing on derived data
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/StaticFor.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

// constants
const float ctauxiPDG = 4.91;     // from PDG
const float ctauomegaPDG = 2.461; // from PDG

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using cascMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascBBs, aod::CascCoreMCLabels>;

struct derivedCascadeAnalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "occupancy axis"};

  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};
  ConfigurableAxis axisMass{"axisMass", {200, 1.222f, 1.422f}, "range of invariant mass, in case of omega take 1.572f, 1.772f"};

  Configurable<bool> isXi{"isXi", 1, "Apply cuts for Xi identification"};
  Configurable<bool> doBefSelCheck{"doBefSelCheck", false, "Fill mass histograms of all candidates before selections"};
  Configurable<bool> doPtDepCutStudy{"doPtDepCutStudy", false, "Fill histogram with a cutting paramer"};
  Configurable<bool> doTriggerTVXEventCut{"doTriggerTVXEventCut", false, "Minimal MB event selection, for MC"};
  Configurable<bool> doTriggerSel8EventCut{"doTriggerSel8EventCut", true, "Standard MB event selection"};
  Configurable<bool> doGoodPVFT0EventCut{"doGoodPVFT0EventCut", true, "check for the PV position diffrence when estimated from tracks and FT0"};
  Configurable<bool> doITSTPCvertexEventCut{"doITSTPCvertexEventCut", true, "checks the presence of at least one ITS-TPC track"};
  Configurable<bool> doSameBunchPileUpEventCut{"doSameBunchPileUpEventCut", true, "removes events  associated with the same \"found-by-T0\" bunch crossing"};
  Configurable<bool> doVertexTOFmatch{"doVertexTOFmatch", false, "Checks wherher at least one of vertex contributors is matched to TOF"};
  Configurable<bool> doVertexTRDmatch{"doVertexTRDmatch", false, "Checks wherher at least one of vertex contributors is matched to TRD"};
  Configurable<bool> doBefSelEventMultCorr{"doBefSelEventMultCorr", false, "Enable histogram of multiplicity correlation before cuts"};
  Configurable<bool> doTFeventCut{"doTFeventCut", false, "Enable TF event Cut"};
  Configurable<bool> doITSFrameBorderCut{"doITSFrameBorderCut", false, "Enable ITSFrame event cut"};
  Configurable<bool> doMultiplicityCorrCut{"doMultiplicityCorrCut", false, "Enable multiplicity vs centrality correlation cut"};
  Configurable<bool> doOccupancyCheck{"doOccupancyCheck", true, ""};
  Configurable<bool> doTimeRangeStandardCut{"doTimeRangeStandardCut", true, "It rejects a given collision if there are other events nearby in |dt|< 10 μs"};
  Configurable<bool> doTimeRangeNarrowCut{"doTimeRangeNarrowCut", false, "It rejects a given collision if there are other events nearby in |dt|< 4 μs"};
  Configurable<int> minOccupancy{"minOccupancy", -1, "Minimal occupancy"};
  Configurable<int> maxOccupancy{"maxOccupancy", -1, "Maximal occupancy"};

  Configurable<int> centMin{"centMin", 0, "Minimal accepted centrality"};
  Configurable<int> centMax{"centMax", 100, "Maximal accepted centrality"};

  Configurable<float> minPt{"minPt", 0.0f, "minPt"};
  Configurable<float> masswin{"masswin", 0.05, "Mass window limit"};
  Configurable<float> lambdaMassWin{"lambdaMassWin", 0.005, "V0 Mass window limit"};
  Configurable<float> rapCut{"rapCut", 0.5, "Rapidity acceptance"};
  Configurable<float> etaDauCut{"etaDauCut", 0.8, "Pseudorapidity acceptance of the cascade daughters"};
  Configurable<float> dcaBaryonToPV{"dcaBaryonToPV", 0.05, "DCA of baryon doughter track To PV"};
  Configurable<float> dcaMesonToPV{"dcaMesonToPV", 0.1, "DCA of meson doughter track To PV"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.04, "DCA Bach To PV"};
  Configurable<double> casccospa{"casccospa", 0.97, "Casc CosPA"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcacascdau{"dcacascdau", 1., "DCA Casc Daughters"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcaV0ToPV{"dcaV0ToPV", 0.06, "DCA V0 To PV"};
  Configurable<float> minRadius{"minRadius", 1.4f, "minRadius"};
  Configurable<float> maxRadius{"maxRadius", 100.0f, "maxRadius"};
  Configurable<float> minV0Radius{"minV0Radius", 1.2f, "V0 transverse decay radius, minimum"};
  Configurable<float> maxV0Radius{"maxV0Radius", 100.0f, "V0 transverse decay radius, maximum"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 5, "N sigma TPC Pion"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 5, "N sigma TPC Proton"};
  Configurable<float> nsigmatpcKa{"nsigmatpcKa", 5, "N sigma TPC Kaon"};
  Configurable<float> nsigmatofPr{"nsigmatofPr", 3, "N sigma TOF Proton"};
  Configurable<float> nsigmatofPion{"nsigmatofPion", 3, "N sigma TOF for Pion from V0"};
  Configurable<float> nsigmatofBachPion{"nsigmatofBachPion", 3, "N sigma TOF for bachelor Pion"};
  Configurable<float> nsigmatofBachKaon{"nsigmatofBachKaon", 3, "N sigma TOF for bachelor Kaon"};
  Configurable<float> bachBaryonCosPA{"bachBaryonCosPA", 0.9999, "Bachelor baryon CosPA"};
  Configurable<float> bachBaryonDCAxyToPV{"bachBaryonDCAxyToPV", 0.08, "DCA bachelor baryon to PV"};
  Configurable<int> mintpccrrows{"mintpccrrows", 50, "min N TPC crossed rows"};
  Configurable<int> dooobrej{"dooobrej", 0, "OOB rejection: 0 no selection, 1 = ITS||TOF, 2 = TOF only for pT > ptthrtof"};
  Configurable<float> ptthrtof{"ptthrtof", 2, "Pt threshold for applying only tof oob rejection"};
  Configurable<float> proplifetime{"proplifetime", 3, "ctau/<ctau>"};
  Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};

  Configurable<bool> doPtDepCosPaCut{"doPtDepCosPaCut", false, "Enable pt dependent cos PA cut"};
  Configurable<bool> doPtDepCascRadiusCut{"doPtDepCascRadiusCut", false, "Enable pt dependent cascade radius cut"};
  Configurable<bool> doPtDepV0RadiusCut{"doPtDepV0RadiusCut", false, "Enable pt dependent V0 radius cut"};
  Configurable<bool> doPtDepV0CosPaCut{"doPtDepV0CosPaCut", false, "Enable pt dependent cos PA cut of the V0 daughter"};
  Configurable<bool> doPtDepDCAcascDauCut{"doPtDepDCAcascDauCut", false, "Enable pt dependent DCA cascade daughter cut"};
  Configurable<bool> doDCAdauToPVCut{"doDCAdauToPVCut", true, "Enable cut DCA daughter track to PV"};
  Configurable<bool> doCascadeCosPaCut{"doCascadeCosPaCut", true, "Enable cos PA cut"};
  Configurable<bool> doV0CosPaCut{"doV0CosPaCut", true, "Enable cos PA cut for the V0 daughter"};
  Configurable<bool> doDCACascadeDauCut{"doDCACascadeDauCut", true, "Enable cut DCA betweenn daughter tracks"};
  Configurable<bool> doDCAV0DauCut{"doDCAV0DauCut", true, "Enable cut DCA betweenn V0 daughter tracks"};
  Configurable<bool> doCascadeRadiusCut{"doCascadeRadiusCut", true, "Enable cut on the cascade radius"};
  Configurable<bool> doV0RadiusCut{"doV0RadiusCut", true, "Enable cut on the V0 radius"};
  Configurable<bool> doDCAV0ToPVCut{"doDCAV0ToPVCut", true, "Enable cut DCA of V0 to PV"};
  Configurable<bool> doNTPCSigmaCut{"doNTPCSigmaCut", false, "Enable cut N sigma TPC"};
  Configurable<bool> doBachelorBaryonCut{"doBachelorBaryonCut", true, "Enable Bachelor-Baryon cut "};
  Configurable<bool> doProperLifeTimeCut{"doProperLifeTimeCut", true, "Enable proper life-time cut "};
  Configurable<bool> doNTOFSigmaProtonCut{"doNTOFSigmaProtonCut", true, "Enable n sigma TOF PID cut for proton from V0"};
  Configurable<bool> doNTOFSigmaV0PionCut{"doNTOFSigmaV0PionCut", false, "Enable n sigma TOF PID cut for pion from V0"};
  Configurable<bool> doNTOFSigmaBachelorCut{"doNTOFSigmaBachelorCut", false, "Enable n sigma TOF PID cut for bachelor track"};
  Configurable<bool> doFillNsigmaTPCHistPionBach{"doFillNsigmaTPCHistPionBach", false, ""};
  Configurable<bool> doFillNsigmaTPCHistKaonBach{"doFillNsigmaTPCHistKaonBach", false, ""};
  Configurable<bool> doFillNsigmaTPCHistProton{"doFillNsigmaTPCHistProton", false, ""};
  Configurable<bool> doFillNsigmaTPCHistV0Pion{"doFillNsigmaTPCHistV0Pion", false, ""};
  Configurable<int> nPtBinsForNsigmaTPC{"nPtBinsForNsigmaTPC", 100, ""};

  Configurable<float> cosPApar0{"cosPApar0", 0.2, "const par for pt dep cosPA cut"};
  Configurable<float> cosPApar1{"cosPApar1", -0.022, "linear par for pt dep cosPA cut"};

  Configurable<float> parCascRadius0{"parCascRadius0", 1.216159, "const par for pt dep radius cut"};
  Configurable<float> parCascRadius1{"parCascRadius1", 0.064462, "linear par for pt dep radius cut"};

  Configurable<float> parV0Radius0{"parV0Radius0", 2.136381, "const par for pt dep V0 radius cut"};
  Configurable<float> parV0Radius1{"parV0Radius1", 0.437074, "linear par for pt dep V0 radius cut"};

  Configurable<float> dcaCacsDauPar0{"dcaCacsDauPar0", 0.8, " par for pt dep DCA cascade daughter cut, p_T < 1 GeV/c"};
  Configurable<float> dcaCacsDauPar1{"dcaCacsDauPar1", 0.5, " par for pt dep DCA cascade daughter cut, 1< p_T < 4 GeV/c"};
  Configurable<float> dcaCacsDauPar2{"dcaCacsDauPar2", 0.2, " par for pt dep DCA cascade daughter cut, p_T > 4 GeV/c"};

  ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles in |y| < 0.5"};

  Service<o2::framework::O2DatabasePDG> pdgDB;

  static constexpr std::string_view Index[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
  static constexpr float centralityIntervals[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  static constexpr std::string_view charge[] = {"Positive", "Negative"};
  static constexpr std::string_view selectionNames[] = {"BachelorBaryonDCA", "DCAV0ToPV", "V0Radius", "CascadeRadius", "DCAV0Daughters", "DCACascDaughters", "V0pa", "CascPA", "DCABachelorToPV", "DCAMesonToPV", "DCABaryonToPV", "CascadeProperLifeTime"};

  // For manual sliceBy
  // Preslice<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

  void init(InitContext const&)
  {
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101, 0, 101}});
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{17, 0, 17}});
    histos.add("hOccupancyVsCentrality", "", kTH2F, {axisOccupancy, {100, 0, 100}});

    histos.add("hEventNchCorrelationAfCuts", "hEventNchCorrelationAfCuts", kTH2F, {{5000, 0, 5000}, {5000, 0, 2500}});
    histos.add("hEventPVcontributorsVsCentrality", "hEventPVcontributorsVsCentrality", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
    histos.add("hEventGlobalTracksVsCentrality", "hEventGlobalTracksVsCentrality", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    if (doBefSelEventMultCorr) {
      histos.add("hEventNchCorrelationBefCuts", "hEventNchCorrelationBefCuts", kTH2F, {{5000, 0, 5000}, {2500, 0, 2500}});
      histos.add("hEventPVcontributorsVsCentralityBefCuts", "hEventPVcontributorsVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
      histos.add("hEventGlobalTracksVsCentralityBefCuts", "hEventGlobalTracksVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    }

    histos.add("hCandidate", "hCandidate", HistType::kTH1D, {{22, -0.5, 21.5}});
    histos.add("hCutValue", "hCutValue", HistType::kTH2D, {{22, -0.5, 21.5}, {300, 0, 3.01}});

    if (doFillNsigmaTPCHistProton) {
      histos.add("hNsigmaProton", "hNsigmaProton", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 6}, {100, 0, 100}});
      histos.add("hNsigmaProtonNeg", "hNsigmaProtonNeg", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 6}, {100, 0, 100}});
    }
    if (doFillNsigmaTPCHistV0Pion) {
      histos.add("hNsigmaPionNeg", "hNsigmaPionNeg", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
      histos.add("hNsigmaPionPos", "hNsigmaPionPos", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
    }
    if (doFillNsigmaTPCHistPionBach) {
      histos.add("hNsigmaPionPosBach", "hNsigmaPionPosBach", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
      histos.add("hNsigmaPionNegBach", "hNsigmaPionNegBach", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
    }
    if (doFillNsigmaTPCHistKaonBach)
      histos.add("hNsigmaKaon", "hNsigmaKaon", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 6}, {100, 0, 100}});

    if (doNTOFSigmaProtonCut)
      histos.add("hNsigmaTOFProton", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
    if (doNTOFSigmaV0PionCut)
      histos.add("hNsigmaTOFV0Pion", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
    if (doNTOFSigmaBachelorCut) {
      if (isXi)
        histos.add("hNsigmaTOFBachelorPion", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
      else
        histos.add("hNsigmaTOFBachelorKaon", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
    }

    TString CutLabel[22] = {"All", "MassWin", "y", "DCACascDau", "DCAV0Dau", "rCasc", "rCascMax", "rV0", "rV0Max", "LambdaMass", "Bach-baryon", "V0CosPA", "CompDecayMass", "DCADauToPV", "EtaDau", "CascCosPA", "DCAV0ToPV", "nSigmaTPCV0Dau", "NTPCrows", "OOBRej", "nSigmaTPCbachelor", "ctau"};
    for (Int_t i = 1; i <= histos.get<TH1>(HIST("hCandidate"))->GetNbinsX(); i++) {
      histos.get<TH1>(HIST("hCandidate"))->GetXaxis()->SetBinLabel(i, CutLabel[i - 1]);
      histos.get<TH2>(HIST("hCutValue"))->GetXaxis()->SetBinLabel(i, CutLabel[i - 1]);
    }

    histos.add("InvMassAfterSel/hNegativeCascade", "hNegativeCascade", HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}});
    histos.add("InvMassAfterSel/hPositiveCascade", "hPositiveCascade", {HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}}});

    if (doOccupancyCheck) {
      histos.add("InvMassAfterSelCent1/hNegativeCascade", "hNegativeCascade", HistType::kTH3F, {axisPt, axisMass, axisOccupancy});
      histos.add("InvMassAfterSelCent1/hPositiveCascade", "hPositiveCascade", HistType::kTH3F, {axisPt, axisMass, axisOccupancy});
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent2/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent3/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent4/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent5/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent6/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent7/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent8/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent9/");
      histos.addClone("InvMassAfterSelCent1/", "InvMassAfterSelCent10/");
    }

    if (doBefSelCheck)
      histos.addClone("InvMassAfterSel/", "InvMassBefSel/");

    if (doprocessCascadesMCrec)
      histos.addClone("InvMassAfterSel/", "InvMassAfterSelMCrecTruth/");

    if (doPtDepCutStudy && !doProperLifeTimeCut) {
      histos.add("PtDepCutStudy/hNegativeCascadeProperLifeTime", "hNegativeCascadeProperLifeTime", HistType::kTH3F, {axisPt, axisMass, {100, 0, 10}});
      histos.add("PtDepCutStudy/hPositiveCascadeProperLifeTime", "hPositiveCascadeProperLifeTime", {HistType::kTH3F, {axisPt, axisMass, {100, 0, 10}}});
    }
    if (doPtDepCutStudy && !doBachelorBaryonCut) {
      histos.add("PtDepCutStudy/hNegativeBachelorBaryonDCA", "hNegativeBachelorBaryonDCA", HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}});
      histos.add("PtDepCutStudy/hPositiveBachelorBaryonDCA", "hPositiveBachelorBaryonDCA", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}}});
    }
    if (doPtDepCutStudy && !doDCAV0ToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCAV0ToPV", "hNegativeDCAV0ToPV", HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}});
      histos.add("PtDepCutStudy/hPositiveDCAV0ToPV", "hPositiveDCAV0ToPV", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}}});
    }
    if (doPtDepCutStudy && !doV0RadiusCut) {
      histos.add("PtDepCutStudy/hNegativeV0Radius", "hNegativeV0Radius", HistType::kTH3F, {axisPt, axisMass, {20, 0, 10}});
      histos.add("PtDepCutStudy/hPositiveV0Radius", "hPositiveV0Radius", {HistType::kTH3F, {axisPt, axisMass, {20, 0, 10}}});
    }
    if (doPtDepCutStudy && !doCascadeRadiusCut) {
      histos.add("PtDepCutStudy/hNegativeCascadeRadius", "hNegativeCascadeRadius", HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}});
      histos.add("PtDepCutStudy/hPositiveCascadeRadius", "hPositiveCascadeRadius", {HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}}});
    }

    if (doPtDepCutStudy && !doDCAV0DauCut) {
      histos.add("PtDepCutStudy/hNegativeDCAV0Daughters", "hNegativeDCAV0Daughters", HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}});
      histos.add("PtDepCutStudy/hPositiveDCAV0Daughters", "hPositiveDCAV0Daughters", {HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}}});
    }

    if (doPtDepCutStudy && !doDCACascadeDauCut) {
      histos.add("PtDepCutStudy/hNegativeDCACascDaughters", "hNegativeDCACascDaughters", {HistType::kTH3F, {axisPt, axisMass, {20, 0, 1}}});
      histos.add("PtDepCutStudy/hPositiveDCACascDaughters", "hPositiveDCACascDaughters", {HistType::kTH3F, {axisPt, axisMass, {20, 0, 1}}});
    }

    if (doPtDepCutStudy && !doV0CosPaCut) {
      histos.add("PtDepCutStudy/hNegativeV0pa", "hNegativeV0pa", HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}});
      histos.add("PtDepCutStudy/hPositiveV0pa", "hPositiveV0pa", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}}});
    }
    if (doPtDepCutStudy && !doDCAdauToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCABachelorToPV", "hNegativeDCABachelorToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hNegativeDCABaryonToPV", "hNegativeDCABaryonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hNegativeDCAMesonToPV", "hNegativeDCAMesonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hPositiveDCABachelorToPV", "hPositiveDCABachelorToPV", {HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}}});
      histos.add("PtDepCutStudy/hPositiveDCABaryonToPV", "hPositiveDCABaryonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hPositiveDCAMesonToPV", "hPositiveDCAMesonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
    }

    if (doPtDepCutStudy && !doCascadeCosPaCut) {
      histos.add("PtDepCutStudy/hNegativeCascPA", "hNegativeCascPA", HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}});
      histos.add("PtDepCutStudy/hPositiveCascPA", "hPositiveCascPA", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}}});
    }

    if (doprocessCascadesMCrec) {
      histos.addClone("PtDepCutStudy/", "PtDepCutStudyMCTruth/");
      histos.add("hNegativeCascadePtForEfficiency", "hNegativeCascadePtForEfficiency", HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}});
      histos.add("hPositiveCascadePtForEfficiency", "hPositiveCascadePtForEfficiency", {HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}}});
    }

    if (doprocessCascadesMCforEff) {
      histos.add("hGenEvents", "", HistType::kTH2F, {{axisNch}, {2, 0, 2}});
      histos.add("hCentralityVsMultMC", "", kTH2F, {{101, 0.0f, 101.0f}, axisNch});

      histos.add("hCentralityVsNcoll_beforeEvSel", "", kTH2F, {{101, 0.0f, 101.0f}, {50, 0.f, 50.f}});
      histos.add("hCentralityVsNcoll_afterEvSel", "", kTH2F, {{101, 0.0f, 101.0f}, {50, 0.f, 50.f}});

      histos.add("h2dGenXiMinus", "h2dGenXiMinus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});
      histos.add("h2dGenXiPlus", "h2dGenXiPlus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});
      histos.add("h2dGenOmegaMinus", "h2dGenOmegaMinus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});
      histos.add("h2dGenOmegaPlus", "h2dGenOmegaPlus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});

      histos.add("h2dGenXiMinusVsMultMC", "h2dGenXiMinusVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("h2dGenXiPlusVsMultMC", "h2dGenXiPlusVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("h2dGenOmegaMinusVsMultMC", "h2dGenOmegaMinusVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("h2dGenOmegaPlusVsMultMC", "h2dGenOmegaPlusVsMultMC", kTH2D, {axisNch, axisPt});
    }
  }
  template <typename TCascade>
  bool IsCosPAAccepted(TCascade casc, float x, float y, float z, bool ptdepcut, bool isCascPa)
  {

    if (ptdepcut) {
      double ptdepCut;
      if (isCascPa)
        ptdepCut = cosPApar0 + cosPApar1 * casc.pt();
      else
        ptdepCut = cosPApar0 + cosPApar1 * casc.pt();
      if (ptdepCut > 0.3 && casc.pt() < 0.5)
        ptdepCut = 0.3;
      if (ptdepCut < 0.012)
        ptdepCut = 0.012;
      if (isCascPa)
        histos.fill(HIST("hCutValue"), 16, TMath::Cos(ptdepCut));
      if (!isCascPa)
        histos.fill(HIST("hCutValue"), 12, TMath::Cos(ptdepCut));
      if (isCascPa && casc.casccosPA(x, y, z) < TMath::Cos(ptdepCut))
        return false;
      if (!isCascPa && casc.v0cosPA(x, y, z) < TMath::Cos(ptdepCut))
        return false;
    } else {
      float cut = casccospa;
      float cutV0 = v0cospa;
      if (isCascPa)
        histos.fill(HIST("hCutValue"), 16, cut);
      if (!isCascPa)
        histos.fill(HIST("hCutValue"), 12, cutV0);
      if (isCascPa && casc.casccosPA(x, y, z) < casccospa)
        return false;
      if (!isCascPa && casc.v0cosPA(x, y, z) < v0cospa)
        return false;
    }

    return true;
  }

  template <typename TCollision>
  bool IsEventAccepted(TCollision coll, bool fillHists)
  {

    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);

    if (doBefSelEventMultCorr) {
      histos.fill(HIST("hEventNchCorrelationBefCuts"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
      histos.fill(HIST("hEventPVcontributorsVsCentralityBefCuts"), coll.centFT0C(), coll.multNTracksPVeta1());
      histos.fill(HIST("hEventGlobalTracksVsCentralityBefCuts"), coll.centFT0C(), coll.multNTracksGlobal());
    }

    if (doTriggerTVXEventCut && !coll.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1.5 /* collisions  after sel*/);

    if (doTriggerSel8EventCut && !coll.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2.5 /* collisions  after sel*/);

    if (TMath::Abs(coll.posZ()) > zVertexCut) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3.5 /* collisions  after sel pvz sel*/);

    if (coll.centFT0C() > centMax || coll.centFT0C() < centMin) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4.5 /* collisions  after centrality sel*/);

    if (doSameBunchPileUpEventCut && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5.5 /* Not same Bunch pile up */);

    if (doGoodPVFT0EventCut && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6.5 /* No large vertexZ difference from tracks and FT0*/);

    if (doITSTPCvertexEventCut && !coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7.5 /* At least one ITS-TPC track in the event*/);

    if (doVertexTOFmatch && !coll.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8.5 /* At least one of vertex contributors is matched to TOF*/);

    if (doVertexTRDmatch && !coll.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9.5 /* At least one of vertex contributors is matched to TRD*/);

    if (doITSFrameBorderCut && !coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 10.5 /* Not at ITS ROF border */);

    if (doTFeventCut && !coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 11.5 /* Not at TF border */);

    if (doMultiplicityCorrCut) {
      if (coll.multNTracksGlobal() < (1343.3 * TMath::Exp(-0.0443259 * coll.centFT0C()) - 50) || coll.multNTracksGlobal() > (2098.9 * TMath::Exp(-0.0332444 * coll.centFT0C())))
        return false;
      if (coll.multNTracksPVeta1() < (3703 * TMath::Exp(-0.0455483 * coll.centFT0C()) - 150) || coll.multNTracksPVeta1() > (4937.33 * TMath::Exp(-0.0372668 * coll.centFT0C()) + 20))
        return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 12.5 /* Remove outlyers */);

    if (doTimeRangeStandardCut && !coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 13.5 /* Rejection of events too close in time */);

    if (doTimeRangeNarrowCut && !coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 14.5 /* No other collision within +/- 4 microseconds */);

    int occupancy = coll.trackOccupancyInTimeRange();
    if (minOccupancy > 0 && occupancy < minOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 15.5 /* Below min occupancy */);
    if (maxOccupancy > 0 && occupancy > maxOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 16.5 /* Above max occupancy */);

    if (fillHists) {
      histos.fill(HIST("hOccupancyVsCentrality"), occupancy, coll.centFT0C());
      histos.fill(HIST("hEventCentrality"), coll.centFT0C());
      histos.fill(HIST("hEventVertexZ"), coll.posZ());
      histos.fill(HIST("hEventNchCorrelationAfCuts"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
      histos.fill(HIST("hEventPVcontributorsVsCentrality"), coll.centFT0C(), coll.multNTracksPVeta1());
      histos.fill(HIST("hEventGlobalTracksVsCentrality"), coll.centFT0C(), coll.multNTracksGlobal());
    }

    return true;
  }

  template <typename TCascade>
  bool IsCascadeCandidateAccepted(TCascade casc, int counter, float /*centrality*/)
  {
    float cut = masswin;
    histos.fill(HIST("hCutValue"), 2, cut);
    cut = rapCut;
    histos.fill(HIST("hCutValue"), 3, cut);
    if (isXi) {
      if (TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) > masswin) {
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (TMath::Abs(casc.yXi()) > rapCut)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) > masswin) {
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (TMath::Abs(casc.yOmega()) > rapCut)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    }

    if (doDCACascadeDauCut) {
      if (doPtDepDCAcascDauCut) {
        float ptDepCut = dcaCacsDauPar0;
        if (casc.pt() > 1 && casc.pt() < 4)
          ptDepCut = dcaCacsDauPar1;
        else if (casc.pt() > 4)
          ptDepCut = dcaCacsDauPar2;
        histos.fill(HIST("hCutValue"), 4, ptDepCut);
        if (casc.dcacascdaughters() > ptDepCut)
          return false;
      } else {
        cut = dcacascdau;
        histos.fill(HIST("hCutValue"), 4, cut);
        if (casc.dcacascdaughters() > dcacascdau)
          return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      ++counter;
    }

    if (doDCAV0DauCut) {
      cut = dcav0dau;
      histos.fill(HIST("hCutValue"), 5, cut);
      if (casc.dcaV0daughters() > dcav0dau)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      ++counter;
    }

    if (doCascadeRadiusCut) {
      if (doPtDepCascRadiusCut) {
        double ptdepminRadius = parCascRadius0 + parCascRadius1 * casc.pt();
        histos.fill(HIST("hCutValue"), 6, ptdepminRadius);
        if (casc.cascradius() < ptdepminRadius)
          return false;
      } else {
        cut = minRadius;
        histos.fill(HIST("hCutValue"), 6, cut);
        if (casc.cascradius() < minRadius)
          return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);

      if (casc.cascradius() > maxRadius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      counter += 2;
    }

    if (doV0RadiusCut) {
      if (doPtDepV0RadiusCut) {
        float cut = parV0Radius0 + casc.pt() * parV0Radius1;
        histos.fill(HIST("hCutValue"), 8, cut);
        if (casc.v0radius() < cut)
          return false;
      } else {
        cut = minV0Radius;
        histos.fill(HIST("hCutValue"), 8, cut);
        if (casc.v0radius() < minV0Radius)
          return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (casc.v0radius() > maxV0Radius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      counter += 2;
    }

    cut = lambdaMassWin;
    histos.fill(HIST("hCutValue"), 10, cut);
    if (TMath::Abs(casc.mLambda() - pdgDB->Mass(3122)) > lambdaMassWin)
      return false;
    histos.fill(HIST("hCandidate"), ++counter);

    cut = bachBaryonDCAxyToPV;
    histos.fill(HIST("hCutValue"), 11, cut);
    if (doBachelorBaryonCut) {
      if ((casc.bachBaryonCosPA() > bachBaryonCosPA || TMath::Abs(casc.bachBaryonDCAxyToPV()) < bachBaryonDCAxyToPV)) { // Bach-baryon selection if required
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      ++counter;
    }

    if (doV0CosPaCut) {
      if (!IsCosPAAccepted(casc, casc.x(), casc.y(), casc.z(), doPtDepV0CosPaCut, false))
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      ++counter;
    }

    cut = rejcomp;
    histos.fill(HIST("hCutValue"), 13, cut);
    if (isXi) {
      if (TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) < rejcomp)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) < rejcomp)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    }

    cut = dcaBachToPV;
    histos.fill(HIST("hCutValue"), 14, cut);
    cut = dcaMesonToPV;
    histos.fill(HIST("hCutValue"), 14, cut);
    cut = dcaBaryonToPV;
    histos.fill(HIST("hCutValue"), 14, cut);
    if (doDCAdauToPVCut) {
      if (TMath::Abs(casc.dcabachtopv()) < dcaBachToPV)
        return false;
      if (casc.sign() > 0 && (TMath::Abs(casc.dcanegtopv()) < dcaBaryonToPV || TMath::Abs(casc.dcapostopv()) < dcaMesonToPV))
        return false;
      if (casc.sign() < 0 && (TMath::Abs(casc.dcapostopv()) < dcaBaryonToPV || TMath::Abs(casc.dcanegtopv()) < dcaMesonToPV))
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      ++counter;
    }

    return true;
  }

  void processCascades(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFNSigmas> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {

    if (!IsEventAccepted(coll, true))
      return;

    for (auto& casc : Cascades) {

      int counter = -1;
      histos.fill(HIST("hCandidate"), ++counter);

      bool isNegative = false;
      bool isPositive = false;
      if (casc.sign() > 0)
        isPositive = true;
      if (casc.sign() < 0)
        isNegative = true;
      if (!isNegative && !isPositive)
        continue;

      double invmass;
      if (isXi)
        invmass = casc.mXi();
      else
        invmass = casc.mOmega();
      // To have trace of how it was before selections

      if (doBefSelCheck) {
        if (isPositive)
          histos.fill(HIST("InvMassBefSel/h") + HIST(charge[0]) + HIST("Cascade"), casc.pt(), invmass, coll.centFT0C());
        if (isNegative)
          histos.fill(HIST("InvMassBefSel/h") + HIST(charge[1]) + HIST("Cascade"), casc.pt(), invmass, coll.centFT0C());
      }

      if (!IsCascadeCandidateAccepted(casc, counter, coll.centFT0C()))
        continue;
      counter += 13;

      auto negExtra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto bachExtra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      auto poseta = RecoDecay::eta(std::array{casc.pxpos(), casc.pypos(), casc.pzpos()});
      auto negeta = RecoDecay::eta(std::array{casc.pxneg(), casc.pyneg(), casc.pzneg()});
      auto bacheta = RecoDecay::eta(std::array{casc.pxbach(), casc.pybach(), casc.pzbach()});

      auto fullMomentumPosDaugh = TMath::Sqrt(TMath::Power(casc.pxpos(), 2) + TMath::Power(casc.pypos(), 2) + TMath::Power(casc.pzpos(), 2));
      auto fullmomentumNegDaugh = TMath::Sqrt(TMath::Power(casc.pxneg(), 2) + TMath::Power(casc.pyneg(), 2) + TMath::Power(casc.pzneg(), 2));
      auto fullmomentumBachelor = TMath::Sqrt(TMath::Power(casc.pxbach(), 2) + TMath::Power(casc.pybach(), 2) + TMath::Power(casc.pzbach(), 2));

      float cut = etaDauCut;
      histos.fill(HIST("hCutValue"), counter + 1, cut);
      if (TMath::Abs(poseta) > etaDauCut || TMath::Abs(negeta) > etaDauCut || TMath::Abs(bacheta) > etaDauCut)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      if (doCascadeCosPaCut) {
        if (!IsCosPAAccepted(casc, coll.posX(), coll.posY(), coll.posZ(), doPtDepCosPaCut, true))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      cut = dcaV0ToPV;
      histos.fill(HIST("hCutValue"), 17, cut);
      if (doDCAV0ToPVCut) {
        if (TMath::Abs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0ToPV)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      if (casc.sign() < 0) {
        if (doFillNsigmaTPCHistProton)
          histos.fill(HIST("hNsigmaProton"), posExtra.tpcNSigmaPr(), fullMomentumPosDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistV0Pion)
          histos.fill(HIST("hNsigmaPionNeg"), negExtra.tpcNSigmaPi(), fullmomentumNegDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && isXi)
          histos.fill(HIST("hNsigmaPionNegBach"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && !isXi)
          histos.fill(HIST("hNsigmaKaon"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());

      } else if (casc.sign() > 0) {
        if (doFillNsigmaTPCHistV0Pion)
          histos.fill(HIST("hNsigmaPionPos"), posExtra.tpcNSigmaPi(), fullMomentumPosDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistProton)
          histos.fill(HIST("hNsigmaProtonNeg"), negExtra.tpcNSigmaPr(), fullmomentumNegDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && isXi)
          histos.fill(HIST("hNsigmaPionPosBach"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && !isXi)
          histos.fill(HIST("hNsigmaKaon"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());
      }

      if (casc.sign() < 0) {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || TMath::Abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      } else if (casc.sign() > 0) {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi || TMath::Abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      }

      if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      bool kHasTOF = (posExtra.hasTOF() || negExtra.hasTOF() || bachExtra.hasTOF());
      bool kHasITS = (posExtra.hasITS() || negExtra.hasITS() || bachExtra.hasITS());
      if (dooobrej == 1) {
        if (!kHasTOF && !kHasITS)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else if (dooobrej == 2) {
        if (!kHasTOF && (casc.pt() > ptthrtof))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      float cascpos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctau = -10;

      cut = ctau;
      histos.fill(HIST("hCutValue"), 22, cut);

      if (posExtra.hasTOF()) {

        if (doNTOFSigmaProtonCut && casc.sign() < 0) {
          histos.fill(HIST("hNsigmaTOFProton"), casc.tofNSigmaXiLaPr(), fullMomentumPosDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPr()) > nsigmatofPr && fullMomentumPosDaugh > 0.6)
            continue;
        }
        if (doNTOFSigmaV0PionCut && casc.sign() > 0) {
          histos.fill(HIST("hNsigmaTOFV0Pion"), casc.tofNSigmaXiLaPi(), fullMomentumPosDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPi()) > nsigmatofPion)
            continue;
        }
      }

      if (negExtra.hasTOF()) {

        if (doNTOFSigmaProtonCut && casc.sign() > 0) {
          histos.fill(HIST("hNsigmaTOFProton"), casc.tofNSigmaXiLaPr(), fullmomentumNegDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPr()) > nsigmatofPr && fullmomentumNegDaugh > 0.6)
            continue;
        }
        if (doNTOFSigmaV0PionCut && casc.sign() < 0) {
          histos.fill(HIST("hNsigmaTOFV0Pion"), casc.tofNSigmaXiLaPi(), fullmomentumNegDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPi()) > nsigmatofPion)
            continue;
        }
      }

      if (isXi) {

        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }

        if (bachExtra.hasTOF() && doNTOFSigmaBachelorCut) {
          histos.fill(HIST("hNsigmaTOFBachelorPion"), casc.tofNSigmaXiPi(), fullmomentumBachelor, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiPi()) > nsigmatofBachPion)
            continue;
        }

        ctau = pdgDB->Mass(3312) * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
        if (doProperLifeTimeCut) {
          if (ctau > proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      } else {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaKa()) > nsigmatpcKa)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }

        if (bachExtra.hasTOF() && doNTOFSigmaBachelorCut) {
          histos.fill(HIST("hNsigmaTOFBachelorKaon"), casc.tofNSigmaOmKa(), TMath::Sqrt(TMath::Power(casc.pxbach(), 2) + TMath::Power(casc.pybach(), 2) + TMath::Power(casc.pzbach(), 2)), coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaOmKa()) > nsigmatofBachKaon)
            continue;
        }

        ctau = pdgDB->Mass(3334) * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
        if (doProperLifeTimeCut) {
          if (ctau > proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      }
      if (isPositive)
        histos.fill(HIST("InvMassAfterSel/h") + HIST(charge[0]) + HIST("Cascade"), casc.pt(), invmass, coll.centFT0C());
      if (isNegative)
        histos.fill(HIST("InvMassAfterSel/h") + HIST(charge[1]) + HIST("Cascade"), casc.pt(), invmass, coll.centFT0C());

      if (doOccupancyCheck) {
        static_for<0, 9>([&](auto i) {
          constexpr int index = i.value;
          if (coll.centFT0C() < centralityIntervals[index + 1] && coll.centFT0C() > centralityIntervals[index]) {
            if (isPositive)
              histos.fill(HIST("InvMassAfterSelCent") + HIST(Index[index]) + HIST("/h") + HIST(charge[0]) + HIST("Cascade"), casc.pt(), invmass, coll.trackOccupancyInTimeRange());
            if (isNegative)
              histos.fill(HIST("InvMassAfterSelCent") + HIST(Index[index]) + HIST("/h") + HIST(charge[1]) + HIST("Cascade"), casc.pt(), invmass, coll.trackOccupancyInTimeRange());
          }
        });
      }

      float dcaMesonToPV = -10;
      float dcaBaryonToPV = -10;
      if (isPositive) {
        dcaMesonToPV = casc.dcanegtopv();
        dcaBaryonToPV = casc.dcapostopv();
      }
      if (isNegative) {
        dcaBaryonToPV = casc.dcanegtopv();
        dcaMesonToPV = casc.dcapostopv();
      }
      double selections[] = {casc.bachBaryonDCAxyToPV(),
                             TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())),
                             casc.v0radius(),
                             casc.cascradius(),
                             casc.dcaV0daughters(),
                             casc.dcacascdaughters(),
                             TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())),
                             casc.dcabachtopv(),
                             dcaMesonToPV,
                             dcaBaryonToPV,
                             ctau};
      bool selectionToBeTested[] = {doBachelorBaryonCut, doDCAV0ToPVCut, doV0RadiusCut, doCascadeRadiusCut, doDCAV0DauCut, doDCACascadeDauCut, doV0CosPaCut, doCascadeCosPaCut, doDCAdauToPVCut, doDCAdauToPVCut, doDCAdauToPVCut, doProperLifeTimeCut};

      if (doPtDepCutStudy) {
        static_for<0, 10>([&](auto i) {
          constexpr int index = i.value;
          if (!selectionToBeTested[index]) {
            if (isPositive)
              histos.fill(HIST("PtDepCutStudy/h") + HIST(charge[0]) + HIST(selectionNames[index]), casc.pt(), invmass, selections[index]);
            if (isNegative)
              histos.fill(HIST("PtDepCutStudy/h") + HIST(charge[1]) + HIST(selectionNames[index]), casc.pt(), invmass, selections[index]);
          }
        });
      }
    }
  }
  void processCascadesMCrec(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>::iterator const& coll, cascMCCandidates const& Cascades, dauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&) //, , , soa::Join<aod::MotherMCParts const&, aod::StraMCCollMults> const& /*mccollisions*/, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
                                                                                                                                                                                                                                                     // soa::Join<aod::CascCollRefs, aod::CascMCCollRefs, aod::CascCores, aod::CascMCCores, aod::CascExtras, aod::CascBBs, aod::CascTOFNSigmas> const& Cascades, soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs> const&)
  {
    if (!IsEventAccepted(coll, true))
      return;

    for (auto& casc : Cascades) {
      float mass = -1;
      if (isXi)
        mass = casc.mXi();
      else
        mass = casc.mOmega();

      int counter = -1;
      histos.fill(HIST("hCandidate"), ++counter);

      bool isNegative = false;
      bool isPositive = false;
      if (casc.sign() > 0)
        isPositive = true;
      if (casc.sign() < 0)
        isNegative = true;
      if (!isNegative && !isPositive)
        continue;
      // To have trace of how it was before selections
      if (doBefSelCheck) {
        if (isPositive)
          histos.fill(HIST("InvMassBefSel/h") + HIST(charge[0]) + HIST("Cascade"), casc.pt(), mass, coll.centFT0C());
        if (isNegative)
          histos.fill(HIST("InvMassBefSel/h") + HIST(charge[1]) + HIST("Cascade"), casc.pt(), mass, coll.centFT0C());
      }

      if (!casc.has_cascMCCore())
        continue;
      auto cascMC = casc.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

      bool isTrueMCCascade = false;
      if (cascMC.isPhysicalPrimary() && ((isXi && std::abs(cascMC.pdgCode()) == 3312) || (!isXi && std::abs(cascMC.pdgCode()) == 3334)))
        isTrueMCCascade = true;

      float ptmc = RecoDecay::sqrtSumOfSquares(cascMC.pxMC(), cascMC.pyMC());

      if (!IsCascadeCandidateAccepted(casc, counter, coll.centFT0C()))
        continue;
      counter += 13;

      auto negExtra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto bachExtra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      auto poseta = RecoDecay::eta(std::array{casc.pxpos(), casc.pypos(), casc.pzpos()});
      auto negeta = RecoDecay::eta(std::array{casc.pxneg(), casc.pyneg(), casc.pzneg()});
      auto bacheta = RecoDecay::eta(std::array{casc.pxbach(), casc.pybach(), casc.pzbach()});

      auto fullMomentumPosDaugh = TMath::Sqrt(TMath::Power(casc.pxpos(), 2) + TMath::Power(casc.pypos(), 2) + TMath::Power(casc.pzpos(), 2));
      auto fullmomentumNegDaugh = TMath::Sqrt(TMath::Power(casc.pxneg(), 2) + TMath::Power(casc.pyneg(), 2) + TMath::Power(casc.pzneg(), 2));
      auto fullmomentumBachelor = TMath::Sqrt(TMath::Power(casc.pxbach(), 2) + TMath::Power(casc.pybach(), 2) + TMath::Power(casc.pzbach(), 2));

      if (TMath::Abs(poseta) > etaDauCut || TMath::Abs(negeta) > etaDauCut || TMath::Abs(bacheta) > etaDauCut)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      if (doCascadeCosPaCut) {
        if (!IsCosPAAccepted(casc, coll.posX(), coll.posY(), coll.posZ(), doPtDepCosPaCut, true))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      if (doDCAV0ToPVCut) {
        if (TMath::Abs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0ToPV)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      if (isNegative) {
        if (doFillNsigmaTPCHistProton)
          histos.fill(HIST("hNsigmaProton"), posExtra.tpcNSigmaPr(), fullMomentumPosDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistV0Pion)
          histos.fill(HIST("hNsigmaPionNeg"), negExtra.tpcNSigmaPi(), fullmomentumNegDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && isXi)
          histos.fill(HIST("hNsigmaPionNegBach"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && !isXi)
          histos.fill(HIST("hNsigmaKaon"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());

        if (doNTPCSigmaCut) {
          if (TMath::Abs(posExtra.tpcNSigmaPr()) > nsigmatpcPr || TMath::Abs(negExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      }
      if (isPositive) {
        if (doFillNsigmaTPCHistV0Pion)
          histos.fill(HIST("hNsigmaPionPos"), posExtra.tpcNSigmaPi(), fullMomentumPosDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistProton)
          histos.fill(HIST("hNsigmaProtonNeg"), negExtra.tpcNSigmaPr(), fullmomentumNegDaugh, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && isXi)
          histos.fill(HIST("hNsigmaPionPosBach"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());
        if (doFillNsigmaTPCHistPionBach && !isXi)
          histos.fill(HIST("hNsigmaKaon"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, coll.centFT0C());

        if (doNTPCSigmaCut) {
          if (TMath::Abs(posExtra.tpcNSigmaPi()) > nsigmatpcPi || TMath::Abs(negExtra.tpcNSigmaPr()) > nsigmatpcPr)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      }

      if (posExtra.tpcCrossedRows() < mintpccrrows || negExtra.tpcCrossedRows() < mintpccrrows || bachExtra.tpcCrossedRows() < mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      bool kHasTOF = (posExtra.hasTOF() || negExtra.hasTOF() || bachExtra.hasTOF());
      bool kHasITS = (posExtra.hasITS() || negExtra.hasITS() || bachExtra.hasITS());
      if (dooobrej == 1) {
        if (!kHasTOF && !kHasITS)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else if (dooobrej == 2) {
        if (!kHasTOF && (casc.pt() > ptthrtof))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      float cascpos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctau = -10;

      if (posExtra.hasTOF()) {
        if (doNTOFSigmaProtonCut && isNegative) {
          histos.fill(HIST("hNsigmaTOFProton"), casc.tofNSigmaXiLaPr(), fullMomentumPosDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPr()) > nsigmatofPr && fullMomentumPosDaugh > 0.6)
            continue;
        }
        if (doNTOFSigmaV0PionCut && isPositive) {
          histos.fill(HIST("hNsigmaTOFV0Pion"), casc.tofNSigmaXiLaPi(), fullMomentumPosDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPi()) > nsigmatofPion)
            continue;
        }
      }

      if (negExtra.hasTOF()) {
        if (doNTOFSigmaProtonCut && isPositive) {
          histos.fill(HIST("hNsigmaTOFProton"), casc.tofNSigmaXiLaPr(), fullmomentumNegDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPr()) > nsigmatofPr && fullmomentumNegDaugh > 0.6)
            continue;
        }
        if (doNTOFSigmaV0PionCut && isNegative) {
          histos.fill(HIST("hNsigmaTOFV0Pion"), casc.tofNSigmaXiLaPi(), fullmomentumNegDaugh, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiLaPi()) > nsigmatofPion)
            continue;
        }
      }

      if (isXi) {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaPi()) > nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }

        if (bachExtra.hasTOF() && doNTOFSigmaBachelorCut) {
          histos.fill(HIST("hNsigmaTOFBachelorPion"), casc.tofNSigmaXiPi(), fullmomentumBachelor, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaXiPi()) > nsigmatofBachPion)
            continue;
        }

        ctau = pdgDB->Mass(3312) * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
      } else {
        if (doNTPCSigmaCut) {
          if (TMath::Abs(bachExtra.tpcNSigmaKa()) > nsigmatpcKa)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }

        if (bachExtra.hasTOF() && doNTOFSigmaBachelorCut) {
          histos.fill(HIST("hNsigmaTOFBachelorKaon"), casc.tofNSigmaOmKa(), fullmomentumBachelor, coll.centFT0C());
          if (TMath::Abs(casc.tofNSigmaOmKa()) > nsigmatofBachKaon)
            continue;
        }

        ctau = pdgDB->Mass(3334) * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
      }

      if (doProperLifeTimeCut) {
        if (ctau > proplifetime)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      if (isPositive)
        histos.fill(HIST("InvMassAfterSel/h") + HIST(charge[0]) + HIST("Cascade"), casc.pt(), mass, coll.centFT0C());
      if (isNegative)
        histos.fill(HIST("InvMassAfterSel/h") + HIST(charge[1]) + HIST("Cascade"), casc.pt(), mass, coll.centFT0C());
      if (isTrueMCCascade) {
        if (isPositive)
          histos.fill(HIST("InvMassAfterSelMCrecTruth/h") + HIST(charge[0]) + HIST("Cascade"), ptmc, mass, coll.centFT0C());
        if (isNegative)
          histos.fill(HIST("InvMassAfterSelMCrecTruth/h") + HIST(charge[1]) + HIST("Cascade"), ptmc, mass, coll.centFT0C());
      }
      if (doOccupancyCheck) {
        static_for<0, 9>([&](auto i) {
          constexpr int index = i.value;
          if (coll.centFT0C() < centralityIntervals[index + 1] && coll.centFT0C() > centralityIntervals[index]) {
            if (isPositive)
              histos.fill(HIST("InvMassAfterSelCent") + HIST(Index[index]) + HIST("/h") + HIST(charge[0]) + HIST("Cascade"), casc.pt(), mass, coll.trackOccupancyInTimeRange());
            if (isNegative)
              histos.fill(HIST("InvMassAfterSelCent") + HIST(Index[index]) + HIST("/h") + HIST(charge[1]) + HIST("Cascade"), casc.pt(), mass, coll.trackOccupancyInTimeRange());
          }
        });
      }
      float dcaMesonToPV = -10;
      float dcaBaryonToPV = -10;
      if (isPositive) {
        dcaMesonToPV = casc.dcanegtopv();
        dcaBaryonToPV = casc.dcapostopv();
      }
      if (isNegative) {
        dcaBaryonToPV = casc.dcanegtopv();
        dcaMesonToPV = casc.dcapostopv();
      }
      double selections[] = {casc.bachBaryonDCAxyToPV(),
                             TMath::Abs(casc.dcav0topv(casc.x(), casc.y(), casc.z())),
                             casc.v0radius(),
                             casc.cascradius(),
                             casc.dcaV0daughters(),
                             casc.dcacascdaughters(),
                             TMath::ACos(casc.v0cosPA(casc.x(), casc.y(), casc.z())),
                             casc.dcabachtopv(),
                             dcaMesonToPV,
                             dcaBaryonToPV,
                             ctau};
      bool selectionToBeTested[] = {doBachelorBaryonCut, doDCAV0ToPVCut, doV0RadiusCut, doCascadeRadiusCut, doDCAV0DauCut, doDCACascadeDauCut, doV0CosPaCut, doCascadeCosPaCut, doDCAdauToPVCut, doDCAdauToPVCut, doDCAdauToPVCut, doProperLifeTimeCut};

      if (doPtDepCutStudy) {
        static_for<0, 10>([&](auto i) {
          constexpr int index = i.value;
          if (!selectionToBeTested[index]) {
            if (isPositive)
              histos.fill(HIST("PtDepCutStudy/h") + HIST(charge[0]) + HIST(selectionNames[index]), casc.pt(), mass, selections[index]);
            if (isNegative)
              histos.fill(HIST("PtDepCutStudy/h") + HIST(charge[1]) + HIST(selectionNames[index]), casc.pt(), mass, selections[index]);
            if (isTrueMCCascade) {
              if (isPositive)
                histos.fill(HIST("PtDepCutStudyMCTruth/h") + HIST(charge[0]) + HIST(selectionNames[index]), ptmc, mass, selections[index]);
              if (isNegative)
                histos.fill(HIST("PtDepCutStudyMCTruth/h") + HIST(charge[1]) + HIST(selectionNames[index]), ptmc, mass, selections[index]);
            }
          }
        });
      }
    }
  }
  void processCascadesMCforEff(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const& Cascades, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels> const& collisions)
  {

    std::vector<int> listBestCollisionIdx = fillGenEventHist(mcCollisions, collisions);

    for (auto const& cascMC : Cascades) {
      if (!cascMC.has_straMCCollision())
        continue;

      if (!cascMC.isPhysicalPrimary())
        continue;

      float ptmc = RecoDecay::sqrtSumOfSquares(cascMC.pxMC(), cascMC.pyMC());
      float ymc = 1e3;
      if (TMath::Abs(cascMC.pdgCode()) == 3312)
        ymc = RecoDecay::y(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}, o2::constants::physics::MassXiMinus);
      else if (TMath::Abs(cascMC.pdgCode()) == 3334)
        ymc = RecoDecay::y(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}, o2::constants::physics::MassOmegaMinus);

      if (TMath::Abs(ymc) > rapCut)
        continue;

      auto mcCollision = cascMC.straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
      float centrality = 100.5f;
      if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
        auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
        centrality = collision.centFT0C();
      }

      if (cascMC.pdgCode() == 3312) {
        histos.fill(HIST("h2dGenXiMinus"), centrality, ptmc);
        histos.fill(HIST("h2dGenXiMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == -3312) {
        histos.fill(HIST("h2dGenXiPlus"), centrality, ptmc);
        histos.fill(HIST("h2dGenXiPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == 3334) {
        histos.fill(HIST("h2dGenOmegaMinus"), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == -3334) {
        histos.fill(HIST("h2dGenOmegaPlus"), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    }
  }
  // ______________________________________________________
  // Simulated processing
  // Fill event information (for event loss estimation) and return the index to the recoed collision associated to a given MC collision.
  std::vector<int> fillGenEventHist(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels> const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 0.5 /* all gen. events*/);

      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      // Check if there is at least one of the reconstructed collisions associated to this MC collision
      // If so, we consider it
      bool atLeastOne = false;
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      float centrality = 100.5f;
      int nCollisions = 0;
      for (auto const& collision : groupedCollisions) {
        if (!IsEventAccepted(collision, false)) {
          continue;
        }

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          bestCollisionIndex = collision.globalIndex();
          centrality = collision.centFT0C();
        }
        nCollisions++;

        atLeastOne = true;
      }
      listBestCollisionIdx[mcCollision.globalIndex()] = bestCollisionIndex;

      histos.fill(HIST("hCentralityVsNcoll_beforeEvSel"), centrality, groupedCollisions.size() + 0.5);
      histos.fill(HIST("hCentralityVsNcoll_afterEvSel"), centrality, nCollisions + 0.5);

      histos.fill(HIST("hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta05());

      if (atLeastOne) {
        histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 1.5 /* at least 1 rec. event*/);
      }
    }
    return listBestCollisionIdx;
  }
  PROCESS_SWITCH(derivedCascadeAnalysis, processCascades, "cascade analysis, run3 data ", true);
  PROCESS_SWITCH(derivedCascadeAnalysis, processCascadesMCrec, "cascade analysis, run3 rec MC", false);
  PROCESS_SWITCH(derivedCascadeAnalysis, processCascadesMCforEff, "cascade analysis, run3 rec MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<derivedCascadeAnalysis>(cfgc)};
}
