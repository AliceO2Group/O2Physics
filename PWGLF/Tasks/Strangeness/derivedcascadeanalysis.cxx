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
/// \file derivedcascadeanalysis.cxx
/// \brief Tasks processing derived data for Cascade analysis in PbPb collisions
/// \author Lucia Anna Tarasovicova (lucia.anna.husova@cern.ch)

#include "MetadataHelper.h"

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "ReconstructionDataFormats/Track.h"

#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <string>
#include <vector>

// constants
const float ctauxiPDG = 4.91;     // from PDG
const float ctauomegaPDG = 2.461; // from PDG

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using CascMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascBBs, aod::CascCoreMCLabels>;

struct Derivedcascadeanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 25000.0f}, "occupancy axis"};
  ConfigurableAxis axisOccupancyFt0{"axisOccupancyFt0", {VARIABLE_WIDTH, 0.0f, 2500.0f, 5000.0f, 7500.0f, 10000.0f, 15000.0f, 20000.0f, 30000.0f, 45000.0f, 60000.0f, 80000.0f, 100000.0f, 250000.0f}, "occupancy axis for the FT0 definition"};
  ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles in |y| < 0.5"};
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};
  ConfigurableAxis axisMass{"axisMass", {200, 1.222f, 1.422f}, "range of invariant mass, in case of omega take 1.572f, 1.772f"};
  ConfigurableAxis axisIR{"axisIR", {510, -1, 50}, "Binning for the interaction rate (kHz)"};

  Configurable<bool> isXi{"isXi", 1, "Apply cuts for Xi identification"};
  Configurable<bool> ispO{"ispO", 0, "Analyse p--O collisions"};
  Configurable<bool> useCentralityFT0M{"useCentralityFT0M", 0, "If true, use centFT0M"};
  Configurable<bool> useCentralityFT0A{"useCentralityFT0A", 0, "If true, use centFT0A"};
  Configurable<bool> useCentralityFT0Cvar1{"useCentralityFT0Cvar1", 0, "If true, use centFT0FT0Cvar1"};
  Configurable<bool> useTrackOccupancyDef{"useTrackOccupancyDef", true, "Use occupancy definition based on the tracks"};
  Configurable<bool> useFT0OccupancyDef{"useFT0OccupancyDef", false, "se occupancy definition based on the FT0 signals"};
  Configurable<float> minPt{"minPt", 0.0f, "minPt"};
  Configurable<int> nPtBinsForNsigmaTPC{"nPtBinsForNsigmaTPC", 100, ""};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  struct : ConfigurableGroup {
    std::string prefix = "qa";
    Configurable<bool> doFillNsigmaTPCHistPionBach{"doFillNsigmaTPCHistPionBach", false, ""};
    Configurable<bool> doFillNsigmaTPCHistKaonBach{"doFillNsigmaTPCHistKaonBach", false, ""};
    Configurable<bool> doFillNsigmaTPCHistProton{"doFillNsigmaTPCHistProton", false, ""};
    Configurable<bool> doFillNsigmaTPCHistV0Pion{"doFillNsigmaTPCHistV0Pion", false, ""};
    Configurable<bool> doBefSelCheck{"doBefSelCheck", false, "Fill mass histograms of all candidates before selections"};
    Configurable<bool> doPtDepCutStudy{"doPtDepCutStudy", false, "Fill histogram with a cutting paramer"};
    Configurable<bool> doBefSelEventMultCorr{"doBefSelEventMultCorr", false, "Enable histogram of multiplicity correlation before cuts"};
    Configurable<bool> doOccupancyCheck{"doOccupancyCheck", true, ""};
    Configurable<bool> doIRCheck{"doIRCheck", true, ""};
    Configurable<bool> doITSTPCmatchingCheck{"doITSTPCmatchingCheck", true, "fill histogram for ITS-TPC matching check"};
  } qaFlags;

  struct : ConfigurableGroup {
    std::string prefix = "evSelectionRun3";
    Configurable<bool> doTFeventCut{"doTFeventCut", false, "Enable TF event Cut"};
    Configurable<bool> doITSFrameBorderCut{"doITSFrameBorderCut", false, "Enable ITSFrame event cut"};
    Configurable<bool> doGoodPVFT0EventCut{"doGoodPVFT0EventCut", true, "check for the PV position diffrence when estimated from tracks and FT0"};
    Configurable<bool> doITSTPCvertexEventCut{"doITSTPCvertexEventCut", true, "checks the presence of at least one ITS-TPC track"};
    Configurable<bool> doSameBunchPileUpEventCut{"doSameBunchPileUpEventCut", true, "removes events  associated with the same \"found-by-T0\" bunch crossing"};
    Configurable<bool> doVertexTOFmatch{"doVertexTOFmatch", false, "Checks wherher at least one of vertex contributors is matched to TOF"};
    Configurable<bool> doVertexTRDmatch{"doVertexTRDmatch", false, "Checks wherher at least one of vertex contributors is matched to TRD"};
    Configurable<bool> doTimeRangeStandardCut{"doTimeRangeStandardCut", true, "It rejects a given collision if there are other events nearby in dtime +/- 2 μs, or mult above some threshold in -4..-2 μs"};
    Configurable<bool> doTimeRangeStrictCut{"doTimeRangeStrictCut", false, "It rejects a given collision if there are other events nearby in |dt|< 10 μs"};
    Configurable<bool> doNoCollInRofStrictCut{"doNoCollInRofStrictCut", false, "Enable an evevnt selection which rejects a collision if there are other events within the same ITS ROF"};
    Configurable<bool> doNoCollInRofStandardCut{"doNoCollInRofStandardCut", true, "Enable an evevnt selection which rejects a collision if there are other events within the same ITS ROF with mult above threshold"};
    Configurable<bool> doMultiplicityCorrCut{"doMultiplicityCorrCut", false, "Enable multiplicity vs centrality correlation cut"};
    Configurable<bool> doITSallLayersCut{"doITSallLayersCut", false, "Enable event selection which rejects collisions when ITS was rebooting."};
    Configurable<int> minOccupancy{"minOccupancy", -1, "Minimal occupancy"};
    Configurable<int> maxOccupancy{"maxOccupancy", -1, "Maximal occupancy"};
    Configurable<float> minOccupancyFT0{"minOccupancyFT0", -1, "Minimal occupancy"};
    Configurable<float> maxOccupancyFT0{"maxOccupancyFT0", -1, "Maximal occupancy"};
    Configurable<float> globalTracksCorrelpar0Low{"globalTracksCorrelpar0Low", 81, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> globalTracksCorrelpar1Low{"globalTracksCorrelpar1Low", -0.0431016, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> globalTracksCorrelpar2Low{"globalTracksCorrelpar2Low", -6, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> globalTracksCorrelpar0High{"globalTracksCorrelpar0High", 226, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> globalTracksCorrelpar1High{"globalTracksCorrelpar1High", -0.0181686, "[0]*exp([1]*centrality)+[2], mean plus 3*sigma"};
    Configurable<float> globalTracksCorrelpar2High{"globalTracksCorrelpar2High", -22, "[0]*exp([1]*centrality)+[2], mean plus 3*sigma"};
    Configurable<float> pvContribCorrelpar0Low{"pvContribCorrelpar0Low", 152, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> pvContribCorrelpar1Low{"pvContribCorrelpar1Low", -0.0431016, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> pvContribCorrelpar2Low{"pvContribCorrelpar2Low", -15.3776, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> pvContribCorrelpar0High{"pvContribCorrelpar0High", 384.861, "[0]*exp([1]*centrality)+[2], mean minus 3*sigma"};
    Configurable<float> pvContribCorrelpar1High{"pvContribCorrelpar1High", -0.0181686, "[0]*exp([1]*centrality)+[2], mean plus 3*sigma"};
    Configurable<float> pvContribCorrelpar2High{"pvContribCorrelpar2High", -39, "[0]*exp([1]*centrality)+[2], mean plus 3*sigma"};
  } eventSelectionRun3Flags;

  struct : ConfigurableGroup {
    std::string prefix = "evSelectionCommon";
    Configurable<bool> doTriggerSel8EventCut{"doTriggerSel8EventCut", true, "Standard MB event selection"};
    Configurable<bool> doTriggerTVXEventCut{"doTriggerTVXEventCut", false, "Minimal MB event selection, for MC"};
    Configurable<bool> doInel0{"doInel0", true, "Enable INEL > 0 selection"};
    Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
    Configurable<int> centMin{"centMin", 0, "Minimal accepted centrality"};
    Configurable<int> centMax{"centMax", 100, "Maximal accepted centrality"};
    Configurable<bool> doInel0MCGen{"doInel0MCGen", true, "Enable INEL > 0 selection for MC gen events"};
    Configurable<bool> applyZVtxSelOnMCPV{"applyZVtxSelOnMCPV", false, "Enable z vertex cut selection on generated events"};
  } eventSelectionCommonFlags;

  struct : ConfigurableGroup {
    std::string prefix = "evSelectionRun2";
    Configurable<bool> doSel7{"doSel7", true, "require sel7 event selection (Run 2 only: event selection decision based on V0A & V0C)"};
    Configurable<bool> doINT7{"doINT7", true, "require INT7 trigger selection (Run 2 only)"};
    Configurable<bool> doIncompleteDAQCut{"doIncompleteDAQCut", true, "reject events with incomplete DAQ (Run 2 only)"};
    Configurable<bool> doConsistentSPDAndTrackVtx{"doConsistentSPDAndTrackVtx", true, "reject events with inconsistent in SPD and Track vertices (Run 2 only)"};
    Configurable<bool> doPileupFromSPDCut{"doPileupFromSPDCut", true, "reject events with pileup according to SPD vertexer (Run 2 only)"};
    Configurable<bool> doV0PFPileupCut{"doV0PFPileupCut", false, "reject events tagged as OOB pileup according to V0 past-future info (Run 2 only)"};
    Configurable<bool> doPileupInMultBinsCut{"doPileupInMultBinsCut", true, "reject events tagged as pileup according to multiplicity-differential pileup checks (Run 2 only)"};
    Configurable<bool> doPileupMVCut{"doPileupMVCut", true, "reject events tagged as pileup according to according to multi-vertexer (Run 2 only)"};
    Configurable<bool> doTPCPileupCut{"doTPCPileupCut", false, "reject events tagged as pileup according to pileup in TPC (Run 2 only)"};
    Configurable<bool> doNoV0MOnVsOffPileup{"doNoV0MOnVsOffPileup", false, "reject events tagged as OOB pileup according to online-vs-offline VOM correlation (Run 2 only)"};
    Configurable<bool> doNoSPDOnVsOffPileup{"doNoSPDOnVsOffPileup", false, "reject events tagged as pileup according to online-vs-offline SPD correlation (Run 2 only)"};
    Configurable<bool> doNoSPDClsVsTklBG{"doNoSPDClsVsTklBG", true, "reject events tagged as beam-gas and pileup according to cluster-vs-tracklet correlation (Run 2 only)"};

    Configurable<bool> useSPDTrackletsCent{"useSPDTrackletsCent", false, "Use SPD tracklets for estimating centrality? If not, use V0M-based centrality (Run 2 only)"};
  } eventSelectionRun2Flags;

  struct : ConfigurableGroup {
    std::string prefix = "candidateSelFlag";
    Configurable<bool> doPtDepCosPaCut{"doPtDepCosPaCut", false, "Enable pt dependent cos PA cut"};
    Configurable<bool> doPtDepCascRadiusCut{"doPtDepCascRadiusCut", false, "Enable pt dependent cascade radius cut"};
    Configurable<bool> doPtDepV0RadiusCut{"doPtDepV0RadiusCut", false, "Enable pt dependent V0 radius cut"};
    Configurable<bool> doPtDepV0CosPaCut{"doPtDepV0CosPaCut", false, "Enable pt dependent cos PA cut of the V0 daughter"};
    Configurable<bool> doPtDepDCAcascDauCut{"doPtDepDCAcascDauCut", false, "Enable pt dependent DCA cascade daughter cut"};
    Configurable<bool> doDCAbachToPVCut{"doDCAbachToPVCut", true, "Enable cut DCA daughter track to PV"};
    Configurable<bool> doDCAmesonToPVCut{"doDCAmesonToPVCut", true, "Enable cut DCA daughter track to PV"};
    Configurable<bool> doDCAbaryonToPVCut{"doDCAbaryonToPVCut", true, "Enable cut DCA daughter track to PV"};
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
    Configurable<int> dooobrej{"dooobrej", 0, "OOB rejection: 0 no selection, 1 = ITS||TOF, 2 = TOF only for pT > ptthrtof"};
    Configurable<bool> doAtLeastOneTrackAB{"doAtLeastOneTrackAB", false, "require that at least one of the daughter tracks is from Afterburner"};
    Configurable<bool> doBachelorITSTracking{"doBachelorITSTracking", false, "require that the bachelor track is from the ITS tracking"};
    Configurable<bool> doAllTracksMinITSClusters{"doAllTracksMinITSClusters", false, "require that all daughter tracks have minimal ITS hits"};
  } candidateSelectionFlags;

  struct : ConfigurableGroup {
    std::string prefix = "candidateSelection";
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
    Configurable<int> mintpccrrows{"mintpccrrows", 70, "min N TPC crossed rows"};
    Configurable<float> cosPApar0{"cosPApar0", 0.2, "const par for pt dep cosPA cut"};
    Configurable<float> cosPApar1{"cosPApar1", -0.022, "linear par for pt dep cosPA cut"};
    Configurable<float> parCascRadius0{"parCascRadius0", 1.216159, "const par for pt dep radius cut"};
    Configurable<float> parCascRadius1{"parCascRadius1", 0.064462, "linear par for pt dep radius cut"};
    Configurable<float> parV0Radius0{"parV0Radius0", 2.136381, "const par for pt dep V0 radius cut"};
    Configurable<float> parV0Radius1{"parV0Radius1", 0.437074, "linear par for pt dep V0 radius cut"};
    Configurable<float> dcaCacsDauPar0{"dcaCacsDauPar0", 0.8, " par for pt dep DCA cascade daughter cut, p_T < 1 GeV/c"};
    Configurable<float> dcaCacsDauPar1{"dcaCacsDauPar1", 0.5, " par for pt dep DCA cascade daughter cut, 1< p_T < 4 GeV/c"};
    Configurable<float> dcaCacsDauPar2{"dcaCacsDauPar2", 0.2, " par for pt dep DCA cascade daughter cut, p_T > 4 GeV/c"};
    Configurable<float> lambdaMassWin{"lambdaMassWin", 0.005, "V0 Mass window limit"};
    Configurable<float> proplifetime{"proplifetime", 3, "ctau/<ctau>"};
    Configurable<float> ptthrtof{"ptthrtof", 2, "Pt threshold for applying only tof oob rejection"};
    Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};
    Configurable<float> masswin{"masswin", 0.05, "Mass window limit"};
    Configurable<float> rapCut{"rapCut", 0.5, "Rapidity acceptance"};
    Configurable<float> minRapCut{"minRapCut", -0.845, "minimal rapidity acceptance in case of p--o"};
    Configurable<float> maxRapCut{"maxRapCut", 0.155, "maximal rapidity acceptance in case of p--o"};
    Configurable<float> etaDauCut{"etaDauCut", 0.8, "Pseudorapidity acceptance of the cascade daughters"};
    Configurable<int> minITSclusters{"minITSclusters", 3, "minimal number of ITS hits for the daughter tracks"};
  } candidateSelectionValues;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  uint16_t selectionCheckMask;
  double selectionCheck;

  static constexpr std::string_view kCentIndex[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
  static constexpr float kCentralityIntervals[11] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
  static constexpr std::string_view kCharge[] = {"Positive", "Negative"};
  static constexpr std::string_view kSelectionNames[] = {"BachelorBaryonDCA", "DCAV0ToPV", "V0Radius", "CascadeRadius", "DCAV0Daughters", "DCACascDaughters", "V0pa", "CascPA", "DCABachelorToPV", "DCAMesonToPV", "DCABaryonToPV", "CascadeProperLifeTime"};

  // For manual sliceBy
  // Preslice<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCentsRun2, aod::StraEvSelsRun2, aod::StraCollLabels>> perMcCollisionRun2 = aod::v0data::straMCCollisionId;

  void init(InitContext const&)
  {
    if ((doprocessCascades || doprocessCascadesMCrec || doprocessCascadesMCforEff) && (doprocessCascadesRun2 || doprocessCascadesMCrecRun2 || doprocessCascadesMCforEffRun2)) {
      LOGF(fatal, "Cannot enable Run2 and Run3 processes at the same time. Please choose one.");
    }

    // setting CCDB service
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    selectionCheck = -1;
    selectionCheckMask = 0;
    if (!candidateSelectionFlags.doBachelorBaryonCut)
      SETBIT(selectionCheckMask, 0);
    if (!candidateSelectionFlags.doDCAV0ToPVCut)
      SETBIT(selectionCheckMask, 1);
    if (!candidateSelectionFlags.doV0RadiusCut)
      SETBIT(selectionCheckMask, 2);
    if (!candidateSelectionFlags.doCascadeRadiusCut)
      SETBIT(selectionCheckMask, 3);
    if (!candidateSelectionFlags.doDCAV0DauCut)
      SETBIT(selectionCheckMask, 4);
    if (!candidateSelectionFlags.doDCACascadeDauCut)
      SETBIT(selectionCheckMask, 5);
    if (!candidateSelectionFlags.doV0CosPaCut)
      SETBIT(selectionCheckMask, 6);
    if (!candidateSelectionFlags.doCascadeCosPaCut)
      SETBIT(selectionCheckMask, 7);
    if (!candidateSelectionFlags.doDCAbachToPVCut)
      SETBIT(selectionCheckMask, 8);
    if (!candidateSelectionFlags.doDCAmesonToPVCut)
      SETBIT(selectionCheckMask, 9);
    if (!candidateSelectionFlags.doDCAbaryonToPVCut)
      SETBIT(selectionCheckMask, 10);
    if (!candidateSelectionFlags.doProperLifeTimeCut)
      SETBIT(selectionCheckMask, 11);

    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});
    histos.add("hEventMultFt0C", "", kTH1F, {{500, 0, 5000}});
    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{101, 0, 101}});

    histos.add("hEventSelection", "hEventSelection", kTH1F, {{22, 0, 22}});
    // TODO adjust labels
    TString eventSelLabelRun3[22] = {"all", "sel8", "TVX", "PV_{z}", "cent", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kIsVertexITSTPC", "kIsVertexTOFmatched", "kIsVertexTRDmatched", "kNoITSROFrameBorder", "kNoTimeFrameBorder", "MultCorrCut", "kNoCollInTimeRangeStrict", "kNoCollInTimeRangeStandard", "min Occup", "mxOccup", "kNoCollInRofStrict", "kNoCollInRofStandard", "kIsGoodITSLayersAll", "occupFt0", "-"};
    TString eventSelLabelRun2[22] = {"all", "sel8", "TVX", "PV_{z}", "cent", "sel7", "kINT7", "kNoIncompleteDAQ", "kNoInconsistentVtx", "kNoPileupFromSPD", "kNoV0PFPileup", "kNoPileupInMultBins", "kNoPileupMV", "kNoPileupTPC", "kNoV0MOnVsOfPileup", "kNoSPDOnVsOfPileup", "kNoSPDClsVsTklBG", "INEL0", "-", "-", "-", "-"};
    for (int i = 1; i <= histos.get<TH1>(HIST("hEventSelection"))->GetNbinsX(); i++) {
      if (doprocessCascades || doprocessCascadesMCrec || doprocessCascadesMCforEff) {
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabelRun3[i - 1]);
      } else {
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(i, eventSelLabelRun2[i - 1]);
      }
    }

    histos.add("hOccupancyVsOccupFt0VsCentrality", "", kTH3F, {axisOccupancy, axisOccupancyFt0, {100, 0, 100}});

    histos.add("hNCrossedRowsNegative", "", kTH1F, {{400, -200, 200}});
    histos.add("hNCrossedRowsPositive", "", kTH1F, {{400, -200, 200}});
    histos.add("hNCrossedRowsBachelor", "", kTH1F, {{400, -200, 200}});

    histos.add("hPseudorapPosDaughter", "", kTH1F, {{50, -1, 1}});
    histos.add("hPseudorapNegDaughter", "", kTH1F, {{50, -1, 1}});
    histos.add("hPseudorapBachelor", "", kTH1F, {{50, -1, 1}});

    histos.add("hEventNchCorrelationAfCuts", "hEventNchCorrelationAfCuts", kTH2F, {{5000, 0, 5000}, {5000, 0, 2500}});
    histos.add("hEventPVcontributorsVsCentrality", "hEventPVcontributorsVsCentrality", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
    histos.add("hEventGlobalTracksVsCentrality", "hEventGlobalTracksVsCentrality", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    if (qaFlags.doBefSelEventMultCorr) {
      histos.add("hEventNchCorrelationBefCuts", "hEventNchCorrelationBefCuts", kTH2F, {{5000, 0, 5000}, {2500, 0, 2500}});
      histos.add("hEventPVcontributorsVsCentralityBefCuts", "hEventPVcontributorsVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {5000, 0, 5000}});
      histos.add("hEventGlobalTracksVsCentralityBefCuts", "hEventGlobalTracksVsCentralityBefCuts", kTH2F, {{100, 0, 100}, {2500, 0, 2500}});
    }

    histos.add("hCandidate", "hCandidate", HistType::kTH1D, {{22, -0.5, 21.5}});
    histos.add("hCutValue", "hCutValue", HistType::kTH2D, {{22, -0.5, 21.5}, {300, 0, 3.1}});

    if (qaFlags.doFillNsigmaTPCHistProton) {
      histos.add("hNsigmaProton", "hNsigmaProton", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 6}, {100, 0, 100}});
      histos.add("hNsigmaProtonNeg", "hNsigmaProtonNeg", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 6}, {100, 0, 100}});
    }
    if (qaFlags.doFillNsigmaTPCHistV0Pion) {
      histos.add("hNsigmaPionNeg", "hNsigmaPionNeg", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
      histos.add("hNsigmaPionPos", "hNsigmaPionPos", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
    }
    if (qaFlags.doFillNsigmaTPCHistPionBach) {
      histos.add("hNsigmaPionPosBach", "hNsigmaPionPosBach", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
      histos.add("hNsigmaPionNegBach", "hNsigmaPionNegBach", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 2}, {100, 0, 100}});
    }
    if (qaFlags.doFillNsigmaTPCHistKaonBach)
      histos.add("hNsigmaKaon", "hNsigmaKaon", HistType::kTH3D, {{280, -7, 7}, {nPtBinsForNsigmaTPC, 0, 6}, {100, 0, 100}});

    if (candidateSelectionFlags.doNTOFSigmaProtonCut)
      histos.add("hNsigmaTOFProton", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
    if (candidateSelectionFlags.doNTOFSigmaV0PionCut)
      histos.add("hNsigmaTOFV0Pion", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
    if (candidateSelectionFlags.doNTOFSigmaBachelorCut) {
      if (isXi)
        histos.add("hNsigmaTOFBachelorPion", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
      else
        histos.add("hNsigmaTOFBachelorKaon", "", HistType::kTH3D, {{70, -7, 7}, {100, 0, 10}, {100, 0, 100}});
    }

    TString cutLabel[22] = {"All", "MassWin", "y", "DCACascDau", "DCAV0Dau", "rCasc", "rCascMax", "rV0", "rV0Max", "LambdaMass", "Bach-baryon", "V0CosPA", "CompDecayMass", "DCADauToPV", "EtaDau", "CascCosPA", "DCAV0ToPV", "nSigmaTPCV0Dau", "NTPCrows", "OOBRej", "nSigmaTPCbachelor", "ctau"};
    for (int i = 1; i <= histos.get<TH1>(HIST("hCandidate"))->GetNbinsX(); i++) {
      histos.get<TH1>(HIST("hCandidate"))->GetXaxis()->SetBinLabel(i, cutLabel[i - 1]);
      histos.get<TH2>(HIST("hCutValue"))->GetXaxis()->SetBinLabel(i, cutLabel[i - 1]);
    }

    histos.add("InvMassAfterSel/hNegativeCascade", "hNegativeCascade", HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}});
    histos.add("InvMassAfterSel/hPositiveCascade", "hPositiveCascade", {HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}}});

    if (qaFlags.doOccupancyCheck) {
      histos.add("InvMassAfterSelCent1/hNegativeCascade", "hNegativeCascade", HistType::kTH3F, {axisPt, axisMass, axisOccupancy});
      histos.add("InvMassAfterSelCent1/hPositiveCascade", "hPositiveCascade", HistType::kTH3F, {axisPt, axisMass, axisOccupancy});
      if (doprocessCascadesMCrec || doprocessCascadesMCrecRun2) {
        histos.add("InvMassAfterSelCent1/hNegativeCascadeMCTruth", "hNegativeCascadeMCTruth", HistType::kTH3F, {axisPt, axisMass, axisOccupancy});
        histos.add("InvMassAfterSelCent1/hPositiveCascadeMCTruth", "hPositiveCascadeMCTruth", HistType::kTH3F, {axisPt, axisMass, axisOccupancy});
      }
      if (!qaFlags.doIRCheck) {
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
    }

    if (qaFlags.doIRCheck) {
      histos.add("InvMassAfterSelCent1/hNegativeCascadeIR", "hNegativeCascadeIR", HistType::kTH3F, {axisPt, axisMass, axisIR});
      histos.add("InvMassAfterSelCent1/hPositiveCascadeIR", "hPositiveCascadeIR", HistType::kTH3F, {axisPt, axisMass, axisIR});
      if (doprocessCascadesMCrec || doprocessCascadesMCrecRun2) {
        histos.add("InvMassAfterSelCent1/hNegativeCascadeMCTruthIR", "hNegativeCascadeMCTruthIR", HistType::kTH3F, {axisPt, axisMass, axisIR});
        histos.add("InvMassAfterSelCent1/hPositiveCascadeMCTruthIR", "hPositiveCascadeMCTruthIR", HistType::kTH3F, {axisPt, axisMass, axisIR});
      }
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

    histos.add("hInteractionRate", "hInteractionRate", kTH1F, {axisIR});
    histos.add("hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2F, {{101, 0.0f, 101.0f}, axisIR});

    if (qaFlags.doBefSelCheck)
      histos.addClone("InvMassAfterSel/", "InvMassBefSel/");

    if (doprocessCascadesMCrec || doprocessCascadesMCrecRun2)
      histos.addClone("InvMassAfterSel/", "InvMassAfterSelMCrecTruth/");

    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doProperLifeTimeCut) {
      histos.add("PtDepCutStudy/hNegativeCascadeProperLifeTime", "hNegativeCascadeProperLifeTime", HistType::kTH3F, {axisPt, axisMass, {100, 0, 10}});
      histos.add("PtDepCutStudy/hPositiveCascadeProperLifeTime", "hPositiveCascadeProperLifeTime", {HistType::kTH3F, {axisPt, axisMass, {100, 0, 10}}});
    }
    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doBachelorBaryonCut) {
      histos.add("PtDepCutStudy/hNegativeBachelorBaryonDCA", "hNegativeBachelorBaryonDCA", HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}});
      histos.add("PtDepCutStudy/hPositiveBachelorBaryonDCA", "hPositiveBachelorBaryonDCA", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}}});
    }
    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doDCAV0ToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCAV0ToPV", "hNegativeDCAV0ToPV", HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}});
      histos.add("PtDepCutStudy/hPositiveDCAV0ToPV", "hPositiveDCAV0ToPV", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 1}}});
    }
    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doV0RadiusCut) {
      histos.add("PtDepCutStudy/hNegativeV0Radius", "hNegativeV0Radius", HistType::kTH3F, {axisPt, axisMass, {20, 0, 10}});
      histos.add("PtDepCutStudy/hPositiveV0Radius", "hPositiveV0Radius", {HistType::kTH3F, {axisPt, axisMass, {20, 0, 10}}});
    }
    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doCascadeRadiusCut) {
      histos.add("PtDepCutStudy/hNegativeCascadeRadius", "hNegativeCascadeRadius", HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}});
      histos.add("PtDepCutStudy/hPositiveCascadeRadius", "hPositiveCascadeRadius", {HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}}});
    }

    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doDCAV0DauCut) {
      histos.add("PtDepCutStudy/hNegativeDCAV0Daughters", "hNegativeDCAV0Daughters", HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}});
      histos.add("PtDepCutStudy/hPositiveDCAV0Daughters", "hPositiveDCAV0Daughters", {HistType::kTH3F, {axisPt, axisMass, {50, 0, 5}}});
    }

    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doDCACascadeDauCut) {
      histos.add("PtDepCutStudy/hNegativeDCACascDaughters", "hNegativeDCACascDaughters", {HistType::kTH3F, {axisPt, axisMass, {20, 0, 1}}});
      histos.add("PtDepCutStudy/hPositiveDCACascDaughters", "hPositiveDCACascDaughters", {HistType::kTH3F, {axisPt, axisMass, {20, 0, 1}}});
    }

    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doV0CosPaCut) {
      histos.add("PtDepCutStudy/hNegativeV0pa", "hNegativeV0pa", HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}});
      histos.add("PtDepCutStudy/hPositiveV0pa", "hPositiveV0pa", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}}});
    }
    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doDCAbachToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCABachelorToPV", "hNegativeDCABachelorToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hPositiveDCABachelorToPV", "hPositiveDCABachelorToPV", {HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}}});
    }
    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doDCAmesonToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCAMesonToPV", "hNegativeDCAMesonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hPositiveDCAMesonToPV", "hPositiveDCAMesonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
    }
    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doDCAbaryonToPVCut) {
      histos.add("PtDepCutStudy/hNegativeDCABaryonToPV", "hNegativeDCABaryonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
      histos.add("PtDepCutStudy/hPositiveDCABaryonToPV", "hPositiveDCABaryonToPV", HistType::kTH3F, {axisPt, axisMass, {50, 0, 0.5}});
    }

    if (qaFlags.doPtDepCutStudy && !candidateSelectionFlags.doCascadeCosPaCut) {
      histos.add("PtDepCutStudy/hNegativeCascPA", "hNegativeCascPA", HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}});
      histos.add("PtDepCutStudy/hPositiveCascPA", "hPositiveCascPA", {HistType::kTH3F, {axisPt, axisMass, {40, 0, 0.4}}});
    }
    if (qaFlags.doITSTPCmatchingCheck) {
      histos.add("histITSTPCmatchPosTrack", "", HistType::kTH3F, {axisPt, {101, 0, 101}, {3, 0, 3}});
      histos.add("histITSTPCmatchNegTrack", "", HistType::kTH3F, {axisPt, {101, 0, 101}, {3, 0, 3}});
      histos.add("histITSTPCmatchBachTrack", "", HistType::kTH3F, {axisPt, {101, 0, 101}, {3, 0, 3}});
    }

    if (doprocessCascadesMCrec || doprocessCascadesMCrecRun2) {
      histos.addClone("PtDepCutStudy/", "PtDepCutStudyMCTruth/");
      histos.add("hNegativeCascadePtForEfficiency", "hNegativeCascadePtForEfficiency", HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}});
      histos.add("hPositiveCascadePtForEfficiency", "hPositiveCascadePtForEfficiency", {HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}}});
      histos.add("hNegativeCascadePtForEfficiencyVsNch", "hNegativeCascadePtForEfficiencyVsNch", HistType::kTH3F, {axisPt, axisMass, axisNch});
      histos.add("hPositiveCascadePtForEfficiencyVsNch", "hPositiveCascadePtForEfficiencyVsNch", {HistType::kTH3F, {axisPt, axisMass, axisNch}});
      if (qaFlags.doBefSelCheck) {
        histos.add("hNegativeCascadePtForEfficiencyBefSel", "hNegativeCascadePtForEfficiencyBefSel", HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}});
        histos.add("hPositiveCascadePtForEfficiencyBefSel", "hPositiveCascadePtForEfficiencyBefSel", {HistType::kTH3F, {axisPt, axisMass, {101, 0, 101}}});
        histos.add("hNegativeCascadePtForEfficiencyVsNchBefSel", "hNegativeCascadePtForEfficiencyVsNchBefSel", HistType::kTH3F, {axisPt, axisMass, axisNch});
        histos.add("hPositiveCascadePtForEfficiencyVsNchBefSel", "hPositiveCascadePtForEfficiencyVsNchBefSel", {HistType::kTH3F, {axisPt, axisMass, axisNch}});
      }
    }

    if (doprocessCascadesMCforEff || doprocessCascadesMCforEffRun2) {
      histos.add("hGenEvents", "", HistType::kTH2F, {{axisNch}, {4, 0, 4}});
      histos.add("hCentralityVsMultMC", "", kTH2F, {{101, 0.0f, 101.0f}, axisNch});
      histos.add("hRecMultVsMultMC", "", kTH2F, {axisNch, axisNch});
      histos.add("hGenMultMCFT0C", "", kTH1F, {{500, 0, 5000}});
      histos.add("hGenMCNParticlesEta10", "", kTH1F, {{500, 0, 5000}});
      histos.add("hCentralityVsIRGen", "", kTH2F, {{101, 0.0f, 101.0f}, axisIR});
      histos.add("hMultMCVsCentralityVsIRGen", "", kTH3F, {axisNch, {101, 0.0f, 101.0f}, axisIR});

      histos.add("hCentralityVsNcoll_beforeEvSel", "", kTH2F, {{101, 0.0f, 101.0f}, {50, 0.f, 50.f}});
      histos.add("hCentralityVsNcoll_afterEvSel", "", kTH2F, {{101, 0.0f, 101.0f}, {50, 0.f, 50.f}});

      if (isXi) {
        histos.add("h2dGenXiMinusEta", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenXiMinusEtaPosDaughter", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenXiMinusEtaNegDaughter", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenXiMinusEtaBach", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenXiMinus", "h2dGenXiMinus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenXiPlus", "h2dGenXiPlus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenXiMinusVsNch", "h2dGenXiMinusVsNch", kTH2D, {axisNch, axisPt});
        histos.add("h2dGenXiPlusVsNch", "h2dGenXiPlusVsNch", kTH2D, {axisNch, axisPt});
        histos.add("h2dGenXiMinusVsCentOccupancy", "h2dGenXiMinusVsCentOccupancy", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisOccupancy});
        histos.add("h2dGenXiPlusVsCentOccupancy", "h2dGenXiPlusVsCentOccupancy", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisOccupancy});
        histos.add("h2dGenXiMinusVsNchVsOccupancy", "h2dGenXiMinusVsNchVsOccupancy", kTH3D, {axisPt, axisNch, axisOccupancy});
        histos.add("h2dGenXiPlusVsNchVsOccupancy", "h2dGenXiPlusVsNchVsOccupancy", kTH3D, {axisPt, axisNch, axisOccupancy});
        histos.add("h2dGenXiMinusVsMultMCVsCentrality", "h2dGenXiMinusVsMultMCVsCentrality", kTH3D, {axisNch, {101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenXiPlusVsMultMCVsCentrality", "h2dGenXiPlusVsMultMCVsCentrality", kTH3D, {axisNch, {101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenXiMinusVsCentIR", "h2dGenXiMinusVsCentIR", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisIR});
        histos.add("h2dGenXiPlusVsCentIR", "h2dGenXiPlusVsCentIR", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisIR});
        histos.add("h2dGenXiMinusVsMultMCVsIR", "h2dGenXiMinusVsMultMCVsIR", kTH3D, {axisNch, axisIR, axisPt});
        histos.add("h2dGenXiPlusVsMultMCVsIR", "h2dGenXiPlusVsMultMCVsIR", kTH3D, {axisNch, axisIR, axisPt});
      } else {
        histos.add("h2dGenOmegaMinusEta", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenOmegaMinusEtaPosDaughter", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenOmegaMinusEtaNegDaughter", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenOmegaMinusEtaBach", "", kTH1F, {{30, -2, 2}});
        histos.add("h2dGenOmegaMinus", "h2dGenOmegaMinus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenOmegaPlus", "h2dGenOmegaPlus", kTH2D, {{101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenOmegaMinusVsNch", "h2dGenOmegaMinusVsNch", kTH2D, {axisNch, axisPt});
        histos.add("h2dGenOmegaPlusVsNch", "h2dGenOmegaPlusVsNch", kTH2D, {axisNch, axisPt});
        histos.add("h2dGenOmegaMinusVsCentOccupancy", "h2dGenOmegaMinusVsCentOccupancy", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisOccupancy});
        histos.add("h2dGenOmegaPlusVsCentOccupancy", "h2dGenOmegaPlusVsCentOccupancy", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisOccupancy});
        histos.add("h2dGenOmegaMinusVsNchVsOccupancy", "h2dGenOmegaMinusVsNchVsOccupancy", kTH3D, {axisPt, axisNch, axisOccupancy});
        histos.add("h2dGenOmegaPlusVsNchVsOccupancy", "h2dGenOmegaPlusVsNchVsOccupancy", kTH3D, {axisPt, axisNch, axisOccupancy});
        histos.add("h2dGenOmegaMinusVsMultMCVsCentrality", "h2dGenOmegaMinusVsMultMCVsCentrality", kTH3D, {axisNch, {101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenOmegaPlusVsMultMCVsCentrality", "h2dGenOmegaPlusVsMultMCVsCentrality", kTH3D, {axisNch, {101, 0.0f, 101.0f}, axisPt});
        histos.add("h2dGenOmegaMinusVsCentIR", "h2dGenOmegaMinusVsCentIR", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisIR});
        histos.add("h2dGenOmegaPlusVsCentIR", "h2dGenOmegaPlusVsCentIR", kTH3D, {axisPt, {101, 0.0f, 101.0f}, axisIR});
        histos.add("h2dGenOmegaMinusVsMultMCVsIR", "h2dGenOmegaMinusVsMultMCVsIR", kTH3D, {axisNch, axisIR, axisPt});
        histos.add("h2dGenOmegaPlusVsMultMCVsIR", "h2dGenOmegaPlusVsMultMCVsIR", kTH3D, {axisNch, axisIR, axisPt});
      }
    }
  }
  // Return slicing output
  template <bool run3, typename TCollisions>
  auto getGroupedCollisions(TCollisions const& collisions, int globalIndex)
  {
    if constexpr (run3) { // check if we are in Run 3
      return collisions.sliceBy(perMcCollision, globalIndex);
    } else { // we are in Run2
      return collisions.sliceBy(perMcCollisionRun2, globalIndex);
    }
  }
  template <typename TCascade>
  bool isCosPAAccepted(TCascade casc, float x, float y, float z, bool ptdepcut, bool isCascPa)
  {

    if (ptdepcut) {
      double ptdepCut;
      if (isCascPa)
        ptdepCut = candidateSelectionValues.cosPApar0 + candidateSelectionValues.cosPApar1 * casc.pt();
      else
        ptdepCut = candidateSelectionValues.cosPApar0 + candidateSelectionValues.cosPApar1 * casc.pt();
      if (ptdepCut > 0.3 && casc.pt() < 0.5)
        ptdepCut = 0.3;
      if (ptdepCut < 0.012)
        ptdepCut = 0.012;
      if (isCascPa)
        histos.fill(HIST("hCutValue"), 15, std::cos(ptdepCut));
      if (!isCascPa)
        histos.fill(HIST("hCutValue"), 11, std::cos(ptdepCut));
      if (isCascPa && casc.casccosPA(x, y, z) < std::cos(ptdepCut))
        return false;
      if (!isCascPa && casc.v0cosPA(x, y, z) < std::cos(ptdepCut))
        return false;
    } else {
      float cut = candidateSelectionValues.casccospa;
      float cutV0 = candidateSelectionValues.v0cospa;
      if (isCascPa)
        histos.fill(HIST("hCutValue"), 15, cut);
      if (!isCascPa)
        histos.fill(HIST("hCutValue"), 11, cutV0);
      if (isCascPa && casc.casccosPA(x, y, z) < candidateSelectionValues.casccospa)
        return false;
      if (!isCascPa && casc.v0cosPA(x, y, z) < candidateSelectionValues.v0cospa)
        return false;
    }

    return true;
  }

  template <typename TCollision>
  bool isEventAccepted(TCollision coll, bool fillHists)
  {

    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);
    if (eventSelectionCommonFlags.doTriggerSel8EventCut && !coll.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1.5 /* collisions  after sel8*/);

    if (eventSelectionCommonFlags.doTriggerTVXEventCut && !coll.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2.5 /* collisions  after TVX cut*/);
    if (std::abs(coll.posZ()) > eventSelectionCommonFlags.zVertexCut) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3.5 /* collisions  after sel pvz sel*/);
    float centrality = -10;
    double interactionRate = -10;
    int occupancy = -2;
    float occupancyFT0 = -2;

    if constexpr (requires { coll.centFT0C(); }) {
      centrality = coll.centFT0C();
      if (useCentralityFT0M)
        centrality = coll.centFT0M();
      if (useCentralityFT0A)
        centrality = coll.centFV0A();
      if (useCentralityFT0Cvar1)
        centrality = coll.centFT0CVariant1();

      if (qaFlags.doBefSelEventMultCorr) {
        histos.fill(HIST("hEventNchCorrelationBefCuts"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
        histos.fill(HIST("hEventPVcontributorsVsCentralityBefCuts"), centrality, coll.multNTracksPVeta1());
        histos.fill(HIST("hEventGlobalTracksVsCentralityBefCuts"), centrality, coll.multNTracksGlobal());
      }

      if (centrality > eventSelectionCommonFlags.centMax || centrality < eventSelectionCommonFlags.centMin) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 4.5 /* collisions  after centrality sel*/);

      if (eventSelectionRun3Flags.doSameBunchPileUpEventCut && !coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 5.5 /* Not same Bunch pile up */);

      if (eventSelectionRun3Flags.doGoodPVFT0EventCut && !coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 6.5 /* No large vertexZ difference from tracks and FT0*/);

      if (eventSelectionRun3Flags.doITSTPCvertexEventCut && !coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 7.5 /* At least one ITS-TPC track in the event*/);

      if (eventSelectionRun3Flags.doVertexTOFmatch && !coll.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 8.5 /* At least one of vertex contributors is matched to TOF*/);

      if (eventSelectionRun3Flags.doVertexTRDmatch && !coll.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 9.5 /* At least one of vertex contributors is matched to TRD*/);

      if (eventSelectionRun3Flags.doITSFrameBorderCut && !coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 10.5 /* Not at ITS ROF border */);

      if (eventSelectionRun3Flags.doTFeventCut && !coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 11.5 /* Not at TF border */);

      if (eventSelectionRun3Flags.doMultiplicityCorrCut) {
        if (coll.multNTracksGlobal() < (eventSelectionRun3Flags.globalTracksCorrelpar0Low * std::exp(eventSelectionRun3Flags.globalTracksCorrelpar1Low * centrality) + eventSelectionRun3Flags.globalTracksCorrelpar2Low) || coll.multNTracksGlobal() > (eventSelectionRun3Flags.globalTracksCorrelpar0High * std::exp(eventSelectionRun3Flags.globalTracksCorrelpar1High * centrality) + eventSelectionRun3Flags.globalTracksCorrelpar2High))
          return false;
        if (coll.multNTracksPVeta1() < (eventSelectionRun3Flags.pvContribCorrelpar0Low * std::exp(eventSelectionRun3Flags.pvContribCorrelpar1Low * centrality) + eventSelectionRun3Flags.pvContribCorrelpar2Low) || coll.multNTracksPVeta1() > (eventSelectionRun3Flags.pvContribCorrelpar0High * std::exp(eventSelectionRun3Flags.pvContribCorrelpar1High * centrality) + eventSelectionRun3Flags.pvContribCorrelpar2High))
          return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 12.5 /* Remove outlyers */);

      if (eventSelectionRun3Flags.doTimeRangeStrictCut && !coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 13.5 /* Rejection of events too close in time, dtime +/- 10 μs */);

      if (eventSelectionRun3Flags.doTimeRangeStandardCut && !coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 14.5 /* No other collision within +/- 2 μs, or mult above some threshold in -4..-2 μs */);

      occupancy = coll.trackOccupancyInTimeRange();
      if (eventSelectionRun3Flags.minOccupancy > 0 && occupancy < eventSelectionRun3Flags.minOccupancy) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 15.5 /* Below min occupancy */);
      if (eventSelectionRun3Flags.maxOccupancy > 0 && occupancy > eventSelectionRun3Flags.maxOccupancy) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 16.5 /* Above max occupancy */);

      if (eventSelectionRun3Flags.doNoCollInRofStrictCut && !coll.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 17.5 /*rejects a collision if there are other events within the same ITS ROF*/);

      if (eventSelectionRun3Flags.doNoCollInRofStandardCut && !coll.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 18.5 /*rejects a collision if there are other events within the same ITS ROF above mult threshold*/);

      if (eventSelectionRun3Flags.doITSallLayersCut && !coll.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 19.5 /*rejects collisions if ITS was in rebooting stage*/);

      occupancyFT0 = coll.ft0cOccupancyInTimeRange();
      if (eventSelectionRun3Flags.minOccupancyFT0 > -1 && occupancyFT0 < eventSelectionRun3Flags.minOccupancyFT0) {
        return false;
      }
      if (eventSelectionRun3Flags.maxOccupancyFT0 > -1 && occupancyFT0 > eventSelectionRun3Flags.maxOccupancyFT0) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 20.5 /* Occupancy FT0 selection */);

      interactionRate = rateFetcher.fetch(ccdb.service, coll.timestamp(), coll.runNumber(), irSource) * 1.e-3;
    } else {
      centrality = eventSelectionRun2Flags.useSPDTrackletsCent ? coll.centRun2SPDTracklets() : coll.centRun2V0M();

      if (centrality > eventSelectionCommonFlags.centMax || centrality < eventSelectionCommonFlags.centMin) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 4.5 /* collisions  after centrality sel*/);

      if (eventSelectionRun2Flags.doSel7 && !coll.sel7()) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 5.5 /* sel7 selection */);

      if (eventSelectionRun2Flags.doINT7 && !coll.alias_bit(kINT7)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 6.5 /* INT7-triggered collisions */);
      if (eventSelectionRun2Flags.doIncompleteDAQCut && !coll.selection_bit(o2::aod::evsel::kNoIncompleteDAQ)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 7.5 /* Complete events according to DAQ flags */);

      if (eventSelectionRun2Flags.doConsistentSPDAndTrackVtx && !coll.selection_bit(o2::aod::evsel::kNoInconsistentVtx)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 8.5 /* No inconsistency in SPD and Track vertices */);
      if (eventSelectionRun2Flags.doPileupFromSPDCut && !coll.selection_bit(o2::aod::evsel::kNoPileupFromSPD)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 9.5 /* No pileup according to SPD vertexer */);
      if (eventSelectionRun2Flags.doV0PFPileupCut && !coll.selection_bit(o2::aod::evsel::kNoV0PFPileup)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 10.5 /* No out-of-bunch pileup according to V0 past-future info */);

      if (eventSelectionRun2Flags.doPileupInMultBinsCut && !coll.selection_bit(o2::aod::evsel::kNoPileupInMultBins)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 11.5 /* No pileup according to multiplicity-differential pileup checks */);

      if (eventSelectionRun2Flags.doPileupMVCut && !coll.selection_bit(o2::aod::evsel::kNoPileupMV)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 12.5 /* No pileup according to multi-vertexer */);

      if (eventSelectionRun2Flags.doTPCPileupCut && !coll.selection_bit(o2::aod::evsel::kNoPileupTPC)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 13.5 /* No pileup in TPC */);

      if (eventSelectionRun2Flags.doNoV0MOnVsOffPileup && !coll.selection_bit(o2::aod::evsel::kNoV0MOnVsOfPileup)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 14.5 /* No out-of-bunch pileup according to online-vs-offline VOM correlation */);

      if (eventSelectionRun2Flags.doNoSPDOnVsOffPileup && !coll.selection_bit(o2::aod::evsel::kNoSPDOnVsOfPileup)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 15.5 /* No out-of-bunch pileup according to online-vs-offline SPD correlation */);

      if (eventSelectionRun2Flags.doNoSPDClsVsTklBG && !coll.selection_bit(o2::aod::evsel::kNoSPDClsVsTklBG)) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 16.5 /* No beam-gas according to cluster-vs-tracklet correlation */);
    }
    if (eventSelectionCommonFlags.doInel0 && coll.multNTracksPVeta1() < 1) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 27.5 /* INEL > 0 selection */);

    if (fillHists) {
      histos.fill(HIST("hInteractionRate"), interactionRate);
      histos.fill(HIST("hCentralityVsInteractionRate"), centrality, interactionRate);
      histos.fill(HIST("hOccupancyVsOccupFt0VsCentrality"), occupancy, occupancyFT0, centrality);
      histos.fill(HIST("hEventCentrality"), centrality);
      histos.fill(HIST("hEventVertexZ"), coll.posZ());
      histos.fill(HIST("hEventMultFt0C"), coll.multFT0C());
      histos.fill(HIST("hEventNchCorrelationAfCuts"), coll.multNTracksPVeta1(), coll.multNTracksGlobal());
      histos.fill(HIST("hEventPVcontributorsVsCentrality"), centrality, coll.multNTracksPVeta1());
      histos.fill(HIST("hEventGlobalTracksVsCentrality"), centrality, coll.multNTracksGlobal());
    }

    return true;
  }

  template <typename TCascade>
  bool isCascadeCandidateAccepted(TCascade casc, int counter, float /*centrality*/)
  {
    float cut = candidateSelectionValues.masswin;
    histos.fill(HIST("hCutValue"), 1, cut);
    cut = candidateSelectionValues.rapCut;
    if (ispO)
      cut = candidateSelectionValues.maxRapCut;
    histos.fill(HIST("hCutValue"), 2, cut);
    if (isXi) {
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > candidateSelectionValues.masswin) {
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (ispO && (casc.yXi() < candidateSelectionValues.minRapCut || casc.yXi() > candidateSelectionValues.maxRapCut))
        return false;
      else if (std::abs(casc.yXi()) > candidateSelectionValues.rapCut)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > candidateSelectionValues.masswin) {
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (ispO && (casc.yOmega() < candidateSelectionValues.minRapCut || casc.yOmega() > candidateSelectionValues.maxRapCut))
        return false;
      else if (std::abs(casc.yOmega()) > candidateSelectionValues.rapCut)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    }

    if (candidateSelectionFlags.doDCACascadeDauCut) {
      if (candidateSelectionFlags.doPtDepDCAcascDauCut) {
        float ptDepCut = candidateSelectionValues.dcaCacsDauPar0;
        if (casc.pt() > 1 && casc.pt() < 4)
          ptDepCut = candidateSelectionValues.dcaCacsDauPar1;
        else if (casc.pt() > 4)
          ptDepCut = candidateSelectionValues.dcaCacsDauPar2;
        histos.fill(HIST("hCutValue"), 3, ptDepCut);
        if (casc.dcacascdaughters() > ptDepCut)
          return false;
      } else {
        cut = candidateSelectionValues.dcacascdau;
        histos.fill(HIST("hCutValue"), 3, cut);
        if (casc.dcacascdaughters() > candidateSelectionValues.dcacascdau)
          return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (qaFlags.doPtDepCutStudy)
        selectionCheck = casc.dcacascdaughters();
      ++counter;
    }

    if (candidateSelectionFlags.doDCAV0DauCut) {
      cut = candidateSelectionValues.dcav0dau;
      histos.fill(HIST("hCutValue"), 4, cut);
      if (casc.dcaV0daughters() > candidateSelectionValues.dcav0dau)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (qaFlags.doPtDepCutStudy)
        selectionCheck = casc.dcaV0daughters();
      ++counter;
    }

    if (candidateSelectionFlags.doCascadeRadiusCut) {
      if (candidateSelectionFlags.doPtDepCascRadiusCut) {
        double ptdepminRadius = candidateSelectionValues.parCascRadius0 + candidateSelectionValues.parCascRadius1 * casc.pt();
        histos.fill(HIST("hCutValue"), 5, ptdepminRadius);
        if (casc.cascradius() < ptdepminRadius)
          return false;
      } else {
        cut = candidateSelectionValues.minRadius;
        histos.fill(HIST("hCutValue"), 5, cut);
        if (casc.cascradius() < candidateSelectionValues.minRadius)
          return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);

      if (casc.cascradius() > candidateSelectionValues.maxRadius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (qaFlags.doPtDepCutStudy)
        selectionCheck = casc.cascradius();
      counter += 2;
    }

    if (candidateSelectionFlags.doV0RadiusCut) {
      if (candidateSelectionFlags.doPtDepV0RadiusCut) {
        float cut = candidateSelectionValues.parV0Radius0 + casc.pt() * candidateSelectionValues.parV0Radius1;
        histos.fill(HIST("hCutValue"), 7, cut);
        if (casc.v0radius() < cut)
          return false;
      } else {
        cut = candidateSelectionValues.minV0Radius;
        histos.fill(HIST("hCutValue"), 7, cut);
        if (casc.v0radius() < candidateSelectionValues.minV0Radius)
          return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
      if (casc.v0radius() > candidateSelectionValues.maxV0Radius)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (qaFlags.doPtDepCutStudy)
        selectionCheck = casc.v0radius();
      counter += 2;
    }

    cut = candidateSelectionValues.lambdaMassWin;
    histos.fill(HIST("hCutValue"), 9, cut);
    if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda0) > candidateSelectionValues.lambdaMassWin)
      return false;
    histos.fill(HIST("hCandidate"), ++counter);

    cut = candidateSelectionValues.bachBaryonDCAxyToPV;
    histos.fill(HIST("hCutValue"), 10, cut);
    if (candidateSelectionFlags.doBachelorBaryonCut) {
      if ((casc.bachBaryonCosPA() > candidateSelectionValues.bachBaryonCosPA || std::abs(casc.bachBaryonDCAxyToPV()) < candidateSelectionValues.bachBaryonDCAxyToPV)) { // Bach-baryon selection if required
        return false;
      }
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (qaFlags.doPtDepCutStudy)
        selectionCheck = casc.bachBaryonDCAxyToPV();
      ++counter;
    }

    if (candidateSelectionFlags.doV0CosPaCut) {
      if (!isCosPAAccepted(casc, casc.x(), casc.y(), casc.z(), candidateSelectionFlags.doPtDepV0CosPaCut, false))
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (qaFlags.doPtDepCutStudy)
        selectionCheck = std::acos(casc.v0cosPA(casc.x(), casc.y(), casc.z()));
      ++counter;
    }

    cut = candidateSelectionValues.rejcomp;
    histos.fill(HIST("hCutValue"), 12, cut);
    if (isXi) {
      if (std::abs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < candidateSelectionValues.rejcomp)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) < candidateSelectionValues.rejcomp)
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    }

    cut = candidateSelectionValues.dcaBachToPV;
    histos.fill(HIST("hCutValue"), 13, cut);
    cut = candidateSelectionValues.dcaMesonToPV;
    histos.fill(HIST("hCutValue"), 13, cut);
    cut = candidateSelectionValues.dcaBaryonToPV;
    histos.fill(HIST("hCutValue"), 13, cut);
    if (candidateSelectionFlags.doDCAbachToPVCut || candidateSelectionFlags.doDCAmesonToPVCut || candidateSelectionFlags.doDCAbaryonToPVCut) {
      if (candidateSelectionFlags.doDCAbachToPVCut && std::abs(casc.dcabachtopv()) < candidateSelectionValues.dcaBachToPV)
        return false;
      if (casc.sign() > 0 && ((std::abs(casc.dcanegtopv()) < candidateSelectionValues.dcaBaryonToPV && candidateSelectionFlags.doDCAbaryonToPVCut) || (std::abs(casc.dcapostopv()) < candidateSelectionValues.dcaMesonToPV && candidateSelectionFlags.doDCAmesonToPVCut)))
        return false;
      if (casc.sign() < 0 && ((std::abs(casc.dcapostopv()) < candidateSelectionValues.dcaBaryonToPV && candidateSelectionFlags.doDCAbaryonToPVCut) || (std::abs(casc.dcanegtopv()) < candidateSelectionValues.dcaMesonToPV && candidateSelectionFlags.doDCAmesonToPVCut)))
        return false;
      histos.fill(HIST("hCandidate"), ++counter);
    } else {
      if (qaFlags.doPtDepCutStudy) {
        if (!candidateSelectionFlags.doDCAbachToPVCut)
          selectionCheck = std::abs(casc.dcabachtopv());
        if (!candidateSelectionFlags.doDCAbaryonToPVCut && casc.sign() > 0)
          selectionCheck = std::abs(casc.dcanegtopv());
        if (!candidateSelectionFlags.doDCAbaryonToPVCut && casc.sign() < 0)
          selectionCheck = std::abs(casc.dcapostopv());
        if (!candidateSelectionFlags.doDCAmesonToPVCut && casc.sign() > 0)
          selectionCheck = std::abs(casc.dcapostopv());
        if (!candidateSelectionFlags.doDCAmesonToPVCut && casc.sign() < 0)
          selectionCheck = std::abs(casc.dcanegtopv());
      }
      ++counter;
    }

    return true;
  }
  void fillHistOccupancyCheck(double pt, double mass, float occup, float centrality, bool isPos)
  {
    static_for<0, 9>([&](auto i) {
      constexpr int In = i.value;
      if (centrality < kCentralityIntervals[In + 1] && centrality > kCentralityIntervals[In]) {
        if (isPos)
          histos.fill(HIST("InvMassAfterSelCent") + HIST(kCentIndex[In]) + HIST("/hPositiveCascade"), pt, mass, occup);
        else
          histos.fill(HIST("InvMassAfterSelCent") + HIST(kCentIndex[In]) + HIST("/hNegativeCascade"), pt, mass, occup);
      }
    });
  }
  void fillHistIRCheckData(double pt, double mass, double ir, float centrality, bool isPos)
  {
    static_for<0, 9>([&](auto i) {
      constexpr int In = i.value;
      if (centrality < kCentralityIntervals[In + 1] && centrality > kCentralityIntervals[In]) {
        if (isPos)
          histos.fill(HIST("InvMassAfterSelCent") + HIST(kCentIndex[In]) + HIST("/hPositiveCascadeIR"), pt, mass, ir);
        else
          histos.fill(HIST("InvMassAfterSelCent") + HIST(kCentIndex[In]) + HIST("/hNegativeCascadeIR"), pt, mass, ir);
      }
    });
  }
  void fillHistIRCheckMC(double pt, double mass, double ir, float centrality, bool isPos)
  {
    static_for<0, 9>([&](auto i) {
      constexpr int In = i.value;
      if (centrality < kCentralityIntervals[In + 1] && centrality > kCentralityIntervals[In]) {
        if (isPos)
          histos.fill(HIST("InvMassAfterSelCent") + HIST(kCentIndex[In]) + HIST("/hPositiveCascadeMCTruthIR"), pt, mass, ir);
        else
          histos.fill(HIST("InvMassAfterSelCent") + HIST(kCentIndex[In]) + HIST("/hNegativeCascadeMCTruthIR"), pt, mass, ir);
      }
    });
  }
  void fillHistPtDepSelectionStudy(double pt, double mass, double sel, uint64_t selectionCheckMask, bool isPos)
  {
    static_for<0, 11>([&](auto i) {
      constexpr int In = i.value;
      if (TESTBIT(selectionCheckMask, In)) {
        if (isPos)
          histos.fill(HIST("PtDepCutStudy/hPositive") + HIST(kSelectionNames[In]), pt, mass, sel);
        else
          histos.fill(HIST("PtDepCutStudy/hNegative") + HIST(kSelectionNames[In]), pt, mass, sel);
      }
    });
  }
  void fillHistPtDepSelectionStudyMC(double pt, double mass, double sel, uint64_t selectionCheckMask, bool isPos)
  {
    static_for<0, 11>([&](auto i) {
      constexpr int In = i.value;
      if (TESTBIT(selectionCheckMask, In)) {
        if (isPos)
          histos.fill(HIST("PtDepCutStudyMCTruth/hPositive") + HIST(kSelectionNames[In]), pt, mass, sel);
        else
          histos.fill(HIST("PtDepCutStudyMCTruth/hNegative") + HIST(kSelectionNames[In]), pt, mass, sel);
      }
    });
  }
  template <typename TCollision, typename TCascade>
  void analyseCascades(TCollision const& coll, TCascade const& Cascades) //(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFNSigmas> const& Cascades, DauTracks const&)//(TCollision const& coll, TCascade const& Cascades)
  {

    if (!isEventAccepted(coll, true))
      return;
    float centrality = -1;
    float nChEta1 = -1;
    float occupancy = -2;

    if constexpr (requires { coll.centFT0C(); }) {
      nChEta1 = coll.multNTracksPVeta1();
      centrality = coll.centFT0C();
      if (useCentralityFT0M)
        centrality = coll.centFT0M();
      if (useCentralityFT0A)
        centrality = coll.centFV0A();
      if (useCentralityFT0Cvar1)
        centrality = coll.centFT0CVariant1();
      occupancy = useTrackOccupancyDef ? coll.trackOccupancyInTimeRange() : coll.ft0cOccupancyInTimeRange();
    } else {
      centrality = eventSelectionRun2Flags.useSPDTrackletsCent ? coll.centRun2SPDTracklets() : coll.centRun2V0M();
    }

    for (const auto& casc : Cascades) {

      int counter = -1;
      histos.fill(HIST("hCandidate"), ++counter);

      double recoPt = casc.pt();

      bool isNegative = false;
      bool isPositive = false;
      if (casc.sign() > 0)
        isPositive = true;
      if (casc.sign() < 0)
        isNegative = true;
      if (!isNegative && !isPositive)
        continue;

      double invmass = -10;
      invmass = isXi ? casc.mXi() : casc.mOmega();
      // To have trace of how it was before selections
      if (qaFlags.doBefSelCheck) {
        if (isPositive)
          histos.fill(HIST("InvMassBefSel/hPositiveCascade"), recoPt, invmass, centrality);
        if (isNegative)
          histos.fill(HIST("InvMassBefSel/hNegativeCascade"), recoPt, invmass, centrality);
      }

      bool isTrueMCCascade = false;
      bool isTrueMCCascadeDecay = false;
      bool isCorrectLambdaDecay = false;
      float ptmc = -1;

      // MC part
      if constexpr (requires { casc.has_cascMCCore(); }) {
        if (!casc.has_cascMCCore())
          continue;
        auto cascMC = casc.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();
        ptmc = RecoDecay::sqrtSumOfSquares(cascMC.pxMC(), cascMC.pyMC());

        if (cascMC.isPhysicalPrimary() && ((isXi && std::abs(cascMC.pdgCode()) == PDG_t::kXiMinus) || (!isXi && std::abs(cascMC.pdgCode()) == PDG_t::kOmegaMinus)))
          isTrueMCCascade = true;
        if (isTrueMCCascade && ((isPositive && cascMC.pdgCodePositive() == PDG_t::kPiPlus && cascMC.pdgCodeNegative() == PDG_t::kProtonBar) || (isNegative && cascMC.pdgCodePositive() == PDG_t::kProton && cascMC.pdgCodeNegative() == PDG_t::kPiMinus)))
          isCorrectLambdaDecay = true;
        if (isTrueMCCascade && isCorrectLambdaDecay && ((isXi && std::abs(cascMC.pdgCodeBachelor()) == PDG_t::kPiPlus) || (!isXi && std::abs(cascMC.pdgCodeBachelor()) == PDG_t::kKPlus)))
          isTrueMCCascadeDecay = true;

        if (qaFlags.doBefSelCheck && isTrueMCCascade) {
          if (isPositive) {
            histos.fill(HIST("InvMassBefSel/hPositiveCascadePtForEfficiencyVsNchBefSel"), ptmc, invmass, nChEta1);
            histos.fill(HIST("InvMassBefSel/hPositiveCascadePtForEfficiencyBefSel"), ptmc, invmass, centrality);
          }
          if (isNegative) {
            histos.fill(HIST("InvMassBefSel/hNegativeCascadePtForEfficiencyVsNchBefSel"), ptmc, invmass, nChEta1);
            histos.fill(HIST("InvMassBefSel/hNegativeCascadePtForEfficiencyBefSel"), ptmc, invmass, centrality);
          }
        }
      }

      if (!isCascadeCandidateAccepted(casc, counter, centrality))
        continue;
      counter += 13;

      auto negExtra = casc.template negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto posExtra = casc.template posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();
      auto bachExtra = casc.template bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>>();

      auto poseta = RecoDecay::eta(std::array{casc.pxpos(), casc.pypos(), casc.pzpos()});
      auto negeta = RecoDecay::eta(std::array{casc.pxneg(), casc.pyneg(), casc.pzneg()});
      auto bacheta = RecoDecay::eta(std::array{casc.pxbach(), casc.pybach(), casc.pzbach()});

      auto fullMomentumPosDaugh = std::sqrt(std::pow(casc.pxpos(), 2) + std::pow(casc.pypos(), 2) + std::pow(casc.pzpos(), 2));
      auto fullmomentumNegDaugh = std::sqrt(std::pow(casc.pxneg(), 2) + std::pow(casc.pyneg(), 2) + std::pow(casc.pzneg(), 2));
      auto fullmomentumBachelor = std::sqrt(std::pow(casc.pxbach(), 2) + std::pow(casc.pybach(), 2) + std::pow(casc.pzbach(), 2));

      float cut = candidateSelectionValues.etaDauCut;
      histos.fill(HIST("hCutValue"), counter + 1, cut);
      if (std::abs(poseta) > candidateSelectionValues.etaDauCut || std::abs(negeta) > candidateSelectionValues.etaDauCut || std::abs(bacheta) > candidateSelectionValues.etaDauCut)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      histos.fill(HIST("hPseudorapPosDaughter"), poseta);
      histos.fill(HIST("hPseudorapNegDaughter"), negeta);
      histos.fill(HIST("hPseudorapBachelor"), bacheta);

      if (candidateSelectionFlags.doCascadeCosPaCut) {
        if (!isCosPAAccepted(casc, coll.posX(), coll.posY(), coll.posZ(), candidateSelectionFlags.doPtDepCosPaCut, true))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        if (qaFlags.doPtDepCutStudy)
          selectionCheck = std::acos(casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()));
        ++counter;
      }

      cut = candidateSelectionValues.dcaV0ToPV;
      histos.fill(HIST("hCutValue"), 16, cut);
      if (candidateSelectionFlags.doDCAV0ToPVCut) {
        if (std::abs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < candidateSelectionValues.dcaV0ToPV)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        if (qaFlags.doPtDepCutStudy)
          selectionCheck = std::abs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()));
        ++counter;
      }

      if (isNegative) {
        if (qaFlags.doFillNsigmaTPCHistProton)
          histos.fill(HIST("hNsigmaProton"), posExtra.tpcNSigmaPr(), fullMomentumPosDaugh, centrality);
        if (qaFlags.doFillNsigmaTPCHistV0Pion)
          histos.fill(HIST("hNsigmaPionNeg"), negExtra.tpcNSigmaPi(), fullmomentumNegDaugh, centrality);
        if (qaFlags.doFillNsigmaTPCHistPionBach && isXi)
          histos.fill(HIST("hNsigmaPionNegBach"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, centrality);
        if (qaFlags.doFillNsigmaTPCHistPionBach && !isXi)
          histos.fill(HIST("hNsigmaKaon"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, centrality);

      } else {
        if (qaFlags.doFillNsigmaTPCHistV0Pion)
          histos.fill(HIST("hNsigmaPionPos"), posExtra.tpcNSigmaPi(), fullMomentumPosDaugh, centrality);
        if (qaFlags.doFillNsigmaTPCHistProton)
          histos.fill(HIST("hNsigmaProtonNeg"), negExtra.tpcNSigmaPr(), fullmomentumNegDaugh, centrality);
        if (qaFlags.doFillNsigmaTPCHistPionBach && isXi)
          histos.fill(HIST("hNsigmaPionPosBach"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, centrality);
        if (qaFlags.doFillNsigmaTPCHistPionBach && !isXi)
          histos.fill(HIST("hNsigmaKaon"), bachExtra.tpcNSigmaPi(), fullmomentumBachelor, centrality);
      }

      if (isNegative) {
        if (candidateSelectionFlags.doNTPCSigmaCut) {
          if (std::abs(posExtra.tpcNSigmaPr()) > candidateSelectionValues.nsigmatpcPr || std::abs(negExtra.tpcNSigmaPi()) > candidateSelectionValues.nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      } else {
        if (candidateSelectionFlags.doNTPCSigmaCut) {
          if (std::abs(posExtra.tpcNSigmaPi()) > candidateSelectionValues.nsigmatpcPi || std::abs(negExtra.tpcNSigmaPr()) > candidateSelectionValues.nsigmatpcPr)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }
      }

      histos.fill(HIST("hNCrossedRowsBachelor"), bachExtra.tpcCrossedRows());
      histos.fill(HIST("hNCrossedRowsNegative"), negExtra.tpcCrossedRows());
      histos.fill(HIST("hNCrossedRowsPositive"), posExtra.tpcCrossedRows());

      if (qaFlags.doITSTPCmatchingCheck) {
        auto ptPosDaugh = std::sqrt(std::pow(casc.pxpos(), 2) + std::pow(casc.pypos(), 2));
        auto ptNegDaugh = std::sqrt(std::pow(casc.pxneg(), 2) + std::pow(casc.pyneg(), 2));
        auto ptBachelor = std::sqrt(std::pow(casc.pxbach(), 2) + std::pow(casc.pybach(), 2));

        histos.fill(HIST("histITSTPCmatchPosTrack"), ptPosDaugh, centrality, 0.5);
        histos.fill(HIST("histITSTPCmatchNegTrack"), ptNegDaugh, centrality, 0.5);
        histos.fill(HIST("histITSTPCmatchBachTrack"), ptBachelor, centrality, 0.5);
        if (candidateSelectionFlags.doAllTracksMinITSClusters) {
          if (posExtra.hasITS() && posExtra.itsNCls() >= candidateSelectionValues.minITSclusters)
            histos.fill(HIST("histITSTPCmatchPosTrack"), ptPosDaugh, centrality, 1.5);
          if (negExtra.hasITS() && negExtra.itsNCls() >= candidateSelectionValues.minITSclusters)
            histos.fill(HIST("histITSTPCmatchNegTrack"), ptNegDaugh, centrality, 1.5);
          if (bachExtra.hasITS() && bachExtra.itsNCls() >= candidateSelectionValues.minITSclusters)
            histos.fill(HIST("histITSTPCmatchBachTrack"), ptBachelor, centrality, 1.5);
          if (posExtra.hasITS() && posExtra.itsNCls() >= candidateSelectionValues.minITSclusters && posExtra.hasTPC() && std::abs(posExtra.tpcCrossedRows()) >= candidateSelectionValues.mintpccrrows)
            histos.fill(HIST("histITSTPCmatchPosTrack"), ptPosDaugh, centrality, 2.5);
          if (negExtra.hasITS() && negExtra.itsNCls() >= candidateSelectionValues.minITSclusters && negExtra.hasTPC() && std::abs(negExtra.tpcCrossedRows()) >= candidateSelectionValues.mintpccrrows)
            histos.fill(HIST("histITSTPCmatchNegTrack"), ptNegDaugh, centrality, 2.5);
          if (bachExtra.hasITS() && bachExtra.itsNCls() >= candidateSelectionValues.minITSclusters && bachExtra.hasTPC() && std::abs(bachExtra.tpcCrossedRows()) >= candidateSelectionValues.mintpccrrows)
            histos.fill(HIST("histITSTPCmatchBachTrack"), ptBachelor, centrality, 2.5);
        } else {
          if (posExtra.hasTPC() && std::abs(posExtra.tpcCrossedRows()) >= candidateSelectionValues.mintpccrrows)
            histos.fill(HIST("histITSTPCmatchPosTrack"), ptPosDaugh, centrality, 2.5);
          if (negExtra.hasTPC() && std::abs(negExtra.tpcCrossedRows()) >= candidateSelectionValues.mintpccrrows)
            histos.fill(HIST("histITSTPCmatchNegTrack"), ptNegDaugh, centrality, 2.5);
          if (bachExtra.hasTPC() && std::abs(bachExtra.tpcCrossedRows()) >= candidateSelectionValues.mintpccrrows)
            histos.fill(HIST("histITSTPCmatchBachTrack"), ptBachelor, centrality, 2.5);
        }
      }

      if (std::abs(posExtra.tpcCrossedRows()) < candidateSelectionValues.mintpccrrows || std::abs(negExtra.tpcCrossedRows()) < candidateSelectionValues.mintpccrrows || std::abs(bachExtra.tpcCrossedRows()) < candidateSelectionValues.mintpccrrows)
        continue;
      histos.fill(HIST("hCandidate"), ++counter);

      bool kHasTOF = (posExtra.hasTOF() || negExtra.hasTOF() || bachExtra.hasTOF());
      bool kHasITS = (posExtra.hasITS() || negExtra.hasITS() || bachExtra.hasITS());
      if (candidateSelectionFlags.dooobrej == 1) {
        if (!kHasTOF && !kHasITS)
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else if (candidateSelectionFlags.dooobrej == 2) {
        if (!kHasTOF && (casc.pt() > candidateSelectionValues.ptthrtof))
          continue;
        histos.fill(HIST("hCandidate"), ++counter);
      } else {
        ++counter;
      }

      float cascpos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      float ctau = -10;

      cut = candidateSelectionValues.proplifetime;
      histos.fill(HIST("hCutValue"), 21, cut);

      if (posExtra.hasTOF()) {

        if (candidateSelectionFlags.doNTOFSigmaProtonCut && casc.sign() < 0) {
          histos.fill(HIST("hNsigmaTOFProton"), casc.tofNSigmaXiLaPr(), fullMomentumPosDaugh, centrality);
          if (std::abs(casc.tofNSigmaXiLaPr()) > candidateSelectionValues.nsigmatofPr && fullMomentumPosDaugh > 0.6)
            continue;
        }
        if (candidateSelectionFlags.doNTOFSigmaV0PionCut && casc.sign() > 0) {
          histos.fill(HIST("hNsigmaTOFV0Pion"), casc.tofNSigmaXiLaPi(), fullMomentumPosDaugh, centrality);
          if (std::abs(casc.tofNSigmaXiLaPi()) > candidateSelectionValues.nsigmatofPion)
            continue;
        }
      }

      if (negExtra.hasTOF()) {

        if (candidateSelectionFlags.doNTOFSigmaProtonCut && casc.sign() > 0) {
          histos.fill(HIST("hNsigmaTOFProton"), casc.tofNSigmaXiLaPr(), fullmomentumNegDaugh, centrality);
          if (std::abs(casc.tofNSigmaXiLaPr()) > candidateSelectionValues.nsigmatofPr && fullmomentumNegDaugh > 0.6)
            continue;
        }
        if (candidateSelectionFlags.doNTOFSigmaV0PionCut && casc.sign() < 0) {
          histos.fill(HIST("hNsigmaTOFV0Pion"), casc.tofNSigmaXiLaPi(), fullmomentumNegDaugh, centrality);
          if (std::abs(casc.tofNSigmaXiLaPi()) > candidateSelectionValues.nsigmatofPion)
            continue;
        }
      }

      if (candidateSelectionFlags.doAtLeastOneTrackAB) {
        if (bachExtra.hasITSTracker() && negExtra.hasITSTracker() && posExtra.hasITSTracker())
          continue;
      }
      if (candidateSelectionFlags.doBachelorITSTracking) {
        if (!bachExtra.hasITSTracker())
          continue;
      }
      if (candidateSelectionFlags.doAllTracksMinITSClusters) {
        if (bachExtra.itsNCls() < candidateSelectionValues.minITSclusters || posExtra.itsNCls() < candidateSelectionValues.minITSclusters || negExtra.itsNCls() < candidateSelectionValues.minITSclusters)
          continue;
      }

      if (isXi) {

        if (candidateSelectionFlags.doNTPCSigmaCut) {
          if (std::abs(bachExtra.tpcNSigmaPi()) > candidateSelectionValues.nsigmatpcPi)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }

        if (bachExtra.hasTOF() && candidateSelectionFlags.doNTOFSigmaBachelorCut) {
          histos.fill(HIST("hNsigmaTOFBachelorPion"), casc.tofNSigmaXiPi(), fullmomentumBachelor, centrality);
          if (std::abs(casc.tofNSigmaXiPi()) > candidateSelectionValues.nsigmatofBachPion)
            continue;
        }

        ctau = o2::constants::physics::MassXiMinus * cascpos / ((cascptotmom + 1e-13) * ctauxiPDG);
        if (candidateSelectionFlags.doProperLifeTimeCut) {
          if (ctau > candidateSelectionValues.proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          if (qaFlags.doPtDepCutStudy)
            selectionCheck = ctau;
          ++counter;
        }
      } else {
        if (candidateSelectionFlags.doNTPCSigmaCut) {
          if (std::abs(bachExtra.tpcNSigmaKa()) > candidateSelectionValues.nsigmatpcKa)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          ++counter;
        }

        if (bachExtra.hasTOF() && candidateSelectionFlags.doNTOFSigmaBachelorCut) {
          histos.fill(HIST("hNsigmaTOFBachelorKaon"), casc.tofNSigmaOmKa(), std::sqrt(std::pow(casc.pxbach(), 2) + std::pow(casc.pybach(), 2) + std::pow(casc.pzbach(), 2)), centrality);
          if (std::abs(casc.tofNSigmaOmKa()) > candidateSelectionValues.nsigmatofBachKaon)
            continue;
        }

        ctau = o2::constants::physics::MassOmegaMinus * cascpos / ((cascptotmom + 1e-13) * ctauomegaPDG);
        if (candidateSelectionFlags.doProperLifeTimeCut) {
          if (ctau > candidateSelectionValues.proplifetime)
            continue;
          histos.fill(HIST("hCandidate"), ++counter);
        } else {
          if (qaFlags.doPtDepCutStudy)
            selectionCheck = ctau;
          ++counter;
        }
      }

      if (qaFlags.doOccupancyCheck) {
        fillHistOccupancyCheck(recoPt, invmass, occupancy, centrality, isPositive);
      }
      if (isPositive) {
        histos.fill(HIST("InvMassAfterSel/hPositiveCascade"), recoPt, invmass, centrality);
        if (isTrueMCCascadeDecay) {
          histos.fill(HIST("hPositiveCascadePtForEfficiency"), ptmc, invmass, centrality);
          histos.fill(HIST("hPositiveCascadePtForEfficiencyVsNch"), ptmc, invmass, nChEta1);
        }
        if (isTrueMCCascadeDecay)
          histos.fill(HIST("InvMassAfterSelMCrecTruth/hPositiveCascade"), ptmc, invmass, centrality);
      } else {
        histos.fill(HIST("InvMassAfterSel/hNegativeCascade"), recoPt, invmass, centrality);
        if (isTrueMCCascadeDecay) {
          histos.fill(HIST("hNegativeCascadePtForEfficiency"), ptmc, invmass, centrality);
          histos.fill(HIST("hNegativeCascadePtForEfficiencyVsNch"), ptmc, invmass, nChEta1);
        }
        if (isTrueMCCascadeDecay)
          histos.fill(HIST("InvMassAfterSelMCrecTruth/hNegativeCascade"), ptmc, invmass, centrality);
      }

      if (qaFlags.doIRCheck) {
        double interactionRate = -2;
        if constexpr (requires { coll.trackOccupancyInTimeRange(); }) {
          interactionRate = rateFetcher.fetch(ccdb.service, coll.timestamp(), coll.runNumber(), irSource) * 1.e-3;
          fillHistIRCheckData(recoPt, invmass, interactionRate, centrality, isPositive);
          if (isTrueMCCascadeDecay)
            fillHistIRCheckMC(ptmc, invmass, interactionRate, centrality, isPositive);
        }
      }

      if (qaFlags.doPtDepCutStudy && selectionCheck != -1) {
        fillHistPtDepSelectionStudy(recoPt, invmass, selectionCheck, selectionCheckMask, isPositive);
        if (isTrueMCCascade)
          fillHistPtDepSelectionStudyMC(ptmc, invmass, selectionCheck, selectionCheckMask, isPositive);
      }
    }
  }
  template <bool run3, typename TMCCollisions, typename TCascMCs, typename TCollisions>
  void analyseCascadesMCforEff(TMCCollisions const& mcCollisions, TCascMCs const& Cascades, TCollisions const& collisions)
  {

    std::vector<int> listBestCollisionIdx = fillGenEventHist<run3>(mcCollisions, collisions);

    for (auto const& cascMC : Cascades) {
      if (!cascMC.has_straMCCollision())
        continue;

      if (!cascMC.isPhysicalPrimary())
        continue;

      float ptmc = RecoDecay::sqrtSumOfSquares(cascMC.pxMC(), cascMC.pyMC());
      float ymc = 1e3;
      if (std::abs(cascMC.pdgCode()) == PDG_t::kXiMinus)
        ymc = RecoDecay::y(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}, o2::constants::physics::MassXiMinus);
      else if (std::abs(cascMC.pdgCode()) == PDG_t::kOmegaMinus)
        ymc = RecoDecay::y(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}, o2::constants::physics::MassOmegaMinus);

      if (ispO && (ymc > candidateSelectionValues.maxRapCut || ymc < candidateSelectionValues.minRapCut))
        continue;
      else if (std::abs(ymc) > candidateSelectionValues.rapCut)
        continue;

      auto mcCollision = cascMC.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
      if (eventSelectionCommonFlags.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelectionCommonFlags.zVertexCut) {
        continue;
      }

      if (eventSelectionCommonFlags.doInel0MCGen && mcCollision.multMCNParticlesEta10() < 1) {
        continue;
      }

      float centrality = 100.5f;
      float occupancy = -2;
      float nChEta1 = -1;
      double intRate = -1;
      if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
        auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
        if constexpr (requires { collision.centFT0C(); }) {
          intRate = rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource) * 1.e-3;
          centrality = collision.centFT0C();
          if (useCentralityFT0M)
            centrality = collision.centFT0M();
          if (useCentralityFT0A)
            centrality = collision.centFV0A();
          if (useCentralityFT0Cvar1)
            centrality = collision.centFT0CVariant1();
          if (useTrackOccupancyDef)
            occupancy = collision.trackOccupancyInTimeRange();
          if (useFT0OccupancyDef)
            occupancy = collision.ft0cOccupancyInTimeRange();
        } else {
          centrality = eventSelectionRun2Flags.useSPDTrackletsCent ? collision.centRun2SPDTracklets() : collision.centRun2V0M();
        }
        nChEta1 = collision.multNTracksPVeta1();
      }

      if (cascMC.pdgCode() == PDG_t::kXiMinus && isXi) {
        histos.fill(HIST("h2dGenXiMinus"), centrality, ptmc);
        histos.fill(HIST("h2dGenXiMinusVsNch"), nChEta1, ptmc);
        histos.fill(HIST("h2dGenXiMinusEta"), RecoDecay::eta(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}));
        histos.fill(HIST("h2dGenXiMinusEtaPosDaughter"), RecoDecay::eta(std::array{cascMC.pxPosMC(), cascMC.pyPosMC(), cascMC.pzPosMC()}));
        histos.fill(HIST("h2dGenXiMinusEtaNegDaughter"), RecoDecay::eta(std::array{cascMC.pxNegMC(), cascMC.pyNegMC(), cascMC.pzNegMC()}));
        histos.fill(HIST("h2dGenXiMinusEtaBach"), RecoDecay::eta(std::array{cascMC.pxBachMC(), cascMC.pyBachMC(), cascMC.pzBachMC()}));
        histos.fill(HIST("h2dGenXiMinusVsMultMCVsCentrality"), mcCollision.multMCNParticlesEta05(), centrality, ptmc);
        histos.fill(HIST("h2dGenXiMinusVsMultMCVsIR"), mcCollision.multMCNParticlesEta05(), intRate, ptmc);
        histos.fill(HIST("h2dGenXiMinusVsCentOccupancy"), ptmc, centrality, occupancy);
        histos.fill(HIST("h2dGenXiMinusVsCentIR"), ptmc, centrality, intRate);
        histos.fill(HIST("h2dGenXiMinusVsNchVsOccupancy"), ptmc, nChEta1, occupancy);
      }
      if (cascMC.pdgCode() == PDG_t::kXiPlusBar && isXi) {
        histos.fill(HIST("h2dGenXiPlus"), centrality, ptmc);
        histos.fill(HIST("h2dGenXiPlusVsNch"), nChEta1, ptmc);
        histos.fill(HIST("h2dGenXiPlusVsMultMCVsCentrality"), mcCollision.multMCNParticlesEta05(), centrality, ptmc);
        histos.fill(HIST("h2dGenXiPlusVsMultMCVsIR"), mcCollision.multMCNParticlesEta05(), intRate, ptmc);
        histos.fill(HIST("h2dGenXiPlusVsCentOccupancy"), ptmc, centrality, occupancy);
        histos.fill(HIST("h2dGenXiPlusVsNchVsOccupancy"), ptmc, nChEta1, occupancy);
        histos.fill(HIST("h2dGenXiPlusVsCentIR"), ptmc, centrality, intRate);
      }
      if (cascMC.pdgCode() == PDG_t::kOmegaMinus && !isXi) {
        histos.fill(HIST("h2dGenOmegaMinus"), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaMinusVsNch"), nChEta1, ptmc);
        histos.fill(HIST("h2dGenOmegaMinusEta"), RecoDecay::eta(std::array{cascMC.pxMC(), cascMC.pyMC(), cascMC.pzMC()}));
        histos.fill(HIST("h2dGenOmegaMinusEtaPosDaughter"), RecoDecay::eta(std::array{cascMC.pxPosMC(), cascMC.pyPosMC(), cascMC.pzPosMC()}));
        histos.fill(HIST("h2dGenOmegaMinusEtaNegDaughter"), RecoDecay::eta(std::array{cascMC.pxNegMC(), cascMC.pyNegMC(), cascMC.pzNegMC()}));
        histos.fill(HIST("h2dGenOmegaMinusEtaBach"), RecoDecay::eta(std::array{cascMC.pxBachMC(), cascMC.pyBachMC(), cascMC.pzBachMC()}));
        histos.fill(HIST("h2dGenOmegaMinusVsMultMCVsCentrality"), mcCollision.multMCNParticlesEta05(), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaMinusVsMultMCVsIR"), mcCollision.multMCNParticlesEta05(), intRate, ptmc);
        histos.fill(HIST("h2dGenOmegaMinusVsCentOccupancy"), ptmc, centrality, occupancy);
        histos.fill(HIST("h2dGenOmegaMinusVsNchVsOccupancy"), ptmc, nChEta1, occupancy);
        histos.fill(HIST("h2dGenOmegaMinusVsCentIR"), ptmc, centrality, intRate);
      }
      if (cascMC.pdgCode() == PDG_t::kOmegaPlusBar && !isXi) {
        histos.fill(HIST("h2dGenOmegaPlus"), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaPlusVsNch"), nChEta1, ptmc);
        histos.fill(HIST("h2dGenOmegaPlusVsMultMCVsCentrality"), mcCollision.multMCNParticlesEta05(), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaPlusVsMultMCVsIR"), mcCollision.multMCNParticlesEta05(), intRate, ptmc);
        histos.fill(HIST("h2dGenOmegaPlusVsCentOccupancy"), ptmc, centrality, occupancy);
        histos.fill(HIST("h2dGenOmegaPlusVsNchVsOccupancy"), ptmc, nChEta1, occupancy);
        histos.fill(HIST("h2dGenOmegaPlusVsCentIR"), ptmc, centrality, intRate);
      }
    }
  }
  // ______________________________________________________
  // Simulated processing
  // Fill event information (for event loss estimation) and return the index to the recoed collision associated to a given MC collision.
  template <bool run3, typename TMCCollisions, typename TCollisions>
  std::vector<int> fillGenEventHist(TMCCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 0.5 /* all gen. events*/);

      if (eventSelectionCommonFlags.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelectionCommonFlags.zVertexCut) {
        continue;
      }
      histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 1.5 /* gen. events with vertex cut*/);

      if (eventSelectionCommonFlags.doInel0MCGen && mcCollision.multMCNParticlesEta10() < 1) {
        continue;
      }
      histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 2.5 /* gen. events with INEL>0t*/);

      auto groupedCollisions = getGroupedCollisions<run3>(collisions, mcCollision.globalIndex());
      // Check if there is at least one of the reconstructed collisions associated to this MC collision
      // If so, we consider it
      bool atLeastOne = false;
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      float centrality = 100.5f;
      int nCollisions = 0;
      float nChEta05 = -1;
      double intRate = -1;
      for (auto const& collision : groupedCollisions) {
        if (!isEventAccepted(collision, false)) {
          continue;
        }

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          bestCollisionIndex = collision.globalIndex();
          if constexpr (requires { collision.centFT0C(); }) {
            centrality = collision.centFT0C();
            intRate = rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource) * 1.e-3;
            if (useCentralityFT0M)
              centrality = collision.centFT0M();
            if (useCentralityFT0A)
              centrality = collision.centFV0A();
            if (useCentralityFT0Cvar1)
              centrality = collision.centFT0CVariant1();
          } else {
            centrality = eventSelectionRun2Flags.useSPDTrackletsCent ? collision.centRun2SPDTracklets() : collision.centRun2V0M();
          }
        }
        nCollisions++;

        atLeastOne = true;
      }
      listBestCollisionIdx[mcCollision.globalIndex()] = bestCollisionIndex;

      histos.fill(HIST("hCentralityVsNcoll_beforeEvSel"), centrality, groupedCollisions.size() + 0.5);
      histos.fill(HIST("hCentralityVsNcoll_afterEvSel"), centrality, nCollisions + 0.5);

      histos.fill(HIST("hCentralityVsIRGen"), centrality, intRate);
      histos.fill(HIST("hMultMCVsCentralityVsIRGen"), mcCollision.multMCNParticlesEta05(), centrality, intRate);

      histos.fill(HIST("hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta05());
      histos.fill(HIST("hRecMultVsMultMC"), nChEta05, mcCollision.multMCNParticlesEta05());

      histos.fill(HIST("hGenMultMCFT0C"), mcCollision.multMCFT0C());
      histos.fill(HIST("hGenMCNParticlesEta10"), mcCollision.multMCNParticlesEta10());

      if (atLeastOne) {
        histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 3.5 /* at least 1 rec. event*/);
      }
    }
    return listBestCollisionIdx;
  }

  void processCascades(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFNSigmas> const& Cascades, DauTracks const&)
  {
    analyseCascades(coll, Cascades);
  }
  void processCascadesRun2(soa::Join<aod::StraCollisions, aod::StraCentsRun2, aod::StraEvSelsRun2, aod::StraStamps>::iterator const& coll, soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFNSigmas> const& Cascades, DauTracks const&)
  {
    analyseCascades(coll, Cascades);
  }

  void processCascadesMCrec(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels, aod::StraStamps>::iterator const& coll, CascMCCandidates const& Cascades, DauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {
    analyseCascades(coll, Cascades);
  }
  void processCascadesMCrecRun2(soa::Join<aod::StraCollisions, aod::StraCentsRun2, aod::StraEvSelsRun2, aod::StraCollLabels, aod::StraStamps>::iterator const& coll, CascMCCandidates const& Cascades, DauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {
    analyseCascades(coll, Cascades);
  }

  void processCascadesMCforEff(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const& Cascades, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels, aod::StraStamps> const& collisions)
  {
    analyseCascadesMCforEff<true>(mcCollisions, Cascades, collisions);
  }
  void processCascadesMCforEffRun2(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const& Cascades, soa::Join<aod::StraCollisions, aod::StraCentsRun2, aod::StraEvSelsRun2, aod::StraCollLabels, aod::StraStamps> const& collisions)
  {
    analyseCascadesMCforEff<false>(mcCollisions, Cascades, collisions);
  }
  PROCESS_SWITCH(Derivedcascadeanalysis, processCascades, "cascade analysis, run3 data ", true);
  PROCESS_SWITCH(Derivedcascadeanalysis, processCascadesRun2, "cascade analysis, run2 data ", false);
  PROCESS_SWITCH(Derivedcascadeanalysis, processCascadesMCrec, "cascade analysis, run3 rec MC", false);
  PROCESS_SWITCH(Derivedcascadeanalysis, processCascadesMCrecRun2, "cascade analysis, run2 rec MC", false);
  PROCESS_SWITCH(Derivedcascadeanalysis, processCascadesMCforEff, "cascade analysis, run3 rec MC", false);
  PROCESS_SWITCH(Derivedcascadeanalysis, processCascadesMCforEffRun2, "cascade analysis, run2 rec MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Derivedcascadeanalysis>(cfgc)};
}
