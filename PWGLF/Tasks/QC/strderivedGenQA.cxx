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
// This code does basic QA of strangeness derived data
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  Strange Derived Generation QA
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <vector>
#include <cstdlib>

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
using std::array;
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using StrCollisionsDatas = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>;
using V0DerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using CascDerivedMCDatas = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascBBs, aod::CascCoreMCLabels>;
using CascDerivedDatas = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascBBs>;

struct strderivedGenQA {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // pack track quality but separte also afterburner
  // dynamic range: 0-31
  enum selection : int { hasTPC = 0,
                         hasITSTracker,
                         hasITSAfterburner,
                         hasTRD,
                         hasTOF };

  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};

  struct : ConfigurableGroup {
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF"};
    Configurable<bool> requireINEL0{"requireINEL0", true, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};

    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};

    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
  } eventSelections;

  struct : ConfigurableGroup {
    Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};

    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

    // Standard 5 topological criteria
    Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // Track quality
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
    Configurable<bool> skipTPConly{"skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
    Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};
    Configurable<bool> rejectPosITSafterburner{"rejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
    Configurable<bool> rejectNegITSafterburner{"rejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};

    // PID (TPC/TOF)
    Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
    Configurable<float> TofPidNsigmaCutLaPr{"TofPidNsigmaCutLaPr", 1e+6, "TofPidNsigmaCutLaPr"};
    Configurable<float> TofPidNsigmaCutLaPi{"TofPidNsigmaCutLaPi", 1e+6, "TofPidNsigmaCutLaPi"};
    Configurable<float> TofPidNsigmaCutK0Pi{"TofPidNsigmaCutK0Pi", 1e+6, "TofPidNsigmaCutK0Pi"};

    // PID (TOF)
    Configurable<float> maxDeltaTimeProton{"maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimePion{"maxDeltaTimePion", 1e+9, "check maximum allowed time"};
  } v0Selections;

  // Axis declarations
  ConfigurableAxis axisPosZ{"axisPosZ", {100, -50.0f, 50.0f}, "Z Position"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};

  // Boolean axes
  ConfigurableAxis axisBool{"axisBool", {2, 0.0f, 2.0f}, "axisBool"};
  ConfigurableAxis axisFt0cOccupancyInTimeRange{"axisFt0cOccupancyInTimeRange", {50, 0, 80000}, "FT0C occupancy"};
  ConfigurableAxis axisTrackOccupancyInTimeRange{"axisTrackOccupancyInTimeRange", {50, 0, 5000}, "Track occupancy"};
  ConfigurableAxis axisMultFT0C{"axisMultFT0C", {1000, 0, 100000}, "FT0C amplitude"};
  ConfigurableAxis axisMultNTracksPVeta1{"axisMultNTracksPVeta1", {200, 0, 6000}, "Mult NTracks PV eta 1"};
  ConfigurableAxis axisMultPVTotalContributors{"axisMultPVTotalContributors", {200, 0, 6000}, "Number of PV Contributors"};
  ConfigurableAxis axisMultAllTracksTPCOnly{"axisMultAllTracksTPCOnly", {200, 0, 6000}, "Mult All Tracks TPC Only"};
  ConfigurableAxis axisMultAllTracksITSTPC{"axisMultAllTracksITSTPC", {200, 0, 6000}, "Mult All Tracks ITS TPC"};
  ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles"};
  ConfigurableAxis axisNumV0sPerColl{"axisNumV0sPerColl", {50000, -0.5f, 49999.5f}, "Num V0s Per Coll"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "Pt Axis"};
  ConfigurableAxis axisPt2{"axisPt2", {VARIABLE_WIDTH, -50.0f, -40.0f, -35.0f, -30.0f, -25.0f, -23.0f, -21.0f, -19.0f, -17.0f, -15.0f, -14.0f, -13.0f, -12.0f, -11.0f, -10.0f, -9.0f, -8.0f, -7.5f, -7.0f, -6.5f, -6.0f, -5.6f, -5.2f, -4.8f, -4.4f, -4.0f, -3.8f, -3.6f, -3.4f, -3.2f, -3.0f, -2.8f, -2.6f, -2.4f, -2.2f, -2.0f, -1.9f, -1.8f, -1.7f, -1.6f, -1.5f, -1.4f, -1.3f, -1.2f, -1.1f, -1.0f, -0.9f, -0.8f, -0.7f, -0.6f, -0.5f, -0.4f, -0.3f, -0.2f, -0.1f, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "Pt Axis"};
  ConfigurableAxis axisAlpha{"axisAlpha", {220, -1.1f, 1.1f}, "Alpha"};
  ConfigurableAxis axisQtarm{"axisQtarm", {220, 0.0f, 0.5f}, "Qtarm"};
  ConfigurableAxis axisCosPA{"axisCosPA", {240, 0.0f, 1.2f}, "CosPA"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA V0 Daughters"};
  ConfigurableAxis axisEta{"axisEta", {100, -3.0f, 3.0f}, "Eta"};
  ConfigurableAxis axisPhi{"axisPhi", {100, 0.0f, TMath::TwoPi()}, "Phi"};
  ConfigurableAxis axisMassGamma{"axisMassGamma", {400, 0.0f, 0.5f}, "Mass Gamma"};
  ConfigurableAxis axisMassLambda{"axisMassLambda", {400, 1.0f, 1.2f}, "Mass Lambda"};
  ConfigurableAxis axisMassK0Short{"axisMassK0Short", {400, 0.4f, 0.6f}, "Mass K0Short"};
  ConfigurableAxis axisV0Type{"axisV0Type", {5, 0.0f, 5.0f}, "V0 Type"};
  ConfigurableAxis axisStraCollisionId{"axisStraCollisionId", {4000, 0.0f, 40000.0f}, "Stra Collision Id"};
  ConfigurableAxis axisGlobalIndex{"axisGlobalIndex", {4000, 0.0f, 40000.0f}, "Global Index"};
  ConfigurableAxis axisNCls{"axisNCls", {8, -0.5, 7.5}, "NCls"};
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisChi2PerNcl{"axisChi2PerNcl", {100, 0, 40}, "Chi2 Per Ncl"};
  ConfigurableAxis axisTPCNSigma{"axisTPCNSigma", {100, -50.0f, 50.0f}, "TPC N Sigma"};
  ConfigurableAxis axisTPCSignal{"axisTPCSignal", {400, -100.0, 300.0}, "TPC Signal"};
  ConfigurableAxis axisTOFNSigma{"axisTOFNSigma", {100, -50.0f, 50.0f}, "TOF N Sigma"};
  ConfigurableAxis axisTOFDeltaT{"axisTOFDeltaT", {200, -1000.0f, +1000.0f}, "TOF Delta T"};
  ConfigurableAxis axisPtResolution{"axisPtResolution", {100, -1.0f, 1.0f}, "Pt Resolution"};
  ConfigurableAxis axisPDGCode{"axisPDGCode", {10001, -5000.5f, +5000.5f}, "PDG Code"};
  ConfigurableAxis axisV0Radius{"axisV0Radius", {400, 0.0f, 200.0f}, "V0 Radius"};
  ConfigurableAxis axisCascRadius{"axisCascRadius", {500, 0.0f, 50.0f}, "Casc Radius"};
  ConfigurableAxis axisDCAToPV{"axisDCAToPV", {500, -50.0f, 50.0f}, "DCA Dau to PV"};
  ConfigurableAxis axisDCAXYCascToPV{"axisDCAXYCascToPV", {1000, 0.0f, 10.0f}, "DCA XY Casc to PV"};
  ConfigurableAxis axisDCAZCascToPV{"axisDCAZCascToPV", {500, -10.0f, 10.0f}, "DCA Z Casc to PV"};
  ConfigurableAxis axisDCAV0ToPV{"axisDCAV0ToPV", {1000, -10.0f, 10.0f}, "DCA V0 to PV"};
  ConfigurableAxis axisDCAV0Dau{"axisDCAV0Dau", {1000, 0.0f, 10.0f}, "DCA V0 Daughters"};
  ConfigurableAxis axisDCACascDau{"axisDCACascDau", {1000, 0.0f, 10.0f}, "DCA Casc Daughters"};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {400, 1.6f, 1.8f}, "Omega Mass"};
  ConfigurableAxis axisXiMass{"axisXiMass", {400, 1.2f, 1.4f}, "Xi Mass"};
  ConfigurableAxis axisTrackProperties{"axisTrackProperties", {32, -0.5, 31.5f}, "Track Properties"};

  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

  void init(InitContext const&)
  {
    // Histogram declarations (can be improved!)
    histos.add("Event/hPosZ", "hPosZ", kTH1F, {axisPosZ});

    // Event Counters
    histos.add("Event/hEventProperties", "hEventProperties", kTH1F, {{20, -0.5f, +18.5f}});
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(11, "kNoCollInTimeRangeStd");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeStrict");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(14, "kNoCollInRofStd");
    histos.get<TH1>(HIST("Event/hEventProperties"))->GetXaxis()->SetBinLabel(15, "kNoCollInRofStrict");

    histos.add("Event/hft0cOccupancyInTimeRange", "hft0cOccupancyInTimeRange", kTH1F, {axisFt0cOccupancyInTimeRange});
    histos.add("Event/htrackOccupancyInTimeRange", "htrackOccupancyInTimeRange", kTH1F, {axisTrackOccupancyInTimeRange});
    histos.add("Event/h2dMultFT0C", "h2dMultFT0C", kTH2F, {axisCentrality, axisMultFT0C});
    histos.add("Event/h2dMultNTracksPVeta1", "h2dMultNTracksPVeta1", kTH2F, {axisCentrality, axisMultNTracksPVeta1});
    histos.add("Event/h2dMultPVTotalContributors", "h2dMultPVTotalContributors", kTH2F, {axisCentrality, axisMultPVTotalContributors});
    histos.add("Event/h2dMultAllTracksTPCOnly", "h2dMultAllTracksTPCOnly", kTH2F, {axisCentrality, axisMultAllTracksTPCOnly});
    histos.add("Event/h2dMultAllTracksITSTPC", "h2dMultAllTracksITSTPC", kTH2F, {axisCentrality, axisMultAllTracksITSTPC});
    histos.add("Event/h2dNumV0sPerColl", "h2dNumV0sPerColl", kTH2F, {axisCentrality, axisNumV0sPerColl});

    histos.add("V0/hpT", "hpT", kTH1F, {axisPt});
    histos.add("V0/h2dArmenterosP", "h2dArmenterosP", kTH2F, {axisAlpha, axisQtarm});
    histos.add("V0/hRadius", "hRadius", kTH1F, {axisV0Radius});
    histos.add("V0/hZ", "hZ", kTH1F, {axisPosZ});
    histos.add("V0/hCosPA", "hCosPA", kTH1F, {axisCosPA});
    histos.add("V0/hdcaDau", "hdcaDau", kTH1F, {axisDCAdau});
    histos.add("V0/hdcaNegtopv", "hdcaNegtopv", kTH1F, {axisDCAToPV});
    histos.add("V0/hdcaPostopv", "hdcaPostopv", kTH1F, {axisDCAToPV});
    histos.add("V0/h2dEtaPhi", "h2dEtaPhi", kTH2F, {axisEta, axisPhi});
    histos.add("V0/hYGamma", "hYGamma", kTH1F, {axisEta});
    histos.add("V0/hYLambda", "hYLambda", kTH1F, {axisEta});
    histos.add("V0/hYK0Short", "hYK0Short", kTH1F, {axisEta});
    histos.add("V0/hMassGamma", "hMassGamma", kTH1F, {axisMassGamma});
    histos.add("V0/hMassLambda", "hMassLambda", kTH1F, {axisMassLambda});
    histos.add("V0/hMassK0Short", "hMassK0Short", kTH1F, {axisMassK0Short});
    histos.add("V0/hV0Type", "hV0Type", kTH1F, {axisV0Type});
    histos.add("V0/h2dV0Indices", "h2dV0Indices", kTH2F, {axisStraCollisionId, axisGlobalIndex});

    histos.add("V0/Track/h2dITSNCls", "h2dITSNCls", kTH2F, {axisPt2, axisNCls});
    histos.add("V0/Track/h2dITSChi2PerNcl", "h2dITSChi2PerNcl", kTH2F, {axisPt2, axisChi2PerNcl});
    histos.add("V0/Track/h2dTPCCrossedRows", "h2dTPCCrossedRows", kTH2F, {axisPt2, axisTPCrows});

    histos.add("V0/Track/h2dPosTrackProperties", "h2dPosTrackProperties", kTH2F, {axisTrackProperties, axisPt});
    histos.add("V0/Track/h3dTrackPropertiesVspT", "h3dTrackPropertiesVspT", kTH3F, {axisTrackProperties, axisTrackProperties, axisPt});
    histos.add("V0/Track/h2dNegTrackProperties", "h2dNegTrackProperties", kTH2F, {axisTrackProperties, axisPt});

    // Add histogram to the list
    histos.add("V0/Track/hTrackCode", "hTrackCode", kTH1F, {axisTrackProperties});

    // Set bin labels for all combinations
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(1, "None");                                     // Code 0
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(2, "TPC");                                      // Code 1
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(3, "ITSTracker");                               // Code 2
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(4, "ITSTracker + TPC");                         // Code 3
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(5, "ITSAfterburner");                           // Code 4
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(6, "ITSAfterburner + TPC");                     // Code 5
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(7, "ITSAfterburner + ITSTracker");              // Code 6
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(8, "ITSAfterburner + ITSTracker + TPC");        // Code 7
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(9, "TRD");                                      // Code 8
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(10, "TRD + TPC");                               // Code 9
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(11, "TRD + ITSTracker");                        // Code 10
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(12, "TRD + ITSTracker + TPC");                  // Code 11
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(13, "TRD + ITSAfterburner");                    // Code 12
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(14, "TRD + ITSAfterburner + TPC");              // Code 13
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(15, "TRD + ITSAfterburner + ITSTracker");       // Code 14
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(16, "TRD + ITSAfterburner + ITSTracker + TPC"); // Code 15
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(17, "TOF");                                     // Code 16
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(18, "TOF + TPC");                               // Code 17
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(19, "TOF + ITSTracker");                        // Code 18
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(20, "TOF + ITSTracker + TPC");                  // Code 19
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(21, "TOF + ITSAfterburner");                    // Code 20
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(22, "TOF + ITSAfterburner + TPC");              // Code 21
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(23, "TOF + ITSAfterburner + ITSTracker");       // Code 22
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(24, "TOF + ITSAfterburner + ITSTracker + TPC"); // Code 23
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(25, "TOF + TRD");                               // Code 24
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(26, "TOF + TRD + TPC");                         // Code 25
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(27, "TOF + TRD + ITSTracker");                  // Code 26
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(28, "TOF + TRD + ITSTracker + TPC");            // Code 27
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(29, "TOF + TRD + ITSAfterburner");              // Code 28
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(30, "TOF + TRD + ITSAfterburner + TPC");        // Code 29
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(31, "TOF + TRD + ITSAfterburner + ITSTracker"); // Code 30
    histos.get<TH1>(HIST("V0/Track/hTrackCode"))->GetXaxis()->SetBinLabel(32, "All");                                     // Code 31

    histos.add("V0/PID/h2dTPCNSigmaEl", "h2dTPCNSigmaEl", kTH2F, {axisPt2, axisTPCNSigma});
    histos.add("V0/PID/h2dTPCNSigmaPr", "h2dTPCNSigmaPr", kTH2F, {axisPt2, axisTPCNSigma});
    histos.add("V0/PID/h2dTPCNSigmaPi", "h2dTPCNSigmaPi", kTH2F, {axisPt2, axisTPCNSigma});
    histos.add("V0/PID/h2dTPCSignal", "h2dTPCSignal", kTH2F, {axisPt2, axisTPCSignal});

    histos.add("V0/PID/h2dTOFNSigmaLaPr", "h2dTOFNSigmaLaPr", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("V0/PID/h2dTOFNSigmaLaPi", "h2dTOFNSigmaLaPi", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("V0/PID/h2dposTOFDeltaTLaPr", "h2dposTOFDeltaTLaPr", kTH2F, {axisPt, axisTOFDeltaT});
    histos.add("V0/PID/h2dnegTOFDeltaTLaPi", "h2dnegTOFDeltaTLaPi", kTH2F, {axisPt, axisTOFDeltaT});
    histos.add("V0/PID/h2dnegTOFDeltaTLaPr", "h2dnegTOFDeltaTLaPr", kTH2F, {axisPt, axisTOFDeltaT});
    histos.add("V0/PID/h2dposTOFDeltaTLaPi", "h2dposTOFDeltaTLaPi", kTH2F, {axisPt, axisTOFDeltaT});
    histos.add("V0/PID/h2dTOFNSigmaALaPr", "h2dTOFNSigmaALaPr", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("V0/PID/h2dTOFNSigmaALaPi", "h2dTOFNSigmaALaPi", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("V0/PID/h2dTOFNSigmaK0PiPlus", "h2dTOFNSigmaK0PiPlus", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("V0/PID/h2dTOFNSigmaK0PiMinus", "h2dTOFNSigmaK0PiMinus", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("V0/PID/h3dTPCVsTOFNSigmaLaPr", "h3dTPCVsTOFNSigmaLaPr", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});
    histos.add("V0/PID/h3dTPCVsTOFNSigmaLaPi", "h3dTPCVsTOFNSigmaLaPi", kTH3F, {axisTPCNSigma, axisTOFNSigma, axisPt});

    histos.add("MCV0/hv0MCCore", "hv0MCCore", kTH1F, {axisBool});
    histos.add("MCV0/h2dPDGV0VsMother", "h2dPDGV0VsMother", kTHnSparseD, {axisPDGCode, axisPDGCode});
    histos.add("MCV0/h2dPDGV0VsPositive", "h2dPDGV0VsPositive", kTHnSparseD, {axisPDGCode, axisPDGCode});
    histos.add("MCV0/h2dPDGV0VsNegative", "h2dPDGV0VsNegative", kTHnSparseD, {axisPDGCode, axisPDGCode});
    histos.add("MCV0/h2dPDGV0VsIsPhysicalPrimary", "h2dPDGV0VsIsPhysicalPrimary", kTH2F, {axisPDGCode, axisBool});
    histos.add("MCV0/h2dArmenterosP", "h2dArmenterosP", kTH2F, {axisAlpha, axisQtarm});
    histos.add("MCV0/Gamma/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCV0/Gamma/h2dMass", "h2dMass", kTH2F, {axisPt, axisMassGamma});
    histos.add("MCV0/Gamma/h2dTPCNSigmaEl", "h2dTPCNSigmaEl", kTH2F, {axisPt2, axisTPCNSigma});
    histos.add("MCV0/Gamma/h2dTPCSignal", "h2dTPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCV0/Gamma/hRadius", "hRadius", kTH1F, {axisV0Radius});
    histos.add("MCV0/Gamma/hCosPA", "hCosPA", kTH1F, {axisCosPA});
    histos.add("MCV0/Gamma/hdcaDau", "hdcaDau", kTH1F, {axisDCAdau});
    histos.add("MCV0/Gamma/hdcaNegtopv", "hdcaNegtopv", kTH1F, {axisDCAToPV});
    histos.add("MCV0/Gamma/hdcaPostopv", "hdcaPostopv", kTH1F, {axisDCAToPV});

    histos.add("MCV0/Lambda/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCV0/Lambda/h2dMass", "h2dMass", kTH2F, {axisPt, axisMassLambda});
    histos.add("MCV0/Lambda/h2dTPCNSigmaPr", "h2dTPCNSigmaPr", kTH2F, {axisPt, axisTPCNSigma});
    histos.add("MCV0/Lambda/h2dTPCNSigmaPi", "h2dTPCNSigmaPi", kTH2F, {axisPt, axisTPCNSigma});
    histos.add("MCV0/Lambda/h2dTPCSignal", "h2dTPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCV0/Lambda/hRadius", "hRadius", kTH1F, {axisV0Radius});
    histos.add("MCV0/Lambda/hCosPA", "hCosPA", kTH1F, {axisCosPA});
    histos.add("MCV0/Lambda/hdcaDau", "hdcaDau", kTH1F, {axisDCAdau});
    histos.add("MCV0/Lambda/hdcaNegtopv", "hdcaNegtopv", kTH1F, {axisDCAToPV});
    histos.add("MCV0/Lambda/hdcaPostopv", "hdcaPostopv", kTH1F, {axisDCAToPV});

    histos.add("MCV0/AntiLambda/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCV0/AntiLambda/h2dMass", "h2dMass", kTH2F, {axisPt, axisMassLambda});
    histos.add("MCV0/AntiLambda/h2dTPCNSigmaPr", "h2dTPCNSigmaPr", kTH2F, {axisPt, axisTPCNSigma});
    histos.add("MCV0/AntiLambda/h2dTPCNSigmaPi", "h2dTPCNSigmaPi", kTH2F, {axisPt, axisTPCNSigma});
    histos.add("MCV0/AntiLambda/h2dTPCSignal", "h2dTPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCV0/AntiLambda/hRadius", "hRadius", kTH1F, {axisV0Radius});
    histos.add("MCV0/AntiLambda/hCosPA", "hCosPA", kTH1F, {axisCosPA});
    histos.add("MCV0/AntiLambda/hdcaDau", "hdcaDau", kTH1F, {axisDCAdau});
    histos.add("MCV0/AntiLambda/hdcaNegtopv", "hdcaNegtopv", kTH1F, {axisDCAToPV});
    histos.add("MCV0/AntiLambda/hdcaPostopv", "hdcaPostopv", kTH1F, {axisDCAToPV});

    histos.add("MCV0/K0Short/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCV0/K0Short/h2dMass", "h2dMass", kTH2F, {axisPt, axisMassK0Short});
    histos.add("MCV0/K0Short/h2dTPCNSigmaPi", "h2dTPCNSigmaPi", kTH2F, {axisPt2, axisTPCNSigma});
    histos.add("MCV0/K0Short/h2dTPCSignal", "h2dTPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCV0/K0Short/hRadius", "hRadius", kTH1F, {axisV0Radius});
    histos.add("MCV0/K0Short/hCosPA", "hCosPA", kTH1F, {axisCosPA});
    histos.add("MCV0/K0Short/hdcaDau", "hdcaDau", kTH1F, {axisDCAdau});
    histos.add("MCV0/K0Short/hdcaNegtopv", "hdcaNegtopv", kTH1F, {axisDCAToPV});
    histos.add("MCV0/K0Short/hdcaPostopv", "hdcaPostopv", kTH1F, {axisDCAToPV});

    histos.add("Casc/Sign", "Sign", kTH1F, {{3, -1.5f, 1.5f}});
    histos.add("Casc/hpT", "hpT", kTH1F, {axisPt});
    histos.add("Casc/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("Casc/hCascRadius", "hCascRadius", kTH1F, {axisCascRadius});
    histos.add("Casc/hV0CosPA", "hV0CosPA", kTH1F, {axisCosPA});
    histos.add("Casc/hCascCosPA", "hCascCosPA", kTH1F, {axisCosPA});
    histos.add("Casc/hDCAPosToPV", "hDCAPosToPV", kTH1F, {axisDCAToPV});
    histos.add("Casc/hDCANegToPV", "hDCANegToPV", kTH1F, {axisDCAToPV});
    histos.add("Casc/hDCABachToPV", "hDCABachToPV", kTH1F, {axisDCAToPV});
    histos.add("Casc/hDCAXYCascToPV", "hDCAXYCascToPV", kTH1F, {axisDCAXYCascToPV});
    histos.add("Casc/hDCAZCascToPV", "hDCAZCascToPV", kTH1F, {axisDCAZCascToPV});
    histos.add("Casc/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
    histos.add("Casc/hDCAV0Dau", "hDCAV0Dau", kTH1F, {axisDCAV0Dau});
    histos.add("Casc/hDCACascDau", "hDCACascDau", kTH1F, {axisDCACascDau});
    histos.add("Casc/hLambdaMass", "hLambdaMass", kTH1F, {axisMassLambda});

    histos.add("Casc/Track/h3dTrackProperties", "h3dTrackProperties", kTH3F, {axisTrackProperties, axisTrackProperties, axisTrackProperties});
    histos.add("Casc/Track/h2dPosTrackProperties", "h2dPosTrackProperties", kTH2F, {axisTrackProperties, axisPt});
    histos.add("Casc/Track/h2dNegTrackProperties", "h2dNegTrackProperties", kTH2F, {axisTrackProperties, axisPt});
    histos.add("Casc/Track/h2dBachTrackProperties", "h2dBachTrackProperties", kTH2F, {axisTrackProperties, axisPt});
    histos.add("Casc/Track/h2dV0ITSChi2PerNcl", "h2dV0ITSChi2PerNcl", kTH2F, {axisPt2, axisChi2PerNcl});
    histos.add("Casc/Track/h2dV0TPCCrossedRows", "h2dV0TPCCrossedRows", kTH2F, {axisPt2, axisTPCrows});
    histos.add("Casc/Track/h2dV0ITSNCls", "h2dV0ITSNCls", kTH2F, {axisPt2, axisNCls});

    histos.add("Casc/PID/h2dV0TPCNSigmaPr", "h2dV0TPCNSigmaPr", kTH2F, {axisPt2, axisTPCNSigma});
    histos.add("Casc/PID/h2dV0TPCNSigmaPi", "h2dV0TPCNSigmaPi", kTH2F, {axisPt2, axisTPCNSigma});
    histos.add("Casc/PID/h2dV0TPCSignal", "h2dV0TPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("Casc/PID/h2dTOFNSigmaXiLaPi", "h2dTOFNSigmaXiLaPi", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("Casc/PID/h2dTOFNSigmaXiLaPr", "h2dTOFNSigmaXiLaPr", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("Casc/PID/h2dTOFNSigmaXiPi", "h2dTOFNSigmaXiPi", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("Casc/PID/h2dTOFNSigmaOmLaPi", "h2dTOFNSigmaOmLaPi", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("Casc/PID/h2dTOFNSigmaOmLaPr", "h2dTOFNSigmaOmLaPr", kTH2F, {axisPt, axisTOFNSigma});
    histos.add("Casc/PID/h2dTOFNSigmaOmKa", "h2dTOFNSigmaOmKa", kTH2F, {axisPt, axisTOFNSigma});

    histos.add("Casc/hMassXiMinus", "hMassXiMinus", kTH1F, {axisXiMass});
    histos.add("Casc/hMassOmegaMinus", "hMassOmegaMinus", kTH1F, {axisOmegaMass});
    histos.add("Casc/hMassXiPlus", "hMassXiPlus", kTH1F, {axisXiMass});
    histos.add("Casc/hMassOmegaPlus", "hMassOmegaPlus", kTH1F, {axisOmegaMass});
    histos.add("Casc/Track/h2dBachITSNCls", "h2dBachITSNCls", kTH2F, {axisPt2, axisNCls});
    histos.add("Casc/Track/h2dBachITSChi2PerNcl", "h2dBachITSChi2PerNcl", kTH2F, {axisPt2, axisChi2PerNcl});
    histos.add("Casc/Track/h2dBachTPCCrossedRows", "h2dBachTPCCrossedRows", kTH2F, {axisPt2, axisTPCrows});
    histos.add("Casc/PID/h2dBachTPCSignal", "h2dBachTPCSignal", kTH2F, {axisPt2, axisTPCSignal});

    histos.add("MCCasc/hcascMCCore", "hcascMCCore", kTH1F, {axisBool});
    histos.add("MCCasc/h2dPDGV0VsMother", "h2dPDGV0VsMother", kTHnSparseD, {axisPDGCode, axisPDGCode});
    histos.add("MCCasc/h2dPDGV0VsPositive", "h2dPDGV0VsPositive", kTHnSparseD, {axisPDGCode, axisPDGCode});
    histos.add("MCCasc/h2dPDGV0VsNegative", "h2dPDGV0VsNegative", kTHnSparseD, {axisPDGCode, axisPDGCode});
    histos.add("MCCasc/h2dPDGV0VsBach", "h2dPDGV0VsBach", kTHnSparseD, {axisPDGCode, axisPDGCode});
    histos.add("MCCasc/h2dPDGV0VsIsPhysicalPrimary", "h2dPDGV0VsIsPhysicalPrimary", kTH2F, {axisPDGCode, axisBool});

    histos.add("MCCasc/XiMinus/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCCasc/XiMinus/h2dMass", "h2dMass", kTH2F, {axisPt, axisXiMass});
    histos.add("MCCasc/XiMinus/h2dV0TPCSignal", "h2dV0TPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCCasc/XiMinus/h2dBachTPCSignal", "h2dBachTPCSignal", kTH2F, {axisPt, axisTPCSignal});
    histos.add("MCCasc/XiMinus/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("MCCasc/XiMinus/hCascRadius", "hCascRadius", kTH1F, {axisCascRadius});
    histos.add("MCCasc/XiMinus/hV0CosPA", "hV0CosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/XiMinus/hCascCosPA", "hCascCosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/XiMinus/hDCAPosToPV", "hDCAPosToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/XiMinus/hDCANegToPV", "hDCANegToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/XiMinus/hDCABachToPV", "hDCABachToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/XiMinus/hDCAXYCascToPV", "hDCAXYCascToPV", kTH1F, {axisDCAXYCascToPV});
    histos.add("MCCasc/XiMinus/hDCAZCascToPV", "hDCAZCascToPV", kTH1F, {axisDCAZCascToPV});
    histos.add("MCCasc/XiMinus/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
    histos.add("MCCasc/XiMinus/hDCAV0Dau", "hDCAV0Dau", kTH1F, {axisDCAV0Dau});
    histos.add("MCCasc/XiMinus/hDCACascDau", "hDCACascDau", kTH1F, {axisDCACascDau});
    histos.add("MCCasc/XiMinus/hLambdaMass", "hLambdaMass", kTH1F, {axisMassLambda});

    histos.add("MCCasc/XiPlus/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCCasc/XiPlus/h2dMass", "h2dMass", kTH2F, {axisPt, axisXiMass});
    histos.add("MCCasc/XiPlus/h2dV0TPCSignal", "h2dV0TPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCCasc/XiPlus/h2dBachTPCSignal", "h2dBachTPCSignal", kTH2F, {axisPt, axisTPCSignal});
    histos.add("MCCasc/XiPlus/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("MCCasc/XiPlus/hCascRadius", "hCascRadius", kTH1F, {axisCascRadius});
    histos.add("MCCasc/XiPlus/hV0CosPA", "hV0CosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/XiPlus/hCascCosPA", "hCascCosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/XiPlus/hDCAPosToPV", "hDCAPosToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/XiPlus/hDCANegToPV", "hDCANegToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/XiPlus/hDCABachToPV", "hDCABachToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/XiPlus/hDCAXYCascToPV", "hDCAXYCascToPV", kTH1F, {axisDCAXYCascToPV});
    histos.add("MCCasc/XiPlus/hDCAZCascToPV", "hDCAZCascToPV", kTH1F, {axisDCAZCascToPV});
    histos.add("MCCasc/XiPlus/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
    histos.add("MCCasc/XiPlus/hDCAV0Dau", "hDCAV0Dau", kTH1F, {axisDCAV0Dau});
    histos.add("MCCasc/XiPlus/hDCACascDau", "hDCACascDau", kTH1F, {axisDCACascDau});
    histos.add("MCCasc/XiPlus/hLambdaMass", "hLambdaMass", kTH1F, {axisMassLambda});

    histos.add("MCCasc/OmegaMinus/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCCasc/OmegaMinus/h2dMass", "h2dMass", kTH2F, {axisPt, axisOmegaMass});
    histos.add("MCCasc/OmegaMinus/h2dV0TPCSignal", "h2dV0TPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCCasc/OmegaMinus/h2dBachTPCSignal", "h2dBachTPCSignal", kTH2F, {axisPt, axisTPCSignal});
    histos.add("MCCasc/OmegaMinus/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("MCCasc/OmegaMinus/hCascRadius", "hCascRadius", kTH1F, {axisCascRadius});
    histos.add("MCCasc/OmegaMinus/hV0CosPA", "hV0CosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/OmegaMinus/hCascCosPA", "hCascCosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/OmegaMinus/hDCAPosToPV", "hDCAPosToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/OmegaMinus/hDCANegToPV", "hDCANegToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/OmegaMinus/hDCABachToPV", "hDCABachToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/OmegaMinus/hDCAXYCascToPV", "hDCAXYCascToPV", kTH1F, {axisDCAXYCascToPV});
    histos.add("MCCasc/OmegaMinus/hDCAZCascToPV", "hDCAZCascToPV", kTH1F, {axisDCAZCascToPV});
    histos.add("MCCasc/OmegaMinus/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
    histos.add("MCCasc/OmegaMinus/hDCAV0Dau", "hDCAV0Dau", kTH1F, {axisDCAV0Dau});
    histos.add("MCCasc/OmegaMinus/hDCACascDau", "hDCACascDau", kTH1F, {axisDCACascDau});
    histos.add("MCCasc/OmegaMinus/hLambdaMass", "hLambdaMass", kTH1F, {axisMassLambda});

    histos.add("MCCasc/OmegaPlus/h2dpTResolution", "h2dpTResolution", kTH2F, {axisPt, axisPtResolution});
    histos.add("MCCasc/OmegaPlus/h2dMass", "h2dMass", kTH2F, {axisPt, axisOmegaMass});
    histos.add("MCCasc/OmegaPlus/h2dV0TPCSignal", "h2dV0TPCSignal", kTH2F, {axisPt2, axisTPCSignal});
    histos.add("MCCasc/OmegaPlus/h2dBachTPCSignal", "h2dBachTPCSignal", kTH2F, {axisPt, axisTPCSignal});
    histos.add("MCCasc/OmegaPlus/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
    histos.add("MCCasc/OmegaPlus/hCascRadius", "hCascRadius", kTH1F, {axisCascRadius});
    histos.add("MCCasc/OmegaPlus/hV0CosPA", "hV0CosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/OmegaPlus/hCascCosPA", "hCascCosPA", kTH1F, {axisCosPA});
    histos.add("MCCasc/OmegaPlus/hDCAPosToPV", "hDCAPosToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/OmegaPlus/hDCANegToPV", "hDCANegToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/OmegaPlus/hDCABachToPV", "hDCABachToPV", kTH1F, {axisDCAToPV});
    histos.add("MCCasc/OmegaPlus/hDCAXYCascToPV", "hDCAXYCascToPV", kTH1F, {axisDCAXYCascToPV});
    histos.add("MCCasc/OmegaPlus/hDCAZCascToPV", "hDCAZCascToPV", kTH1F, {axisDCAZCascToPV});
    histos.add("MCCasc/OmegaPlus/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
    histos.add("MCCasc/OmegaPlus/hDCAV0Dau", "hDCAV0Dau", kTH1F, {axisDCAV0Dau});
    histos.add("MCCasc/OmegaPlus/hDCACascDau", "hDCACascDau", kTH1F, {axisDCACascDau});
    histos.add("MCCasc/OmegaPlus/hLambdaMass", "hLambdaMass", kTH1F, {axisMassLambda});

    // MC Generated level
    histos.add("GenMC/hGenEvents", "hGenEvents", kTH2F, {{axisNch}, {2, -0.5f, +1.5f}});
    histos.get<TH2>(HIST("GenMC/hGenEvents"))->GetYaxis()->SetBinLabel(1, "All gen. events");
    histos.get<TH2>(HIST("GenMC/hGenEvents"))->GetYaxis()->SetBinLabel(2, "Gen. with at least 1 rec. events");
    histos.add("GenMC/hGenEventCentrality", "hGenEventCentrality", kTH1F, {{101, 0.0f, 101.0f}});
    histos.add("GenMC/hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", kTH2F, {axisCentrality, {50, -0.5f, 49.5f}});
    histos.add("GenMC/hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", kTH2F, {axisCentrality, {50, -0.5f, 49.5f}});
    histos.add("GenMC/hCentralityVsMultMC", "hCentralityVsMultMC", kTH2F, {{101, 0.0f, 101.0f}, axisNch});
    histos.add("GenMC/h2dGenGamma", "h2dGenGamma", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenK0Short", "h2dGenK0Short", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenLambda", "h2dGenLambda", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenAntiLambda", "h2dGenAntiLambda", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenXiMinus", "h2dGenXiMinus", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenXiPlus", "h2dGenXiPlus", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenOmegaMinus", "h2dGenOmegaMinus", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenOmegaPlus", "h2dGenOmegaPlus", kTH2D, {axisCentrality, axisPt});
    histos.add("GenMC/h2dGenK0ShortVsMultMC_RecoedEvt", "h2dGenK0ShortVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenLambdaVsMultMC_RecoedEvt", "h2dGenLambdaVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenAntiLambdaVsMultMC_RecoedEvt", "h2dGenAntiLambdaVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenXiMinusVsMultMC_RecoedEvt", "h2dGenXiMinusVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenXiPlusVsMultMC_RecoedEvt", "h2dGenXiPlusVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenOmegaMinusVsMultMC_RecoedEvt", "h2dGenOmegaMinusVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenOmegaPlusVsMultMC_RecoedEvt", "h2dGenOmegaPlusVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenGammaVsMultMC", "h2dGenGammaVsMultMC", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenK0ShortVsMultMC", "h2dGenK0ShortVsMultMC", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenLambdaVsMultMC", "h2dGenLambdaVsMultMC", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenAntiLambdaVsMultMC", "h2dGenAntiLambdaVsMultMC", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenXiMinusVsMultMC", "h2dGenXiMinusVsMultMC", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenXiPlusVsMultMC", "h2dGenXiPlusVsMultMC", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenOmegaMinusVsMultMC", "h2dGenOmegaMinusVsMultMC", kTH2D, {axisNch, axisPt});
    histos.add("GenMC/h2dGenOmegaPlusVsMultMC", "h2dGenOmegaPlusVsMultMC", kTH2D, {axisNch, axisPt});
  }

  template <typename TCollision>
  bool IsEventAccepted(TCollision collision)
  // check whether the collision passes our collision selections
  {
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (doPPAnalysis) { // we are in pp
      if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
        return false;
      }
      if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
        return false;
      }
    } else { // we are in Pb-Pb
      float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
        return false;
      }
      if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
        return false;
      }
    }

    return true;
  }

  // ______________________________________________________
  // Simulated processing
  // Return the list of indices to the recoed collision associated to a given MC collision.
  std::vector<int> getListOfRecoCollIndices(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels> const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      // Find the collision with the biggest nbr of PV contributors
      // Follows what was done here: https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/mcCollsExtra.cxx#L93
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      for (auto const& collision : groupedCollisions) {
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
  void fillGeneratedEventProperties(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels> const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("GenMC/hGenEvents"), mcCollision.multMCNParticlesEta05(), 0 /* all gen. events*/);

      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      // Check if there is at least one of the reconstructed collisions associated to this MC collision
      // If so, we consider it
      bool atLeastOne = false;
      int biggestNContribs = -1;
      float centrality = 100.5f;
      int nCollisions = 0;
      for (auto const& collision : groupedCollisions) {

        if (!IsEventAccepted(collision)) {
          continue;
        }

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
        }
        nCollisions++;

        atLeastOne = true;
      }

      histos.fill(HIST("GenMC/hCentralityVsNcoll_beforeEvSel"), centrality, groupedCollisions.size());
      histos.fill(HIST("GenMC/hCentralityVsNcoll_afterEvSel"), centrality, nCollisions);
      histos.fill(HIST("GenMC/hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta05());

      if (atLeastOne) {
        histos.fill(HIST("GenMC/hGenEvents"), mcCollision.multMCNParticlesEta05(), 1 /* at least 1 rec. event*/);
        histos.fill(HIST("GenMC/hGenEventCentrality"), centrality);
      }
    }
    return;
  }

  void processDerivedV0s(StrCollisionsDatas::iterator const& coll, V0DerivedDatas const& V0s, dauTracks const&)
  {
    // Event Level
    float centrality = coll.centFT0C();
    histos.fill(HIST("Event/hPosZ"), coll.posZ());
    histos.fill(HIST("Event/hEventProperties"), 0. /* all collisions */);

    if (coll.sel8())
      histos.fill(HIST("Event/hEventProperties"), 1.);
    if (coll.selection_bit(aod::evsel::kIsTriggerTVX))
      histos.fill(HIST("Event/hEventProperties"), 2.);
    if (coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      histos.fill(HIST("Event/hEventProperties"), 3.);
    if (coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      histos.fill(HIST("Event/hEventProperties"), 4.);
    if (coll.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      histos.fill(HIST("Event/hEventProperties"), 5.);
    if (coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      histos.fill(HIST("Event/hEventProperties"), 6.);
    if (coll.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      histos.fill(HIST("Event/hEventProperties"), 7.);
    if (coll.selection_bit(o2::aod::evsel::kIsVertexTRDmatched))
      histos.fill(HIST("Event/hEventProperties"), 8.);
    if (coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      histos.fill(HIST("Event/hEventProperties"), 9.);
    if (coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      histos.fill(HIST("Event/hEventProperties"), 10.);
    if (coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict))
      histos.fill(HIST("Event/hEventProperties"), 11.);
    if (coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow))
      histos.fill(HIST("Event/hEventProperties"), 12.);
    if (coll.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
      histos.fill(HIST("Event/hEventProperties"), 13.);
    if (coll.selection_bit(o2::aod::evsel::kNoCollInRofStrict))
      histos.fill(HIST("Event/hEventProperties"), 14.);

    histos.fill(HIST("Event/hft0cOccupancyInTimeRange"), coll.ft0cOccupancyInTimeRange());
    histos.fill(HIST("Event/htrackOccupancyInTimeRange"), coll.trackOccupancyInTimeRange());
    histos.fill(HIST("Event/h2dMultFT0C"), centrality, coll.multFT0C());
    histos.fill(HIST("Event/h2dMultNTracksPVeta1"), centrality, coll.multNTracksPVeta1());
    histos.fill(HIST("Event/h2dMultPVTotalContributors"), centrality, coll.multPVTotalContributors());
    histos.fill(HIST("Event/h2dMultAllTracksTPCOnly"), centrality, coll.multAllTracksTPCOnly());
    histos.fill(HIST("Event/h2dMultAllTracksITSTPC"), centrality, coll.multAllTracksITSTPC());
    histos.fill(HIST("Event/h2dNumV0sPerColl"), centrality, V0s.size());

    for (auto const& v0 : V0s) {

      // V0-Level
      float V0Y_Gamma = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);
      float V0Y_Lambda = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassLambda);
      float V0Y_K0Short = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassK0Short);

      float pT = v0.pt();
      histos.fill(HIST("V0/hpT"), pT);
      histos.fill(HIST("V0/h2dArmenterosP"), v0.alpha(), v0.qtarm());
      histos.fill(HIST("V0/hRadius"), v0.v0radius());
      histos.fill(HIST("V0/hZ"), v0.z());
      histos.fill(HIST("V0/hCosPA"), v0.v0cosPA());
      histos.fill(HIST("V0/hdcaDau"), v0.dcaV0daughters());
      histos.fill(HIST("V0/hdcaNegtopv"), v0.dcanegtopv());
      histos.fill(HIST("V0/hdcaPostopv"), v0.dcapostopv());
      histos.fill(HIST("V0/h2dEtaPhi"), v0.eta(), RecoDecay::phi(v0.px(), v0.py()));
      histos.fill(HIST("V0/hYGamma"), V0Y_Gamma);
      histos.fill(HIST("V0/hYLambda"), V0Y_Lambda);
      histos.fill(HIST("V0/hYK0Short"), V0Y_K0Short);
      histos.fill(HIST("V0/hMassGamma"), v0.mGamma());
      histos.fill(HIST("V0/hMassLambda"), v0.mLambda());
      histos.fill(HIST("V0/hMassK0Short"), v0.mK0Short());
      histos.fill(HIST("V0/hV0Type"), v0.v0Type());
      histos.fill(HIST("V0/h2dV0Indices"), v0.straCollisionId(), coll.globalIndex()); // cross-check index correctness

      // Track-level
      auto posTrack = v0.template posTrackExtra_as<dauTracks>();
      auto negTrack = v0.template negTrackExtra_as<dauTracks>();

      uint8_t positiveTrackCode = ((uint8_t(posTrack.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrack.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrack.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrack.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrack.hasTOF()) << hasTOF));

      uint8_t negativeTrackCode = ((uint8_t(negTrack.hasTPC()) << hasTPC) |
                                   (uint8_t(negTrack.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(negTrack.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(negTrack.hasTRD()) << hasTRD) |
                                   (uint8_t(negTrack.hasTOF()) << hasTOF));

      histos.fill(HIST("V0/Track/h2dITSNCls"), v0.positivept(), posTrack.itsNCls());
      histos.fill(HIST("V0/Track/h2dITSNCls"), -1 * v0.negativept(), negTrack.itsNCls());
      histos.fill(HIST("V0/Track/h2dITSChi2PerNcl"), v0.positivept(), posTrack.itsChi2PerNcl());
      histos.fill(HIST("V0/Track/h2dITSChi2PerNcl"), -1 * v0.negativept(), negTrack.itsChi2PerNcl());
      histos.fill(HIST("V0/Track/h2dTPCCrossedRows"), v0.positivept(), posTrack.tpcCrossedRows());
      histos.fill(HIST("V0/Track/h2dTPCCrossedRows"), -1 * v0.negativept(), negTrack.tpcCrossedRows());
      histos.fill(HIST("V0/Track/hTrackCode"), positiveTrackCode);                                    // pos track info
      histos.fill(HIST("V0/Track/hTrackCode"), negativeTrackCode);                                    // neg track info
      histos.fill(HIST("V0/Track/h3dTrackPropertiesVspT"), positiveTrackCode, negativeTrackCode, pT); // tracking complete info
      histos.fill(HIST("V0/Track/h2dPosTrackProperties"), positiveTrackCode, v0.positivept());        // pos track info
      histos.fill(HIST("V0/Track/h2dNegTrackProperties"), negativeTrackCode, v0.negativept());        // neg track info

      // PID (TPC)
      histos.fill(HIST("V0/PID/h2dTPCNSigmaEl"), v0.positivept(), posTrack.tpcNSigmaEl());
      histos.fill(HIST("V0/PID/h2dTPCNSigmaEl"), -1 * v0.negativept(), negTrack.tpcNSigmaEl());
      histos.fill(HIST("V0/PID/h2dTPCNSigmaPr"), v0.positivept(), posTrack.tpcNSigmaPr());
      histos.fill(HIST("V0/PID/h2dTPCNSigmaPr"), -1 * v0.negativept(), negTrack.tpcNSigmaPr());
      histos.fill(HIST("V0/PID/h2dTPCNSigmaPi"), v0.positivept(), posTrack.tpcNSigmaPi());
      histos.fill(HIST("V0/PID/h2dTPCNSigmaPi"), -1 * v0.negativept(), negTrack.tpcNSigmaPi());
      histos.fill(HIST("V0/PID/h2dTPCSignal"), v0.positivept(), posTrack.tpcSignal());
      histos.fill(HIST("V0/PID/h2dTPCSignal"), -1 * v0.negativept(), negTrack.tpcSignal());

      // PID (TOF)
      histos.fill(HIST("V0/PID/h2dTOFNSigmaLaPr"), pT, v0.tofNSigmaLaPr());
      histos.fill(HIST("V0/PID/h2dTOFNSigmaLaPi"), pT, v0.tofNSigmaLaPi());
      histos.fill(HIST("V0/PID/h2dposTOFDeltaTLaPr"), pT, v0.posTOFDeltaTLaPr());
      histos.fill(HIST("V0/PID/h2dnegTOFDeltaTLaPi"), pT, v0.negTOFDeltaTLaPi());
      histos.fill(HIST("V0/PID/h2dnegTOFDeltaTLaPr"), pT, v0.negTOFDeltaTLaPr());
      histos.fill(HIST("V0/PID/h2dposTOFDeltaTLaPi"), pT, v0.posTOFDeltaTLaPi());
      histos.fill(HIST("V0/PID/h2dTOFNSigmaALaPr"), pT, v0.tofNSigmaALaPr());
      histos.fill(HIST("V0/PID/h2dTOFNSigmaALaPi"), pT, v0.tofNSigmaALaPi());

      histos.fill(HIST("V0/PID/h2dTOFNSigmaK0PiPlus"), pT, v0.tofNSigmaK0PiPlus());
      histos.fill(HIST("V0/PID/h2dTOFNSigmaK0PiMinus"), pT, v0.tofNSigmaK0PiMinus());

      // PID TPC + TOF
      histos.fill(HIST("V0/PID/h3dTPCVsTOFNSigmaLaPr"), posTrack.tpcNSigmaPr(), v0.tofNSigmaLaPr(), v0.positivept());
      histos.fill(HIST("V0/PID/h3dTPCVsTOFNSigmaLaPi"), negTrack.tpcNSigmaPi(), v0.tofNSigmaLaPi(), v0.negativept());
    }
  }

  void processMCDerivedV0s(V0DerivedMCDatas const& V0s, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    for (auto const& v0 : V0s) {
      histos.fill(HIST("MCV0/hv0MCCore"), v0.has_v0MCCore());
      if (!v0.has_v0MCCore())
        continue;

      auto v0MC = v0.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

      // General
      histos.fill(HIST("MCV0/h2dPDGV0VsMother"), v0MC.pdgCode(), v0MC.pdgCodeMother());
      histos.fill(HIST("MCV0/h2dPDGV0VsPositive"), v0MC.pdgCode(), v0MC.pdgCodePositive());
      histos.fill(HIST("MCV0/h2dPDGV0VsNegative"), v0MC.pdgCode(), v0MC.pdgCodeNegative());
      histos.fill(HIST("MCV0/h2dPDGV0VsIsPhysicalPrimary"), v0MC.pdgCode(), v0MC.isPhysicalPrimary());

      // Track-level
      auto posTrack = v0.template posTrackExtra_as<dauTracks>();
      auto negTrack = v0.template negTrackExtra_as<dauTracks>();

      // Specific analysis by species:
      if (v0MC.pdgCode() == 22) { // IsGamma
        histos.fill(HIST("MCV0/h2dArmenterosP"), v0.alpha(), v0.qtarm());
        histos.fill(HIST("MCV0/Gamma/h2dpTResolution"), v0.pt(), v0.pt() - v0MC.ptMC());
        histos.fill(HIST("MCV0/Gamma/h2dMass"), v0.pt(), v0.mGamma());
        histos.fill(HIST("MCV0/Gamma/h2dTPCNSigmaEl"), v0.positivept(), posTrack.tpcNSigmaEl());
        histos.fill(HIST("MCV0/Gamma/h2dTPCNSigmaEl"), -1 * v0.negativept(), negTrack.tpcNSigmaEl());
        histos.fill(HIST("MCV0/Gamma/h2dTPCSignal"), v0.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCV0/Gamma/h2dTPCSignal"), -1 * v0.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCV0/Gamma/hRadius"), v0.v0radius());
        histos.fill(HIST("MCV0/Gamma/hCosPA"), v0.v0cosPA());
        histos.fill(HIST("MCV0/Gamma/hdcaDau"), v0.dcaV0daughters());
        histos.fill(HIST("MCV0/Gamma/hdcaNegtopv"), v0.dcanegtopv());
        histos.fill(HIST("MCV0/Gamma/hdcaPostopv"), v0.dcapostopv());
      }
      if (v0MC.pdgCode() == 3122) { // IsLambda
        histos.fill(HIST("MCV0/h2dArmenterosP"), v0.alpha(), v0.qtarm());
        histos.fill(HIST("MCV0/Lambda/h2dpTResolution"), v0.pt(), v0.pt() - v0MC.ptMC());
        histos.fill(HIST("MCV0/Lambda/h2dMass"), v0.pt(), v0.mLambda());
        histos.fill(HIST("MCV0/Lambda/h2dTPCNSigmaPr"), v0.positivept(), posTrack.tpcNSigmaPr());
        histos.fill(HIST("MCV0/Lambda/h2dTPCNSigmaPi"), v0.negativept(), negTrack.tpcNSigmaPi());
        histos.fill(HIST("MCV0/Lambda/h2dTPCSignal"), v0.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCV0/Lambda/h2dTPCSignal"), -1 * v0.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCV0/Lambda/hRadius"), v0.v0radius());
        histos.fill(HIST("MCV0/Lambda/hCosPA"), v0.v0cosPA());
        histos.fill(HIST("MCV0/Lambda/hdcaDau"), v0.dcaV0daughters());
        histos.fill(HIST("MCV0/Lambda/hdcaNegtopv"), v0.dcanegtopv());
        histos.fill(HIST("MCV0/Lambda/hdcaPostopv"), v0.dcapostopv());
      }
      if (v0MC.pdgCode() == -3122) { // IsAntiLambda
        histos.fill(HIST("MCV0/h2dArmenterosP"), v0.alpha(), v0.qtarm());
        histos.fill(HIST("MCV0/AntiLambda/h2dpTResolution"), v0.pt(), v0.pt() - v0MC.ptMC());
        histos.fill(HIST("MCV0/AntiLambda/h2dMass"), v0.pt(), v0.mAntiLambda());
        histos.fill(HIST("MCV0/AntiLambda/h2dTPCNSigmaPr"), v0.negativept(), negTrack.tpcNSigmaPr());
        histos.fill(HIST("MCV0/AntiLambda/h2dTPCNSigmaPi"), v0.positivept(), posTrack.tpcNSigmaPi());
        histos.fill(HIST("MCV0/AntiLambda/h2dTPCSignal"), v0.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCV0/AntiLambda/h2dTPCSignal"), -1 * v0.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCV0/AntiLambda/hRadius"), v0.v0radius());
        histos.fill(HIST("MCV0/AntiLambda/hCosPA"), v0.v0cosPA());
        histos.fill(HIST("MCV0/AntiLambda/hdcaDau"), v0.dcaV0daughters());
        histos.fill(HIST("MCV0/AntiLambda/hdcaNegtopv"), v0.dcanegtopv());
        histos.fill(HIST("MCV0/AntiLambda/hdcaPostopv"), v0.dcapostopv());
      }
      if (v0MC.pdgCode() == 310) { // IsK0Short
        histos.fill(HIST("MCV0/h2dArmenterosP"), v0.alpha(), v0.qtarm());
        histos.fill(HIST("MCV0/K0Short/h2dpTResolution"), v0.pt(), v0.pt() - v0MC.ptMC());
        histos.fill(HIST("MCV0/K0Short/h2dMass"), v0.pt(), v0.mK0Short());
        histos.fill(HIST("MCV0/K0Short/h2dTPCNSigmaPi"), v0.positivept(), posTrack.tpcNSigmaPi());
        histos.fill(HIST("MCV0/K0Short/h2dTPCNSigmaPi"), -1 * v0.negativept(), negTrack.tpcNSigmaPi());
        histos.fill(HIST("MCV0/K0Short/h2dTPCSignal"), v0.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCV0/K0Short/h2dTPCSignal"), -1 * v0.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCV0/K0Short/hRadius"), v0.v0radius());
        histos.fill(HIST("MCV0/K0Short/hCosPA"), v0.v0cosPA());
        histos.fill(HIST("MCV0/K0Short/hdcaDau"), v0.dcaV0daughters());
        histos.fill(HIST("MCV0/K0Short/hdcaNegtopv"), v0.dcanegtopv());
        histos.fill(HIST("MCV0/K0Short/hdcaPostopv"), v0.dcapostopv());
      }
    }
  }

  void processDerivedCascades(StrCollisionsDatas::iterator const& coll, CascDerivedDatas const& Cascades, dauTracks const&)
  {
    for (auto& casc : Cascades) {
      // Cascade level
      float pT = casc.pt();
      histos.fill(HIST("Casc/Sign"), casc.sign());
      histos.fill(HIST("Casc/hpT"), pT);
      histos.fill(HIST("Casc/hV0Radius"), casc.v0radius());
      histos.fill(HIST("Casc/hCascRadius"), casc.cascradius());
      histos.fill(HIST("Casc/hV0CosPA"), casc.v0cosPA(casc.x(), casc.y(), casc.z()));
      histos.fill(HIST("Casc/hCascCosPA"), casc.casccosPA(casc.x(), casc.y(), casc.z()));
      histos.fill(HIST("Casc/hDCAPosToPV"), casc.dcapostopv());
      histos.fill(HIST("Casc/hDCANegToPV"), casc.dcanegtopv());
      histos.fill(HIST("Casc/hDCABachToPV"), casc.dcabachtopv());
      histos.fill(HIST("Casc/hDCAXYCascToPV"), casc.dcaXYCascToPV());
      histos.fill(HIST("Casc/hDCAZCascToPV"), casc.dcaZCascToPV());
      histos.fill(HIST("Casc/hDCAV0ToPV"), casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()));
      histos.fill(HIST("Casc/hDCAV0Dau"), casc.dcaV0daughters());
      histos.fill(HIST("Casc/hDCACascDau"), casc.dcacascdaughters());
      histos.fill(HIST("Casc/hLambdaMass"), casc.mLambda());

      // Track level
      auto negTrack = casc.template negTrackExtra_as<dauTracks>();
      auto posTrack = casc.template posTrackExtra_as<dauTracks>();
      auto bachTrack = casc.template bachTrackExtra_as<dauTracks>();

      uint8_t positiveTrackCode = ((uint8_t(posTrack.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrack.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrack.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrack.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrack.hasTOF()) << hasTOF));

      uint8_t negativeTrackCode = ((uint8_t(negTrack.hasTPC()) << hasTPC) |
                                   (uint8_t(negTrack.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(negTrack.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(negTrack.hasTRD()) << hasTRD) |
                                   (uint8_t(negTrack.hasTOF()) << hasTOF));

      uint8_t bachTrackCode = ((uint8_t(bachTrack.hasTPC()) << hasTPC) |
                               (uint8_t(bachTrack.hasITSTracker()) << hasITSTracker) |
                               (uint8_t(bachTrack.hasITSAfterburner()) << hasITSAfterburner) |
                               (uint8_t(bachTrack.hasTRD()) << hasTRD) |
                               (uint8_t(bachTrack.hasTOF()) << hasTOF));

      histos.fill(HIST("Casc/Track/h3dTrackProperties"), positiveTrackCode, negativeTrackCode, bachTrackCode); // complete tracking info
      histos.fill(HIST("Casc/Track/h2dPosTrackProperties"), positiveTrackCode, casc.positivept());             // positive track info
      histos.fill(HIST("Casc/Track/h2dNegTrackProperties"), negativeTrackCode, casc.negativept());             // negative track info
      histos.fill(HIST("Casc/Track/h2dBachTrackProperties"), bachTrackCode, casc.bachelorpt());                // bach track info
      histos.fill(HIST("Casc/Track/h2dV0ITSChi2PerNcl"), casc.positivept(), posTrack.itsChi2PerNcl());
      histos.fill(HIST("Casc/Track/h2dV0ITSChi2PerNcl"), -1 * casc.negativept(), negTrack.itsChi2PerNcl());
      histos.fill(HIST("Casc/Track/h2dV0TPCCrossedRows"), casc.positivept(), posTrack.tpcCrossedRows());
      histos.fill(HIST("Casc/Track/h2dV0TPCCrossedRows"), -1 * casc.negativept(), negTrack.tpcCrossedRows());
      histos.fill(HIST("Casc/Track/h2dV0ITSNCls"), casc.positivept(), posTrack.itsNCls());
      histos.fill(HIST("Casc/Track/h2dV0ITSNCls"), -1 * casc.negativept(), negTrack.itsNCls());

      // PID (TPC)
      histos.fill(HIST("Casc/PID/h2dV0TPCNSigmaPr"), casc.positivept(), posTrack.tpcNSigmaPr());
      histos.fill(HIST("Casc/PID/h2dV0TPCNSigmaPr"), -1 * casc.negativept(), negTrack.tpcNSigmaPr());
      histos.fill(HIST("Casc/PID/h2dV0TPCNSigmaPi"), casc.positivept(), posTrack.tpcNSigmaPi());
      histos.fill(HIST("Casc/PID/h2dV0TPCNSigmaPi"), -1 * casc.negativept(), negTrack.tpcNSigmaPi());
      histos.fill(HIST("Casc/PID/h2dV0TPCSignal"), casc.positivept(), posTrack.tpcSignal());
      histos.fill(HIST("Casc/PID/h2dV0TPCSignal"), -1 * casc.negativept(), negTrack.tpcSignal());

      // PID (TOF)
      histos.fill(HIST("Casc/PID/h2dTOFNSigmaXiLaPi"), pT, casc.tofNSigmaXiLaPi()); //! meson track NSigma from pion <- lambda <- xi expectation
      histos.fill(HIST("Casc/PID/h2dTOFNSigmaXiLaPr"), pT, casc.tofNSigmaXiLaPr()); //! baryon track NSigma from proton <- lambda <- xi expectation
      histos.fill(HIST("Casc/PID/h2dTOFNSigmaXiPi"), pT, casc.tofNSigmaXiPi());     //! bachelor track NSigma from pion <- xi expectation
      histos.fill(HIST("Casc/PID/h2dTOFNSigmaOmLaPi"), pT, casc.tofNSigmaOmLaPi()); //! meson track NSigma from pion <- lambda <- om expectation
      histos.fill(HIST("Casc/PID/h2dTOFNSigmaOmLaPr"), pT, casc.tofNSigmaOmLaPr()); //! baryon track NSigma from proton <- lambda <- om expectation
      histos.fill(HIST("Casc/PID/h2dTOFNSigmaOmKa"), pT, casc.tofNSigmaOmKa());     //! bachelor track NSigma from kaon <- om expectation

      // By particle species
      if (casc.sign() < 0) {
        histos.fill(HIST("Casc/hMassXiMinus"), casc.mXi());
        histos.fill(HIST("Casc/hMassOmegaMinus"), casc.mOmega());
        histos.fill(HIST("Casc/Track/h2dBachITSNCls"), -1 * casc.bachelorpt(), bachTrack.itsNCls());
        histos.fill(HIST("Casc/Track/h2dBachITSChi2PerNcl"), -1 * casc.bachelorpt(), bachTrack.itsChi2PerNcl());
        histos.fill(HIST("Casc/Track/h2dBachTPCCrossedRows"), -1 * casc.bachelorpt(), bachTrack.tpcCrossedRows());
        histos.fill(HIST("Casc/PID/h2dBachTPCSignal"), -1 * casc.bachelorpt(), bachTrack.tpcSignal());
      } else {
        histos.fill(HIST("Casc/hMassXiPlus"), casc.mXi());
        histos.fill(HIST("Casc/hMassOmegaPlus"), casc.mOmega());
        histos.fill(HIST("Casc/Track/h2dBachITSNCls"), casc.bachelorpt(), bachTrack.itsNCls());
        histos.fill(HIST("Casc/Track/h2dBachITSChi2PerNcl"), casc.bachelorpt(), bachTrack.itsChi2PerNcl());
        histos.fill(HIST("Casc/Track/h2dBachTPCCrossedRows"), casc.bachelorpt(), bachTrack.tpcCrossedRows());
        histos.fill(HIST("Casc/PID/h2dBachTPCSignal"), casc.bachelorpt(), bachTrack.tpcSignal());
      }
    }
  }

  void processMCDerivedCascades(StrCollisionsDatas::iterator const& coll, CascDerivedMCDatas const& Cascades, dauTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {
    for (auto& casc : Cascades) {

      float pT = casc.pt();
      histos.fill(HIST("MCCasc/hcascMCCore"), casc.has_cascMCCore());
      if (!casc.has_cascMCCore())
        continue;
      auto cascMC = casc.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

      // General
      histos.fill(HIST("MCCasc/h2dPDGV0VsMother"), cascMC.pdgCode(), cascMC.pdgCodeMother());
      histos.fill(HIST("MCCasc/h2dPDGV0VsPositive"), cascMC.pdgCode(), cascMC.pdgCodePositive());
      histos.fill(HIST("MCCasc/h2dPDGV0VsNegative"), cascMC.pdgCode(), cascMC.pdgCodeNegative());
      histos.fill(HIST("MCCasc/h2dPDGV0VsBach"), cascMC.pdgCode(), cascMC.pdgCodeBachelor());
      histos.fill(HIST("MCCasc/h2dPDGV0VsIsPhysicalPrimary"), cascMC.pdgCode(), cascMC.isPhysicalPrimary());

      // Track level
      auto negTrack = casc.template negTrackExtra_as<dauTracks>();
      auto posTrack = casc.template posTrackExtra_as<dauTracks>();
      auto bachTrack = casc.template bachTrackExtra_as<dauTracks>();

      // Specific analysis by species:
      if (cascMC.pdgCode() == 3312) { // XiMinus
        histos.fill(HIST("MCCasc/XiMinus/h2dpTResolution"), pT, pT - cascMC.ptMC());
        histos.fill(HIST("MCCasc/XiMinus/h2dMass"), pT, casc.mXi());
        histos.fill(HIST("MCCasc/XiMinus/h2dV0TPCSignal"), casc.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCCasc/XiMinus/h2dV0TPCSignal"), -1 * casc.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCCasc/XiMinus/h2dBachTPCSignal"), casc.bachelorpt(), bachTrack.tpcSignal());
        histos.fill(HIST("MCCasc/XiMinus/hV0Radius"), casc.v0radius());
        histos.fill(HIST("MCCasc/XiMinus/hCascRadius"), casc.cascradius());
        histos.fill(HIST("MCCasc/XiMinus/hV0CosPA"), casc.v0cosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/XiMinus/hCascCosPA"), casc.casccosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/XiMinus/hDCAPosToPV"), casc.dcapostopv());
        histos.fill(HIST("MCCasc/XiMinus/hDCANegToPV"), casc.dcanegtopv());
        histos.fill(HIST("MCCasc/XiMinus/hDCABachToPV"), casc.dcabachtopv());
        histos.fill(HIST("MCCasc/XiMinus/hDCAXYCascToPV"), casc.dcaXYCascToPV());
        histos.fill(HIST("MCCasc/XiMinus/hDCAZCascToPV"), casc.dcaZCascToPV());
        histos.fill(HIST("MCCasc/XiMinus/hDCAV0ToPV"), casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()));
        histos.fill(HIST("MCCasc/XiMinus/hDCAV0Dau"), casc.dcaV0daughters());
        histos.fill(HIST("MCCasc/XiMinus/hDCACascDau"), casc.dcacascdaughters());
        histos.fill(HIST("MCCasc/XiMinus/hLambdaMass"), casc.mLambda());
      }
      if (cascMC.pdgCode() == -3312) { // XiPlus
        histos.fill(HIST("MCCasc/XiPlus/h2dpTResolution"), pT, pT - cascMC.ptMC());
        histos.fill(HIST("MCCasc/XiPlus/h2dMass"), pT, casc.mXi());
        histos.fill(HIST("MCCasc/XiPlus/h2dV0TPCSignal"), casc.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCCasc/XiPlus/h2dV0TPCSignal"), -1 * casc.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCCasc/XiPlus/h2dBachTPCSignal"), casc.bachelorpt(), bachTrack.tpcSignal());
        histos.fill(HIST("MCCasc/XiPlus/hV0Radius"), casc.v0radius());
        histos.fill(HIST("MCCasc/XiPlus/hCascRadius"), casc.cascradius());
        histos.fill(HIST("MCCasc/XiPlus/hV0CosPA"), casc.v0cosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/XiPlus/hCascCosPA"), casc.casccosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/XiPlus/hDCAPosToPV"), casc.dcapostopv());
        histos.fill(HIST("MCCasc/XiPlus/hDCANegToPV"), casc.dcanegtopv());
        histos.fill(HIST("MCCasc/XiPlus/hDCABachToPV"), casc.dcabachtopv());
        histos.fill(HIST("MCCasc/XiPlus/hDCAXYCascToPV"), casc.dcaXYCascToPV());
        histos.fill(HIST("MCCasc/XiPlus/hDCAZCascToPV"), casc.dcaZCascToPV());
        histos.fill(HIST("MCCasc/XiPlus/hDCAV0ToPV"), casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()));
        histos.fill(HIST("MCCasc/XiPlus/hDCAV0Dau"), casc.dcaV0daughters());
        histos.fill(HIST("MCCasc/XiPlus/hDCACascDau"), casc.dcacascdaughters());
        histos.fill(HIST("MCCasc/XiPlus/hLambdaMass"), casc.mLambda());
      }
      if (cascMC.pdgCode() == 3334) { // OmegaMinus
        histos.fill(HIST("MCCasc/OmegaMinus/h2dpTResolution"), pT, pT - cascMC.ptMC());
        histos.fill(HIST("MCCasc/OmegaMinus/h2dMass"), pT, casc.mOmega());
        histos.fill(HIST("MCCasc/OmegaMinus/h2dV0TPCSignal"), casc.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCCasc/OmegaMinus/h2dV0TPCSignal"), -1 * casc.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCCasc/OmegaMinus/h2dBachTPCSignal"), casc.bachelorpt(), bachTrack.tpcSignal());
        histos.fill(HIST("MCCasc/OmegaMinus/hV0Radius"), casc.v0radius());
        histos.fill(HIST("MCCasc/OmegaMinus/hCascRadius"), casc.cascradius());
        histos.fill(HIST("MCCasc/OmegaMinus/hV0CosPA"), casc.v0cosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/OmegaMinus/hCascCosPA"), casc.casccosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/OmegaMinus/hDCAPosToPV"), casc.dcapostopv());
        histos.fill(HIST("MCCasc/OmegaMinus/hDCANegToPV"), casc.dcanegtopv());
        histos.fill(HIST("MCCasc/OmegaMinus/hDCABachToPV"), casc.dcabachtopv());
        histos.fill(HIST("MCCasc/OmegaMinus/hDCAXYCascToPV"), casc.dcaXYCascToPV());
        histos.fill(HIST("MCCasc/OmegaMinus/hDCAZCascToPV"), casc.dcaZCascToPV());
        histos.fill(HIST("MCCasc/OmegaMinus/hDCAV0ToPV"), casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()));
        histos.fill(HIST("MCCasc/OmegaMinus/hDCAV0Dau"), casc.dcaV0daughters());
        histos.fill(HIST("MCCasc/OmegaMinus/hDCACascDau"), casc.dcacascdaughters());
        histos.fill(HIST("MCCasc/OmegaMinus/hLambdaMass"), casc.mLambda());
      }
      if (cascMC.pdgCode() == -3334) { // OmegaPlus
        histos.fill(HIST("MCCasc/OmegaPlus/h2dpTResolution"), pT, pT - cascMC.ptMC());
        histos.fill(HIST("MCCasc/OmegaPlus/h2dMass"), pT, casc.mOmega());
        histos.fill(HIST("MCCasc/OmegaPlus/h2dV0TPCSignal"), casc.positivept(), posTrack.tpcSignal());
        histos.fill(HIST("MCCasc/OmegaPlus/h2dV0TPCSignal"), -1 * casc.negativept(), negTrack.tpcSignal());
        histos.fill(HIST("MCCasc/OmegaPlus/h2dBachTPCSignal"), casc.bachelorpt(), bachTrack.tpcSignal());
        histos.fill(HIST("MCCasc/OmegaPlus/hV0Radius"), casc.v0radius());
        histos.fill(HIST("MCCasc/OmegaPlus/hCascRadius"), casc.cascradius());
        histos.fill(HIST("MCCasc/OmegaPlus/hV0CosPA"), casc.v0cosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/OmegaPlus/hCascCosPA"), casc.casccosPA(casc.x(), casc.y(), casc.z()));
        histos.fill(HIST("MCCasc/OmegaPlus/hDCAPosToPV"), casc.dcapostopv());
        histos.fill(HIST("MCCasc/OmegaPlus/hDCANegToPV"), casc.dcanegtopv());
        histos.fill(HIST("MCCasc/OmegaPlus/hDCABachToPV"), casc.dcabachtopv());
        histos.fill(HIST("MCCasc/OmegaPlus/hDCAXYCascToPV"), casc.dcaXYCascToPV());
        histos.fill(HIST("MCCasc/OmegaPlus/hDCAZCascToPV"), casc.dcaZCascToPV());
        histos.fill(HIST("MCCasc/OmegaPlus/hDCAV0ToPV"), casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()));
        histos.fill(HIST("MCCasc/OmegaPlus/hDCAV0Dau"), casc.dcaV0daughters());
        histos.fill(HIST("MCCasc/OmegaPlus/hDCACascDau"), casc.dcacascdaughters());
        histos.fill(HIST("MCCasc/OmegaPlus/hLambdaMass"), casc.mLambda());
      }
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void processGenerated(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const& V0MCCores, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const& CascMCCores, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels> const& collisions)
  {
    fillGeneratedEventProperties(mcCollisions, collisions);
    std::vector<int> listBestCollisionIdx = getListOfRecoCollIndices(mcCollisions, collisions);
    for (auto const& v0MC : V0MCCores) {
      if (!v0MC.has_straMCCollision())
        continue;

      if (!v0MC.isPhysicalPrimary())
        continue;

      float ptmc = v0MC.ptMC();
      float ymc = 1e3;
      if (v0MC.pdgCode() == 310)
        ymc = v0MC.rapidityMC(0);
      else if (TMath::Abs(v0MC.pdgCode()) == 3122)
        ymc = v0MC.rapidityMC(1);

      if (TMath::Abs(ymc) > v0Selections.rapidityCut)
        continue;

      auto mcCollision = v0MC.straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
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
        float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();

        if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
          continue;
        }
        if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
          continue;
        }

        if (v0MC.pdgCode() == 310) {
          histos.fill(HIST("GenMC/h2dGenK0ShortVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (v0MC.pdgCode() == 3122) {
          histos.fill(HIST("GenMC/h2dGenLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (v0MC.pdgCode() == -3122) {
          histos.fill(HIST("GenMC/h2dGenAntiLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
      }
      if (v0MC.pdgCode() == 22) {
        histos.fill(HIST("GenMC/h2dGenGamma"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenGammaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (v0MC.pdgCode() == 310) {
        histos.fill(HIST("GenMC/h2dGenK0Short"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenK0ShortVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (v0MC.pdgCode() == 3122) {
        histos.fill(HIST("GenMC/h2dGenLambda"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (v0MC.pdgCode() == -3122) {
        histos.fill(HIST("GenMC/h2dGenAntiLambda"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenAntiLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    }

    for (auto const& cascMC : CascMCCores) {
      if (!cascMC.has_straMCCollision())
        continue;

      if (!cascMC.isPhysicalPrimary())
        continue;

      float ptmc = cascMC.ptMC();
      float ymc = 1e3;
      if (TMath::Abs(cascMC.pdgCode()) == 3312)
        ymc = cascMC.rapidityMC(0);
      else if (TMath::Abs(cascMC.pdgCode()) == 3334)
        ymc = cascMC.rapidityMC(2);

      if (TMath::Abs(ymc) > v0Selections.rapidityCut)
        continue;

      auto mcCollision = cascMC.straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
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
        float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();

        if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
          continue;
        }
        if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
          continue;
        }

        if (cascMC.pdgCode() == 3312) {
          histos.fill(HIST("GenMC/h2dGenXiMinusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (cascMC.pdgCode() == -3312) {
          histos.fill(HIST("GenMC/h2dGenXiPlusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (cascMC.pdgCode() == 3334) {
          histos.fill(HIST("GenMC/h2dGenOmegaMinusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (cascMC.pdgCode() == -3334) {
          histos.fill(HIST("GenMC/h2dGenOmegaPlusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
      }

      if (cascMC.pdgCode() == 3312) {
        histos.fill(HIST("GenMC/h2dGenXiMinus"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenXiMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == -3312) {
        histos.fill(HIST("GenMC/h2dGenXiPlus"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenXiPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == 3334) {
        histos.fill(HIST("GenMC/h2dGenOmegaMinus"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenOmegaMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (cascMC.pdgCode() == -3334) {
        histos.fill(HIST("GenMC/h2dGenOmegaPlus"), centrality, ptmc);
        histos.fill(HIST("GenMC/h2dGenOmegaPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    }
  }

  PROCESS_SWITCH(strderivedGenQA, processDerivedV0s, "Process derived data", true);
  PROCESS_SWITCH(strderivedGenQA, processMCDerivedV0s, "Process derived data", false);
  PROCESS_SWITCH(strderivedGenQA, processDerivedCascades, "Process derived data", true);
  PROCESS_SWITCH(strderivedGenQA, processMCDerivedCascades, "Process derived data", false);
  PROCESS_SWITCH(strderivedGenQA, processGenerated, "process MC generated", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strderivedGenQA>(cfgc)};
}
