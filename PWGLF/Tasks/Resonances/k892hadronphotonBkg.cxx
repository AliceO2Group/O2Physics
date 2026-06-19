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
// This is a ask that computes the same-event rotational and the
// mixed-event combinatorial backgrounds for the K*(892) -> K0S + gamma analysis.
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  K892 hadron-photon background task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    oussama.benchikhi@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <TH1.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0StandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;

struct k892hadronphotonBkg {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;
  TRandom3 rotRng{12345}; // struct member; fixed seed for reproducibility across grid jobs

  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};

  // For ML Selection
  Configurable<bool> useMLScores{"useMLScores", false, "use ML scores to select candidates"};

  // Interaction-rate retrieval (used by the event selection)
  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<bool> fIRCrashOnNull{"fIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  struct : ConfigurableGroup {
    std::string prefix = "kstarBkgConfig";
    Configurable<bool> doSameEvtRotation{"doSameEvtRotation", false, "Same-event rotational background"};
    Configurable<bool> doEvtMixing{"doEvtMixing", false, "Mixed-event background"};
    Configurable<int> nMix{"nMix", 5, "Number of mixed events"};
    Configurable<int> deltaCollision{"deltaCollision", 25, "Min |Δ globalIndex| for mixing"};
    Configurable<float> kstarMaxOPAngle{"kstarMaxOPAngle", 7.f, "Max opening angle (rad)"};
    Configurable<float> kstarMaxRap{"kstarMaxRap", 0.5f, "Max |y(K*)|"};
    Configurable<int> nBkgRot{"nBkgRot", 3, "Rotations per pair (rotational bkg)"};
    Configurable<int> rotationalCut{"rotationalCut", 10, "theta band: [pi - pi/cut, pi + pi/cut]"};
  } kstarBkgConfig;

  ConfigurableAxis axisVertexMixBkg{"axisVertexMixBkg", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "z-vertex bins for mixing"};
  ConfigurableAxis axisCentralityMixBkg{"axisCentralityMixBkg", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "centrality bins for mixing"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> fUseEventSelection{"fUseEventSelection", false, "Apply event selection cuts"};
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
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF"};
    Configurable<bool> requireINEL0{"requireINEL0", true, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};
    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};
    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    // fast check on interaction rate
    Configurable<float> minIR{"minIR", -1, "minimum IR collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR collisions"};
  } eventSelections;

  //// Photon criteria:
  struct : ConfigurableGroup {
    std::string prefix = "photonSelections"; // JSON group name
    Configurable<float> gammaMLThreshold{"gammaMLThreshold", 0.1, "Decision Threshold value to select gammas"};
    Configurable<int> photonv0TypeSel{"photonv0TypeSel", 7, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> photonMinDCADauToPv{"photonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
    Configurable<float> photonMaxDCAV0Dau{"photonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
    Configurable<int> photonMinTPCCrossedRows{"photonMinTPCCrossedRows", 30, "Min daughter TPC Crossed Rows"};
    Configurable<float> photonMinTPCNSigmas{"photonMinTPCNSigmas", -7, "Min TPC NSigmas for daughters"};
    Configurable<float> photonMaxTPCNSigmas{"photonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
    Configurable<float> photonMinRapidity{"photonMinRapidity", -0.5, "v0 min rapidity"};
    Configurable<float> photonMaxRapidity{"photonMaxRapidity", 0.5, "v0 max rapidity"};
    Configurable<float> photonDauEtaMin{"photonDauEtaMin", -0.8, "Min pseudorapidity of daughter tracks"};
    Configurable<float> photonDauEtaMax{"photonDauEtaMax", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<float> photonMinRadius{"photonMinRadius", 3.0, "Min photon conversion radius (cm)"};
    Configurable<float> photonMaxRadius{"photonMaxRadius", 115, "Max photon conversion radius (cm)"};
    Configurable<float> photonMinZ{"photonMinZ", -240, "Min photon conversion point z value (cm)"};
    Configurable<float> photonMaxZ{"photonMaxZ", 240, "Max photon conversion point z value (cm)"};
    Configurable<float> photonMaxQt{"photonMaxQt", 0.08, "Max photon qt value (AP plot) (GeV/c)"};
    Configurable<float> photonMaxAlpha{"photonMaxAlpha", 1.0, "Max photon alpha absolute value (AP plot)"};
    Configurable<float> photonMinV0cospa{"photonMinV0cospa", 0.80, "Min V0 CosPA"};
    Configurable<float> photonMaxMass{"photonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};
    Configurable<float> photonPhiMin1{"photonPhiMin1", -1, "Phi min value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> photonPhiMax1{"photonPhiMax1", -1, "Phi max value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> photonPhiMin2{"photonPhiMin2", -1, "Phi max value to reject photons, region 2 (leave negative if no selection desired)"};
    Configurable<float> photonPhiMax2{"photonPhiMax2", -1, "Phi min value to reject photons, region 2 (leave negative if no selection desired)"};
  } photonSelections;

  // KShort criteria:
  struct : ConfigurableGroup {
    std::string prefix = "kshortSelections"; // JSON group name
    Configurable<float> kshortMLThreshold{"kshortMLThreshold", 0.1, "Decision Threshold value to select kshorts"};
    Configurable<float> kshortMinDCANegToPv{"kshortMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> kshortMinDCAPosToPv{"kshortMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> kshortMaxDCAV0Dau{"kshortMaxDCAV0Dau", 2.5, "Max DCA V0 Daughters (cm)"};
    Configurable<float> kshortMinv0radius{"kshortMinv0radius", 0.0, "Min V0 radius (cm)"};
    Configurable<float> kshortMaxv0radius{"kshortMaxv0radius", 40, "Max V0 radius (cm)"};
    Configurable<float> kshortMinv0cospa{"kshortMinv0cospa", 0.95, "Min V0 CosPA"};
    Configurable<float> kshortMaxLifeTime{"kshortMaxLifeTime", 20, "Max lifetime"};
    Configurable<float> kshortWindow{"kshortWindow", 0.015, "Mass window around expected (in GeV/c2). Leave negative to disable"};
    Configurable<float> kshortMinRapidity{"kshortMinRapidity", -0.5, "v0 min rapidity"};
    Configurable<float> kshortMaxRapidity{"kshortMaxRapidity", 0.5, "v0 max rapidity"};
    Configurable<float> kshortDauEtaMin{"kshortDauEtaMin", -0.8, "Min pseudorapidity of daughter tracks"};
    Configurable<float> kshortDauEtaMax{"kshortDauEtaMax", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<float> kshortMinZ{"kshortMinZ", -240, "Min kshort decay point z value (cm)"};
    Configurable<float> kshortMaxZ{"kshortMaxZ", 240, "Max kshort decay point z value (cm)"};
    Configurable<int> kshortMinTPCCrossedRows{"kshortMinTPCCrossedRows", 50, "Min daughter TPC Crossed Rows"};
    Configurable<int> kshortMinITSclusters{"kshortMinITSclusters", 1, "minimum ITS clusters"};
    Configurable<bool> kshortRejectPosITSafterburner{"kshortRejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
    Configurable<bool> kshortRejectNegITSafterburner{"kshortRejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};
    Configurable<float> kshortArmenterosCoefficient{"kshortArmenterosCoefficient", 0.2, "Armenteros-Podolanski coefficient to reject lambdas"};
  } kshortSelections;

  struct : ConfigurableGroup {
    // base properties
    std::string prefix = "axisConfig"; // JSON group name
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
    ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
    ConfigurableAxis axisKStarMass{"axisKStarMass", {500, 0.6f, 1.6f}, "M_{K^{*}} (GeV/c^{2})"};
    ConfigurableAxis axisIRBinning{"axisIRBinning", {151, -10, 1500}, "Binning for the interaction rate (kHz)"};
  } axisConfig;

  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    histos.add("hEventCentrality", "hEventCentrality", kTH1D, {axisConfig.axisCentrality});

    if (eventSelections.fUseEventSelection) {
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
        histos.add("GeneralQA/hInteractionRate", "hInteractionRate", kTH1D, {axisConfig.axisIRBinning});
        histos.add("GeneralQA/hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2D, {axisConfig.axisCentrality, axisConfig.axisIRBinning});
      }
    }

    if (kstarBkgConfig.doSameEvtRotation || kstarBkgConfig.doEvtMixing) {
      histos.add("KStarBkg/hDeltaCollision", "hDeltaCollision", kTH1D, {{2000, -1000.f, 1000.f}});
      histos.add("KStarBkg/h2dCentralityCollPair", "h2dCentralityCollPair", kTH2D, {axisConfig.axisCentrality, axisConfig.axisCentrality});
    }
    if (kstarBkgConfig.doSameEvtRotation) {
      histos.add("KStarBkg/h2dRotKStarMassVsPt", "h2dRotKStarMassVsPt", kTH2D, {axisConfig.axisKStarMass, axisConfig.axisPt});
      histos.add("KStarBkg/h3dRotKStarMassVsPt", "h3dRotKStarMassVsPt", kTH3D, {axisConfig.axisCentrality, axisConfig.axisPt, axisConfig.axisKStarMass});
      histos.add("KStarBkg/h3dRotKStarPtVsOPAngle", "h3dRotKStarPtVsOPAngle", kTH3D, {{140, 0.f, 7.f}, axisConfig.axisPt, axisConfig.axisKStarMass});
    }
    if (kstarBkgConfig.doEvtMixing) {
      histos.add("KStarBkg/h2dMixedKStarMassVsPt", "h2dMixedKStarMassVsPt", kTH2D, {axisConfig.axisKStarMass, axisConfig.axisPt});
      histos.add("KStarBkg/h3dMixedKStarMassVsPt", "h3dMixedKStarMassVsPt", kTH3D, {axisConfig.axisCentrality, axisConfig.axisPt, axisConfig.axisKStarMass});
      histos.add("KStarBkg/h3dMixedKStarPtVsOPAngle", "h3dMixedKStarPtVsOPAngle", kTH3D, {{140, 0.f, 7.f}, axisConfig.axisPt, axisConfig.axisKStarMass});
    }
  }

  //_______________________________________________
  // Event selection (identical to the builder)
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
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 19 /* Above max IR */);

    // Fill centrality histogram after event selection
    if (fillHists)
      histos.fill(HIST("hEventCentrality"), centrality);

    return true;
  }

  //_______________________________________________
  // Process v0 photon candidate (data only, no QA fills)
  template <typename TV0Object>
  bool processPhotonCandidate(TV0Object const& gamma)
  {
    // V0 type selection
    if (gamma.v0Type() != photonSelections.photonv0TypeSel && photonSelections.photonv0TypeSel > -1)
      return false;

    float photonY = RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma);

    if (useMLScores) {
      if (gamma.gammaBDTScore() <= photonSelections.gammaMLThreshold)
        return false;

    } else {
      // Standard selection
      // Gamma basic selection criteria:
      if ((gamma.mGamma() < 0) || (gamma.mGamma() > photonSelections.photonMaxMass))
        return false;

      if ((photonY < photonSelections.photonMinRapidity) || (photonY > photonSelections.photonMaxRapidity))
        return false;

      if (gamma.negativeeta() < photonSelections.photonDauEtaMin || gamma.negativeeta() > photonSelections.photonDauEtaMax)
        return false;

      if (gamma.positiveeta() < photonSelections.photonDauEtaMin || gamma.positiveeta() > photonSelections.photonDauEtaMax)
        return false;

      if ((TMath::Abs(gamma.dcapostopv()) < photonSelections.photonMinDCADauToPv) || (TMath::Abs(gamma.dcanegtopv()) < photonSelections.photonMinDCADauToPv))
        return false;

      if (TMath::Abs(gamma.dcaV0daughters()) > photonSelections.photonMaxDCAV0Dau)
        return false;

      if ((gamma.v0radius() < photonSelections.photonMinRadius) || (gamma.v0radius() > photonSelections.photonMaxRadius))
        return false;

      if ((gamma.z() < photonSelections.photonMinZ) || (gamma.z() > photonSelections.photonMaxZ))
        return false;

      if (gamma.v0cosPA() < photonSelections.photonMinV0cospa)
        return false;

      float photonPhi = RecoDecay::phi(gamma.px(), gamma.py());
      if ((((photonPhi > photonSelections.photonPhiMin1) && (photonPhi < photonSelections.photonPhiMax1)) || ((photonPhi > photonSelections.photonPhiMin2) && (photonPhi < photonSelections.photonPhiMax2))) && ((photonSelections.photonPhiMin1 != -1) && (photonSelections.photonPhiMax1 != -1) && (photonSelections.photonPhiMin2 != -1) && (photonSelections.photonPhiMax2 != -1)))
        return false;

      if (gamma.qtarm() > photonSelections.photonMaxQt)
        return false;

      if (TMath::Abs(gamma.alpha()) > photonSelections.photonMaxAlpha)
        return false;

      auto posTrackGamma = gamma.template posTrackExtra_as<dauTracks>();
      auto negTrackGamma = gamma.template negTrackExtra_as<dauTracks>();

      if ((posTrackGamma.tpcCrossedRows() < photonSelections.photonMinTPCCrossedRows) || (negTrackGamma.tpcCrossedRows() < photonSelections.photonMinTPCCrossedRows))
        return false;

      if (((posTrackGamma.tpcNSigmaEl() < photonSelections.photonMinTPCNSigmas) || (posTrackGamma.tpcNSigmaEl() > photonSelections.photonMaxTPCNSigmas)))
        return false;

      if (((negTrackGamma.tpcNSigmaEl() < photonSelections.photonMinTPCNSigmas) || (negTrackGamma.tpcNSigmaEl() > photonSelections.photonMaxTPCNSigmas)))
        return false;
    }

    return true;
  }

  //_______________________________________________
  // Process K0Short candidate (data only, no QA fills)
  template <typename TV0Object, typename TCollision>
  bool processKShortCandidate(TV0Object const& kshort, TCollision const& collision)
  {
    // V0 type selection
    if (kshort.v0Type() != 1)
      return false;

    if (useMLScores) {
      // if (kshort.k0ShortBDTScore() <= kshortSelections.kshortMLThreshold)
      return false;

    } else {
      // KShort basic selection criteria:
      if ((TMath::Abs(kshort.mK0Short() - o2::constants::physics::MassK0Short) > kshortSelections.kshortWindow) && kshortSelections.kshortWindow > 0)
        return false;

      if ((kshort.yK0Short() < kshortSelections.kshortMinRapidity) || (kshort.yK0Short() > kshortSelections.kshortMaxRapidity))
        return false;

      if ((kshort.negativeeta() < kshortSelections.kshortDauEtaMin) || (kshort.negativeeta() > kshortSelections.kshortDauEtaMax))
        return false;

      if ((kshort.positiveeta() < kshortSelections.kshortDauEtaMin) || (kshort.positiveeta() > kshortSelections.kshortDauEtaMax))
        return false;

      if ((TMath::Abs(kshort.dcapostopv()) < kshortSelections.kshortMinDCAPosToPv) || (TMath::Abs(kshort.dcanegtopv()) < kshortSelections.kshortMinDCANegToPv))
        return false;

      if ((kshort.v0radius() < kshortSelections.kshortMinv0radius) || (kshort.v0radius() > kshortSelections.kshortMaxv0radius))
        return false;

      if ((kshort.z() < kshortSelections.kshortMinZ) || (kshort.z() > kshortSelections.kshortMaxZ))
        return false;

      if (TMath::Abs(kshort.dcaV0daughters()) > kshortSelections.kshortMaxDCAV0Dau)
        return false;

      if (kshort.qtarm() < kshortSelections.kshortArmenterosCoefficient * TMath::Abs(kshort.alpha()))
        return false;

      if (kshort.v0cosPA() < kshortSelections.kshortMinv0cospa)
        return false;

      auto posTrackKShort = kshort.template posTrackExtra_as<dauTracks>();
      auto negTrackKShort = kshort.template negTrackExtra_as<dauTracks>();

      if ((posTrackKShort.tpcCrossedRows() < kshortSelections.kshortMinTPCCrossedRows) || (negTrackKShort.tpcCrossedRows() < kshortSelections.kshortMinTPCCrossedRows))
        return false;

      // MinITSCls
      bool posIsFromAfterburner = posTrackKShort.itsChi2PerNcl() < 0;
      bool negIsFromAfterburner = negTrackKShort.itsChi2PerNcl() < 0;

      if (posTrackKShort.itsNCls() < kshortSelections.kshortMinITSclusters && (!kshortSelections.kshortRejectPosITSafterburner || posIsFromAfterburner))
        return false;
      if (negTrackKShort.itsNCls() < kshortSelections.kshortMinITSclusters && (!kshortSelections.kshortRejectNegITSafterburner || negIsFromAfterburner))
        return false;

      float fKShortLifeTime = kshort.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
      if (fKShortLifeTime > kshortSelections.kshortMaxLifeTime)
        return false;
    }
    return true;
  }

  //_______________________________________________
  // Compute same-event rotational background for K* within a single collision
  template <typename TCollision, typename TV0s>
  void calculateRotBackground(TCollision const& coll,
                              std::vector<int> const& photonIndices,
                              std::vector<int> const& kshortIndices,
                              TV0s const& fullV0s)
  {
    if (photonIndices.empty() || kshortIndices.empty())
      return;

    const float centrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();
    for (int kIdx : kshortIndices) {
      auto kshort = fullV0s.rawIteratorAt(kIdx);

      for (int pIdx : photonIndices) {
        auto photon = fullV0s.rawIteratorAt(pIdx);

        // photon as a massless 4-vector
        ROOT::Math::PtEtaPhiMVector pGamma(photon.pt(),
                                           photon.eta(),
                                           photon.phi(),
                                           0.0);

        for (int irot = 0; irot < kstarBkgConfig.nBkgRot; ++irot) {
          float theta = rotRng.Uniform(o2::constants::math::PI - o2::constants::math::PI / kstarBkgConfig.rotationalCut,
                                       o2::constants::math::PI + o2::constants::math::PI / kstarBkgConfig.rotationalCut);

          ROOT::Math::PtEtaPhiMVector kRot(kshort.pt(), kshort.eta(), kshort.phi() + theta, o2::constants::physics::MassK0Short);

          auto kstar = pGamma + kRot;

          float rapidity = RecoDecay::y(std::array{static_cast<float>(kstar.Px()),
                                                   static_cast<float>(kstar.Py()),
                                                   static_cast<float>(kstar.Pz())},
                                        o2::constants::physics::MassK0Star892);
          if (std::abs(rapidity) > kstarBkgConfig.kstarMaxRap)
            continue;

          // Opening angle between photon and rotated K0s (QA only, not used as a cut)
          double cosOA = pGamma.Vect().Dot(kRot.Vect()) / (pGamma.P() * kRot.P());
          double openAngle = std::acos(cosOA);

          histos.fill(HIST("KStarBkg/h2dRotKStarMassVsPt"), kstar.M(), kstar.Pt());
          histos.fill(HIST("KStarBkg/h3dRotKStarMassVsPt"), centrality, kstar.Pt(), kstar.M());
          histos.fill(HIST("KStarBkg/h3dRotKStarPtVsOPAngle"), openAngle, kstar.Pt(), kstar.M());
        }
      }
    }
  }

  //_______________________________________________
  // Centrality helper for the background (keeps the builder's semantics)
  template <typename TCollision>
  float getCentralityRun3Bkg(TCollision const& collision)
  {
    return doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
  }

  //_______________________________________________
  // Main: same-event rotation + event mixing for K* background (data only)
  using BkgBinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  template <typename TCollisions, typename TV0s>
  void calculateKStarBkg(TCollisions const& collisions, TV0s const& fullV0s)
  {
    // Per-collision pools of selected photon and K0s V0 indices
    std::vector<std::vector<int>> photonPool(collisions.size());
    std::vector<std::vector<int>> kshortPool(collisions.size());

    // V0 grouping by straCollisionId
    std::vector<std::vector<int>> v0grouped(collisions.size());
    for (const auto& v0 : fullV0s) {
      v0grouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }

    // ── Pass 1: populate pools using single-particle selections ──
    for (const auto& coll : collisions) {

      if (eventSelections.fUseEventSelection && !isEventAccepted(coll, true))
        continue;

      for (size_t i = 0; i < v0grouped[coll.globalIndex()].size(); i++) {
        auto v0 = fullV0s.rawIteratorAt(v0grouped[coll.globalIndex()][i]);

        if (processPhotonCandidate(v0))
          photonPool[coll.globalIndex()].push_back(v0.globalIndex());

        if (processKShortCandidate(v0, coll))
          kshortPool[coll.globalIndex()].push_back(v0.globalIndex());
      }

      // Same-event rotational background
      if (kstarBkgConfig.doSameEvtRotation) {
        calculateRotBackground(coll,
                               photonPool[coll.globalIndex()],
                               kshortPool[coll.globalIndex()],
                               fullV0s);
      }
    }

    // Event Mixing
    if (!kstarBkgConfig.doEvtMixing)
      return;

    // Build the mixing binning locally: a struct member initialized from a
    // ConfigurableAxis captures the default bins at task construction time (before
    // the framework applies JSON overrides), silently ignoring user configuration.
    BkgBinningType bkgColBinning{{axisVertexMixBkg, axisCentralityMixBkg}, true};

    for (auto& [coll1, coll2] : selfCombinations(bkgColBinning, kstarBkgConfig.nMix, -1,
                                                 collisions, collisions)) {
      if (coll1.globalIndex() == coll2.globalIndex())
        continue;

      histos.fill(HIST("KStarBkg/hDeltaCollision"),
                  coll1.globalIndex() - coll2.globalIndex());
      histos.fill(HIST("KStarBkg/h2dCentralityCollPair"),
                  getCentralityRun3Bkg(coll1), getCentralityRun3Bkg(coll2));

      if (std::abs(static_cast<int64_t>(coll1.globalIndex()) - static_cast<int64_t>(coll2.globalIndex())) < kstarBkgConfig.deltaCollision)
        continue;

      auto const& photons1 = photonPool[coll1.globalIndex()];
      auto const& kshorts1 = kshortPool[coll1.globalIndex()];
      auto const& photons2 = photonPool[coll2.globalIndex()];
      auto const& kshorts2 = kshortPool[coll2.globalIndex()];

      // K0s(coll1) × γ(coll2)
      if (!kshorts1.empty() && !photons2.empty()) {
        for (int kIdx : kshorts1) {
          auto kshort = fullV0s.rawIteratorAt(kIdx);
          float kP = std::hypot(kshort.px(), kshort.py(), kshort.pz());
          ROOT::Math::PxPyPzEVector fourMomKShort(
            kshort.px(), kshort.py(), kshort.pz(),
            std::sqrt(kP * kP +
                      o2::constants::physics::MassK0Short *
                        o2::constants::physics::MassK0Short));

          for (int pIdx : photons2) {
            auto photon = fullV0s.rawIteratorAt(pIdx);
            float pP = std::hypot(photon.px(), photon.py(), photon.pz());
            ROOT::Math::PxPyPzEVector fourMomPhoton(
              photon.px(), photon.py(), photon.pz(), pP);

            auto fourMomKStar = fourMomPhoton + fourMomKShort;

            double cosOA = fourMomPhoton.Vect().Dot(fourMomKShort.Vect()) /
                           (fourMomPhoton.P() * fourMomKShort.P());
            double openAngle = std::acos(cosOA);
            float mass = fourMomKStar.M();
            float pt = fourMomKStar.Pt();
            // Rapidity computed under the K*(892) mass hypothesis (NOT the actual pair
            // invariant mass) to match the same-event rotational background and the
            // buildKStar signal selection, so the rapidity acceptance is identical for
            // signal and all backgrounds.
            float rapidity = RecoDecay::y(std::array{static_cast<float>(fourMomKStar.Px()),
                                                     static_cast<float>(fourMomKStar.Py()),
                                                     static_cast<float>(fourMomKStar.Pz())},
                                          o2::constants::physics::MassK0Star892);

            if (openAngle > kstarBkgConfig.kstarMaxOPAngle)
              continue;
            if (std::abs(rapidity) > kstarBkgConfig.kstarMaxRap)
              continue;

            histos.fill(HIST("KStarBkg/h2dMixedKStarMassVsPt"), mass, pt);
            histos.fill(HIST("KStarBkg/h3dMixedKStarMassVsPt"),
                        getCentralityRun3Bkg(coll1), pt, mass);
            histos.fill(HIST("KStarBkg/h3dMixedKStarPtVsOPAngle"),
                        openAngle, pt, mass);
          }
        }
      }

      // γ(coll1) × K0s(coll2)
      if (!photons1.empty() && !kshorts2.empty()) {
        for (int pIdx : photons1) {
          auto photon = fullV0s.rawIteratorAt(pIdx);
          float pP = std::hypot(photon.px(), photon.py(), photon.pz());
          ROOT::Math::PxPyPzEVector fourMomPhoton(
            photon.px(), photon.py(), photon.pz(), pP);

          for (int kIdx : kshorts2) {
            auto kshort = fullV0s.rawIteratorAt(kIdx);
            float kP = std::hypot(kshort.px(), kshort.py(), kshort.pz());
            ROOT::Math::PxPyPzEVector fourMomKShort(
              kshort.px(), kshort.py(), kshort.pz(),
              std::sqrt(kP * kP +
                        o2::constants::physics::MassK0Short *
                          o2::constants::physics::MassK0Short));

            auto fourMomKStar = fourMomPhoton + fourMomKShort;

            double cosOA = fourMomPhoton.Vect().Dot(fourMomKShort.Vect()) /
                           (fourMomPhoton.P() * fourMomKShort.P());
            double openAngle = std::acos(cosOA);
            float mass = fourMomKStar.M();
            float pt = fourMomKStar.Pt();
            // Rapidity computed under the K*(892) mass hypothesis (NOT the actual pair
            // invariant mass) to match the same-event rotational background and the
            // buildKStar signal selection, so the rapidity acceptance is identical for
            // signal and all backgrounds.
            float rapidity = RecoDecay::y(std::array{static_cast<float>(fourMomKStar.Px()),
                                                     static_cast<float>(fourMomKStar.Py()),
                                                     static_cast<float>(fourMomKStar.Pz())},
                                          o2::constants::physics::MassK0Star892);

            if (openAngle > kstarBkgConfig.kstarMaxOPAngle)
              continue;
            if (std::abs(rapidity) > kstarBkgConfig.kstarMaxRap)
              continue;

            histos.fill(HIST("KStarBkg/h2dMixedKStarMassVsPt"), mass, pt);
            histos.fill(HIST("KStarBkg/h3dMixedKStarMassVsPt"),
                        getCentralityRun3Bkg(coll1), pt, mass);
            histos.fill(HIST("KStarBkg/h3dMixedKStarPtVsOPAngle"),
                        openAngle, pt, mass);
          }
        }
      }
    }
  }

  //_______________________________________________
  // Data process: same-event rotational + mixed-event K* background
  void processKStarBkg(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraEvSelExtras> const& collisions,
                       V0StandardDerivedDatas const& fullV0s,
                       dauTracks const&)
  {
    calculateKStarBkg(collisions, fullV0s);
  }

  PROCESS_SWITCH(k892hadronphotonBkg, processKStarBkg, "Compute K* same-event rotational and mixed-event background (data)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892hadronphotonBkg>(cfgc)};
}
