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
/// \file strangenessderivedbinnedinfo.cxx
/// \brief analysis task producing V0/cascade info in binned format
///
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & SMI
//
// ================
//
// This code loops over V0Cores and CascCores tables and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    romain.schotter@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

// constants
const float ctauXiPDG = 4.91;     // Xi PDG lifetime
const float ctauOmegaPDG = 2.461; // Omega PDG lifetime

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using namespace o2::aod::rctsel;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>;
using CascadeCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs>;

struct strangenessderivedbinnedinfo {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> analyseK0Short{"analyseK0Short", true, "process K0Short-like candidates"};
  Configurable<bool> analyseLambda{"analyseLambda", false, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", false, "process AntiLambda-like candidates"};
  Configurable<bool> analyseXi{"analyseXi", false, "process Xi-like candidates"};
  Configurable<bool> analyseOmega{"analyseOmega", false, "process Omega-like candidates"};
  Configurable<bool> isPP{"isPP", true, "If running on pp collision, switch it on true"};

  // for running over skimmed dataset
  Configurable<bool> doPPAnalysis{"doPPAnalysis", false, "if in pp, set to true"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "If running over skimmed data, switch it on true"};
  Configurable<std::string> cfgSkimmedTrigger{"cfgSkimmedTrigger", "fDoubleXi,fTripleXi,fQuadrupleXi", "(std::string) Comma separated list of triggers of interest"};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  struct : ConfigurableGroup {
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track (Run 3 only)"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF (Run 3 only)"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD (Run 3 only)"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF (Run 3 only)"};
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

  static constexpr float defaultSqrtScalingParameters[1][4] = {{0.1, 0.1, 0, 128}};

  // preselection options
  struct : ConfigurableGroup {
    std::string prefix = "encodingOpts";
    Configurable<bool> useSqrtEncodingForOccupancy{"useSqrtEncodingForOccupancy", false, "Store sqrt(occupancy) instead of occupancy"};
    Configurable<bool> useSqrtEncodingForRadius{"useSqrtEncodingForRadius", false, "Store sqrt(radius) instead of radius"};
    Configurable<bool> useSqrtScalingForEncodingPt{"useSqrtScalingForEncodingPt", false, "Store sqrt scaling(pT) instead of pT"};
    Configurable<LabeledArray<float>> sqrtScalingParameters{"sqrtScalingParameters", {defaultSqrtScalingParameters[0], 4, {"sigma0", "sigma1", "clampMin", "clampMax"}}, "Sqrt scaling parameters"};
  } encodingOpts;

  struct : ConfigurableGroup {
    Configurable<int> v0TypeSelection{"v0Selections.v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};

    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"v0Selections.rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"v0Selections.daughterEtaCut", 0.8, "max eta for daughters"};

    // Standard 6 topological criteria
    Configurable<float> v0cospa{"v0Selections.v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"v0Selections.dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcav0topv{"v0Selections.dcav0topv", .05, "min DCA V0 to PV (cm)"};
    Configurable<float> dcapiontopv{"v0Selections.dcapiontopv", .05, "min DCA Pion To PV (cm)"};
    Configurable<float> dcaprotontopv{"v0Selections.dcaprotontopv", .05, "min DCA Proton To PV (cm)"};
    Configurable<float> v0radius{"v0Selections.v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0Selections.v0radiusMax", 1E5, "maximum V0 radius (cm)"};

    // invariant mass selection
    Configurable<float> v0MassWindow{"v0Selections.v0MassWindow", 0.008, "#Lambda mass (GeV/#it{c}^{2})"};
    Configurable<float> compMassRejection{"v0Selections.compMassRejection", 0.008, "Competing mass rejection (GeV/#it{c}^{2})"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"v0Selections.armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // Track quality
    Configurable<int> minTPCrows{"v0Selections.minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"v0Selections.minITSclusters", -1, "minimum ITS clusters"};
    Configurable<bool> requireTPConly{"v0Selections.requireTPConly", false, "require V0s comprised of at least one TPC only prong"};
    Configurable<bool> requirePosITSonly{"v0Selections.requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireNegITSonly{"v0Selections.requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};

    // PID (TPC)
    Configurable<float> tpcPidNsigmaCut{"v0Selections.tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
  } v0Selections;

  struct : ConfigurableGroup {
    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"cascSelections.rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"cascSelections.daughterEtaCut", 0.8, "max eta for daughters"};

    // Standard 6 topological criteria on V0
    Configurable<float> v0cospa{"cascSelections.v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"cascSelections.dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcav0topv{"cascSelections.dcav0topv", .05, "min DCA V0 to PV (cm)"};
    Configurable<float> dcapiontopv{"cascSelections.dcapiontopv", .05, "min DCA Pion To PV (cm)"};
    Configurable<float> dcaprotontopv{"cascSelections.dcaprotontopv", .05, "min DCA Proton To PV (cm)"};
    Configurable<float> v0radius{"cascSelections.v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"cascSelections.v0radiusMax", 1E5, "maximum V0 radius (cm)"};

    // Standard 6 topological criteria on cascades
    Configurable<float> casccospa{"cascSelections.casccospa", 0.97, "min Cascade CosPA"};
    Configurable<float> dcacascdau{"cascSelections.dcacascdau", 1.0, "max DCA Cascade Daughters (cm)"};
    Configurable<float> dcaxybachbaryontopv{"cascSelections.dcaxybachbaryontopv", -1, "DCAxy Bachelor-Baryon to PV (cm)"};
    Configurable<float> bachbaryoncospa{"cascSelections.bachbaryoncospa", -1, "Bachelor-Baryon CosPA"};
    Configurable<float> dcabachtopv{"cascSelections.dcabachtopv", .05, "min DCA Bachelor To PV (cm)"};
    Configurable<float> cascradius{"cascSelections.cascradius", 0.5, "minimum Cascade radius (cm)"};
    Configurable<float> cascradiusMax{"cascSelections.cascradiusMax", 1E5, "maximum Cascade radius (cm)"};
    Configurable<float> cascProperLifeTime{"cascSelections.cascProperLifeTime", 3, "maximum lifetime (ctau)"};

    // invariant mass selection
    Configurable<float> v0MassWindow{"cascSelections.v0MassWindow", 0.008, "#Lambda mass (GeV/#it{c}^{2})"};
    Configurable<float> cascMassWindow{"cascSelections.cascMassWindow", 0.008, "#Lambda mass (GeV/#it{c}^{2})"};
    Configurable<float> compMassRejection{"cascSelections.compMassRejection", 0.008, "Competing mass rejection (GeV/#it{c}^{2})"};

    // Track quality
    Configurable<int> minTPCrows{"cascSelections.minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"cascSelections.minITSclusters", -1, "minimum ITS clusters"};
    Configurable<bool> skipTPConly{"cascSelections.skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
    Configurable<bool> requireBachITSonly{"cascSelections.requireBachITSonly", false, "require that bachelor track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requirePosITSonly{"cascSelections.requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireNegITSonly{"cascSelections.requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};

    // PID (TPC)
    Configurable<float> tpcPidNsigmaCut{"cascSelections.tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
  } cascSelections;

  struct : ConfigurableGroup {
    std::string prefix = "rctConfigurations"; // JSON group name
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "", "Which detector condition requirements? (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "Include ZDC flags in the bit selection (for Pb-Pb only)"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } rctConfigurations;

  RCTFlagsChecker rctFlagsChecker{rctConfigurations.cfgRCTLabel.value};

  // CCDB options
  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdbConfigurations.ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"ccdbConfigurations.grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"ccdbConfigurations.grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"ccdbConfigurations.lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"ccdbConfigurations.geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> mVtxPath{"ccdbConfigurations.mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } ccdbConfigurations;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  ctpRateFetcher rateFetcher;
  int mRunNumber;
  std::map<std::string, std::string> metadata;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f}, "Centrality"};
  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 1000.0f, 3000.0f, 10000.0f, 30000.0f}, "Occupancy"};

  // topological variable QA axes
  ConfigurableAxis axisMass{"axisV0Mass", {25, -0.05f, 0.05f}, "Invariant mass (GeV/#it{c}^{2})"};
  ConfigurableAxis axisPhi{"axisPhi", {36, 0.0f, constants::math::TwoPI}, "#varphi (rad)"};
  ConfigurableAxis axisEta{"axisEta", {10, -1.0f, 1.0f}, "Pseudo-rapidity #eta"};
  ConfigurableAxis axisRadius{"axisRadius", {10, 0.0f, 250.0f}, "Decay radius (cm)"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 4.0f, 5.0f, 7.0f, 9.0f, 11.0f, 15.0f, 30.0f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisAlphaV0{"axisAlphaV0", {1, -1.0f, 1.f}, "V0 #alpha Armenteros"};
  ConfigurableAxis axisPtArmV0{"axisPtArmV0", {1, 0.0f, 10.f}, "V0 #it{p}_{T} Armenteros"};

  // PDG database
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Sqrt scaling function
  // Author: Marian Ivanov
  int codeSqrtScaling(float val, float sigma0, float sigma1, int clampMin, int clampMax)
  {
    float code_f = std::asinh((sigma1 * val) / sigma0) / sigma0;
    return std::clamp(static_cast<int>(std::round(code_f)), clampMin, clampMax);
  }

  // Function to decode the sqrt scaling
  // Author: Marian Ivanov
  float decodeSqrtScaling(int code, float sigma0, float sigma1)
  {
    float code_f = static_cast<float>(code);
    return (sigma0 / sigma1) * std::sinh(sigma0 * code_f);
  }

  void init(InitContext const&)
  {
    if (analyseK0Short + analyseLambda + analyseAntiLambda + analyseXi + analyseOmega > 1) {
      LOGF(fatal, "Cannot enable several particles at the same time. Please choose one.");
    }

    // Initialise the RCTFlagsChecker
    rctFlagsChecker.init(rctConfigurations.cfgRCTLabel.value, rctConfigurations.cfgCheckZDC, rctConfigurations.cfgTreatLimitedAcceptanceAsBad);

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
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(21, "RCT flags");

    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, +100.0f}});
    histos.add("hEventOccupancy", "hEventOccupancy", kTH1F, {axisOccupancy});

    histos.add("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc", "h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc", kTHnSparseF, {axisMass, axisPt, axisPhi, axisEta, axisPtArmV0, axisAlphaV0, axisRadius, axisCentrality, axisOccupancy});
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(0)->SetName("Mass");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(1)->SetName("Pt");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(2)->SetName("Phi");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(3)->SetName("Eta");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(4)->SetName("V0PtArm");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(5)->SetName("V0Alpha");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(6)->SetName("Radius");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(7)->SetName("Centrality");
    histos.get<THnSparse>(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"))->GetAxis(8)->SetName("Occupancy");

    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    // inspect histogram sizes, please
    histos.print();
  }

  template <typename TCollision> // TCollision should be of the type: soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator or so
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }

    mRunNumber = collision.runNumber();
    if (cfgSkimmedProcessing) {
      ccdb->setURL(ccdbConfigurations.ccdburl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);

      zorro.initCCDB(ccdb.service, collision.runNumber(), collision.timestamp(), cfgSkimmedTrigger.value);
      zorro.populateHistRegistry(histos, collision.runNumber());
    }
  }

  template <typename TCollision>
  bool isEventAccepted(TCollision collision, bool fillHists)
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
    double interactionRate = (eventSelections.minIR >= 0 || eventSelections.maxIR >= 0) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource) * 1.e-3 : -1;
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

    if (!rctConfigurations.cfgRCTLabel.value.empty() && !rctFlagsChecker(collision)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 20 /* Pass CBT condition */);

    return true;
  }

  template <typename TCollision>
  void fillEventHistograms(TCollision collision, float& centrality, float& occupancy)
  {
    if (isPP) { //
      centrality = collision.centFT0M();
    } else {
      centrality = collision.centFT0C();
      occupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (encodingOpts.useSqrtEncodingForOccupancy)
        occupancy = std::sqrt(occupancy);
    }
    histos.fill(HIST("hEventCentrality"), centrality);
    histos.fill(HIST("hEventOccupancy"), occupancy);

    return;
  }

  template <typename TV0, typename TCollision>
  bool isV0Selected(TV0 v0, TCollision collision, float rapidity)
  // precalculate this information so that a check is one mask operation, not many
  {
    //
    // Base topological variables
    //

    // v0 radius min/max selections
    if (v0.v0radius() < v0Selections.v0radius)
      return false;
    if (v0.v0radius() > v0Selections.v0radiusMax)
      return false;
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (analyseK0Short) {
      if (std::fabs(v0.dcapostopv()) < v0Selections.dcapiontopv)
        return false;
      if (std::fabs(v0.dcanegtopv()) < v0Selections.dcapiontopv)
        return false;
    }
    if (analyseLambda) {
      if (std::fabs(v0.dcapostopv()) < v0Selections.dcaprotontopv)
        return false;
      if (std::fabs(v0.dcanegtopv()) < v0Selections.dcapiontopv)
        return false;
    }
    if (analyseAntiLambda) {
      if (std::fabs(v0.dcapostopv()) < v0Selections.dcapiontopv)
        return false;
      if (std::fabs(v0.dcanegtopv()) < v0Selections.dcaprotontopv)
        return false;
    }
    // V0 cosine of pointing angle
    if (v0.v0cosPA() < v0Selections.v0cospa)
      return false;
    // DCA between v0 daughters
    if (v0.dcaV0daughters() > v0Selections.dcav0dau)
      return false;
    // DCA V0 to prim vtx
    if (v0.dcav0topv() < v0Selections.dcav0topv)
      return false;

    //
    // rapidity
    //
    if (std::fabs(rapidity) > v0Selections.rapidityCut)
      return false;

    //
    // competing mass rejection
    //
    if ((analyseLambda || analyseAntiLambda) && std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0Selections.compMassRejection)
      return false;
    if (analyseK0Short && std::fabs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0Selections.compMassRejection)
      return false;

    auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

    //
    // ITS quality flags
    //
    if (posTrackExtra.itsNCls() < v0Selections.minITSclusters)
      return false;
    if (negTrackExtra.itsNCls() < v0Selections.minITSclusters)
      return false;

    //
    // TPC quality flags
    //
    if (posTrackExtra.tpcCrossedRows() < v0Selections.minTPCrows)
      return false;
    if (negTrackExtra.tpcCrossedRows() < v0Selections.minTPCrows)
      return false;

    //
    // TPC PID
    //
    if (analyseK0Short) {
      if (std::fabs(posTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut)
        return false;
    }
    if (analyseLambda) {
      if (std::fabs(posTrackExtra.tpcNSigmaPr()) > v0Selections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut)
        return false;
    }
    if (analyseAntiLambda) {
      if (std::fabs(posTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPr()) > v0Selections.tpcPidNsigmaCut)
        return false;
    }

    //
    // ITS only tag
    if (v0Selections.requirePosITSonly && posTrackExtra.tpcCrossedRows() > 0)
      return false;
    if (v0Selections.requireNegITSonly && negTrackExtra.tpcCrossedRows() > 0)
      return false;

    //
    // TPC only tag
    if (v0Selections.requireTPConly && posTrackExtra.detectorMap() != o2::aod::track::TPC)
      return false;
    if (v0Selections.requireTPConly && negTrackExtra.detectorMap() != o2::aod::track::TPC)
      return false;

    //
    // proper lifetime
    if ((analyseLambda || analyseAntiLambda) &&
        v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 > lifetimecut->get("lifetimecutLambda"))
      return false;
    if (analyseK0Short &&
        v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short > lifetimecut->get("lifetimecutK0S"))
      return false;

    //
    // armenteros
    if (v0Selections.armPodCut > 1e-4 && v0.qtarm() * v0Selections.armPodCut < std::fabs(v0.alpha()))
      return false;

    return true;
  }

  template <typename TCascade, typename TCollision>
  bool isCascadeSelected(TCascade casc, TCollision collision, float rapidity)
  // precalculate this information so that a check is one mask operation, not many
  {
    //
    // Base topological variables
    //

    // v0 radius min/max selections
    if (casc.v0radius() < cascSelections.v0radius)
      return false;
    if (casc.v0radius() > cascSelections.v0radiusMax)
      return false;
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (casc.sign() < 0) { // Xi- or Omega- --> positive/negative daughter = proton/pion
      if (std::fabs(casc.dcapostopv()) < cascSelections.dcaprotontopv)
        return false;
      if (std::fabs(casc.dcanegtopv()) < cascSelections.dcapiontopv)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (std::fabs(casc.dcapostopv()) < cascSelections.dcapiontopv)
        return false;
      if (std::fabs(casc.dcanegtopv()) < cascSelections.dcaprotontopv)
        return false;
    }
    // V0 cosine of pointing angle
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.v0cospa)
      return false;
    // DCA between v0 daughters
    if (casc.dcaV0daughters() > cascSelections.dcav0dau)
      return false;
    // DCA V0 to prim vtx
    if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.dcav0topv)
      return false;

    // casc radius min/max selections
    if (casc.cascradius() < cascSelections.cascradius)
      return false;
    if (casc.cascradius() > cascSelections.cascradiusMax)
      return false;
    // DCA bachelor selection
    if (std::fabs(casc.dcabachtopv()) < cascSelections.dcabachtopv)
      return false;
    // Bachelor-baryon cosPA selection
    if (casc.bachBaryonCosPA() < cascSelections.bachbaryoncospa)
      return false;
    // DCA bachelor-baryon selection
    if (std::fabs(casc.bachBaryonDCAxyToPV()) < cascSelections.dcaxybachbaryontopv)
      return false;
    // casc cosine of pointing angle
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.casccospa)
      return false;
    // DCA between casc daughters
    if (casc.dcacascdaughters() > cascSelections.dcacascdau)
      return false;

    //
    // rapidity
    //
    if (std::fabs(rapidity) > cascSelections.rapidityCut)
      return false;

    //
    // competing mass rejection
    //
    if (analyseXi && std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < cascSelections.compMassRejection)
      return false;
    if (analyseOmega && std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) < cascSelections.compMassRejection)
      return false;

    auto bachTrackExtra = casc.template bachTrackExtra_as<DauTracks>();
    auto posTrackExtra = casc.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = casc.template negTrackExtra_as<DauTracks>();

    //
    // ITS quality flags
    //
    if (bachTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;
    if (posTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;
    if (negTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;

    //
    // TPC quality flags
    //
    if (bachTrackExtra.tpcCrossedRows() < cascSelections.minTPCrows)
      return false;
    if (posTrackExtra.tpcCrossedRows() < cascSelections.minTPCrows)
      return false;
    if (negTrackExtra.tpcCrossedRows() < cascSelections.minTPCrows)
      return false;

    //
    // TPC PID
    //
    if (analyseXi && std::fabs(bachTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
      return false;
    if (analyseOmega && std::fabs(bachTrackExtra.tpcNSigmaKa()) > cascSelections.tpcPidNsigmaCut)
      return false;
    if (casc.sign() < 0) { // Xi- or Omega- --> positive/negative daughter = proton/pion
      if (std::fabs(posTrackExtra.tpcNSigmaPr()) > cascSelections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (std::fabs(posTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPr()) > cascSelections.tpcPidNsigmaCut)
        return false;
    }

    //
    // proper lifetime
    float distOverTotMom = std::sqrt(std::pow(casc.x() - collision.posX(), 2) + std::pow(casc.y() - collision.posY(), 2) + std::pow(casc.z() - collision.posZ(), 2)) / (casc.p() + 1E-10);
    if (analyseXi && distOverTotMom * o2::constants::physics::MassXiMinus / ctauXiPDG > cascSelections.cascProperLifeTime)
      return false;
    if (analyseOmega && distOverTotMom * o2::constants::physics::MassOmegaMinus / ctauOmegaPDG > cascSelections.cascProperLifeTime)
      return false;

    return true;
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  void process(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator const& collision, V0Candidates const& fullV0s, CascadeCandidates const& fullCascades, DauTracks const&)
  {
    // Fire up CCDB
    if (cfgSkimmedProcessing) {
      initCCDB(collision);
    }

    if (!isEventAccepted(collision, true)) {
      return;
    }

    if (cfgSkimmedProcessing) {
      zorro.isSelected(collision.globalBC()); /// Just let Zorro do the accounting
    }

    float centrality = -1;
    float occupancy = -1;
    fillEventHistograms(collision, centrality, occupancy);

    // __________________________________________
    // perform main analysis
    //
    if (analyseK0Short || analyseLambda || analyseAntiLambda) { // Look at V0s
      for (const auto& v0 : fullV0s) {
        if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut || std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        if (v0.v0Type() != v0Selections.v0TypeSelection && v0Selections.v0TypeSelection > -1)
          continue; // skip V0s that are not standard

        float pT = encodingOpts.useSqrtScalingForEncodingPt ? codeSqrtScaling(v0.pt(), encodingOpts.sqrtScalingParameters->get("sigma0"), encodingOpts.sqrtScalingParameters->get("sigma1"), encodingOpts.sqrtScalingParameters->get("clampMin"), encodingOpts.sqrtScalingParameters->get("clampMax")) : v0.pt();
        float decayRadius = encodingOpts.useSqrtEncodingForRadius ? std::sqrt(v0.v0radius()) : v0.v0radius();

        if (analyseK0Short && isV0Selected(v0, collision, v0.yK0Short())) {
          histos.fill(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"), v0.mK0Short() - o2::constants::physics::MassK0Short, pT, v0.phi(), v0.eta(), v0.qtarm(), v0.alpha(), decayRadius, centrality, occupancy);
        }
        if (analyseLambda && isV0Selected(v0, collision, v0.yLambda())) {
          histos.fill(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"), v0.mLambda() - o2::constants::physics::MassLambda0, pT, v0.phi(), v0.eta(), v0.qtarm(), v0.alpha(), decayRadius, centrality, occupancy);
        }
        if (analyseAntiLambda && isV0Selected(v0, collision, v0.yLambda())) {
          histos.fill(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"), v0.mAntiLambda() - o2::constants::physics::MassLambda0, pT, v0.phi(), v0.eta(), v0.qtarm(), v0.alpha(), decayRadius, centrality, occupancy);
        }
      } // end v0 loop
    }

    if (analyseXi || analyseOmega) { // Look at Cascades
      for (const auto& cascade : fullCascades) {
        if (std::abs(cascade.negativeeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.positiveeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.bacheloreta()) > cascSelections.daughterEtaCut)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        float pT = encodingOpts.useSqrtScalingForEncodingPt ? codeSqrtScaling(cascade.pt(), encodingOpts.sqrtScalingParameters->get("sigma0"), encodingOpts.sqrtScalingParameters->get("sigma1"), encodingOpts.sqrtScalingParameters->get("clampMin"), encodingOpts.sqrtScalingParameters->get("clampMax")) : cascade.pt();
        float decayRadius = encodingOpts.useSqrtEncodingForRadius ? std::sqrt(cascade.cascradius()) : cascade.cascradius();

        if (analyseXi && isCascadeSelected(cascade, collision, cascade.yXi())) {
          histos.fill(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"), cascade.m(1) - o2::constants::physics::MassXiMinus, pT, cascade.phi(), cascade.eta(), 0., 0., decayRadius, centrality, occupancy);
        }
        if (analyseOmega && isCascadeSelected(cascade, collision, cascade.yOmega())) {
          histos.fill(HIST("h9dMassPtPhiEtaPtArmV0AlphaV0RadiusCentOcc"), cascade.m(2) - o2::constants::physics::MassOmegaMinus, pT, cascade.phi(), cascade.eta(), 0., 0., decayRadius, centrality, occupancy);
        }
      } // end cascade loop
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessderivedbinnedinfo>(cfgc)};
}
