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
// V0 analysis task
// ================
//
// This code loops over a V0Cores table and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
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
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/Core/SGSelector.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

// constants
const float ctauXiPDG = 4.91;     // from PDG
const float ctauOmegaPDG = 2.461; // from PDG

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using dauMCTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackMCIds, aod::DauTrackTPCPIDs>;
using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0K0ShortMLScores>;
// using v0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0MCCores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0MCCollRefs>;
using v0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0K0ShortMLScores>;

using cascadeCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascXiMLScores, aod::CascOmMLScores>;
using cascadeMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascXiMLScores, aod::CascOmMLScores, aod::CascCoreMCLabels>;

// simple checkers, but ensure 64 bit integers
#define bitset(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

struct quarkoniaToHyperons {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> isPP{"isPP", true, "If running on pp collision, switch it on true"};

  Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
  Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
  Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", true, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
  Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};

  Configurable<bool> buildLaLaBarPairs{"buildLaLaBarPairs", false, "Build Lambda antiLambda from charmonia decay"};
  Configurable<bool> buildXiXiBarPairs{"buildXiXiBarPairs", false, "Build Xi antiXi from charmonia decay"};
  Configurable<bool> buildOmOmBarPairs{"buildOmOmBarPairs", false, "Build Omega antiOmega from charmonia decay"};

  // fast check on occupancy
  Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
  Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};

  // rapidity cut on the hyperon-antiHyperon pair
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity cut on the hyp-antiHyp pair"};

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
    Configurable<bool> skipTPConly{"v0Selections.skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
    Configurable<bool> requirePosITSonly{"v0Selections.requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireNegITSonly{"v0Selections.requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};

    // PID (TPC/TOF)
    Configurable<float> TpcPidNsigmaCut{"v0Selections.TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
    Configurable<float> TofPidNsigmaCutLaPr{"v0Selections.TofPidNsigmaCutLaPr", 1e+6, "TofPidNsigmaCutLaPr"};
    Configurable<float> TofPidNsigmaCutLaPi{"v0Selections.TofPidNsigmaCutLaPi", 1e+6, "TofPidNsigmaCutLaPi"};
    Configurable<float> TofPidNsigmaCutK0Pi{"v0Selections.TofPidNsigmaCutK0Pi", 1e+6, "TofPidNsigmaCutK0Pi"};

    // PID (TOF)
    Configurable<float> maxDeltaTimeProton{"v0Selections.maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimePion{"v0Selections.maxDeltaTimePion", 1e+9, "check maximum allowed time"};
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

    // PID (TPC/TOF)
    Configurable<float> TpcPidNsigmaCut{"cascSelections.TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};
    Configurable<float> TofPidNsigmaCutLaPr{"cascSelections.TofPidNsigmaCutLaPr", 1e+6, "TofPidNsigmaCutLaPr"};
    Configurable<float> TofPidNsigmaCutLaPi{"cascSelections.TofPidNsigmaCutLaPi", 1e+6, "TofPidNsigmaCutLaPi"};
    Configurable<float> TofPidNsigmaCutXiPi{"cascSelections.TofPidNsigmaCutXiPi", 1e+6, "TofPidNsigmaCutXiPi"};
    Configurable<float> TofPidNsigmaCutOmKa{"cascSelections.TofPidNsigmaCutOmKa", 1e+6, "TofPidNsigmaCutOmKa"};

    // PID (TOF)
    Configurable<float> maxDeltaTimeProton{"cascSelections.maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimePion{"cascSelections.maxDeltaTimePion", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimeKaon{"cascSelections.maxDeltaTimeKaon", 1e+9, "check maximum allowed time"};
  } cascSelections;

  Configurable<bool> qaCentrality{"qaCentrality", false, "qa centrality flag: check base raw values"};

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};

  // UPC selections
  SGSelector sgSelector;
  struct : ConfigurableGroup {
    Configurable<float> FV0cut{"upcCuts.FV0cut", 100., "FV0A threshold"};
    Configurable<float> FT0Acut{"upcCuts.FT0Acut", 200., "FT0A threshold"};
    Configurable<float> FT0Ccut{"upcCuts.FT0Ccut", 100., "FT0C threshold"};
    Configurable<float> ZDCcut{"upcCuts.ZDCcut", 10., "ZDC threshold"};
    // Configurable<float> gapSel{"upcCuts.gapSel", 2, "Gap selection"};
  } upcCuts;

  // Machine learning evaluation for pre-selection and corresponding information generation
  o2::ml::OnnxModel mlCustomModelK0Short;
  o2::ml::OnnxModel mlCustomModelLambda;
  o2::ml::OnnxModel mlCustomModelAntiLambda;
  o2::ml::OnnxModel mlCustomModelGamma;

  struct : ConfigurableGroup {
    // ML classifiers: master flags to control whether we should use custom ML classifiers or the scores in the derived data
    Configurable<bool> useK0ShortScores{"mlConfigurations.useK0ShortScores", false, "use ML scores to select K0Short"};
    Configurable<bool> useLambdaScores{"mlConfigurations.useLambdaScores", false, "use ML scores to select Lambda"};
    Configurable<bool> useAntiLambdaScores{"mlConfigurations.useAntiLambdaScores", false, "use ML scores to select AntiLambda"};

    Configurable<bool> calculateK0ShortScores{"mlConfigurations.calculateK0ShortScores", false, "calculate K0Short ML scores"};
    Configurable<bool> calculateLambdaScores{"mlConfigurations.calculateLambdaScores", false, "calculate Lambda ML scores"};
    Configurable<bool> calculateAntiLambdaScores{"mlConfigurations.calculateAntiLambdaScores", false, "calculate AntiLambda ML scores"};

    // ML input for ML calculation
    Configurable<std::string> customModelPathCCDB{"mlConfigurations.customModelPathCCDB", "", "Custom ML Model path in CCDB"};
    Configurable<int64_t> timestampCCDB{"mlConfigurations.timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadCustomModelsFromCCDB{"mlConfigurations.loadCustomModelsFromCCDB", false, "Flag to enable or disable the loading of custom models from CCDB"};
    Configurable<bool> enableOptimizations{"mlConfigurations.enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

    // Local paths for test purposes
    Configurable<std::string> localModelPathLambda{"mlConfigurations.localModelPathLambda", "Lambda_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
    Configurable<std::string> localModelPathAntiLambda{"mlConfigurations.localModelPathAntiLambda", "AntiLambda_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
    Configurable<std::string> localModelPathK0Short{"mlConfigurations.localModelPathK0Short", "KZeroShort_BDTModel.onnx", "(std::string) Path to the local .onnx file."};

    // Thresholds for choosing to populate V0Cores tables with pre-selections
    Configurable<float> thresholdLambda{"mlConfigurations.thresholdLambda", -1.0f, "Threshold to keep Lambda candidates"};
    Configurable<float> thresholdAntiLambda{"mlConfigurations.thresholdAntiLambda", -1.0f, "Threshold to keep AntiLambda candidates"};
    Configurable<float> thresholdK0Short{"mlConfigurations.thresholdK0Short", -1.0f, "Threshold to keep K0Short candidates"};
  } mlConfigurations;

  // CCDB options
  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdbConfigurations.ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"ccdbConfigurations.grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"ccdbConfigurations.grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"ccdbConfigurations.lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"ccdbConfigurations.geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> mVtxPath{"ccdbConfigurations.mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } ccdbConfigurations;

  o2::ccdb::CcdbApi ccdbApi;
  int mRunNumber;
  std::map<std::string, std::string> metadata;

  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.2f, 1.4f, 1.6f, 1.8f, 2.0f, 2.4f, 2.8f, 3.2f, 3.6f, 4.0f, 4.8f, 5.6f, 6.5f, 7.5f, 9.0f, 11.0f, 13.0f, 15.0f, 19.0f, 23.0f, 30.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisQuarkoniumMass{"axisQuarkoniumMass", {500, 2.600f, 4.000f}, "M (hyp. #bar{hyp.} ) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};
  ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles"};

  ConfigurableAxis axisRawCentrality{"axisRawCentrality", {VARIABLE_WIDTH, 0.000f, 52.320f, 75.400f, 95.719f, 115.364f, 135.211f, 155.791f, 177.504f, 200.686f, 225.641f, 252.645f, 281.906f, 313.850f, 348.302f, 385.732f, 426.307f, 470.146f, 517.555f, 568.899f, 624.177f, 684.021f, 748.734f, 818.078f, 892.577f, 973.087f, 1058.789f, 1150.915f, 1249.319f, 1354.279f, 1465.979f, 1584.790f, 1710.778f, 1844.863f, 1985.746f, 2134.643f, 2291.610f, 2456.943f, 2630.653f, 2813.959f, 3006.631f, 3207.229f, 3417.641f, 3637.318f, 3865.785f, 4104.997f, 4354.938f, 4615.786f, 4885.335f, 5166.555f, 5458.021f, 5762.584f, 6077.881f, 6406.834f, 6746.435f, 7097.958f, 7462.579f, 7839.165f, 8231.629f, 8635.640f, 9052.000f, 9484.268f, 9929.111f, 10389.350f, 10862.059f, 11352.185f, 11856.823f, 12380.371f, 12920.401f, 13476.971f, 14053.087f, 14646.190f, 15258.426f, 15890.617f, 16544.433f, 17218.024f, 17913.465f, 18631.374f, 19374.983f, 20136.700f, 20927.783f, 21746.796f, 22590.880f, 23465.734f, 24372.274f, 25314.351f, 26290.488f, 27300.899f, 28347.512f, 29436.133f, 30567.840f, 31746.818f, 32982.664f, 34276.329f, 35624.859f, 37042.588f, 38546.609f, 40139.742f, 41837.980f, 43679.429f, 45892.130f, 400000.000f}, "raw centrality signal"}; // for QA

  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};

  // topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {20, 0.0f, 1.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {20, 0.0f, 2.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAV0ToPV{"axisDCAV0ToPV", {20, 0.0f, 2.0f}, "DCA (cm)"};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {20, 0.0f, 2.0f}, "pointing angle (rad)"};
  ConfigurableAxis axisRadius{"axisRadius", {20, 0.0f, 60.0f}, "Decay radius (cm)"};
  ConfigurableAxis axisProperLifeTime{"axisV0ProperLifeTime", {100, 0.0f, 50.0f}, "ProperLifeTime 2D radius (cm)"};
  ConfigurableAxis axisMassWindow{"axisMassWindow", {40, -0.020f, 0.020f}, "Inv. mass - PDG mass (GeV/#it{c}^{2})"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {500, 0.400f, 0.600f}, "K0Short mass (GeV/#it{c}^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {500, 1.098f, 1.198f}, "Lambda mass (GeV/#it{c}^{2})"};
  ConfigurableAxis axisXiMass{"axisXiMass", {500, 1.318f, 1.370f}, "Xi mass (GeV/#it{c}^{2})"};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {500, 1.670f, 1.675f}, "Omega mass (GeV/#it{c}^{2})"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};

  // Track quality axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisITSclus{"axisITSclus", {7, 0.0f, 7.0f}, "N ITS Clusters"};

  // UPC axes
  ConfigurableAxis axisSelGap{"axisSelGap", {4, -1.5, 2.5}, "Gap side"};

  // PDG database
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // For manual sliceBy
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

  enum selection : uint64_t { selCosPA = 0,
                              selRadius,
                              selRadiusMax,
                              selDCANegToPV,
                              selDCAPosToPV,
                              selDCAV0ToPV,
                              selDCAV0Dau,
                              selK0ShortRapidity,
                              selLambdaRapidity,
                              selK0ShortMassWindow,
                              selLambdaMassWindow,
                              selAntiLambdaMassWindow,
                              selK0ShortMassRejection,
                              selLambdaMassRejection,
                              selTPCPIDPositivePion,
                              selTPCPIDNegativePion,
                              selTPCPIDPositiveProton,
                              selTPCPIDNegativeProton,
                              selTOFDeltaTPositiveProtonLambda,
                              selTOFDeltaTPositivePionLambda,
                              selTOFDeltaTPositivePionK0Short,
                              selTOFDeltaTNegativeProtonLambda,
                              selTOFDeltaTNegativePionLambda,
                              selTOFDeltaTNegativePionK0Short,
                              selTOFNSigmaPositiveProtonLambda, // Nsigma
                              selTOFNSigmaPositivePionLambda,   // Nsigma
                              selTOFNSigmaPositivePionK0Short,  // Nsigma
                              selTOFNSigmaNegativeProtonLambda, // Nsigma
                              selTOFNSigmaNegativePionLambda,   // Nsigma
                              selTOFNSigmaNegativePionK0Short,  // Nsigma
                              selK0ShortCTau,
                              selLambdaCTau,
                              selK0ShortArmenteros,
                              selPosGoodTPCTrack, // at least min # TPC rows
                              selNegGoodTPCTrack, // at least min # TPC rows
                              selPosGoodITSTrack, // at least min # ITS clusters
                              selNegGoodITSTrack, // at least min # ITS clusters
                              selPosItsOnly,
                              selNegItsOnly,
                              selPosNotTPCOnly,
                              selNegNotTPCOnly,
                              selConsiderK0Short,    // for mc tagging
                              selConsiderLambda,     // for mc tagging
                              selConsiderAntiLambda, // for mc tagging
                              selPhysPrimK0Short,    // for mc tagging
                              selPhysPrimLambda,     // for mc tagging
                              selPhysPrimAntiLambda, // for mc tagging
  };

  uint64_t maskTopological;
  uint64_t maskTopoNoV0Radius;
  uint64_t maskTopoNoDCANegToPV;
  uint64_t maskTopoNoDCAPosToPV;
  uint64_t maskTopoNoCosPA;
  uint64_t maskTopoNoDCAV0Dau;
  uint64_t maskTopoNoDCAV0ToPV;
  uint64_t maskTrackProperties;

  uint64_t maskK0ShortSpecific;
  uint64_t maskLambdaSpecific;
  uint64_t maskAntiLambdaSpecific;

  uint64_t maskSelectionK0Short;
  uint64_t maskSelectionLambda;
  uint64_t maskSelectionAntiLambda;

  uint64_t secondaryMaskSelectionLambda;
  uint64_t secondaryMaskSelectionAntiLambda;

  void init(InitContext const&)
  {
    // initialise bit masks
    maskTopological = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0ToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoV0Radius = (uint64_t(1) << selCosPA) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0ToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCANegToPV = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0ToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCAPosToPV = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAV0ToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoCosPA = (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0ToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCAV0Dau = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0ToPV) | (uint64_t(1) << selRadiusMax);
    maskTopoNoDCAV0ToPV = (uint64_t(1) << selCosPA) | (uint64_t(1) << selRadius) | (uint64_t(1) << selDCANegToPV) | (uint64_t(1) << selDCAPosToPV) | (uint64_t(1) << selDCAV0Dau) | (uint64_t(1) << selRadiusMax);

    maskK0ShortSpecific = (uint64_t(1) << selK0ShortRapidity) | (uint64_t(1) << selK0ShortCTau) | (uint64_t(1) << selK0ShortArmenteros) | (uint64_t(1) << selConsiderK0Short) | (uint64_t(1) << selK0ShortMassWindow) | (uint64_t(1) << selLambdaMassRejection);
    maskLambdaSpecific = (uint64_t(1) << selLambdaRapidity) | (uint64_t(1) << selLambdaCTau) | (uint64_t(1) << selConsiderLambda) | (uint64_t(1) << selLambdaMassWindow) | (uint64_t(1) << selK0ShortMassRejection);
    maskAntiLambdaSpecific = (uint64_t(1) << selLambdaRapidity) | (uint64_t(1) << selLambdaCTau) | (uint64_t(1) << selConsiderAntiLambda) | (uint64_t(1) << selAntiLambdaMassWindow) | (uint64_t(1) << selK0ShortMassRejection);

    // ask for specific TPC/TOF PID selections
    maskTrackProperties = 0;
    if (v0Selections.requirePosITSonly) {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selPosItsOnly) | (uint64_t(1) << selPosGoodITSTrack);
    } else {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selPosGoodTPCTrack) | (uint64_t(1) << selPosGoodITSTrack);
      // TPC signal is available: ask for positive track PID
      if (v0Selections.TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTPCPIDPositivePion);
        maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTPCPIDPositiveProton);
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTPCPIDPositivePion);
      }
      // TOF PID
      if (v0Selections.TofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTOFNSigmaPositivePionK0Short) | (uint64_t(1) << selTOFDeltaTPositivePionK0Short);
      if (v0Selections.TofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
        maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTOFNSigmaPositiveProtonLambda) | (uint64_t(1) << selTOFDeltaTPositiveProtonLambda);
      if (v0Selections.TofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTOFNSigmaPositivePionLambda) | (uint64_t(1) << selTOFDeltaTPositivePionLambda);
    }
    if (v0Selections.requireNegITSonly) {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selNegItsOnly) | (uint64_t(1) << selNegGoodITSTrack);
    } else {
      maskTrackProperties = maskTrackProperties | (uint64_t(1) << selNegGoodTPCTrack) | (uint64_t(1) << selNegGoodITSTrack);
      // TPC signal is available: ask for negative track PID
      if (v0Selections.TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTPCPIDNegativePion);
        maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTPCPIDNegativePion);
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTPCPIDNegativeProton);
      }
      // TOF PID
      if (v0Selections.TofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selTOFNSigmaNegativePionK0Short) | (uint64_t(1) << selTOFDeltaTNegativePionK0Short);
      if (v0Selections.TofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
        maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selTOFNSigmaNegativePionLambda) | (uint64_t(1) << selTOFDeltaTNegativePionLambda);
      if (v0Selections.TofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selTOFNSigmaNegativeProtonLambda) | (uint64_t(1) << selTOFDeltaTNegativeProtonLambda);
    }

    if (v0Selections.skipTPConly) {
      maskK0ShortSpecific = maskK0ShortSpecific | (uint64_t(1) << selPosNotTPCOnly) | (uint64_t(1) << selNegNotTPCOnly);
      maskLambdaSpecific = maskLambdaSpecific | (uint64_t(1) << selPosNotTPCOnly) | (uint64_t(1) << selNegNotTPCOnly);
      maskAntiLambdaSpecific = maskAntiLambdaSpecific | (uint64_t(1) << selPosNotTPCOnly) | (uint64_t(1) << selNegNotTPCOnly);
    }

    // Primary particle selection, central to analysis
    maskSelectionK0Short = maskTopological | maskTrackProperties | maskK0ShortSpecific;
    maskSelectionLambda = maskTopological | maskTrackProperties | maskLambdaSpecific;
    maskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific;

    // No primary requirement for feeddown matrix
    secondaryMaskSelectionLambda = maskTopological | maskTrackProperties | maskLambdaSpecific;
    secondaryMaskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific;

    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{20, -0.5f, +19.5f}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(11, "kNoCollInTimeRangeStd");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(13, "Below min occup.");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(14, "Above max occup.");

    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, +100.0f}});
    histos.add("hCentralityVsNch", "hCentralityVsNch", kTH2F, {axisCentrality, axisNch});

    histos.add("hEventOccupancy", "hEventOccupancy", kTH1F, {axisOccupancy});
    histos.add("hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH2F, {axisCentrality, axisOccupancy});

    if (!isPP) {
      histos.add("hGapSide", "Gap side; Entries", kTH1F, {{5, -0.5, 4.5}});
      histos.add("hSelGapSide", "Selected gap side; Entries", kTH1F, {axisSelGap});
      histos.add("hEventCentralityVsSelGapSide", ";Centrality (%); Selected gap side", kTH2F, {{100, 0.0f, +100.0f}, axisSelGap});
    }

    // for QA and test purposes
    auto hRawCentrality = histos.add<TH1>("hRawCentrality", "hRawCentrality", kTH1F, {axisRawCentrality});

    for (int ii = 1; ii < 101; ii++) {
      float value = 100.5f - static_cast<float>(ii);
      hRawCentrality->SetBinContent(ii, value);
    }

    // histograms versus mass
    if (buildLaLaBarPairs) {
      histos.add("LaLaBar/h3dMassLaLabar", "h3dMassLaLabar", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      if (!isPP) {
        // Non-UPC info
        histos.add("LaLaBar/h3dMassLaLabarHadronic", "h3dMassLaLabarHadronic", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        // UPC info
        histos.add("LaLaBar/h3dMassLaLabarSGA", "h3dMassLaLabarSGA", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dMassLaLabarSGC", "h3dMassLaLabarSGC", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dMassLaLabarDG", "h3dMassLaLabarDG", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      }
      histos.add("LaLaBar/h2dNbrOfK0ShortVsCentrality", "h2dNbrOfK0ShortVsCentrality", kTH2F, {axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("LaLaBar/h2dNbrOfLambdaVsCentrality", "h2dNbrOfLambdaVsCentrality", kTH2F, {axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("LaLaBar/h2dNbrOfAntiLambdaVsCentrality", "h2dNbrOfAntiLambdaVsCentrality", kTH2F, {axisCentrality, {10, -0.5f, 9.5f}});
      // QA plot
      // Candidates after Lambda selections
      histos.add("LaLaBar/Lambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("LaLaBar/Lambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("LaLaBar/Lambda/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axisDCAdau});
      histos.add("LaLaBar/Lambda/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
      histos.add("LaLaBar/Lambda/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axisPointingAngle});
      histos.add("LaLaBar/Lambda/hV0Radius", "hV0Radius", kTH1F, {axisRadius});
      histos.add("LaLaBar/Lambda/hV0DecayLength", "hDecayLength", kTH1F, {axisProperLifeTime});
      histos.add("LaLaBar/Lambda/hV0InvMassWindow", "hInvMassWindow", kTH1F, {axisMassWindow});
      histos.add("LaLaBar/Lambda/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axisLambdaMass, axisK0Mass});
      histos.add("LaLaBar/Lambda/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("LaLaBar/Lambda/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("LaLaBar/Lambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("LaLaBar/Lambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});

      // Candidates after AntiLambda selections
      histos.add("LaLaBar/AntiLambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("LaLaBar/AntiLambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("LaLaBar/AntiLambda/hDCAV0Daughters", "hDCADaughters", kTH1F, {axisDCAdau});
      histos.add("LaLaBar/AntiLambda/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
      histos.add("LaLaBar/AntiLambda/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axisPointingAngle});
      histos.add("LaLaBar/AntiLambda/hV0Radius", "hV0Radius", kTH1F, {axisRadius});
      histos.add("LaLaBar/AntiLambda/hV0DecayLength", "hDecayLength", kTH1F, {axisProperLifeTime});
      histos.add("LaLaBar/AntiLambda/hV0InvMassWindow", "hInvMassWindow", kTH1F, {axisMassWindow});
      histos.add("LaLaBar/AntiLambda/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axisLambdaMass, axisK0Mass});
      histos.add("LaLaBar/AntiLambda/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("LaLaBar/AntiLambda/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("LaLaBar/AntiLambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("LaLaBar/AntiLambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      if (doMCAssociation) {
        histos.add("LaLaBar/h3dInvMassTrueEtaC1S", "h3dInvMassTrueEtaC1S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueJPsi", "h3dInvMassTrueJPsi", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueChiC0", "h3dInvMassTrueChiC0", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueChiC1", "h3dInvMassTrueChiC1", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueHC", "h3dInvMassTrueHC", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueChiC2", "h3dInvMassTrueChiC2", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueEtaC2S", "h3dInvMassTrueEtaC2S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTruePsi2S", "h3dInvMassTruePsi2S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      }
    }
    if (buildXiXiBarPairs) {
      histos.add("XiXiBar/h3dMassXiXibar", "h3dMassXiXibar", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      if (!isPP) {
        // Non-UPC info
        histos.add("XiXiBar/h3dMassXiXibarHadronic", "h3dMassXiXibarHadronic", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        // UPC info
        histos.add("XiXiBar/h3dMassXiXibarSGA", "h3dMassXiXibarSGA", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dMassXiXibarSGC", "h3dMassXiXibarSGC", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dMassXiXibarDG", "h3dMassXiXibarDG", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      }
      histos.add("XiXiBar/h2dNbrOfXiVsCentrality", "h2dNbrOfXiVsCentrality", kTH2F, {axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("XiXiBar/h2dNbrOfAntiXiVsCentrality", "h2dNbrOfAntiXiVsCentrality", kTH2F, {axisCentrality, {10, -0.5f, 9.5f}});
      // QA plot
      // Candidates after Xi selections
      histos.add("XiXiBar/Xi/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("XiXiBar/Xi/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("XiXiBar/Xi/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("XiXiBar/Xi/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axisDCAdau});
      histos.add("XiXiBar/Xi/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axisDCAdau});
      histos.add("XiXiBar/Xi/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
      histos.add("XiXiBar/Xi/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axisPointingAngle});
      histos.add("XiXiBar/Xi/hV0Radius", "hV0Radius", kTH1F, {axisRadius});
      histos.add("XiXiBar/Xi/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axisPointingAngle});
      histos.add("XiXiBar/Xi/hCascRadius", "hCascRadius", kTH1F, {axisRadius});
      histos.add("XiXiBar/Xi/hCascDecayLength", "hCascDecayLength", kTH1F, {axisProperLifeTime});
      histos.add("XiXiBar/Xi/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axisMassWindow});
      histos.add("XiXiBar/Xi/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axisMassWindow});
      histos.add("XiXiBar/Xi/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axisXiMass, axisOmegaMass});
      histos.add("XiXiBar/Xi/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("XiXiBar/Xi/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("XiXiBar/Xi/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("XiXiBar/Xi/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("XiXiBar/Xi/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("XiXiBar/Xi/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      // Candidates after AntiXi selections
      histos.add("XiXiBar/AntiXi/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("XiXiBar/AntiXi/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("XiXiBar/AntiXi/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("XiXiBar/AntiXi/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axisDCAdau});
      histos.add("XiXiBar/AntiXi/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axisDCAdau});
      histos.add("XiXiBar/AntiXi/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
      histos.add("XiXiBar/AntiXi/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axisPointingAngle});
      histos.add("XiXiBar/AntiXi/hV0Radius", "hV0Radius", kTH1F, {axisRadius});
      histos.add("XiXiBar/AntiXi/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axisPointingAngle});
      histos.add("XiXiBar/AntiXi/hCascRadius", "hCascRadius", kTH1F, {axisRadius});
      histos.add("XiXiBar/AntiXi/hCascDecayLength", "hCascDecayLength", kTH1F, {axisProperLifeTime});
      histos.add("XiXiBar/AntiXi/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axisMassWindow});
      histos.add("XiXiBar/AntiXi/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axisMassWindow});
      histos.add("XiXiBar/AntiXi/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axisXiMass, axisOmegaMass});
      histos.add("XiXiBar/AntiXi/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("XiXiBar/AntiXi/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("XiXiBar/AntiXi/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("XiXiBar/AntiXi/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("XiXiBar/AntiXi/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("XiXiBar/AntiXi/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      if (doMCAssociation) {
        histos.add("XiXiBar/h3dInvMassTrueEtaC1S", "h3dInvMassTrueEtaC1S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueJPsi", "h3dInvMassTrueJPsi", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueChiC0", "h3dInvMassTrueChiC0", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueChiC1", "h3dInvMassTrueChiC1", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueHC", "h3dInvMassTrueHC", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueChiC2", "h3dInvMassTrueChiC2", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueEtaC2S", "h3dInvMassTrueEtaC2S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTruePsi2S", "h3dInvMassTruePsi2S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      }
    }
    if (buildOmOmBarPairs) {
      histos.add("OmOmBar/h3dMassOmOmbar", "h3dMassOmOmbar", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      if (!isPP) {
        // Non-UPC info
        histos.add("OmOmBar/h3dMassOmOmbarHadronic", "h3dMassOmOmbarHadronic", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        // UPC info
        histos.add("OmOmBar/h3dMassOmOmbarSGA", "h3dMassOmOmbarSGA", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("OmOmBar/h3dMassOmOmbarSGC", "h3dMassOmOmbarSGC", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("OmOmBar/h3dMassOmOmbarDG", "h3dMassOmOmbarDG", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      }
      histos.add("OmOmBar/h2dNbrOfOmegaVsCentrality", "h2dNbrOfOmegaVsCentrality", kTH2F, {axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("OmOmBar/h2dNbrOfAntiOmegaVsCentrality", "h2dNbrOfAntiOmegaVsCentrality", kTH2F, {axisCentrality, {10, -0.5f, 9.5f}});
      // QA plot
      // Candidates after Omega selections
      histos.add("OmOmBar/Omega/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("OmOmBar/Omega/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("OmOmBar/Omega/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("OmOmBar/Omega/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axisDCAdau});
      histos.add("OmOmBar/Omega/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axisDCAdau});
      histos.add("OmOmBar/Omega/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
      histos.add("OmOmBar/Omega/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axisPointingAngle});
      histos.add("OmOmBar/Omega/hV0Radius", "hV0Radius", kTH1F, {axisRadius});
      histos.add("OmOmBar/Omega/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axisPointingAngle});
      histos.add("OmOmBar/Omega/hCascRadius", "hCascRadius", kTH1F, {axisRadius});
      histos.add("OmOmBar/Omega/hCascDecayLength", "hCascDecayLength", kTH1F, {axisProperLifeTime});
      histos.add("OmOmBar/Omega/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axisMassWindow});
      histos.add("OmOmBar/Omega/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axisMassWindow});
      histos.add("OmOmBar/Omega/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axisXiMass, axisOmegaMass});
      histos.add("OmOmBar/Omega/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("OmOmBar/Omega/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("OmOmBar/Omega/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("OmOmBar/Omega/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("OmOmBar/Omega/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("OmOmBar/Omega/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      // Candidates after AntiOmega selections
      histos.add("OmOmBar/AntiOmega/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("OmOmBar/AntiOmega/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("OmOmBar/AntiOmega/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add("OmOmBar/AntiOmega/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axisDCAdau});
      histos.add("OmOmBar/AntiOmega/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axisDCAdau});
      histos.add("OmOmBar/AntiOmega/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axisDCAV0ToPV});
      histos.add("OmOmBar/AntiOmega/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axisPointingAngle});
      histos.add("OmOmBar/AntiOmega/hV0Radius", "hV0Radius", kTH1F, {axisRadius});
      histos.add("OmOmBar/AntiOmega/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axisPointingAngle});
      histos.add("OmOmBar/AntiOmega/hCascRadius", "hCascRadius", kTH1F, {axisRadius});
      histos.add("OmOmBar/AntiOmega/hCascDecayLength", "hCascDecayLength", kTH1F, {axisProperLifeTime});
      histos.add("OmOmBar/AntiOmega/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axisMassWindow});
      histos.add("OmOmBar/AntiOmega/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axisMassWindow});
      histos.add("OmOmBar/AntiOmega/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axisXiMass, axisOmegaMass});
      histos.add("OmOmBar/AntiOmega/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("OmOmBar/AntiOmega/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("OmOmBar/AntiOmega/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axisNsigmaTPC});
      histos.add("OmOmBar/AntiOmega/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("OmOmBar/AntiOmega/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add("OmOmBar/AntiOmega/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      if (doMCAssociation) {
        histos.add("OmOmBar/h3dInvMassTrueEtaC2S", "h3dInvMassTrueEtaC2S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
        histos.add("OmOmBar/h3dInvMassTruePsi2S", "h3dInvMassTruePsi2S", kTH3F, {axisCentrality, axisPt, axisQuarkoniumMass});
      }
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

    // machine learning initialization if requested
    if (mlConfigurations.calculateK0ShortScores ||
        mlConfigurations.calculateLambdaScores ||
        mlConfigurations.calculateAntiLambdaScores) {
      int64_t timeStampML = collision.timestamp();
      if (mlConfigurations.timestampCCDB.value != -1)
        timeStampML = mlConfigurations.timestampCCDB.value;
      LoadMachines(timeStampML);
    }
  }

  // function to load models for ML-based classifiers
  void LoadMachines(int64_t timeStampML)
  {
    if (mlConfigurations.loadCustomModelsFromCCDB) {
      ccdbApi.init(ccdbConfigurations.ccdburl);
      LOG(info) << "Fetching models for timestamp: " << timeStampML;

      if (mlConfigurations.calculateLambdaScores) {
        bool retrieveSuccessLambda = ccdbApi.retrieveBlob(mlConfigurations.customModelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathLambda.value);
        if (retrieveSuccessLambda) {
          mlCustomModelLambda.initModel(mlConfigurations.localModelPathLambda.value, mlConfigurations.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Lambda model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (mlConfigurations.calculateAntiLambdaScores) {
        bool retrieveSuccessAntiLambda = ccdbApi.retrieveBlob(mlConfigurations.customModelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathAntiLambda.value);
        if (retrieveSuccessAntiLambda) {
          mlCustomModelAntiLambda.initModel(mlConfigurations.localModelPathAntiLambda.value, mlConfigurations.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the AntiLambda model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (mlConfigurations.calculateK0ShortScores) {
        bool retrieveSuccessKZeroShort = ccdbApi.retrieveBlob(mlConfigurations.customModelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathK0Short.value);
        if (retrieveSuccessKZeroShort) {
          mlCustomModelK0Short.initModel(mlConfigurations.localModelPathK0Short.value, mlConfigurations.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the K0Short model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }
    } else {
      if (mlConfigurations.calculateLambdaScores)
        mlCustomModelLambda.initModel(mlConfigurations.localModelPathLambda.value, mlConfigurations.enableOptimizations.value);
      if (mlConfigurations.calculateAntiLambdaScores)
        mlCustomModelAntiLambda.initModel(mlConfigurations.localModelPathAntiLambda.value, mlConfigurations.enableOptimizations.value);
      if (mlConfigurations.calculateK0ShortScores)
        mlCustomModelK0Short.initModel(mlConfigurations.localModelPathK0Short.value, mlConfigurations.enableOptimizations.value);
    }
    LOG(info) << "ML Models loaded.";
  }

  template <typename TCollision>
  bool IsEventAccepted(TCollision collision, bool fillHists)
  // check whether the collision passes our collision selections
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

    if (std::abs(collision.posZ()) > 10.f) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2 /* vertex-Z selected */);

    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);

    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);

    if (requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5 /* Contains at least one ITS-TPC track */);

    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6 /* PV position consistency check */);

    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7 /* PV with at least one contributor matched with TOF */);

    if (requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TRD */);

    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9 /* Not at same bunch pile-up */);

    if (requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 10 /* No other collision within +/- 10 microseconds */);

    if (requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 4 microseconds */);

    if (minOccupancy > 0 && collision.trackOccupancyInTimeRange() < minOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 12 /* Below min occupancy */);
    if (maxOccupancy > 0 && collision.trackOccupancyInTimeRange() > maxOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 13 /* Above max occupancy */);

    return true;
  }

  template <typename TCollision>
  void fillEventHistograms(TCollision collision, float& centrality, int& selGapSide)
  {
    if (isPP) { //
      centrality = collision.centFT0M();

      if (qaCentrality) {
        auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
        centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(collision.multFT0A() + collision.multFT0C()));
      }
    } else {
      centrality = collision.centFT0C();

      if (qaCentrality) {
        auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
        centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(collision.multFT0C()));
      }
    }

    // in case we want to push the analysis to Pb-Pb UPC
    int gapSide = collision.gapSide();
    if (!isPP) {
      // -1 --> Hadronic
      // 0 --> Single Gap - A side
      // 1 --> Single Gap - C side
      // 2 --> Double Gap - both A & C sides
      selGapSide = sgSelector.trueGap(collision, upcCuts.FV0cut, upcCuts.FT0Acut, upcCuts.FT0Ccut, upcCuts.ZDCcut);
      histos.fill(HIST("hGapSide"), gapSide);
      histos.fill(HIST("hSelGapSide"), selGapSide);
      histos.fill(HIST("hEventCentralityVsSelGapSide"), centrality, selGapSide <= 2 ? selGapSide : -1);
    }

    histos.fill(HIST("hEventCentrality"), centrality);

    histos.fill(HIST("hCentralityVsNch"), centrality, collision.multNTracksPVeta1());

    histos.fill(HIST("hEventOccupancy"), collision.trackOccupancyInTimeRange());
    histos.fill(HIST("hCentralityVsOccupancy"), centrality, collision.trackOccupancyInTimeRange());

    return;
  }

  template <typename TV0, typename TCollision>
  uint64_t computeReconstructionBitmap(TV0 v0, TCollision collision, float rapidityLambda, float rapidityK0Short, float /*pT*/)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;

    //
    // Base topological variables
    //

    // v0 radius min/max selections
    if (v0.v0radius() > v0Selections.v0radius)
      bitset(bitMap, selRadius);
    if (v0.v0radius() < v0Selections.v0radiusMax)
      bitset(bitMap, selRadiusMax);
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (TMath::Abs(v0.dcapostopv()) > v0Selections.dcaprotontopv &&
        TMath::Abs(v0.dcanegtopv()) > v0Selections.dcapiontopv) {
      bitset(bitMap, selDCAPosToPV);
      bitset(bitMap, selDCANegToPV);
    } else if (TMath::Abs(v0.dcapostopv()) > v0Selections.dcapiontopv &&
               TMath::Abs(v0.dcanegtopv()) > v0Selections.dcaprotontopv) {
      bitset(bitMap, selDCAPosToPV);
      bitset(bitMap, selDCANegToPV);
    }
    // V0 cosine of pointing angle
    if (v0.v0cosPA() > v0Selections.v0cospa)
      bitset(bitMap, selCosPA);
    // DCA between v0 daughters
    if (v0.dcaV0daughters() < v0Selections.dcav0dau)
      bitset(bitMap, selDCAV0Dau);
    // DCA V0 to prim vtx
    if (v0.dcav0topv() > v0Selections.dcav0topv)
      bitset(bitMap, selDCAV0ToPV);

    //
    // rapidity
    //
    if (TMath::Abs(rapidityLambda) < v0Selections.rapidityCut)
      bitset(bitMap, selLambdaRapidity);
    if (TMath::Abs(rapidityK0Short) < v0Selections.rapidityCut)
      bitset(bitMap, selK0ShortRapidity);

    //
    // invariant mass window
    //
    if (TMath::Abs(v0.mK0Short() - pdgDB->Mass(310)) < v0Selections.v0MassWindow)
      bitset(bitMap, selK0ShortMassWindow);
    if (TMath::Abs(v0.mLambda() - pdgDB->Mass(3122)) < v0Selections.v0MassWindow)
      bitset(bitMap, selLambdaMassWindow);
    if (TMath::Abs(v0.mAntiLambda() - pdgDB->Mass(3122)) < v0Selections.v0MassWindow)
      bitset(bitMap, selAntiLambdaMassWindow);

    //
    // competing mass rejection
    //
    if (TMath::Abs(v0.mK0Short() - pdgDB->Mass(310)) > v0Selections.compMassRejection)
      bitset(bitMap, selK0ShortMassRejection);
    if (TMath::Abs(v0.mLambda() - pdgDB->Mass(3122)) > v0Selections.compMassRejection)
      bitset(bitMap, selLambdaMassRejection);

    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

    //
    // ITS quality flags
    //
    if (posTrackExtra.itsNCls() >= v0Selections.minITSclusters)
      bitset(bitMap, selPosGoodITSTrack);
    if (negTrackExtra.itsNCls() >= v0Selections.minITSclusters)
      bitset(bitMap, selNegGoodITSTrack);

    //
    // TPC quality flags
    //
    if (posTrackExtra.tpcCrossedRows() >= v0Selections.minTPCrows)
      bitset(bitMap, selPosGoodTPCTrack);
    if (negTrackExtra.tpcCrossedRows() >= v0Selections.minTPCrows)
      bitset(bitMap, selNegGoodTPCTrack);

    //
    // TPC PID
    //
    if (fabs(posTrackExtra.tpcNSigmaPi()) < v0Selections.TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDPositivePion);
    if (fabs(posTrackExtra.tpcNSigmaPr()) < v0Selections.TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDPositiveProton);
    if (fabs(negTrackExtra.tpcNSigmaPi()) < v0Selections.TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDNegativePion);
    if (fabs(negTrackExtra.tpcNSigmaPr()) < v0Selections.TpcPidNsigmaCut)
      bitset(bitMap, selTPCPIDNegativeProton);

    //
    // TOF PID in DeltaT
    // Positive track
    if (fabs(v0.posTOFDeltaTLaPr()) < v0Selections.maxDeltaTimeProton)
      bitset(bitMap, selTOFDeltaTPositiveProtonLambda);
    if (fabs(v0.posTOFDeltaTLaPi()) < v0Selections.maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTPositivePionLambda);
    if (fabs(v0.posTOFDeltaTK0Pi()) < v0Selections.maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTPositivePionK0Short);
    // Negative track
    if (fabs(v0.negTOFDeltaTLaPr()) < v0Selections.maxDeltaTimeProton)
      bitset(bitMap, selTOFDeltaTNegativeProtonLambda);
    if (fabs(v0.negTOFDeltaTLaPi()) < v0Selections.maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTNegativePionLambda);
    if (fabs(v0.negTOFDeltaTK0Pi()) < v0Selections.maxDeltaTimePion)
      bitset(bitMap, selTOFDeltaTNegativePionK0Short);

    //
    // TOF PID in NSigma
    // Positive track
    if (fabs(v0.tofNSigmaLaPr()) < v0Selections.TofPidNsigmaCutLaPr)
      bitset(bitMap, selTOFNSigmaPositiveProtonLambda);
    if (fabs(v0.tofNSigmaALaPi()) < v0Selections.TofPidNsigmaCutLaPi)
      bitset(bitMap, selTOFNSigmaPositivePionLambda);
    if (fabs(v0.tofNSigmaK0PiPlus()) < v0Selections.TofPidNsigmaCutK0Pi)
      bitset(bitMap, selTOFNSigmaPositivePionK0Short);
    // Negative track
    if (fabs(v0.tofNSigmaALaPr()) < v0Selections.TofPidNsigmaCutLaPr)
      bitset(bitMap, selTOFNSigmaNegativeProtonLambda);
    if (fabs(v0.tofNSigmaLaPi()) < v0Selections.TofPidNsigmaCutLaPi)
      bitset(bitMap, selTOFNSigmaNegativePionLambda);
    if (fabs(v0.tofNSigmaK0PiMinus()) < v0Selections.TofPidNsigmaCutK0Pi)
      bitset(bitMap, selTOFNSigmaNegativePionK0Short);

    //
    // ITS only tag
    if (posTrackExtra.tpcCrossedRows() < 1)
      bitset(bitMap, selPosItsOnly);
    if (negTrackExtra.tpcCrossedRows() < 1)
      bitset(bitMap, selNegItsOnly);

    //
    // TPC only tag
    if (posTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitset(bitMap, selPosNotTPCOnly);
    if (negTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitset(bitMap, selNegNotTPCOnly);

    //
    // proper lifetime
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda"))
      bitset(bitMap, selLambdaCTau);
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S"))
      bitset(bitMap, selK0ShortCTau);

    //
    // armenteros
    if (v0.qtarm() * v0Selections.armPodCut > TMath::Abs(v0.alpha()) || v0Selections.armPodCut < 1e-4)
      bitset(bitMap, selK0ShortArmenteros);

    return bitMap;
  }

  template <typename TCascade, typename TCollision>
  bool isCascadeSelected(TCascade casc, TCollision collision, float rapidity, bool isXi)
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
      if (TMath::Abs(casc.dcapostopv()) < cascSelections.dcaprotontopv)
        return false;
      if (TMath::Abs(casc.dcanegtopv()) < cascSelections.dcapiontopv)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (TMath::Abs(casc.dcapostopv()) < cascSelections.dcapiontopv)
        return false;
      if (TMath::Abs(casc.dcanegtopv()) < cascSelections.dcaprotontopv)
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
    if (TMath::Abs(casc.dcabachtopv()) < cascSelections.dcabachtopv)
      return false;
    // Bachelor-baryon cosPA selection
    if (casc.bachBaryonCosPA() < cascSelections.bachbaryoncospa)
      return false;
    // DCA bachelor-baryon selection
    if (TMath::Abs(casc.bachBaryonDCAxyToPV()) < cascSelections.dcaxybachbaryontopv)
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
    if (TMath::Abs(rapidity) > cascSelections.rapidityCut)
      return false;

    //
    // invariant mass window
    //
    if (TMath::Abs(casc.mLambda() - pdgDB->Mass(3122)) > cascSelections.v0MassWindow)
      return false;
    if (isXi && TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) > cascSelections.cascMassWindow)
      return false;
    if (!isXi && TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) > cascSelections.cascMassWindow)
      return false;

    //
    // competing mass rejection
    //
    if (isXi && TMath::Abs(casc.mOmega() - pdgDB->Mass(3334)) < cascSelections.compMassRejection)
      return false;
    if (!isXi && TMath::Abs(casc.mXi() - pdgDB->Mass(3312)) < cascSelections.compMassRejection)
      return false;

    auto bachTrackExtra = casc.template bachTrackExtra_as<dauTracks>();
    auto posTrackExtra = casc.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = casc.template negTrackExtra_as<dauTracks>();

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
    if (isXi && fabs(bachTrackExtra.tpcNSigmaPi()) > cascSelections.TpcPidNsigmaCut)
      return false;
    if (!isXi && fabs(bachTrackExtra.tpcNSigmaKa()) > cascSelections.TpcPidNsigmaCut)
      return false;
    if (casc.sign() < 0) { // Xi- or Omega- --> positive/negative daughter = proton/pion
      if (fabs(posTrackExtra.tpcNSigmaPr()) > cascSelections.TpcPidNsigmaCut)
        return false;
      if (fabs(negTrackExtra.tpcNSigmaPi()) > cascSelections.TpcPidNsigmaCut)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (fabs(posTrackExtra.tpcNSigmaPi()) > cascSelections.TpcPidNsigmaCut)
        return false;
      if (fabs(negTrackExtra.tpcNSigmaPr()) > cascSelections.TpcPidNsigmaCut)
        return false;
    }

    //
    // TOF PID in DeltaT
    // Bachelor track
    if (bachTrackExtra.hasTOF()) {
      if (isXi && fabs(casc.bachTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
        return false;
      if (!isXi && fabs(casc.bachTOFDeltaTOmKa()) > cascSelections.maxDeltaTimeKaon)
        return false;
    }
    // Positive track
    if (posTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> positive daughter = proton
        if (isXi && fabs(casc.posTOFDeltaTXiPr()) > cascSelections.maxDeltaTimeProton)
          return false;
        if (!isXi && fabs(casc.posTOFDeltaTOmPr()) > cascSelections.maxDeltaTimeProton)
          return false;
      } else { // Xi+ or Omega+ --> positive daughter = pion
        if (isXi && fabs(casc.posTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
          return false;
        if (!isXi && fabs(casc.posTOFDeltaTOmPi()) > cascSelections.maxDeltaTimePion)
          return false;
      }
    }
    // Negative track
    if (negTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> negative daughter = pion
        if (isXi && fabs(casc.negTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
          return false;
        if (!isXi && fabs(casc.negTOFDeltaTOmPi()) > cascSelections.maxDeltaTimePion)
          return false;
      } else { // Xi+ or Omega+ --> negative daughter = proton
        if (isXi && fabs(casc.negTOFDeltaTXiPr()) > cascSelections.maxDeltaTimeProton)
          return false;
        if (!isXi && fabs(casc.negTOFDeltaTOmPr()) > cascSelections.maxDeltaTimeProton)
          return false;
      }
    }

    //
    // TOF PID in NSigma
    // Bachelor track
    if (bachTrackExtra.hasTOF()) {
      if (isXi && fabs(casc.tofNSigmaXiPi()) > cascSelections.TofPidNsigmaCutXiPi)
        return false;
      if (!isXi && fabs(casc.tofNSigmaOmKa()) > cascSelections.TofPidNsigmaCutOmKa)
        return false;
    }
    // Positive track
    if (posTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> positive daughter = proton
        if (isXi && fabs(casc.tofNSigmaXiLaPr()) > cascSelections.TofPidNsigmaCutLaPr)
          return false;
        if (!isXi && fabs(casc.tofNSigmaOmLaPr()) > cascSelections.TofPidNsigmaCutLaPr)
          return false;
      } else { // Xi+ or Omega+ --> positive daughter = pion
        if (isXi && fabs(casc.tofNSigmaXiLaPi()) > cascSelections.TofPidNsigmaCutLaPi)
          return false;
        if (!isXi && fabs(casc.tofNSigmaOmLaPi()) > cascSelections.TofPidNsigmaCutLaPi)
          return false;
      }
    }
    // Negative track
    if (negTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> negative daughter = pion
        if (isXi && fabs(casc.tofNSigmaXiLaPr()) > cascSelections.TofPidNsigmaCutLaPi)
          return false;
        if (!isXi && fabs(casc.tofNSigmaOmLaPr()) > cascSelections.TofPidNsigmaCutLaPi)
          return false;
      } else { // Xi+ or Omega+ --> negative daughter = proton
        if (isXi && fabs(casc.tofNSigmaXiLaPi()) > cascSelections.TofPidNsigmaCutLaPr)
          return false;
        if (!isXi && fabs(casc.tofNSigmaOmLaPi()) > cascSelections.TofPidNsigmaCutLaPr)
          return false;
      }
    }

    //
    // proper lifetime
    float distOverTotMom = std::sqrt(std::pow(casc.x() - collision.posX(), 2) + std::pow(casc.y() - collision.posY(), 2) + std::pow(casc.z() - collision.posZ(), 2)) / (casc.p() + 1E-10);
    if (isXi && distOverTotMom * o2::constants::physics::MassXiMinus / ctauXiPDG > cascSelections.cascProperLifeTime)
      return false;
    if (!isXi && distOverTotMom * o2::constants::physics::MassOmegaMinus / ctauOmegaPDG > cascSelections.cascProperLifeTime)
      return false;

    //
    // MC association (if asked)
    if (doMCAssociation) {
      if constexpr (requires { casc.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>(); }) { // check if MC information is available
        auto cascMC = casc.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

        if (isXi) {
          if (casc.sign() < 0) {
            if (cascMC.pdgCode() != 3312 || cascMC.pdgCodePositive() != 2212 || cascMC.pdgCodeNegative() != -211 || cascMC.pdgCodeBachelor() != -211)
              return false;
          } else {
            if (cascMC.pdgCode() != -3312 || cascMC.pdgCodePositive() != 211 || cascMC.pdgCodeNegative() != -2212 || cascMC.pdgCodeBachelor() != 211)
              return false;
          }
        } else {
          if (casc.sign() < 0) {
            if (cascMC.pdgCode() != 3334 || cascMC.pdgCodePositive() != 2212 || cascMC.pdgCodeNegative() != -211 || cascMC.pdgCodeBachelor() != -321)
              return false;
          } else {
            if (cascMC.pdgCode() != -3334 || cascMC.pdgCodePositive() != 211 || cascMC.pdgCodeNegative() != -2212 || cascMC.pdgCodeBachelor() != 321)
              return false;
          }
        }
      }
    }

    return true;
  }

  template <typename TV0>
  uint64_t computeMCAssociation(TV0 v0)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    // check for specific particle species

    if (v0.pdgCode() == 310 && v0.pdgCodePositive() == 211 && v0.pdgCodeNegative() == -211) {
      bitset(bitMap, selConsiderK0Short);
      if (v0.isPhysicalPrimary())
        bitset(bitMap, selPhysPrimK0Short);
    }
    if (v0.pdgCode() == 3122 && v0.pdgCodePositive() == 2212 && v0.pdgCodeNegative() == -211) {
      bitset(bitMap, selConsiderLambda);
      if (v0.isPhysicalPrimary())
        bitset(bitMap, selPhysPrimLambda);
    }
    if (v0.pdgCode() == -3122 && v0.pdgCodePositive() == 211 && v0.pdgCodeNegative() == -2212) {
      bitset(bitMap, selConsiderAntiLambda);
      if (v0.isPhysicalPrimary())
        bitset(bitMap, selPhysPrimAntiLambda);
    }
    return bitMap;
  }

  bool verifyMask(uint64_t bitmap, uint64_t mask)
  {
    return (bitmap & mask) == mask;
  }

  template <typename TV0>
  void analyseV0Candidate(TV0 v0, float pt, float /*centrality*/, uint64_t selMap, std::vector<bool>& selK0ShortIndices, std::vector<bool>& selLambdaIndices, std::vector<bool>& selAntiLambdaIndices, int v0TableOffset)
  // precalculate this information so that a check is one mask operation, not many
  {
    bool passK0ShortSelections = false;
    bool passLambdaSelections = false;
    bool passAntiLambdaSelections = false;

    // machine learning is on, go for calculation of thresholds
    // FIXME THIS NEEDS ADJUSTING
    std::vector<float> inputFeatures{pt, 0.0f, 0.0f, v0.v0radius(), v0.v0cosPA(), v0.dcaV0daughters(), v0.dcapostopv(), v0.dcanegtopv()};

    if (mlConfigurations.useK0ShortScores) {
      float k0shortScore = -1;
      if (mlConfigurations.calculateK0ShortScores) {
        // evaluate machine-learning scores
        float* k0shortProbability = mlCustomModelK0Short.evalModel(inputFeatures);
        k0shortScore = k0shortProbability[1];
      } else {
        k0shortScore = v0.k0ShortBDTScore();
      }
      if (k0shortScore > mlConfigurations.thresholdK0Short.value) {
        passK0ShortSelections = true;
      }
    } else {
      passK0ShortSelections = verifyMask(selMap, maskSelectionK0Short);
    }
    if (mlConfigurations.useLambdaScores) {
      float lambdaScore = -1;
      if (mlConfigurations.calculateLambdaScores) {
        // evaluate machine-learning scores
        float* lambdaProbability = mlCustomModelLambda.evalModel(inputFeatures);
        lambdaScore = lambdaProbability[1];
      } else {
        lambdaScore = v0.lambdaBDTScore();
      }
      if (lambdaScore > mlConfigurations.thresholdK0Short.value) {
        passLambdaSelections = true;
      }
    } else {
      passLambdaSelections = verifyMask(selMap, maskSelectionLambda);
    }
    if (mlConfigurations.useLambdaScores) {
      float antiLambdaScore = -1;
      if (mlConfigurations.calculateAntiLambdaScores) {
        // evaluate machine-learning scores
        float* antilambdaProbability = mlCustomModelAntiLambda.evalModel(inputFeatures);
        antiLambdaScore = antilambdaProbability[1];
      } else {
        antiLambdaScore = v0.antiLambdaBDTScore();
      }
      if (antiLambdaScore > mlConfigurations.thresholdK0Short.value) {
        passAntiLambdaSelections = true;
      }
    } else {
      passAntiLambdaSelections = verifyMask(selMap, maskSelectionAntiLambda);
    }

    // need local index because of the grouping of collisions
    selK0ShortIndices[v0.globalIndex() - v0TableOffset] = passK0ShortSelections;
    selLambdaIndices[v0.globalIndex() - v0TableOffset] = passLambdaSelections;
    selAntiLambdaIndices[v0.globalIndex() - v0TableOffset] = passAntiLambdaSelections;
  }

  template <typename TCollision, typename THyperon>
  void fillQAplot(TCollision collision, THyperon hyperon, THyperon antiHyperon, int type)
  { // fill QA information about hyperon - antihyperon pair
    if (type == 0) {
      if constexpr (requires { hyperon.mK0Short(); antiHyperon.mK0Short(); }) { // check if v0 information is available
        auto posTrackExtraHyperon = hyperon.template posTrackExtra_as<dauTracks>();
        auto negTrackExtraHyperon = hyperon.template negTrackExtra_as<dauTracks>();

        auto posTrackExtraAntiHyperon = antiHyperon.template posTrackExtra_as<dauTracks>();
        auto negTrackExtraAntiHyperon = antiHyperon.template negTrackExtra_as<dauTracks>();

        float hyperonDecayLength = std::sqrt(std::pow(hyperon.x() - collision.posX(), 2) + std::pow(hyperon.y() - collision.posY(), 2) + std::pow(hyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassLambda0 / (hyperon.p() + 1E-10);
        float antiHyperonDecayLength = std::sqrt(std::pow(antiHyperon.x() - collision.posX(), 2) + std::pow(antiHyperon.y() - collision.posY(), 2) + std::pow(antiHyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassLambda0 / (antiHyperon.p() + 1E-10);

        // Candidates after Xi selections
        histos.fill(HIST("LaLaBar/Lambda/hPosDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("LaLaBar/Lambda/hNegDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("LaLaBar/Lambda/hDCAV0Daughters"), hyperon.dcaV0daughters());
        histos.fill(HIST("LaLaBar/Lambda/hDCAV0ToPV"), hyperon.dcav0topv());
        histos.fill(HIST("LaLaBar/Lambda/hV0PointingAngle"), hyperon.v0cosPA());
        histos.fill(HIST("LaLaBar/Lambda/hV0Radius"), hyperon.v0radius());
        histos.fill(HIST("LaLaBar/Lambda/hV0DecayLength"), hyperonDecayLength);
        histos.fill(HIST("LaLaBar/Lambda/hV0InvMassWindow"), hyperon.mLambda() - pdgDB->Mass(3122));
        histos.fill(HIST("LaLaBar/Lambda/h2dCompetingMassRej"), hyperon.mLambda(), hyperon.mK0Short());
        histos.fill(HIST("LaLaBar/Lambda/hPosTPCNsigma"), posTrackExtraHyperon.tpcNSigmaPr());
        histos.fill(HIST("LaLaBar/Lambda/hNegTPCNsigma"), negTrackExtraHyperon.tpcNSigmaPi());
        histos.fill(HIST("LaLaBar/Lambda/h2dPositiveITSvsTPCpts"), posTrackExtraHyperon.tpcCrossedRows(), posTrackExtraHyperon.itsNCls());
        histos.fill(HIST("LaLaBar/Lambda/h2dNegativeITSvsTPCpts"), negTrackExtraHyperon.tpcCrossedRows(), negTrackExtraHyperon.itsNCls());
        // Candidates after AntiXi selections
        histos.fill(HIST("LaLaBar/AntiLambda/hPosDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("LaLaBar/AntiLambda/hNegDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("LaLaBar/AntiLambda/hDCAV0Daughters"), antiHyperon.dcaV0daughters());
        histos.fill(HIST("LaLaBar/AntiLambda/hDCAV0ToPV"), antiHyperon.dcav0topv());
        histos.fill(HIST("LaLaBar/AntiLambda/hV0PointingAngle"), antiHyperon.v0cosPA());
        histos.fill(HIST("LaLaBar/AntiLambda/hV0Radius"), antiHyperon.v0radius());
        histos.fill(HIST("LaLaBar/AntiLambda/hV0DecayLength"), antiHyperonDecayLength);
        histos.fill(HIST("LaLaBar/AntiLambda/hV0InvMassWindow"), antiHyperon.mLambda() - pdgDB->Mass(3122));
        histos.fill(HIST("LaLaBar/AntiLambda/h2dCompetingMassRej"), antiHyperon.mLambda(), antiHyperon.mK0Short());
        histos.fill(HIST("LaLaBar/AntiLambda/hPosTPCNsigma"), posTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("LaLaBar/AntiLambda/hNegTPCNsigma"), negTrackExtraAntiHyperon.tpcNSigmaPr());
        histos.fill(HIST("LaLaBar/AntiLambda/h2dPositiveITSvsTPCpts"), posTrackExtraAntiHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("LaLaBar/AntiLambda/h2dNegativeITSvsTPCpts"), negTrackExtraAntiHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.itsNCls());
      }
    }
    if (type == 1) {
      if constexpr (requires { hyperon.dcabachtopv(); antiHyperon.dcabachtopv(); }) { // check if Cascade information is available
        auto bachTrackExtraHyperon = hyperon.template bachTrackExtra_as<dauTracks>();
        auto posTrackExtraHyperon = hyperon.template posTrackExtra_as<dauTracks>();
        auto negTrackExtraHyperon = hyperon.template negTrackExtra_as<dauTracks>();

        auto bachTrackExtraAntiHyperon = antiHyperon.template bachTrackExtra_as<dauTracks>();
        auto posTrackExtraAntiHyperon = antiHyperon.template posTrackExtra_as<dauTracks>();
        auto negTrackExtraAntiHyperon = antiHyperon.template negTrackExtra_as<dauTracks>();

        float hyperonDecayLength = std::sqrt(std::pow(hyperon.x() - collision.posX(), 2) + std::pow(hyperon.y() - collision.posY(), 2) + std::pow(hyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassXiMinus / (hyperon.p() + 1E-10);
        float antiHyperonDecayLength = std::sqrt(std::pow(antiHyperon.x() - collision.posX(), 2) + std::pow(antiHyperon.y() - collision.posY(), 2) + std::pow(antiHyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassXiMinus / (antiHyperon.p() + 1E-10);

        // Candidates after Xi selections
        histos.fill(HIST("XiXiBar/Xi/hBachDCAToPV"), hyperon.dcabachtopv());
        histos.fill(HIST("XiXiBar/Xi/hPosDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("XiXiBar/Xi/hNegDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("XiXiBar/Xi/hDCACascDaughters"), hyperon.dcacascdaughters());
        histos.fill(HIST("XiXiBar/Xi/hDCAV0Daughters"), hyperon.dcaV0daughters());
        histos.fill(HIST("XiXiBar/Xi/hDCAV0ToPV"), hyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/Xi/hV0PointingAngle"), hyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/Xi/hV0Radius"), hyperon.v0radius());
        histos.fill(HIST("XiXiBar/Xi/hCascPointingAngle"), hyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/Xi/hCascRadius"), hyperon.cascradius());
        histos.fill(HIST("XiXiBar/Xi/hCascDecayLength"), hyperonDecayLength);
        histos.fill(HIST("XiXiBar/Xi/hV0InvMassWindow"), hyperon.mLambda() - pdgDB->Mass(3122));
        histos.fill(HIST("XiXiBar/Xi/hCascInvMassWindow"), hyperon.mXi() - pdgDB->Mass(3312));
        histos.fill(HIST("XiXiBar/Xi/h2dCompetingMassRej"), hyperon.mXi(), hyperon.mOmega());
        histos.fill(HIST("XiXiBar/Xi/hBachTPCNsigma"), bachTrackExtraHyperon.tpcNSigmaPi());
        histos.fill(HIST("XiXiBar/Xi/hPosTPCNsigma"), posTrackExtraHyperon.tpcNSigmaPr());
        histos.fill(HIST("XiXiBar/Xi/hNegTPCNsigma"), negTrackExtraHyperon.tpcNSigmaPi());
        histos.fill(HIST("XiXiBar/Xi/h2dBachelorITSvsTPCpts"), bachTrackExtraHyperon.tpcCrossedRows(), bachTrackExtraHyperon.itsNCls());
        histos.fill(HIST("XiXiBar/Xi/h2dPositiveITSvsTPCpts"), posTrackExtraHyperon.tpcCrossedRows(), posTrackExtraHyperon.itsNCls());
        histos.fill(HIST("XiXiBar/Xi/h2dNegativeITSvsTPCpts"), negTrackExtraHyperon.tpcCrossedRows(), negTrackExtraHyperon.itsNCls());
        // Candidates after AntiXi selections
        histos.fill(HIST("XiXiBar/AntiXi/hBachDCAToPV"), antiHyperon.dcabachtopv());
        histos.fill(HIST("XiXiBar/AntiXi/hPosDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("XiXiBar/AntiXi/hNegDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("XiXiBar/AntiXi/hDCACascDaughters"), antiHyperon.dcacascdaughters());
        histos.fill(HIST("XiXiBar/AntiXi/hDCAV0Daughters"), antiHyperon.dcaV0daughters());
        histos.fill(HIST("XiXiBar/AntiXi/hDCAV0ToPV"), antiHyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/AntiXi/hV0PointingAngle"), antiHyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/AntiXi/hV0Radius"), antiHyperon.v0radius());
        histos.fill(HIST("XiXiBar/AntiXi/hCascPointingAngle"), antiHyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/AntiXi/hCascRadius"), antiHyperon.cascradius());
        histos.fill(HIST("XiXiBar/AntiXi/hCascDecayLength"), antiHyperonDecayLength);
        histos.fill(HIST("XiXiBar/AntiXi/hV0InvMassWindow"), antiHyperon.mLambda() - pdgDB->Mass(3122));
        histos.fill(HIST("XiXiBar/AntiXi/hCascInvMassWindow"), antiHyperon.mXi() - pdgDB->Mass(3312));
        histos.fill(HIST("XiXiBar/AntiXi/h2dCompetingMassRej"), antiHyperon.mXi(), antiHyperon.mOmega());
        histos.fill(HIST("XiXiBar/AntiXi/hBachTPCNsigma"), bachTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("XiXiBar/AntiXi/hPosTPCNsigma"), posTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("XiXiBar/AntiXi/hNegTPCNsigma"), negTrackExtraAntiHyperon.tpcNSigmaPr());
        histos.fill(HIST("XiXiBar/AntiXi/h2dBachelorITSvsTPCpts"), bachTrackExtraAntiHyperon.tpcCrossedRows(), bachTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("XiXiBar/AntiXi/h2dPositiveITSvsTPCpts"), posTrackExtraAntiHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("XiXiBar/AntiXi/h2dNegativeITSvsTPCpts"), negTrackExtraAntiHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.itsNCls());
      }
    }
    if (type == 2) {
      if constexpr (requires { hyperon.dcabachtopv(); antiHyperon.dcabachtopv(); }) { // check if Cascade information is available
        auto bachTrackExtraHyperon = hyperon.template bachTrackExtra_as<dauTracks>();
        auto posTrackExtraHyperon = hyperon.template posTrackExtra_as<dauTracks>();
        auto negTrackExtraHyperon = hyperon.template negTrackExtra_as<dauTracks>();

        auto bachTrackExtraAntiHyperon = antiHyperon.template bachTrackExtra_as<dauTracks>();
        auto posTrackExtraAntiHyperon = antiHyperon.template posTrackExtra_as<dauTracks>();
        auto negTrackExtraAntiHyperon = antiHyperon.template negTrackExtra_as<dauTracks>();

        float hyperonDecayLength = std::sqrt(std::pow(hyperon.x() - collision.posX(), 2) + std::pow(hyperon.y() - collision.posY(), 2) + std::pow(hyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassOmegaMinus / (hyperon.p() + 1E-10);
        float antiHyperonDecayLength = std::sqrt(std::pow(antiHyperon.x() - collision.posX(), 2) + std::pow(antiHyperon.y() - collision.posY(), 2) + std::pow(antiHyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassOmegaMinus / (antiHyperon.p() + 1E-10);

        // Candidates after Omega selections
        histos.fill(HIST("OmOmBar/Omega/hBachDCAToPV"), hyperon.dcabachtopv());
        histos.fill(HIST("OmOmBar/Omega/hPosDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("OmOmBar/Omega/hNegDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("OmOmBar/Omega/hDCACascDaughters"), hyperon.dcacascdaughters());
        histos.fill(HIST("OmOmBar/Omega/hDCAV0Daughters"), hyperon.dcaV0daughters());
        histos.fill(HIST("OmOmBar/Omega/hDCAV0ToPV"), hyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/Omega/hV0PointingAngle"), hyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/Omega/hV0Radius"), hyperon.v0radius());
        histos.fill(HIST("OmOmBar/Omega/hCascPointingAngle"), hyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/Omega/hCascRadius"), hyperon.cascradius());
        histos.fill(HIST("OmOmBar/Omega/hCascDecayLength"), hyperonDecayLength);
        histos.fill(HIST("OmOmBar/Omega/hV0InvMassWindow"), hyperon.mLambda() - pdgDB->Mass(3122));
        histos.fill(HIST("OmOmBar/Omega/hCascInvMassWindow"), hyperon.mOmega() - pdgDB->Mass(3334));
        histos.fill(HIST("OmOmBar/Omega/h2dCompetingMassRej"), hyperon.mXi(), hyperon.mOmega());
        histos.fill(HIST("OmOmBar/Omega/hBachTPCNsigma"), bachTrackExtraHyperon.tpcNSigmaKa());
        histos.fill(HIST("OmOmBar/Omega/hPosTPCNsigma"), posTrackExtraHyperon.tpcNSigmaPr());
        histos.fill(HIST("OmOmBar/Omega/hNegTPCNsigma"), negTrackExtraHyperon.tpcNSigmaPi());
        histos.fill(HIST("OmOmBar/Omega/h2dBachelorITSvsTPCpts"), bachTrackExtraHyperon.tpcCrossedRows(), bachTrackExtraHyperon.itsNCls());
        histos.fill(HIST("OmOmBar/Omega/h2dPositiveITSvsTPCpts"), posTrackExtraHyperon.tpcCrossedRows(), posTrackExtraHyperon.itsNCls());
        histos.fill(HIST("OmOmBar/Omega/h2dNegativeITSvsTPCpts"), negTrackExtraHyperon.tpcCrossedRows(), negTrackExtraHyperon.itsNCls());
        // Candidates after AntiOmega selections
        histos.fill(HIST("OmOmBar/AntiOmega/hBachDCAToPV"), antiHyperon.dcabachtopv());
        histos.fill(HIST("OmOmBar/AntiOmega/hPosDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("OmOmBar/AntiOmega/hNegDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("OmOmBar/AntiOmega/hDCACascDaughters"), antiHyperon.dcacascdaughters());
        histos.fill(HIST("OmOmBar/AntiOmega/hDCAV0Daughters"), antiHyperon.dcaV0daughters());
        histos.fill(HIST("OmOmBar/AntiOmega/hDCAV0ToPV"), antiHyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/AntiOmega/hV0PointingAngle"), antiHyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/AntiOmega/hV0Radius"), antiHyperon.v0radius());
        histos.fill(HIST("OmOmBar/AntiOmega/hCascPointingAngle"), antiHyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/AntiOmega/hCascRadius"), antiHyperon.cascradius());
        histos.fill(HIST("OmOmBar/AntiOmega/hCascDecayLength"), antiHyperonDecayLength);
        histos.fill(HIST("OmOmBar/AntiOmega/hV0InvMassWindow"), antiHyperon.mLambda() - pdgDB->Mass(3122));
        histos.fill(HIST("OmOmBar/AntiOmega/hCascInvMassWindow"), antiHyperon.mOmega() - pdgDB->Mass(3334));
        histos.fill(HIST("OmOmBar/AntiOmega/h2dCompetingMassRej"), antiHyperon.mXi(), antiHyperon.mOmega());
        histos.fill(HIST("OmOmBar/AntiOmega/hBachTPCNsigma"), bachTrackExtraAntiHyperon.tpcNSigmaKa());
        histos.fill(HIST("OmOmBar/AntiOmega/hPosTPCNsigma"), posTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("OmOmBar/AntiOmega/hNegTPCNsigma"), negTrackExtraAntiHyperon.tpcNSigmaPr());
        histos.fill(HIST("OmOmBar/AntiOmega/h2dBachelorITSvsTPCpts"), bachTrackExtraAntiHyperon.tpcCrossedRows(), bachTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("OmOmBar/AntiOmega/h2dPositiveITSvsTPCpts"), posTrackExtraAntiHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("OmOmBar/AntiOmega/h2dNegativeITSvsTPCpts"), negTrackExtraAntiHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.itsNCls());
      }
    }
  }

  template <typename TCollision, typename THyperon>
  void analyseHyperonPairCandidate(TCollision collision, THyperon hyperon, THyperon antiHyperon, float centrality, uint8_t gapSide, int type)
  // fill information related to the quarkonium mother
  // type = 0 (Lambda), 1 (Xi), 2 (Omega)
  {
    float pt = RecoDecay::pt(hyperon.px() + antiHyperon.px(), hyperon.py() + antiHyperon.py());

    float invmass = -1;
    if (type == 0)
      invmass = RecoDecay::m(std::array{std::array{hyperon.px(), hyperon.py(), hyperon.pz()}, std::array{antiHyperon.px(), antiHyperon.py(), antiHyperon.pz()}}, std::array{o2::constants::physics::MassLambda0, o2::constants::physics::MassLambda0Bar});
    if (type == 1)
      invmass = RecoDecay::m(std::array{std::array{hyperon.px(), hyperon.py(), hyperon.pz()}, std::array{antiHyperon.px(), antiHyperon.py(), antiHyperon.pz()}}, std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassXiPlusBar});
    if (type == 2)
      invmass = RecoDecay::m(std::array{std::array{hyperon.px(), hyperon.py(), hyperon.pz()}, std::array{antiHyperon.px(), antiHyperon.py(), antiHyperon.pz()}}, std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassOmegaPlusBar});

    float rapidity = RecoDecay::y(std::array{hyperon.px() + antiHyperon.px(), hyperon.py() + antiHyperon.py(), hyperon.pz() + antiHyperon.pz()}, invmass);

    // rapidity cut on the quarkonium mother
    if (!doMCAssociation && TMath::Abs(rapidity) > rapidityCut)
      return;

    // fillV0sInfo(lambda, antiLambda, centrality);

    // __________________________________________
    // main analysis
    if (type == 0) {
      if (doMCAssociation) {
        if constexpr (requires { hyperon.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>(); }) { // check if MC information is available
          auto hyperonMC = hyperon.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
          auto antiHyperonMC = antiHyperon.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

          if (hyperonMC.pdgCodeMother() != antiHyperonMC.pdgCodeMother()) {
            return;
          }

          float ptmc = RecoDecay::pt(hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC());
          float rapiditymc = RecoDecay::y(std::array{hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC(), hyperonMC.pzMC() + antiHyperonMC.pzMC()}, pdgDB->Mass(hyperonMC.pdgCodeMother()));

          if (TMath::Abs(rapiditymc) > rapidityCut)
            return;

          if (hyperonMC.pdgCodeMother() == 441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // EtaC(1S)
            histos.fill(HIST("LaLaBar/h3dInvMassTrueEtaC1S"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // J/psi
            histos.fill(HIST("LaLaBar/h3dInvMassTrueJPsi"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 10441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC0
            histos.fill(HIST("LaLaBar/h3dInvMassTrueChiC0"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 20443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC1
            histos.fill(HIST("LaLaBar/h3dInvMassTrueChiC1"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 10443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // hC
            histos.fill(HIST("LaLaBar/h3dInvMassTrueHC"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 445 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC2
            histos.fill(HIST("LaLaBar/h3dInvMassTrueChiC2"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 100441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // EtaC(2S)
            histos.fill(HIST("LaLaBar/h3dInvMassTrueEtaC2S"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 100443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // Psi(2S)
            histos.fill(HIST("LaLaBar/h3dInvMassTruePsi2S"), centrality, ptmc, invmass);
          }
        }
      }

      histos.fill(HIST("LaLaBar/h3dMassLaLabar"), centrality, pt, invmass);
      if (!isPP) { // in case of PbPb data
        if (gapSide == 0)
          histos.fill(HIST("LaLaBar/h3dMassLaLabarSGA"), centrality, pt, invmass);
        else if (gapSide == 1)
          histos.fill(HIST("LaLaBar/h3dMassLaLabarSGC"), centrality, pt, invmass);
        else if (gapSide == 2)
          histos.fill(HIST("LaLaBar/h3dMassLaLabarDG"), centrality, pt, invmass);
        else
          histos.fill(HIST("LaLaBar/h3dMassLaLabarHadronic"), centrality, pt, invmass);
      }
      fillQAplot(collision, hyperon, antiHyperon, 0);
    }
    if (type == 1) {
      if (doMCAssociation) {
        if constexpr (requires { hyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>(); }) { // check if MC information is available
          auto hyperonMC = hyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();
          auto antiHyperonMC = antiHyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

          if (hyperonMC.pdgCodeMother() != antiHyperonMC.pdgCodeMother()) {
            return;
          }

          float ptmc = RecoDecay::pt(hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC());
          float rapiditymc = RecoDecay::y(std::array{hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC(), hyperonMC.pzMC() + antiHyperonMC.pzMC()}, pdgDB->Mass(hyperonMC.pdgCodeMother()));

          if (TMath::Abs(rapiditymc) > rapidityCut)
            return;

          if (hyperonMC.pdgCodeMother() == 441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // EtaC(1S)
            histos.fill(HIST("XiXiBar/h3dInvMassTrueEtaC1S"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // J/psi
            histos.fill(HIST("XiXiBar/h3dInvMassTrueJPsi"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 10441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC0
            histos.fill(HIST("XiXiBar/h3dInvMassTrueChiC0"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 20443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC1
            histos.fill(HIST("XiXiBar/h3dInvMassTrueChiC1"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 10443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // hC
            histos.fill(HIST("XiXiBar/h3dInvMassTrueHC"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 445 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC2
            histos.fill(HIST("XiXiBar/h3dInvMassTrueChiC2"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 100441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // EtaC(2S)
            histos.fill(HIST("XiXiBar/h3dInvMassTrueEtaC2S"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 100443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // Psi(2S)
            histos.fill(HIST("XiXiBar/h3dInvMassTruePsi2S"), centrality, ptmc, invmass);
          }
        }
      }

      histos.fill(HIST("XiXiBar/h3dMassXiXibar"), centrality, pt, invmass);
      if (!isPP) { // in case of PbPb data
        if (gapSide == 0)
          histos.fill(HIST("XiXiBar/h3dMassXiXibarSGA"), centrality, pt, invmass);
        else if (gapSide == 1)
          histos.fill(HIST("XiXiBar/h3dMassXiXibarSGC"), centrality, pt, invmass);
        else if (gapSide == 2)
          histos.fill(HIST("XiXiBar/h3dMassXiXibarDG"), centrality, pt, invmass);
        else
          histos.fill(HIST("XiXiBar/h3dMassXiXibarHadronic"), centrality, pt, invmass);
      }
      fillQAplot(collision, hyperon, antiHyperon, 1);
    }
    if (type == 2) {
      if (doMCAssociation) {
        if constexpr (requires { hyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>(); }) { // check if MC information is available
          auto hyperonMC = hyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();
          auto antiHyperonMC = antiHyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

          if (hyperonMC.pdgCodeMother() != antiHyperonMC.pdgCodeMother()) {
            return;
          }

          float ptmc = RecoDecay::pt(hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC());
          float rapiditymc = RecoDecay::y(std::array{hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC(), hyperonMC.pzMC() + antiHyperonMC.pzMC()}, pdgDB->Mass(hyperonMC.pdgCodeMother()));

          if (TMath::Abs(rapiditymc) > rapidityCut)
            return;

          if (hyperonMC.pdgCodeMother() == 100441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // EtaC(2S)
            histos.fill(HIST("OmOmBar/h3dInvMassTrueEtaC2S"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 100443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // Psi(2S)
            histos.fill(HIST("OmOmBar/h3dInvMassTruePsi2S"), centrality, ptmc, invmass);
          }
        }
      }

      histos.fill(HIST("OmOmBar/h3dMassOmOmbar"), centrality, pt, invmass);
      if (!isPP) { // in case of PbPb data
        if (gapSide == 0)
          histos.fill(HIST("OmOmBar/h3dMassOmOmbarSGA"), centrality, pt, invmass);
        else if (gapSide == 1)
          histos.fill(HIST("OmOmBar/h3dMassOmOmbarSGC"), centrality, pt, invmass);
        else if (gapSide == 2)
          histos.fill(HIST("OmOmBar/h3dMassOmOmbarDG"), centrality, pt, invmass);
        else
          histos.fill(HIST("OmOmBar/h3dMassOmOmbarHadronic"), centrality, pt, invmass);
      }
      fillQAplot(collision, hyperon, antiHyperon, 2);
    }
  }

  // function to check that the hyperon and antihyperon have different daughter tracks
  template <typename THyperon>
  bool checkTrackIndices(THyperon hyperon, THyperon antiHyperon)
  {
    if constexpr (requires { hyperon.template bachTrackExtra_as<dauTracks>(); }) { // cascade case: check if bachelor information is available
      // check that bachelor track from hyperon is different from daughter tracks of antiHyperon
      if (hyperon.bachTrackExtraId() == antiHyperon.bachTrackExtraId() ||
          hyperon.bachTrackExtraId() == antiHyperon.posTrackExtraId() ||
          hyperon.bachTrackExtraId() == antiHyperon.negTrackExtraId())
        return false;
      // check that positive track from hyperon is different from daughter tracks of antiHyperon
      if (hyperon.posTrackExtraId() == antiHyperon.bachTrackExtraId() ||
          hyperon.posTrackExtraId() == antiHyperon.posTrackExtraId() ||
          hyperon.posTrackExtraId() == antiHyperon.negTrackExtraId())
        return false;
      // check that negative track from hyperon is different from daughter tracks of antiHyperon
      if (hyperon.negTrackExtraId() == antiHyperon.bachTrackExtraId() ||
          hyperon.negTrackExtraId() == antiHyperon.posTrackExtraId() ||
          hyperon.negTrackExtraId() == antiHyperon.negTrackExtraId())
        return false;
    } else { // v0 case
      // check that positive track from hyperon is different from daughter tracks of antiHyperon
      if (hyperon.posTrackExtraId() == antiHyperon.posTrackExtraId() ||
          hyperon.posTrackExtraId() == antiHyperon.negTrackExtraId())
        return false;
      // check that negative track from hyperon is different from daughter tracks of antiHyperon
      if (hyperon.negTrackExtraId() == antiHyperon.posTrackExtraId() ||
          hyperon.negTrackExtraId() == antiHyperon.negTrackExtraId())
        return false;
    }
    return true;
  }

  template <typename TCollision, typename THyperons>
  void buildHyperonAntiHyperonPairs(TCollision const& collision, THyperons const& fullHyperons, std::vector<bool> selHypIndices, std::vector<bool> selAntiHypIndices, float centrality, uint8_t gapSide, int type)
  {
    // 1st loop over all v0s/cascades
    for (auto& hyperon : fullHyperons) {
      // select only v0s matching Lambda selections
      if (!selHypIndices[hyperon.globalIndex() - fullHyperons.offset()]) { // local index needed due to collisions grouping
        continue;
      }

      // 2nd loop over all v0s/cascade
      for (auto& antiHyperon : fullHyperons) {
        // select only v0s matching Anti-Lambda selections
        if (!selAntiHypIndices[antiHyperon.globalIndex() - fullHyperons.offset()]) { // local index needed due to collisions grouping
          continue;
        }

        // check we don't look at the same v0s/cascades
        if (hyperon.globalIndex() == antiHyperon.globalIndex()) {
          continue;
        }

        // check that the two hyperons have different daughter tracks
        if (!checkTrackIndices(hyperon, antiHyperon)) {
          continue;
        }

        // form V0 pairs and fill histograms
        analyseHyperonPairCandidate(collision, hyperon, antiHyperon, centrality, gapSide, type);
      } // end antiHyperon loop
    } // end hyperon loop

    return;
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator const& collision, v0Candidates const& fullV0s, cascadeCandidates const& fullCascades, dauTracks const&)
  {
    // Fire up CCDB
    if ((mlConfigurations.useK0ShortScores && mlConfigurations.calculateK0ShortScores) ||
        (mlConfigurations.useLambdaScores && mlConfigurations.calculateLambdaScores) ||
        (mlConfigurations.useAntiLambdaScores && mlConfigurations.calculateAntiLambdaScores)) {
      initCCDB(collision);
    }

    if (!IsEventAccepted(collision, true)) {
      return;
    }

    float centrality = -1;
    int selGapSide = -1; // only useful in case one wants to use this task in Pb-Pb UPC
    fillEventHistograms(collision, centrality, selGapSide);

    // __________________________________________
    // perform main analysis
    //
    if (buildLaLaBarPairs) { // Look at V0s
      std::vector<bool> selK0ShortIndices(fullV0s.size());
      std::vector<bool> selLambdaIndices(fullV0s.size());
      std::vector<bool> selAntiLambdaIndices(fullV0s.size());
      for (auto& v0 : fullV0s) {
        if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut || std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        if (v0.v0Type() != v0Selections.v0TypeSelection && v0Selections.v0TypeSelection > -1)
          continue; // skip V0s that are not standard

        uint64_t selMap = computeReconstructionBitmap(v0, collision, v0.yLambda(), v0.yK0Short(), v0.pt());

        // consider for histograms for all species
        selMap = selMap | (uint64_t(1) << selConsiderK0Short) | (uint64_t(1) << selConsiderLambda) | (uint64_t(1) << selConsiderAntiLambda);
        selMap = selMap | (uint64_t(1) << selPhysPrimK0Short) | (uint64_t(1) << selPhysPrimLambda) | (uint64_t(1) << selPhysPrimAntiLambda);

        analyseV0Candidate(v0, v0.pt(), centrality, selMap, selK0ShortIndices, selLambdaIndices, selAntiLambdaIndices, fullV0s.offset());
      } // end v0 loop

      // count the number of K0s, Lambda and AntiLambdas passsing the selections
      int nK0Shorts = std::count(selK0ShortIndices.begin(), selK0ShortIndices.end(), true);
      int nLambdas = std::count(selLambdaIndices.begin(), selLambdaIndices.end(), true);
      int nAntiLambdas = std::count(selAntiLambdaIndices.begin(), selAntiLambdaIndices.end(), true);

      // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
      histos.fill(HIST("LaLaBar/h2dNbrOfK0ShortVsCentrality"), centrality, nK0Shorts);
      histos.fill(HIST("LaLaBar/h2dNbrOfLambdaVsCentrality"), centrality, nLambdas);
      histos.fill(HIST("LaLaBar/h2dNbrOfAntiLambdaVsCentrality"), centrality, nAntiLambdas);

      // Check the number of Lambdas and antiLambdas
      // needs at least 1 of each
      if (nLambdas >= 1 && nAntiLambdas >= 1) {
        buildHyperonAntiHyperonPairs(collision, fullV0s, selLambdaIndices, selAntiLambdaIndices, centrality, selGapSide, 0);
      }
    }

    if (buildXiXiBarPairs || buildOmOmBarPairs) { // Look at Cascades
      std::vector<bool> selXiIndices(fullCascades.size());
      std::vector<bool> selAntiXiIndices(fullCascades.size());
      std::vector<bool> selOmIndices(fullCascades.size());
      std::vector<bool> selAntiOmIndices(fullCascades.size());
      for (auto& cascade : fullCascades) {
        if (std::abs(cascade.negativeeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.positiveeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.bacheloreta()) > cascSelections.daughterEtaCut)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        if (buildXiXiBarPairs) {
          if (cascade.sign() < 0) {
            selXiIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, cascade.yXi(), true);
          } else {
            selAntiXiIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, cascade.yXi(), true);
          }
        }
        if (buildOmOmBarPairs) {
          if (cascade.sign() < 0) {
            selOmIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, cascade.yOmega(), false);
          } else {
            selAntiOmIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, cascade.yOmega(), false);
          }
        }
      } // end cascade loop

      // count the number of Xi and antiXi passsing the selections
      int nXis = std::count(selXiIndices.begin(), selXiIndices.end(), true);
      int nAntiXis = std::count(selAntiXiIndices.begin(), selAntiXiIndices.end(), true);
      int nOmegas = std::count(selOmIndices.begin(), selOmIndices.end(), true);
      int nAntiOmegas = std::count(selAntiOmIndices.begin(), selAntiOmIndices.end(), true);

      // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
      if (buildXiXiBarPairs) {
        histos.fill(HIST("XiXiBar/h2dNbrOfXiVsCentrality"), centrality, nXis);
        histos.fill(HIST("XiXiBar/h2dNbrOfAntiXiVsCentrality"), centrality, nAntiXis);

        // Check the number of Lambdas and antiLambdas
        // needs at least 1 of each
        if (nXis >= 1 && nAntiXis >= 1) {
          buildHyperonAntiHyperonPairs(collision, fullCascades, selXiIndices, selAntiXiIndices, centrality, selGapSide, 1);
        }
      }
      if (buildOmOmBarPairs) {
        histos.fill(HIST("OmOmBar/h2dNbrOfOmegaVsCentrality"), centrality, nOmegas);
        histos.fill(HIST("OmOmBar/h2dNbrOfAntiOmegaVsCentrality"), centrality, nAntiOmegas);

        // Check the number of Lambdas and antiLambdas
        // needs at least 1 of each
        if (nOmegas >= 1 && nAntiOmegas >= 1) {
          buildHyperonAntiHyperonPairs(collision, fullCascades, selOmIndices, selAntiOmIndices, centrality, selGapSide, 2);
        }
      }
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels>::iterator const& collision, v0MCCandidates const& fullV0s, cascadeMCCandidates const& fullCascades, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& /*mccollisions*/, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {
    // Fire up CCDB
    if ((mlConfigurations.useK0ShortScores && mlConfigurations.calculateK0ShortScores) ||
        (mlConfigurations.useLambdaScores && mlConfigurations.calculateLambdaScores) ||
        (mlConfigurations.useAntiLambdaScores && mlConfigurations.calculateAntiLambdaScores)) {
      initCCDB(collision);
    }

    if (!IsEventAccepted(collision, true)) {
      return;
    }

    float centrality = -1;
    int selGapSide = -1; // only useful in case one wants to use this task in Pb-Pb UPC
    fillEventHistograms(collision, centrality, selGapSide);

    // __________________________________________
    // perform main analysis
    if (buildLaLaBarPairs) { // Look at V0s
      std::vector<bool> selK0ShortIndices(fullV0s.size());
      std::vector<bool> selLambdaIndices(fullV0s.size());
      std::vector<bool> selAntiLambdaIndices(fullV0s.size());
      for (auto& v0 : fullV0s) {
        if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut || std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        if (!v0.has_v0MCCore())
          continue;

        auto v0MC = v0.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

        float ptmc = RecoDecay::sqrtSumOfSquares(v0MC.pxPosMC() + v0MC.pxNegMC(), v0MC.pyPosMC() + v0MC.pyNegMC());
        float ymc = 1e-3;
        if (v0MC.pdgCode() == 310)
          ymc = RecoDecay::y(std::array{v0MC.pxPosMC() + v0MC.pxNegMC(), v0MC.pyPosMC() + v0MC.pyNegMC(), v0MC.pzPosMC() + v0MC.pzNegMC()}, o2::constants::physics::MassKaonNeutral);
        else if (TMath::Abs(v0MC.pdgCode()) == 3122)
          ymc = RecoDecay::y(std::array{v0MC.pxPosMC() + v0MC.pxNegMC(), v0MC.pyPosMC() + v0MC.pyNegMC(), v0MC.pzPosMC() + v0MC.pzNegMC()}, o2::constants::physics::MassLambda);

        uint64_t selMap = computeReconstructionBitmap(v0, collision, ymc, ymc, ptmc);
        selMap = selMap | computeMCAssociation(v0MC);

        // consider only associated candidates if asked to do so, disregard association
        if (!doMCAssociation) {
          selMap = selMap | (uint64_t(1) << selConsiderK0Short) | (uint64_t(1) << selConsiderLambda) | (uint64_t(1) << selConsiderAntiLambda);
          selMap = selMap | (uint64_t(1) << selPhysPrimK0Short) | (uint64_t(1) << selPhysPrimLambda) | (uint64_t(1) << selPhysPrimAntiLambda);
        }

        analyseV0Candidate(v0, ptmc, centrality, selMap, selK0ShortIndices, selLambdaIndices, selAntiLambdaIndices, fullV0s.offset());
      } // end v0 loop

      /// count the number of K0s, Lambda and AntiLambdas passsing the selections
      int nK0Shorts = std::count(selK0ShortIndices.begin(), selK0ShortIndices.end(), true);
      int nLambdas = std::count(selLambdaIndices.begin(), selLambdaIndices.end(), true);
      int nAntiLambdas = std::count(selAntiLambdaIndices.begin(), selAntiLambdaIndices.end(), true);

      // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
      histos.fill(HIST("LaLaBar/h2dNbrOfK0ShortVsCentrality"), centrality, nK0Shorts);
      histos.fill(HIST("LaLaBar/h2dNbrOfLambdaVsCentrality"), centrality, nLambdas);
      histos.fill(HIST("LaLaBar/h2dNbrOfAntiLambdaVsCentrality"), centrality, nAntiLambdas);

      if (nLambdas >= 1 && nAntiLambdas >= 1) {
        buildHyperonAntiHyperonPairs(collision, fullV0s, selLambdaIndices, selAntiLambdaIndices, centrality, selGapSide, 0);
      }
    }

    if (buildXiXiBarPairs || buildOmOmBarPairs) { // Look at Cascades
      std::vector<bool> selXiIndices(fullCascades.size());
      std::vector<bool> selAntiXiIndices(fullCascades.size());
      std::vector<bool> selOmIndices(fullCascades.size());
      std::vector<bool> selAntiOmIndices(fullCascades.size());
      for (auto& cascade : fullCascades) {
        if (std::abs(cascade.negativeeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.positiveeta()) > cascSelections.daughterEtaCut ||
            std::abs(cascade.bacheloreta()) > cascSelections.daughterEtaCut)
          continue; // remove acceptance that's badly reproduced by MC / superfluous in future

        if (!cascade.has_cascMCCore())
          continue;

        auto cascadeMC = cascade.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

        // float ptmc = RecoDecay::sqrtSumOfSquares(cascadeMC.pxMC(), cascadeMC.pyMC());
        float ymc = 1e-3;
        if (TMath::Abs(cascadeMC.pdgCode()) == 3312)
          ymc = RecoDecay::y(std::array{cascadeMC.pxMC(), cascadeMC.pyMC(), cascadeMC.pzMC()}, o2::constants::physics::MassXiMinus);
        else if (TMath::Abs(cascadeMC.pdgCode()) == 3334)
          ymc = RecoDecay::y(std::array{cascadeMC.pxMC(), cascadeMC.pyMC(), cascadeMC.pzMC()}, o2::constants::physics::MassOmegaMinus);

        if (buildXiXiBarPairs) {
          if (cascade.sign() < 0) {
            selXiIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, ymc, true);
          } else {
            selAntiXiIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, ymc, true);
          }
        }
        if (buildOmOmBarPairs) {
          if (cascade.sign() < 0) {
            selOmIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, ymc, false);
          } else {
            selAntiOmIndices[cascade.globalIndex() - fullCascades.offset()] = isCascadeSelected(cascade, collision, ymc, false);
          }
        }
      } // end cascade loop

      // count the number of Xi and antiXi passsing the selections
      int nXis = std::count(selXiIndices.begin(), selXiIndices.end(), true);
      int nAntiXis = std::count(selAntiXiIndices.begin(), selAntiXiIndices.end(), true);
      int nOmegas = std::count(selOmIndices.begin(), selOmIndices.end(), true);
      int nAntiOmegas = std::count(selAntiOmIndices.begin(), selAntiOmIndices.end(), true);

      // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
      if (buildXiXiBarPairs) {
        histos.fill(HIST("XiXiBar/h2dNbrOfXiVsCentrality"), centrality, nXis);
        histos.fill(HIST("XiXiBar/h2dNbrOfAntiXiVsCentrality"), centrality, nAntiXis);

        // Check the number of Lambdas and antiLambdas
        // needs at least 1 of each
        if (nXis >= 1 && nAntiXis >= 1) {
          buildHyperonAntiHyperonPairs(collision, fullCascades, selXiIndices, selAntiXiIndices, centrality, selGapSide, 1);
        }
      }
      if (buildOmOmBarPairs) {
        histos.fill(HIST("OmOmBar/h2dNbrOfOmegaVsCentrality"), centrality, nOmegas);
        histos.fill(HIST("OmOmBar/h2dNbrOfAntiOmegaVsCentrality"), centrality, nAntiOmegas);

        // Check the number of Lambdas and antiLambdas
        // needs at least 1 of each
        if (nOmegas >= 1 && nAntiOmegas >= 1) {
          buildHyperonAntiHyperonPairs(collision, fullCascades, selOmIndices, selAntiOmIndices, centrality, selGapSide, 2);
        }
      }
    }
  }

  PROCESS_SWITCH(quarkoniaToHyperons, processRealData, "process as if real data", true);
  PROCESS_SWITCH(quarkoniaToHyperons, processMonteCarlo, "process as if MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<quarkoniaToHyperons>(cfgc)};
}
