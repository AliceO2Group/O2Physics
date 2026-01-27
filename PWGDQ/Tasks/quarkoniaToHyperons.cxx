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
/// \file quarkoniaToHyperons.cxx
/// \brief quarkonia --> hyperon antihyperon analysis task
///
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences & SMI
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & SMI
//
// V0 analysis task
// ================
//
// This code loops over a V0Cores table and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCHelpers.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector3D.h"
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

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using DauMCTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackMCIds, aod::DauTrackTPCPIDs>;
using V0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0K0ShortMLScores>;
// using V0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0MCCores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0MCCollRefs>;
using V0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0K0ShortMLScores>;

using CascadeCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascXiMLScores, aod::CascOmMLScores>;
using CascadeMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascXiMLScores, aod::CascOmMLScores, aod::CascCoreMCLabels>;

// simple checkers, but ensure 64 bit integers
#define BITSET(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define BITCHECK(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

struct QuarkoniaToHyperons {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> isPP{"isPP", true, "If running on pp collision, switch it on true"};
  Configurable<bool> doQA{"doQA", false, "Produce additional QA histograms?"};
  Configurable<bool> doPairPropagationToDCA{"doPairPropagationToDCA", false, "Propagate the hyperon pair to their DCA?"};

  // for running over skimmed dataset
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "If running over skimmed data, switch it on true"};
  Configurable<std::string> cfgSkimmedTrigger{"cfgSkimmedTrigger", "fDoubleXi,fTripleXi,fQuadrupleXi", "(std::string) Comma separated list of triggers of interest"};

  // Custom grouping
  std::vector<std::vector<int>> v0sGrouped;
  std::vector<std::vector<int>> cascadesGrouped;

  // vector of selected V0/cascade indices
  std::vector<int> selK0ShortIndices;
  std::vector<int> selLambdaIndices;
  std::vector<int> selAntiLambdaIndices;
  std::vector<int> selXiIndices;
  std::vector<int> selAntiXiIndices;
  std::vector<int> selOmIndices;
  std::vector<int> selAntiOmIndices;

  // switch on/off event selections
  Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
  Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
  Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
  Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", true, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
  Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};

  Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};

  Configurable<bool> buildK0sK0sPairs{"buildK0sK0sPairs", false, "Build K0s K0s from charmonia decay"};
  Configurable<bool> buildLaLaBarPairs{"buildLaLaBarPairs", false, "Build Lambda antiLambda from charmonia decay"};
  Configurable<bool> buildXiXiBarPairs{"buildXiXiBarPairs", false, "Build Xi antiXi from charmonia decay"};
  Configurable<bool> buildOmOmBarPairs{"buildOmOmBarPairs", false, "Build Omega antiOmega from charmonia decay"};

  Configurable<bool> buildSameSignPairs{"buildSameSignPairs", false, "If true: build same-sign pairs, otherwise consider only opposite-sign pairs"};

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
    Configurable<float> tpcPidNsigmaCut{"v0Selections.tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
    Configurable<float> tofPidNsigmaCutLaPr{"v0Selections.tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
    Configurable<float> tofPidNsigmaCutLaPi{"v0Selections.tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
    Configurable<float> tofPidNsigmaCutK0Pi{"v0Selections.tofPidNsigmaCutK0Pi", 1e+6, "tofPidNsigmaCutK0Pi"};

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
    Configurable<float> tpcPidNsigmaCut{"cascSelections.tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
    Configurable<float> tofPidNsigmaCutLaPr{"cascSelections.tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
    Configurable<float> tofPidNsigmaCutLaPi{"cascSelections.tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
    Configurable<float> tofPidNsigmaCutXiPi{"cascSelections.tofPidNsigmaCutXiPi", 1e+6, "tofPidNsigmaCutXiPi"};
    Configurable<float> tofPidNsigmaCutOmKa{"cascSelections.tofPidNsigmaCutOmKa", 1e+6, "tofPidNsigmaCutOmKa"};

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
    Configurable<float> fv0Cut{"upcCuts.fv0Cut", 100., "FV0A threshold"};
    Configurable<float> ft0aCut{"upcCuts.ft0aCut", 200., "FT0A threshold"};
    Configurable<float> ft0cCut{"upcCuts.ft0cCut", 100., "FT0C threshold"};
    Configurable<float> zdcCut{"upcCuts.zdcCut", 10., "ZDC threshold"};
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

    // manual
    Configurable<bool> useCustomMagField{"ccdbConfigurations.useCustomMagField", false, "Use custom magnetic field value"};
    Configurable<bool> useCustomRunNumber{"ccdbConfigurations.useCustomRunNumber", false, "Use custom run number"};
    Configurable<float> customMagField{"ccdbConfigurations.customMagField", 5.0f, "Manually set magnetic field"};
    Configurable<int> customRunNumber{"ccdbConfigurations.customRunNumber", 544122, "Manually set the run number "};
  } ccdbConfigurations;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  int mRunNumber;
  std::map<std::string, std::string> metadata;
  float magField;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  struct : ConfigurableGroup {
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

    // AP plot axes
    ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
    ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

    // Track quality axes
    ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
    ConfigurableAxis axisITSclus{"axisITSclus", {7, 0.0f, 7.0f}, "N ITS Clusters"};

    // UPC axes
    ConfigurableAxis axisSelGap{"axisSelGap", {4, -1.5, 2.5}, "Gap side"};

    ConfigurableAxis axisHypPairRadius3D{"axisHypPairRadius3D", {500, 0, 2}, "Hyperon pair 3D distance to PV (cm)"};
    ConfigurableAxis axisHypPairRadius2D{"axisHypPairRadius2D", {500, 0, 1}, "Hyperon pair XY distance to PV (cm)"};
    ConfigurableAxis axisHypPairZ{"axisHypPairZ", {500, -2, 2}, "Hyperon pair longitudinal distance to PV (cm)"};
    ConfigurableAxis axisDCAHypPair{"axisDCAHypPair", {500, 0, 2}, "DCA between the hyperon pair (cm)"};
    ConfigurableAxis axisHypPairCosPA{"axisHypPairCosPA", {500, 0, 1}, "Hyperon pair cosine of pointing angle"};
    ConfigurableAxis axisHypPairOpAngle{"axisHypPairOpAngle", {360, 0, o2::constants::math::TwoPI}, "Hyperon pair momentum opening angle (rad)"};
    ConfigurableAxis axisHypPairEta{"axisHypPairEta", {20, -2, 2}, "Hyperon pair pseudo-rapidity"};
    ConfigurableAxis axisHypPairPhi{"axisHypPairPhi", {180, 0.0f, constants::math::TwoPI}, "Hyperon pair azimuthal angle (rad)"};
  } axes;

  o2::base::MatLayerCylSet* lut; // material LUT for DCA fitter
  o2::vertexing::DCAFitterN<2> fitter;

  // helper object
  o2::pwglf::strangenessBuilderHelper straHelper;

  // PDG database
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // For manual sliceBy
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

  enum Selection : uint64_t { selCosPA = 0,
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
    maskTopological = (static_cast<uint64_t>(1) << selCosPA) | (static_cast<uint64_t>(1) << selRadius) | (static_cast<uint64_t>(1) << selDCANegToPV) | (static_cast<uint64_t>(1) << selDCAPosToPV) | (static_cast<uint64_t>(1) << selDCAV0ToPV) | (static_cast<uint64_t>(1) << selDCAV0Dau) | (static_cast<uint64_t>(1) << selRadiusMax);
    maskTopoNoV0Radius = (static_cast<uint64_t>(1) << selCosPA) | (static_cast<uint64_t>(1) << selDCANegToPV) | (static_cast<uint64_t>(1) << selDCAPosToPV) | (static_cast<uint64_t>(1) << selDCAV0ToPV) | (static_cast<uint64_t>(1) << selDCAV0Dau) | (static_cast<uint64_t>(1) << selRadiusMax);
    maskTopoNoDCANegToPV = (static_cast<uint64_t>(1) << selCosPA) | (static_cast<uint64_t>(1) << selRadius) | (static_cast<uint64_t>(1) << selDCAPosToPV) | (static_cast<uint64_t>(1) << selDCAV0ToPV) | (static_cast<uint64_t>(1) << selDCAV0Dau) | (static_cast<uint64_t>(1) << selRadiusMax);
    maskTopoNoDCAPosToPV = (static_cast<uint64_t>(1) << selCosPA) | (static_cast<uint64_t>(1) << selRadius) | (static_cast<uint64_t>(1) << selDCANegToPV) | (static_cast<uint64_t>(1) << selDCAV0ToPV) | (static_cast<uint64_t>(1) << selDCAV0Dau) | (static_cast<uint64_t>(1) << selRadiusMax);
    maskTopoNoCosPA = (static_cast<uint64_t>(1) << selRadius) | (static_cast<uint64_t>(1) << selDCANegToPV) | (static_cast<uint64_t>(1) << selDCAPosToPV) | (static_cast<uint64_t>(1) << selDCAV0ToPV) | (static_cast<uint64_t>(1) << selDCAV0Dau) | (static_cast<uint64_t>(1) << selRadiusMax);
    maskTopoNoDCAV0Dau = (static_cast<uint64_t>(1) << selCosPA) | (static_cast<uint64_t>(1) << selRadius) | (static_cast<uint64_t>(1) << selDCANegToPV) | (static_cast<uint64_t>(1) << selDCAPosToPV) | (static_cast<uint64_t>(1) << selDCAV0ToPV) | (static_cast<uint64_t>(1) << selRadiusMax);
    maskTopoNoDCAV0ToPV = (static_cast<uint64_t>(1) << selCosPA) | (static_cast<uint64_t>(1) << selRadius) | (static_cast<uint64_t>(1) << selDCANegToPV) | (static_cast<uint64_t>(1) << selDCAPosToPV) | (static_cast<uint64_t>(1) << selDCAV0Dau) | (static_cast<uint64_t>(1) << selRadiusMax);

    maskK0ShortSpecific = (static_cast<uint64_t>(1) << selK0ShortRapidity) | (static_cast<uint64_t>(1) << selK0ShortCTau) | (static_cast<uint64_t>(1) << selK0ShortArmenteros) | (static_cast<uint64_t>(1) << selConsiderK0Short) | (static_cast<uint64_t>(1) << selK0ShortMassWindow) | (static_cast<uint64_t>(1) << selLambdaMassRejection);
    maskLambdaSpecific = (static_cast<uint64_t>(1) << selLambdaRapidity) | (static_cast<uint64_t>(1) << selLambdaCTau) | (static_cast<uint64_t>(1) << selConsiderLambda) | (static_cast<uint64_t>(1) << selLambdaMassWindow) | (static_cast<uint64_t>(1) << selK0ShortMassRejection);
    maskAntiLambdaSpecific = (static_cast<uint64_t>(1) << selLambdaRapidity) | (static_cast<uint64_t>(1) << selLambdaCTau) | (static_cast<uint64_t>(1) << selConsiderAntiLambda) | (static_cast<uint64_t>(1) << selAntiLambdaMassWindow) | (static_cast<uint64_t>(1) << selK0ShortMassRejection);

    // ask for specific TPC/TOF PID selections
    maskTrackProperties = 0;
    if (v0Selections.requirePosITSonly) {
      maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << selPosItsOnly) | (static_cast<uint64_t>(1) << selPosGoodITSTrack);
    } else {
      maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << selPosGoodTPCTrack) | (static_cast<uint64_t>(1) << selPosGoodITSTrack);
      // TPC signal is available: ask for positive track PID
      if (v0Selections.tpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << selTPCPIDPositivePion);
        maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << selTPCPIDPositiveProton);
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << selTPCPIDPositivePion);
      }
      // TOF PID
      if (v0Selections.tofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << selTOFNSigmaPositivePionK0Short) | (static_cast<uint64_t>(1) << selTOFDeltaTPositivePionK0Short);
      if (v0Selections.tofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
        maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << selTOFNSigmaPositiveProtonLambda) | (static_cast<uint64_t>(1) << selTOFDeltaTPositiveProtonLambda);
      if (v0Selections.tofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << selTOFNSigmaPositivePionLambda) | (static_cast<uint64_t>(1) << selTOFDeltaTPositivePionLambda);
    }
    if (v0Selections.requireNegITSonly) {
      maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << selNegItsOnly) | (static_cast<uint64_t>(1) << selNegGoodITSTrack);
    } else {
      maskTrackProperties = maskTrackProperties | (static_cast<uint64_t>(1) << selNegGoodTPCTrack) | (static_cast<uint64_t>(1) << selNegGoodITSTrack);
      // TPC signal is available: ask for negative track PID
      if (v0Selections.tpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << selTPCPIDNegativePion);
        maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << selTPCPIDNegativePion);
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << selTPCPIDNegativeProton);
      }
      // TOF PID
      if (v0Selections.tofPidNsigmaCutK0Pi < 1e+5) // safeguard for no cut
        maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << selTOFNSigmaNegativePionK0Short) | (static_cast<uint64_t>(1) << selTOFDeltaTNegativePionK0Short);
      if (v0Selections.tofPidNsigmaCutLaPi < 1e+5) // safeguard for no cut
        maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << selTOFNSigmaNegativePionLambda) | (static_cast<uint64_t>(1) << selTOFDeltaTNegativePionLambda);
      if (v0Selections.tofPidNsigmaCutLaPr < 1e+5) // safeguard for no cut
        maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << selTOFNSigmaNegativeProtonLambda) | (static_cast<uint64_t>(1) << selTOFDeltaTNegativeProtonLambda);
    }

    if (v0Selections.skipTPConly) {
      maskK0ShortSpecific = maskK0ShortSpecific | (static_cast<uint64_t>(1) << selPosNotTPCOnly) | (static_cast<uint64_t>(1) << selNegNotTPCOnly);
      maskLambdaSpecific = maskLambdaSpecific | (static_cast<uint64_t>(1) << selPosNotTPCOnly) | (static_cast<uint64_t>(1) << selNegNotTPCOnly);
      maskAntiLambdaSpecific = maskAntiLambdaSpecific | (static_cast<uint64_t>(1) << selPosNotTPCOnly) | (static_cast<uint64_t>(1) << selNegNotTPCOnly);
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
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(14, "Below min occup.");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(15, "Above max occup.");

    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, +100.0f}});
    histos.add("hCentralityVsNch", "hCentralityVsNch", kTH2F, {axes.axisCentrality, axes.axisNch});

    histos.add("hEventOccupancy", "hEventOccupancy", kTH1F, {axes.axisOccupancy});
    histos.add("hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH2F, {axes.axisCentrality, axes.axisOccupancy});

    if (!isPP) {
      histos.add("hGapSide", "Gap side; Entries", kTH1F, {{5, -0.5, 4.5}});
      histos.add("hSelGapSide", "Selected gap side; Entries", kTH1F, {axes.axisSelGap});
      histos.add("hEventCentralityVsSelGapSide", ";Centrality (%); Selected gap side", kTH2F, {{100, 0.0f, +100.0f}, axes.axisSelGap});
    }

    // for QA and test purposes
    auto hRawCentrality = histos.add<TH1>("hRawCentrality", "hRawCentrality", kTH1F, {axes.axisRawCentrality});

    for (int ii = 1; ii < 101; ii++) {
      float value = 100.5f - static_cast<float>(ii);
      hRawCentrality->SetBinContent(ii, value);
    }

    // histograms versus mass
    if (buildK0sK0sPairs) {
      histos.add("K0sK0s/h3dMassK0sK0s", "h3dMassK0sK0s", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      if (!isPP) {
        // Non-UPC info
        histos.add("K0sK0s/h3dMassK0sK0sHadronic", "h3dMassK0sK0sHadronic", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        // UPC info
        histos.add("K0sK0s/h3dMassK0sK0sSGA", "h3dMassK0sK0sSGA", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dMassK0sK0sSGC", "h3dMassK0sK0sSGC", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dMassK0sK0sDG", "h3dMassK0sK0sDG", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      histos.add("K0sK0s/h2dNbrOfK0ShortVsCentrality", "h2dNbrOfK0ShortVsCentrality", kTH2F, {axes.axisCentrality, {10, -0.5f, 9.5f}});
      // QA plot
      // Candidates after K0s selections
      histos.add("K0sK0s/K0s/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("K0sK0s/K0s/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("K0sK0s/K0s/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axes.axisDCAdau});
      histos.add("K0sK0s/K0s/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axes.axisDCAV0ToPV});
      histos.add("K0sK0s/K0s/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("K0sK0s/K0s/hV0Radius", "hV0Radius", kTH1F, {axes.axisRadius});
      histos.add("K0sK0s/K0s/hV0DecayLength", "hDecayLength", kTH1F, {axes.axisProperLifeTime});
      histos.add("K0sK0s/K0s/hV0InvMassWindow", "hInvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("K0sK0s/K0s/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axes.axisLambdaMass, axes.axisK0Mass});
      histos.add("K0sK0s/K0s/h2dArmenteros", "h2dArmenteros", kTH2F, {axes.axisAPAlpha, axes.axisAPQt});
      histos.add("K0sK0s/K0s/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("K0sK0s/K0s/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("K0sK0s/K0s/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("K0sK0s/K0s/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      if (doMCAssociation) {
        histos.add("K0sK0s/h3dInvMassTrueEtaC1S", "h3dInvMassTrueEtaC1S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dInvMassTrueJPsi", "h3dInvMassTrueJPsi", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dInvMassTrueChiC0", "h3dInvMassTrueChiC0", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dInvMassTrueChiC1", "h3dInvMassTrueChiC1", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dInvMassTrueHC", "h3dInvMassTrueHC", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dInvMassTrueChiC2", "h3dInvMassTrueChiC2", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dInvMassTrueEtaC2S", "h3dInvMassTrueEtaC2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("K0sK0s/h3dInvMassTruePsi2S", "h3dInvMassTruePsi2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      if (doQA) {
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPosDCAToPV", "h3dMassK0sK0sVsPosDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsNegDCAToPV", "h3dMassK0sK0sVsNegDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsDCAV0Daughters", "h3dMassK0sK0sVsDCAV0Daughters", kTH3F, {axes.axisDCAdau, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsDCAV0ToPV", "h3dMassK0sK0sVsDCAV0ToPV", kTH3F, {axes.axisDCAV0ToPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsV0PointingAngle", "h3dMassK0sK0sVsV0PointingAngle", kTH3F, {axes.axisPointingAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsV0Radius", "h3dMassK0sK0sVsV0Radius", kTH3F, {axes.axisRadius, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsDecayLength", "h3dMassK0sK0sVsDecayLength", kTH3F, {axes.axisProperLifeTime, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsInvMassWindow", "h3dMassK0sK0sVsInvMassWindow", kTH3F, {axes.axisMassWindow, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPosTPCNsigma", "h3dMassK0sK0sVsPosTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsNegTPCNsigma", "h3dMassK0sK0sVsNegTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPosITSclusters", "h3dMassK0sK0sVsPosITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsNegITSclusters", "h3dMassK0sK0sVsNegITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPosNbrCrossedRows", "h3dMassK0sK0sVsPosNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsNegNbrCrossedRows", "h3dMassK0sK0sVsNegNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPairRadius3D", "h3dMassK0sK0sVsPairRadius3D", kTH3F, {axes.axisHypPairRadius3D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPairRadius2D", "h3dMassK0sK0sVsPairRadius2D", kTH3F, {axes.axisHypPairRadius2D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPairZ", "h3dMassK0sK0sVsPairZ", kTH3F, {axes.axisHypPairZ, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsDCAPair", "h3dMassK0sK0sVsDCAPair", kTH3F, {axes.axisDCAHypPair, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPairCosPA", "h3dMassK0sK0sVsPairCosPA", kTH3F, {axes.axisHypPairCosPA, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPairOpAngle", "h3dMassK0sK0sVsPairOpAngle", kTH3F, {axes.axisHypPairOpAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPairEta", "h3dMassK0sK0sVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dMassK0sK0sVsPairPhi", "h3dMassK0sK0sVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/K0sK0s/h3dDeltaEtaK0sK0sVsPairEta", "h3dDeltaEtaK0sK0sVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisHypPairEta});
        histos.add("QA/K0sK0s/h3dDeltaPhiK0sK0sVsPairPhi", "h3dDeltaPhiK0sK0sVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisHypPairPhi});
      }
    }
    if (buildLaLaBarPairs) {
      histos.add("LaLaBar/h3dMassLaLabar", "h3dMassLaLabar", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      if (!isPP) {
        // Non-UPC info
        histos.add("LaLaBar/h3dMassLaLabarHadronic", "h3dMassLaLabarHadronic", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        // UPC info
        histos.add("LaLaBar/h3dMassLaLabarSGA", "h3dMassLaLabarSGA", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dMassLaLabarSGC", "h3dMassLaLabarSGC", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dMassLaLabarDG", "h3dMassLaLabarDG", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      histos.add("LaLaBar/h2dNbrOfLambdaVsCentrality", "h2dNbrOfLambdaVsCentrality", kTH2F, {axes.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("LaLaBar/h2dNbrOfAntiLambdaVsCentrality", "h2dNbrOfAntiLambdaVsCentrality", kTH2F, {axes.axisCentrality, {10, -0.5f, 9.5f}});
      // QA plot
      // Candidates after Lambda selections
      histos.add("LaLaBar/Lambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("LaLaBar/Lambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("LaLaBar/Lambda/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axes.axisDCAdau});
      histos.add("LaLaBar/Lambda/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axes.axisDCAV0ToPV});
      histos.add("LaLaBar/Lambda/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("LaLaBar/Lambda/hV0Radius", "hV0Radius", kTH1F, {axes.axisRadius});
      histos.add("LaLaBar/Lambda/hV0DecayLength", "hDecayLength", kTH1F, {axes.axisProperLifeTime});
      histos.add("LaLaBar/Lambda/hV0InvMassWindow", "hInvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("LaLaBar/Lambda/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axes.axisLambdaMass, axes.axisK0Mass});
      histos.add("LaLaBar/Lambda/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("LaLaBar/Lambda/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("LaLaBar/Lambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("LaLaBar/Lambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});

      // Candidates after AntiLambda selections
      histos.add("LaLaBar/AntiLambda/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("LaLaBar/AntiLambda/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("LaLaBar/AntiLambda/hDCAV0Daughters", "hDCADaughters", kTH1F, {axes.axisDCAdau});
      histos.add("LaLaBar/AntiLambda/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axes.axisDCAV0ToPV});
      histos.add("LaLaBar/AntiLambda/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("LaLaBar/AntiLambda/hV0Radius", "hV0Radius", kTH1F, {axes.axisRadius});
      histos.add("LaLaBar/AntiLambda/hV0DecayLength", "hDecayLength", kTH1F, {axes.axisProperLifeTime});
      histos.add("LaLaBar/AntiLambda/hV0InvMassWindow", "hInvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("LaLaBar/AntiLambda/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axes.axisLambdaMass, axes.axisK0Mass});
      histos.add("LaLaBar/AntiLambda/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("LaLaBar/AntiLambda/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("LaLaBar/AntiLambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("LaLaBar/AntiLambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      if (doMCAssociation) {
        histos.add("LaLaBar/h3dInvMassTrueEtaC1S", "h3dInvMassTrueEtaC1S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueJPsi", "h3dInvMassTrueJPsi", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueChiC0", "h3dInvMassTrueChiC0", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueChiC1", "h3dInvMassTrueChiC1", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueHC", "h3dInvMassTrueHC", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueChiC2", "h3dInvMassTrueChiC2", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTrueEtaC2S", "h3dInvMassTrueEtaC2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("LaLaBar/h3dInvMassTruePsi2S", "h3dInvMassTruePsi2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      if (doQA) {
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPosDCAToPV", "h3dMassLaLaBarVsPosDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsNegDCAToPV", "h3dMassLaLaBarVsNegDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsDCAV0Daughters", "h3dMassLaLaBarVsDCAV0Daughters", kTH3F, {axes.axisDCAdau, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsDCAV0ToPV", "h3dMassLaLaBarVsDCAV0ToPV", kTH3F, {axes.axisDCAV0ToPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsV0PointingAngle", "h3dMassLaLaBarVsV0PointingAngle", kTH3F, {axes.axisPointingAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsV0Radius", "h3dMassLaLaBarVsV0Radius", kTH3F, {axes.axisRadius, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsDecayLength", "h3dMassLaLaBarVsDecayLength", kTH3F, {axes.axisProperLifeTime, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsInvMassWindow", "h3dMassLaLaBarVsInvMassWindow", kTH3F, {axes.axisMassWindow, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPosTPCNsigma", "h3dMassLaLaBarVsPosTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsNegTPCNsigma", "h3dMassLaLaBarVsNegTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPosITSclusters", "h3dMassLaLaBarVsPosITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsNegITSclusters", "h3dMassLaLaBarVsNegITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPosNbrCrossedRows", "h3dMassLaLaBarVsPosNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsNegNbrCrossedRows", "h3dMassLaLaBarVsNegNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPairRadius3D", "h3dMassLaLaBarVsPairRadius3D", kTH3F, {axes.axisHypPairRadius3D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPairRadius2D", "h3dMassLaLaBarVsPairRadius2D", kTH3F, {axes.axisHypPairRadius2D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPairZ", "h3dMassLaLaBarVsPairZ", kTH3F, {axes.axisHypPairZ, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsDCAPair", "h3dMassLaLaBarVsDCAPair", kTH3F, {axes.axisDCAHypPair, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPairCosPA", "h3dMassLaLaBarVsPairCosPA", kTH3F, {axes.axisHypPairCosPA, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPairOpAngle", "h3dMassLaLaBarVsPairOpAngle", kTH3F, {axes.axisHypPairOpAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPairEta", "h3dMassLaLaBarVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dMassLaLaBarVsPairPhi", "h3dMassLaLaBarVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dDeltaEtaLaLaBarVsPairEta", "h3dDeltaEtaLaLaBarVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisHypPairEta});
        histos.add("QA/LaLaBar/h3dDeltaPhiLaLaBarVsPairPhi", "h3dDeltaPhiLaLaBarVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisHypPairPhi});
      }
    }
    if (buildXiXiBarPairs) {
      histos.add("XiXiBar/h3dMassXiXibar", "h3dMassXiXibar", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      if (!isPP) {
        // Non-UPC info
        histos.add("XiXiBar/h3dMassXiXibarHadronic", "h3dMassXiXibarHadronic", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        // UPC info
        histos.add("XiXiBar/h3dMassXiXibarSGA", "h3dMassXiXibarSGA", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dMassXiXibarSGC", "h3dMassXiXibarSGC", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dMassXiXibarDG", "h3dMassXiXibarDG", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      histos.add("XiXiBar/h2dNbrOfXiVsCentrality", "h2dNbrOfXiVsCentrality", kTH2F, {axes.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("XiXiBar/h2dNbrOfAntiXiVsCentrality", "h2dNbrOfAntiXiVsCentrality", kTH2F, {axes.axisCentrality, {10, -0.5f, 9.5f}});
      // QA plot
      // Candidates after Xi selections
      histos.add("XiXiBar/Xi/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("XiXiBar/Xi/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("XiXiBar/Xi/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("XiXiBar/Xi/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axes.axisDCAdau});
      histos.add("XiXiBar/Xi/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axes.axisDCAdau});
      histos.add("XiXiBar/Xi/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axes.axisDCAV0ToPV});
      histos.add("XiXiBar/Xi/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("XiXiBar/Xi/hV0Radius", "hV0Radius", kTH1F, {axes.axisRadius});
      histos.add("XiXiBar/Xi/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("XiXiBar/Xi/hCascRadius", "hCascRadius", kTH1F, {axes.axisRadius});
      histos.add("XiXiBar/Xi/hCascDecayLength", "hCascDecayLength", kTH1F, {axes.axisProperLifeTime});
      histos.add("XiXiBar/Xi/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("XiXiBar/Xi/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("XiXiBar/Xi/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axes.axisXiMass, axes.axisOmegaMass});
      histos.add("XiXiBar/Xi/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("XiXiBar/Xi/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("XiXiBar/Xi/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("XiXiBar/Xi/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("XiXiBar/Xi/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("XiXiBar/Xi/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      // Candidates after AntiXi selections
      histos.add("XiXiBar/AntiXi/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("XiXiBar/AntiXi/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("XiXiBar/AntiXi/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("XiXiBar/AntiXi/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axes.axisDCAdau});
      histos.add("XiXiBar/AntiXi/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axes.axisDCAdau});
      histos.add("XiXiBar/AntiXi/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axes.axisDCAV0ToPV});
      histos.add("XiXiBar/AntiXi/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("XiXiBar/AntiXi/hV0Radius", "hV0Radius", kTH1F, {axes.axisRadius});
      histos.add("XiXiBar/AntiXi/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("XiXiBar/AntiXi/hCascRadius", "hCascRadius", kTH1F, {axes.axisRadius});
      histos.add("XiXiBar/AntiXi/hCascDecayLength", "hCascDecayLength", kTH1F, {axes.axisProperLifeTime});
      histos.add("XiXiBar/AntiXi/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("XiXiBar/AntiXi/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("XiXiBar/AntiXi/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axes.axisXiMass, axes.axisOmegaMass});
      histos.add("XiXiBar/AntiXi/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("XiXiBar/AntiXi/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("XiXiBar/AntiXi/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("XiXiBar/AntiXi/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("XiXiBar/AntiXi/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("XiXiBar/AntiXi/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      if (doMCAssociation) {
        histos.add("XiXiBar/h3dInvMassTrueEtaC1S", "h3dInvMassTrueEtaC1S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueJPsi", "h3dInvMassTrueJPsi", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueChiC0", "h3dInvMassTrueChiC0", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueChiC1", "h3dInvMassTrueChiC1", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueHC", "h3dInvMassTrueHC", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueChiC2", "h3dInvMassTrueChiC2", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTrueEtaC2S", "h3dInvMassTrueEtaC2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("XiXiBar/h3dInvMassTruePsi2S", "h3dInvMassTruePsi2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      if (doQA) {
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsBachDCAToPV", "h3dMassXiXiBarVsBachDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPosDCAToPV", "h3dMassXiXiBarVsPosDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsNegDCAToPV", "h3dMassXiXiBarVsNegDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsDCACascDaughters", "h3dMassXiXiBarVsDCACascDaughters", kTH3F, {axes.axisDCAdau, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsDCAV0Daughters", "h3dMassXiXiBarVsDCAV0Daughters", kTH3F, {axes.axisDCAdau, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsDCAV0ToPV", "h3dMassXiXiBarVsDCAV0ToPV", kTH3F, {axes.axisDCAV0ToPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsV0PointingAngle", "h3dMassXiXiBarVsV0PointingAngle", kTH3F, {axes.axisPointingAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsV0Radius", "h3dMassXiXiBarVsV0Radius", kTH3F, {axes.axisRadius, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsCascPointingAngle", "h3dMassXiXiBarVsCascPointingAngle", kTH3F, {axes.axisPointingAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsCascRadius", "h3dMassXiXiBarVsCascRadius", kTH3F, {axes.axisRadius, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsDecayLength", "h3dMassXiXiBarVsDecayLength", kTH3F, {axes.axisProperLifeTime, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsV0InvMassWindow", "h3dMassXiXiBarVsV0InvMassWindow", kTH3F, {axes.axisMassWindow, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsCascInvMassWindow", "h3dMassXiXiBarVsCascInvMassWindow", kTH3F, {axes.axisMassWindow, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsBachTPCNsigma", "h3dMassXiXiBarVsBachTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPosTPCNsigma", "h3dMassXiXiBarVsPosTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsNegTPCNsigma", "h3dMassXiXiBarVsNegTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsBachITSclusters", "h3dMassXiXiBarVsBachITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPosITSclusters", "h3dMassXiXiBarVsPosITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsNegITSclusters", "h3dMassXiXiBarVsNegITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsBachNbrCrossedRows", "h3dMassXiXiBarVsBachNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPosNbrCrossedRows", "h3dMassXiXiBarVsPosNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsNegNbrCrossedRows", "h3dMassXiXiBarVsNegNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPairRadius3D", "h3dMassXiXiBarVsPairRadius3D", kTH3F, {axes.axisHypPairRadius3D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPairRadius2D", "h3dMassXiXiBarVsPairRadius2D", kTH3F, {axes.axisHypPairRadius2D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPairZ", "h3dMassXiXiBarVsPairZ", kTH3F, {axes.axisHypPairZ, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsDCAPair", "h3dMassXiXiBarVsDCAPair", kTH3F, {axes.axisDCAHypPair, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPairCosPA", "h3dMassXiXiBarVsPairCosPA", kTH3F, {axes.axisHypPairCosPA, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPairOpAngle", "h3dMassXiXiBarVsPairOpAngle", kTH3F, {axes.axisHypPairOpAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPairEta", "h3dMassXiXiBarVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/XiXiBar/h3dMassXiXiBarVsPairPhi", "h3dMassXiXiBarVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dDeltaEtaXiXiBarVsPairEta", "h3dDeltaEtaXiXiBarVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisHypPairEta});
        histos.add("QA/LaLaBar/h3dDeltaPhiXiXiBarVsPairPhi", "h3dDeltaPhiXiXiBarVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisHypPairPhi});
      }
    }
    if (buildOmOmBarPairs) {
      histos.add("OmOmBar/h3dMassOmOmbar", "h3dMassOmOmbar", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      if (!isPP) {
        // Non-UPC info
        histos.add("OmOmBar/h3dMassOmOmbarHadronic", "h3dMassOmOmbarHadronic", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        // UPC info
        histos.add("OmOmBar/h3dMassOmOmbarSGA", "h3dMassOmOmbarSGA", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("OmOmBar/h3dMassOmOmbarSGC", "h3dMassOmOmbarSGC", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("OmOmBar/h3dMassOmOmbarDG", "h3dMassOmOmbarDG", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      histos.add("OmOmBar/h2dNbrOfOmegaVsCentrality", "h2dNbrOfOmegaVsCentrality", kTH2F, {axes.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("OmOmBar/h2dNbrOfAntiOmegaVsCentrality", "h2dNbrOfAntiOmegaVsCentrality", kTH2F, {axes.axisCentrality, {10, -0.5f, 9.5f}});
      // QA plot
      // Candidates after Omega selections
      histos.add("OmOmBar/Omega/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("OmOmBar/Omega/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("OmOmBar/Omega/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("OmOmBar/Omega/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axes.axisDCAdau});
      histos.add("OmOmBar/Omega/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axes.axisDCAdau});
      histos.add("OmOmBar/Omega/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axes.axisDCAV0ToPV});
      histos.add("OmOmBar/Omega/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("OmOmBar/Omega/hV0Radius", "hV0Radius", kTH1F, {axes.axisRadius});
      histos.add("OmOmBar/Omega/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("OmOmBar/Omega/hCascRadius", "hCascRadius", kTH1F, {axes.axisRadius});
      histos.add("OmOmBar/Omega/hCascDecayLength", "hCascDecayLength", kTH1F, {axes.axisProperLifeTime});
      histos.add("OmOmBar/Omega/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("OmOmBar/Omega/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("OmOmBar/Omega/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axes.axisXiMass, axes.axisOmegaMass});
      histos.add("OmOmBar/Omega/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("OmOmBar/Omega/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("OmOmBar/Omega/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("OmOmBar/Omega/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("OmOmBar/Omega/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("OmOmBar/Omega/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      // Candidates after AntiOmega selections
      histos.add("OmOmBar/AntiOmega/hBachDCAToPV", "hBachDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("OmOmBar/AntiOmega/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("OmOmBar/AntiOmega/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axes.axisDCAtoPV});
      histos.add("OmOmBar/AntiOmega/hDCACascDaughters", "hDCACascDaughters", kTH1F, {axes.axisDCAdau});
      histos.add("OmOmBar/AntiOmega/hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {axes.axisDCAdau});
      histos.add("OmOmBar/AntiOmega/hDCAV0ToPV", "hDCAV0ToPV", kTH1F, {axes.axisDCAV0ToPV});
      histos.add("OmOmBar/AntiOmega/hV0PointingAngle", "hV0PointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("OmOmBar/AntiOmega/hV0Radius", "hV0Radius", kTH1F, {axes.axisRadius});
      histos.add("OmOmBar/AntiOmega/hCascPointingAngle", "hCascPointingAngle", kTH1F, {axes.axisPointingAngle});
      histos.add("OmOmBar/AntiOmega/hCascRadius", "hCascRadius", kTH1F, {axes.axisRadius});
      histos.add("OmOmBar/AntiOmega/hCascDecayLength", "hCascDecayLength", kTH1F, {axes.axisProperLifeTime});
      histos.add("OmOmBar/AntiOmega/hV0InvMassWindow", "hV0InvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("OmOmBar/AntiOmega/hCascInvMassWindow", "hCascInvMassWindow", kTH1F, {axes.axisMassWindow});
      histos.add("OmOmBar/AntiOmega/h2dCompetingMassRej", "h2dCompetingMassRej", kTH2F, {axes.axisXiMass, axes.axisOmegaMass});
      histos.add("OmOmBar/AntiOmega/hBachTPCNsigma", "hBachTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("OmOmBar/AntiOmega/hPosTPCNsigma", "hPosTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("OmOmBar/AntiOmega/hNegTPCNsigma", "hNegTPCNsigma", kTH1F, {axes.axisNsigmaTPC});
      histos.add("OmOmBar/AntiOmega/h2dBachelorITSvsTPCpts", "h2dBachelorITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("OmOmBar/AntiOmega/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      histos.add("OmOmBar/AntiOmega/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axes.axisTPCrows, axes.axisITSclus});
      if (doMCAssociation) {
        histos.add("OmOmBar/h3dInvMassTrueEtaC2S", "h3dInvMassTrueEtaC2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("OmOmBar/h3dInvMassTruePsi2S", "h3dInvMassTruePsi2S", kTH3F, {axes.axisCentrality, axes.axisPt, axes.axisQuarkoniumMass});
      }
      if (doQA) {
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsBachDCAToPV", "h3dMassOmOmBarVsBachDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPosDCAToPV", "h3dMassOmOmBarVsPosDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsNegDCAToPV", "h3dMassOmOmBarVsNegDCAToPV", kTH3F, {axes.axisDCAtoPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsDCACascDaughters", "h3dMassOmOmBarVsDCACascDaughters", kTH3F, {axes.axisDCAdau, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsDCAV0Daughters", "h3dMassOmOmBarVsDCAV0Daughters", kTH3F, {axes.axisDCAdau, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsDCAV0ToPV", "h3dMassOmOmBarVsDCAV0ToPV", kTH3F, {axes.axisDCAV0ToPV, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsV0PointingAngle", "h3dMassOmOmBarVsV0PointingAngle", kTH3F, {axes.axisPointingAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsV0Radius", "h3dMassOmOmBarVsV0Radius", kTH3F, {axes.axisRadius, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsCascPointingAngle", "h3dMassOmOmBarVsCascPointingAngle", kTH3F, {axes.axisPointingAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsCascRadius", "h3dMassOmOmBarVsCascRadius", kTH3F, {axes.axisRadius, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsDecayLength", "h3dMassOmOmBarVsDecayLength", kTH3F, {axes.axisProperLifeTime, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsV0InvMassWindow", "h3dMassOmOmBarVsV0InvMassWindow", kTH3F, {axes.axisMassWindow, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsCascInvMassWindow", "h3dMassOmOmBarVsCascInvMassWindow", kTH3F, {axes.axisMassWindow, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsBachTPCNsigma", "h3dMassOmOmBarVsBachTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPosTPCNsigma", "h3dMassOmOmBarVsPosTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsNegTPCNsigma", "h3dMassOmOmBarVsNegTPCNsigma", kTH3F, {axes.axisNsigmaTPC, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsBachITSclusters", "h3dMassOmOmBarVsBachITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPosITSclusters", "h3dMassOmOmBarVsPosITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsNegITSclusters", "h3dMassOmOmBarVsNegITSclusters", kTH3F, {axes.axisITSclus, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsBachNbrCrossedRows", "h3dMassOmOmBarVsBachNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPosNbrCrossedRows", "h3dMassOmOmBarVsPosNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsNegNbrCrossedRows", "h3dMassOmOmBarVsNegNbrCrossedRows", kTH3F, {axes.axisTPCrows, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPairRadius3D", "h3dMassOmOmBarVsPairRadius3D", kTH3F, {axes.axisHypPairRadius3D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPairRadius2D", "h3dMassOmOmBarVsPairRadius2D", kTH3F, {axes.axisHypPairRadius2D, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPairZ", "h3dMassOmOmBarVsPairZ", kTH3F, {axes.axisHypPairZ, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsDCAPair", "h3dMassOmOmBarVsDCAPair", kTH3F, {axes.axisDCAHypPair, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPairCosPA", "h3dMassOmOmBarVsPairCosPA", kTH3F, {axes.axisHypPairCosPA, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPairOpAngle", "h3dMassOmOmBarVsPairOpAngle", kTH3F, {axes.axisHypPairOpAngle, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPairEta", "h3dMassOmOmBarVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/OmOmBar/h3dMassOmOmBarVsPairPhi", "h3dMassOmOmBarVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisQuarkoniumMass});
        histos.add("QA/LaLaBar/h3dDeltaEtaOmOmBarVsPairEta", "h3dDeltaEtaOmOmBarVsPairEta", kTH3F, {axes.axisHypPairEta, axes.axisPt, axes.axisHypPairEta});
        histos.add("QA/LaLaBar/h3dDeltaPhiOmOmBarVsPairPhi", "h3dDeltaPhiOmOmBarVsPairPhi", kTH3F, {axes.axisHypPairPhi, axes.axisPt, axes.axisHypPairPhi});
      }
    }

    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    // standards hardcoded in builder ...
    // ...but can be changed easily since fitter is public
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxDXYIni(4.0f);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
    // LUT has to be loaded later
    lut = nullptr;
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);

    // mag field has to be set later
    fitter.setBz(-999.9f); // will NOT make sense if not changed

    // set V0 parameters in the helper
    straHelper.v0selections.minCrossedRows = v0Selections.minTPCrows;
    straHelper.v0selections.dcanegtopv = std::min(v0Selections.dcapiontopv.value, v0Selections.dcaprotontopv.value);
    straHelper.v0selections.dcapostopv = std::min(v0Selections.dcapiontopv.value, v0Selections.dcaprotontopv.value);
    straHelper.v0selections.v0cospa = v0Selections.v0cospa;
    straHelper.v0selections.dcav0dau = v0Selections.dcav0dau;
    straHelper.v0selections.v0radius = v0Selections.v0radius;
    straHelper.v0selections.maxDaughterEta = v0Selections.daughterEtaCut;

    ccdb->setURL(ccdbConfigurations.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // inspect histogram sizes, please
    histos.print();
  }

  template <typename TCollision> // TCollision should be of the type: soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator or so
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber() || (ccdbConfigurations.useCustomRunNumber && mRunNumber == ccdbConfigurations.customRunNumber)) {
      return;
    }

    mRunNumber = ccdbConfigurations.useCustomRunNumber ? ccdbConfigurations.customRunNumber : collision.runNumber();

    if (doPairPropagationToDCA) {
      // In case override, don't proceed, please - no CCDB access required
      if (ccdbConfigurations.useCustomMagField) {
        magField = ccdbConfigurations.customMagField;
        o2::parameters::GRPMagField grpmag;
        if (fabs(magField) > 1e-5) {
          grpmag.setL3Current(30000.f / (magField / 5.0f));
        }
        o2::base::Propagator::initFieldFromGRP(&grpmag);
      } else {
        o2::parameters::GRPObject* grpo = ccdb->getForRun<o2::parameters::GRPObject>(ccdbConfigurations.grpPath, mRunNumber);
        o2::parameters::GRPMagField* grpmag = 0x0;
        if (grpo) {
          o2::base::Propagator::initFieldFromGRP(grpo);
          // Fetch magnetic field from ccdb for current collision
          magField = grpo->getNominalL3Field();
          LOG(info) << "Retrieved GRP for run " << mRunNumber << " with magnetic field of " << magField << " kZG";
        } else {
          grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, mRunNumber);
          if (!grpmag) {
            LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField and " << ccdbConfigurations.grpPath << " of object GRPObject for run " << mRunNumber;
          }
          o2::base::Propagator::initFieldFromGRP(grpmag);
          // Fetch magnetic field from ccdb for current collision
          magField = std::lround(5.f * grpmag->getL3Current() / 30000.f);
          LOG(info) << "Retrieved GRP for run " << mRunNumber << " with magnetic field of " << magField << " kZG";
        }
      }

      // load matLUT for this timestamp
      if (!lut) {
        LOG(info) << "Loading material look-up table for timestamp: " << mRunNumber;
        lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->template getForRun<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath.value, mRunNumber));
        straHelper.lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->template getForRun<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath.value, mRunNumber));
      } else {
        LOG(info) << "Material look-up table already in place. Not reloading.";
      }
      LOG(info) << "Setting global propagator material propagation LUT";
      o2::base::Propagator::Instance()->setMatLUT(lut);
      o2::base::Propagator::Instance()->setMatLUT(straHelper.lut);

      fitter.setBz(magField);
      straHelper.fitter.setBz(magField);
    }

    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, collision.runNumber(), collision.timestamp(), cfgSkimmedTrigger.value);
      zorro.populateHistRegistry(histos, collision.runNumber());
    }

    // machine learning initialization if requested
    if (mlConfigurations.calculateK0ShortScores ||
        mlConfigurations.calculateLambdaScores ||
        mlConfigurations.calculateAntiLambdaScores) {
      int64_t timeStampML = collision.timestamp();
      if (mlConfigurations.timestampCCDB.value != -1)
        timeStampML = mlConfigurations.timestampCCDB.value;
      loadMachines(timeStampML);
    }
  }

  // function to load models for ML-based classifiers
  void loadMachines(int64_t timeStampML)
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

  // Taken from https://github.com/AliceO2Group/O2Physics/blob/master/PWGLF/TableProducer/Strangeness/sigma0builder.cxx#L319
  // Thanks Gianni!
  // ______________________________________________________
  // Struct to store V0Pair properties
  struct PairTopoInfo {
    float X = -999.f;
    float Y = -999.f;
    float Z = -999.f;
    std::array<float, 3> hyperonMomentum = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> antiHyperonMomentum = {0.0f, 0.0f, 0.0f};
    float DCADau = -999.f;
    float CosPA = -1.f;
    float OpAngle = -999.f;
    float Eta() const
    {
      return RecoDecay::eta(std::array{hyperonMomentum[0] + antiHyperonMomentum[0], hyperonMomentum[1] + antiHyperonMomentum[1], hyperonMomentum[2] + antiHyperonMomentum[2]});
    }
    float Phi() const
    {
      return RecoDecay::phi(std::array{hyperonMomentum[0] + antiHyperonMomentum[0], hyperonMomentum[1] + antiHyperonMomentum[1]});
    }
  };

  template <typename TV0>
  PairTopoInfo propagateV0PairToDCA(float pvX, float pvY, float pvZ, TV0 const& v01, TV0 const& v02)
  {
    PairTopoInfo pairInfo;

    // Positions
    ROOT::Math::XYZVector v01position(v01.x(), v01.y(), v01.z());
    ROOT::Math::XYZVector v02position(v02.x(), v02.y(), v02.z());

    // Momenta
    ROOT::Math::XYZVector v01momentum(v01.px(), v01.py(), v01.pz());
    ROOT::Math::XYZVector v02momentum(v02.px(), v02.py(), v02.pz());

    // Momenta (normalized)
    ROOT::Math::XYZVector v01momentumNorm(v01.px() / v01.p(), v01.py() / v01.p(), v01.pz() / v01.p());
    ROOT::Math::XYZVector v02momentumNorm(v02.px() / v02.p(), v02.py() / v02.p(), v02.pz() / v02.p());

    // DCADau calculation (using full momenta for precision)
    ROOT::Math::XYZVector posdiff = v02position - v01position;
    ROOT::Math::XYZVector cross = v01momentum.Cross(v02momentum);

    float d = 1.0f - TMath::Power(v01momentumNorm.Dot(v02momentumNorm), 2);
    float t = posdiff.Dot(v01momentumNorm - v01momentumNorm.Dot(v02momentumNorm) * v02momentumNorm) / d;
    float s = -posdiff.Dot(v02momentumNorm - v01momentumNorm.Dot(v02momentumNorm) * v01momentumNorm) / d;

    ROOT::Math::XYZVector pointOn1 = v01position + t * v01momentumNorm;
    ROOT::Math::XYZVector pointOn2 = v02position + s * v02momentumNorm;
    ROOT::Math::XYZVector PCA = 0.5 * (pointOn1 + pointOn2);

    pairInfo.hyperonMomentum[0] = v01.px();
    pairInfo.hyperonMomentum[1] = v01.py();
    pairInfo.hyperonMomentum[2] = v01.pz();
    pairInfo.antiHyperonMomentum[0] = v02.px();
    pairInfo.antiHyperonMomentum[1] = v02.py();
    pairInfo.antiHyperonMomentum[2] = v02.pz();

    // Calculate properties and fill struct
    pairInfo.DCADau = (cross.Mag2() > 0) ? std::abs(posdiff.Dot(cross)) / cross.R() : 999.f;

    pairInfo.CosPA = RecoDecay::cpa(
      std::array{pvX, pvY, pvZ},
      std::array{PCA.X(), PCA.Y(), PCA.Z()},
      std::array{v01.px() + v02.px(),
                 v01.py() + v02.py(),
                 v01.pz() + v02.pz()});

    TVector3 hyp1Momentum(pairInfo.hyperonMomentum[0], pairInfo.hyperonMomentum[1], pairInfo.hyperonMomentum[2]);
    TVector3 hyp2Momentum(pairInfo.antiHyperonMomentum[0], pairInfo.antiHyperonMomentum[1], pairInfo.antiHyperonMomentum[2]);
    pairInfo.OpAngle = hyp1Momentum.Angle(hyp2Momentum);

    if (d < 1e-5f) {                                 // Parallel or nearly parallel lines
      pairInfo.X = pairInfo.Y = pairInfo.Z = -999.f; // should we use another dummy value? Perhaps 999.f?
      return pairInfo;
    }

    pairInfo.X = PCA.X();
    pairInfo.Y = PCA.Y();
    pairInfo.Z = PCA.Z();

    return pairInfo;
  }

  template <typename TCascade>
  PairTopoInfo propagateCascPairToDCA(float pvX, float pvY, float pvZ, TCascade const& casc1, TCascade const& casc2)
  {
    PairTopoInfo pairInfo;

    const std::array<float, 3> vtxCasc1 = {casc1.x(), casc1.y(), casc1.z()};
    const std::array<float, 3> vtxCasc2 = {casc2.x(), casc2.y(), casc2.z()};

    const std::array<float, 3> momCasc1 = {casc1.px(), casc1.py(), casc1.pz()};
    const std::array<float, 3> momCasc2 = {casc2.px(), casc2.py(), casc2.pz()};

    const std::array<float, 21> covCasc1 = {999.};
    const std::array<float, 21> covCasc2 = {999.};

    o2::track::TrackParCov cascTrack1(vtxCasc1, momCasc1, covCasc1, casc1.sign(), true);
    cascTrack1.setPID(o2::track::PID::XiMinus);
    // cascTrack1.setPID(o2::track::PID::OmegaMinus);
    o2::track::TrackParCov cascTrack2(vtxCasc2, momCasc2, covCasc2, casc2.sign(), true);
    cascTrack2.setPID(o2::track::PID::XiMinus);
    // cascTrack2.setPID(o2::track::PID::OmegaMinus);

    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(cascTrack1, cascTrack2);
    } catch (...) {
      return pairInfo;
    }
    if (nCand == 0) {
      return pairInfo;
    }

    fitter.propagateTracksToVertex(); // propagate e and K to D vertex
    if (!fitter.isPropagateTracksToVertexDone()) {
      return pairInfo;
    }

    auto lCasc1Track = fitter.getTrack(0);
    auto lCasc2Track = fitter.getTrack(1);

    lCasc1Track.getPxPyPzGlo(pairInfo.hyperonMomentum);
    lCasc2Track.getPxPyPzGlo(pairInfo.antiHyperonMomentum);

    // DCA between cascade daughters
    pairInfo.DCADau = std::sqrt(fitter.getChi2AtPCACandidate());

    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    pairInfo.X = vtx[0];
    pairInfo.Y = vtx[1];
    pairInfo.Z = vtx[2];

    pairInfo.CosPA = RecoDecay::cpa(
      std::array{pvX, pvY, pvZ},
      std::array{vtx[0], vtx[1], vtx[2]},
      std::array{pairInfo.hyperonMomentum[0] + pairInfo.antiHyperonMomentum[0],
                 pairInfo.hyperonMomentum[1] + pairInfo.antiHyperonMomentum[1],
                 pairInfo.hyperonMomentum[2] + pairInfo.antiHyperonMomentum[2]});

    // Momenta
    TVector3 casc1Momentum(pairInfo.hyperonMomentum[0], pairInfo.hyperonMomentum[1], pairInfo.hyperonMomentum[2]);
    TVector3 casc2Momentum(pairInfo.antiHyperonMomentum[0], pairInfo.antiHyperonMomentum[1], pairInfo.antiHyperonMomentum[2]);

    pairInfo.OpAngle = casc1Momentum.Angle(casc2Momentum);

    return pairInfo;
  }

  template <typename TCollision>
  bool isEventAccepted(TCollision collision, bool fillHists)
  // check whether the collision passes our collision selections
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

    if (requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);

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

    if (std::abs(collision.posZ()) > maxZVtxPosition) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);

    if (requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);

    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);

    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);

    if (requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);

    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);

    if (requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 10 microseconds */);

    if (requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 4 microseconds */);

    if (minOccupancy > 0 && collision.trackOccupancyInTimeRange() < minOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 13 /* Below min occupancy */);
    if (maxOccupancy > 0 && collision.trackOccupancyInTimeRange() > maxOccupancy) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 14 /* Above max occupancy */);

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
      selGapSide = sgSelector.trueGap(collision, upcCuts.fv0Cut, upcCuts.ft0aCut, upcCuts.ft0cCut, upcCuts.zdcCut);
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
      BITSET(bitMap, selRadius);
    if (v0.v0radius() < v0Selections.v0radiusMax)
      BITSET(bitMap, selRadiusMax);
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (std::fabs(v0.dcapostopv()) > v0Selections.dcaprotontopv &&
        std::fabs(v0.dcanegtopv()) > v0Selections.dcapiontopv) {
      BITSET(bitMap, selDCAPosToPV);
      BITSET(bitMap, selDCANegToPV);
    } else if (std::fabs(v0.dcapostopv()) > v0Selections.dcapiontopv &&
               std::fabs(v0.dcanegtopv()) > v0Selections.dcaprotontopv) {
      BITSET(bitMap, selDCAPosToPV);
      BITSET(bitMap, selDCANegToPV);
    }
    // V0 cosine of pointing angle
    if (v0.v0cosPA() > v0Selections.v0cospa)
      BITSET(bitMap, selCosPA);
    // DCA between v0 daughters
    if (v0.dcaV0daughters() < v0Selections.dcav0dau)
      BITSET(bitMap, selDCAV0Dau);
    // DCA V0 to prim vtx
    if (v0.dcav0topv() < v0Selections.dcav0topv)
      BITSET(bitMap, selDCAV0ToPV);

    //
    // rapidity
    //
    if (std::fabs(rapidityLambda) < v0Selections.rapidityCut)
      BITSET(bitMap, selLambdaRapidity);
    if (std::fabs(rapidityK0Short) < v0Selections.rapidityCut)
      BITSET(bitMap, selK0ShortRapidity);

    //
    // invariant mass window
    //
    if (std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0Selections.v0MassWindow)
      BITSET(bitMap, selK0ShortMassWindow);
    if (std::fabs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0Selections.v0MassWindow)
      BITSET(bitMap, selLambdaMassWindow);
    if (std::fabs(v0.mAntiLambda() - o2::constants::physics::MassLambda0Bar) < v0Selections.v0MassWindow)
      BITSET(bitMap, selAntiLambdaMassWindow);

    //
    // competing mass rejection
    //
    if (std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) > v0Selections.compMassRejection)
      BITSET(bitMap, selK0ShortMassRejection);
    if (std::fabs(v0.mLambda() - o2::constants::physics::MassLambda0) > v0Selections.compMassRejection)
      BITSET(bitMap, selLambdaMassRejection);

    auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

    //
    // ITS quality flags
    //
    if (posTrackExtra.itsNCls() >= v0Selections.minITSclusters)
      BITSET(bitMap, selPosGoodITSTrack);
    if (negTrackExtra.itsNCls() >= v0Selections.minITSclusters)
      BITSET(bitMap, selNegGoodITSTrack);

    //
    // TPC quality flags
    //
    if (posTrackExtra.tpcCrossedRows() >= v0Selections.minTPCrows)
      BITSET(bitMap, selPosGoodTPCTrack);
    if (negTrackExtra.tpcCrossedRows() >= v0Selections.minTPCrows)
      BITSET(bitMap, selNegGoodTPCTrack);

    //
    // TPC PID
    //
    if (std::fabs(posTrackExtra.tpcNSigmaPi()) < v0Selections.tpcPidNsigmaCut)
      BITSET(bitMap, selTPCPIDPositivePion);
    if (std::fabs(posTrackExtra.tpcNSigmaPr()) < v0Selections.tpcPidNsigmaCut)
      BITSET(bitMap, selTPCPIDPositiveProton);
    if (std::fabs(negTrackExtra.tpcNSigmaPi()) < v0Selections.tpcPidNsigmaCut)
      BITSET(bitMap, selTPCPIDNegativePion);
    if (std::fabs(negTrackExtra.tpcNSigmaPr()) < v0Selections.tpcPidNsigmaCut)
      BITSET(bitMap, selTPCPIDNegativeProton);

    //
    // TOF PID in DeltaT
    // Positive track
    if (!posTrackExtra.hasTOF() || std::fabs(v0.posTOFDeltaTLaPr()) < v0Selections.maxDeltaTimeProton)
      BITSET(bitMap, selTOFDeltaTPositiveProtonLambda);
    if (!posTrackExtra.hasTOF() || std::fabs(v0.posTOFDeltaTLaPi()) < v0Selections.maxDeltaTimePion)
      BITSET(bitMap, selTOFDeltaTPositivePionLambda);
    if (!posTrackExtra.hasTOF() || std::fabs(v0.posTOFDeltaTK0Pi()) < v0Selections.maxDeltaTimePion)
      BITSET(bitMap, selTOFDeltaTPositivePionK0Short);
    // Negative track
    if (!negTrackExtra.hasTOF() || std::fabs(v0.negTOFDeltaTLaPr()) < v0Selections.maxDeltaTimeProton)
      BITSET(bitMap, selTOFDeltaTNegativeProtonLambda);
    if (!negTrackExtra.hasTOF() || std::fabs(v0.negTOFDeltaTLaPi()) < v0Selections.maxDeltaTimePion)
      BITSET(bitMap, selTOFDeltaTNegativePionLambda);
    if (!negTrackExtra.hasTOF() || std::fabs(v0.negTOFDeltaTK0Pi()) < v0Selections.maxDeltaTimePion)
      BITSET(bitMap, selTOFDeltaTNegativePionK0Short);

    //
    // TOF PID in NSigma
    // Positive track
    if (!posTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaLaPr()) < v0Selections.tofPidNsigmaCutLaPr)
      BITSET(bitMap, selTOFNSigmaPositiveProtonLambda);
    if (!posTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaALaPi()) < v0Selections.tofPidNsigmaCutLaPi)
      BITSET(bitMap, selTOFNSigmaPositivePionLambda);
    if (!posTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaK0PiPlus()) < v0Selections.tofPidNsigmaCutK0Pi)
      BITSET(bitMap, selTOFNSigmaPositivePionK0Short);
    // Negative track
    if (!negTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaALaPr()) < v0Selections.tofPidNsigmaCutLaPr)
      BITSET(bitMap, selTOFNSigmaNegativeProtonLambda);
    if (!negTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaLaPi()) < v0Selections.tofPidNsigmaCutLaPi)
      BITSET(bitMap, selTOFNSigmaNegativePionLambda);
    if (!negTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaK0PiMinus()) < v0Selections.tofPidNsigmaCutK0Pi)
      BITSET(bitMap, selTOFNSigmaNegativePionK0Short);

    //
    // ITS only tag
    if (posTrackExtra.tpcCrossedRows() < 1)
      BITSET(bitMap, selPosItsOnly);
    if (negTrackExtra.tpcCrossedRows() < 1)
      BITSET(bitMap, selNegItsOnly);

    //
    // TPC only tag
    if (posTrackExtra.detectorMap() != o2::aod::track::TPC)
      BITSET(bitMap, selPosNotTPCOnly);
    if (negTrackExtra.detectorMap() != o2::aod::track::TPC)
      BITSET(bitMap, selNegNotTPCOnly);

    //
    // proper lifetime
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda"))
      BITSET(bitMap, selLambdaCTau);
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S"))
      BITSET(bitMap, selK0ShortCTau);

    //
    // armenteros
    if (v0.qtarm() * v0Selections.armPodCut > std::fabs(v0.alpha()) || v0Selections.armPodCut < 1e-4)
      BITSET(bitMap, selK0ShortArmenteros);

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
    // invariant mass window
    //
    if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cascSelections.v0MassWindow)
      return false;
    if (isXi && std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) > cascSelections.cascMassWindow)
      return false;
    if (!isXi && std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > cascSelections.cascMassWindow)
      return false;

    //
    // competing mass rejection
    //
    if (isXi && std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < cascSelections.compMassRejection)
      return false;
    if (!isXi && std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) < cascSelections.compMassRejection)
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
    if (isXi && std::fabs(bachTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
      return false;
    if (!isXi && std::fabs(bachTrackExtra.tpcNSigmaKa()) > cascSelections.tpcPidNsigmaCut)
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
    // TOF PID in DeltaT
    // Bachelor track
    if (bachTrackExtra.hasTOF()) {
      if (isXi && std::fabs(casc.bachTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
        return false;
      if (!isXi && std::fabs(casc.bachTOFDeltaTOmKa()) > cascSelections.maxDeltaTimeKaon)
        return false;
    }
    // Positive track
    if (posTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> positive daughter = proton
        if (isXi && std::fabs(casc.posTOFDeltaTXiPr()) > cascSelections.maxDeltaTimeProton)
          return false;
        if (!isXi && std::fabs(casc.posTOFDeltaTOmPr()) > cascSelections.maxDeltaTimeProton)
          return false;
      } else { // Xi+ or Omega+ --> positive daughter = pion
        if (isXi && std::fabs(casc.posTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
          return false;
        if (!isXi && std::fabs(casc.posTOFDeltaTOmPi()) > cascSelections.maxDeltaTimePion)
          return false;
      }
    }
    // Negative track
    if (negTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> negative daughter = pion
        if (isXi && std::fabs(casc.negTOFDeltaTXiPi()) > cascSelections.maxDeltaTimePion)
          return false;
        if (!isXi && std::fabs(casc.negTOFDeltaTOmPi()) > cascSelections.maxDeltaTimePion)
          return false;
      } else { // Xi+ or Omega+ --> negative daughter = proton
        if (isXi && std::fabs(casc.negTOFDeltaTXiPr()) > cascSelections.maxDeltaTimeProton)
          return false;
        if (!isXi && std::fabs(casc.negTOFDeltaTOmPr()) > cascSelections.maxDeltaTimeProton)
          return false;
      }
    }

    //
    // TOF PID in NSigma
    // Bachelor track
    if (bachTrackExtra.hasTOF()) {
      if (isXi && std::fabs(casc.tofNSigmaXiPi()) > cascSelections.tofPidNsigmaCutXiPi)
        return false;
      if (!isXi && std::fabs(casc.tofNSigmaOmKa()) > cascSelections.tofPidNsigmaCutOmKa)
        return false;
    }
    // Positive track
    if (posTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> positive daughter = proton
        if (isXi && std::fabs(casc.tofNSigmaXiLaPr()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPr()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
      } else { // Xi+ or Omega+ --> positive daughter = pion
        if (isXi && std::fabs(casc.tofNSigmaXiLaPi()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPi()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
      }
    }
    // Negative track
    if (negTrackExtra.hasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> negative daughter = pion
        if (isXi && std::fabs(casc.tofNSigmaXiLaPr()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPr()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
      } else { // Xi+ or Omega+ --> negative daughter = proton
        if (isXi && std::fabs(casc.tofNSigmaXiLaPi()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
        if (!isXi && std::fabs(casc.tofNSigmaOmLaPi()) > cascSelections.tofPidNsigmaCutLaPr)
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
            if (cascMC.pdgCode() != PDG_t::kXiMinus || cascMC.pdgCodePositive() != PDG_t::kProton || cascMC.pdgCodeNegative() != PDG_t::kPiMinus || cascMC.pdgCodeBachelor() != -PDG_t::kPiMinus)
              return false;
          } else {
            if (cascMC.pdgCode() != PDG_t::kXiPlusBar || cascMC.pdgCodePositive() != PDG_t::kPiPlus || cascMC.pdgCodeNegative() != PDG_t::kProtonBar || cascMC.pdgCodeBachelor() != PDG_t::kPiPlus)
              return false;
          }
        } else {
          if (casc.sign() < 0) {
            if (cascMC.pdgCode() != PDG_t::kOmegaMinus || cascMC.pdgCodePositive() != PDG_t::kProton || cascMC.pdgCodeNegative() != PDG_t::kPiMinus || cascMC.pdgCodeBachelor() != PDG_t::kKMinus)
              return false;
          } else {
            if (cascMC.pdgCode() != PDG_t::kOmegaPlusBar || cascMC.pdgCodePositive() != PDG_t::kPiPlus || cascMC.pdgCodeNegative() != PDG_t::kProtonBar || cascMC.pdgCodeBachelor() != PDG_t::kKPlus)
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

    if (v0.pdgCode() == PDG_t::kK0Short && v0.pdgCodePositive() == PDG_t::kPiPlus && v0.pdgCodeNegative() == PDG_t::kPiMinus) {
      BITSET(bitMap, selConsiderK0Short);
      if (v0.isPhysicalPrimary())
        BITSET(bitMap, selPhysPrimK0Short);
    }
    if (v0.pdgCode() == PDG_t::kLambda0 && v0.pdgCodePositive() == PDG_t::kProton && v0.pdgCodeNegative() == PDG_t::kPiMinus) {
      BITSET(bitMap, selConsiderLambda);
      if (v0.isPhysicalPrimary())
        BITSET(bitMap, selPhysPrimLambda);
    }
    if (v0.pdgCode() == PDG_t::kLambda0Bar && v0.pdgCodePositive() == PDG_t::kPiPlus && v0.pdgCodeNegative() == PDG_t::kProtonBar) {
      BITSET(bitMap, selConsiderAntiLambda);
      if (v0.isPhysicalPrimary())
        BITSET(bitMap, selPhysPrimAntiLambda);
    }
    return bitMap;
  }

  bool verifyMask(uint64_t bitmap, uint64_t mask)
  {
    return (bitmap & mask) == mask;
  }

  template <typename TV0>
  void analyseV0Candidate(TV0 v0, float pt, uint64_t selMap, std::vector<int>& selK0ShortIndices, std::vector<int>& selLambdaIndices, std::vector<int>& selAntiLambdaIndices /*, int v0TableOffset*/)
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
    if (passK0ShortSelections)
      selK0ShortIndices.push_back(v0.globalIndex());
    if (passLambdaSelections)
      selLambdaIndices.push_back(v0.globalIndex());
    if (passAntiLambdaSelections)
      selAntiLambdaIndices.push_back(v0.globalIndex());
  }

  template <typename TCollision, typename THyperon>
  void fillQAplot(TCollision collision, PairTopoInfo pair, THyperon hyperon, THyperon antiHyperon, float pt, float invmass, int type)
  { // fill QA information about hyperon - antihyperon pair
    if (type == 0) {
      if constexpr (requires { hyperon.mK0Short(); antiHyperon.mK0Short(); }) { // check if v0 information is available
        auto posTrackExtraHyperon = hyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraHyperon = hyperon.template negTrackExtra_as<DauTracks>();

        auto posTrackExtraAntiHyperon = antiHyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraAntiHyperon = antiHyperon.template negTrackExtra_as<DauTracks>();

        float hyperonDecayLength = std::sqrt(std::pow(hyperon.x() - collision.posX(), 2) + std::pow(hyperon.y() - collision.posY(), 2) + std::pow(hyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassKaonNeutral / (hyperon.p() + 1E-10);
        float antiHyperonDecayLength = std::sqrt(std::pow(antiHyperon.x() - collision.posX(), 2) + std::pow(antiHyperon.y() - collision.posY(), 2) + std::pow(antiHyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassKaonNeutral / (antiHyperon.p() + 1E-10);

        // Candidates after K0s selections
        histos.fill(HIST("K0sK0s/K0s/hPosDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("K0sK0s/K0s/hNegDCAToPV"), hyperon.dcanegtopv());
        histos.fill(HIST("K0sK0s/K0s/hDCAV0Daughters"), hyperon.dcaV0daughters());
        histos.fill(HIST("K0sK0s/K0s/hDCAV0ToPV"), hyperon.dcav0topv());
        histos.fill(HIST("K0sK0s/K0s/hV0PointingAngle"), hyperon.v0cosPA());
        histos.fill(HIST("K0sK0s/K0s/hV0Radius"), hyperon.v0radius());
        histos.fill(HIST("K0sK0s/K0s/hV0DecayLength"), hyperonDecayLength);
        histos.fill(HIST("K0sK0s/K0s/hV0InvMassWindow"), hyperon.mK0Short() - o2::constants::physics::MassK0Short);
        histos.fill(HIST("K0sK0s/K0s/h2dCompetingMassRej"), hyperon.mLambda(), hyperon.mK0Short());
        histos.fill(HIST("K0sK0s/K0s/h2dArmenteros"), hyperon.alpha(), hyperon.qtarm()); // cross-check
        histos.fill(HIST("K0sK0s/K0s/hPosTPCNsigma"), posTrackExtraHyperon.tpcNSigmaPi());
        histos.fill(HIST("K0sK0s/K0s/hNegTPCNsigma"), negTrackExtraHyperon.tpcNSigmaPi());
        histos.fill(HIST("K0sK0s/K0s/h2dPositiveITSvsTPCpts"), posTrackExtraHyperon.tpcCrossedRows(), posTrackExtraHyperon.itsNCls());
        histos.fill(HIST("K0sK0s/K0s/h2dNegativeITSvsTPCpts"), negTrackExtraHyperon.tpcCrossedRows(), negTrackExtraHyperon.itsNCls());
        // Candidates after K0s selections
        histos.fill(HIST("K0sK0s/K0s/hPosDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("K0sK0s/K0s/hNegDCAToPV"), antiHyperon.dcanegtopv());
        histos.fill(HIST("K0sK0s/K0s/hDCAV0Daughters"), antiHyperon.dcaV0daughters());
        histos.fill(HIST("K0sK0s/K0s/hDCAV0ToPV"), antiHyperon.dcav0topv());
        histos.fill(HIST("K0sK0s/K0s/hV0PointingAngle"), antiHyperon.v0cosPA());
        histos.fill(HIST("K0sK0s/K0s/hV0Radius"), antiHyperon.v0radius());
        histos.fill(HIST("K0sK0s/K0s/hV0DecayLength"), antiHyperonDecayLength);
        histos.fill(HIST("K0sK0s/K0s/hV0InvMassWindow"), antiHyperon.mK0Short() - o2::constants::physics::MassK0Short);
        histos.fill(HIST("K0sK0s/K0s/h2dCompetingMassRej"), antiHyperon.mLambda(), antiHyperon.mK0Short());
        histos.fill(HIST("K0sK0s/K0s/h2dArmenteros"), antiHyperon.alpha(), antiHyperon.qtarm()); // cross-check
        histos.fill(HIST("K0sK0s/K0s/hPosTPCNsigma"), posTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("K0sK0s/K0s/hNegTPCNsigma"), negTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("K0sK0s/K0s/h2dPositiveITSvsTPCpts"), posTrackExtraAntiHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("K0sK0s/K0s/h2dNegativeITSvsTPCpts"), negTrackExtraAntiHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.itsNCls());

        if (doQA) {
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPosDCAToPV"), std::min(hyperon.dcapostopv(), antiHyperon.dcapostopv()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsNegDCAToPV"), std::min(hyperon.dcanegtopv(), antiHyperon.dcanegtopv()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsDCAV0Daughters"), std::max(hyperon.dcaV0daughters(), antiHyperon.dcaV0daughters()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsDCAV0ToPV"), std::max(hyperon.dcav0topv(), antiHyperon.dcav0topv()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsV0PointingAngle"), std::min(hyperon.v0cosPA(), antiHyperon.v0cosPA()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsV0Radius"), std::min(hyperon.v0radius(), antiHyperon.v0radius()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsDecayLength"), std::max(hyperonDecayLength, antiHyperonDecayLength), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsInvMassWindow"), std::max(std::abs(hyperon.mK0Short() - o2::constants::physics::MassK0Short), std::abs(antiHyperon.mK0Short() - o2::constants::physics::MassK0Short)), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPosTPCNsigma"), std::max(posTrackExtraHyperon.tpcNSigmaPi(), posTrackExtraAntiHyperon.tpcNSigmaPi()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsNegTPCNsigma"), std::max(negTrackExtraHyperon.tpcNSigmaPi(), negTrackExtraAntiHyperon.tpcNSigmaPi()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPosITSclusters"), std::min(posTrackExtraHyperon.itsNCls(), posTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsNegITSclusters"), std::min(negTrackExtraHyperon.itsNCls(), negTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPosNbrCrossedRows"), std::min(posTrackExtraHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsNegNbrCrossedRows"), std::min(negTrackExtraHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPairRadius3D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y, pair.Z), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPairRadius2D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPairZ"), pair.Z, pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsDCAPair"), pair.DCADau, pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPairCosPA"), pair.CosPA, pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPairOpAngle"), pair.OpAngle, pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPairEta"), pair.Eta(), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dMassK0sK0sVsPairPhi"), pair.Phi(), pt, invmass);
          histos.fill(HIST("QA/K0sK0s/h3dDeltaEtaK0sK0sVsPairEta"), pair.Eta(), pt, hyperon.eta() - antiHyperon.eta());
          histos.fill(HIST("QA/K0sK0s/h3dDeltaPhiK0sK0sVsPairPhi"), pair.Phi(), pt, std::abs(hyperon.phi() - antiHyperon.phi()));
        }
      }
    }
    if (type == 1) {
      if constexpr (requires { hyperon.mK0Short(); antiHyperon.mK0Short(); }) { // check if v0 information is available
        auto posTrackExtraHyperon = hyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraHyperon = hyperon.template negTrackExtra_as<DauTracks>();

        auto posTrackExtraAntiHyperon = antiHyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraAntiHyperon = antiHyperon.template negTrackExtra_as<DauTracks>();

        float hyperonDecayLength = std::sqrt(std::pow(hyperon.x() - collision.posX(), 2) + std::pow(hyperon.y() - collision.posY(), 2) + std::pow(hyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassLambda0 / (hyperon.p() + 1E-10);
        float antiHyperonDecayLength = std::sqrt(std::pow(antiHyperon.x() - collision.posX(), 2) + std::pow(antiHyperon.y() - collision.posY(), 2) + std::pow(antiHyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassLambda0 / (antiHyperon.p() + 1E-10);

        // Candidates after Lambda selections
        histos.fill(HIST("LaLaBar/Lambda/hPosDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("LaLaBar/Lambda/hNegDCAToPV"), hyperon.dcanegtopv());
        histos.fill(HIST("LaLaBar/Lambda/hDCAV0Daughters"), hyperon.dcaV0daughters());
        histos.fill(HIST("LaLaBar/Lambda/hDCAV0ToPV"), hyperon.dcav0topv());
        histos.fill(HIST("LaLaBar/Lambda/hV0PointingAngle"), hyperon.v0cosPA());
        histos.fill(HIST("LaLaBar/Lambda/hV0Radius"), hyperon.v0radius());
        histos.fill(HIST("LaLaBar/Lambda/hV0DecayLength"), hyperonDecayLength);
        histos.fill(HIST("LaLaBar/Lambda/hV0InvMassWindow"), hyperon.mLambda() - o2::constants::physics::MassLambda0);
        histos.fill(HIST("LaLaBar/Lambda/h2dCompetingMassRej"), hyperon.mLambda(), hyperon.mK0Short());
        histos.fill(HIST("LaLaBar/Lambda/hPosTPCNsigma"), posTrackExtraHyperon.tpcNSigmaPr());
        histos.fill(HIST("LaLaBar/Lambda/hNegTPCNsigma"), negTrackExtraHyperon.tpcNSigmaPi());
        histos.fill(HIST("LaLaBar/Lambda/h2dPositiveITSvsTPCpts"), posTrackExtraHyperon.tpcCrossedRows(), posTrackExtraHyperon.itsNCls());
        histos.fill(HIST("LaLaBar/Lambda/h2dNegativeITSvsTPCpts"), negTrackExtraHyperon.tpcCrossedRows(), negTrackExtraHyperon.itsNCls());
        // Candidates after AntiLambda selections
        histos.fill(HIST("LaLaBar/AntiLambda/hPosDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("LaLaBar/AntiLambda/hNegDCAToPV"), antiHyperon.dcapostopv());
        histos.fill(HIST("LaLaBar/AntiLambda/hDCAV0Daughters"), antiHyperon.dcaV0daughters());
        histos.fill(HIST("LaLaBar/AntiLambda/hDCAV0ToPV"), antiHyperon.dcav0topv());
        histos.fill(HIST("LaLaBar/AntiLambda/hV0PointingAngle"), antiHyperon.v0cosPA());
        histos.fill(HIST("LaLaBar/AntiLambda/hV0Radius"), antiHyperon.v0radius());
        histos.fill(HIST("LaLaBar/AntiLambda/hV0DecayLength"), antiHyperonDecayLength);
        histos.fill(HIST("LaLaBar/AntiLambda/hV0InvMassWindow"), antiHyperon.mAntiLambda() - o2::constants::physics::MassLambda0);
        histos.fill(HIST("LaLaBar/AntiLambda/h2dCompetingMassRej"), antiHyperon.mAntiLambda(), antiHyperon.mK0Short());
        histos.fill(HIST("LaLaBar/AntiLambda/hPosTPCNsigma"), posTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("LaLaBar/AntiLambda/hNegTPCNsigma"), negTrackExtraAntiHyperon.tpcNSigmaPr());
        histos.fill(HIST("LaLaBar/AntiLambda/h2dPositiveITSvsTPCpts"), posTrackExtraAntiHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("LaLaBar/AntiLambda/h2dNegativeITSvsTPCpts"), negTrackExtraAntiHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.itsNCls());

        if (doQA) {
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPosDCAToPV"), std::min(hyperon.dcapostopv(), antiHyperon.dcanegtopv()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsNegDCAToPV"), std::min(hyperon.dcanegtopv(), antiHyperon.dcapostopv()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsDCAV0Daughters"), std::max(hyperon.dcaV0daughters(), antiHyperon.dcaV0daughters()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsDCAV0ToPV"), std::max(hyperon.dcav0topv(), antiHyperon.dcav0topv()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsV0PointingAngle"), std::min(hyperon.v0cosPA(), antiHyperon.v0cosPA()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsV0Radius"), std::min(hyperon.v0radius(), antiHyperon.v0radius()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsDecayLength"), std::max(hyperonDecayLength, antiHyperonDecayLength), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsInvMassWindow"), std::max(std::abs(hyperon.mLambda() - o2::constants::physics::MassLambda0), std::abs(antiHyperon.mAntiLambda() - o2::constants::physics::MassLambda0)), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPosTPCNsigma"), std::max(posTrackExtraHyperon.tpcNSigmaPr(), negTrackExtraAntiHyperon.tpcNSigmaPr()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsNegTPCNsigma"), std::max(negTrackExtraHyperon.tpcNSigmaPi(), posTrackExtraAntiHyperon.tpcNSigmaPi()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPosITSclusters"), std::min(posTrackExtraHyperon.itsNCls(), negTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsNegITSclusters"), std::min(negTrackExtraHyperon.itsNCls(), posTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPosNbrCrossedRows"), std::min(posTrackExtraHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsNegNbrCrossedRows"), std::min(negTrackExtraHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPairRadius3D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y, pair.Z), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPairRadius2D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPairZ"), pair.Z, pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsDCAPair"), pair.DCADau, pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPairCosPA"), pair.CosPA, pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPairOpAngle"), pair.OpAngle, pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPairEta"), pair.Eta(), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dMassLaLaBarVsPairPhi"), pair.Phi(), pt, invmass);
          histos.fill(HIST("QA/LaLaBar/h3dDeltaEtaLaLaBarVsPairEta"), pair.Eta(), pt, hyperon.eta() - antiHyperon.eta());
          histos.fill(HIST("QA/LaLaBar/h3dDeltaPhiLaLaBarVsPairPhi"), pair.Phi(), pt, std::abs(hyperon.phi() - antiHyperon.phi()));
        }
      }
    }
    if (type == 2) {
      if constexpr (requires { hyperon.dcabachtopv(); antiHyperon.dcabachtopv(); }) { // check if Cascade information is available
        auto bachTrackExtraHyperon = hyperon.template bachTrackExtra_as<DauTracks>();
        auto posTrackExtraHyperon = hyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraHyperon = hyperon.template negTrackExtra_as<DauTracks>();

        auto bachTrackExtraAntiHyperon = antiHyperon.template bachTrackExtra_as<DauTracks>();
        auto posTrackExtraAntiHyperon = antiHyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraAntiHyperon = antiHyperon.template negTrackExtra_as<DauTracks>();

        float hyperonDecayLength = std::sqrt(std::pow(hyperon.x() - collision.posX(), 2) + std::pow(hyperon.y() - collision.posY(), 2) + std::pow(hyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassXiMinus / (hyperon.p() + 1E-10);
        float antiHyperonDecayLength = std::sqrt(std::pow(antiHyperon.x() - collision.posX(), 2) + std::pow(antiHyperon.y() - collision.posY(), 2) + std::pow(antiHyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassXiMinus / (antiHyperon.p() + 1E-10);

        // Candidates after Xi selections
        histos.fill(HIST("XiXiBar/Xi/hBachDCAToPV"), hyperon.dcabachtopv());
        histos.fill(HIST("XiXiBar/Xi/hPosDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("XiXiBar/Xi/hNegDCAToPV"), hyperon.dcanegtopv());
        histos.fill(HIST("XiXiBar/Xi/hDCACascDaughters"), hyperon.dcacascdaughters());
        histos.fill(HIST("XiXiBar/Xi/hDCAV0Daughters"), hyperon.dcaV0daughters());
        histos.fill(HIST("XiXiBar/Xi/hDCAV0ToPV"), hyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/Xi/hV0PointingAngle"), hyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/Xi/hV0Radius"), hyperon.v0radius());
        histos.fill(HIST("XiXiBar/Xi/hCascPointingAngle"), hyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/Xi/hCascRadius"), hyperon.cascradius());
        histos.fill(HIST("XiXiBar/Xi/hCascDecayLength"), hyperonDecayLength);
        histos.fill(HIST("XiXiBar/Xi/hV0InvMassWindow"), hyperon.mLambda() - o2::constants::physics::MassLambda0);
        histos.fill(HIST("XiXiBar/Xi/hCascInvMassWindow"), hyperon.mXi() - o2::constants::physics::MassXiMinus);
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
        histos.fill(HIST("XiXiBar/AntiXi/hNegDCAToPV"), antiHyperon.dcanegtopv());
        histos.fill(HIST("XiXiBar/AntiXi/hDCACascDaughters"), antiHyperon.dcacascdaughters());
        histos.fill(HIST("XiXiBar/AntiXi/hDCAV0Daughters"), antiHyperon.dcaV0daughters());
        histos.fill(HIST("XiXiBar/AntiXi/hDCAV0ToPV"), antiHyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/AntiXi/hV0PointingAngle"), antiHyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/AntiXi/hV0Radius"), antiHyperon.v0radius());
        histos.fill(HIST("XiXiBar/AntiXi/hCascPointingAngle"), antiHyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("XiXiBar/AntiXi/hCascRadius"), antiHyperon.cascradius());
        histos.fill(HIST("XiXiBar/AntiXi/hCascDecayLength"), antiHyperonDecayLength);
        histos.fill(HIST("XiXiBar/AntiXi/hV0InvMassWindow"), antiHyperon.mLambda() - o2::constants::physics::MassLambda0);
        histos.fill(HIST("XiXiBar/AntiXi/hCascInvMassWindow"), antiHyperon.mXi() - o2::constants::physics::MassXiMinus);
        histos.fill(HIST("XiXiBar/AntiXi/h2dCompetingMassRej"), antiHyperon.mXi(), antiHyperon.mOmega());
        histos.fill(HIST("XiXiBar/AntiXi/hBachTPCNsigma"), bachTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("XiXiBar/AntiXi/hPosTPCNsigma"), posTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("XiXiBar/AntiXi/hNegTPCNsigma"), negTrackExtraAntiHyperon.tpcNSigmaPr());
        histos.fill(HIST("XiXiBar/AntiXi/h2dBachelorITSvsTPCpts"), bachTrackExtraAntiHyperon.tpcCrossedRows(), bachTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("XiXiBar/AntiXi/h2dPositiveITSvsTPCpts"), posTrackExtraAntiHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("XiXiBar/AntiXi/h2dNegativeITSvsTPCpts"), negTrackExtraAntiHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.itsNCls());

        if (doQA) {
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsBachDCAToPV"), std::min(hyperon.dcabachtopv(), antiHyperon.dcabachtopv()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPosDCAToPV"), std::min(hyperon.dcapostopv(), antiHyperon.dcanegtopv()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsNegDCAToPV"), std::min(hyperon.dcanegtopv(), antiHyperon.dcapostopv()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsDCACascDaughters"), std::max(hyperon.dcacascdaughters(), antiHyperon.dcacascdaughters()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsDCAV0Daughters"), std::max(hyperon.dcaV0daughters(), antiHyperon.dcaV0daughters()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsDCAV0ToPV"), std::max(hyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), antiHyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ())), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsV0PointingAngle"), std::min(hyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), antiHyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ())), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsV0Radius"), std::min(hyperon.v0radius(), antiHyperon.v0radius()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsCascPointingAngle"), std::min(hyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()), antiHyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ())), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsCascRadius"), std::min(hyperon.cascradius(), antiHyperon.cascradius()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsDecayLength"), std::max(hyperonDecayLength, antiHyperonDecayLength), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsV0InvMassWindow"), std::max(std::abs(hyperon.mLambda() - o2::constants::physics::MassLambda0), std::abs(antiHyperon.mLambda() - o2::constants::physics::MassLambda0)), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsCascInvMassWindow"), std::max(std::abs(hyperon.mXi() - o2::constants::physics::MassXiMinus), std::abs(antiHyperon.mXi() - o2::constants::physics::MassXiMinus)), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsBachTPCNsigma"), std::max(bachTrackExtraHyperon.tpcNSigmaPi(), bachTrackExtraAntiHyperon.tpcNSigmaPi()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPosTPCNsigma"), std::max(posTrackExtraHyperon.tpcNSigmaPr(), negTrackExtraAntiHyperon.tpcNSigmaPr()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsNegTPCNsigma"), std::max(negTrackExtraHyperon.tpcNSigmaPi(), posTrackExtraAntiHyperon.tpcNSigmaPi()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsBachITSclusters"), std::min(bachTrackExtraHyperon.itsNCls(), bachTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPosITSclusters"), std::min(posTrackExtraHyperon.itsNCls(), negTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsNegITSclusters"), std::min(negTrackExtraHyperon.itsNCls(), posTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsBachNbrCrossedRows"), std::min(bachTrackExtraHyperon.tpcCrossedRows(), bachTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPosNbrCrossedRows"), std::min(posTrackExtraHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsNegNbrCrossedRows"), std::min(negTrackExtraHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPairRadius3D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y, pair.Z), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPairRadius2D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPairZ"), pair.Z, pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsDCAPair"), pair.DCADau, pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPairCosPA"), pair.CosPA, pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPairOpAngle"), pair.OpAngle, pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPairEta"), pair.Eta(), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dMassXiXiBarVsPairPhi"), pair.Phi(), pt, invmass);
          histos.fill(HIST("QA/XiXiBar/h3dDeltaEtaXiXiBarVsPairEta"), pair.Eta(), pt, hyperon.eta() - antiHyperon.eta());
          histos.fill(HIST("QA/XiXiBar/h3dDeltaPhiXiXiBarVsPairPhi"), pair.Phi(), pt, std::abs(hyperon.phi() - antiHyperon.phi()));
        }
      }
    }
    if (type == 3) {
      if constexpr (requires { hyperon.dcabachtopv(); antiHyperon.dcabachtopv(); }) { // check if Cascade information is available
        auto bachTrackExtraHyperon = hyperon.template bachTrackExtra_as<DauTracks>();
        auto posTrackExtraHyperon = hyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraHyperon = hyperon.template negTrackExtra_as<DauTracks>();

        auto bachTrackExtraAntiHyperon = antiHyperon.template bachTrackExtra_as<DauTracks>();
        auto posTrackExtraAntiHyperon = antiHyperon.template posTrackExtra_as<DauTracks>();
        auto negTrackExtraAntiHyperon = antiHyperon.template negTrackExtra_as<DauTracks>();

        float hyperonDecayLength = std::sqrt(std::pow(hyperon.x() - collision.posX(), 2) + std::pow(hyperon.y() - collision.posY(), 2) + std::pow(hyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassOmegaMinus / (hyperon.p() + 1E-10);
        float antiHyperonDecayLength = std::sqrt(std::pow(antiHyperon.x() - collision.posX(), 2) + std::pow(antiHyperon.y() - collision.posY(), 2) + std::pow(antiHyperon.z() - collision.posZ(), 2)) * o2::constants::physics::MassOmegaMinus / (antiHyperon.p() + 1E-10);

        // Candidates after Omega selections
        histos.fill(HIST("OmOmBar/Omega/hBachDCAToPV"), hyperon.dcabachtopv());
        histos.fill(HIST("OmOmBar/Omega/hPosDCAToPV"), hyperon.dcapostopv());
        histos.fill(HIST("OmOmBar/Omega/hNegDCAToPV"), hyperon.dcanegtopv());
        histos.fill(HIST("OmOmBar/Omega/hDCACascDaughters"), hyperon.dcacascdaughters());
        histos.fill(HIST("OmOmBar/Omega/hDCAV0Daughters"), hyperon.dcaV0daughters());
        histos.fill(HIST("OmOmBar/Omega/hDCAV0ToPV"), hyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/Omega/hV0PointingAngle"), hyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/Omega/hV0Radius"), hyperon.v0radius());
        histos.fill(HIST("OmOmBar/Omega/hCascPointingAngle"), hyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/Omega/hCascRadius"), hyperon.cascradius());
        histos.fill(HIST("OmOmBar/Omega/hCascDecayLength"), hyperonDecayLength);
        histos.fill(HIST("OmOmBar/Omega/hV0InvMassWindow"), hyperon.mLambda() - o2::constants::physics::MassLambda0);
        histos.fill(HIST("OmOmBar/Omega/hCascInvMassWindow"), hyperon.mOmega() - o2::constants::physics::MassOmegaMinus);
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
        histos.fill(HIST("OmOmBar/AntiOmega/hNegDCAToPV"), antiHyperon.dcanegtopv());
        histos.fill(HIST("OmOmBar/AntiOmega/hDCACascDaughters"), antiHyperon.dcacascdaughters());
        histos.fill(HIST("OmOmBar/AntiOmega/hDCAV0Daughters"), antiHyperon.dcaV0daughters());
        histos.fill(HIST("OmOmBar/AntiOmega/hDCAV0ToPV"), antiHyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/AntiOmega/hV0PointingAngle"), antiHyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/AntiOmega/hV0Radius"), antiHyperon.v0radius());
        histos.fill(HIST("OmOmBar/AntiOmega/hCascPointingAngle"), antiHyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        histos.fill(HIST("OmOmBar/AntiOmega/hCascRadius"), antiHyperon.cascradius());
        histos.fill(HIST("OmOmBar/AntiOmega/hCascDecayLength"), antiHyperonDecayLength);
        histos.fill(HIST("OmOmBar/AntiOmega/hV0InvMassWindow"), antiHyperon.mLambda() - o2::constants::physics::MassLambda0);
        histos.fill(HIST("OmOmBar/AntiOmega/hCascInvMassWindow"), antiHyperon.mOmega() - o2::constants::physics::MassOmegaMinus);
        histos.fill(HIST("OmOmBar/AntiOmega/h2dCompetingMassRej"), antiHyperon.mXi(), antiHyperon.mOmega());
        histos.fill(HIST("OmOmBar/AntiOmega/hBachTPCNsigma"), bachTrackExtraAntiHyperon.tpcNSigmaKa());
        histos.fill(HIST("OmOmBar/AntiOmega/hPosTPCNsigma"), posTrackExtraAntiHyperon.tpcNSigmaPi());
        histos.fill(HIST("OmOmBar/AntiOmega/hNegTPCNsigma"), negTrackExtraAntiHyperon.tpcNSigmaPr());
        histos.fill(HIST("OmOmBar/AntiOmega/h2dBachelorITSvsTPCpts"), bachTrackExtraAntiHyperon.tpcCrossedRows(), bachTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("OmOmBar/AntiOmega/h2dPositiveITSvsTPCpts"), posTrackExtraAntiHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.itsNCls());
        histos.fill(HIST("OmOmBar/AntiOmega/h2dNegativeITSvsTPCpts"), negTrackExtraAntiHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.itsNCls());

        if (doQA) {
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsBachDCAToPV"), std::min(hyperon.dcabachtopv(), antiHyperon.dcabachtopv()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPosDCAToPV"), std::min(hyperon.dcapostopv(), antiHyperon.dcanegtopv()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsNegDCAToPV"), std::min(hyperon.dcanegtopv(), antiHyperon.dcapostopv()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsDCACascDaughters"), std::max(hyperon.dcacascdaughters(), antiHyperon.dcacascdaughters()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsDCAV0Daughters"), std::max(hyperon.dcaV0daughters(), antiHyperon.dcaV0daughters()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsDCAV0ToPV"), std::max(hyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), antiHyperon.dcav0topv(collision.posX(), collision.posY(), collision.posZ())), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsV0PointingAngle"), std::min(hyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), antiHyperon.v0cosPA(collision.posX(), collision.posY(), collision.posZ())), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsV0Radius"), std::min(hyperon.v0radius(), antiHyperon.v0radius()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsCascPointingAngle"), std::min(hyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ()), antiHyperon.casccosPA(collision.posX(), collision.posY(), collision.posZ())), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsCascRadius"), std::min(hyperon.cascradius(), antiHyperon.cascradius()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsDecayLength"), std::max(hyperonDecayLength, antiHyperonDecayLength), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsV0InvMassWindow"), std::max(std::abs(hyperon.mLambda() - o2::constants::physics::MassLambda0), std::abs(antiHyperon.mLambda() - o2::constants::physics::MassLambda0)), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsCascInvMassWindow"), std::max(std::abs(hyperon.mOmega() - o2::constants::physics::MassOmegaMinus), std::abs(antiHyperon.mOmega() - o2::constants::physics::MassOmegaMinus)), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsBachTPCNsigma"), std::max(bachTrackExtraHyperon.tpcNSigmaPi(), bachTrackExtraAntiHyperon.tpcNSigmaPi()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPosTPCNsigma"), std::max(posTrackExtraHyperon.tpcNSigmaPr(), negTrackExtraAntiHyperon.tpcNSigmaPr()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsNegTPCNsigma"), std::max(negTrackExtraHyperon.tpcNSigmaPi(), posTrackExtraAntiHyperon.tpcNSigmaPi()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsBachITSclusters"), std::min(bachTrackExtraHyperon.itsNCls(), bachTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPosITSclusters"), std::min(posTrackExtraHyperon.itsNCls(), negTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsNegITSclusters"), std::min(negTrackExtraHyperon.itsNCls(), posTrackExtraAntiHyperon.itsNCls()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsBachNbrCrossedRows"), std::min(bachTrackExtraHyperon.tpcCrossedRows(), bachTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPosNbrCrossedRows"), std::min(posTrackExtraHyperon.tpcCrossedRows(), negTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsNegNbrCrossedRows"), std::min(negTrackExtraHyperon.tpcCrossedRows(), posTrackExtraAntiHyperon.tpcCrossedRows()), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPairRadius3D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y, pair.Z), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPairRadius2D"), RecoDecay::sqrtSumOfSquares(pair.X, pair.Y), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPairZ"), pair.Z, pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsDCAPair"), pair.DCADau, pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPairCosPA"), pair.CosPA, pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPairOpAngle"), pair.OpAngle, pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPairEta"), pair.Eta(), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dMassOmOmBarVsPairPhi"), pair.Phi(), pt, invmass);
          histos.fill(HIST("QA/OmOmBar/h3dDeltaEtaOmOmBarVsPairEta"), pair.Eta(), pt, hyperon.eta() - antiHyperon.eta());
          histos.fill(HIST("QA/OmOmBar/h3dDeltaPhiOmOmBarVsPairPhi"), pair.Phi(), pt, std::abs(hyperon.phi() - antiHyperon.phi()));
        }
      }
    }
  }

  template <typename TCollision, typename THyperon>
  void analyseHyperonPairCandidate(TCollision collision, THyperon hyperon, THyperon antiHyperon, float centrality, uint8_t gapSide, int type)
  // fill information related to the quarkonium mother
  // type = 0 (Lambda), 1 (Xi), 2 (Omega)
  {
    PairTopoInfo pair;
    if (doPairPropagationToDCA) {
      if (type == 0 || type == 1) {
        if constexpr (requires { hyperon.mK0Short(); antiHyperon.mK0Short(); }) {
          pair = propagateV0PairToDCA(collision.posX(), collision.posY(), collision.posZ(), hyperon, antiHyperon);
        }
      } else {
        if constexpr (requires { hyperon.dcabachtopv(); antiHyperon.dcabachtopv(); }) { // check if Cascade information is available
          pair = propagateCascPairToDCA(collision.posX(), collision.posY(), collision.posZ(), hyperon, antiHyperon);
        }
      }
    } else {
      pair.hyperonMomentum[0] = hyperon.px();
      pair.hyperonMomentum[1] = hyperon.py();
      pair.hyperonMomentum[2] = hyperon.pz();
      pair.antiHyperonMomentum[0] = antiHyperon.px();
      pair.antiHyperonMomentum[1] = antiHyperon.py();
      pair.antiHyperonMomentum[2] = antiHyperon.pz();
    }

    float pt = RecoDecay::pt(pair.hyperonMomentum[0] + pair.antiHyperonMomentum[0], pair.hyperonMomentum[1] + pair.antiHyperonMomentum[1]);

    float invmass = -1;
    if (type == 0)
      invmass = RecoDecay::m(std::array{pair.hyperonMomentum, pair.antiHyperonMomentum}, std::array{o2::constants::physics::MassKaonNeutral, o2::constants::physics::MassKaonNeutral});
    if (type == 1)
      invmass = RecoDecay::m(std::array{pair.hyperonMomentum, pair.antiHyperonMomentum}, std::array{o2::constants::physics::MassLambda0, o2::constants::physics::MassLambda0Bar});
    if (type == 2)
      invmass = RecoDecay::m(std::array{pair.hyperonMomentum, pair.antiHyperonMomentum}, std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassXiPlusBar});
    if (type == 3)
      invmass = RecoDecay::m(std::array{pair.hyperonMomentum, pair.antiHyperonMomentum}, std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassOmegaPlusBar});

    float rapidity = RecoDecay::y(std::array{pair.hyperonMomentum[0] + pair.antiHyperonMomentum[0], pair.hyperonMomentum[1] + pair.antiHyperonMomentum[1], pair.hyperonMomentum[2] + pair.antiHyperonMomentum[2]}, invmass);

    // rapidity cut on the quarkonium mother
    if (!doMCAssociation && std::fabs(rapidity) > rapidityCut)
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

          if (std::fabs(rapiditymc) > rapidityCut)
            return;

          if (hyperonMC.pdgCodeMother() == 441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // EtaC(1S)
            histos.fill(HIST("K0sK0s/h3dInvMassTrueEtaC1S"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // J/psi
            histos.fill(HIST("K0sK0s/h3dInvMassTrueJPsi"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 10441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC0
            histos.fill(HIST("K0sK0s/h3dInvMassTrueChiC0"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 20443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC1
            histos.fill(HIST("K0sK0s/h3dInvMassTrueChiC1"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 10443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // hC
            histos.fill(HIST("K0sK0s/h3dInvMassTrueHC"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 445 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // ChiC2
            histos.fill(HIST("K0sK0s/h3dInvMassTrueChiC2"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 100441 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // EtaC(2S)
            histos.fill(HIST("K0sK0s/h3dInvMassTrueEtaC2S"), centrality, ptmc, invmass);
          }
          if (hyperonMC.pdgCodeMother() == 100443 && hyperonMC.pdgCodeMother() == antiHyperonMC.pdgCodeMother()) { // Psi(2S)
            histos.fill(HIST("K0sK0s/h3dInvMassTruePsi2S"), centrality, ptmc, invmass);
          }
        }
      }

      histos.fill(HIST("K0sK0s/h3dMassK0sK0s"), centrality, pt, invmass);
      if (!isPP) { // in case of PbPb data
        if (gapSide == 0)
          histos.fill(HIST("K0sK0s/h3dMassK0sK0sSGA"), centrality, pt, invmass);
        else if (gapSide == 1)
          histos.fill(HIST("K0sK0s/h3dMassK0sK0sSGC"), centrality, pt, invmass);
        else if (gapSide == 2)
          histos.fill(HIST("K0sK0s/h3dMassK0sK0sDG"), centrality, pt, invmass);
        else
          histos.fill(HIST("K0sK0s/h3dMassK0sK0sHadronic"), centrality, pt, invmass);
      }
      fillQAplot(collision, pair, hyperon, antiHyperon, pt, invmass, type);
    }
    if (type == 1) {
      if (doMCAssociation) {
        if constexpr (requires { hyperon.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>(); }) { // check if MC information is available
          auto hyperonMC = hyperon.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
          auto antiHyperonMC = antiHyperon.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

          if (hyperonMC.pdgCodeMother() != antiHyperonMC.pdgCodeMother()) {
            return;
          }

          float ptmc = RecoDecay::pt(hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC());
          float rapiditymc = RecoDecay::y(std::array{hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC(), hyperonMC.pzMC() + antiHyperonMC.pzMC()}, pdgDB->Mass(hyperonMC.pdgCodeMother()));

          if (std::fabs(rapiditymc) > rapidityCut)
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
      fillQAplot(collision, pair, hyperon, antiHyperon, pt, invmass, type);
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

          if (std::fabs(rapiditymc) > rapidityCut)
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
      fillQAplot(collision, pair, hyperon, antiHyperon, pt, invmass, type);
    }
    if (type == 3) {
      if (doMCAssociation) {
        if constexpr (requires { hyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>(); }) { // check if MC information is available
          auto hyperonMC = hyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();
          auto antiHyperonMC = antiHyperon.template cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

          if (hyperonMC.pdgCodeMother() != antiHyperonMC.pdgCodeMother()) {
            return;
          }

          float ptmc = RecoDecay::pt(hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC());
          float rapiditymc = RecoDecay::y(std::array{hyperonMC.pxMC() + antiHyperonMC.pxMC(), hyperonMC.pyMC() + antiHyperonMC.pyMC(), hyperonMC.pzMC() + antiHyperonMC.pzMC()}, pdgDB->Mass(hyperonMC.pdgCodeMother()));

          if (std::fabs(rapiditymc) > rapidityCut)
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
      fillQAplot(collision, pair, hyperon, antiHyperon, pt, invmass, type);
    }
  }

  // function to check that the hyperon and antihyperon have different daughter tracks
  template <typename THyperon>
  bool checkTrackIndices(THyperon hyperon, THyperon antiHyperon)
  {
    if constexpr (requires { hyperon.template bachTrackExtra_as<DauTracks>(); }) { // cascade case: check if bachelor information is available
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
  void buildHyperonAntiHyperonPairs(TCollision const& collision, THyperons const& fullHyperons, std::vector<int> selHypIndices, std::vector<int> selAntiHypIndices, float centrality, uint8_t gapSide, int type)
  {
    // 1st loop over all v0s/cascades
    for (std::size_t iHyp = 0; iHyp < selHypIndices.size(); iHyp++) {
      auto hyperon = fullHyperons.rawIteratorAt(selHypIndices[iHyp]);

      // 2nd loop over all v0s/cascade
      for (std::size_t iAntiHyp = 0; iAntiHyp < selAntiHypIndices.size(); iAntiHyp++) {
        // check we don't look at the same v0s/cascades
        if (selHypIndices[iHyp] == selAntiHypIndices[iAntiHyp]) {
          continue;
        }

        auto antiHyperon = fullHyperons.rawIteratorAt(selAntiHypIndices[iAntiHyp]);
        // check that the two hyperons have different daughter tracks
        if (!checkTrackIndices(hyperon, antiHyperon)) {
          continue;
        }

        // form V0 pairs and fill histograms
        analyseHyperonPairCandidate(collision, hyperon, antiHyperon, centrality, gapSide, type);
      } // end antiHyperon loop
    } // end hyperon loop

    // for (const auto& hyperon : fullHyperons) {
    //   // select only v0s matching Lambda selections
    //   if (!selHypIndices[hyperon.globalIndex() /*- fullHyperons.offset()*/]) { // local index needed due to collisions grouping
    //     continue;
    //   }

    //   // 2nd loop over all v0s/cascade
    //   for (const auto& antiHyperon : fullHyperons) {
    //     // select only v0s matching Anti-Lambda selections
    //     if (!selAntiHypIndices[antiHyperon.globalIndex() /*- fullHyperons.offset()*/]) { // local index needed due to collisions grouping
    //       continue;
    //     }

    //     // check we don't look at the same v0s/cascades
    //     if (hyperon.globalIndex() == antiHyperon.globalIndex()) {
    //       continue;
    //     }

    //     // check that the two hyperons have different daughter tracks
    //     if (!checkTrackIndices(hyperon, antiHyperon)) {
    //       continue;
    //     }

    //     // form V0 pairs and fill histograms
    //     analyseHyperonPairCandidate(collision, hyperon, antiHyperon, centrality, gapSide, type);
    //   } // end antiHyperon loop
    // } // end hyperon loop

    return;
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0Candidates const& fullV0s, CascadeCandidates const& fullCascades, DauTracks const&)
  {
    // Custom grouping
    v0sGrouped.clear();
    cascadesGrouped.clear();
    v0sGrouped.resize(collisions.size());
    cascadesGrouped.resize(collisions.size());

    for (const auto& v0 : fullV0s) {
      v0sGrouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }
    for (const auto& cascade : fullCascades) {
      cascadesGrouped[cascade.straCollisionId()].push_back(cascade.globalIndex());
    }

    for (const auto& collision : collisions) {
      // Fire up CCDB
      if (cfgSkimmedProcessing || doPairPropagationToDCA ||
          (mlConfigurations.useK0ShortScores && mlConfigurations.calculateK0ShortScores) ||
          (mlConfigurations.useLambdaScores && mlConfigurations.calculateLambdaScores) ||
          (mlConfigurations.useAntiLambdaScores && mlConfigurations.calculateAntiLambdaScores)) {
        initCCDB(collision);
      }

      if (!isEventAccepted(collision, true)) {
        continue;
      }

      if (cfgSkimmedProcessing) {
        zorro.isSelected(collision.globalBC()); /// Just let Zorro do the accounting
      }

      float centrality = -1;
      int selGapSide = -1; // only useful in case one wants to use this task in Pb-Pb UPC
      fillEventHistograms(collision, centrality, selGapSide);

      // __________________________________________
      // perform main analysis
      //
      if (buildK0sK0sPairs || buildLaLaBarPairs) { // Look at V0s
        std::size_t nV0sThisColl = v0sGrouped[collision.globalIndex()].size();
        selK0ShortIndices.clear();
        selLambdaIndices.clear();
        selAntiLambdaIndices.clear();
        for (std::size_t i = 0; i < nV0sThisColl; i++) {
          auto v0 = fullV0s.rawIteratorAt(v0sGrouped[collision.globalIndex()][i]);

          if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut || std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
            continue; // remove acceptance that's badly reproduced by MC / superfluous in future

          if (v0.v0Type() != v0Selections.v0TypeSelection && v0Selections.v0TypeSelection > -1)
            continue; // skip V0s that are not standard

          uint64_t selMap = computeReconstructionBitmap(v0, collision, v0.yLambda(), v0.yK0Short(), v0.pt());

          // consider for histograms for all species
          selMap = selMap | (static_cast<uint64_t>(1) << selConsiderK0Short) | (static_cast<uint64_t>(1) << selConsiderLambda) | (static_cast<uint64_t>(1) << selConsiderAntiLambda);
          selMap = selMap | (static_cast<uint64_t>(1) << selPhysPrimK0Short) | (static_cast<uint64_t>(1) << selPhysPrimLambda) | (static_cast<uint64_t>(1) << selPhysPrimAntiLambda);

          analyseV0Candidate(v0, v0.pt(), selMap, selK0ShortIndices, selLambdaIndices, selAntiLambdaIndices /*, fullV0s.offset()*/);
        } // end v0 loop

        // count the number of K0s, Lambda and AntiLambdas passsing the selections
        std::size_t nK0Shorts = selK0ShortIndices.size();
        std::size_t nLambdas = selLambdaIndices.size();
        std::size_t nAntiLambdas = selAntiLambdaIndices.size();

        if (buildK0sK0sPairs) {
          // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
          histos.fill(HIST("K0sK0s/h2dNbrOfK0ShortVsCentrality"), centrality, nK0Shorts);

          // Check the number of K0Short
          // needs at least 2 to form K0s-K0s pairs
          if (nK0Shorts >= 2) { // consider K0s K0s pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selK0ShortIndices, selK0ShortIndices, centrality, selGapSide, 0);
          }
        }

        if (buildLaLaBarPairs) {
          // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
          histos.fill(HIST("LaLaBar/h2dNbrOfLambdaVsCentrality"), centrality, nLambdas);
          histos.fill(HIST("LaLaBar/h2dNbrOfAntiLambdaVsCentrality"), centrality, nAntiLambdas);

          // Check the number of Lambdas and antiLambdas
          // needs at least 1 of each
          if (!buildSameSignPairs && nLambdas >= 1 && nAntiLambdas >= 1) { // consider Lambda antiLambda pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selLambdaIndices, selAntiLambdaIndices, centrality, selGapSide, 1);
          }
          if (buildSameSignPairs && nLambdas > 1) { // consider Lambda Lambda pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selLambdaIndices, selLambdaIndices, centrality, selGapSide, 1);
          }
          if (buildSameSignPairs && nAntiLambdas > 1) { // consider antiLambda antiLambda pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selAntiLambdaIndices, selAntiLambdaIndices, centrality, selGapSide, 1);
          }
        }
      }

      if (buildXiXiBarPairs || buildOmOmBarPairs) { // Look at Cascades
        std::size_t nCascadesThisColl = cascadesGrouped[collision.globalIndex()].size();

        selXiIndices.clear();
        selAntiXiIndices.clear();
        selOmIndices.clear();
        selAntiOmIndices.clear();
        for (std::size_t i = 0; i < nCascadesThisColl; i++) {
          auto cascade = fullCascades.rawIteratorAt(cascadesGrouped[collision.globalIndex()][i]);

          if (std::abs(cascade.negativeeta()) > cascSelections.daughterEtaCut ||
              std::abs(cascade.positiveeta()) > cascSelections.daughterEtaCut ||
              std::abs(cascade.bacheloreta()) > cascSelections.daughterEtaCut)
            continue; // remove acceptance that's badly reproduced by MC / superfluous in future

          if (buildXiXiBarPairs) {
            if (isCascadeSelected(cascade, collision, cascade.yXi(), true)) {
              if (cascade.sign() < 0) {
                selXiIndices.push_back(cascade.globalIndex());
              } else {
                selAntiXiIndices.push_back(cascade.globalIndex());
              }
            }
          }
          if (buildOmOmBarPairs) {
            if (isCascadeSelected(cascade, collision, cascade.yOmega(), false)) {
              if (cascade.sign() < 0) {
                selOmIndices.push_back(cascade.globalIndex());
              } else {
                selAntiOmIndices.push_back(cascade.globalIndex());
              }
            }
          }
        } // end cascade loop

        // count the number of Xi and antiXi passsing the selections
        std::size_t nXis = selXiIndices.size();
        std::size_t nAntiXis = selAntiXiIndices.size();
        std::size_t nOmegas = selOmIndices.size();
        std::size_t nAntiOmegas = selAntiOmIndices.size();

        // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
        if (buildXiXiBarPairs) {
          histos.fill(HIST("XiXiBar/h2dNbrOfXiVsCentrality"), centrality, nXis);
          histos.fill(HIST("XiXiBar/h2dNbrOfAntiXiVsCentrality"), centrality, nAntiXis);

          // Check the number of Lambdas and antiLambdas
          // needs at least 1 of each
          if (!buildSameSignPairs && nXis >= 1 && nAntiXis >= 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selXiIndices, selAntiXiIndices, centrality, selGapSide, 2);
          }
          if (buildSameSignPairs && nXis > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selXiIndices, selXiIndices, centrality, selGapSide, 2);
          }
          if (buildSameSignPairs && nAntiXis > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selAntiXiIndices, selAntiXiIndices, centrality, selGapSide, 2);
          }
        }
        if (buildOmOmBarPairs) {
          histos.fill(HIST("OmOmBar/h2dNbrOfOmegaVsCentrality"), centrality, nOmegas);
          histos.fill(HIST("OmOmBar/h2dNbrOfAntiOmegaVsCentrality"), centrality, nAntiOmegas);

          // Check the number of Lambdas and antiLambdas
          // needs at least 1 of each
          if (!buildSameSignPairs && nOmegas >= 1 && nAntiOmegas >= 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selOmIndices, selAntiOmIndices, centrality, selGapSide, 3);
          }
          if (buildSameSignPairs && nOmegas > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selOmIndices, selOmIndices, centrality, selGapSide, 3);
          }
          if (buildSameSignPairs && nAntiOmegas > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selAntiOmIndices, selAntiOmIndices, centrality, selGapSide, 3);
          }
        }
      }
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0MCCandidates const& fullV0s, CascadeMCCandidates const& fullCascades, DauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& /*mccollisions*/, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&)
  {
    // Custom grouping
    v0sGrouped.clear();
    cascadesGrouped.clear();
    v0sGrouped.resize(collisions.size());
    cascadesGrouped.resize(collisions.size());

    for (const auto& v0 : fullV0s) {
      v0sGrouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }
    for (const auto& cascade : fullCascades) {
      cascadesGrouped[cascade.straCollisionId()].push_back(cascade.globalIndex());
    }

    for (const auto& collision : collisions) {
      // Fire up CCDB
      if (cfgSkimmedProcessing || doPairPropagationToDCA ||
          (mlConfigurations.useK0ShortScores && mlConfigurations.calculateK0ShortScores) ||
          (mlConfigurations.useLambdaScores && mlConfigurations.calculateLambdaScores) ||
          (mlConfigurations.useAntiLambdaScores && mlConfigurations.calculateAntiLambdaScores)) {
        initCCDB(collision);
      }

      if (!isEventAccepted(collision, true)) {
        continue;
      }

      if (cfgSkimmedProcessing) {
        zorro.isSelected(collision.globalBC()); /// Just let Zorro do the accounting
      }

      float centrality = -1;
      int selGapSide = -1; // only useful in case one wants to use this task in Pb-Pb UPC
      fillEventHistograms(collision, centrality, selGapSide);

      // __________________________________________
      // perform main analysis
      if (buildK0sK0sPairs || buildLaLaBarPairs) { // Look at V0s
        std::size_t nV0sThisColl = v0sGrouped[collision.globalIndex()].size();
        selK0ShortIndices.clear();
        selLambdaIndices.clear();
        selAntiLambdaIndices.clear();

        for (std::size_t i = 0; i < nV0sThisColl; i++) {
          auto v0 = fullV0s.rawIteratorAt(v0sGrouped[collision.globalIndex()][i]);

          if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut || std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
            continue; // remove acceptance that's badly reproduced by MC / superfluous in future

          if (!v0.has_v0MCCore())
            continue;

          auto v0MC = v0.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

          float ptmc = RecoDecay::sqrtSumOfSquares(v0MC.pxPosMC() + v0MC.pxNegMC(), v0MC.pyPosMC() + v0MC.pyNegMC());
          float ymc = 1e-3;
          if (v0MC.pdgCode() == PDG_t::kK0Short)
            ymc = RecoDecay::y(std::array{v0MC.pxPosMC() + v0MC.pxNegMC(), v0MC.pyPosMC() + v0MC.pyNegMC(), v0MC.pzPosMC() + v0MC.pzNegMC()}, o2::constants::physics::MassKaonNeutral);
          else if (std::abs(v0MC.pdgCode()) == PDG_t::kLambda0)
            ymc = RecoDecay::y(std::array{v0MC.pxPosMC() + v0MC.pxNegMC(), v0MC.pyPosMC() + v0MC.pyNegMC(), v0MC.pzPosMC() + v0MC.pzNegMC()}, o2::constants::physics::MassLambda);

          uint64_t selMap = computeReconstructionBitmap(v0, collision, ymc, ymc, ptmc);
          selMap = selMap | computeMCAssociation(v0MC);

          // selMap |= maskTopological | maskTrackProperties | maskLambdaSpecific;
          // selMap |= maskTopological | maskTrackProperties | maskAntiLambdaSpecific;

          // consider only associated candidates if asked to do so, disregard association
          if (!doMCAssociation) {
            selMap = selMap | (static_cast<uint64_t>(1) << selConsiderK0Short) | (static_cast<uint64_t>(1) << selConsiderLambda) | (static_cast<uint64_t>(1) << selConsiderAntiLambda);
            selMap = selMap | (static_cast<uint64_t>(1) << selPhysPrimK0Short) | (static_cast<uint64_t>(1) << selPhysPrimLambda) | (static_cast<uint64_t>(1) << selPhysPrimAntiLambda);
          }

          analyseV0Candidate(v0, ptmc, selMap, selK0ShortIndices, selLambdaIndices, selAntiLambdaIndices /*, fullV0s.offset()*/);
        } // end v0 loop

        /// count the number of K0s, Lambda and AntiLambdas passsing the selections
        std::size_t nK0Shorts = selK0ShortIndices.size();
        std::size_t nLambdas = selLambdaIndices.size();
        std::size_t nAntiLambdas = selAntiLambdaIndices.size();

        if (buildK0sK0sPairs) {
          // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
          histos.fill(HIST("K0sK0s/h2dNbrOfK0ShortVsCentrality"), centrality, nK0Shorts);

          // Check the number of K0Short
          // needs at least 2 to form K0s-K0s pairs
          if (nK0Shorts >= 2) { // consider K0s K0s pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selK0ShortIndices, selK0ShortIndices, centrality, selGapSide, 0);
          }
        }

        if (buildLaLaBarPairs) {
          // fill the histograms with the number of reconstructed Lambda/antiLambda per collision
          histos.fill(HIST("LaLaBar/h2dNbrOfLambdaVsCentrality"), centrality, nLambdas);
          histos.fill(HIST("LaLaBar/h2dNbrOfAntiLambdaVsCentrality"), centrality, nAntiLambdas);

          if (!buildSameSignPairs && nLambdas >= 1 && nAntiLambdas >= 1) { // consider Lambda antiLambda pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selLambdaIndices, selAntiLambdaIndices, centrality, selGapSide, 1);
          }
          if (buildSameSignPairs && nLambdas > 1) { // consider Lambda Lambda pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selLambdaIndices, selLambdaIndices, centrality, selGapSide, 1);
          }
          if (buildSameSignPairs && nAntiLambdas > 1) { // consider antiLambda antiLambda pairs
            buildHyperonAntiHyperonPairs(collision, fullV0s, selAntiLambdaIndices, selAntiLambdaIndices, centrality, selGapSide, 1);
          }
        }
      }

      if (buildXiXiBarPairs || buildOmOmBarPairs) { // Look at Cascades
        std::size_t nCascadesThisColl = cascadesGrouped[collision.globalIndex()].size();

        selXiIndices.clear();
        selAntiXiIndices.clear();
        selOmIndices.clear();
        selAntiOmIndices.clear();
        for (std::size_t i = 0; i < nCascadesThisColl; i++) {
          auto cascade = fullCascades.rawIteratorAt(cascadesGrouped[collision.globalIndex()][i]);

          if (std::abs(cascade.negativeeta()) > cascSelections.daughterEtaCut ||
              std::abs(cascade.positiveeta()) > cascSelections.daughterEtaCut ||
              std::abs(cascade.bacheloreta()) > cascSelections.daughterEtaCut)
            continue; // remove acceptance that's badly reproduced by MC / superfluous in future

          if (!cascade.has_cascMCCore())
            continue;

          auto cascadeMC = cascade.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();

          float ymc = 1e-3;
          if (std::abs(cascadeMC.pdgCode()) == PDG_t::kXiMinus)
            ymc = RecoDecay::y(std::array{cascadeMC.pxMC(), cascadeMC.pyMC(), cascadeMC.pzMC()}, o2::constants::physics::MassXiMinus);
          else if (std::abs(cascadeMC.pdgCode()) == PDG_t::kOmegaMinus)
            ymc = RecoDecay::y(std::array{cascadeMC.pxMC(), cascadeMC.pyMC(), cascadeMC.pzMC()}, o2::constants::physics::MassOmegaMinus);

          if (buildXiXiBarPairs) {
            if (isCascadeSelected(cascade, collision, ymc, true)) {
              if (cascade.sign() < 0) {
                selXiIndices.push_back(cascade.globalIndex());
              } else {
                selAntiXiIndices.push_back(cascade.globalIndex());
              }
            }
          }
          if (buildOmOmBarPairs) {
            if (isCascadeSelected(cascade, collision, ymc, false)) {
              if (cascade.sign() < 0) {
                selOmIndices.push_back(cascade.globalIndex());
              } else {
                selAntiOmIndices.push_back(cascade.globalIndex());
              }
            }
          }
        } // end cascade loop

        // count the number of Xi and antiXi passsing the selections
        std::size_t nXis = selXiIndices.size();
        std::size_t nAntiXis = selAntiXiIndices.size();
        std::size_t nOmegas = selOmIndices.size();
        std::size_t nAntiOmegas = selAntiOmIndices.size();

        // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
        if (buildXiXiBarPairs) {
          histos.fill(HIST("XiXiBar/h2dNbrOfXiVsCentrality"), centrality, nXis);
          histos.fill(HIST("XiXiBar/h2dNbrOfAntiXiVsCentrality"), centrality, nAntiXis);

          // Check the number of Lambdas and antiLambdas
          // needs at least 1 of each
          if (!buildSameSignPairs && nXis >= 1 && nAntiXis >= 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selXiIndices, selAntiXiIndices, centrality, selGapSide, 2);
          }
          if (buildSameSignPairs && nXis > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selXiIndices, selXiIndices, centrality, selGapSide, 2);
          }
          if (buildSameSignPairs && nAntiXis > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selAntiXiIndices, selAntiXiIndices, centrality, selGapSide, 2);
          }
        }
        if (buildOmOmBarPairs) {
          histos.fill(HIST("OmOmBar/h2dNbrOfOmegaVsCentrality"), centrality, nOmegas);
          histos.fill(HIST("OmOmBar/h2dNbrOfAntiOmegaVsCentrality"), centrality, nAntiOmegas);

          // Check the number of Lambdas and antiLambdas
          // needs at least 1 of each
          if (!buildSameSignPairs && nOmegas >= 1 && nAntiOmegas >= 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selOmIndices, selAntiOmIndices, centrality, selGapSide, 3);
          }
          if (buildSameSignPairs && nOmegas > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selOmIndices, selOmIndices, centrality, selGapSide, 3);
          }
          if (buildSameSignPairs && nAntiOmegas > 1) {
            buildHyperonAntiHyperonPairs(collision, fullCascades, selAntiOmIndices, selAntiOmIndices, centrality, selGapSide, 3);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(QuarkoniaToHyperons, processRealData, "process as if real data", true);
  PROCESS_SWITCH(QuarkoniaToHyperons, processMonteCarlo, "process as if MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QuarkoniaToHyperons>(cfgc)};
}
