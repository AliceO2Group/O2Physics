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
/// \file forwardlambdakzeroanalysis.cxx
/// \brief V0s (K0s, Lambda and antiLambda) analysis task using derived data
///
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences & MBI
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences & MBI
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

#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/FwdDCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using namespace o2::aod::rctsel;

using Vec3D = ROOT::Math::SVector<double, 3>;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

// simple checkers, but ensure 64 bit integers
#define BITSET(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define BITCHECK(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

enum CentEstimator {
  kCentFT0C = 0,
  kCentFT0M,
  kCentFT0CVariant1,
  kCentMFT,
  kCentNGlobal,
  kCentFV0A
};

struct forwardlambdakzeroanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> analyseK0Short{"analyseK0Short", true, "process K0Short-like candidates"};
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", true, "process AntiLambda-like candidates"};
  Configurable<bool> analyseD0{"analyseD0", true, "process D0-like candidates"};
  Configurable<bool> analyseAntiD0{"analyseAntiD0", true, "process AntiD0-like candidates"};
  Configurable<bool> calculateFeeddownMatrix{"calculateFeeddownMatrix", true, "fill feeddown matrix if MC"};

  Configurable<bool> doPPAnalysis{"doPPAnalysis", false, "if in pp, set to true"};
  Configurable<std::string> irSource{"irSource", "", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0C, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFT0CVariant1, 3:CentMFT, 4:CentNGlobal, 5:CentFV0A)"};
  Configurable<bool> doUPCanalysis{"doUPCanalysis", true, "Study V0s in hadronic and UPC collisions"};

  Configurable<bool> doEventQA{"doEventQA", false, "do event QA histograms"};
  Configurable<bool> doCompleteTopoQA{"doCompleteTopoQA", false, "do topological variable QA histograms"};
  Configurable<bool> doPlainTopoQA{"doPlainTopoQA", true, "do simple 1D QA of candidates"};

  // for MC
  Configurable<bool> doTreatPiToMuon{"doTreatPiToMuon", false, "Take pi decay into muon into account in MC"};
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
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

  static constexpr float DefaultLifetimeCuts[1][3] = {{20., 30., 20.}};

  struct : ConfigurableGroup {
    std::string prefix = "v0Selections"; // JSON group name

    // Selection criteria: acceptance
    Configurable<float> rapidityMin{"rapidityMin", 0.5, "min rapidity"};
    Configurable<float> rapidityMax{"rapidityMax", 0.5, "max rapidity"};
    Configurable<float> daughterEtaCutMin{"daughterEtaCutMin", 0.8, "min eta for daughters"};
    Configurable<float> daughterEtaCutMax{"daughterEtaCutMax", 0.8, "max eta for daughters"};
    Configurable<float> minTrackPt{"minTrackPt", 0.8, "min track pT (GeV/c)"};

    // Standard 5 topological criteria
    Configurable<float> v0CosPA{"v0CosPA", 0.97, "min V0 CosPA"};
    Configurable<float> v0OpAngle{"v0OpAngle", 0., "min V0 opening angle"};
    Configurable<float> dcaV0Dau{"dcaV0Dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcaNegToPVxy{"dcaNegToPVxy", .05, "min DCAxy Neg To PV (cm)"};
    Configurable<float> dcaPosToPVxy{"dcaPosToPVxy", .05, "min DCAxy Pos To PV (cm)"};
    Configurable<float> dcaNegToPVz{"dcaNegToPVz", .05, "min DCAz Neg To PV (cm)"};
    Configurable<float> dcaPosToPVz{"dcaPosToPVz", .05, "min DCAz Pos To PV (cm)"};
    Configurable<float> v0RadiusMin{"v0RadiusMin", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0RadiusMax{"v0RadiusMax", 1e+09, "maximum V0 radius (cm)"};
    Configurable<float> v0Zmin{"v0Zmin", 1.2, "minimum V0 z (cm)"};
    Configurable<float> v0Zmax{"v0Zmax", 1e+09, "maximum V0 z (cm)"};
    Configurable<float> minPseudolifetime{"minPseudolifetime", -1e+09, "minimum V0 pseudo-proper lifetime (cm)"};
    Configurable<float> maxPseudolifetime{"maxPseudolifetime", 1e+09, "maximum V0 pseudo-proper lifetime (cm)"};
    Configurable<LabeledArray<float>> lifetimeCut{"lifetimeCut", {DefaultLifetimeCuts[0], 3, {"lifetimecutD0", "lifetimecutLambda", "lifetimecutK0S"}}, "lifetimeCut"};

    // invariant mass selection
    Configurable<float> compMassRejectionK0Short{"compMassRejectionK0Short", -1, "Competing K^{0}_{S} mass rejection (GeV/#it{c}^{2})"};
    Configurable<float> compMassRejectionLambda{"compMassRejectionLambda", -1, "Competing Lambda mass rejection (GeV/#it{c}^{2})"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // Track quality
    Configurable<int> minMFTclusters{"minMFTclusters", -1, "minimum MFT clusters"};
    Configurable<float> maxMFTchi2PerNcls{"maxMFTchi2PerNcls", 1e+09, "maximum MFT chi2 per clusters"};
    Configurable<float> maxMFTchi2{"maxMFTchi2", 1e+09, "max. MFT chi2"};
  } v0Selections;

  struct : ConfigurableGroup {
    std::string prefix = "rctConfigurations"; // JSON group name
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "", "Which detector condition requirements? (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "Include ZDC flags in the bit selection (for Pb-Pb only)"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } rctConfigurations;

  RCTFlagsChecker rctFlagsChecker{rctConfigurations.cfgRCTLabel.value};

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdbConfigurations"; // JSON group name
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

    // manual
    Configurable<bool> useCustomMagField{"useCustomMagField", false, "Use custom magnetic field value"};
    Configurable<float> customMagField{"customMagField", 5.0f, "Manually set magnetic field"};
  } ccdbConfigurations;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;
  int mRunNumber;
  float magField;
  std::map<std::string, std::string> metadata;
  o2::parameters::GRPMagField* grpmag = nullptr;

  // CCDB options
  struct : ConfigurableGroup {
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
    ConfigurableAxis axisPtXi{"axisPtXi", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for feeddown from Xi"};
    ConfigurableAxis axisPtCoarse{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
    ConfigurableAxis axisPtResol{"axisPtResol", {100, -1.0f, 1.0f}, "Axis for momentum resolution (GeV/c)"};
    ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, "Axis for K0S invariant mass (GeV/c2)"};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, "Axis for Lambda invariant mass (GeV/c2)"};
    ConfigurableAxis axisD0Mass{"axisD0Mass", {200, 1.75f, 1.95f}, "Axis for D0 invariant mass (GeV/c2)"};
    ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Centrality"};
    ConfigurableAxis axisCentralityFine{"axisCentralityFine", {101, 0.0f, 101.0f}, "Centrality"};
    ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles"};
    ConfigurableAxis axisIRBinning{"axisIRBinning", {500, 0, 50}, "Binning for the interaction rate (kHz)"};
    ConfigurableAxis axisMultFT0M{"axisMultFT0M", {500, 0.0f, +100000.0f}, "Multiplicity FT0M"};
    ConfigurableAxis axisMultFT0C{"axisMultFT0C", {500, 0.0f, +10000.0f}, "Multiplicity FT0C"};
    ConfigurableAxis axisMultFV0A{"axisMultFV0A", {500, 0.0f, +100000.0f}, "Multiplicity FV0A"};

    ConfigurableAxis axisRawCentrality{"axisRawCentrality", {VARIABLE_WIDTH, 0.000f, 52.320f, 75.400f, 95.719f, 115.364f, 135.211f, 155.791f, 177.504f, 200.686f, 225.641f, 252.645f, 281.906f, 313.850f, 348.302f, 385.732f, 426.307f, 470.146f, 517.555f, 568.899f, 624.177f, 684.021f, 748.734f, 818.078f, 892.577f, 973.087f, 1058.789f, 1150.915f, 1249.319f, 1354.279f, 1465.979f, 1584.790f, 1710.778f, 1844.863f, 1985.746f, 2134.643f, 2291.610f, 2456.943f, 2630.653f, 2813.959f, 3006.631f, 3207.229f, 3417.641f, 3637.318f, 3865.785f, 4104.997f, 4354.938f, 4615.786f, 4885.335f, 5166.555f, 5458.021f, 5762.584f, 6077.881f, 6406.834f, 6746.435f, 7097.958f, 7462.579f, 7839.165f, 8231.629f, 8635.640f, 9052.000f, 9484.268f, 9929.111f, 10389.350f, 10862.059f, 11352.185f, 11856.823f, 12380.371f, 12920.401f, 13476.971f, 14053.087f, 14646.190f, 15258.426f, 15890.617f, 16544.433f, 17218.024f, 17913.465f, 18631.374f, 19374.983f, 20136.700f, 20927.783f, 21746.796f, 22590.880f, 23465.734f, 24372.274f, 25314.351f, 26290.488f, 27300.899f, 28347.512f, 29436.133f, 30567.840f, 31746.818f, 32982.664f, 34276.329f, 35624.859f, 37042.588f, 38546.609f, 40139.742f, 41837.980f, 43679.429f, 45892.130f, 400000.000f}, "raw centrality signal"}; // for QA

    ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};

    // topological variable QA axes
    ConfigurableAxis axisDCAtoPVxy{"axisDCAtoPVxy", {100, 0.0f, 5.0f}, "DCAxy (cm)"};
    ConfigurableAxis axisDCAtoPVz{"axisDCAtoPVz", {100, -5.0f, 5.0f}, "DCAz (cm)"};
    ConfigurableAxis axisDCAdau{"axisDCAdau", {20, 0.0f, 2.0f}, "DCA (cm)"};
    ConfigurableAxis axisCosPA{"axisCosPA", {1000, 0.0f, 1.0f}, "cos(PA)"};
    ConfigurableAxis axisOpeningAngle{"axisOpeningAngle", {1000, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}, "Opening angle (rad)"};
    ConfigurableAxis axisV0Radius{"axisV0Radius", {20, 0.0f, 60.0f}, "V0 2D radius (cm)"};
    ConfigurableAxis axisV0Z{"axisV0Z", {100, -50.0f, 0.0f}, "V0 Z (cm)"};
    ConfigurableAxis axisV0Rapidity{"axisV0Rapidity", {100, -4.0f, -1.0f}, "Rapidity"};
    ConfigurableAxis axisV0ProperLifeTime{"axisV0ProperLifeTime", {100, 0.0f, 50.0f}, "Proper lifetime (cm)"};
    ConfigurableAxis axisV0PseudoProperLifeTime{"axisV0PseudoProperLifeTime", {100, -50.0f, 50.0f}, "Pseudo-proper lifetime (cm)"};

    ConfigurableAxis axisMFTchi2{"axisMFTchi2", {1000, 0.0f, 1000.0f}, "#chi^{2}"};
    ConfigurableAxis axisMFTchi2NDF{"axisMFTchi2NDF", {100, 0.0f, 100.0f}, "#chi^{2} per MFT clusters"};
    ConfigurableAxis axisMFTclus{"axisMFTclus", {10, 0.0f, 10.0f}, "N MFT Clusters"};

    // UPC axes
    ConfigurableAxis axisSelGap{"axisSelGap", {4, -1.5, 2.5}, "Gap side"};

    // AP plot axes
    ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
    ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

    // MC coll assoc QA axis
    ConfigurableAxis axisMonteCarloNch{"axisMonteCarloNch", {300, 0.0f, 3000.0f}, "N_{ch} MC"};
  } axisConfigurations;

  // UPC selections
  SGSelector sgSelector;
  struct : ConfigurableGroup {
    std::string prefix = "upcCuts"; // JSON group name
    Configurable<float> fv0Cut{"fv0Cut", 100., "FV0A threshold"};
    Configurable<float> ft0Acut{"ft0Acut", 200., "FT0A threshold"};
    Configurable<float> ft0Ccut{"ft0Ccut", 100., "FT0C threshold"};
    Configurable<float> zdcCut{"zdcCut", 10., "ZDC threshold"};
    // Configurable<float> gapSel{"gapSel", 2, "Gap selection"};
  } upcCuts;

  // DCA fitter configurations
  struct : ConfigurableGroup {
    std::string prefix = "fitterConfigurations"; // JSON group name
    Configurable<bool> propagateToPCA{"propagateToPCA", true, "Propagate to PCA?"};
    Configurable<float> minParamChange{"minParamChange", 4., "Stop minimization iterations if largest change of any X is smaller than this."};
    Configurable<float> minRelChi2Change{"minRelChi2Change", 0.9, "Stop iterations is chi2/chi2old > this"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Use abs. distance minimization rather than chi2"};
    Configurable<bool> useLUTMatCorr{"useLUTMatCorr", true, "Use material LUT correction, instead of TGeo material correction or no material correction."};
    Configurable<bool> useTGeoMatCorr{"useTGeoMatCorr", false, "Use material TGeo correction, instead of LUT material correction or no material correction."};
  } fitterConfigurations;

  o2::base::MatLayerCylSet* lut; // material LUT for DCA fitter
  o2::vertexing::FwdDCAFitterN<2> fitter;

  // Taken from https://github.com/AliceO2Group/O2Physics/blob/master/PWGLF/TableProducer/Strangeness/sigma0builder.cxx#L319
  // Thanks Gianni!
  // ______________________________________________________
  // Struct to store V0Pair properties
  struct PairTopoInfo {
    float X = -999.f;
    float Y = -999.f;
    float Z = -999.f;
    float Radius = -999.f;
    std::array<float, 3> positiveMomentum = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> negativeMomentum = {0.0f, 0.0f, 0.0f};
    float dcaPosToPVxy = -999.f;
    float dcaNegToPVxy = -999.f;
    float dcaPosToPVz = -999.f;
    float dcaNegToPVz = -999.f;
    float pTot = -999.f;
    float pT = -999.f;
    float pZ = -999.f;
    float DistOverTotMom = -1.f;
    float ZdistOverPz = -1.f;
    float DCADau = -999.f;
    float CosPA = -1.f;
    float OpAngle = -999.f;
    float mK0s = -999.f;
    float mLambda = -999.f;
    float mAntiLambda = -999.f;
    float mD0 = -999.f;
    float mAntiD0 = -999.f;
    float QtArm = -999.f;
    float AlphaArm = -999.f;
    float rapidityK0s = -999.f;
    float rapidityLambda = -999.f;
    float rapidityD0 = -999.f;
    int posNclusters = -1;
    int negNclusters = -1;
    float posChi2 = 999.f;
    float negChi2 = 999.f;
    float posChi2PerNclus = 999.f;
    float negChi2PerNclus = 999.f;

    // MC specific information
    int label = -1;
    int motherLabel = -1;
    float xMc = -999.f;
    float yMc = -999.f;
    float zMc = -999.f;
    std::array<float, 3> momentumMc = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> positiveMomentumMc = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> negativeMomentumMc = {0.0f, 0.0f, 0.0f};
    std::array<float, 3> momentumMotherMc = {0.0f, 0.0f, 0.0f};
    int processPositive = -1;
    int processNegative = -1;
    float pTotMc = -999.f;
    float pTMc = -999.f;
    float pZMc = -999.f;
    float rapMc = -999.f;
    int pdgCodePositive = 0;
    int pdgCodeNegative = 0;
    int pdgCode = 0;
    int pdgCodeMother = 0;
    int mcCollision = -1;
    bool isPhysicalPrimary = false;
    bool motherIsPhysicalPrimary = false;
  };

  // For manual sliceBy
  // Preslice<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels>> perMcCollision = o2::aod::mccollisionlabel::mcCollisionId;

  enum Selection : uint64_t { selCosPA = 0,
                              selOpAngle,
                              selRadius,
                              selRadiusMax,
                              selZmin,
                              selZmax,
                              selDCANegToPVxy,
                              selDCAPosToPVxy,
                              selDCANegToPVz,
                              selDCAPosToPVz,
                              selDCAV0Dau,
                              selK0ShortRapidityMin,
                              selK0ShortRapidityMax,
                              selLambdaRapidityMin,
                              selLambdaRapidityMax,
                              selD0RapidityMin,
                              selD0RapidityMax,
                              selK0ShortMassRejection,
                              selLambdaMassRejection,
                              selK0ShortCTau,
                              selLambdaCTau,
                              selD0CTau,
                              selK0ShortPseudoLifetimeMin,
                              selK0ShortPseudoLifetimeMax,
                              selLambdaPseudoLifetimeMin,
                              selLambdaPseudoLifetimeMax,
                              selD0PseudoLifetimeMin,
                              selD0PseudoLifetimeMax,
                              selK0ShortArmenteros,
                              selPosGoodMFTTrack,
                              selNegGoodMFTTrack,
                              selConsiderK0Short,    // for mc tagging
                              selConsiderLambda,     // for mc tagging
                              selConsiderAntiLambda, // for mc tagging
                              selConsiderD0,         // for mc tagging
                              selConsiderAntiD0,     // for mc tagging
                              selPhysPrimK0Short,    // for mc tagging
                              selPhysPrimLambda,     // for mc tagging
                              selPhysPrimAntiLambda, // for mc tagging
                              selPhysPrimD0,         // for mc tagging
                              selPhysPrimAntiD0,     // for mc tagging
  };

  uint64_t maskTopological;
  uint64_t maskTrackProperties;

  uint64_t maskK0ShortSpecific;
  uint64_t maskLambdaSpecific;
  uint64_t maskAntiLambdaSpecific;
  uint64_t maskD0Specific;
  uint64_t maskAntiD0Specific;

  uint64_t maskSelectionK0Short;
  uint64_t maskSelectionLambda;
  uint64_t maskSelectionAntiLambda;
  uint64_t maskSelectionD0;
  uint64_t maskSelectionAntiD0;

  uint64_t secondaryMaskSelectionLambda;
  uint64_t secondaryMaskSelectionAntiLambda;

  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL(ccdbConfigurations.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    // initialise bit masks
    // Mask with all topologic selections
    maskTopological = 0;
    BITSET(maskTopological, selCosPA);
    BITSET(maskTopological, selRadius);
    BITSET(maskTopological, selRadiusMax);
    BITSET(maskTopological, selDCANegToPVxy);
    BITSET(maskTopological, selDCAPosToPVxy);
    BITSET(maskTopological, selDCANegToPVz);
    BITSET(maskTopological, selDCAPosToPVz);
    BITSET(maskTopological, selDCAV0Dau);
    BITSET(maskTopological, selZmin);
    BITSET(maskTopological, selZmax);
    BITSET(maskTopological, selOpAngle);

    // Mask for specifically selecting K0Short
    maskK0ShortSpecific = 0;
    BITSET(maskK0ShortSpecific, selK0ShortRapidityMin);
    BITSET(maskK0ShortSpecific, selK0ShortRapidityMax);
    BITSET(maskK0ShortSpecific, selK0ShortCTau);
    BITSET(maskK0ShortSpecific, selK0ShortPseudoLifetimeMin);
    BITSET(maskK0ShortSpecific, selK0ShortPseudoLifetimeMax);
    BITSET(maskK0ShortSpecific, selK0ShortArmenteros);
    BITSET(maskK0ShortSpecific, selConsiderK0Short);
    BITSET(maskK0ShortSpecific, selLambdaMassRejection);
    // Mask for specifically selecting Lambda
    maskLambdaSpecific = 0;
    BITSET(maskLambdaSpecific, selLambdaRapidityMin);
    BITSET(maskLambdaSpecific, selLambdaRapidityMax);
    BITSET(maskLambdaSpecific, selLambdaCTau);
    BITSET(maskLambdaSpecific, selLambdaPseudoLifetimeMin);
    BITSET(maskLambdaSpecific, selLambdaPseudoLifetimeMax);
    BITSET(maskLambdaSpecific, selConsiderLambda);
    BITSET(maskLambdaSpecific, selK0ShortMassRejection);
    // Mask for specifically selecting AntiLambda
    maskAntiLambdaSpecific = 0;
    BITSET(maskAntiLambdaSpecific, selLambdaRapidityMin);
    BITSET(maskAntiLambdaSpecific, selLambdaRapidityMax);
    BITSET(maskAntiLambdaSpecific, selLambdaCTau);
    BITSET(maskAntiLambdaSpecific, selLambdaPseudoLifetimeMin);
    BITSET(maskAntiLambdaSpecific, selLambdaPseudoLifetimeMax);
    BITSET(maskAntiLambdaSpecific, selConsiderAntiLambda);
    BITSET(maskAntiLambdaSpecific, selK0ShortMassRejection);
    // Mask for specifically selecting D0
    maskD0Specific = 0;
    BITSET(maskD0Specific, selD0RapidityMin);
    BITSET(maskD0Specific, selD0RapidityMax);
    BITSET(maskD0Specific, selD0CTau);
    BITSET(maskD0Specific, selD0PseudoLifetimeMin);
    BITSET(maskD0Specific, selD0PseudoLifetimeMax);
    BITSET(maskD0Specific, selConsiderD0);
    BITSET(maskD0Specific, selK0ShortMassRejection);
    BITSET(maskD0Specific, selLambdaMassRejection);
    // Mask for specifically selecting D0
    maskAntiD0Specific = 0;
    BITSET(maskAntiD0Specific, selD0RapidityMin);
    BITSET(maskAntiD0Specific, selD0RapidityMax);
    BITSET(maskAntiD0Specific, selD0CTau);
    BITSET(maskAntiD0Specific, selD0PseudoLifetimeMin);
    BITSET(maskAntiD0Specific, selD0PseudoLifetimeMax);
    BITSET(maskAntiD0Specific, selConsiderAntiD0);
    BITSET(maskAntiD0Specific, selK0ShortMassRejection);
    BITSET(maskAntiD0Specific, selLambdaMassRejection);

    // ask for specific TPC/TOF PID selections
    maskTrackProperties = 0;
    BITSET(maskTrackProperties, selPosGoodMFTTrack);
    BITSET(maskTrackProperties, selNegGoodMFTTrack);

    // Primary particle selection, central to analysis
    maskSelectionK0Short = maskTopological | maskTrackProperties | maskK0ShortSpecific;
    maskSelectionLambda = maskTopological | maskTrackProperties | maskLambdaSpecific;
    maskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific;
    maskSelectionD0 = maskTopological | maskTrackProperties | maskD0Specific;
    maskSelectionAntiD0 = maskTopological | maskTrackProperties | maskAntiD0Specific;

    BITSET(maskSelectionK0Short, selPhysPrimK0Short);
    BITSET(maskSelectionLambda, selPhysPrimLambda);
    BITSET(maskSelectionAntiLambda, selPhysPrimAntiLambda);
    BITSET(maskSelectionD0, selPhysPrimD0);
    BITSET(maskSelectionAntiD0, selPhysPrimAntiD0);

    // No primary requirement for feeddown matrix
    secondaryMaskSelectionLambda = maskTopological | maskTrackProperties | maskLambdaSpecific;
    secondaryMaskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific;

    // Initialise the RCTFlagsChecker
    rctFlagsChecker.init(rctConfigurations.cfgRCTLabel.value, rctConfigurations.cfgCheckZDC, rctConfigurations.cfgTreatLimitedAcceptanceAsBad);

    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1D, {{23, -0.5f, +22.5f}});
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
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(17, "INEL>0");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(18, "INEL>1");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(19, "Below min occup.");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(20, "Above max occup.");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(21, "Below min IR");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(22, "Above max IR");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(23, "RCT flags");

    histos.add("hEventCentrality", "hEventCentrality", kTH1D, {axisConfigurations.axisCentralityFine});
    histos.add("hCentralityVsNch", "hCentralityVsNch", kTH2D, {axisConfigurations.axisCentralityFine, axisConfigurations.axisNch});

    histos.add("hEventPVz", "hEventPVz", kTH1D, {{100, -20.0f, +20.0f}});
    histos.add("hCentralityVsPVz", "hCentralityVsPVz", kTH2D, {axisConfigurations.axisCentralityFine, {100, -20.0f, +20.0f}});
    if (doprocessGenerated) {
      histos.add("hEventPVzMC", "hEventPVzMC", kTH1D, {{100, -20.0f, +20.0f}});
      histos.add("hCentralityVsPVzMC", "hCentralityVsPVzMC", kTH2D, {axisConfigurations.axisCentralityFine, {100, -20.0f, +20.0f}});
    }

    histos.add("hEventOccupancy", "hEventOccupancy", kTH1D, {axisConfigurations.axisOccupancy});
    histos.add("hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH2D, {axisConfigurations.axisCentralityFine, axisConfigurations.axisOccupancy});

    if (doUPCanalysis) {
      histos.add("hGapSide", "Gap side; Entries", kTH1D, {{5, -0.5, 4.5}});
      histos.add("hSelGapSide", "Selected gap side; Entries", kTH1D, {axisConfigurations.axisSelGap});
      histos.add("hEventCentralityVsSelGapSide", ";Centrality (%); Selected gap side", kTH2D, {axisConfigurations.axisCentralityFine, axisConfigurations.axisSelGap});
    }

    histos.add("hInteractionRate", "hInteractionRate", kTH1D, {axisConfigurations.axisIRBinning});
    histos.add("hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2D, {axisConfigurations.axisCentralityFine, axisConfigurations.axisIRBinning});

    histos.add("hInteractionRateVsOccupancy", "hInteractionRateVsOccupancy", kTH2D, {axisConfigurations.axisIRBinning, axisConfigurations.axisOccupancy});

    auto hSelectionV0s = histos.add<TH1>("GeneralQA/hSelectionV0s", "hSelectionV0s", kTH1D, {{static_cast<int>(selPhysPrimAntiD0) + 3, -0.5f, static_cast<double>(selPhysPrimAntiD0) + 2.5f}});
    hSelectionV0s->GetXaxis()->SetBinLabel(1, "All");
    hSelectionV0s->GetXaxis()->SetBinLabel(selCosPA + 2, "cosPA");
    hSelectionV0s->GetXaxis()->SetBinLabel(selOpAngle + 2, "Op. angle");
    hSelectionV0s->GetXaxis()->SetBinLabel(selRadius + 2, "Radius min.");
    hSelectionV0s->GetXaxis()->SetBinLabel(selRadiusMax + 2, "Radius max.");
    hSelectionV0s->GetXaxis()->SetBinLabel(selZmin + 2, "Z min.");
    hSelectionV0s->GetXaxis()->SetBinLabel(selZmax + 2, "Z max.");
    hSelectionV0s->GetXaxis()->SetBinLabel(selDCANegToPVxy + 2, "DCAxy neg. to PV");
    hSelectionV0s->GetXaxis()->SetBinLabel(selDCAPosToPVxy + 2, "DCAxy pos. to PV");
    hSelectionV0s->GetXaxis()->SetBinLabel(selDCANegToPVz + 2, "DCAz neg. to PV");
    hSelectionV0s->GetXaxis()->SetBinLabel(selDCAPosToPVz + 2, "DCAz pos. to PV");
    hSelectionV0s->GetXaxis()->SetBinLabel(selDCAV0Dau + 2, "DCA V0 dau.");
    hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortRapidityMin + 2, "K^{0}_{S} rap. min");
    hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortRapidityMax + 2, "K^{0}_{S} rap. max");
    hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaRapidityMin + 2, "#Lambda rap. min");
    hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaRapidityMax + 2, "#Lambda rap. max");
    hSelectionV0s->GetXaxis()->SetBinLabel(selD0RapidityMin + 2, "D^{0} rap. min");
    hSelectionV0s->GetXaxis()->SetBinLabel(selD0RapidityMax + 2, "D^{0} rap. max");
    hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortMassRejection + 2, "K^{0}_{S} mass rej.");
    hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaMassRejection + 2, "#Lambda mass rej.");
    hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortCTau + 2, "K^{0}_{S} lifetime");
    hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaCTau + 2, "#Lambda lifetime");
    hSelectionV0s->GetXaxis()->SetBinLabel(selD0CTau + 2, "D^{0} lifetime");
    hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortPseudoLifetimeMin + 2, "K^{0}_{S} pseudo-time min");
    hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortPseudoLifetimeMax + 2, "K^{0}_{S} pseudo-time max");
    hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaPseudoLifetimeMin + 2, "#Lambda pseudo-time min");
    hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaPseudoLifetimeMax + 2, "#Lambda pseudo-time max");
    hSelectionV0s->GetXaxis()->SetBinLabel(selD0PseudoLifetimeMin + 2, "D^{0} pseudo-time min");
    hSelectionV0s->GetXaxis()->SetBinLabel(selD0PseudoLifetimeMax + 2, "D^{0} pseudo-time max");
    hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortArmenteros + 2, "Arm. pod. cut");
    hSelectionV0s->GetXaxis()->SetBinLabel(selPosGoodMFTTrack + 2, "Pos. good MFT track");
    hSelectionV0s->GetXaxis()->SetBinLabel(selNegGoodMFTTrack + 2, "Neg. good MFT track");
    hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderK0Short + 2, "True K^{0}_{S}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderLambda + 2, "True #Lambda");
    hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderAntiLambda + 2, "True #bar{#Lambda}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderD0 + 2, "True D^{0}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderAntiD0 + 2, "True #bar{D}^{0}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimK0Short + 2, "Phys. prim. K^{0}_{S}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimLambda + 2, "Phys. prim. #Lambda");
    hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimAntiLambda + 2, "Phys. prim. #bar{#Lambda}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimD0 + 2, "Phys. prim. D^{0}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimAntiD0 + 2, "Phys. prim. #bar{D}^{0}");
    hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimAntiD0 + 3, "Cand. selected");

    // histograms versus mass
    if (analyseK0Short) {
      histos.add("h2dNbrOfK0ShortVsCentrality", "h2dNbrOfK0ShortVsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("h3dMassK0Short", "h3dMassK0Short", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisK0Mass});
      if (doUPCanalysis) {
        // Non-UPC info
        histos.add("h3dMassK0ShortHadronic", "h3dMassK0ShortHadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisK0Mass});
        // UPC info
        histos.add("h3dMassK0ShortSGA", "h3dMassK0ShortSGA", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisK0Mass});
        histos.add("h3dMassK0ShortSGC", "h3dMassK0ShortSGC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisK0Mass});
        histos.add("h3dMassK0ShortDG", "h3dMassK0ShortDG", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisK0Mass});
      }
    }
    if (analyseLambda) {
      histos.add("h2dNbrOfLambdaVsCentrality", "h2dNbrOfLambdaVsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("h3dMassLambda", "h3dMassLambda", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      if (doUPCanalysis) {
        // Non-UPC info
        histos.add("h3dMassLambdaHadronic", "h3dMassLambdaHadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        // UPC info
        histos.add("h3dMassLambdaSGA", "h3dMassLambdaSGA", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        histos.add("h3dMassLambdaSGC", "h3dMassLambdaSGC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        histos.add("h3dMassLambdaDG", "h3dMassLambdaDG", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      }
    }
    if (analyseAntiLambda) {
      histos.add("h2dNbrOfAntiLambdaVsCentrality", "h2dNbrOfAntiLambdaVsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      if (doUPCanalysis) {
        // Non-UPC info
        histos.add("h3dMassAntiLambdaHadronic", "h3dMassAntiLambdaHadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        // UPC info
        histos.add("h3dMassAntiLambdaSGA", "h3dMassAntiLambdaSGA", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        histos.add("h3dMassAntiLambdaSGC", "h3dMassAntiLambdaSGC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        histos.add("h3dMassAntiLambdaDG", "h3dMassAntiLambdaDG", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      }
    }
    if (analyseD0) {
      histos.add("h2dNbrOfD0VsCentrality", "h2dNbrOfD0VsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("h3dMassD0", "h3dMassD0", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
      if (doUPCanalysis) {
        // Non-UPC info
        histos.add("h3dMassD0Hadronic", "h3dMassD0Hadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
        // UPC info
        histos.add("h3dMassD0SGA", "h3dMassD0SGA", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
        histos.add("h3dMassD0SGC", "h3dMassD0SGC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
        histos.add("h3dMassD0DG", "h3dMassD0DG", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
      }
    }
    if (analyseAntiD0) {
      histos.add("h2dNbrOfAntiD0VsCentrality", "h2dNbrOfAntiD0VsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("h3dMassAntiD0", "h3dMassAntiD0", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
      if (doUPCanalysis) {
        // Non-UPC info
        histos.add("h3dMassAntiD0Hadronic", "h3dMassAntiLambdaHadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
        // UPC info
        histos.add("h3dMassAntiD0SGA", "h3dMassAntiD0SGA", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
        histos.add("h3dMassAntiD0SGC", "h3dMassAntiD0SGC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
        histos.add("h3dMassAntiD0DG", "h3dMassAntiD0DG", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisD0Mass});
      }
    }

    if (analyseLambda && calculateFeeddownMatrix && doprocessMonteCarlo) {
      histos.add("h3dLambdaFeeddown", "h3dLambdaFeeddown", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisPtXi});
      histos.add("h3dLambdaFeeddownFromXi0", "h3dLambdaFeeddownFromXi0", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisPtXi});
    }
    if (analyseAntiLambda && calculateFeeddownMatrix && doprocessMonteCarlo) {
      histos.add("h3dAntiLambdaFeeddown", "h3dAntiLambdaFeeddown", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisPtXi});
      histos.add("h3dAntiLambdaFeeddownFromXi0", "h3dAntiLambdaFeeddownFromXi0", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisPtXi});
    }

    if (analyseK0Short)
      histos.add("hMassK0Short", "hMassK0Short", kTH1D, {axisConfigurations.axisK0Mass});
    if (analyseLambda)
      histos.add("hMassLambda", "hMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
    if (analyseAntiLambda)
      histos.add("hMassAntiLambda", "hMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
    if (analyseD0)
      histos.add("hMassD0", "hMassD0", kTH1D, {axisConfigurations.axisD0Mass});
    if (analyseAntiD0)
      histos.add("hMassAntiD0", "hMassAntiD0", kTH1D, {axisConfigurations.axisD0Mass});

    if (doPlainTopoQA) {
      // All candidates received
      histos.add("hPosDCAToPVxy", "hPosDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
      histos.add("hNegDCAToPVxy", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVz});
      histos.add("hPosDCAToPVz", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
      histos.add("hNegDCAToPVz", "hNegDCAToPVz", kTH1D, {axisConfigurations.axisDCAtoPVz});
      histos.add("hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
      histos.add("hCosPA", "hCosPA", kTH1D, {axisConfigurations.axisCosPA});
      histos.add("hOpeningAngle", "hOpeningAngle", kTH1D, {axisConfigurations.axisOpeningAngle});
      histos.add("hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
      histos.add("hV0Z", "hV0Z", kTH1D, {axisConfigurations.axisV0Z});
      histos.add("hV0Rapidity", "hV0Rapidity", kTH1D, {axisConfigurations.axisV0Rapidity});
      histos.add("hV0LifetimeK0s", "hV0LifetimeK0s", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
      histos.add("hV0LifetimeLambda", "hV0LifetimeLambda", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
      histos.add("hV0LifetimeD0", "hV0LifetimeD0", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
      histos.add("hV0PseudoLifetimeK0s", "hV0PseudoLifetimeK0s", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
      histos.add("hV0PseudoLifetimeLambda", "hV0PseudoLifetimeLambda", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
      histos.add("hV0PseudoLifetimeD0", "hV0PseudoLifetimeD0", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
      histos.add("hV0InvMassK0s", "hV0InvMassK0s", kTH1D, {axisConfigurations.axisK0Mass});
      histos.add("hV0InvMassLambda", "hV0InvMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
      histos.add("hV0InvMassAntiLambda", "hV0InvMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
      histos.add("hV0InvMassD0", "hV0InvMassD0", kTH1D, {axisConfigurations.axisD0Mass});
      histos.add("hV0InvMassAntiD0", "hV0InvMassAntiD0", kTH1D, {axisConfigurations.axisD0Mass});

      histos.add("hPositiveMFTcls", "hPositiveMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
      histos.add("hNegativeMFTcls", "hNegativeMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
      histos.add("hPositiveMFTchi2PerNcls", "hPositiveMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
      histos.add("hNegativeMFTchi2PerNcls", "hNegativeMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
      histos.add("hPositiveMFTchi2", "hPositiveMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
      histos.add("hNegativeMFTchi2", "hNegativeMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
      if (doprocessMonteCarlo) {
        histos.add("hPositiveMomResolution", "hPositiveMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
        histos.add("hNegativeMomResolution", "hNegativeMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
      }
      if (analyseK0Short) {
        histos.add("K0Short/hPosDCAToPVxy", "hPosDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("K0Short/hNegDCAToPVxy", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("K0Short/hPosDCAToPVz", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("K0Short/hNegDCAToPVz", "hNegDCAToPVz", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("K0Short/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("K0Short/hCosPA", "hCosPA", kTH1D, {axisConfigurations.axisCosPA});
        histos.add("K0Short/hOpeningAngle", "hOpeningAngle", kTH1D, {axisConfigurations.axisOpeningAngle});
        histos.add("K0Short/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("K0Short/hV0Z", "hV0Z", kTH1D, {axisConfigurations.axisV0Z});
        histos.add("K0Short/hV0Rapidity", "hV0Rapidity", kTH1D, {axisConfigurations.axisV0Rapidity});
        histos.add("K0Short/hV0LifetimeK0s", "hV0LifetimeK0s", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("K0Short/hV0LifetimeLambda", "hV0LifetimeLambda", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("K0Short/hV0LifetimeD0", "hV0LifetimeD0", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("K0Short/hV0PseudoLifetimeK0s", "hV0PseudoLifetimeK0s", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("K0Short/hV0PseudoLifetimeLambda", "hV0PseudoLifetimeLambda", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("K0Short/hV0PseudoLifetimeD0", "hV0PseudoLifetimeD0", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("K0Short/hV0InvMassK0s", "hV0InvMassK0s", kTH1D, {axisConfigurations.axisK0Mass});
        histos.add("K0Short/hV0InvMassLambda", "hV0InvMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("K0Short/hV0InvMassAntiLambda", "hV0InvMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("K0Short/hV0InvMassD0", "hV0InvMassD0", kTH1D, {axisConfigurations.axisD0Mass});
        histos.add("K0Short/hV0InvMassAntiD0", "hV0InvMassAntiD0", kTH1D, {axisConfigurations.axisD0Mass});

        histos.add("K0Short/hPositiveMFTcls", "hPositiveMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("K0Short/hNegativeMFTcls", "hNegativeMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("K0Short/hPositiveMFTchi2PerNcls", "hPositiveMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("K0Short/hNegativeMFTchi2PerNcls", "hNegativeMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("K0Short/hPositiveMFTchi2", "hPositiveMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        histos.add("K0Short/hNegativeMFTchi2", "hNegativeMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        if (doprocessMonteCarlo) {
          histos.add("K0Short/hPositiveMomResolution", "hPositiveMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
          histos.add("K0Short/hNegativeMomResolution", "hNegativeMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
        }
      }
      if (analyseLambda) {
        histos.add("Lambda/hPosDCAToPVxy", "hPosDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("Lambda/hNegDCAToPVxy", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("Lambda/hPosDCAToPVz", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("Lambda/hNegDCAToPVz", "hNegDCAToPVz", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("Lambda/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("Lambda/hCosPA", "hCosPA", kTH1D, {axisConfigurations.axisCosPA});
        histos.add("Lambda/hOpeningAngle", "hOpeningAngle", kTH1D, {axisConfigurations.axisOpeningAngle});
        histos.add("Lambda/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("Lambda/hV0Z", "hV0Z", kTH1D, {axisConfigurations.axisV0Z});
        histos.add("Lambda/hV0Rapidity", "hV0Rapidity", kTH1D, {axisConfigurations.axisV0Rapidity});
        histos.add("Lambda/hV0LifetimeK0s", "hV0LifetimeK0s", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("Lambda/hV0LifetimeLambda", "hV0LifetimeLambda", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("Lambda/hV0LifetimeD0", "hV0LifetimeD0", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("Lambda/hV0PseudoLifetimeK0s", "hV0PseudoLifetimeK0s", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("Lambda/hV0PseudoLifetimeLambda", "hV0PseudoLifetimeLambda", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("Lambda/hV0PseudoLifetimeD0", "hV0PseudoLifetimeD0", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("Lambda/hV0InvMassK0s", "hV0InvMassK0s", kTH1D, {axisConfigurations.axisK0Mass});
        histos.add("Lambda/hV0InvMassLambda", "hV0InvMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("Lambda/hV0InvMassAntiLambda", "hV0InvMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("Lambda/hV0InvMassD0", "hV0InvMassD0", kTH1D, {axisConfigurations.axisD0Mass});
        histos.add("Lambda/hV0InvMassAntiD0", "hV0InvMassAntiD0", kTH1D, {axisConfigurations.axisD0Mass});

        histos.add("Lambda/hPositiveMFTcls", "hPositiveMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("Lambda/hNegativeMFTcls", "hNegativeMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("Lambda/hPositiveMFTchi2PerNcls", "hPositiveMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("Lambda/hNegativeMFTchi2PerNcls", "hNegativeMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("Lambda/hPositiveMFTchi2", "hPositiveMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        histos.add("Lambda/hNegativeMFTchi2", "hNegativeMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        if (doprocessMonteCarlo) {
          histos.add("Lambda/hPositiveMomResolution", "hPositiveMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
          histos.add("Lambda/hNegativeMomResolution", "hNegativeMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
        }
      }
      if (analyseAntiLambda) {
        histos.add("AntiLambda/hPosDCAToPVxy", "hPosDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("AntiLambda/hNegDCAToPVxy", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("AntiLambda/hPosDCAToPVz", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("AntiLambda/hNegDCAToPVz", "hNegDCAToPVz", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("AntiLambda/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("AntiLambda/hCosPA", "hCosPA", kTH1D, {axisConfigurations.axisCosPA});
        histos.add("AntiLambda/hOpeningAngle", "hOpeningAngle", kTH1D, {axisConfigurations.axisOpeningAngle});
        histos.add("AntiLambda/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("AntiLambda/hV0Z", "hV0Z", kTH1D, {axisConfigurations.axisV0Z});
        histos.add("AntiLambda/hV0Rapidity", "hV0Rapidity", kTH1D, {axisConfigurations.axisV0Rapidity});
        histos.add("AntiLambda/hV0LifetimeK0s", "hV0LifetimeK0s", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("AntiLambda/hV0LifetimeLambda", "hV0LifetimeLambda", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("AntiLambda/hV0LifetimeD0", "hV0LifetimeD0", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("AntiLambda/hV0PseudoLifetimeK0s", "hV0PseudoLifetimeK0s", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("AntiLambda/hV0PseudoLifetimeLambda", "hV0PseudoLifetimeLambda", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("AntiLambda/hV0PseudoLifetimeD0", "hV0PseudoLifetimeD0", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("AntiLambda/hV0InvMassK0s", "hV0InvMassK0s", kTH1D, {axisConfigurations.axisK0Mass});
        histos.add("AntiLambda/hV0InvMassLambda", "hV0InvMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("AntiLambda/hV0InvMassAntiLambda", "hV0InvMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("AntiLambda/hV0InvMassD0", "hV0InvMassD0", kTH1D, {axisConfigurations.axisD0Mass});
        histos.add("AntiLambda/hV0InvMassAntiD0", "hV0InvMassAntiD0", kTH1D, {axisConfigurations.axisD0Mass});

        histos.add("AntiLambda/hPositiveMFTcls", "hPositiveMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("AntiLambda/hNegativeMFTcls", "hNegativeMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("AntiLambda/hPositiveMFTchi2PerNcls", "hPositiveMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("AntiLambda/hNegativeMFTchi2PerNcls", "hNegativeMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("AntiLambda/hPositiveMFTchi2", "hPositiveMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        histos.add("AntiLambda/hNegativeMFTchi2", "hNegativeMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        if (doprocessMonteCarlo) {
          histos.add("AntiLambda/hPositiveMomResolution", "hPositiveMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
          histos.add("AntiLambda/hNegativeMomResolution", "hNegativeMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
        }
      }
      if (analyseD0) {
        histos.add("D0/hPosDCAToPVxy", "hPosDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("D0/hNegDCAToPVxy", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("D0/hPosDCAToPVz", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("D0/hNegDCAToPVz", "hNegDCAToPVz", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("D0/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("D0/hCosPA", "hCosPA", kTH1D, {axisConfigurations.axisCosPA});
        histos.add("D0/hOpeningAngle", "hOpeningAngle", kTH1D, {axisConfigurations.axisOpeningAngle});
        histos.add("D0/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("D0/hV0Z", "hV0Z", kTH1D, {axisConfigurations.axisV0Z});
        histos.add("D0/hV0Rapidity", "hV0Rapidity", kTH1D, {axisConfigurations.axisV0Rapidity});
        histos.add("D0/hV0LifetimeK0s", "hV0LifetimeK0s", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("D0/hV0LifetimeLambda", "hV0LifetimeLambda", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("D0/hV0LifetimeD0", "hV0LifetimeD0", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("D0/hV0PseudoLifetimeK0s", "hV0PseudoLifetimeK0s", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("D0/hV0PseudoLifetimeLambda", "hV0PseudoLifetimeLambda", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("D0/hV0PseudoLifetimeD0", "hV0PseudoLifetimeD0", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("D0/hV0InvMassK0s", "hV0InvMassK0s", kTH1D, {axisConfigurations.axisK0Mass});
        histos.add("D0/hV0InvMassLambda", "hV0InvMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("D0/hV0InvMassAntiLambda", "hV0InvMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("D0/hV0InvMassD0", "hV0InvMassD0", kTH1D, {axisConfigurations.axisD0Mass});
        histos.add("D0/hV0InvMassAntiD0", "hV0InvMassAntiD0", kTH1D, {axisConfigurations.axisD0Mass});

        histos.add("D0/hPositiveMFTcls", "hPositiveMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("D0/hNegativeMFTcls", "hNegativeMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("D0/hPositiveMFTchi2PerNcls", "hPositiveMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("D0/hNegativeMFTchi2PerNcls", "hNegativeMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("D0/hPositiveMFTchi2", "hPositiveMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        histos.add("D0/hNegativeMFTchi2", "hNegativeMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        if (doprocessMonteCarlo) {
          histos.add("D0/hPositiveMomResolution", "hPositiveMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
          histos.add("D0/hNegativeMomResolution", "hNegativeMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
        }
      }
      if (analyseAntiD0) {
        histos.add("AntiD0/hPosDCAToPVxy", "hPosDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("AntiD0/hNegDCAToPVxy", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("AntiD0/hPosDCAToPVz", "hNegDCAToPVxy", kTH1D, {axisConfigurations.axisDCAtoPVxy});
        histos.add("AntiD0/hNegDCAToPVz", "hNegDCAToPVz", kTH1D, {axisConfigurations.axisDCAtoPVz});
        histos.add("AntiD0/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("AntiD0/hCosPA", "hCosPA", kTH1D, {axisConfigurations.axisCosPA});
        histos.add("AntiD0/hOpeningAngle", "hOpeningAngle", kTH1D, {axisConfigurations.axisOpeningAngle});
        histos.add("AntiD0/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("AntiD0/hV0Z", "hV0Z", kTH1D, {axisConfigurations.axisV0Z});
        histos.add("AntiD0/hV0Rapidity", "hV0Rapidity", kTH1D, {axisConfigurations.axisV0Rapidity});
        histos.add("AntiD0/hV0LifetimeK0s", "hV0LifetimeK0s", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("AntiD0/hV0LifetimeLambda", "hV0LifetimeLambda", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("AntiD0/hV0LifetimeD0", "hV0LifetimeD0", kTH1D, {axisConfigurations.axisV0ProperLifeTime});
        histos.add("AntiD0/hV0PseudoLifetimeK0s", "hV0PseudoLifetimeK0s", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("AntiD0/hV0PseudoLifetimeLambda", "hV0PseudoLifetimeLambda", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("AntiD0/hV0PseudoLifetimeD0", "hV0PseudoLifetimeD0", kTH1D, {axisConfigurations.axisV0PseudoProperLifeTime});
        histos.add("AntiD0/hV0InvMassK0s", "hV0InvMassK0s", kTH1D, {axisConfigurations.axisK0Mass});
        histos.add("AntiD0/hV0InvMassLambda", "hV0InvMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("AntiD0/hV0InvMassAntiLambda", "hV0InvMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        histos.add("AntiD0/hV0InvMassD0", "hV0InvMassD0", kTH1D, {axisConfigurations.axisD0Mass});
        histos.add("AntiD0/hV0InvMassAntiD0", "hV0InvMassAntiD0", kTH1D, {axisConfigurations.axisD0Mass});

        histos.add("AntiD0/hPositiveMFTcls", "hPositiveMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("AntiD0/hNegativeMFTcls", "hNegativeMFTcls", kTH1D, {axisConfigurations.axisMFTclus});
        histos.add("AntiD0/hPositiveMFTchi2PerNcls", "hPositiveMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("AntiD0/hNegativeMFTchi2PerNcls", "hNegativeMFTchi2PerNcls", kTH1D, {axisConfigurations.axisMFTchi2NDF});
        histos.add("AntiD0/hPositiveMFTchi2", "hPositiveMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        histos.add("AntiD0/hNegativeMFTchi2", "hNegativeMFTchi2", kTH1D, {axisConfigurations.axisMFTchi2});
        if (doprocessMonteCarlo) {
          histos.add("AntiD0/hPositiveMomResolution", "hPositiveMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
          histos.add("AntiD0/hNegativeMomResolution", "hNegativeMomResolution", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPtResol});
        }
      }
    }

    // Check if doing the right thing in AP space please
    histos.add("GeneralQA/h2dArmenterosAll", "h2dArmenterosAll", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosK0sSelected", "h2dArmenterosK0sSelected", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosLambdaSelected", "h2dArmenterosLambdaSelected", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosD0Selected", "h2dArmenterosD0Selected", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});

    // Creation of histograms: MC generated
    if (doprocessGenerated) {
      histos.add("hGenEvents", "hGenEvents", kTH2D, {{axisConfigurations.axisNch}, {2, -0.5f, +1.5f}});
      histos.get<TH2>(HIST("hGenEvents"))->GetYaxis()->SetBinLabel(1, "All gen. events");
      histos.get<TH2>(HIST("hGenEvents"))->GetYaxis()->SetBinLabel(2, "Gen. with at least 1 rec. events");
      histos.add("hGenEventCentrality", "hGenEventCentrality", kTH1D, {{101, 0.0f, 101.0f}});

      histos.add("hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", kTH2D, {axisConfigurations.axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", kTH2D, {axisConfigurations.axisCentrality, {50, -0.5f, 49.5f}});

      histos.add("hCentralityVsMultMC", "hCentralityVsMultMC", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisNch});

      histos.add("h2dGenK0Short", "h2dGenK0Short;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});
      histos.add("h2dGenLambda", "h2dGenLambda;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});
      histos.add("h2dGenAntiLambda", "h2dGenAntiLambda;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});
      histos.add("h2dGenXiMinus", "h2dGenXiMinus;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});
      histos.add("h2dGenXiPlus", "h2dGenXiPlus;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});
      histos.add("h2dGenOmegaMinus", "h2dGenOmegaMinus;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});
      histos.add("h2dGenOmegaPlus", "h2dGenOmegaPlus;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});

      histos.add("h2dGenD0", "h2dGenD0;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});
      histos.add("h2dGenAntiD0", "h2dGenAntiD0;Centrality (%); #it{p}_{T} (GeV/#it{c})", kTH2D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt});

      histos.add("h2dGenK0ShortVsMultMC_RecoedEvt", "h2dGenK0ShortVsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenLambdaVsMultMC_RecoedEvt", "h2dGenLambdaVsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenAntiLambdaVsMultMC_RecoedEvt", "h2dGenAntiLambdaVsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenXiMinusVsMultMC_RecoedEvt", "h2dGenXiMinusVsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenXiPlusVsMultMC_RecoedEvt", "h2dGenXiPlusVsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenOmegaMinusVsMultMC_RecoedEvt", "h2dGenOmegaMinusVsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenOmegaPlusVsMultMC_RecoedEvt", "h2dGenOmegaPlusVsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});

      histos.add("h2dGenD0VsMultMC_RecoedEvt", "h2dGenD0VsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenAntiD0VsMultMC_RecoedEvt", "h2dGenAntiD0VsMultMC_RecoedEvt", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});

      histos.add("h2dGenK0ShortVsMultMC", "h2dGenK0ShortVsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenLambdaVsMultMC", "h2dGenLambdaVsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenAntiLambdaVsMultMC", "h2dGenAntiLambdaVsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenXiMinusVsMultMC", "h2dGenXiMinusVsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenXiPlusVsMultMC", "h2dGenXiPlusVsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenOmegaMinusVsMultMC", "h2dGenOmegaMinusVsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenOmegaPlusVsMultMC", "h2dGenOmegaPlusVsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});

      histos.add("h2dGenD0VsMultMC", "h2dGenD0VsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
      histos.add("h2dGenAntiD0VsMultMC", "h2dGenAntiD0VsMultMC", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisPt});
    }

    if (fitterConfigurations.useLUTMatCorr && fitterConfigurations.useTGeoMatCorr) {
      LOG(fatal) << "Cannot run with both useLUTMatCorr = on and useTGeoMatCorr = on. Please check your configuration!";
    }

    // standards hardcoded in builder ...
    // ...but can be changed easily since fitter is public
    fitter.setPropagateToPCA(fitterConfigurations.propagateToPCA);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(fitterConfigurations.minParamChange);     // 1e-3 for DCAfitter ; 4.0 for FwdDCAfitter
    fitter.setMinRelChi2Change(fitterConfigurations.minRelChi2Change); // 0.9 for DCAfitter ; 1e-3 for FwdDCAfitter
    // fitter.setMaxDZIni(1e9);
    // fitter.setMaxDXYIni(4.0f);
    // fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(fitterConfigurations.useAbsDCA);
    // fitter.setWeightedFinalPCA(false);
    fitter.setTGeoMat(fitterConfigurations.useTGeoMatCorr);
    // LUT has to be loaded later
    lut = nullptr;
    // fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);

    // mag field has to be set later
    fitter.setBz(-999.9f); // will NOT make sense if not changed

    // inspect histogram sizes, please
    histos.print();
  }

  // ______________________________________________________
  // Return centrality estimate for a given collision.
  // If takeMcCentrality is enabled, the centrality is taken from the MC collision; otherwise it is taken
  // from the reconstructed collision. Returns -1 if no corresponding centrality estimator is found or if no MC collision is associated to the recoed collision.
  template <typename TCollision>
  auto getCentralityRun3(TCollision const& collision)
  {
    // Helper lambda to extract centrality from any object exposing the cent* columns
    auto extractCentrality = [this](auto const& coll) -> float {
      switch (centralityEstimator) {
        case kCentFT0C:
          return coll.centFT0C();
        case kCentFT0M:
          return coll.centFT0M();
        case kCentFT0CVariant1:
          return coll.centFT0CVariant1();
        // case kCentMFT:          return coll.centMFT();
        case kCentNGlobal:
          return coll.centNGlobal();
        case kCentFV0A:
          return coll.centFV0A();
        default:
          return -1.f;
      }
    };
    return extractCentrality(collision);
  }

  // ______________________________________________________
  // Return slicing output
  template <typename TCollisions>
  auto getGroupedCollisions(TCollisions const& collisions, int globalIndex)
  {
    return collisions.sliceBy(perMcCollision, globalIndex);
  }

  template <typename TBCs, typename TCollision>
  void initCCDB(TBCs const& bcs, TCollision collision)
  {
    auto bc = collision.template bc_as<TBCs>();
    if (!bcs.size()) {
      LOGF(warn, "No BC found, skipping this DF.");
      return; // signal to skip this DF
    }

    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();
    // Fetching magnetic field if requested
    // In case override, don't proceed, please - no CCDB access required
    if (ccdbConfigurations.useCustomMagField) {
      magField = ccdbConfigurations.customMagField;
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
    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << mRunNumber;
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->template getForRun<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath.value, mRunNumber));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }
    LOG(info) << "Setting global propagator material propagation LUT";
    o2::base::Propagator::Instance()->setMatLUT(lut);

    if (fitterConfigurations.useLUTMatCorr) {
      fitter.setMatLUT(lut);
    }
    fitter.setBz(magField);
  }

  template <typename TV0>
  uint64_t computeReconstructionBitmap(TV0 v0, float rapK0s, float rapLambda, float rapD0)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    //
    // Base topological variables
    //
    // v0 radius min/max selections
    if (v0.Radius > v0Selections.v0RadiusMin) {
      BITSET(bitMap, selRadius);
    }
    if (v0.Radius < v0Selections.v0RadiusMax) {
      BITSET(bitMap, selRadiusMax);
    }
    // v0 radius min/max selections
    if (v0.Z > v0Selections.v0Zmin) {
      BITSET(bitMap, selZmin);
    }
    if (v0.Z < v0Selections.v0Zmax) {
      BITSET(bitMap, selZmax);
    }
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (std::fabs(v0.dcaPosToPVxy) > v0Selections.dcaPosToPVxy) {
      BITSET(bitMap, selDCAPosToPVxy);
    }
    if (std::fabs(v0.dcaNegToPVxy) > v0Selections.dcaNegToPVxy) {
      BITSET(bitMap, selDCANegToPVxy);
    }
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (std::fabs(v0.dcaPosToPVz) > v0Selections.dcaPosToPVz) {
      BITSET(bitMap, selDCAPosToPVz);
    }
    if (std::fabs(v0.dcaNegToPVz) > v0Selections.dcaNegToPVz) {
      BITSET(bitMap, selDCANegToPVz);
    }
    // V0 cosine of pointing angle
    if (v0.CosPA > v0Selections.v0CosPA) {
      BITSET(bitMap, selCosPA);
    }
    // V0 opening angle
    if (v0.OpAngle > v0Selections.v0OpAngle) {
      BITSET(bitMap, selOpAngle);
    }
    // DCA between v0 daughters
    if (v0.DCADau < v0Selections.dcaV0Dau) {
      BITSET(bitMap, selDCAV0Dau);
    }

    //
    // rapidity
    //
    if (rapK0s > v0Selections.rapidityMin) {
      BITSET(bitMap, selK0ShortRapidityMin);
    }
    if (rapK0s < v0Selections.rapidityMax) {
      BITSET(bitMap, selK0ShortRapidityMax);
    }
    if (rapLambda > v0Selections.rapidityMin) {
      BITSET(bitMap, selLambdaRapidityMin);
    }
    if (rapLambda < v0Selections.rapidityMax) {
      BITSET(bitMap, selLambdaRapidityMax);
    }
    if (rapD0 > v0Selections.rapidityMin) {
      BITSET(bitMap, selD0RapidityMin);
    }
    if (rapD0 < v0Selections.rapidityMax) {
      BITSET(bitMap, selD0RapidityMax);
    }

    //
    // competing mass rejection
    //
    if ((analyseD0 || analyseAntiD0) &&
        (std::fabs(v0.mK0s - o2::constants::physics::MassK0Short) > v0Selections.compMassRejectionK0Short)) {
      BITSET(bitMap, selK0ShortMassRejection);
    }
    if ((analyseD0 || analyseAntiD0) &&
        (std::fabs(v0.mLambda - o2::constants::physics::MassLambda0) > v0Selections.compMassRejectionLambda)) {
      BITSET(bitMap, selLambdaMassRejection);
    }
    if ((analyseLambda || analyseAntiLambda) && std::fabs(v0.mK0s - o2::constants::physics::MassK0Short) > v0Selections.compMassRejectionK0Short) {
      BITSET(bitMap, selK0ShortMassRejection);
    }
    if (analyseK0Short && std::fabs(v0.mLambda - o2::constants::physics::MassLambda0) > v0Selections.compMassRejectionLambda) {
      BITSET(bitMap, selLambdaMassRejection);
    }

    //
    // MFT quality flags
    //
    if (v0.posNclusters > v0Selections.minMFTclusters &&
        v0.posChi2 < v0Selections.maxMFTchi2 &&
        v0.posChi2PerNclus < v0Selections.maxMFTchi2PerNcls) {
      BITSET(bitMap, selPosGoodMFTTrack);
    }
    if (v0.negNclusters > v0Selections.minMFTclusters &&
        v0.negChi2 < v0Selections.maxMFTchi2 &&
        v0.negChi2PerNclus < v0Selections.maxMFTchi2PerNcls) {
      BITSET(bitMap, selNegGoodMFTTrack);
    }

    //
    // proper lifetime
    if ((analyseD0 || analyseAntiD0) && v0.DistOverTotMom * o2::constants::physics::MassD0 < v0Selections.lifetimeCut->get("lifetimecutD0")) {
      BITSET(bitMap, selD0CTau);
    }
    if ((analyseLambda || analyseAntiLambda) && v0.DistOverTotMom * o2::constants::physics::MassLambda0 < v0Selections.lifetimeCut->get("lifetimecutLambda")) {
      BITSET(bitMap, selLambdaCTau);
    }
    if (analyseK0Short && v0.DistOverTotMom * o2::constants::physics::MassK0Short < v0Selections.lifetimeCut->get("lifetimecutK0S")) {
      BITSET(bitMap, selK0ShortCTau);
    }

    //
    // Pseudo-proper lifetime
    if ((analyseD0 || analyseAntiD0) && v0.ZdistOverPz * o2::constants::physics::MassD0 > v0Selections.minPseudolifetime) {
      BITSET(bitMap, selD0PseudoLifetimeMin);
    }
    if ((analyseD0 || analyseAntiD0) && v0.ZdistOverPz * o2::constants::physics::MassD0 < v0Selections.maxPseudolifetime) {
      BITSET(bitMap, selD0PseudoLifetimeMax);
    }
    if ((analyseLambda || analyseAntiLambda) && v0.ZdistOverPz * o2::constants::physics::MassLambda0 > v0Selections.minPseudolifetime) {
      BITSET(bitMap, selLambdaPseudoLifetimeMin);
    }
    if ((analyseLambda || analyseAntiLambda) && v0.ZdistOverPz * o2::constants::physics::MassLambda0 < v0Selections.maxPseudolifetime) {
      BITSET(bitMap, selLambdaPseudoLifetimeMax);
    }
    if (analyseK0Short && v0.ZdistOverPz * o2::constants::physics::MassK0Short > v0Selections.minPseudolifetime) {
      BITSET(bitMap, selK0ShortPseudoLifetimeMin);
    }
    if (analyseK0Short && v0.ZdistOverPz * o2::constants::physics::MassK0Short < v0Selections.maxPseudolifetime) {
      BITSET(bitMap, selK0ShortPseudoLifetimeMax);
    }

    //
    // armenteros
    if (v0Selections.armPodCut < 1e-4 || v0.QtArm * v0Selections.armPodCut > std::abs(v0.AlphaArm)) {
      BITSET(bitMap, selK0ShortArmenteros);
    }

    return bitMap;
  }

  template <typename TV0>
  uint64_t computeMCAssociation(TV0 v0)
  // precalculate this information so that a check is one mask operation, not many
  {
    uint64_t bitMap = 0;
    bool isPositiveProton = v0.pdgCodePositive == PDG_t::kProton;
    bool isPositiveKaon = v0.pdgCodePositive == PDG_t::kKPlus;
    bool isPositivePion = v0.pdgCodePositive == PDG_t::kPiPlus || (doTreatPiToMuon && v0.pdgCodePositive == PDG_t::kMuonPlus);
    bool isNegativeProton = v0.pdgCodeNegative == PDG_t::kProtonBar;
    bool isNegativeKaon = v0.pdgCodeNegative == PDG_t::kKMinus;
    bool isNegativePion = v0.pdgCodeNegative == PDG_t::kPiMinus || (doTreatPiToMuon && v0.pdgCodeNegative == PDG_t::kMuonMinus);

    if (v0.pdgCode == PDG_t::kK0Short && isPositivePion && isNegativePion) {
      BITSET(bitMap, selConsiderK0Short);
      if (v0.isPhysicalPrimary)
        BITSET(bitMap, selPhysPrimK0Short);
    }
    if (v0.pdgCode == PDG_t::kLambda0 && isPositiveProton && isNegativePion) {
      BITSET(bitMap, selConsiderLambda);
      if (v0.isPhysicalPrimary)
        BITSET(bitMap, selPhysPrimLambda);
    }
    if (v0.pdgCode == PDG_t::kLambda0Bar && isPositivePion && isNegativeProton) {
      BITSET(bitMap, selConsiderAntiLambda);
      if (v0.isPhysicalPrimary)
        BITSET(bitMap, selPhysPrimAntiLambda);
    }
    if (v0.pdgCode == o2::constants::physics::kD0 && isPositivePion && isNegativeKaon) {
      BITSET(bitMap, selConsiderD0);
      if (v0.isPhysicalPrimary)
        BITSET(bitMap, selPhysPrimD0);
    }
    if (v0.pdgCode == o2::constants::physics::kD0Bar && isPositiveKaon && isNegativePion) {
      BITSET(bitMap, selConsiderAntiD0);
      if (v0.isPhysicalPrimary)
        BITSET(bitMap, selPhysPrimAntiD0);
    }
    return bitMap;
  }

  bool verifyMask(uint64_t bitmap, uint64_t mask)
  {
    return (bitmap & mask) == mask;
  }

  template <typename TV0>
  void analyseCandidate(TV0 v0, float pt, float centrality, uint64_t selMap, uint8_t gapSide, int& nK0Shorts, int& nLambdas, int& nAntiLambdas, int& nD0s, int& nAntiD0s)
  // precalculate this information so that a check is one mask operation, not many
  {
    bool passK0ShortSelections = false;
    bool passLambdaSelections = false;
    bool passAntiLambdaSelections = false;
    bool passD0Selections = false;
    bool passAntiD0Selections = false;

    // machine learning is on, go for calculation of thresholds
    // FIXME THIS NEEDS ADJUSTING
    passK0ShortSelections = verifyMask(selMap, maskSelectionK0Short);
    passLambdaSelections = verifyMask(selMap, maskSelectionLambda);
    passAntiLambdaSelections = verifyMask(selMap, maskSelectionAntiLambda);
    passD0Selections = verifyMask(selMap, maskSelectionD0);
    passAntiD0Selections = verifyMask(selMap, maskSelectionAntiD0);

    double invMassK0Short = v0.mK0s;
    double invMassLambda = v0.mLambda;
    double invMassAntiLambda = v0.mAntiLambda;
    double invMassD0 = v0.mD0;
    double invMassAntiD0 = v0.mAntiD0;

    // __________________________________________
    // fill with no selection if plain QA requested
    if (doPlainTopoQA) {
      histos.fill(HIST("hPosDCAToPVxy"), v0.dcaPosToPVxy);
      histos.fill(HIST("hNegDCAToPVxy"), v0.dcaNegToPVxy);
      histos.fill(HIST("hPosDCAToPVz"), v0.dcaPosToPVz);
      histos.fill(HIST("hNegDCAToPVz"), v0.dcaNegToPVz);
      histos.fill(HIST("hDCADaughters"), v0.DCADau);
      histos.fill(HIST("hCosPA"), v0.CosPA);
      histos.fill(HIST("hOpeningAngle"), v0.OpAngle);
      histos.fill(HIST("hV0Radius"), v0.Radius);
      histos.fill(HIST("hV0Z"), v0.Z);
      if (analyseK0Short) {
        histos.fill(HIST("hV0Rapidity"), v0.rapidityK0s);
      } else if (analyseLambda || analyseAntiLambda) {
        histos.fill(HIST("hV0Rapidity"), v0.rapidityLambda);
      } else if (analyseD0 || analyseAntiD0) {
        histos.fill(HIST("hV0Rapidity"), v0.rapidityD0);
      }
      histos.fill(HIST("hV0LifetimeK0s"), v0.DistOverTotMom * o2::constants::physics::MassK0Short);
      histos.fill(HIST("hV0LifetimeLambda"), v0.DistOverTotMom * o2::constants::physics::MassLambda);
      histos.fill(HIST("hV0LifetimeD0"), v0.DistOverTotMom * o2::constants::physics::MassD0);
      histos.fill(HIST("hV0PseudoLifetimeK0s"), v0.ZdistOverPz * o2::constants::physics::MassK0Short);
      histos.fill(HIST("hV0PseudoLifetimeLambda"), v0.ZdistOverPz * o2::constants::physics::MassLambda);
      histos.fill(HIST("hV0PseudoLifetimeD0"), v0.ZdistOverPz * o2::constants::physics::MassD0);
      histos.fill(HIST("hV0InvMassK0s"), invMassK0Short);
      histos.fill(HIST("hV0InvMassLambda"), invMassLambda);
      histos.fill(HIST("hV0InvMassAntiLambda"), invMassAntiLambda);
      histos.fill(HIST("hV0InvMassD0"), invMassD0);
      histos.fill(HIST("hV0InvMassAntiD0"), invMassAntiD0);

      histos.fill(HIST("hPositiveMFTcls"), v0.posNclusters);
      histos.fill(HIST("hNegativeMFTcls"), v0.negNclusters);
      histos.fill(HIST("hPositiveMFTchi2PerNcls"), v0.posChi2PerNclus);
      histos.fill(HIST("hNegativeMFTchi2PerNcls"), v0.negChi2PerNclus);
      histos.fill(HIST("hPositiveMFTchi2"), v0.posChi2);
      histos.fill(HIST("hNegativeMFTchi2"), v0.negChi2);
      if (doprocessMonteCarlo) {
        histos.fill(HIST("hPositiveMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
        histos.fill(HIST("hNegativeMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
      }
    }

    // Fill first bin: all candidates
    histos.fill(HIST("GeneralQA/hSelectionV0s"), 0);
    // Loop over all bits in the enum and fill if passed
    for (uint64_t i = 0; i <= selPhysPrimAntiD0; i++) {
      if (BITCHECK(selMap, i)) {
        histos.fill(HIST("GeneralQA/hSelectionV0s"), i + 1); // +1 because bin 0 = "All"
      }
    }

    // __________________________________________
    // main analysis
    if (passK0ShortSelections && analyseK0Short) {
      histos.fill(HIST("GeneralQA/hSelectionV0s"), selPhysPrimAntiD0 + 2);            //
      histos.fill(HIST("GeneralQA/h2dArmenterosK0sSelected"), v0.AlphaArm, v0.QtArm); // cross-check
      histos.fill(HIST("h3dMassK0Short"), centrality, pt, invMassK0Short);
      if (doUPCanalysis) {
        if (gapSide == 0)
          histos.fill(HIST("h3dMassK0ShortSGA"), centrality, pt, invMassK0Short);
        else if (gapSide == 1)
          histos.fill(HIST("h3dMassK0ShortSGC"), centrality, pt, invMassK0Short);
        else if (gapSide == 2)
          histos.fill(HIST("h3dMassK0ShortDG"), centrality, pt, invMassK0Short);
        else
          histos.fill(HIST("h3dMassK0ShortHadronic"), centrality, pt, invMassK0Short);
      }
      histos.fill(HIST("hMassK0Short"), invMassK0Short);
      if (doPlainTopoQA) {
        histos.fill(HIST("K0Short/hPosDCAToPVxy"), v0.dcaPosToPVxy);
        histos.fill(HIST("K0Short/hNegDCAToPVxy"), v0.dcaNegToPVxy);
        histos.fill(HIST("K0Short/hPosDCAToPVz"), v0.dcaPosToPVz);
        histos.fill(HIST("K0Short/hNegDCAToPVz"), v0.dcaNegToPVz);
        histos.fill(HIST("K0Short/hDCADaughters"), v0.DCADau);
        histos.fill(HIST("K0Short/hCosPA"), v0.CosPA);
        histos.fill(HIST("K0Short/hOpeningAngle"), v0.OpAngle);
        histos.fill(HIST("K0Short/hV0Radius"), v0.Radius);
        histos.fill(HIST("K0Short/hV0Z"), v0.Z);
        histos.fill(HIST("K0Short/hV0Rapidity"), v0.rapidityK0s);
        histos.fill(HIST("K0Short/hV0LifetimeK0s"), v0.DistOverTotMom * o2::constants::physics::MassK0Short);
        histos.fill(HIST("K0Short/hV0LifetimeLambda"), v0.DistOverTotMom * o2::constants::physics::MassLambda);
        histos.fill(HIST("K0Short/hV0LifetimeD0"), v0.DistOverTotMom * o2::constants::physics::MassD0);
        histos.fill(HIST("K0Short/hV0PseudoLifetimeK0s"), v0.ZdistOverPz * o2::constants::physics::MassK0Short);
        histos.fill(HIST("K0Short/hV0PseudoLifetimeLambda"), v0.ZdistOverPz * o2::constants::physics::MassLambda);
        histos.fill(HIST("K0Short/hV0PseudoLifetimeD0"), v0.ZdistOverPz * o2::constants::physics::MassD0);
        histos.fill(HIST("K0Short/hV0InvMassK0s"), invMassK0Short);
        histos.fill(HIST("K0Short/hV0InvMassLambda"), invMassLambda);
        histos.fill(HIST("K0Short/hV0InvMassAntiLambda"), invMassAntiLambda);
        histos.fill(HIST("K0Short/hV0InvMassD0"), invMassD0);
        histos.fill(HIST("K0Short/hV0InvMassAntiD0"), invMassAntiD0);

        histos.fill(HIST("K0Short/hPositiveMFTcls"), v0.posNclusters);
        histos.fill(HIST("K0Short/hNegativeMFTcls"), v0.negNclusters);
        histos.fill(HIST("K0Short/hPositiveMFTchi2PerNcls"), v0.posChi2PerNclus);
        histos.fill(HIST("K0Short/hNegativeMFTchi2PerNcls"), v0.negChi2PerNclus);
        histos.fill(HIST("K0Short/hPositiveMFTchi2"), v0.posChi2);
        histos.fill(HIST("K0Short/hNegativeMFTchi2"), v0.negChi2);
        if (doprocessMonteCarlo) {
          histos.fill(HIST("K0Short/hPositiveMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
          histos.fill(HIST("K0Short/hNegativeMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
        }
      }
      nK0Shorts++;
    }
    if (passLambdaSelections && analyseLambda) {
      histos.fill(HIST("GeneralQA/hSelectionV0s"), selPhysPrimAntiD0 + 2);               //
      histos.fill(HIST("GeneralQA/h2dArmenterosLambdaSelected"), v0.AlphaArm, v0.QtArm); // cross-check
      histos.fill(HIST("h3dMassLambda"), centrality, pt, invMassLambda);
      if (doUPCanalysis) {
        if (gapSide == 0)
          histos.fill(HIST("h3dMassLambdaSGA"), centrality, pt, invMassLambda);
        else if (gapSide == 1)
          histos.fill(HIST("h3dMassLambdaSGC"), centrality, pt, invMassLambda);
        else if (gapSide == 2)
          histos.fill(HIST("h3dMassLambdaDG"), centrality, pt, invMassLambda);
        else
          histos.fill(HIST("h3dMassLambdaHadronic"), centrality, pt, invMassLambda);
      }
      histos.fill(HIST("hMassLambda"), invMassLambda);
      if (doPlainTopoQA) {
        histos.fill(HIST("Lambda/hPosDCAToPVxy"), v0.dcaPosToPVxy);
        histos.fill(HIST("Lambda/hNegDCAToPVxy"), v0.dcaNegToPVxy);
        histos.fill(HIST("Lambda/hPosDCAToPVz"), v0.dcaPosToPVz);
        histos.fill(HIST("Lambda/hNegDCAToPVz"), v0.dcaNegToPVz);
        histos.fill(HIST("Lambda/hDCADaughters"), v0.DCADau);
        histos.fill(HIST("Lambda/hCosPA"), v0.CosPA);
        histos.fill(HIST("Lambda/hOpeningAngle"), v0.OpAngle);
        histos.fill(HIST("Lambda/hV0Radius"), v0.Radius);
        histos.fill(HIST("Lambda/hV0Z"), v0.Z);
        histos.fill(HIST("Lambda/hV0Rapidity"), v0.rapidityLambda);
        histos.fill(HIST("Lambda/hV0LifetimeK0s"), v0.DistOverTotMom * o2::constants::physics::MassK0Short);
        histos.fill(HIST("Lambda/hV0LifetimeLambda"), v0.DistOverTotMom * o2::constants::physics::MassLambda);
        histos.fill(HIST("Lambda/hV0LifetimeD0"), v0.DistOverTotMom * o2::constants::physics::MassD0);
        histos.fill(HIST("Lambda/hV0PseudoLifetimeK0s"), v0.ZdistOverPz * o2::constants::physics::MassK0Short);
        histos.fill(HIST("Lambda/hV0PseudoLifetimeLambda"), v0.ZdistOverPz * o2::constants::physics::MassLambda);
        histos.fill(HIST("Lambda/hV0PseudoLifetimeD0"), v0.ZdistOverPz * o2::constants::physics::MassD0);
        histos.fill(HIST("Lambda/hV0InvMassK0s"), invMassK0Short);
        histos.fill(HIST("Lambda/hV0InvMassLambda"), invMassLambda);
        histos.fill(HIST("Lambda/hV0InvMassAntiLambda"), invMassAntiLambda);
        histos.fill(HIST("Lambda/hV0InvMassD0"), invMassD0);
        histos.fill(HIST("Lambda/hV0InvMassAntiD0"), invMassAntiD0);

        histos.fill(HIST("Lambda/hPositiveMFTcls"), v0.posNclusters);
        histos.fill(HIST("Lambda/hNegativeMFTcls"), v0.negNclusters);
        histos.fill(HIST("Lambda/hPositiveMFTchi2PerNcls"), v0.posChi2PerNclus);
        histos.fill(HIST("Lambda/hNegativeMFTchi2PerNcls"), v0.negChi2PerNclus);
        histos.fill(HIST("Lambda/hPositiveMFTchi2"), v0.posChi2);
        histos.fill(HIST("Lambda/hNegativeMFTchi2"), v0.negChi2);
        if (doprocessMonteCarlo) {
          histos.fill(HIST("Lambda/hPositiveMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
          histos.fill(HIST("Lambda/hNegativeMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
        }
      }
      nLambdas++;
    }
    if (passAntiLambdaSelections && analyseAntiLambda) {
      histos.fill(HIST("GeneralQA/hSelectionV0s"), selPhysPrimAntiD0 + 2);               //
      histos.fill(HIST("GeneralQA/h2dArmenterosLambdaSelected"), v0.AlphaArm, v0.QtArm); // cross-check
      histos.fill(HIST("h3dMassAntiLambda"), centrality, pt, invMassAntiLambda);
      if (doUPCanalysis) {
        if (gapSide == 0)
          histos.fill(HIST("h3dMassAntiLambdaSGA"), centrality, pt, invMassAntiLambda);
        else if (gapSide == 1)
          histos.fill(HIST("h3dMassAntiLambdaSGC"), centrality, pt, invMassAntiLambda);
        else if (gapSide == 2)
          histos.fill(HIST("h3dMassAntiLambdaDG"), centrality, pt, invMassAntiLambda);
        else
          histos.fill(HIST("h3dMassAntiLambdaHadronic"), centrality, pt, invMassAntiLambda);
      }
      histos.fill(HIST("hMassAntiLambda"), invMassAntiLambda);
      if (doPlainTopoQA) {
        histos.fill(HIST("AntiLambda/hPosDCAToPVxy"), v0.dcaPosToPVxy);
        histos.fill(HIST("AntiLambda/hNegDCAToPVxy"), v0.dcaNegToPVxy);
        histos.fill(HIST("AntiLambda/hPosDCAToPVz"), v0.dcaPosToPVz);
        histos.fill(HIST("AntiLambda/hNegDCAToPVz"), v0.dcaNegToPVz);
        histos.fill(HIST("AntiLambda/hDCADaughters"), v0.DCADau);
        histos.fill(HIST("AntiLambda/hCosPA"), v0.CosPA);
        histos.fill(HIST("AntiLambda/hOpeningAngle"), v0.OpAngle);
        histos.fill(HIST("AntiLambda/hV0Radius"), v0.Radius);
        histos.fill(HIST("AntiLambda/hV0Z"), v0.Z);
        histos.fill(HIST("AntiLambda/hV0Rapidity"), v0.rapidityLambda);
        histos.fill(HIST("AntiLambda/hV0LifetimeK0s"), v0.DistOverTotMom * o2::constants::physics::MassK0Short);
        histos.fill(HIST("AntiLambda/hV0LifetimeLambda"), v0.DistOverTotMom * o2::constants::physics::MassLambda);
        histos.fill(HIST("AntiLambda/hV0LifetimeD0"), v0.DistOverTotMom * o2::constants::physics::MassD0);
        histos.fill(HIST("AntiLambda/hV0PseudoLifetimeK0s"), v0.ZdistOverPz * o2::constants::physics::MassK0Short);
        histos.fill(HIST("AntiLambda/hV0PseudoLifetimeLambda"), v0.ZdistOverPz * o2::constants::physics::MassLambda);
        histos.fill(HIST("AntiLambda/hV0PseudoLifetimeD0"), v0.ZdistOverPz * o2::constants::physics::MassD0);
        histos.fill(HIST("AntiLambda/hV0InvMassK0s"), invMassK0Short);
        histos.fill(HIST("AntiLambda/hV0InvMassLambda"), invMassLambda);
        histos.fill(HIST("AntiLambda/hV0InvMassAntiLambda"), invMassAntiLambda);
        histos.fill(HIST("AntiLambda/hV0InvMassD0"), invMassD0);
        histos.fill(HIST("AntiLambda/hV0InvMassAntiD0"), invMassAntiD0);

        histos.fill(HIST("AntiLambda/hPositiveMFTcls"), v0.posNclusters);
        histos.fill(HIST("AntiLambda/hNegativeMFTcls"), v0.negNclusters);
        histos.fill(HIST("AntiLambda/hPositiveMFTchi2PerNcls"), v0.posChi2PerNclus);
        histos.fill(HIST("AntiLambda/hNegativeMFTchi2PerNcls"), v0.negChi2PerNclus);
        histos.fill(HIST("AntiLambda/hPositiveMFTchi2"), v0.posChi2);
        histos.fill(HIST("AntiLambda/hNegativeMFTchi2"), v0.negChi2);
        if (doprocessMonteCarlo) {
          histos.fill(HIST("AntiLambda/hPositiveMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
          histos.fill(HIST("AntiLambda/hNegativeMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
        }
      }
      nAntiLambdas++;
    }
    if (passD0Selections && analyseD0) {
      histos.fill(HIST("GeneralQA/hSelectionV0s"), selPhysPrimAntiD0 + 2);           //
      histos.fill(HIST("GeneralQA/h2dArmenterosD0Selected"), v0.AlphaArm, v0.QtArm); // cross-check
      histos.fill(HIST("h3dMassD0"), centrality, pt, invMassD0);
      if (doUPCanalysis) {
        if (gapSide == 0)
          histos.fill(HIST("h3dMassD0SGA"), centrality, pt, invMassD0);
        else if (gapSide == 1)
          histos.fill(HIST("h3dMassD0SGC"), centrality, pt, invMassD0);
        else if (gapSide == 2)
          histos.fill(HIST("h3dMassD0DG"), centrality, pt, invMassD0);
        else
          histos.fill(HIST("h3dMassD0Hadronic"), centrality, pt, invMassD0);
      }
      histos.fill(HIST("hMassD0"), invMassD0);
      if (doPlainTopoQA) {
        histos.fill(HIST("D0/hPosDCAToPVxy"), v0.dcaPosToPVxy);
        histos.fill(HIST("D0/hNegDCAToPVxy"), v0.dcaNegToPVxy);
        histos.fill(HIST("D0/hPosDCAToPVz"), v0.dcaPosToPVz);
        histos.fill(HIST("D0/hNegDCAToPVz"), v0.dcaNegToPVz);
        histos.fill(HIST("D0/hDCADaughters"), v0.DCADau);
        histos.fill(HIST("D0/hCosPA"), v0.CosPA);
        histos.fill(HIST("D0/hOpeningAngle"), v0.OpAngle);
        histos.fill(HIST("D0/hV0Radius"), v0.Radius);
        histos.fill(HIST("D0/hV0Z"), v0.Z);
        histos.fill(HIST("D0/hV0Rapidity"), v0.rapidityD0);
        histos.fill(HIST("D0/hV0LifetimeK0s"), v0.DistOverTotMom * o2::constants::physics::MassK0Short);
        histos.fill(HIST("D0/hV0LifetimeLambda"), v0.DistOverTotMom * o2::constants::physics::MassLambda);
        histos.fill(HIST("D0/hV0LifetimeD0"), v0.DistOverTotMom * o2::constants::physics::MassD0);
        histos.fill(HIST("D0/hV0PseudoLifetimeK0s"), v0.ZdistOverPz * o2::constants::physics::MassK0Short);
        histos.fill(HIST("D0/hV0PseudoLifetimeLambda"), v0.ZdistOverPz * o2::constants::physics::MassLambda);
        histos.fill(HIST("D0/hV0PseudoLifetimeD0"), v0.ZdistOverPz * o2::constants::physics::MassD0);
        histos.fill(HIST("D0/hV0InvMassK0s"), invMassK0Short);
        histos.fill(HIST("D0/hV0InvMassLambda"), invMassLambda);
        histos.fill(HIST("D0/hV0InvMassAntiLambda"), invMassAntiLambda);
        histos.fill(HIST("D0/hV0InvMassD0"), invMassD0);
        histos.fill(HIST("D0/hV0InvMassAntiD0"), invMassAntiD0);

        histos.fill(HIST("D0/hPositiveMFTcls"), v0.posNclusters);
        histos.fill(HIST("D0/hNegativeMFTcls"), v0.negNclusters);
        histos.fill(HIST("D0/hPositiveMFTchi2PerNcls"), v0.posChi2PerNclus);
        histos.fill(HIST("D0/hNegativeMFTchi2PerNcls"), v0.negChi2PerNclus);
        histos.fill(HIST("D0/hPositiveMFTchi2"), v0.posChi2);
        histos.fill(HIST("D0/hNegativeMFTchi2"), v0.negChi2);
        if (doprocessMonteCarlo) {
          histos.fill(HIST("D0/hPositiveMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
          histos.fill(HIST("D0/hNegativeMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
        }
      }
      nD0s++;
    }
    if (passAntiD0Selections && analyseAntiD0) {
      histos.fill(HIST("GeneralQA/hSelectionV0s"), selPhysPrimAntiD0 + 2);           //
      histos.fill(HIST("GeneralQA/h2dArmenterosD0Selected"), v0.AlphaArm, v0.QtArm); // cross-check
      histos.fill(HIST("h3dMassAntiD0"), centrality, pt, invMassAntiD0);
      if (doUPCanalysis) {
        if (gapSide == 0)
          histos.fill(HIST("h3dMassAntiD0SGA"), centrality, pt, invMassAntiD0);
        else if (gapSide == 1)
          histos.fill(HIST("h3dMassAntiD0SGC"), centrality, pt, invMassAntiD0);
        else if (gapSide == 2)
          histos.fill(HIST("h3dMassAntiD0DG"), centrality, pt, invMassAntiD0);
        else
          histos.fill(HIST("h3dMassAntiD0Hadronic"), centrality, pt, invMassAntiD0);
      }
      histos.fill(HIST("hMassAntiD0"), invMassAntiD0);
      if (doPlainTopoQA) {
        histos.fill(HIST("AntiD0/hPosDCAToPVxy"), v0.dcaPosToPVxy);
        histos.fill(HIST("AntiD0/hNegDCAToPVxy"), v0.dcaNegToPVxy);
        histos.fill(HIST("AntiD0/hPosDCAToPVz"), v0.dcaPosToPVz);
        histos.fill(HIST("AntiD0/hNegDCAToPVz"), v0.dcaNegToPVz);
        histos.fill(HIST("AntiD0/hDCADaughters"), v0.DCADau);
        histos.fill(HIST("AntiD0/hCosPA"), v0.CosPA);
        histos.fill(HIST("AntiD0/hOpeningAngle"), v0.OpAngle);
        histos.fill(HIST("AntiD0/hV0Radius"), v0.Radius);
        histos.fill(HIST("AntiD0/hV0Z"), v0.Z);
        histos.fill(HIST("AntiD0/hV0Rapidity"), v0.rapidityD0);
        histos.fill(HIST("AntiD0/hV0LifetimeK0s"), v0.DistOverTotMom * o2::constants::physics::MassK0Short);
        histos.fill(HIST("AntiD0/hV0LifetimeLambda"), v0.DistOverTotMom * o2::constants::physics::MassLambda);
        histos.fill(HIST("AntiD0/hV0LifetimeD0"), v0.DistOverTotMom * o2::constants::physics::MassD0);
        histos.fill(HIST("AntiD0/hV0PseudoLifetimeK0s"), v0.ZdistOverPz * o2::constants::physics::MassK0Short);
        histos.fill(HIST("AntiD0/hV0PseudoLifetimeLambda"), v0.ZdistOverPz * o2::constants::physics::MassLambda);
        histos.fill(HIST("AntiD0/hV0PseudoLifetimeD0"), v0.ZdistOverPz * o2::constants::physics::MassD0);
        histos.fill(HIST("AntiD0/hV0InvMassK0s"), invMassK0Short);
        histos.fill(HIST("AntiD0/hV0InvMassLambda"), invMassLambda);
        histos.fill(HIST("AntiD0/hV0InvMassAntiLambda"), invMassAntiLambda);
        histos.fill(HIST("AntiD0/hV0InvMassD0"), invMassD0);
        histos.fill(HIST("AntiD0/hV0InvMassAntiD0"), invMassAntiD0);

        histos.fill(HIST("AntiD0/hPositiveMFTcls"), v0.posNclusters);
        histos.fill(HIST("AntiD0/hNegativeMFTcls"), v0.negNclusters);
        histos.fill(HIST("AntiD0/hPositiveMFTchi2PerNcls"), v0.posChi2PerNclus);
        histos.fill(HIST("AntiD0/hNegativeMFTchi2PerNcls"), v0.negChi2PerNclus);
        histos.fill(HIST("AntiD0/hPositiveMFTchi2"), v0.posChi2);
        histos.fill(HIST("AntiD0/hNegativeMFTchi2"), v0.negChi2);
        if (doprocessMonteCarlo) {
          histos.fill(HIST("AntiD0/hPositiveMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
          histos.fill(HIST("AntiD0/hNegativeMomResolution"), v0.pTMc, v0.pTMc - v0.pT);
        }
      }
      nAntiD0s++;
    }
  }

  template <typename TV0>
  void fillFeeddownMatrix(TV0 v0, float pt, float centrality, uint64_t selMap)
  // fill feeddown matrix for Lambdas or AntiLambdas
  // fixme: a potential improvement would be to consider mass windows for the l/al
  {
    if (v0.motherLabel < 0)
      return; // does not have mother particle in record, skip

    float rapidityXi = 999.f;
    if (std::abs(v0.pdgCodeMother) == PDG_t::kXiMinus)
      rapidityXi = RecoDecay::y(std::array{v0.momentumMotherMc[0], v0.momentumMotherMc[1], v0.momentumMotherMc[2]}, o2::constants::physics::MassXiMinus);
    if (std::abs(v0.pdgCodeMother) == o2::constants::physics::Pdg::kXi0)
      rapidityXi = RecoDecay::y(std::array{v0.momentumMotherMc[0], v0.momentumMotherMc[1], v0.momentumMotherMc[2]}, o2::constants::physics::MassXi0);

    if (std::fabs(rapidityXi) > 0.5f)
      return; // not a valid mother rapidity (PDG selection is later)

    // __________________________________________
    if (verifyMask(selMap, secondaryMaskSelectionLambda) && analyseLambda) {
      if (v0.motherIsPhysicalPrimary) {
        if (v0.pdgCodeMother == PDG_t::kXiMinus) {
          histos.fill(HIST("h3dLambdaFeeddown"), centrality, pt, std::hypot(v0.momentumMotherMc[0], v0.momentumMotherMc[1]));
        }
        if (v0.pdgCodeMother == PDG_t::kXiMinus || v0.pdgCodeMother == o2::constants::physics::Pdg::kXi0) {
          histos.fill(HIST("h3dLambdaFeeddownFromXi0"), centrality, pt, std::hypot(v0.momentumMotherMc[0], v0.momentumMotherMc[1]));
        }
      }
    }
    if (verifyMask(selMap, secondaryMaskSelectionAntiLambda) && analyseAntiLambda) {
      if (v0.motherIsPhysicalPrimary) {
        if (v0.pdgCodeMother == PDG_t::kXiPlusBar) {
          histos.fill(HIST("h3dAntiLambdaFeeddown"), centrality, pt, std::hypot(v0.momentumMotherMc[0], v0.momentumMotherMc[1]));
        }
        if (v0.pdgCodeMother == PDG_t::kXiPlusBar || v0.pdgCodeMother == -o2::constants::physics::Pdg::kXi0) {
          histos.fill(HIST("h3dAntiLambdaFeeddownFromXi0"), centrality, pt, std::hypot(v0.momentumMotherMc[0], v0.momentumMotherMc[1]));
        }
      }
    }
  }

  template <typename TCollision, typename TBCs>
  bool isEventAccepted(TCollision const& collision, TBCs const& /*bcs*/, bool fillHists)
  // check whether the collision passes our collision selections
  {
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    }
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    }

    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    }

    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
    }

    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
    }

    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
    }

    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
    }

    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
    }

    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
    }

    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
    }

    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
    }

    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);
    }

    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
    }

    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
    }

    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
    }

    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
    }

    if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 16 /* INEL > 0 */);
    }

    if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 17 /* INEL > 1 */);
    }

    float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 18 /* Below min occupancy */);
    }

    if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 19 /* Above max occupancy */);
    }

    // Fetch interaction rate only if required (in order to limit ccdb calls)
    auto bc = collision.template bc_as<TBCs>();
    double interactionRate = (eventSelections.minIR >= 0 || eventSelections.maxIR >= 0) ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3 : -1;
    if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 20 /* Below min IR */);
    }

    if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 21 /* Above max IR */);
    }

    if (!rctConfigurations.cfgRCTLabel.value.empty() && !rctFlagsChecker(collision)) {
      return false;
    }
    if (fillHists) {
      histos.fill(HIST("hEventSelection"), 22 /* Pass CBT condition */);
    }
    return true;
  }

  // ______________________________________________________
  // Simulated processing
  // Return the list of indices to the recoed collision associated to a given MC collision.
  template <typename TMCollisions, typename TCollisions, typename TBCs>
  std::vector<int> getListOfRecoCollIndices(TMCollisions const& mcCollisions, TCollisions const& collisions, TBCs const& bcs)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      auto groupedCollisions = getGroupedCollisions(collisions, mcCollision.globalIndex());
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      for (auto const& collision : groupedCollisions) {
        // consider event selections in the recoed <-> gen collision association, for the denominator (or numerator) of the efficiency (or signal loss)?
        if (eventSelections.useEvtSelInDenomEff) {
          if (!isEventAccepted(collision, bcs, false)) {
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
  // Reconstructed data processing
  // Fill reconstructed event information
  // Return centrality, occupancy, interaction rate, gap side and selGapside via reference-passing in arguments
  template <typename TCollision, typename TBCs>
  void fillReconstructedEventProperties(TCollision const& collision, TBCs const& /*bcs*/, float& centrality, float& collisionOccupancy, double& interactionRate, int& gapSide, int& selGapSide)
  {
    centrality = getCentralityRun3(collision);
    collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    // Fetch interaction rate only if required (in order to limit ccdb calls)
    auto bc = collision.template bc_as<TBCs>();
    interactionRate = !irSource.value.empty() ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3 : -1;

    if (doUPCanalysis) {
      // gap side
      // gapSide = collision.gapSide();
      gapSide = -1;
      // -1 --> Hadronic
      // 0 --> Single Gap - A side
      // 1 --> Single Gap - C side
      // 2 --> Double Gap - both A & C sides
      // selGapSide = sgSelector.trueGap(collision, upcCuts.fv0Cut, upcCuts.ft0Acut, upcCuts.ft0Ccut, upcCuts.zdcCut);

      histos.fill(HIST("hGapSide"), gapSide);
      histos.fill(HIST("hSelGapSide"), selGapSide);
      histos.fill(HIST("hEventCentralityVsSelGapSide"), centrality, selGapSide <= 2 ? selGapSide : -1);
    }

    histos.fill(HIST("hEventCentrality"), centrality);

    histos.fill(HIST("hCentralityVsNch"), centrality, collision.multNTracksPVeta1());
    if (doEventQA) {
      histos.fill(HIST("hCentralityVsNGlobal"), centrality, collision.multNTracksGlobal());
      histos.fill(HIST("hEventCentVsMultFT0M"), collision.centFT0M(), collision.multFT0A() + collision.multFT0C());
      histos.fill(HIST("hEventCentVsMultFT0C"), collision.centFT0C(), collision.multFT0C());
      histos.fill(HIST("hEventCentVsMultNGlobal"), collision.centNGlobal(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventCentVsMultFV0A"), collision.centFV0A(), collision.multFV0A());
      histos.fill(HIST("hEventMultFT0MvsMultNGlobal"), collision.multFT0A() + collision.multFT0C(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventMultFT0CvsMultNGlobal"), collision.multFT0C(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventMultFV0AvsMultNGlobal"), collision.multFV0A(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventMultPVvsMultNGlobal"), collision.multNTracksPVeta1(), collision.multNTracksGlobal());
      histos.fill(HIST("hEventMultFT0CvsMultFV0A"), collision.multFT0C(), collision.multFV0A());
    }

    histos.fill(HIST("hCentralityVsPVz"), centrality, collision.posZ());
    histos.fill(HIST("hEventPVz"), collision.posZ());

    histos.fill(HIST("hEventOccupancy"), collisionOccupancy);
    histos.fill(HIST("hCentralityVsOccupancy"), centrality, collisionOccupancy);

    histos.fill(HIST("hInteractionRate"), interactionRate);
    histos.fill(HIST("hCentralityVsInteractionRate"), centrality, interactionRate);

    histos.fill(HIST("hInteractionRateVsOccupancy"), interactionRate, collisionOccupancy);
    return;
  }

  // ______________________________________________________
  // Simulated processing
  // Fill generated event information (for event loss/splitting estimation)
  template <typename TMCCollisions, typename TCollisions, typename TBCs>
  void fillGeneratedEventProperties(TMCCollisions const& mcCollisions, TCollisions const& collisions, TBCs const& bcs)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      // Apply selections on MC collisions
      if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
        continue;
      }
      if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
        continue;
      }

      if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
        continue;
      }

      histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 0 /* all gen. events*/);

      auto groupedCollisions = getGroupedCollisions(collisions, mcCollision.globalIndex());
      // Check if there is at least one of the reconstructed collisions associated to this MC collision
      // If so, we consider it
      bool atLeastOne = false;
      int biggestNContribs = -1;
      float centrality = 100.5f;
      int nCollisions = 0;
      for (auto const& collision : groupedCollisions) {

        if (!isEventAccepted(collision, bcs, false)) {
          continue;
        }

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          centrality = getCentralityRun3(collision);
        }
        nCollisions++;

        atLeastOne = true;
      }

      histos.fill(HIST("hCentralityVsNcoll_beforeEvSel"), centrality, groupedCollisions.size());
      histos.fill(HIST("hCentralityVsNcoll_afterEvSel"), centrality, nCollisions);

      histos.fill(HIST("hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta05());
      histos.fill(HIST("hCentralityVsPVzMC"), centrality, mcCollision.posZ());
      histos.fill(HIST("hEventPVzMC"), mcCollision.posZ());

      if (atLeastOne) {
        histos.fill(HIST("hGenEvents"), mcCollision.multMCNParticlesEta05(), 1 /* at least 1 rec. event*/);

        histos.fill(HIST("hGenEventCentrality"), centrality);
      }
    }
    return;
  }

  // BUILDER PART: largely inspired of the strangenessBuilderModule
  // FIXME: to be put in a dedicated table producer if it works
  //__________________________________________________
  // MC kink handling
  template <typename mcpart>
  int getOriginatingParticle(mcpart const& part, int& indexForPositionOfDecay, bool treatPiToMuDecays)
  {
    int returnValue = -1;
    if (part.has_mothers()) {
      auto const& motherList = part.template mothers_as<aod::McParticles>();
      if (motherList.size() == 1) {
        for (const auto& mother : motherList) {
          if (std::abs(part.pdgCode()) == PDG_t::kMuonMinus && treatPiToMuDecays) {
            // muon decay, de-ref mother twice
            if (mother.has_mothers()) {
              auto grandMotherList = mother.template mothers_as<aod::McParticles>();
              if (grandMotherList.size() == 1) {
                for (const auto& grandMother : grandMotherList) {
                  returnValue = grandMother.globalIndex();
                  indexForPositionOfDecay = mother.globalIndex(); // for V0 decay position: grab muon
                }
              }
            }
          } else {
            returnValue = mother.globalIndex();
            indexForPositionOfDecay = part.globalIndex();
          }
        }
      }
    }
    return returnValue;
  }

  //___________________________________________________________________
  // Taken from https://github.com/AliceO2Group/AliceO2/blob/be5a2bb6c1be6757b7a496f9e90d470304d98fe7/Common/DCAFitter/include/DCAFitter/FwdDCAFitterN.h#L1281
  // Re-adapted for this task
  bool propagateToVtx(o2::track::TrackParCovFwd& t, const std::array<float, 3>& PV, const std::array<float, 2>& PVcov) const
  {
    // propagate track to vertex including MCS effects if material budget included, simple propagation to Z otherwise
    float x2x0 = 0;
    auto mb = lut->getMatBudget(t.getX(), t.getY(), t.getZ(), PV[0], PV[1], PV[2]);
    x2x0 = static_cast<float>(mb.meanX2X0);
    return t.propagateToVtxhelixWithMCS(PV[2], {PV[0], PV[1]}, PVcov, magField, x2x0);
  }

  template <typename TCollision, typename TTracks, typename TMFTTracks, typename TMCParticles>
  std::vector<PairTopoInfo> buildV0s(TCollision const& collision, TMFTTracks const& /*mftTracks*/, TTracks const& besttracks, TMCParticles const& mcParticles)
  {
    std::vector<PairTopoInfo> v0;
    for (const auto& [amft1, amft2] : combinations(besttracks, besttracks)) {
      auto mftPositive = amft1.template mfttrack_as<TMFTTracks>();
      auto mftNegative = amft2.template mfttrack_as<TMFTTracks>();

      if (mftPositive.eta() < v0Selections.daughterEtaCutMin)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future
      if (mftPositive.eta() > v0Selections.daughterEtaCutMax)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future
      if (mftNegative.eta() < v0Selections.daughterEtaCutMin)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future
      if (mftNegative.eta() > v0Selections.daughterEtaCutMax)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      if (mftPositive.pt() < v0Selections.minTrackPt)
        continue;
      if (mftNegative.pt() < v0Selections.minTrackPt)
        continue;

      // Consider only tracks of opposite charges
      if (mftPositive.sign() < 0)
        continue;
      if (mftNegative.sign() > 0)
        continue;

      SMatrix5 tpars1(mftPositive.x(), mftPositive.y(), mftPositive.phi(), mftPositive.tgl(), mftPositive.signed1Pt());
      std::vector<double> v1;
      SMatrix55 tcovs1(v1.begin(), v1.end());
      o2::track::TrackParCovFwd pars1{mftPositive.z(), tpars1, tcovs1, mftPositive.chi2()};
      o2::track::TrackParCovFwd pars1Copy{mftPositive.z(), tpars1, tcovs1, mftPositive.chi2()};

      SMatrix5 tpars2(mftNegative.x(), mftNegative.y(), mftNegative.phi(), mftNegative.tgl(), mftNegative.signed1Pt());
      std::vector<double> v2;
      SMatrix55 tcovs2(v2.begin(), v2.end());
      o2::track::TrackParCovFwd pars2{mftNegative.z(), tpars2, tcovs2, mftNegative.chi2()};
      o2::track::TrackParCovFwd pars2Copy{mftNegative.z(), tpars2, tcovs2, mftNegative.chi2()};

      propagateToVtx(pars1Copy, std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{collision.covXX(), collision.covYY()});
      propagateToVtx(pars2Copy, std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{collision.covXX(), collision.covYY()});
      float dcaPosToPVx = pars1Copy.getX() - collision.posX();
      float dcaPosToPVy = pars1Copy.getY() - collision.posY();
      float dcaPosToPVz = pars1Copy.getZ() - collision.posZ();

      float dcaNegToPVx = pars2Copy.getX() - collision.posX();
      float dcaNegToPVy = pars2Copy.getY() - collision.posY();
      float dcaNegToPVz = pars2Copy.getZ() - collision.posZ();

      // Move close to minima
      int nCand = 0;
      try {
        nCand = fitter.process(pars1, pars2);
      } catch (...) {
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      fitter.FwdpropagateTracksToVertex(); // propagate e and K to D vertex
      if (!fitter.isPropagateTracksToVertexDone()) {
        continue;
      }

      auto lTrack1 = fitter.getTrack(0);
      auto lTrack2 = fitter.getTrack(1);

      PairTopoInfo pairInfo;
      // DCA between cascade daughters
      pairInfo.DCADau = std::sqrt(fitter.getChi2AtPCACandidate());

      // get decay vertex coordinates
      Vec3D vtx = fitter.getPCACandidate();
      pairInfo.X = vtx[0];
      pairInfo.Y = vtx[1];
      pairInfo.Z = vtx[2];

      // get daughter DCA to PV
      pairInfo.dcaPosToPVxy = std::sqrt(dcaPosToPVx * dcaPosToPVx + dcaPosToPVy * dcaPosToPVy);
      pairInfo.dcaNegToPVxy = std::sqrt(dcaNegToPVx * dcaNegToPVx + dcaNegToPVy * dcaNegToPVy);
      pairInfo.dcaPosToPVz = dcaPosToPVz;
      pairInfo.dcaNegToPVz = dcaNegToPVz;

      // get daughter momenta
      pairInfo.positiveMomentum[0] = lTrack1.getPx();
      pairInfo.positiveMomentum[1] = lTrack1.getPy();
      pairInfo.positiveMomentum[2] = lTrack1.getPz();
      pairInfo.negativeMomentum[0] = lTrack2.getPx();
      pairInfo.negativeMomentum[1] = lTrack2.getPy();
      pairInfo.negativeMomentum[2] = lTrack2.getPz();

      pairInfo.CosPA = RecoDecay::cpa(
        std::array{collision.posX(), collision.posY(), collision.posZ()},
        std::array{vtx[0], vtx[1], vtx[2]},
        std::array{pairInfo.positiveMomentum[0] + pairInfo.negativeMomentum[0],
                   pairInfo.positiveMomentum[1] + pairInfo.negativeMomentum[1],
                   pairInfo.positiveMomentum[2] + pairInfo.negativeMomentum[2]});

      // Momenta
      TVector3 track1Momentum(pairInfo.positiveMomentum[0], pairInfo.positiveMomentum[1], pairInfo.positiveMomentum[2]);
      TVector3 track2Momentum(pairInfo.negativeMomentum[0], pairInfo.negativeMomentum[1], pairInfo.negativeMomentum[2]);

      pairInfo.OpAngle = track1Momentum.Angle(track2Momentum);

      // Radius
      pairInfo.Radius = std::sqrt(vtx[0] * vtx[0] + vtx[1] * vtx[1]);

      // Dist over tot mom.
      float px = pairInfo.positiveMomentum[0] + pairInfo.negativeMomentum[0];
      float py = pairInfo.positiveMomentum[1] + pairInfo.negativeMomentum[1];
      float pz = pairInfo.positiveMomentum[2] + pairInfo.negativeMomentum[2];
      pairInfo.DistOverTotMom = std::sqrt(vtx[0] * vtx[0] + vtx[1] * vtx[1] + vtx[2] * vtx[2]) / std::sqrt(px * px + py * py + pz * pz);

      // Z dist over pz
      pairInfo.ZdistOverPz = vtx[2] / pz;

      // V0 Momenta
      pairInfo.pT = std::sqrt(px * px + py * py);
      pairInfo.pTot = std::sqrt(px * px + py * py + pz * pz);
      pairInfo.pZ = pz;

      // Armenteros-Podolanski variables
      float momTot2 = RecoDecay::p2(px, py, pz);
      float dp = RecoDecay::dotProd(std::array{pairInfo.negativeMomentum[0], pairInfo.negativeMomentum[1], pairInfo.negativeMomentum[2]}, std::array{px, py, pz});
      pairInfo.QtArm = std::sqrt(RecoDecay::p2(pairInfo.negativeMomentum[0], pairInfo.negativeMomentum[1], pairInfo.negativeMomentum[2]) - dp * dp / momTot2); // qtarm

      float momTot = RecoDecay::p(px, py, pz);
      float lQlNeg = RecoDecay::dotProd(std::array{pairInfo.negativeMomentum[0], pairInfo.negativeMomentum[1], pairInfo.negativeMomentum[2]}, std::array{px, py, pz}) / momTot;
      float lQlPos = RecoDecay::dotProd(std::array{pairInfo.positiveMomentum[0], pairInfo.positiveMomentum[1], pairInfo.positiveMomentum[2]}, std::array{px, py, pz}) / momTot;
      pairInfo.AlphaArm = (lQlPos - lQlNeg) / (lQlPos + lQlNeg); // alphav0

      // Inv. mass
      pairInfo.mK0s = RecoDecay::m(std::array{pairInfo.positiveMomentum, pairInfo.negativeMomentum}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassPiMinus});
      pairInfo.mLambda = RecoDecay::m(std::array{pairInfo.positiveMomentum, pairInfo.negativeMomentum}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPiMinus});
      pairInfo.mAntiLambda = RecoDecay::m(std::array{pairInfo.positiveMomentum, pairInfo.negativeMomentum}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassProtonBar});
      pairInfo.mD0 = RecoDecay::m(std::array{pairInfo.positiveMomentum, pairInfo.negativeMomentum}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiMinus});
      pairInfo.mAntiD0 = RecoDecay::m(std::array{pairInfo.positiveMomentum, pairInfo.negativeMomentum}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKMinus});

      // Rapidity
      pairInfo.rapidityK0s = RecoDecay::y(std::array{pairInfo.positiveMomentum[0] + pairInfo.negativeMomentum[0], pairInfo.positiveMomentum[1] + pairInfo.negativeMomentum[1], pairInfo.positiveMomentum[2] + pairInfo.negativeMomentum[2]}, o2::constants::physics::MassKaonNeutral);
      pairInfo.rapidityLambda = RecoDecay::y(std::array{pairInfo.positiveMomentum[0] + pairInfo.negativeMomentum[0], pairInfo.positiveMomentum[1] + pairInfo.negativeMomentum[1], pairInfo.positiveMomentum[2] + pairInfo.negativeMomentum[2]}, o2::constants::physics::MassLambda);
      pairInfo.rapidityD0 = RecoDecay::y(std::array{pairInfo.positiveMomentum[0] + pairInfo.negativeMomentum[0], pairInfo.positiveMomentum[1] + pairInfo.negativeMomentum[1], pairInfo.positiveMomentum[2] + pairInfo.negativeMomentum[2]}, o2::constants::physics::MassD0);

      //
      pairInfo.posNclusters = mftPositive.nClusters();
      pairInfo.posChi2 = mftPositive.chi2();
      pairInfo.posChi2PerNclus = mftPositive.chi2() / mftPositive.nClusters();

      pairInfo.negNclusters = mftNegative.nClusters();
      pairInfo.negChi2 = mftNegative.chi2();
      pairInfo.negChi2PerNclus = mftNegative.chi2() / mftNegative.nClusters();

      //_________________________________________________________
      // MC handling part
      if constexpr (requires { mftNegative.mcParticleId(); mftPositive.mcParticleId(); }) {
        // Association check
        // There might be smarter ways of doing this in the future
        if (mftNegative.has_mcParticle() && mftPositive.has_mcParticle()) {
          auto lMCNegTrack = mftNegative.template mcParticle_as<aod::McParticles>();
          auto lMCPosTrack = mftPositive.template mcParticle_as<aod::McParticles>();

          pairInfo.pdgCodePositive = lMCPosTrack.pdgCode();
          pairInfo.pdgCodeNegative = lMCNegTrack.pdgCode();
          pairInfo.processPositive = lMCPosTrack.getProcess();
          pairInfo.processNegative = lMCNegTrack.getProcess();
          pairInfo.positiveMomentumMc[0] = lMCPosTrack.px();
          pairInfo.positiveMomentumMc[1] = lMCPosTrack.py();
          pairInfo.positiveMomentumMc[2] = lMCPosTrack.pz();
          pairInfo.negativeMomentumMc[0] = lMCNegTrack.px();
          pairInfo.negativeMomentumMc[1] = lMCNegTrack.py();
          pairInfo.negativeMomentumMc[2] = lMCNegTrack.pz();

          // check for pi -> mu + antineutrino decay
          // if present, de-reference original V0 correctly and provide label to original object
          // NOTA BENE: the prong info will still correspond to a muon, treat carefully!
          int negOriginating = -1, posOriginating = -1, particleForDecayPositionIdx = -1;
          negOriginating = getOriginatingParticle(lMCNegTrack, particleForDecayPositionIdx, doTreatPiToMuon);
          posOriginating = getOriginatingParticle(lMCPosTrack, particleForDecayPositionIdx, doTreatPiToMuon);

          if (negOriginating > -1 && negOriginating == posOriginating) {
            auto originatingV0 = mcParticles.rawIteratorAt(negOriginating);
            auto particleForDecayPosition = mcParticles.rawIteratorAt(particleForDecayPositionIdx);

            pairInfo.label = originatingV0.globalIndex();
            pairInfo.xMc = particleForDecayPosition.vx();
            pairInfo.yMc = particleForDecayPosition.vy();
            pairInfo.zMc = particleForDecayPosition.vz();

            if (originatingV0.has_mcCollision()) {
              pairInfo.mcCollision = originatingV0.mcCollisionId(); // save this reference, please
            }

            // acquire information
            pairInfo.pdgCode = originatingV0.pdgCode();
            pairInfo.isPhysicalPrimary = originatingV0.isPhysicalPrimary();
            pairInfo.momentumMc[0] = originatingV0.px();
            pairInfo.momentumMc[1] = originatingV0.py();
            pairInfo.momentumMc[2] = originatingV0.pz();
            pairInfo.pTotMc = std::sqrt(pairInfo.momentumMc[0] * pairInfo.momentumMc[0] + pairInfo.momentumMc[1] * pairInfo.momentumMc[1] + pairInfo.momentumMc[2] * pairInfo.momentumMc[2]);
            pairInfo.pTMc = std::sqrt(pairInfo.momentumMc[0] * pairInfo.momentumMc[0] + pairInfo.momentumMc[1] * pairInfo.momentumMc[1]);
            pairInfo.pZMc = pairInfo.momentumMc[2];

            if (pairInfo.pdgCode == PDG_t::kK0Short)
              pairInfo.rapMc = RecoDecay::y(std::array{pairInfo.momentumMc[0], pairInfo.momentumMc[1], pairInfo.momentumMc[2]}, o2::constants::physics::MassKaonNeutral);
            else if (std::abs(pairInfo.pdgCode) == PDG_t::kLambda0)
              pairInfo.rapMc = RecoDecay::y(std::array{pairInfo.momentumMc[0], pairInfo.momentumMc[1], pairInfo.momentumMc[2]}, o2::constants::physics::MassLambda);
            else if (std::abs(pairInfo.pdgCode) == o2::constants::physics::kD0)
              pairInfo.rapMc = RecoDecay::y(std::array{pairInfo.momentumMc[0], pairInfo.momentumMc[1], pairInfo.momentumMc[2]}, o2::constants::physics::MassD0);

            if (originatingV0.has_mothers()) {
              for (const auto& lV0Mother : originatingV0.template mothers_as<aod::McParticles>()) {
                pairInfo.pdgCodeMother = lV0Mother.pdgCode();
                pairInfo.motherLabel = lV0Mother.globalIndex();
                pairInfo.momentumMotherMc[0] = lV0Mother.px();
                pairInfo.momentumMotherMc[1] = lV0Mother.py();
                pairInfo.momentumMotherMc[2] = lV0Mother.pz();
                pairInfo.motherIsPhysicalPrimary = lV0Mother.isPhysicalPrimary();
              }
            }
          }
        } // end association check
      } // end MC handling

      v0.push_back(pairInfo);
    } // end v0 loop
    return v0;
  }

  // ______________________________________________________
  // Real data processing - no MC subscription
  template <typename TCollision, typename TMFTTracks, typename TTracks, typename TBCs>
  void analyzeRecoedV0sInRealData(TCollision const& collision, TMFTTracks const& mftTracks, TTracks const& besttracks, TBCs const& bcs)
  {
    // Fire up CCDB
    initCCDB(bcs, collision);

    if (!isEventAccepted(collision, bcs, true)) {
      return;
    }

    float centrality = -1;
    float collisionOccupancy = -2; // -1 already taken for the case where occupancy cannot be evaluated
    double interactionRate = -1;
    // gap side
    int gapSide = -1;
    int selGapSide = -1; // -1 --> Hadronic ; 0 --> Single Gap - A side ; 1 --> Single Gap - C side ; 2 --> Double Gap - both A & C sides
    // Fill recoed event properties
    fillReconstructedEventProperties(collision, bcs, centrality, collisionOccupancy, interactionRate, gapSide, selGapSide);

    histos.fill(HIST("hInteractionRateVsOccupancy"), interactionRate, collisionOccupancy);

    // __________________________________________
    // perform main analysis
    int nK0Shorts = 0;
    int nLambdas = 0;
    int nAntiLambdas = 0;
    int nD0s = 0;
    int nAntiD0s = 0;
    std::vector<PairTopoInfo> V0s = buildV0s(collision, mftTracks, besttracks, static_cast<TObject*>(nullptr));
    for (const auto& v0 : V0s) {
      // fill AP plot for all V0s
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.AlphaArm, v0.QtArm);

      uint64_t selMap = computeReconstructionBitmap(v0, v0.rapidityK0s, v0.rapidityLambda, v0.rapidityD0);

      // consider for histograms for all species
      BITSET(selMap, selConsiderK0Short);
      BITSET(selMap, selConsiderLambda);
      BITSET(selMap, selConsiderAntiLambda);
      BITSET(selMap, selConsiderD0);
      BITSET(selMap, selConsiderAntiD0);

      BITSET(selMap, selPhysPrimK0Short);
      BITSET(selMap, selPhysPrimLambda);
      BITSET(selMap, selPhysPrimAntiLambda);
      BITSET(selMap, selPhysPrimD0);
      BITSET(selMap, selPhysPrimAntiD0);

      analyseCandidate(v0, v0.pT, centrality, selMap, selGapSide, nK0Shorts, nLambdas, nAntiLambdas, nD0s, nAntiD0s);
    } // end v0 loop

    // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
    if (analyseK0Short) {
      histos.fill(HIST("h2dNbrOfK0ShortVsCentrality"), centrality, nK0Shorts);
    }
    if (analyseLambda) {
      histos.fill(HIST("h2dNbrOfLambdaVsCentrality"), centrality, nLambdas);
    }
    if (analyseAntiLambda) {
      histos.fill(HIST("h2dNbrOfAntiLambdaVsCentrality"), centrality, nAntiLambdas);
    }
    if (analyseD0) {
      histos.fill(HIST("h2dNbrOfD0VsCentrality"), centrality, nD0s);
    }
    if (analyseAntiD0) {
      histos.fill(HIST("h2dNbrOfAntiD0VsCentrality"), centrality, nAntiD0s);
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  template <typename TCollision, typename TMFTTracks, typename TTracks, typename TBCs, typename TMCParticles>
  void analyzeRecoedV0sInMonteCarlo(TCollision const& collision, TMFTTracks const& mftTracks, TTracks const& besttracks, TBCs const& bcs, TMCParticles const& mcParticles)
  {
    // Fire up CCDB
    initCCDB(bcs, collision);

    if (!isEventAccepted(collision, bcs, true)) {
      return;
    }

    float centrality = -1;
    float collisionOccupancy = -2; // -1 already taken for the case where occupancy cannot be evaluated
    double interactionRate = -1;
    // gap side
    int gapSide = -1;
    int selGapSide = -1; // -1 --> Hadronic ; 0 --> Single Gap - A side ; 1 --> Single Gap - C side ; 2 --> Double Gap - both A & C sides
    // Fill recoed event properties
    fillReconstructedEventProperties(collision, bcs, centrality, collisionOccupancy, interactionRate, gapSide, selGapSide);

    histos.fill(HIST("hInteractionRateVsOccupancy"), interactionRate, collisionOccupancy);

    // __________________________________________
    // perform main analysis
    int nK0Shorts = 0;
    int nLambdas = 0;
    int nAntiLambdas = 0;
    int nD0s = 0;
    int nAntiD0s = 0;
    std::vector<PairTopoInfo> V0s = buildV0s(collision, mftTracks, besttracks, mcParticles);
    for (const auto& v0 : V0s) {
      if (v0.label < 0) {
        continue;
      }

      // fill AP plot for all V0s
      histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.AlphaArm, v0.QtArm);

      uint64_t selMap = computeReconstructionBitmap(v0, v0.rapMc, v0.rapMc, v0.rapMc);
      selMap = selMap | computeMCAssociation(v0);

      // feeddown matrix always with association
      if (calculateFeeddownMatrix)
        fillFeeddownMatrix(v0, v0.pTMc, centrality, selMap);

      // consider only associated candidates if asked to do so, disregard association
      if (!doMCAssociation) {
        BITSET(selMap, selConsiderK0Short);
        BITSET(selMap, selConsiderLambda);
        BITSET(selMap, selConsiderAntiLambda);
        BITSET(selMap, selConsiderD0);
        BITSET(selMap, selConsiderAntiD0);

        BITSET(selMap, selPhysPrimK0Short);
        BITSET(selMap, selPhysPrimLambda);
        BITSET(selMap, selPhysPrimAntiLambda);
        BITSET(selMap, selPhysPrimD0);
        BITSET(selMap, selPhysPrimAntiD0);
      }

      analyseCandidate(v0, v0.pTMc, centrality, selMap, selGapSide, nK0Shorts, nLambdas, nAntiLambdas, nD0s, nAntiD0s);
    } // end v0 loop

    // fill the histograms with the number of reconstructed K0s/Lambda/antiLambda per collision
    if (analyseK0Short) {
      histos.fill(HIST("h2dNbrOfK0ShortVsCentrality"), centrality, nK0Shorts);
    }
    if (analyseLambda) {
      histos.fill(HIST("h2dNbrOfLambdaVsCentrality"), centrality, nLambdas);
    }
    if (analyseAntiLambda) {
      histos.fill(HIST("h2dNbrOfAntiLambdaVsCentrality"), centrality, nAntiLambdas);
    }
    if (analyseD0) {
      histos.fill(HIST("h2dNbrOfD0VsCentrality"), centrality, nD0s);
    }
    if (analyseAntiD0) {
      histos.fill(HIST("h2dNbrOfAntiD0VsCentrality"), centrality, nAntiD0s);
    }
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  template <typename TMCCollisions, typename TMcParticles, typename TCollisions, typename TBCs>
  void analyzeGeneratedV0s(TMCCollisions const& mcCollisions, TMcParticles const& mcParticles, TCollisions const& collisions, TBCs const& bcs)
  {
    fillGeneratedEventProperties(mcCollisions, collisions, bcs);
    std::vector<int> listBestCollisionIdx = getListOfRecoCollIndices(mcCollisions, collisions, bcs);
    for (auto const& mcParticle : mcParticles) {
      if (!mcParticle.has_mcCollision())
        continue;

      auto mcCollision = mcParticle.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
      if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
        continue;
      }
      if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
        continue;
      }

      if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
        continue;
      }

      if (!mcParticle.isPhysicalPrimary())
        continue;

      float ptmc = mcParticle.pt();
      float ymc = 1e3;
      if (mcParticle.pdgCode() == PDG_t::kK0Short)
        ymc = RecoDecay::y(std::array{mcParticle.px(), mcParticle.py(), mcParticle.pz()}, o2::constants::physics::MassKaonNeutral);
      else if (std::abs(mcParticle.pdgCode()) == PDG_t::kLambda0)
        ymc = RecoDecay::y(std::array{mcParticle.px(), mcParticle.py(), mcParticle.pz()}, o2::constants::physics::MassLambda);
      else if (std::abs(mcParticle.pdgCode()) == o2::constants::physics::kD0)
        ymc = RecoDecay::y(std::array{mcParticle.px(), mcParticle.py(), mcParticle.pz()}, o2::constants::physics::MassD0);
      else if (std::abs(mcParticle.pdgCode()) == PDG_t::kXiMinus)
        ymc = RecoDecay::y(std::array{mcParticle.px(), mcParticle.py(), mcParticle.pz()}, o2::constants::physics::MassXiMinus);
      else if (std::abs(mcParticle.pdgCode()) == PDG_t::kOmegaMinus)
        ymc = RecoDecay::y(std::array{mcParticle.px(), mcParticle.py(), mcParticle.pz()}, o2::constants::physics::MassOmegaMinus);

      if (ymc > v0Selections.rapidityMax)
        continue;
      if (ymc < v0Selections.rapidityMin)
        continue;

      float centrality = 100.5f;
      if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
        auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);
        centrality = getCentralityRun3(collision);

        if (mcParticle.pdgCode() == PDG_t::kK0Short) {
          histos.fill(HIST("h2dGenK0ShortVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0) {
          histos.fill(HIST("h2dGenLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (mcParticle.pdgCode() == PDG_t::kLambda0Bar) {
          histos.fill(HIST("h2dGenAntiLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiMinus) {
          histos.fill(HIST("h2dGenXiMinusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (mcParticle.pdgCode() == PDG_t::kXiPlusBar) {
          histos.fill(HIST("h2dGenXiPlusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaMinus) {
          histos.fill(HIST("h2dGenOmegaMinusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (mcParticle.pdgCode() == PDG_t::kOmegaPlusBar) {
          histos.fill(HIST("h2dGenOmegaPlusVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
      }

      if (mcParticle.pdgCode() == PDG_t::kK0Short) {
        histos.fill(HIST("h2dGenK0Short"), centrality, ptmc);
        histos.fill(HIST("h2dGenK0ShortVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (mcParticle.pdgCode() == PDG_t::kLambda0) {
        histos.fill(HIST("h2dGenLambda"), centrality, ptmc);
        histos.fill(HIST("h2dGenLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (mcParticle.pdgCode() == PDG_t::kLambda0Bar) {
        histos.fill(HIST("h2dGenAntiLambda"), centrality, ptmc);
        histos.fill(HIST("h2dGenAntiLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (mcParticle.pdgCode() == PDG_t::kXiMinus) {
        histos.fill(HIST("h2dGenXiMinus"), centrality, ptmc);
        histos.fill(HIST("h2dGenXiMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (mcParticle.pdgCode() == PDG_t::kXiPlusBar) {
        histos.fill(HIST("h2dGenXiPlus"), centrality, ptmc);
        histos.fill(HIST("h2dGenXiPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (mcParticle.pdgCode() == PDG_t::kOmegaMinus) {
        histos.fill(HIST("h2dGenOmegaMinus"), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaMinusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (mcParticle.pdgCode() == PDG_t::kOmegaPlusBar) {
        histos.fill(HIST("h2dGenOmegaPlus"), centrality, ptmc);
        histos.fill(HIST("h2dGenOmegaPlusVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    }
  }

  // ______________________________________________________
  // Real data processing in Run 3 - no MC subscription
  void processRealData(soa::Join<aod::Collisions, aod::EvSels, aod::MultsGlobal, aod::FT0Mults, aod::FV0Mults, aod::PVMults, aod::MultsExtra, aod::CentNGlobals, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0CVariant1s>::iterator const& collision,
                       aod::MFTTracks const& tracks,
                       soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                       aod::BCsWithTimestamps const& bcs)
  {
    analyzeRecoedV0sInRealData(collision, tracks, besttracks, bcs);
  }

  // ______________________________________________________
  // Simulated processing in Run 3 (subscribes to MC information too)
  void processMonteCarlo(soa::Join<aod::Collisions, aod::EvSels, aod::MultsGlobal, aod::FT0Mults, aod::FV0Mults, aod::PVMults, aod::MultsExtra, aod::CentNGlobals, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& tracks,
                         soa::SmallGroups<soa::Join<aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>> const& besttracks,
                         aod::BCsWithTimestamps const& bcs,
                         soa::Join<aod::McCollisions, aod::MultsExtraMC> const& /*mccollisions*/,
                         aod::McParticles const& mcParticles)
  {
    analyzeRecoedV0sInMonteCarlo(collision, tracks, besttracks, bcs, mcParticles);
  }

  // ______________________________________________________
  // Simulated processing in Run 3 (subscribes to MC information too)
  void processGenerated(soa::Join<aod::McCollisions, aod::MultsExtraMC> const& mcCollisions, aod::McParticles const& mcParticles, soa::Join<aod::Collisions, aod::EvSels, aod::MultsGlobal, aod::FT0Mults, aod::FV0Mults, aod::PVMults, aod::MultsExtra, aod::CentNGlobals, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::McCollisionLabels> const& collisions, aod::BCsWithTimestamps const& bcs)
  {
    analyzeGeneratedV0s(mcCollisions, mcParticles, collisions, bcs);
  }

  PROCESS_SWITCH(forwardlambdakzeroanalysis, processRealData, "process as if real data", true);
  PROCESS_SWITCH(forwardlambdakzeroanalysis, processMonteCarlo, "process as if MC", false);
  PROCESS_SWITCH(forwardlambdakzeroanalysis, processGenerated, "process MC generated", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<forwardlambdakzeroanalysis>(cfgc)};
}
