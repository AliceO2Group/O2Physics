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

/// \file lambda1520Analysispo.cxx
/// \brief Task for Lambda(1520) resonance reconstruction via proton-kaon invariant mass analysis for pO system
///
/// \author Durgesh Bhatt <durgesh.bhatt@cern.ch>
/// \author Narayan khadanga <narayan.khadanga@cern.ch>
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr float ProtonMass = o2::constants::physics::MassProton;
constexpr float KaonMass = o2::constants::physics::MassKaonCharged;
constexpr float LambdaMass = o2::constants::physics::MassLambda0;
constexpr float XiMass = o2::constants::physics::MassXiMinus;
constexpr float OmegaMass = o2::constants::physics::MassOmegaMinus;
constexpr float Lambda1520Mass = o2::constants::physics::MassLambda1520;
constexpr float MinAbsCharge = 1e-3f;
constexpr float InvalidCentrality = -999.0f;
const std::array<float, 2> kProtonKaonMasses{ProtonMass, KaonMass};
} // namespace

namespace lambda1520analysispo_analysis
{
static constexpr int PdgProton{2212};
static constexpr int PdgKaon{321};
static constexpr int PdgLambda0{3122};
static constexpr int PdgXiMinus{3312};
static constexpr int PdgOmegaMinus{3334};

enum CentEstimator : int {
  kFT0M = 0,
  kFT0A = 1,
  kFT0C = 2,
  kFV0A = 3
};

} // namespace lambda1520analysispo_analysis

struct Lambda1520Analysispo {

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdgDB{};

  struct : ConfigurableGroup {
    Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted |z-vertex| range [cm]"};
    Configurable<bool> cfgUseNoITSROFrameBorderCut{"cfgUseNoITSROFrameBorderCut", false, "Apply kNoITSROFrameBorder selection"};
    Configurable<bool> cfgUseSel8{"cfgUseSel8", false, "Apply sel8 selection"};
    Configurable<bool> cfgRequireRecoINELgt0{"cfgRequireRecoINELgt0", true, "Require reconstructed INEL>0 event class (applied in ALL processes)"};
    Configurable<bool> cfgUseNoSameBunchPileupCut{"cfgUseNoSameBunchPileupCut", false, "Apply kNoSameBunchPileup event selection"};
    Configurable<bool> cfgUseGoodZvtxFT0vsPVCut{"cfgUseGoodZvtxFT0vsPVCut", false, "Apply kIsGoodZvtxFT0vsPV event selection"};
    Configurable<bool> cfgUseGoodITSLayersAll{"cfgUseGoodITSLayersAll", false, "Apply kIsGoodITSLayersAll event selection"};
    Configurable<bool> applyOccupancyInTimeRangeCut{"applyOccupancyInTimeRangeCut", false, "Apply occupancy-in-time-range cut"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", 0, "Minimum track occupancy in time range"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 9999, "Maximum track occupancy in time range"};
    Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator: 0=FT0M 1=FT0A 2=FT0C 3=FV0A"};
    Configurable<bool> cfgUseMCTruthCentrality{"cfgUseMCTruthCentrality", false, "Use MC truth centrality for generated Lambda(1520) histograms (otherwise the centrality of the accepted reconstructed collision is used when there is one)"};
    Configurable<bool> cfgApplyMCTruthINELgt0{"cfgApplyMCTruthINELgt0", true, "Apply truth-level INEL>0 selection to MC generated events (event factor, signal loss, generated QA)"};
  } evSel;

  struct : ConfigurableGroup {
    Configurable<float> cfgCutPt{"cfgCutPt", 0.15f, "Minimum pT of daughter track [GeV/c]"};
    Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Maximum |eta| of daughter track"};
    Configurable<float> cfgMinY{"cfgMinY", -0.5f, "Minimum rapidity of the Lambda(1520) (reconstructed, generated and signal-loss particles)"};
    Configurable<float> cfgMaxY{"cfgMaxY", 0.5f, "Maximum rapidity of the Lambda(1520) (reconstructed, generated and signal-loss particles)"};
    Configurable<int> cfgMinTPCCrossedRows{"cfgMinTPCCrossedRows", 70, "Minimum TPC crossed pad rows"};
    Configurable<int> cfgMinTPCClusters{"cfgMinTPCClusters", 120, "Minimum TPC clusters found"};
    Configurable<int> cfgMinITSClusters{"cfgMinITSClusters", 5, "Minimum ITS clusters"};
    Configurable<int> cfgMaxTPCSharedClusters{"cfgMaxTPCSharedClusters", 0, "Maximum TPC shared clusters"};
    Configurable<float> cfgMinCrossedRowsOverFindable{"cfgMinCrossedRowsOverFindable", 0.8f, "Minimum ratio crossedRows/findableClusters"};
    Configurable<float> cfgMaxTPCChi2NCl{"cfgMaxTPCChi2NCl", 4.0f, "Maximum TPC chi2/NCl"};
    Configurable<float> cfgMaxITSChi2NCl{"cfgMaxITSChi2NCl", 36.0f, "Maximum ITS chi2/NCl"};
    Configurable<bool> applyTPCCrossedRowsCut{"applyTPCCrossedRowsCut", true, "Enable/disable minimum TPC crossed rows cut"};
    Configurable<bool> applyTPCFoundClusterCut{"applyTPCFoundClusterCut", false, "Enable/disable minimum TPC clusters found cut"};
    Configurable<bool> applyITSClusterCut{"applyITSClusterCut", false, "Enable/disable minimum ITS clusters cut"};
    Configurable<bool> applyTPCSharedClusterCut{"applyTPCSharedClusterCut", false, "Enable/disable maximum TPC shared clusters cut"};
    Configurable<bool> applyTPCCrossedRowsOverFindableCut{"applyTPCCrossedRowsOverFindableCut", false, "Enable/disable crossed rows / findable clusters cut"};
    Configurable<bool> applyTPCChi2NClCut{"applyTPCChi2NClCut", false, "Enable/disable maximum TPC chi2/NCl cut"};
    Configurable<bool> applyITSChi2NClCut{"applyITSChi2NClCut", false, "Enable/disable maximum ITS chi2/NCl cut"};
    Configurable<bool> requirePrimaryTrack{"requirePrimaryTrack", true, "Require isPrimaryTrack flag"};
    Configurable<bool> requireGlobalTrackWoDCA{"requireGlobalTrackWoDCA", true, "Require isGlobalTrackWoDCA flag"};
    Configurable<bool> requireGlobalTrack{"requireGlobalTrack", false, "Require isGlobalTrack flag"};
    Configurable<bool> requirePVContributor{"requirePVContributor", false, "Require PV-contributor flag"};
  } trackCuts;

  struct : ConfigurableGroup {
    Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "Maximum |DCAz| for all tracks [cm]"};
    Configurable<std::vector<float>> protonDCAPtBinEdges{"protonDCAPtBinEdges", {0.0f, 1000.f}, "Proton pT bin edges for DCAxy cut [GeV/c]"};
    Configurable<std::vector<float>> protonMaxDCAxyPerPtBin{"protonMaxDCAxyPerPtBin", {0.1f}, "Max |DCAxy| for proton per pT bin [cm]"};
    Configurable<std::vector<float>> kaonDCAPtBinEdges{"kaonDCAPtBinEdges", {0.0f, 1000.f}, "Kaon pT bin edges for DCAxy cut [GeV/c]"};
    Configurable<std::vector<float>> kaonMaxDCAxyPerPtBin{"kaonMaxDCAxyPerPtBin", {0.1f}, "Max |DCAxy| for kaon per pT bin [cm]"};
    Configurable<bool> isApplyPtDepDCAxyCut{"isApplyPtDepDCAxyCut", false, "Additionally apply the continuous pT-dependent DCAxy cut"};
    Configurable<float> dcaXYConst{"dcaXYConst", 0.004f, "Constant term of pT-dependent DCAxy cut [cm] (7-sigma: 0.004)"};
    Configurable<float> dcaXYPtCoeff{"dcaXYPtCoeff", 0.013f, "pT-dependent coefficient of DCAxy cut [cm] (7-sigma: 0.013)"};
    Configurable<float> dcaXYPtExponent{"dcaXYPtExponent", 1.0f, "pT exponent of DCAxy cut"};
  } dcaCuts;

  struct : ConfigurableGroup {
    Configurable<bool> useTPCOnlyPIDAtAllMomenta{"useTPCOnlyPIDAtAllMomenta", false, "Bypass TOF: TPC-only PID (|nSigma| < max TPC nSigma) at all momenta"};
    Configurable<bool> applyKaonPionRejection{"applyKaonPionRejection", false, "Reject pion-like kaon candidates (no TOF, pT window, small pion nSigma)"};
    Configurable<float> pionRejectionMinPt{"pionRejectionMinPt", 1.0f, "Pion rejection: minimum pT [GeV/c]"};
    Configurable<float> pionRejectionMaxPt{"pionRejectionMaxPt", 2.0f, "Pion rejection: maximum pT [GeV/c]"};
    Configurable<float> pionRejectionMaxTPCNSigma{"pionRejectionMaxTPCNSigma", 2.0f, "Pion rejection: |TPC nSigma pion| below this is rejected"};
  } pidShared;

  struct : ConfigurableGroup {
    Configurable<bool> requireTOFForProton{"requireTOFForProton", false, "Require TOF signal for proton candidates"};
    Configurable<float> maxTPCNSigmaProton{"maxTPCNSigmaProton", 2.0f, "Maximum |TPC nSigma| (used only with useTPCOnlyPIDAtAllMomenta)"};
    Configurable<float> protonThresholdMomentum{"protonThresholdMomentum", 0.9f, "Proton momentum p [GeV/c] for tracks WITHOUT TOF: below it TPC-only PID, at/above it the track is rejected (tracks with TOF always use the TPC+TOF cut)"};
    Configurable<std::vector<float>> protonTPCPIDMomentumBins{"protonTPCPIDMomentumBins", {0.f, 0.9f}, "momentum edges p for TPC-only proton PID [GeV/c]"};
    Configurable<std::vector<float>> protonTPCNSigmaCutPerBin{"protonTPCNSigmaCutPerBin", {2.f}, "Max |TPC nSigma| per momentum bin, TPC-only proton PID"};
    Configurable<std::vector<float>> protonCombinedPIDMomentumBins{"protonCombinedPIDMomentumBins", {0.f, 999.f}, "momentum edges p for circular TPC+TOF proton PID [GeV/c]"};
    Configurable<std::vector<float>> protonCombinedNSigmaCutPerBin{"protonCombinedNSigmaCutPerBin", {2.f}, "Max combined nSigma radius per momentum bin, circular cut"};
  } protonPID;

  struct : ConfigurableGroup {
    Configurable<bool> requireTOFForKaon{"requireTOFForKaon", false, "Require TOF signal for kaon candidates"};
    Configurable<float> maxTPCNSigmaKaon{"maxTPCNSigmaKaon", 2.0f, "Maximum |TPC nSigma| (used only with useTPCOnlyPIDAtAllMomenta)"};
    Configurable<float> kaonThresholdMomentum{"kaonThresholdMomentum", 0.5f, "Kaon momentum p [GeV/c] for tracks WITHOUT TOF: below it TPC-only PID, at/above it the track is rejected (tracks with TOF always use the TPC+TOF cut)"};
    Configurable<std::vector<float>> kaonTPCPIDMomentumBins{"kaonTPCPIDMomentumBins", {0.f, 0.5f}, "momentum edges p for TPC-only kaon PID [GeV/c]"};
    Configurable<std::vector<float>> kaonTPCNSigmaCutPerBin{"kaonTPCNSigmaCutPerBin", {2.f}, "Max |TPC nSigma| per momentum bin, TPC-only kaon PID"};
    Configurable<std::vector<float>> kaonCombinedPIDMomentumBins{"kaonCombinedPIDMomentumBins", {0.f, 999.f}, "momentum edges p for circular TPC+TOF kaon PID [GeV/c]"};
    Configurable<std::vector<float>> kaonCombinedNSigmaCutPerBin{"kaonCombinedNSigmaCutPerBin", {2.f}, "Max combined nSigma radius per momentum bin, circular cut"};
  } kaonPID;

  struct : ConfigurableGroup {
    Configurable<bool> applyDaughterDeltaEtaPhiCuts{"applyDaughterDeltaEtaPhiCuts", false, "Apply |deltaEta| and |deltaPhi| cuts between daughters"};
    Configurable<float> maxDaughterDeltaEta{"maxDaughterDeltaEta", 0.7f, "Maximum |deltaEta| between daughters"};
    Configurable<float> maxDaughterDeltaPhi{"maxDaughterDeltaPhi", 1.5f, "Maximum |deltaPhi| between daughters"};
    Configurable<bool> applyDeepAngleCut{"applyDeepAngleCut", false, "Apply minimum opening-angle cut (removes split-track background)"};
    Configurable<float> deepAngleCutValue{"deepAngleCutValue", 0.04f, "Minimum opening angle between proton and kaon [rad]"};
    Configurable<bool> applyOpeningAngleVeto{"applyOpeningAngleVeto", false, "Reject pairs whose opening angle lies INSIDE [veto min, veto max]"};
    Configurable<float> openingAngleVetoMin{"openingAngleVetoMin", 1.4f, "Opening-angle veto window: lower edge [rad]"};
    Configurable<float> openingAngleVetoMax{"openingAngleVetoMax", 2.4f, "Opening-angle veto window: upper edge [rad]"};
    Configurable<bool> applyPairPtMassUpperLimits{"applyPairPtMassUpperLimits", false, "Exclude pairs with pT or mass above the limits below (overflow bins)"};
    Configurable<float> cfgMaxPairPt{"cfgMaxPairPt", 10.0f, "Maximum pair pT [GeV/c]"};
    Configurable<float> cfgMaxPairMass{"cfgMaxPairMass", 3.0f, "Maximum pair mass [GeV/c^2]"};
  } pairCuts;

  struct : ConfigurableGroup {
    Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 10, "Number of mixed events per signal event"};
    Configurable<bool> cfgFillMixLikeSign{"cfgFillMixLikeSign", false, "Also fill like-sign pairs in mixed events"};
  } mixingCfg;

  struct : ConfigurableGroup {
    Configurable<bool> enableRotationalBackground{"enableRotationalBackground", false, "Compute rotational background (unlike-sign pairs only)"};
    Configurable<int> numberOfRotations{"numberOfRotations", 10, "Number of rotations per pair (evenly spaced, weight 1/N)"};
    Configurable<float> rotationAngleWindow{"rotationAngleWindow", 10.f, "Rotation angles within PI +/- PI/rotationAngleWindow"};
    Configurable<bool> rotateProtonInsteadOfKaon{"rotateProtonInsteadOfKaon", true, "Rotate the proton phi (false = rotate the kaon phi)"};
  } rotBkg;

  struct : ConfigurableGroup {
    Configurable<int> cfgLambda1520PdgCode{"cfgLambda1520PdgCode", static_cast<int>(o2::constants::physics::Pdg::kLambda1520_Py), "PDG code of the Lambda(1520) in the MC sample (change it for injected samples that use a different code). Sign is ignored: +code = Lambda(1520), -code = anti-Lambda(1520)"};
    Configurable<bool> applyGeneratedEtaCut{"applyGeneratedEtaCut", false, "Apply |eta| cut to generated particles and to truth-matched reconstructed pairs"};
    Configurable<float> generatedEtaMax{"generatedEtaMax", 0.5f, "MC |eta| cut"};
    Configurable<bool> applyGeneratedDaughterEtaCut{"applyGeneratedDaughterEtaCut", true, "Require both generated daughters within |eta| < trackCuts.cfgCutEta"};
    Configurable<bool> applyMtScaleWeight{"applyMtScaleWeight", false, "Weight mT-scaled signal-loss histograms by mtScaledPt/pt"};
    Configurable<bool> useGeneratedMultiplicityCentrality{"useGeneratedMultiplicityCentrality", false, "Derive MC centrality from generated multiplicity"};
    Configurable<std::vector<float>> generatedMultiplicityLowerEdges{"generatedMultiplicityLowerEdges", std::vector<float>{500, 300, 200, 120, 80, 50, 30, 10, 0}, "Generated multiplicity lower cuts"};
    Configurable<std::vector<float>> centralityPerGeneratedMultiplicityClass{"centralityPerGeneratedMultiplicityClass", std::vector<float>{2.5f, 7.5f, 15.0f, 25.0f, 35.0f, 45.0f, 55.0f, 75.0f, 95.0f}, "Centrality bin centres"};
    Configurable<float> generatedMultiplicityEtaMax{"generatedMultiplicityEtaMax", 0.5f, "Max eta for generated mid-rapidity multiplicity"};
    Configurable<float> generatedMultiplicityEtaMin{"generatedMultiplicityEtaMin", -0.5f, "Min eta for generated mid-rapidity multiplicity"};
  } mcCfg;

  struct : ConfigurableGroup {
    Configurable<bool> cfgFillTrackQA{"cfgFillTrackQA", true, "Fill track QA histograms (before/after selection)"};
    Configurable<bool> cfgFillPairQA{"cfgFillPairQA", false, "Fill pair QA histograms (opening angle, deltaEta, deltaPhi)"};
    Configurable<bool> cfgFillMultQA{"cfgFillMultQA", false, "Fill centrality vs track multiplicity correlation"};
    Configurable<bool> cfgFillMixQA{"cfgFillMixQA", false, "Fill mixed-event pool QA histograms"};
  } qa;

  struct : ConfigurableGroup {
    Configurable<int> referenceFrameDaughter{"referenceFrameDaughter", 0, "Reference daughter for cos(theta*): 0 = proton, 1 = kaon"};
    Configurable<bool> activateHelicityFrame{"activateHelicityFrame", false, "Compute cos(theta*) in the Helicity frame (quantization axis = Lambda(1520) flight direction)"};
    Configurable<bool> activateCollinsSoperFrame{"activateCollinsSoperFrame", false, "Compute cos(theta*) in the Collins-Soper frame (beam bisector)"};
    Configurable<bool> activateProductionFrame{"activateProductionFrame", false, "Compute cos(theta*) wrt the production-plane normal"};
    Configurable<bool> activateBeamAxisFrame{"activateBeamAxisFrame", false, "Compute cos(theta*) in the beam-axis (Gottfried-Jackson) frame"};
    Configurable<bool> activateRandomFrame{"activateRandomFrame", false, "Compute cos(theta*) wrt a random axis (null/control test, should be flat)"};
    Configurable<bool> fillRotatedBackgroundCosThetaStar{"fillRotatedBackgroundCosThetaStar", false, "Also fill rotated-background cos(theta*) histograms"};
    Configurable<bool> fillMixedEventCosThetaStar{"fillMixedEventCosThetaStar", false, "Also fill mixed-event-background cos(theta*) histograms"};
    Configurable<float> beamEnergy{"beamEnergy", 13600.0f, "sqrt(s) beam energy in GeV (used for Collins-Soper/beam-axis/production frames)"};
    ConfigurableAxis cosThetaBins{"cosThetaBins", {20, -1.0f, 1.0f}, "cos(theta*) binning"};
  } spinAlignment;

  struct : ConfigurableGroup {
    ConfigurableAxis cfgPtAxis{"cfgPtAxis", {VARIABLE_WIDTH, 0.0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.75, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 4.6, 4.8, 4.9, 5.0, 5.5, 5.6, 6.0, 6.4, 6.5, 7.0, 7.2, 8.0, 9.0, 9.5, 9.6, 10.0, 11.0, 11.5, 12.0, 13.0, 14.0, 14.4, 15.0, 16.0, 18.0, 19.2, 20.0}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis cfgMassAxis{"cfgMassAxis", {70, 1.3f, 2.5f}, "Invariant mass (GeV/#it{c}^{2})"};
    ConfigurableAxis cfgCentAxis{"cfgCentAxis", {VARIABLE_WIDTH, 0.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f}, "Centrality (%)"};
    ConfigurableAxis cfgRapAxis{"cfgRapAxis", {20, -1.0, 1.0}, "Rapidity y"};
    ConfigurableAxis cfgVtxAxis{"cfgVtxAxis", {VARIABLE_WIDTH, -10.0f, -8.0f, -6.0f, -4.0f, -2.0f, 0.0f, 2.0f, 4.0f, 6.0f, 8.0f, 10.0f}, "Vertex z [cm] (histograms and event-mixing bins)"};
  } axes;

  struct DaughterTrackCandidate {
    float px = 0.f, py = 0.f, pz = 0.f;
    float pt = 0.f, eta = 0.f, phi = 0.f;
    float mom = 0.f;
    float dcaXY = 0.f, dcaZ = 0.f;
    int sign = 0;
    bool hasTOF = false;
    float tpcNSigmaPr = 0.f, tpcNSigmaKa = 0.f;
    float tofNSigmaPr = -999.f, tofNSigmaKa = -999.f;
    int64_t globalIndex = 0;
    bool hasMcLabel = false;
    int mcPdgCode = 0;
    int64_t mcLambda1520MotherIndex = -1;
    int mcLambda1520MotherPdgCode = 0;
  };

  struct SpeciesPIDSettings {
    std::vector<float> tpcOnlyBins, tpcOnlyCuts;
    std::vector<float> combinedBins, combinedCuts;
    float maxTPCNSigma = 2.f;
    float thresholdMomentum = 0.f;
  };

  struct DecayKinematics {
    std::array<float, 3> motherMom{};
    float motherE = 0.f;
    std::array<float, 3> protonMom{};
    float protonE = 0.f;
    std::array<float, 3> kaonMom{};
    float kaonE = 0.f;
  };

  struct AcceptedRecoCollision {
    bool found = false;
    float centrality = 0.f;
  };

  struct NoCutFlowTag {
  };
  struct DataCutFlowTag {
  };
  struct McRecoCutFlowTag {
  };

  enum class CosThetaStarSample {
    Signal,
    LikeSign,
    Rotated,
    Mixed,
    LikeSignMixed,
    TruthMatched
  };

  SpeciesPIDSettings mProtonPIDCuts;
  SpeciesPIDSettings mKaonPIDCuts;
  std::vector<float> mProtonDCAPtEdges;
  std::vector<float> mProtonMaxDCAxy;
  std::vector<float> mKaonDCAPtEdges;
  std::vector<float> mKaonMaxDCAxy;
  int mLambda1520AbsPdgCode = 0;
  std::vector<DaughterTrackCandidate> mProtonCandidates;
  std::vector<DaughterTrackCandidate> mKaonCandidates;

  float mBeam1Energy = 0.f;
  std::array<float, 3> mBeam1Momentum{};
  float mBeam2Energy = 0.f;
  std::array<float, 3> mBeam2Momentum{};
  std::mt19937 mRandomEngine{std::random_device{}()};
  std::uniform_real_distribution<double> mUniform01{0., 1.};

  template <typename HistPtr>
  static void setHistogramBinLabels(HistPtr const& hist, std::vector<const char*> const& labels)
  {
    for (size_t i = 0; i < labels.size(); ++i) {
      hist->GetXaxis()->SetBinLabel(static_cast<int>(i) + 1, labels[i]);
    }
  }

  bool isSpinAlignmentActive() const
  {
    return spinAlignment.activateHelicityFrame.value ||
           spinAlignment.activateCollinsSoperFrame.value ||
           spinAlignment.activateProductionFrame.value ||
           spinAlignment.activateBeamAxisFrame.value ||
           spinAlignment.activateRandomFrame.value;
  }

  void validateConfiguration()
  {
    auto checkEdgesMatchCuts = [](auto const& cuts, auto const& bins, const char* label) {
      const auto& cutsVec = cuts.value;
      const auto& binsVec = bins.value;
      if (cutsVec.size() != binsVec.size() - 1) {
        LOG(fatal) << "PID bin-size mismatch for: " << label
                   << " (cuts=" << cutsVec.size()
                   << ", bins=" << binsVec.size() << ")";
      }
    };
    checkEdgesMatchCuts(protonPID.protonTPCNSigmaCutPerBin, protonPID.protonTPCPIDMomentumBins, "proton TPC-only PID");
    checkEdgesMatchCuts(protonPID.protonCombinedNSigmaCutPerBin, protonPID.protonCombinedPIDMomentumBins, "proton combined PID");
    checkEdgesMatchCuts(kaonPID.kaonTPCNSigmaCutPerBin, kaonPID.kaonTPCPIDMomentumBins, "kaon TPC-only PID");
    checkEdgesMatchCuts(kaonPID.kaonCombinedNSigmaCutPerBin, kaonPID.kaonCombinedPIDMomentumBins, "kaon combined PID");
    checkEdgesMatchCuts(dcaCuts.protonMaxDCAxyPerPtBin, dcaCuts.protonDCAPtBinEdges, "proton DCAxy pT bins");
    checkEdgesMatchCuts(dcaCuts.kaonMaxDCAxyPerPtBin, dcaCuts.kaonDCAPtBinEdges, "kaon DCAxy pT bins");

    if (mcCfg.useGeneratedMultiplicityCentrality &&
        mcCfg.generatedMultiplicityLowerEdges.value.size() != mcCfg.centralityPerGeneratedMultiplicityClass.value.size()) {
      LOG(fatal) << "generatedMultiplicityLowerEdges and centralityPerGeneratedMultiplicityClass must have the same size";
    }

    const int nActiveFrames = static_cast<int>(spinAlignment.activateHelicityFrame) +
                              static_cast<int>(spinAlignment.activateCollinsSoperFrame) +
                              static_cast<int>(spinAlignment.activateProductionFrame) +
                              static_cast<int>(spinAlignment.activateBeamAxisFrame) +
                              static_cast<int>(spinAlignment.activateRandomFrame);
    if (nActiveFrames > 1) {
      LOG(fatal) << "Multiple spin-alignment frames activated simultaneously. "
                 << "Please activate only ONE frame at a time to prevent data mixing in the cos(theta*) THnSparse.";
    }

    if (doprocessMCGenSpinAlignment && !isSpinAlignmentActive()) {
      LOG(fatal) << "processMCGenSpinAlignment requires at least one spin-alignment frame to be active";
    }
  }

  void cacheCutSettings()
  {
    mLambda1520AbsPdgCode = std::abs(mcCfg.cfgLambda1520PdgCode.value);
    if (mLambda1520AbsPdgCode == 0) {
      LOG(fatal) << "cfgLambda1520PdgCode must be non-zero";
    }
    if (doprocessMC || doprocessMCGen || doprocessSignalLoss) {
      LOG(info) << "Lambda(1520) PDG code used for MC truth: " << mLambda1520AbsPdgCode;
    }

    mProtonDCAPtEdges = dcaCuts.protonDCAPtBinEdges.value;
    mProtonMaxDCAxy = dcaCuts.protonMaxDCAxyPerPtBin.value;
    mKaonDCAPtEdges = dcaCuts.kaonDCAPtBinEdges.value;
    mKaonMaxDCAxy = dcaCuts.kaonMaxDCAxyPerPtBin.value;

    mProtonPIDCuts.tpcOnlyBins = protonPID.protonTPCPIDMomentumBins.value;
    mProtonPIDCuts.tpcOnlyCuts = protonPID.protonTPCNSigmaCutPerBin.value;
    mProtonPIDCuts.combinedBins = protonPID.protonCombinedPIDMomentumBins.value;
    mProtonPIDCuts.combinedCuts = protonPID.protonCombinedNSigmaCutPerBin.value;
    mProtonPIDCuts.maxTPCNSigma = protonPID.maxTPCNSigmaProton;
    mProtonPIDCuts.thresholdMomentum = protonPID.protonThresholdMomentum;

    mKaonPIDCuts.tpcOnlyBins = kaonPID.kaonTPCPIDMomentumBins.value;
    mKaonPIDCuts.tpcOnlyCuts = kaonPID.kaonTPCNSigmaCutPerBin.value;
    mKaonPIDCuts.combinedBins = kaonPID.kaonCombinedPIDMomentumBins.value;
    mKaonPIDCuts.combinedCuts = kaonPID.kaonCombinedNSigmaCutPerBin.value;
    mKaonPIDCuts.maxTPCNSigma = kaonPID.maxTPCNSigmaKaon;
    mKaonPIDCuts.thresholdMomentum = kaonPID.kaonThresholdMomentum;

    auto checkTPCOnlyMomentumCoverage = [](std::vector<float> const& edges, float threshold, const char* label) {
      if (threshold > 0.f && (edges.empty() || edges.front() > 0.f || edges.back() < threshold)) {
        LOG(warning) << label << ": TPC-only momentum bins do not cover [0, " << threshold
                     << ") GeV/c - tracks WITHOUT TOF in the uncovered range are rejected";
      }
    };
    if (!pidShared.useTPCOnlyPIDAtAllMomenta) {
      checkTPCOnlyMomentumCoverage(mProtonPIDCuts.tpcOnlyBins, mProtonPIDCuts.thresholdMomentum, "proton PID");
      checkTPCOnlyMomentumCoverage(mKaonPIDCuts.tpcOnlyBins, mKaonPIDCuts.thresholdMomentum, "kaon PID");
    }
  }

  void computeBeamKinematics()
  {
    const auto sqrtS = static_cast<double>(spinAlignment.beamEnergy);
    const double beamMomentum = std::sqrt(sqrtS * sqrtS / 4. - static_cast<double>(ProtonMass) * static_cast<double>(ProtonMass));
    mBeam1Energy = static_cast<float>(sqrtS * 0.5);
    mBeam1Momentum = {0.f, 0.f, static_cast<float>(-beamMomentum)};
    mBeam2Energy = static_cast<float>(sqrtS * 0.5);
    mBeam2Momentum = {0.f, 0.f, static_cast<float>(beamMomentum)};
  }

  void init(InitContext const&)
  {
    validateConfiguration();
    cacheCutSettings();
    computeBeamKinematics();

    const AxisSpec ptAxis{axes.cfgPtAxis, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec massAxis{axes.cfgMassAxis, "M_{inv} (GeV/#it{c}^{2})"};
    const AxisSpec centAxis{axes.cfgCentAxis, "Centrality (%)"};
    const AxisSpec rapAxis{axes.cfgRapAxis, "Rapidity y"};
    const AxisSpec vtxAxis{axes.cfgVtxAxis, "Vertex z [cm]"};
    const AxisSpec cosThetaAxis{spinAlignment.cosThetaBins, "cos#theta*"};
    const AxisSpec nSigmaTPCaxis{100, -5., 5., "n#sigma^{TPC}"};
    const AxisSpec nSigmaTOFaxis{100, -5., 5., "n#sigma^{TOF}"};
    const AxisSpec momentumAxis{100, 0., 10., "p (GeV/#it{c})"};
    const AxisSpec ptForPIDAxis{100, 0., 10., "p_{T} (GeV/#it{c})"};
    const AxisSpec dcaXYaxis{200, -0.2f, 0.2f, "DCA_{xy} (cm)"};
    const AxisSpec dcaZaxis{200, -1.2f, 1.2f, "DCA_{z} (cm)"};
    const AxisSpec tpcClusAxis{200, 0, 200, "N_{clusters}^{TPC}"};
    const AxisSpec tpcRowsAxis{200, 0, 200, "Crossed rows^{TPC}"};
    const AxisSpec tpcSignalAxis{500, 0, 500, "TPC signal (a.u.)"};
    const AxisSpec occupancyAxis{200, 0, 10000, "Track occupancy [-40, 100] ns"};
    const AxisSpec etaAxis{150, -1.5, 1.5, "#eta"};
    const AxisSpec phiAxis{72, 0., o2::constants::math::TwoPI, "#varphi (rad)"};
    const AxisSpec openAngleAxis{100, 0., o2::constants::math::PI, "Opening angle (rad)"};
    const AxisSpec deltaEtaAxis{160, 0., 1.6, "|#Delta#eta|"};
    const AxisSpec deltaPhiAxis{100, 0., o2::constants::math::PI, "|#Delta#varphi| (rad)"};
    const AxisSpec massResolutionAxis{200, -0.05f, 0.05f, "M_{reco} - M_{gen} (GeV/c^{2})"};
    const std::vector<const char*> eventCutLabels{"All", "Vz", "kNoITSROFrameBorder", "sel8", "INEL>0", "Occupancy", "NoSameBunchPileup", "GoodZvtxFT0vsPV", "GoodITSLayersAll", "Accepted"};

    if (doprocessSameEvent) {
      histos.add("CutFlow/Data/hEventSelection", "Data event cut flow", kTH1F, {{10, -0.5f, 9.5f}});
      setHistogramBinLabels(histos.get<TH1>(HIST("CutFlow/Data/hEventSelection")), eventCutLabels);

      histos.add("EventQA/Data/hVertexZ", "Vertex z; z (cm)", kTH1F, {{400, -20., 20.}});
      histos.add("EventQA/Data/hPVContributors", "PV contributors; N", kTH1F, {{2001, -0.5f, 2000.5f}});
      histos.add("EventQA/Data/hCentrality", "Centrality", kTH1F, {centAxis});
      histos.add("EventQA/Data/hOccupancy", "Occupancy in time range", kTH1F, {occupancyAxis});
      histos.add("EventQA/Data/hCentralityVsVertexZ", "Centrality vs vertex z", kTH2F, {vtxAxis, centAxis});
      histos.add("EventQA/Data/hCentralityVsOccupancy", "Centrality vs occupancy", kTH2F, {occupancyAxis, centAxis});

      if (qa.cfgFillMultQA) {
        histos.add("EventQA/Data/hCentralityVsFilteredTracksVsPVTracks", "Centrality vs filtered tracks vs PV tracks", kTHnSparseF, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Filtered tracks"}, {500, 0, 5000, "PV tracks"}});
      }

      if (qa.cfgFillTrackQA) {
        histos.add("TrackQA/Data/PreSelection/Proton/hNSigmaTPCvsP", "TPC n#sigma proton vs p (before PID)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/PreSelection/Proton/hNSigmaTOFvsP", "TOF n#sigma proton vs p (before PID)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/PreSelection/Proton/hNSigmaTOFvsNSigmaTPC", "TOF vs TPC n#sigma proton (before PID)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/PreSelection/Kaon/hNSigmaTPCvsP", "TPC n#sigma kaon vs p (before PID)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/PreSelection/Kaon/hNSigmaTOFvsP", "TOF n#sigma kaon vs p (before PID)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/PreSelection/Kaon/hNSigmaTOFvsNSigmaTPC", "TOF vs TPC n#sigma kaon (before PID)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});

        histos.add("TrackQA/Data/AllTracks/hEta", "Reco dN/d#eta; #eta; Tracks", kTH1F, {{50, -1.0, 1.0}});
        histos.add("TrackQA/Data/AllTracks/hPt", "Reco pT; p_{T} (GeV/c); Tracks", kTH1F, {ptForPIDAxis});
        histos.add("TrackQA/Data/AllTracks/hPhi", "Reco #varphi; #varphi (rad); Tracks", kTH1F, {phiAxis});
        histos.add("TrackQA/Data/AllTracks/hEtaPhi", "Reco #eta vs #varphi; #eta; #varphi", kTH2F, {etaAxis, phiAxis});

        histos.add("TrackQA/Data/Selected/Proton/hDCAxyVsPt", "Proton DCA_{xy} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaXYaxis});
        histos.add("TrackQA/Data/Selected/Proton/hDCAzVsPt", "Proton DCA_{z} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaZaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTPCvsP", "Proton TPC n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTPCvsPt", "Proton TPC n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTPCvsCentrality", "Proton TPC n#sigma vs centrality", kTH2F, {centAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTOFvsP", "Proton TOF n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTOFvsPt", "Proton TOF n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTOFvsCentrality", "Proton TOF n#sigma vs centrality", kTH2F, {centAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTOFvsNSigmaTPC", "Proton TOF vs TPC n#sigma (after cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTPCPionVsPt", "Proton track: TPC n#sigma pion contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Proton/hNSigmaTPCKaonVsPt", "Proton track: TPC n#sigma kaon contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Proton/hTPCSignalVsInnerP", "Proton TPC signal vs p (after cuts)", kTH2F, {momentumAxis, tpcSignalAxis});
        histos.add("TrackQA/Data/Selected/Proton/hTPCClustersVsPhi", "Proton TPC clusters vs #varphi", kTH2F, {tpcClusAxis, phiAxis});
        histos.add("TrackQA/Data/Selected/Proton/hTPCCrossedRowsVsPt", "Proton TPC crossed rows vs p_{T}", kTH2F, {ptForPIDAxis, tpcRowsAxis});
        histos.add("TrackQA/Data/Selected/Proton/hTPCFoundClustersVsPt", "Proton TPC clusters found vs p_{T}", kTH2F, {ptForPIDAxis, tpcClusAxis});
        histos.add("TrackQA/Data/Selected/Proton/hPt", "Proton p_{T}; p_{T} (GeV/c); Tracks", kTH1F, {ptForPIDAxis});
        histos.add("TrackQA/Data/Selected/Proton/hEta", "Proton #eta; #eta; Tracks", kTH1F, {etaAxis});

        histos.add("TrackQA/Data/Selected/Kaon/hDCAxyVsPt", "Kaon DCA_{xy} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaXYaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hDCAzVsPt", "Kaon DCA_{z} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaZaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTPCvsP", "Kaon TPC n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTPCvsPt", "Kaon TPC n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTPCvsCentrality", "Kaon TPC n#sigma vs centrality", kTH2F, {centAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsP", "Kaon TOF n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsPt", "Kaon TOF n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsCentrality", "Kaon TOF n#sigma vs centrality", kTH2F, {centAxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsNSigmaTPC", "Kaon TOF vs TPC n#sigma (after cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTPCPionVsPt", "Kaon track: TPC n#sigma pion contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hNSigmaTPCProtonVsPt", "Kaon track: TPC n#sigma proton contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/Data/Selected/Kaon/hTPCSignalVsInnerP", "Kaon TPC signal vs p (after cuts)", kTH2F, {momentumAxis, tpcSignalAxis});
        histos.add("TrackQA/Data/Selected/Kaon/hTPCClustersVsPhi", "Kaon TPC clusters vs #varphi", kTH2F, {tpcClusAxis, phiAxis});
        histos.add("TrackQA/Data/Selected/Kaon/hTPCCrossedRowsVsPt", "Kaon TPC crossed rows vs p_{T}", kTH2F, {ptForPIDAxis, tpcRowsAxis});
        histos.add("TrackQA/Data/Selected/Kaon/hTPCFoundClustersVsPt", "Kaon TPC clusters found vs p_{T}", kTH2F, {ptForPIDAxis, tpcClusAxis});
        histos.add("TrackQA/Data/Selected/Kaon/hPt", "Kaon p_{T}; p_{T} (GeV/c); Tracks", kTH1F, {ptForPIDAxis});
        histos.add("TrackQA/Data/Selected/Kaon/hEta", "Kaon #eta; #eta; Tracks", kTH1F, {etaAxis});
      }

      if (qa.cfgFillPairQA) {
        histos.add("PairQA/Data/hOpeningAngleBeforeCuts", "Opening angle (before pair cuts)", kTH1F, {openAngleAxis});
        histos.add("PairQA/Data/hOpeningAngleAfterCuts", "Opening angle (after pair cuts)", kTH1F, {openAngleAxis});
        histos.add("PairQA/Data/hDeltaEtaBeforeCuts", "|#Delta#eta| (before pair cuts)", kTH1F, {deltaEtaAxis});
        histos.add("PairQA/Data/hDeltaEtaAfterCuts", "|#Delta#eta| (after pair cuts)", kTH1F, {deltaEtaAxis});
        histos.add("PairQA/Data/hDeltaPhiBeforeCuts", "|#Delta#varphi| (before pair cuts)", kTH1F, {deltaPhiAxis});
        histos.add("PairQA/Data/hDeltaPhiAfterCuts", "|#Delta#varphi| (after pair cuts)", kTH1F, {deltaPhiAxis});
      }

      if (!isSpinAlignmentActive()) {
        histos.add("InvMass/Data/SameEvent/hLambda1520", "THnSparse #Lambda(1520) same-event", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
        histos.add("InvMass/Data/SameEvent/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) same-event", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
        histos.add("InvMass/Data/SameEvent/hLikeSignProtonPlusKaonPlus", "THnSparse like-sign p K^{+}", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
        histos.add("InvMass/Data/SameEvent/hLikeSignProtonMinusKaonMinus", "THnSparse like-sign #bar{p} K^{-}", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});

        if (rotBkg.enableRotationalBackground) {
          histos.add("InvMass/Data/RotatedBackground/hLambda1520", "THnSparse #Lambda(1520) rotational background", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
          histos.add("InvMass/Data/RotatedBackground/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) rotational background", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
          histos.get<THnSparse>(HIST("InvMass/Data/RotatedBackground/hLambda1520"))->Sumw2();
          histos.get<THnSparse>(HIST("InvMass/Data/RotatedBackground/hAntiLambda1520"))->Sumw2();
        }
      }

      if (isSpinAlignmentActive()) {
        histos.add("CosThetaStar/Data/SameEvent/hLambda1520", "cos#theta* #Lambda(1520) unlike-sign", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/Data/SameEvent/hAntiLambda1520", "cos#theta* #bar{#Lambda}(1520) unlike-sign", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/Data/SameEvent/hLikeSignProtonPlusKaonPlus", "cos#theta* like-sign p K^{+}", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/Data/SameEvent/hLikeSignProtonMinusKaonMinus", "cos#theta* like-sign #bar{p} K^{-}", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});

        if (spinAlignment.fillRotatedBackgroundCosThetaStar && rotBkg.enableRotationalBackground) {
          histos.add("CosThetaStar/Data/RotatedBackground/hLambda1520", "cos#theta* rotated bkg #Lambda(1520)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
          histos.add("CosThetaStar/Data/RotatedBackground/hAntiLambda1520", "cos#theta* rotated bkg #bar{#Lambda}(1520)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        }
      }
    }

    if (doprocessMixedEvent) {
      if (!isSpinAlignmentActive()) {
        histos.add("InvMass/Data/MixedEvent/hLambda1520", "THnSparse #Lambda(1520) event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
        histos.add("InvMass/Data/MixedEvent/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});

        if (mixingCfg.cfgFillMixLikeSign) {
          histos.add("InvMass/Data/MixedEvent/hLikeSignProtonPlusKaonPlus", "THnSparse like-sign p K^{+} event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
          histos.add("InvMass/Data/MixedEvent/hLikeSignProtonMinusKaonMinus", "THnSparse like-sign #bar{p} K^{-} event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
        }

        if (qa.cfgFillMixQA) {
          histos.add("MixingQA/Data/hPoolVertexZ", "Mixed event pool: vertex z; z (cm); Pairs", kTH1F, {vtxAxis});
          histos.add("MixingQA/Data/hPoolCentrality", "Mixed event pool: centrality; Centrality (%); Pairs", kTH1F, {centAxis});
          histos.add("MixingQA/Data/hPoolVertexZVsCentrality", "Mixed event pool: vertex z vs centrality", kTH2F, {vtxAxis, centAxis});
          histos.add("MixingQA/Data/hCollisionPairs", "Number of mixed-event collision pairs processed", kTH1F, {{1, 0.5f, 1.5f}});
        }
      }

      if (isSpinAlignmentActive() && spinAlignment.fillMixedEventCosThetaStar) {
        histos.add("CosThetaStar/Data/MixedEvent/hLambda1520", "cos#theta* mixed-event #Lambda(1520)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/Data/MixedEvent/hAntiLambda1520", "cos#theta* mixed-event #bar{#Lambda}(1520)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/Data/MixedEvent/hLikeSignProtonPlusKaonPlus", "cos#theta* mixed-event like-sign p K^{+}", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/Data/MixedEvent/hLikeSignProtonMinusKaonMinus", "cos#theta* mixed-event like-sign #bar{p} K^{-}", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
      }
    }

    if (doprocessMC) {
      histos.add("CutFlow/MCReco/hEventSelection", "MC reconstructed event cut flow", kTH1F, {{10, -0.5f, 9.5f}});
      setHistogramBinLabels(histos.get<TH1>(HIST("CutFlow/MCReco/hEventSelection")), eventCutLabels);

      histos.add("EventQA/MCReco/hAcceptedEventsVsCentrality", "Reconstructed INEL>0 events (Nrec, MC reco)", kTH1F, {centAxis});
      histos.add("EventQA/MCReco/hVertexZ", "MC reco vertex z; z (cm)", kTH1F, {{400, -20., 20.}});
      histos.add("EventQA/MCReco/hPVContributors", "MC reco PV contributors; N", kTH1F, {{2001, -0.5f, 2000.5f}});

      histos.add("InvMass/MCReco/AllPairs/hLambda1520", "THnSparse #Lambda(1520) MC reconstructed (pre-truth-matching)", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("InvMass/MCReco/AllPairs/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) MC reconstructed (pre-truth-matching)", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});

      histos.add("InvMass/MCReco/TruthMatched/hLambda1520", "THnSparse #Lambda(1520) MC reconstructed (truth-matched)", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("InvMass/MCReco/TruthMatched/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) MC reconstructed (truth-matched)", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});

      histos.add("MCRecoQA/DaughterPt/hProtonFromLambda1520", "Rec proton from #Lambda(1520)", kTH1F, {ptAxis});
      histos.add("MCRecoQA/DaughterPt/hProtonFromAntiLambda1520", "Rec proton from #bar{#Lambda}(1520)", kTH1F, {ptAxis});
      histos.add("MCRecoQA/DaughterPt/hKaonFromLambda1520", "Rec kaon from #Lambda(1520)", kTH1F, {ptAxis});
      histos.add("MCRecoQA/DaughterPt/hKaonFromAntiLambda1520", "Rec kaon from #bar{#Lambda}(1520)", kTH1F, {ptAxis});

      histos.add("TrackQA/MCReco/AllTracks/hEta", "MC Reco dN/d#eta; #eta; Tracks", kTH1F, {{50, -1.0, 1.0}});
      histos.add("TrackQA/MCReco/AllTracks/hPt", "MC Reco pT; p_{T} (GeV/c); Tracks", kTH1F, {ptForPIDAxis});
      histos.add("TrackQA/MCReco/AllTracks/hPhi", "MC Reco #varphi; #varphi (rad); Tracks", kTH1F, {phiAxis});
      histos.add("TrackQA/MCReco/AllTracks/hEtaPhi", "MC Reco #eta vs #varphi; #eta; #varphi", kTH2F, {etaAxis, phiAxis});
      histos.add("PairQA/MCReco/hTruthMatchedRapidityVsPt", "Truth-matched reconstructed #Lambda(1520): y vs p_{T}", kTH2F, {ptAxis, rapAxis});

      if (qa.cfgFillTrackQA) {
        histos.add("TrackQA/MCReco/PreSelection/Proton/hNSigmaTPCvsP", "MC reco: TPC n#sigma proton vs p (before PID)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/PreSelection/Proton/hNSigmaTOFvsP", "MC reco: TOF n#sigma proton vs p (before PID)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/PreSelection/Proton/hNSigmaTOFvsNSigmaTPC", "MC reco: TOF vs TPC n#sigma proton (before PID)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/PreSelection/Kaon/hNSigmaTPCvsP", "MC reco: TPC n#sigma kaon vs p (before PID)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/PreSelection/Kaon/hNSigmaTOFvsP", "MC reco: TOF n#sigma kaon vs p (before PID)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/PreSelection/Kaon/hNSigmaTOFvsNSigmaTPC", "MC reco: TOF vs TPC n#sigma kaon (before PID)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});

        histos.add("TrackQA/MCReco/Selected/Proton/hDCAxyVsPt", "MC reco: Proton DCA_{xy} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaXYaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hDCAzVsPt", "MC reco: Proton DCA_{z} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaZaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTPCvsP", "MC reco: Proton TPC n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTPCvsPt", "MC reco: Proton TPC n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTPCvsCentrality", "MC reco: Proton TPC n#sigma vs centrality", kTH2F, {centAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsP", "MC reco: Proton TOF n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsPt", "MC reco: Proton TOF n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsCentrality", "MC reco: Proton TOF n#sigma vs centrality", kTH2F, {centAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsNSigmaTPC", "MC reco: Proton TOF vs TPC n#sigma (after cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTPCPionVsPt", "MC reco: Proton track: TPC n#sigma pion contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hNSigmaTPCKaonVsPt", "MC reco: Proton track: TPC n#sigma kaon contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hTPCSignalVsInnerP", "MC reco: Proton TPC signal vs p (after cuts)", kTH2F, {momentumAxis, tpcSignalAxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hTPCClustersVsPhi", "MC reco: Proton TPC clusters vs #varphi", kTH2F, {tpcClusAxis, phiAxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hTPCCrossedRowsVsPt", "MC reco: Proton TPC crossed rows vs p_{T}", kTH2F, {ptForPIDAxis, tpcRowsAxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hTPCFoundClustersVsPt", "MC reco: Proton TPC clusters found vs p_{T}", kTH2F, {ptForPIDAxis, tpcClusAxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hPt", "MC reco: Proton p_{T}; p_{T} (GeV/c); Tracks", kTH1F, {ptForPIDAxis});
        histos.add("TrackQA/MCReco/Selected/Proton/hEta", "MC reco: Proton #eta; #eta; Tracks", kTH1F, {etaAxis});

        histos.add("TrackQA/MCReco/Selected/Kaon/hDCAxyVsPt", "MC reco: Kaon DCA_{xy} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaXYaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hDCAzVsPt", "MC reco: Kaon DCA_{z} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaZaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCvsP", "MC reco: Kaon TPC n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCvsPt", "MC reco: Kaon TPC n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCvsCentrality", "MC reco: Kaon TPC n#sigma vs centrality", kTH2F, {centAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsP", "MC reco: Kaon TOF n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsPt", "MC reco: Kaon TOF n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsCentrality", "MC reco: Kaon TOF n#sigma vs centrality", kTH2F, {centAxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsNSigmaTPC", "MC reco: Kaon TOF vs TPC n#sigma (after cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCPionVsPt", "MC reco: Kaon track: TPC n#sigma pion contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCProtonVsPt", "MC reco: Kaon track: TPC n#sigma proton contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hTPCSignalVsInnerP", "MC reco: Kaon TPC signal vs p (after cuts)", kTH2F, {momentumAxis, tpcSignalAxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hTPCClustersVsPhi", "MC reco: Kaon TPC clusters vs #varphi", kTH2F, {tpcClusAxis, phiAxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hTPCCrossedRowsVsPt", "MC reco: Kaon TPC crossed rows vs p_{T}", kTH2F, {ptForPIDAxis, tpcRowsAxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hTPCFoundClustersVsPt", "MC reco: Kaon TPC clusters found vs p_{T}", kTH2F, {ptForPIDAxis, tpcClusAxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hPt", "MC reco: Kaon p_{T}; p_{T} (GeV/c); Tracks", kTH1F, {ptForPIDAxis});
        histos.add("TrackQA/MCReco/Selected/Kaon/hEta", "MC reco: Kaon #eta; #eta; Tracks", kTH1F, {etaAxis});
      }

      histos.add("TrackQA/MCReco/TruthMatched/Proton/hDCAxyVsPt", "Truth-matched proton DCA_{xy} vs p_{T}", kTH2F, {ptForPIDAxis, dcaXYaxis});
      histos.add("TrackQA/MCReco/TruthMatched/Proton/hDCAzVsPt", "Truth-matched proton DCA_{z} vs p_{T}", kTH2F, {ptForPIDAxis, dcaZaxis});
      histos.add("TrackQA/MCReco/TruthMatched/Proton/hNSigmaTPCvsPt", "Truth-matched proton TPC n#sigma vs p_{T}", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
      histos.add("TrackQA/MCReco/TruthMatched/Proton/hNSigmaTOFvsPt", "Truth-matched proton TOF n#sigma vs p_{T}", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});
      histos.add("TrackQA/MCReco/TruthMatched/Kaon/hDCAxyVsPt", "Truth-matched kaon DCA_{xy} vs p_{T}", kTH2F, {ptForPIDAxis, dcaXYaxis});
      histos.add("TrackQA/MCReco/TruthMatched/Kaon/hDCAzVsPt", "Truth-matched kaon DCA_{z} vs p_{T}", kTH2F, {ptForPIDAxis, dcaZaxis});
      histos.add("TrackQA/MCReco/TruthMatched/Kaon/hNSigmaTPCvsPt", "Truth-matched kaon TPC n#sigma vs p_{T}", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
      histos.add("TrackQA/MCReco/TruthMatched/Kaon/hNSigmaTOFvsPt", "Truth-matched kaon TOF n#sigma vs p_{T}", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});

      histos.add("CutFlow/MCReco/hGeneratedToReconstructed", "Generated-to-reconstructed #Lambda(1520) cut flow (efficiency QA)", kTH1F, {{19, -0.5f, 18.5f}});
      setHistogramBinLabels(histos.get<TH1>(HIST("CutFlow/MCReco/hGeneratedToReconstructed")),
                            {"Generated Lambda(1520) (event accepted)", "Lambda(1520) -> p + K", "Rapidity accepted",
                             "Proton acceptance", "Kaon acceptance", "Both daughters acceptance",
                             "Proton reconstructed", "Kaon reconstructed", "Both daughters reconstructed",
                             "Proton track quality", "Kaon track quality", "Both track quality",
                             "Proton PID", "Kaon PID", "Both PID",
                             "Proton DCA", "Kaon DCA", "Both DCA",
                             "Pair cuts (final reconstructed)"});

      histos.add("CutFlow/MCReco/hTruthMatching", "MC truth-matching cut flow (efficiency QA)", kTH1F, {{6, -0.5f, 5.5f}});
      setHistogramBinLabels(histos.get<TH1>(HIST("CutFlow/MCReco/hTruthMatching")),
                            {"Reco pair", "Proton has MC label", "Kaon has MC label", "Proton is proton", "Kaon is kaon", "Common Lambda(1520) mother"});

      histos.add("InvMass/MCReco/MassResolution/hLambda1520", "(M_{reco} - M_{gen}) truth-matched #Lambda(1520)", kTHnSparseF, {massResolutionAxis, ptAxis, centAxis});
      histos.add("InvMass/MCReco/MassResolution/hAntiLambda1520", "(M_{reco} - M_{gen}) truth-matched #bar{#Lambda}(1520)", kTHnSparseF, {massResolutionAxis, ptAxis, centAxis});

      if (isSpinAlignmentActive()) {
        histos.add("CosThetaStar/MCReco/TruthMatched/hLambda1520", "cos#theta* truth-matched #Lambda(1520)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/MCReco/TruthMatched/hAntiLambda1520", "cos#theta* truth-matched #bar{#Lambda}(1520)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/MCReco/MassResolution/hLambda1520", "(M_{reco} - M_{gen}) vs cos#theta* #Lambda(1520)", kTHnSparseF, {massResolutionAxis, ptAxis, centAxis, cosThetaAxis});
        histos.add("CosThetaStar/MCReco/MassResolution/hAntiLambda1520", "(M_{reco} - M_{gen}) vs cos#theta* #bar{#Lambda}(1520)", kTHnSparseF, {massResolutionAxis, ptAxis, centAxis, cosThetaAxis});
      }
    }

    if (doprocessMCEventMixing) {
      histos.add("InvMass/MCReco/MixedEvent/hLambda1520", "THnSparse #Lambda(1520) MC reconstructed event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("InvMass/MCReco/MixedEvent/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) MC reconstructed event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("MixingQA/MCReco/hCollisionPairs", "Number of MC mixed-event collision pairs processed", kTH1F, {{1, 0.5f, 1.5f}});
      histos.add("MixingQA/MCReco/hCollisionPairsVsCentrality", "MC mixed-event collision pairs vs centrality", kTH1F, {centAxis});
    }

    if (doprocessMCGen) {
      histos.add("CutFlow/MCGen/hEventSelection", "MC generated event cut flow", kTH1F, {{6, -0.5f, 5.5f}});
      setHistogramBinLabels(histos.get<TH1>(HIST("CutFlow/MCGen/hEventSelection")),
                            {"All", "Vz", "INEL>0", "Centrality", "Final generated event", "Associated reco event accepted"});

      histos.add("InvMass/MCGen/Generated/hLambda1520", "THnSparse #Lambda(1520) generated", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("InvMass/MCGen/Generated/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) generated", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});

      histos.add("InvMass/MCGen/GeneratedInAcceptedEvents/hLambda1520", "THnSparse #Lambda(1520) generated, event-accepted", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("InvMass/MCGen/GeneratedInAcceptedEvents/hAntiLambda1520", "THnSparse #bar{#Lambda}(1520) generated, event-accepted", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});

      histos.add("MCGenQA/DaughterPt/hProtonFromLambda1520", "Gen proton from #Lambda(1520)", kTH1F, {ptAxis});
      histos.add("MCGenQA/DaughterPt/hProtonFromAntiLambda1520", "Gen proton from #bar{#Lambda}(1520)", kTH1F, {ptAxis});
      histos.add("MCGenQA/DaughterPt/hKaonFromLambda1520", "Gen kaon from #Lambda(1520)", kTH1F, {ptAxis});
      histos.add("MCGenQA/DaughterPt/hKaonFromAntiLambda1520", "Gen kaon from #bar{#Lambda}(1520)", kTH1F, {ptAxis});

      histos.add("CutFlow/MCGen/hGeneratedLambda1520Selection", "Generated #Lambda(1520) cut flow", kTH1F, {{8, -0.5f, 7.5f}});
      setHistogramBinLabels(histos.get<TH1>(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection")),
                            {"All particles", "producedByGenerator", "Lambda(1520) PDG", "Decay p + K", "Rapidity", "Daughter eta", "Eta", "Final"});

      histos.add("MCGenQA/Kinematics/hEtaVsPtBeforeCuts", "Generated #Lambda(1520) #eta-p_{T} (before kinematic cuts)", kTH2F, {etaAxis, ptAxis});
      histos.add("MCGenQA/Kinematics/hPhiVsRapidityBeforeCuts", "Generated #Lambda(1520) #varphi-y (before kinematic cuts)", kTH2F, {phiAxis, rapAxis});
      histos.add("MCGenQA/Kinematics/hEtaVsPtAfterCuts", "Generated #Lambda(1520) #eta-p_{T} (after kinematic cuts)", kTH2F, {etaAxis, ptAxis});
      histos.add("MCGenQA/Kinematics/hPhiVsRapidityAfterCuts", "Generated #Lambda(1520) #varphi-y (after kinematic cuts)", kTH2F, {phiAxis, rapAxis});
      histos.add("MCGenQA/Kinematics/hRapidityVsPtAccepted", "Generated #Lambda(1520) y vs p_{T} (final)", kTH2F, {ptAxis, rapAxis});
    }

    if (doprocessMCGenSpinAlignment) {
      histos.add("CosThetaStar/MCGen/Generated/hLambda1520", "cos#theta* generated #Lambda(1520) (spin-alignment process)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
      histos.add("CosThetaStar/MCGen/Generated/hAntiLambda1520", "cos#theta* generated #bar{#Lambda}(1520) (spin-alignment process)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
      histos.add("CosThetaStar/MCGen/GeneratedInAcceptedEvents/hLambda1520", "cos#theta* generated #Lambda(1520), event-accepted (spin-alignment process)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
      histos.add("CosThetaStar/MCGen/GeneratedInAcceptedEvents/hAntiLambda1520", "cos#theta* generated #bar{#Lambda}(1520), event-accepted (spin-alignment process)", kTHnSparseF, {massAxis, ptAxis, centAxis, cosThetaAxis});
    }

    if (doprocessEventFactor) {
      histos.add("EventLoss/hAcceptedRecoEventsVsTruthCentrality", "Accepted reconstructed events of generated INEL>0 events (Nrec, MC-truth centrality)", kTH1F, {centAxis});
      histos.add("EventLoss/hGeneratedEventsVsTruthCentrality", "Generated events passing truth |Zvtx| + truth INEL>0 cuts", kTH1F, {centAxis});
      histos.add("EventLoss/hGeneratedAcceptedEventsVsTruthCentrality", "Generated events with >=1 accepted reconstructed collision (truth centrality)", kTH1F, {centAxis});
      histos.add("EventLoss/hGeneratedAcceptedEventsVsRecoCentrality", "Generated events with >=1 accepted reconstructed collision (reconstructed centrality)", kTH1F, {centAxis});
      histos.add("EventLoss/hRecoCollisionMultiplicity", "Number of reconstructed collisions per generated collision", kTH1F, {{21, -0.5f, 20.5f}});
      histos.add("CutFlow/EventLoss/hGeneratedEventSelection", "Event Loss / Event Splitting: generated-event cut flow", kTH1F, {{6, -0.5f, 5.5f}});
      setHistogramBinLabels(histos.get<TH1>(HIST("CutFlow/EventLoss/hGeneratedEventSelection")),
                            {"All", "Vz", "INEL>0", "HasRecoColl", "RecoCollAccepted", "EventAccepted"});
    }

    if (doprocessSignalLoss) {
      histos.add("SignalLoss/Lambda1520/hGeneratedPtVsCentrality", "Generated #Lambda(1520), generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/AntiLambda1520/hGeneratedPtVsCentrality", "Generated #bar{#Lambda}(1520), generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Lambda1520/hAcceptedEventPtVsCentrality", "Generated #Lambda(1520), accepted events", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/AntiLambda1520/hAcceptedEventPtVsCentrality", "Generated #bar{#Lambda}(1520), accepted events", kTH2F, {ptAxis, centAxis});

      histos.add("SignalLoss/Kaon/hGeneratedTruePtVsCentrality", "Kaon, true pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Kaon/hGeneratedMtScaledPtVsCentrality", "Kaon, mT-scaled pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Kaon/hAcceptedEventTruePtVsCentrality", "Kaon, true pT, accepted events", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Kaon/hAcceptedEventMtScaledPtVsCentrality", "Kaon, mT-scaled pT, accepted events", kTH2F, {ptAxis, centAxis});

      histos.add("SignalLoss/Proton/hGeneratedTruePtVsCentrality", "Proton, true pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Proton/hGeneratedMtScaledPtVsCentrality", "Proton, mT-scaled pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Proton/hAcceptedEventTruePtVsCentrality", "Proton, true pT, accepted events", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Proton/hAcceptedEventMtScaledPtVsCentrality", "Proton, mT-scaled pT, accepted events", kTH2F, {ptAxis, centAxis});

      histos.add("SignalLoss/Lambda0/hGeneratedTruePtVsCentrality", "Lambda, true pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Lambda0/hGeneratedMtScaledPtVsCentrality", "Lambda, mT-scaled pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Lambda0/hAcceptedEventTruePtVsCentrality", "Lambda, true pT, accepted events", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Lambda0/hAcceptedEventMtScaledPtVsCentrality", "Lambda, mT-scaled pT, accepted events", kTH2F, {ptAxis, centAxis});

      histos.add("SignalLoss/Xi/hGeneratedTruePtVsCentrality", "Xi, true pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Xi/hGeneratedMtScaledPtVsCentrality", "Xi, mT-scaled pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Xi/hAcceptedEventTruePtVsCentrality", "Xi, true pT, accepted events", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Xi/hAcceptedEventMtScaledPtVsCentrality", "Xi, mT-scaled pT, accepted events", kTH2F, {ptAxis, centAxis});

      histos.add("SignalLoss/Omega/hGeneratedTruePtVsCentrality", "Omega, true pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Omega/hGeneratedMtScaledPtVsCentrality", "Omega, mT-scaled pT, generated", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Omega/hAcceptedEventTruePtVsCentrality", "Omega, true pT, accepted events", kTH2F, {ptAxis, centAxis});
      histos.add("SignalLoss/Omega/hAcceptedEventMtScaledPtVsCentrality", "Omega, mT-scaled pT, accepted events", kTH2F, {ptAxis, centAxis});
    }
  }

  template <typename CollisionType>
  float getCentrality(CollisionType const& collision)
  {
    switch (evSel.cfgCentralityEstimator) {
      case lambda1520analysispo_analysis::kFT0M:
        return collision.centFT0M();
      case lambda1520analysispo_analysis::kFT0A:
        return collision.centFT0A();
      case lambda1520analysispo_analysis::kFT0C:
        return collision.centFT0C();
      case lambda1520analysispo_analysis::kFV0A:
        return collision.centFV0A();
      default:
        return collision.centFT0M();
    }
  }

  template <typename CollisionType>
  bool isRecoInelGt0(CollisionType const& collision)
  {
    return collision.isInelGt0();
  }

  template <typename McParticleType>
  bool isGeneratedChargedPrimary(McParticleType const& mcPart)
  {
    if (!mcPart.isPhysicalPrimary()) {
      return false;
    }
    auto* pdgParticle = pdgDB->GetParticle(mcPart.pdgCode());
    if (!pdgParticle) {
      return false;
    }
    return std::abs(pdgParticle->Charge()) > MinAbsCharge;
  }

  template <typename McParticleTableType>
  bool isTruthInelGt0(McParticleTableType const& mcParticlesInCollision)
  {
    for (auto const& mcPart : mcParticlesInCollision) {
      if (std::abs(mcPart.eta()) < 1.f && isGeneratedChargedPrimary(mcPart)) {
        return true;
      }
    }
    return false;
  }

  template <typename McPartsT>
  int countGeneratedChargedPrimaries(McPartsT const& partsThisMc)
  {
    int nChargedPrimaries = 0;
    for (auto const& part : partsThisMc) {
      if (!part.isPhysicalPrimary()) {
        continue;
      }
      if (part.eta() > mcCfg.generatedMultiplicityEtaMax || part.eta() < mcCfg.generatedMultiplicityEtaMin) {
        continue;
      }
      auto* pdgParticle = pdgDB->GetParticle(part.pdgCode());
      if (!pdgParticle || std::abs(pdgParticle->Charge()) < MinAbsCharge) {
        continue;
      }
      ++nChargedPrimaries;
    }
    return nChargedPrimaries;
  }

  float getCentralityFromGeneratedMultiplicity(int nChargedPrimaries)
  {
    const auto& cuts = mcCfg.generatedMultiplicityLowerEdges.value;
    const auto& centres = mcCfg.centralityPerGeneratedMultiplicityClass.value;
    for (size_t i = 0; i < cuts.size(); ++i) {
      if (nChargedPrimaries >= cuts[i]) {
        return centres[i];
      }
    }
    return InvalidCentrality;
  }

  template <typename McCollisionType, typename McPartsT>
  float getTruthCentrality(McCollisionType const& mcCollision, McPartsT const& partsThisMc)
  {
    if (mcCfg.useGeneratedMultiplicityCentrality) {
      return getCentralityFromGeneratedMultiplicity(countGeneratedChargedPrimaries(partsThisMc));
    }
    return mcCollision.centFT0M();
  }

  template <typename Tag>
  void fillEventSelectionCutFlowBin(float bin)
  {
    if constexpr (std::is_same_v<Tag, DataCutFlowTag>) {
      histos.fill(HIST("CutFlow/Data/hEventSelection"), bin);
    } else if constexpr (std::is_same_v<Tag, McRecoCutFlowTag>) {
      histos.fill(HIST("CutFlow/MCReco/hEventSelection"), bin);
    }
  }

  template <typename Tag, typename CollisionType>
  bool passesEventSelectionImpl(CollisionType const& collision)
  {
    fillEventSelectionCutFlowBin<Tag>(0.f);

    if (std::abs(collision.posZ()) > evSel.cfgCutVertex) {
      return false;
    }
    fillEventSelectionCutFlowBin<Tag>(1.f);

    if (evSel.cfgUseNoITSROFrameBorderCut && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    fillEventSelectionCutFlowBin<Tag>(2.f);

    if (evSel.cfgUseSel8 && !collision.sel8()) {
      return false;
    }
    fillEventSelectionCutFlowBin<Tag>(3.f);

    if (evSel.cfgRequireRecoINELgt0 && !isRecoInelGt0(collision)) {
      return false;
    }
    fillEventSelectionCutFlowBin<Tag>(4.f);

    if (evSel.applyOccupancyInTimeRangeCut) {
      const int trackOccupancy = collision.trackOccupancyInTimeRange();
      if (trackOccupancy < evSel.cfgOccupancyMin || trackOccupancy > evSel.cfgOccupancyMax) {
        return false;
      }
    }
    fillEventSelectionCutFlowBin<Tag>(5.f);

    if (evSel.cfgUseNoSameBunchPileupCut && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    fillEventSelectionCutFlowBin<Tag>(6.f);

    if (evSel.cfgUseGoodZvtxFT0vsPVCut && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    fillEventSelectionCutFlowBin<Tag>(7.f);

    if (evSel.cfgUseGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    fillEventSelectionCutFlowBin<Tag>(8.f);

    fillEventSelectionCutFlowBin<Tag>(9.f);
    return true;
  }

  template <typename CollisionType>
  bool passesEventSelection(CollisionType const& collision)
  {
    return passesEventSelectionImpl<NoCutFlowTag>(collision);
  }

  template <typename CollisionGroup>
  AcceptedRecoCollision findFirstAcceptedRecoCollision(CollisionGroup const& collisions)
  {
    for (auto const& collision : collisions) {
      if (passesEventSelection(collision)) {
        return {true, getCentrality(collision)};
      }
    }
    return {};
  }

  template <typename TrackType>
  bool passesBasicTrackSelection(TrackType const& track)
  {
    if (trackCuts.applyTPCCrossedRowsCut &&
        track.tpcNClsCrossedRows() < trackCuts.cfgMinTPCCrossedRows) {
      return false;
    }
    if (trackCuts.applyTPCFoundClusterCut &&
        track.tpcNClsFound() < trackCuts.cfgMinTPCClusters) {
      return false;
    }
    if (trackCuts.applyITSClusterCut &&
        track.itsNCls() < trackCuts.cfgMinITSClusters) {
      return false;
    }
    if (trackCuts.applyTPCSharedClusterCut &&
        track.tpcNClsShared() > trackCuts.cfgMaxTPCSharedClusters) {
      return false;
    }
    if (trackCuts.applyTPCCrossedRowsOverFindableCut &&
        track.tpcNClsCrossedRows() < trackCuts.cfgMinCrossedRowsOverFindable * track.tpcNClsFindable()) {
      return false;
    }
    if (trackCuts.applyTPCChi2NClCut &&
        track.tpcChi2NCl() > trackCuts.cfgMaxTPCChi2NCl) {
      return false;
    }
    if (trackCuts.applyITSChi2NClCut &&
        track.itsChi2NCl() > trackCuts.cfgMaxITSChi2NCl) {
      return false;
    }
    if (trackCuts.requirePrimaryTrack && !track.isPrimaryTrack()) {
      return false;
    }
    if (trackCuts.requireGlobalTrackWoDCA && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (trackCuts.requireGlobalTrack && !track.isGlobalTrack()) {
      return false;
    }
    if (trackCuts.requirePVContributor && !track.isPVContributor()) {
      return false;
    }
    return true;
  }

  static int findBinIndex(std::vector<float> const& edges, float x)
  {
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
      if (x >= edges[i] && x < edges[i + 1]) {
        return static_cast<int>(i);
      }
    }
    return -1;
  }

  template <typename TrackType>
  bool passesDCAxyDCAzSelection(TrackType const& track, std::vector<float> const& ptEdges, std::vector<float> const& maxDCAxyPerBin)
  {
    if (std::abs(track.dcaZ()) >= dcaCuts.cfgCutDCAz) {
      return false;
    }
    const float pt = track.pt();
    const float dcaXY = std::abs(track.dcaXY());
    const int bin = findBinIndex(ptEdges, pt);
    if (bin < 0 || dcaXY >= maxDCAxyPerBin[bin]) {
      return false;
    }
    if (dcaCuts.isApplyPtDepDCAxyCut) {
      const float dcaConst = dcaCuts.dcaXYConst;
      const float dcaCoeff = dcaCuts.dcaXYPtCoeff;
      const float dcaExp = dcaCuts.dcaXYPtExponent;
      const float maxPtDepDCAxy = dcaConst + dcaCoeff / std::pow(pt, dcaExp);
      if (dcaXY > maxPtDepDCAxy) {
        return false;
      }
    }
    return true;
  }

  template <typename TrackType>
  bool passesProtonDCASelection(TrackType const& track)
  {
    return passesDCAxyDCAzSelection(track, mProtonDCAPtEdges, mProtonMaxDCAxy);
  }

  template <typename TrackType>
  bool passesKaonDCASelection(TrackType const& track)
  {
    return passesDCAxyDCAzSelection(track, mKaonDCAPtEdges, mKaonMaxDCAxy);
  }

  bool passesMomentumDependentPID(float p, bool hasTOF, float tpcNSigma, float tofNSigma, SpeciesPIDSettings const& cuts)
  {
    const float absTPC = std::abs(tpcNSigma);
    if (pidShared.useTPCOnlyPIDAtAllMomenta) {
      return absTPC < cuts.maxTPCNSigma;
    }
    if (!hasTOF) {
      if (p >= cuts.thresholdMomentum) {
        return false;
      }
      const int bin = findBinIndex(cuts.tpcOnlyBins, p);
      return bin >= 0 && absTPC < cuts.tpcOnlyCuts[bin];
    }
    const float absTOF = std::abs(tofNSigma);
    const int bin = findBinIndex(cuts.combinedBins, p);
    return bin >= 0 && (absTPC * absTPC + absTOF * absTOF) < cuts.combinedCuts[bin] * cuts.combinedCuts[bin];
  }

  template <typename TrackType>
  bool passesProtonPID(TrackType const& track, float totalMomentum)
  {
    if (protonPID.requireTOFForProton && !track.hasTOF()) {
      return false;
    }
    const bool hasTOF = track.hasTOF();
    const float tofNSigma = hasTOF ? track.tofNSigmaPr() : -999.f;
    return passesMomentumDependentPID(totalMomentum, hasTOF, track.tpcNSigmaPr(), tofNSigma, mProtonPIDCuts);
  }

  template <typename TrackType>
  bool isPionContaminatedKaonCandidate(TrackType const& track)
  {
    return track.pt() > pidShared.pionRejectionMinPt &&
           track.pt() < pidShared.pionRejectionMaxPt &&
           !track.hasTOF() &&
           std::abs(track.tpcNSigmaPi()) < pidShared.pionRejectionMaxTPCNSigma;
  }

  template <typename TrackType>
  bool passesKaonPID(TrackType const& track, float totalMomentum)
  {
    if (kaonPID.requireTOFForKaon && !track.hasTOF()) {
      return false;
    }
    if (pidShared.applyKaonPionRejection && isPionContaminatedKaonCandidate(track)) {
      return false;
    }
    const bool hasTOF = track.hasTOF();
    const float tofNSigma = hasTOF ? track.tofNSigmaKa() : -999.f;
    return passesMomentumDependentPID(totalMomentum, hasTOF, track.tpcNSigmaKa(), tofNSigma, mKaonPIDCuts);
  }

  template <typename TrackType>
  DaughterTrackCandidate makeTrackCandidate(TrackType const& track)
  {
    DaughterTrackCandidate cand{};
    cand.px = track.px();
    cand.py = track.py();
    cand.pz = track.pz();
    cand.pt = track.pt();
    cand.eta = track.eta();
    cand.phi = track.phi();
    cand.mom = RecoDecay::p(track.px(), track.py(), track.pz());
    cand.dcaXY = track.dcaXY();
    cand.dcaZ = track.dcaZ();
    cand.sign = track.sign();
    cand.hasTOF = track.hasTOF();
    cand.tpcNSigmaPr = track.tpcNSigmaPr();
    cand.tpcNSigmaKa = track.tpcNSigmaKa();
    if (cand.hasTOF) {
      cand.tofNSigmaPr = track.tofNSigmaPr();
      cand.tofNSigmaKa = track.tofNSigmaKa();
    }
    cand.globalIndex = track.globalIndex();
    return cand;
  }

  template <typename TrackType>
  void attachMCTruthToCandidate(DaughterTrackCandidate& cand, TrackType const& track)
  {
    if (!track.has_mcParticle()) {
      return;
    }
    auto const mcPart = track.mcParticle();
    cand.hasMcLabel = true;
    cand.mcPdgCode = mcPart.pdgCode();
    if (!mcPart.has_mothers()) {
      return;
    }
    for (auto const& mother : mcPart.template mothers_as<aod::McParticles>()) {
      if (std::abs(mother.pdgCode()) == mLambda1520AbsPdgCode) {
        cand.mcLambda1520MotherIndex = mother.globalIndex();
        cand.mcLambda1520MotherPdgCode = mother.pdgCode();
        break;
      }
    }
  }

  template <bool IsMC, typename TrackCollection>
  void buildProtonKaonCandidatePools(TrackCollection const& tracks, std::vector<DaughterTrackCandidate>& protons, std::vector<DaughterTrackCandidate>& kaons)
  {
    protons.clear();
    kaons.clear();
    for (auto const& track : tracks) {
      if (!passesBasicTrackSelection(track)) {
        continue;
      }
      const float mom = RecoDecay::p(track.px(), track.py(), track.pz());
      const bool isProton = passesProtonPID(track, mom) && passesProtonDCASelection(track);
      const bool isKaon = passesKaonPID(track, mom) && passesKaonDCASelection(track);
      if (!isProton && !isKaon) {
        continue;
      }
      DaughterTrackCandidate cand = makeTrackCandidate(track);
      if constexpr (IsMC) {
        attachMCTruthToCandidate(cand, track);
      }
      if (isProton) {
        protons.push_back(cand);
      }
      if (isKaon) {
        kaons.push_back(cand);
      }
    }
  }

  void computeDaughterPairGeometry(DaughterTrackCandidate const& pr, DaughterTrackCandidate const& ka, float& deltaEta, float& deltaPhi, float& openingAngle)
  {
    deltaEta = std::abs(pr.eta - ka.eta);
    deltaPhi = std::abs(pr.phi - ka.phi);
    if (deltaPhi > o2::constants::math::PI) {
      deltaPhi = o2::constants::math::TwoPI - deltaPhi;
    }
    const float dot = pr.px * ka.px + pr.py * ka.py + pr.pz * ka.pz;
    openingAngle = std::acos(std::clamp(dot / (pr.mom * ka.mom), -1.f, 1.f));
  }

  bool passesDaughterAngularCuts(DaughterTrackCandidate const& pr, DaughterTrackCandidate const& ka)
  {
    if (!pairCuts.applyDaughterDeltaEtaPhiCuts && !pairCuts.applyDeepAngleCut && !pairCuts.applyOpeningAngleVeto) {
      return true;
    }
    float deltaEta = 0.f, deltaPhi = 0.f, openingAngle = 0.f;
    computeDaughterPairGeometry(pr, ka, deltaEta, deltaPhi, openingAngle);
    if (pairCuts.applyDaughterDeltaEtaPhiCuts) {
      const float maxDaughterDeltaEtaValue = pairCuts.maxDaughterDeltaEta;
      const float maxDaughterDeltaPhiValue = pairCuts.maxDaughterDeltaPhi;
      if (deltaEta >= maxDaughterDeltaEtaValue || deltaPhi >= maxDaughterDeltaPhiValue) {
        return false;
      }
    }
    if (pairCuts.applyDeepAngleCut) {
      const float minOpeningAngle = pairCuts.deepAngleCutValue;
      if (openingAngle < minOpeningAngle) {
        return false;
      }
    }
    if (pairCuts.applyOpeningAngleVeto) {
      const float openingAngleVetoLow = pairCuts.openingAngleVetoMin;
      const float openingAngleVetoHigh = pairCuts.openingAngleVetoMax;
      if (openingAngle > openingAngleVetoLow && openingAngle < openingAngleVetoHigh) {
        return false;
      }
    }
    return true;
  }

  bool buildProtonKaonPair(DaughterTrackCandidate const& pr, DaughterTrackCandidate const& ka, float& pairMass, float& pairPt, float& pairY)
  {
    if (!passesDaughterAngularCuts(pr, ka)) {
      return false;
    }
    const std::array<std::array<float, 3>, 2> momenta = {
      std::array<float, 3>{pr.px, pr.py, pr.pz},
      std::array<float, 3>{ka.px, ka.py, ka.pz}};
    pairMass = RecoDecay::m(momenta, kProtonKaonMasses);
    pairPt = RecoDecay::pt(std::array{pr.px + ka.px, pr.py + ka.py});
    pairY = RecoDecay::y(std::array{pr.px + ka.px, pr.py + ka.py, pr.pz + ka.pz}, pairMass);

    const float minRapidity = trackCuts.cfgMinY;
    const float maxRapidity = trackCuts.cfgMaxY;
    if (pairY < minRapidity || pairY > maxRapidity) {
      return false;
    }
    if (pairCuts.applyPairPtMassUpperLimits) {
      const float maxPairPt = pairCuts.cfgMaxPairPt;
      const float maxPairMass = pairCuts.cfgMaxPairMass;
      if (pairPt >= maxPairPt || pairMass >= maxPairMass) {
        return false;
      }
    }
    return true;
  }

  void fillPairQABeforeCuts(DaughterTrackCandidate const& pr, DaughterTrackCandidate const& ka)
  {
    float deltaEta = 0.f, deltaPhi = 0.f, openingAngle = 0.f;
    computeDaughterPairGeometry(pr, ka, deltaEta, deltaPhi, openingAngle);
    histos.fill(HIST("PairQA/Data/hOpeningAngleBeforeCuts"), openingAngle);
    histos.fill(HIST("PairQA/Data/hDeltaEtaBeforeCuts"), deltaEta);
    histos.fill(HIST("PairQA/Data/hDeltaPhiBeforeCuts"), deltaPhi);
  }

  void fillPairQAAfterCuts(DaughterTrackCandidate const& pr, DaughterTrackCandidate const& ka)
  {
    float deltaEta = 0.f, deltaPhi = 0.f, openingAngle = 0.f;
    computeDaughterPairGeometry(pr, ka, deltaEta, deltaPhi, openingAngle);
    histos.fill(HIST("PairQA/Data/hOpeningAngleAfterCuts"), openingAngle);
    histos.fill(HIST("PairQA/Data/hDeltaEtaAfterCuts"), deltaEta);
    histos.fill(HIST("PairQA/Data/hDeltaPhiAfterCuts"), deltaPhi);
  }

  static float vecDot(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  static std::array<float, 3> vecSub(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }

  static std::array<float, 3> vecCross(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
  }

  static constexpr float VecUnitEpsilon = 1e-9f;
  static std::array<float, 3> vecUnit(const std::array<float, 3>& a)
  {
    float mag = std::sqrt(vecDot(a, a));
    if (mag < VecUnitEpsilon) {
      return {0.f, 0.f, 0.f};
    }
    return {a[0] / mag, a[1] / mag, a[2] / mag};
  }

  static constexpr float BoostMomentumEpsilon = 1e-6f;
  static std::array<float, 3> boostToRestFrame(float motherE, const std::array<float, 3>& motherP,
                                               float dE, const std::array<float, 3>& dP)
  {
    float pMag2 = vecDot(motherP, motherP);
    float pMag = std::sqrt(pMag2);
    if (pMag < BoostMomentumEpsilon) {
      return dP;
    }

    std::array<float, 3> nHat = {motherP[0] / pMag, motherP[1] / pMag, motherP[2] / pMag};
    float beta = pMag / motherE;
    float motherMass2 = motherE * motherE - pMag2;
    float gamma = motherMass2 > 0.f ? motherE / std::sqrt(motherMass2) : 1.f;

    float pParallel = vecDot(dP, nHat);
    std::array<float, 3> pPerp = {dP[0] - pParallel * nHat[0], dP[1] - pParallel * nHat[1], dP[2] - pParallel * nHat[2]};
    float pParallelBoosted = gamma * (pParallel - beta * dE);

    return {pPerp[0] + pParallelBoosted * nHat[0], pPerp[1] + pParallelBoosted * nHat[1], pPerp[2] + pParallelBoosted * nHat[2]};
  }

  static DecayKinematics makeDecayKinematics(std::array<float, 3> const& protonMom, std::array<float, 3> const& kaonMom)
  {
    DecayKinematics decay;
    decay.protonMom = protonMom;
    decay.kaonMom = kaonMom;
    decay.motherMom = {protonMom[0] + kaonMom[0], protonMom[1] + kaonMom[1], protonMom[2] + kaonMom[2]};
    decay.protonE = std::sqrt(vecDot(protonMom, protonMom) + ProtonMass * ProtonMass);
    decay.kaonE = std::sqrt(vecDot(kaonMom, kaonMom) + KaonMass * KaonMass);
    decay.motherE = decay.protonE + decay.kaonE;
    return decay;
  }

  static DecayKinematics makeDecayKinematics(DaughterTrackCandidate const& pr, DaughterTrackCandidate const& ka)
  {
    return makeDecayKinematics(std::array<float, 3>{pr.px, pr.py, pr.pz}, std::array<float, 3>{ka.px, ka.py, ka.pz});
  }

  template <typename McParticleType>
  static DecayKinematics makeGeneratedDecayKinematics(McParticleType const& mother, McParticleType const& proton, McParticleType const& kaon)
  {
    DecayKinematics decay;
    decay.motherMom = {mother.px(), mother.py(), mother.pz()};
    decay.motherE = mother.e();
    decay.protonMom = {proton.px(), proton.py(), proton.pz()};
    decay.protonE = proton.e();
    decay.kaonMom = {kaon.px(), kaon.py(), kaon.pz()};
    decay.kaonE = kaon.e();
    return decay;
  }

  template <typename McParticleType>
  static float computeGeneratedMass(McParticleType const& particle)
  {
    const float energy = particle.e();
    const float px = particle.px();
    const float py = particle.py();
    const float pz = particle.pz();
    return std::sqrt(std::max(0.f, energy * energy - px * px - py * py - pz * pz));
  }

  float computeCosThetaStar(DecayKinematics const& decay)
  {
    const std::array<float, 3> protonCM = boostToRestFrame(decay.motherE, decay.motherMom, decay.protonE, decay.protonMom);
    const std::array<float, 3> kaonCM = boostToRestFrame(decay.motherE, decay.motherMom, decay.kaonE, decay.kaonMom);
    const std::array<float, 3> referenceCM = vecUnit(spinAlignment.referenceFrameDaughter == 0 ? protonCM : kaonCM);

    if (spinAlignment.activateHelicityFrame) {
      return vecDot(referenceCM, vecUnit(decay.motherMom));
    }

    if (spinAlignment.activateCollinsSoperFrame || spinAlignment.activateProductionFrame || spinAlignment.activateBeamAxisFrame) {
      const std::array<float, 3> beam1CM = vecUnit(boostToRestFrame(decay.motherE, decay.motherMom, mBeam1Energy, mBeam1Momentum));
      const std::array<float, 3> beam2CM = vecUnit(boostToRestFrame(decay.motherE, decay.motherMom, mBeam2Energy, mBeam2Momentum));
      if (spinAlignment.activateCollinsSoperFrame) {
        return vecDot(referenceCM, vecUnit(vecSub(beam1CM, beam2CM)));
      }
      if (spinAlignment.activateProductionFrame) {
        return vecDot(referenceCM, vecUnit(vecCross(beam1CM, vecUnit(decay.motherMom))));
      }
      return vecDot(referenceCM, beam1CM);
    }

    if (spinAlignment.activateRandomFrame) {
      const double randomCos = mUniform01(mRandomEngine) * 2. - 1.;
      const double randomPhi = mUniform01(mRandomEngine) * o2::constants::math::TwoPI;
      const double randomSin = std::sqrt(1. - randomCos * randomCos);
      const std::array<float, 3> randomAxis = {static_cast<float>(randomSin * std::cos(randomPhi)),
                                               static_cast<float>(randomSin * std::sin(randomPhi)),
                                               static_cast<float>(randomCos)};
      return vecDot(referenceCM, randomAxis);
    }

    return 0.f;
  }

  template <CosThetaStarSample Sample>
  void fillCosThetaStarHistograms(float candidateMass, float candidatePt, DecayKinematics const& decay, float centrality, bool isParticle)
  {
    if (!isSpinAlignmentActive()) {
      return;
    }

    const float cosTheta = computeCosThetaStar(decay);

    if constexpr (Sample == CosThetaStarSample::Signal) {
      if (isParticle) {
        histos.fill(HIST("CosThetaStar/Data/SameEvent/hLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      } else {
        histos.fill(HIST("CosThetaStar/Data/SameEvent/hAntiLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      }
    } else if constexpr (Sample == CosThetaStarSample::LikeSign) {
      if (isParticle) {
        histos.fill(HIST("CosThetaStar/Data/SameEvent/hLikeSignProtonPlusKaonPlus"), candidateMass, candidatePt, centrality, cosTheta);
      } else {
        histos.fill(HIST("CosThetaStar/Data/SameEvent/hLikeSignProtonMinusKaonMinus"), candidateMass, candidatePt, centrality, cosTheta);
      }
    } else if constexpr (Sample == CosThetaStarSample::Rotated) {
      if (isParticle) {
        histos.fill(HIST("CosThetaStar/Data/RotatedBackground/hLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      } else {
        histos.fill(HIST("CosThetaStar/Data/RotatedBackground/hAntiLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      }
    } else if constexpr (Sample == CosThetaStarSample::Mixed) {
      if (isParticle) {
        histos.fill(HIST("CosThetaStar/Data/MixedEvent/hLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      } else {
        histos.fill(HIST("CosThetaStar/Data/MixedEvent/hAntiLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      }
    } else if constexpr (Sample == CosThetaStarSample::LikeSignMixed) {
      if (isParticle) {
        histos.fill(HIST("CosThetaStar/Data/MixedEvent/hLikeSignProtonPlusKaonPlus"), candidateMass, candidatePt, centrality, cosTheta);
      } else {
        histos.fill(HIST("CosThetaStar/Data/MixedEvent/hLikeSignProtonMinusKaonMinus"), candidateMass, candidatePt, centrality, cosTheta);
      }
    } else if constexpr (Sample == CosThetaStarSample::TruthMatched) {
      if (isParticle) {
        histos.fill(HIST("CosThetaStar/MCReco/TruthMatched/hLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      } else {
        histos.fill(HIST("CosThetaStar/MCReco/TruthMatched/hAntiLambda1520"), candidateMass, candidatePt, centrality, cosTheta);
      }
    }
  }

  void fillInvariantMassSameEvent(int protonSign, int kaonSign, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (isSpinAlignmentActive()) {
      return;
    }
    if (protonSign * kaonSign < 0) {
      if (protonSign > 0) {
        histos.fill(HIST("InvMass/Data/SameEvent/hLambda1520"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("InvMass/Data/SameEvent/hAntiLambda1520"), pairMass, pairPt, centrality, rapidity);
      }
    } else {
      if (protonSign > 0) {
        histos.fill(HIST("InvMass/Data/SameEvent/hLikeSignProtonPlusKaonPlus"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("InvMass/Data/SameEvent/hLikeSignProtonMinusKaonMinus"), pairMass, pairPt, centrality, rapidity);
      }
    }
  }

  void fillInvariantMassMixedEvent(int protonSign, int kaonSign, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (isSpinAlignmentActive()) {
      return;
    }
    if (protonSign * kaonSign < 0) {
      if (protonSign > 0) {
        histos.fill(HIST("InvMass/Data/MixedEvent/hLambda1520"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("InvMass/Data/MixedEvent/hAntiLambda1520"), pairMass, pairPt, centrality, rapidity);
      }
    } else if (mixingCfg.cfgFillMixLikeSign) {
      if (protonSign > 0) {
        histos.fill(HIST("InvMass/Data/MixedEvent/hLikeSignProtonPlusKaonPlus"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("InvMass/Data/MixedEvent/hLikeSignProtonMinusKaonMinus"), pairMass, pairPt, centrality, rapidity);
      }
    }
  }

  void fillInvariantMassMCRecoAllPairs(int protonSign, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (protonSign > 0) {
      histos.fill(HIST("InvMass/MCReco/AllPairs/hLambda1520"), pairMass, pairPt, centrality, rapidity);
    } else {
      histos.fill(HIST("InvMass/MCReco/AllPairs/hAntiLambda1520"), pairMass, pairPt, centrality, rapidity);
    }
  }

  void fillInvariantMassMCRecoTruthMatched(int motherPdg, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (motherPdg > 0) {
      histos.fill(HIST("InvMass/MCReco/TruthMatched/hLambda1520"), pairMass, pairPt, centrality, rapidity);
    } else {
      histos.fill(HIST("InvMass/MCReco/TruthMatched/hAntiLambda1520"), pairMass, pairPt, centrality, rapidity);
    }
  }

  void fillInvariantMassMCRecoMixedEvent(int protonSign, int kaonSign, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (protonSign * kaonSign >= 0) {
      return;
    }
    if (protonSign > 0) {
      histos.fill(HIST("InvMass/MCReco/MixedEvent/hLambda1520"), pairMass, pairPt, centrality, rapidity);
    } else {
      histos.fill(HIST("InvMass/MCReco/MixedEvent/hAntiLambda1520"), pairMass, pairPt, centrality, rapidity);
    }
  }

  void fillRotatedPairBackground(DaughterTrackCandidate const& pr, DaughterTrackCandidate const& ka, float centrality)
  {
    const int nRot = rotBkg.numberOfRotations;
    if (nRot <= 0) {
      return;
    }
    const float weight = 1.f / static_cast<float>(nRot);
    const float rotationHalfWindow = o2::constants::math::PI / static_cast<float>(rotBkg.rotationAngleWindow);
    const bool rotateProtonInsteadOfKaon = rotBkg.rotateProtonInsteadOfKaon;
    const DaughterTrackCandidate& rotated = rotateProtonInsteadOfKaon ? pr : ka;

    for (int iRot = 0; iRot < nRot; ++iRot) {
      float rotationAngle = o2::constants::math::PI;
      if (nRot > 1) {
        rotationAngle = (o2::constants::math::PI - rotationHalfWindow) +
                        static_cast<float>(iRot) * (2.f * rotationHalfWindow / static_cast<float>(nRot - 1));
      }
      const float rotatedPhi = RecoDecay::constrainAngle(rotated.phi + rotationAngle, 0.f);
      const float rotatedPx = rotated.pt * std::cos(rotatedPhi);
      const float rotatedPy = rotated.pt * std::sin(rotatedPhi);

      const std::array<float, 3> momPr = rotateProtonInsteadOfKaon ? std::array<float, 3>{rotatedPx, rotatedPy, pr.pz}
                                                                   : std::array<float, 3>{pr.px, pr.py, pr.pz};
      const std::array<float, 3> momKa = rotateProtonInsteadOfKaon ? std::array<float, 3>{ka.px, ka.py, ka.pz}
                                                                   : std::array<float, 3>{rotatedPx, rotatedPy, ka.pz};
      const std::array<std::array<float, 3>, 2> rotMomenta = {momPr, momKa};

      const float rotatedMass = RecoDecay::m(rotMomenta, kProtonKaonMasses);
      const float pairPxSum = momPr[0] + momKa[0];
      const float pairPySum = momPr[1] + momKa[1];
      const float pairPzSum = momPr[2] + momKa[2];
      const float rotatedPt = RecoDecay::pt(std::array{pairPxSum, pairPySum});
      const float rotatedRapidity = RecoDecay::y(std::array{pairPxSum, pairPySum, pairPzSum}, rotatedMass);

      const float minRapidity = trackCuts.cfgMinY;
      const float maxRapidity = trackCuts.cfgMaxY;
      if (rotatedRapidity < minRapidity || rotatedRapidity > maxRapidity) {
        continue;
      }
      if (pairCuts.applyPairPtMassUpperLimits) {
        const float maxPairPt = pairCuts.cfgMaxPairPt;
        const float maxPairMass = pairCuts.cfgMaxPairMass;
        if (rotatedPt >= maxPairPt || rotatedMass >= maxPairMass) {
          continue;
        }
      }

      if (!isSpinAlignmentActive()) {
        if (pr.sign > 0) {
          histos.fill(HIST("InvMass/Data/RotatedBackground/hLambda1520"), rotatedMass, rotatedPt, centrality, rotatedRapidity, weight);
        } else {
          histos.fill(HIST("InvMass/Data/RotatedBackground/hAntiLambda1520"), rotatedMass, rotatedPt, centrality, rotatedRapidity, weight);
        }
      }

      if (isSpinAlignmentActive() && spinAlignment.fillRotatedBackgroundCosThetaStar) {
        fillCosThetaStarHistograms<CosThetaStarSample::Rotated>(
          rotatedMass, rotatedPt, makeDecayKinematics(momPr, momKa), centrality, pr.sign > 0);
      }
    }
  }

  template <bool IsMixed, bool IsMC = false>
  void fillInvariantMassFromCandidatePools(
    std::vector<DaughterTrackCandidate> const& protonPool,
    std::vector<DaughterTrackCandidate> const& kaonPool,
    float centrality,
    bool fillPairQA = false)
  {
    for (auto const& pr : protonPool) {
      for (auto const& ka : kaonPool) {

        if constexpr (!IsMixed) {
          if (pr.globalIndex == ka.globalIndex) {
            continue;
          }
          if (fillPairQA) {
            fillPairQABeforeCuts(pr, ka);
          }
        }

        float pairMass = 0.f, pairPt = 0.f, pairY = 0.f;
        if (!buildProtonKaonPair(pr, ka, pairMass, pairPt, pairY)) {
          continue;
        }

        if constexpr (!IsMixed) {
          if (fillPairQA) {
            fillPairQAAfterCuts(pr, ka);
          }
        }

        if constexpr (IsMixed) {
          if constexpr (IsMC) {
            fillInvariantMassMCRecoMixedEvent(pr.sign, ka.sign, pairPt, pairMass, centrality, pairY);
          } else {
            fillInvariantMassMixedEvent(pr.sign, ka.sign, pairPt, pairMass, centrality, pairY);

            if (isSpinAlignmentActive() && spinAlignment.fillMixedEventCosThetaStar) {
              const DecayKinematics decay = makeDecayKinematics(pr, ka);
              if (pr.sign * ka.sign < 0) {
                fillCosThetaStarHistograms<CosThetaStarSample::Mixed>(pairMass, pairPt, decay, centrality, pr.sign > 0);
              } else if (mixingCfg.cfgFillMixLikeSign) {
                fillCosThetaStarHistograms<CosThetaStarSample::LikeSignMixed>(pairMass, pairPt, decay, centrality, pr.sign > 0);
              }
            }
          }
        } else {
          fillInvariantMassSameEvent(pr.sign, ka.sign, pairPt, pairMass, centrality, pairY);

          if (isSpinAlignmentActive()) {
            const DecayKinematics decay = makeDecayKinematics(pr, ka);
            if (pr.sign * ka.sign < 0) {
              fillCosThetaStarHistograms<CosThetaStarSample::Signal>(pairMass, pairPt, decay, centrality, pr.sign > 0);
            } else {
              fillCosThetaStarHistograms<CosThetaStarSample::LikeSign>(pairMass, pairPt, decay, centrality, pr.sign > 0);
            }
          }

          if (rotBkg.enableRotationalBackground && pr.sign * ka.sign < 0) {
            fillRotatedPairBackground(pr, ka, centrality);
          }
        }
      }
    }
  }

  template <typename TrackType>
  void fillSelectedProtonQA(TrackType const& track, float totalMomentum, float centrality)
  {
    const float pt = track.pt();
    const float mom = totalMomentum;
    const float tpcNSigPr = track.tpcNSigmaPr();
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hDCAxyVsPt"), pt, track.dcaXY());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hDCAzVsPt"), pt, track.dcaZ());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTPCvsP"), mom, tpcNSigPr);
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTPCvsPt"), pt, tpcNSigPr);
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTPCvsCentrality"), centrality, tpcNSigPr);
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTPCPionVsPt"), pt, track.tpcNSigmaPi());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTPCKaonVsPt"), pt, track.tpcNSigmaKa());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hTPCSignalVsInnerP"), track.tpcInnerParam(), track.tpcSignal());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hTPCClustersVsPhi"), track.tpcNClsFound(), track.phi());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hTPCCrossedRowsVsPt"), pt, track.tpcNClsCrossedRows());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hTPCFoundClustersVsPt"), pt, track.tpcNClsFound());
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hPt"), pt);
    histos.fill(HIST("TrackQA/Data/Selected/Proton/hEta"), track.eta());
    if (track.hasTOF()) {
      const float tofNSigPr = track.tofNSigmaPr();
      histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTOFvsP"), mom, tofNSigPr);
      histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTOFvsPt"), pt, tofNSigPr);
      histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTOFvsCentrality"), centrality, tofNSigPr);
      histos.fill(HIST("TrackQA/Data/Selected/Proton/hNSigmaTOFvsNSigmaTPC"), tpcNSigPr, tofNSigPr);
    }
  }

  template <typename TrackType>
  void fillSelectedKaonQA(TrackType const& track, float totalMomentum, float centrality)
  {
    const float pt = track.pt();
    const float mom = totalMomentum;
    const float tpcNSigKa = track.tpcNSigmaKa();
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hDCAxyVsPt"), pt, track.dcaXY());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hDCAzVsPt"), pt, track.dcaZ());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTPCvsP"), mom, tpcNSigKa);
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTPCvsPt"), pt, tpcNSigKa);
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTPCvsCentrality"), centrality, tpcNSigKa);
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTPCPionVsPt"), pt, track.tpcNSigmaPi());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTPCProtonVsPt"), pt, track.tpcNSigmaPr());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hTPCSignalVsInnerP"), track.tpcInnerParam(), track.tpcSignal());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hTPCClustersVsPhi"), track.tpcNClsFound(), track.phi());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hTPCCrossedRowsVsPt"), pt, track.tpcNClsCrossedRows());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hTPCFoundClustersVsPt"), pt, track.tpcNClsFound());
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hPt"), pt);
    histos.fill(HIST("TrackQA/Data/Selected/Kaon/hEta"), track.eta());
    if (track.hasTOF()) {
      const float tofNSigKa = track.tofNSigmaKa();
      histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsP"), mom, tofNSigKa);
      histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsPt"), pt, tofNSigKa);
      histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsCentrality"), centrality, tofNSigKa);
      histos.fill(HIST("TrackQA/Data/Selected/Kaon/hNSigmaTOFvsNSigmaTPC"), tpcNSigKa, tofNSigKa);
    }
  }

  template <typename TrackType>
  void fillSelectedProtonQAMC(TrackType const& track, float totalMomentum, float centrality)
  {
    const float pt = track.pt();
    const float mom = totalMomentum;
    const float tpcNSigPr = track.tpcNSigmaPr();
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hDCAxyVsPt"), pt, track.dcaXY());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hDCAzVsPt"), pt, track.dcaZ());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTPCvsP"), mom, tpcNSigPr);
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTPCvsPt"), pt, tpcNSigPr);
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTPCvsCentrality"), centrality, tpcNSigPr);
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTPCPionVsPt"), pt, track.tpcNSigmaPi());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTPCKaonVsPt"), pt, track.tpcNSigmaKa());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hTPCSignalVsInnerP"), track.tpcInnerParam(), track.tpcSignal());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hTPCClustersVsPhi"), track.tpcNClsFound(), track.phi());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hTPCCrossedRowsVsPt"), pt, track.tpcNClsCrossedRows());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hTPCFoundClustersVsPt"), pt, track.tpcNClsFound());
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hPt"), pt);
    histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hEta"), track.eta());
    if (track.hasTOF()) {
      const float tofNSigPr = track.tofNSigmaPr();
      histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsP"), mom, tofNSigPr);
      histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsPt"), pt, tofNSigPr);
      histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsCentrality"), centrality, tofNSigPr);
      histos.fill(HIST("TrackQA/MCReco/Selected/Proton/hNSigmaTOFvsNSigmaTPC"), tpcNSigPr, tofNSigPr);
    }
  }

  template <typename TrackType>
  void fillSelectedKaonQAMC(TrackType const& track, float totalMomentum, float centrality)
  {
    const float pt = track.pt();
    const float mom = totalMomentum;
    const float tpcNSigKa = track.tpcNSigmaKa();
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hDCAxyVsPt"), pt, track.dcaXY());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hDCAzVsPt"), pt, track.dcaZ());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCvsP"), mom, tpcNSigKa);
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCvsPt"), pt, tpcNSigKa);
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCvsCentrality"), centrality, tpcNSigKa);
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCPionVsPt"), pt, track.tpcNSigmaPi());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTPCProtonVsPt"), pt, track.tpcNSigmaPr());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hTPCSignalVsInnerP"), track.tpcInnerParam(), track.tpcSignal());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hTPCClustersVsPhi"), track.tpcNClsFound(), track.phi());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hTPCCrossedRowsVsPt"), pt, track.tpcNClsCrossedRows());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hTPCFoundClustersVsPt"), pt, track.tpcNClsFound());
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hPt"), pt);
    histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hEta"), track.eta());
    if (track.hasTOF()) {
      const float tofNSigKa = track.tofNSigmaKa();
      histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsP"), mom, tofNSigKa);
      histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsPt"), pt, tofNSigKa);
      histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsCentrality"), centrality, tofNSigKa);
      histos.fill(HIST("TrackQA/MCReco/Selected/Kaon/hNSigmaTOFvsNSigmaTPC"), tpcNSigKa, tofNSigKa);
    }
  }

  template <typename McPart>
  bool findProtonKaonDaughters(McPart const& mother, McPart& proton, McPart& kaon)
  {
    if (!mother.has_daughters()) {
      return false;
    }
    bool hasProtonDaughter = false, hasKaonDaughter = false;
    for (auto const& d : mother.template daughters_as<aod::McParticles>()) {
      const int absPdgCode = std::abs(d.pdgCode());
      if (absPdgCode == lambda1520analysispo_analysis::PdgProton) {
        hasProtonDaughter = true;
        proton = d;
      } else if (absPdgCode == lambda1520analysispo_analysis::PdgKaon) {
        hasKaonDaughter = true;
        kaon = d;
      }
    }
    return hasProtonDaughter && hasKaonDaughter;
  }

  struct GeneratedDaughterMatch {
    bool found = false;
    bool trackQuality = false;
    bool pid = false;
    bool dca = false;
    DaughterTrackCandidate cand{};
  };

  template <typename TracksSlice, typename McPart>
  GeneratedDaughterMatch findRecoTrackOfGeneratedDaughter(TracksSlice const& tracksSlice, McPart const& mcDaughter, bool isProton)
  {
    GeneratedDaughterMatch result;
    for (auto const& trk : tracksSlice) {
      if (!trk.has_mcParticle() || trk.mcParticle().globalIndex() != mcDaughter.globalIndex()) {
        continue;
      }
      result.found = true;
      result.trackQuality = passesBasicTrackSelection(trk);
      const float mom = RecoDecay::p(trk.px(), trk.py(), trk.pz());
      if (isProton) {
        result.pid = passesProtonPID(trk, mom);
        result.dca = passesProtonDCASelection(trk);
      } else {
        result.pid = passesKaonPID(trk, mom);
        result.dca = passesKaonDCASelection(trk);
      }
      result.cand = makeTrackCandidate(trk);
      break;
    }
    return result;
  }

  template <typename TracksSlice>
  void fillGeneratedToReconstructedCutFlow(TracksSlice const& tracksOfCollision, aod::McParticles const& mcParticles, int64_t mcCollisionId)
  {
    auto genParts = mcParticles.sliceBy(perMcCollision, mcCollisionId);
    for (auto const& gen : genParts) {
      if (!gen.producedByGenerator() || std::abs(gen.pdgCode()) != mLambda1520AbsPdgCode) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 0.f);

      auto genProton = gen;
      auto genKaon = gen;
      if (!findProtonKaonDaughters(gen, genProton, genKaon)) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 1.f);

      const float minRapidity = trackCuts.cfgMinY;
      const float maxRapidity = trackCuts.cfgMaxY;
      if (gen.y() < minRapidity || gen.y() > maxRapidity) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 2.f);

      const float maxTrackEta = trackCuts.cfgCutEta;
      const float minTrackPt = trackCuts.cfgCutPt;
      const bool protonInAcceptance = std::abs(genProton.eta()) < maxTrackEta && genProton.pt() > minTrackPt;
      const bool kaonInAcceptance = std::abs(genKaon.eta()) < maxTrackEta && genKaon.pt() > minTrackPt;
      if (protonInAcceptance) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 3.f);
      }
      if (kaonInAcceptance) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 4.f);
      }
      if (!protonInAcceptance || !kaonInAcceptance) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 5.f);

      const auto protonMatch = findRecoTrackOfGeneratedDaughter(tracksOfCollision, genProton, true);
      const auto kaonMatch = findRecoTrackOfGeneratedDaughter(tracksOfCollision, genKaon, false);
      if (protonMatch.found) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 6.f);
      }
      if (kaonMatch.found) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 7.f);
      }
      if (!protonMatch.found || !kaonMatch.found) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 8.f);

      if (protonMatch.trackQuality) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 9.f);
      }
      if (kaonMatch.trackQuality) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 10.f);
      }
      if (!protonMatch.trackQuality || !kaonMatch.trackQuality) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 11.f);

      if (protonMatch.pid) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 12.f);
      }
      if (kaonMatch.pid) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 13.f);
      }
      if (!protonMatch.pid || !kaonMatch.pid) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 14.f);

      if (protonMatch.dca) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 15.f);
      }
      if (kaonMatch.dca) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 16.f);
      }
      if (!protonMatch.dca || !kaonMatch.dca) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 17.f);

      float pairMass = 0.f, pairPt = 0.f, pairY = 0.f;
      if (buildProtonKaonPair(protonMatch.cand, kaonMatch.cand, pairMass, pairPt, pairY)) {
        histos.fill(HIST("CutFlow/MCReco/hGeneratedToReconstructed"), 18.f);
      }
    }
  }

  static float computeMtScaledPt(float pt, float referenceMass)
  {
    const float mtScaledPtSquared = pt * pt + referenceMass * referenceMass - Lambda1520Mass * Lambda1520Mass;
    return (mtScaledPtSquared > 0.f) ? std::sqrt(mtScaledPtSquared) : -1.f;
  }

  template <bool IsAcceptedEvent, typename McPartType>
  void fillSignalLossProxySpecies(McPartType const& part, float centrality)
  {
    const int absPdgCode = std::abs(part.pdgCode());
    const float pt = part.pt();
    if (pt <= 0.f) {
      return;
    }

    if (absPdgCode == lambda1520analysispo_analysis::PdgKaon) {
      const float mtScaledPt = computeMtScaledPt(pt, KaonMass);
      const float weight = (mcCfg.applyMtScaleWeight && mtScaledPt > 0.f) ? mtScaledPt / pt : 1.f;
      if constexpr (IsAcceptedEvent) {
        histos.fill(HIST("SignalLoss/Kaon/hAcceptedEventTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Kaon/hAcceptedEventMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      } else {
        histos.fill(HIST("SignalLoss/Kaon/hGeneratedTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Kaon/hGeneratedMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      }
    } else if (absPdgCode == lambda1520analysispo_analysis::PdgProton) {
      const float mtScaledPt = computeMtScaledPt(pt, ProtonMass);
      const float weight = (mcCfg.applyMtScaleWeight && mtScaledPt > 0.f) ? mtScaledPt / pt : 1.f;
      if constexpr (IsAcceptedEvent) {
        histos.fill(HIST("SignalLoss/Proton/hAcceptedEventTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Proton/hAcceptedEventMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      } else {
        histos.fill(HIST("SignalLoss/Proton/hGeneratedTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Proton/hGeneratedMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      }
    } else if (absPdgCode == lambda1520analysispo_analysis::PdgLambda0) {
      const float mtScaledPt = computeMtScaledPt(pt, LambdaMass);
      const float weight = (mcCfg.applyMtScaleWeight && mtScaledPt > 0.f) ? mtScaledPt / pt : 1.f;
      if constexpr (IsAcceptedEvent) {
        histos.fill(HIST("SignalLoss/Lambda0/hAcceptedEventTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Lambda0/hAcceptedEventMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      } else {
        histos.fill(HIST("SignalLoss/Lambda0/hGeneratedTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Lambda0/hGeneratedMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      }
    } else if (absPdgCode == lambda1520analysispo_analysis::PdgXiMinus) {
      const float mtScaledPt = computeMtScaledPt(pt, XiMass);
      const float weight = (mcCfg.applyMtScaleWeight && mtScaledPt > 0.f) ? mtScaledPt / pt : 1.f;
      if constexpr (IsAcceptedEvent) {
        histos.fill(HIST("SignalLoss/Xi/hAcceptedEventTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Xi/hAcceptedEventMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      } else {
        histos.fill(HIST("SignalLoss/Xi/hGeneratedTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Xi/hGeneratedMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      }
    } else if (absPdgCode == lambda1520analysispo_analysis::PdgOmegaMinus) {
      const float mtScaledPt = computeMtScaledPt(pt, OmegaMass);
      const float weight = (mcCfg.applyMtScaleWeight && mtScaledPt > 0.f) ? mtScaledPt / pt : 1.f;
      if constexpr (IsAcceptedEvent) {
        histos.fill(HIST("SignalLoss/Omega/hAcceptedEventTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Omega/hAcceptedEventMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      } else {
        histos.fill(HIST("SignalLoss/Omega/hGeneratedTruePtVsCentrality"), pt, centrality);
        if (mtScaledPt > 0.f) {
          histos.fill(HIST("SignalLoss/Omega/hGeneratedMtScaledPtVsCentrality"), mtScaledPt, centrality, weight);
        }
      }
    }
  }

  Filter acceptanceFilter = (nabs(aod::track::eta) < trackCuts.cfgCutEta && nabs(aod::track::pt) > trackCuts.cfgCutPt);

  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::Mults>;
  using TrackBase = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension,
                              aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                              aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using TrackCandidates = soa::Filtered<TrackBase>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<TrackBase, aod::McTrackLabels>>;
  using EventCandidatesMC = soa::Join<EventCandidates, aod::McCollisionLabels>;
  using McCollisionsCent = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Preslice<TrackCandidates> perCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> perColMC = aod::track::collisionId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  using BinningTypeFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;

  void processSameEvent(EventCandidates const& collisions,
                        TrackCandidates const& tracks,
                        aod::BCs const&)
  {
    for (auto const& collision : collisions) {
      if (!passesEventSelectionImpl<DataCutFlowTag>(collision)) {
        continue;
      }
      const float centrality = getCentrality(collision);
      const int occupancy = collision.trackOccupancyInTimeRange();

      histos.fill(HIST("EventQA/Data/hVertexZ"), collision.posZ());
      histos.fill(HIST("EventQA/Data/hPVContributors"), collision.numContrib());
      histos.fill(HIST("EventQA/Data/hCentrality"), centrality);
      histos.fill(HIST("EventQA/Data/hOccupancy"), occupancy);
      histos.fill(HIST("EventQA/Data/hCentralityVsVertexZ"), collision.posZ(), centrality);
      histos.fill(HIST("EventQA/Data/hCentralityVsOccupancy"), occupancy, centrality);

      const int64_t collisionIndex = collision.globalIndex();
      auto tracksOfCollision = tracks.sliceBy(perCol, collisionIndex);
      tracksOfCollision.bindExternalIndices(&tracks);

      if (qa.cfgFillMultQA) {
        histos.fill(HIST("EventQA/Data/hCentralityVsFilteredTracksVsPVTracks"), centrality, tracksOfCollision.size(), collision.multNTracksPV());
      }

      if (qa.cfgFillTrackQA) {
        for (auto const& track : tracksOfCollision) {
          const float mom = RecoDecay::p(track.px(), track.py(), track.pz());
          histos.fill(HIST("TrackQA/Data/PreSelection/Proton/hNSigmaTPCvsP"), mom, track.tpcNSigmaPr());
          histos.fill(HIST("TrackQA/Data/PreSelection/Kaon/hNSigmaTPCvsP"), mom, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            histos.fill(HIST("TrackQA/Data/PreSelection/Proton/hNSigmaTOFvsP"), mom, track.tofNSigmaPr());
            histos.fill(HIST("TrackQA/Data/PreSelection/Proton/hNSigmaTOFvsNSigmaTPC"), track.tpcNSigmaPr(), track.tofNSigmaPr());
            histos.fill(HIST("TrackQA/Data/PreSelection/Kaon/hNSigmaTOFvsP"), mom, track.tofNSigmaKa());
            histos.fill(HIST("TrackQA/Data/PreSelection/Kaon/hNSigmaTOFvsNSigmaTPC"), track.tpcNSigmaKa(), track.tofNSigmaKa());
          }

          histos.fill(HIST("TrackQA/Data/AllTracks/hEta"), track.eta());
          histos.fill(HIST("TrackQA/Data/AllTracks/hPt"), track.pt());
          histos.fill(HIST("TrackQA/Data/AllTracks/hPhi"), track.phi());
          histos.fill(HIST("TrackQA/Data/AllTracks/hEtaPhi"), track.eta(), track.phi());

          if (!passesBasicTrackSelection(track)) {
            continue;
          }
          if (passesProtonPID(track, mom) && passesProtonDCASelection(track)) {
            fillSelectedProtonQA(track, mom, centrality);
          }
          if (passesKaonPID(track, mom) && passesKaonDCASelection(track)) {
            fillSelectedKaonQA(track, mom, centrality);
          }
        }
      }

      buildProtonKaonCandidatePools<false>(tracksOfCollision, mProtonCandidates, mKaonCandidates);
      fillInvariantMassFromCandidatePools<false>(mProtonCandidates, mKaonCandidates, centrality, qa.cfgFillPairQA);
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processSameEvent, "Process same event (data)", true);

  template <bool IsMC, typename PairType>
  void loopOverMixedEventPairs(PairType& mixingPair)
  {
    for (auto const& [c1, tracks1, c2, tracks2] : mixingPair) {
      if (!passesEventSelection(c1) || !passesEventSelection(c2)) {
        continue;
      }
      const float centrality = getCentrality(c1);

      if constexpr (IsMC) {
        histos.fill(HIST("MixingQA/MCReco/hCollisionPairs"), 1.f);
        histos.fill(HIST("MixingQA/MCReco/hCollisionPairsVsCentrality"), centrality);
      } else {
        if (qa.cfgFillMixQA) {
          histos.fill(HIST("MixingQA/Data/hPoolVertexZ"), c1.posZ());
          histos.fill(HIST("MixingQA/Data/hPoolCentrality"), centrality);
          histos.fill(HIST("MixingQA/Data/hPoolVertexZVsCentrality"), c1.posZ(), centrality);
          histos.fill(HIST("MixingQA/Data/hCollisionPairs"), 1.f);
        }
      }

      std::vector<DaughterTrackCandidate> protonPool1, kaonPool1, protonPool2, kaonPool2;
      buildProtonKaonCandidatePools<false>(tracks1, protonPool1, kaonPool1);
      buildProtonKaonCandidatePools<false>(tracks2, protonPool2, kaonPool2);
      fillInvariantMassFromCandidatePools<true, IsMC>(protonPool1, kaonPool2, centrality);
      fillInvariantMassFromCandidatePools<true, IsMC>(protonPool2, kaonPool1, centrality);
    }
  }

  template <bool IsMC, typename BinningType>
  void runMixingWithBinning(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    BinningType binning{{axes.cfgVtxAxis, axes.cfgCentAxis}, true};
    const int nMixedEvents = mixingCfg.cfgNoMixedEvents;
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<EventCandidates, TrackCandidates, BinningType> mixingPair{binning, nMixedEvents, -1, collisions, tracksTuple, &cache};
    loopOverMixedEventPairs<IsMC>(mixingPair);
  }

  template <bool IsMC>
  void runMixedEventForSelectedEstimator(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    switch (evSel.cfgCentralityEstimator) {
      case lambda1520analysispo_analysis::kFT0M:
        runMixingWithBinning<IsMC, BinningTypeFT0M>(collisions, tracks);
        break;
      case lambda1520analysispo_analysis::kFT0A:
        runMixingWithBinning<IsMC, BinningTypeFT0A>(collisions, tracks);
        break;
      case lambda1520analysispo_analysis::kFT0C:
        runMixingWithBinning<IsMC, BinningTypeFT0C>(collisions, tracks);
        break;
      case lambda1520analysispo_analysis::kFV0A:
        runMixingWithBinning<IsMC, BinningTypeFV0A>(collisions, tracks);
        break;
      default:
        runMixingWithBinning<IsMC, BinningTypeFT0M>(collisions, tracks);
        break;
    }
  }

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    runMixedEventForSelectedEstimator<false>(collisions, tracks);
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMixedEvent, "Process mixed event (data)", true);

  void processMCEventMixing(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    runMixedEventForSelectedEstimator<true>(collisions, tracks);
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMCEventMixing, "Process MC reconstructed event mixing (reuses the mixing pairs; no truth matching)", false);

  void processMC(EventCandidatesMC const& collisions,
                 aod::BCs const&,
                 TrackCandidatesMC const& tracks,
                 aod::McParticles const& mcParticles)
  {
    std::unordered_set<int64_t> mcCollisionsWithEfficiencyQA;

    for (auto const& collision : collisions) {
      if (!passesEventSelectionImpl<McRecoCutFlowTag>(collision)) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      const float centrality = getCentrality(collision);

      histos.fill(HIST("EventQA/MCReco/hAcceptedEventsVsCentrality"), centrality);
      histos.fill(HIST("EventQA/MCReco/hVertexZ"), collision.posZ());
      histos.fill(HIST("EventQA/MCReco/hPVContributors"), collision.numContrib());

      const int64_t collisionIndex = collision.globalIndex();
      auto tracksOfCollision = tracks.sliceBy(perColMC, collisionIndex);
      tracksOfCollision.bindExternalIndices(&tracks);
      tracksOfCollision.bindExternalIndices(&mcParticles);

      for (auto const& track : tracksOfCollision) {
        histos.fill(HIST("TrackQA/MCReco/AllTracks/hEta"), track.eta());
        histos.fill(HIST("TrackQA/MCReco/AllTracks/hPt"), track.pt());
        histos.fill(HIST("TrackQA/MCReco/AllTracks/hPhi"), track.phi());
        histos.fill(HIST("TrackQA/MCReco/AllTracks/hEtaPhi"), track.eta(), track.phi());

        if (qa.cfgFillTrackQA) {
          const float mom = RecoDecay::p(track.px(), track.py(), track.pz());
          histos.fill(HIST("TrackQA/MCReco/PreSelection/Proton/hNSigmaTPCvsP"), mom, track.tpcNSigmaPr());
          histos.fill(HIST("TrackQA/MCReco/PreSelection/Kaon/hNSigmaTPCvsP"), mom, track.tpcNSigmaKa());
          if (track.hasTOF()) {
            histos.fill(HIST("TrackQA/MCReco/PreSelection/Proton/hNSigmaTOFvsP"), mom, track.tofNSigmaPr());
            histos.fill(HIST("TrackQA/MCReco/PreSelection/Proton/hNSigmaTOFvsNSigmaTPC"), track.tpcNSigmaPr(), track.tofNSigmaPr());
            histos.fill(HIST("TrackQA/MCReco/PreSelection/Kaon/hNSigmaTOFvsP"), mom, track.tofNSigmaKa());
            histos.fill(HIST("TrackQA/MCReco/PreSelection/Kaon/hNSigmaTOFvsNSigmaTPC"), track.tpcNSigmaKa(), track.tofNSigmaKa());
          }
          if (passesBasicTrackSelection(track)) {
            if (passesProtonPID(track, mom) && passesProtonDCASelection(track)) {
              fillSelectedProtonQAMC(track, mom, centrality);
            }
            if (passesKaonPID(track, mom) && passesKaonDCASelection(track)) {
              fillSelectedKaonQAMC(track, mom, centrality);
            }
          }
        }
      }

      buildProtonKaonCandidatePools<true>(tracksOfCollision, mProtonCandidates, mKaonCandidates);

      for (auto const& pr : mProtonCandidates) {
        if (pr.mcLambda1520MotherIndex >= 0 && std::abs(pr.mcPdgCode) == lambda1520analysispo_analysis::PdgProton) {
          if (pr.mcLambda1520MotherPdgCode > 0) {
            histos.fill(HIST("MCRecoQA/DaughterPt/hProtonFromLambda1520"), pr.pt);
          } else {
            histos.fill(HIST("MCRecoQA/DaughterPt/hProtonFromAntiLambda1520"), pr.pt);
          }
        }
      }
      for (auto const& ka : mKaonCandidates) {
        if (ka.mcLambda1520MotherIndex >= 0 && std::abs(ka.mcPdgCode) == lambda1520analysispo_analysis::PdgKaon) {
          if (ka.mcLambda1520MotherPdgCode > 0) {
            histos.fill(HIST("MCRecoQA/DaughterPt/hKaonFromLambda1520"), ka.pt);
          } else {
            histos.fill(HIST("MCRecoQA/DaughterPt/hKaonFromAntiLambda1520"), ka.pt);
          }
        }
      }

      for (auto const& pr : mProtonCandidates) {
        for (auto const& ka : mKaonCandidates) {
          if (pr.globalIndex == ka.globalIndex) {
            continue;
          }
          if (pr.sign * ka.sign >= 0) {
            continue;
          }
          float pairMass = 0.f, pairPt = 0.f, pairY = 0.f;
          if (!buildProtonKaonPair(pr, ka, pairMass, pairPt, pairY)) {
            continue;
          }

          fillInvariantMassMCRecoAllPairs(pr.sign, pairPt, pairMass, centrality, pairY);

          histos.fill(HIST("CutFlow/MCReco/hTruthMatching"), 0.f);
          if (pr.hasMcLabel) {
            histos.fill(HIST("CutFlow/MCReco/hTruthMatching"), 1.f);
          }
          if (ka.hasMcLabel) {
            histos.fill(HIST("CutFlow/MCReco/hTruthMatching"), 2.f);
          }
          if (!pr.hasMcLabel || !ka.hasMcLabel) {
            continue;
          }
          if (std::abs(pr.mcPdgCode) != lambda1520analysispo_analysis::PdgProton) {
            continue;
          }
          histos.fill(HIST("CutFlow/MCReco/hTruthMatching"), 3.f);
          if (std::abs(ka.mcPdgCode) != lambda1520analysispo_analysis::PdgKaon) {
            continue;
          }
          histos.fill(HIST("CutFlow/MCReco/hTruthMatching"), 4.f);
          if (pr.mcLambda1520MotherIndex < 0 || pr.mcLambda1520MotherIndex != ka.mcLambda1520MotherIndex) {
            continue;
          }
          histos.fill(HIST("CutFlow/MCReco/hTruthMatching"), 5.f);

          if (mcCfg.applyGeneratedEtaCut) {
            const float pairEta = RecoDecay::eta(std::array{pr.px + ka.px, pr.py + ka.py, pr.pz + ka.pz});
            if (std::abs(pairEta) > mcCfg.generatedEtaMax) {
              continue;
            }
          }

          const bool motherIsParticle = pr.mcLambda1520MotherPdgCode > 0;
          fillInvariantMassMCRecoTruthMatched(pr.mcLambda1520MotherPdgCode, pairPt, pairMass, centrality, pairY);

          const float massResolution = pairMass - computeGeneratedMass(mcParticles.iteratorAt(pr.mcLambda1520MotherIndex));
          if (motherIsParticle) {
            histos.fill(HIST("InvMass/MCReco/MassResolution/hLambda1520"), massResolution, pairPt, centrality);
          } else {
            histos.fill(HIST("InvMass/MCReco/MassResolution/hAntiLambda1520"), massResolution, pairPt, centrality);
          }

          if (isSpinAlignmentActive()) {
            const DecayKinematics decay = makeDecayKinematics(pr, ka);
            fillCosThetaStarHistograms<CosThetaStarSample::TruthMatched>(pairMass, pairPt, decay, centrality, motherIsParticle);

            const float cosTheta = computeCosThetaStar(decay);
            if (motherIsParticle) {
              histos.fill(HIST("CosThetaStar/MCReco/MassResolution/hLambda1520"), massResolution, pairPt, centrality, cosTheta);
            } else {
              histos.fill(HIST("CosThetaStar/MCReco/MassResolution/hAntiLambda1520"), massResolution, pairPt, centrality, cosTheta);
            }
          }

          histos.fill(HIST("PairQA/MCReco/hTruthMatchedRapidityVsPt"), pairPt, pairY);
          histos.fill(HIST("TrackQA/MCReco/TruthMatched/Proton/hDCAxyVsPt"), pr.pt, pr.dcaXY);
          histos.fill(HIST("TrackQA/MCReco/TruthMatched/Proton/hDCAzVsPt"), pr.pt, pr.dcaZ);
          histos.fill(HIST("TrackQA/MCReco/TruthMatched/Proton/hNSigmaTPCvsPt"), pr.pt, pr.tpcNSigmaPr);
          if (pr.hasTOF) {
            histos.fill(HIST("TrackQA/MCReco/TruthMatched/Proton/hNSigmaTOFvsPt"), pr.pt, pr.tofNSigmaPr);
          }
          histos.fill(HIST("TrackQA/MCReco/TruthMatched/Kaon/hDCAxyVsPt"), ka.pt, ka.dcaXY);
          histos.fill(HIST("TrackQA/MCReco/TruthMatched/Kaon/hDCAzVsPt"), ka.pt, ka.dcaZ);
          histos.fill(HIST("TrackQA/MCReco/TruthMatched/Kaon/hNSigmaTPCvsPt"), ka.pt, ka.tpcNSigmaKa);
          if (ka.hasTOF) {
            histos.fill(HIST("TrackQA/MCReco/TruthMatched/Kaon/hNSigmaTOFvsPt"), ka.pt, ka.tofNSigmaKa);
          }
        }
      }

      if (mcCollisionsWithEfficiencyQA.insert(collision.mcCollisionId()).second) {
        fillGeneratedToReconstructedCutFlow(tracksOfCollision, mcParticles, collision.mcCollisionId());
      }
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMC, "Process MC reconstructed (truth matching)", false);

  void processMCGen(McCollisionsCent::iterator const& mcCollision,
                    aod::McParticles const& mcParticles,
                    soa::SmallGroups<EventCandidatesMC> const& collisions)
  {
    histos.fill(HIST("CutFlow/MCGen/hEventSelection"), 0.f);

    if (std::abs(mcCollision.posZ()) > evSel.cfgCutVertex) {
      return;
    }
    histos.fill(HIST("CutFlow/MCGen/hEventSelection"), 1.f);

    if (evSel.cfgApplyMCTruthINELgt0 && !isTruthInelGt0(mcParticles)) {
      return;
    }
    histos.fill(HIST("CutFlow/MCGen/hEventSelection"), 2.f);

    float genCentrality = getTruthCentrality(mcCollision, mcParticles);
    if (genCentrality == InvalidCentrality) {
      return;
    }
    histos.fill(HIST("CutFlow/MCGen/hEventSelection"), 3.f);
    histos.fill(HIST("CutFlow/MCGen/hEventSelection"), 4.f);

    const AcceptedRecoCollision acceptedReco = findFirstAcceptedRecoCollision(collisions);
    const bool hasAcceptedReco = acceptedReco.found;
    if (hasAcceptedReco) {
      if (!evSel.cfgUseMCTruthCentrality) {
        genCentrality = acceptedReco.centrality;
      }
      histos.fill(HIST("CutFlow/MCGen/hEventSelection"), 5.f);
    }

    for (auto const& mcParticle : mcParticles) {
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 0.f);

      if (!mcParticle.producedByGenerator()) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 1.f);

      const int pdg = mcParticle.pdgCode();
      if (std::abs(pdg) != mLambda1520AbsPdgCode) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 2.f);

      auto genProton = mcParticle;
      auto genKaon = mcParticle;
      if (!findProtonKaonDaughters(mcParticle, genProton, genKaon)) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 3.f);

      const float genPt = mcParticle.pt();
      const float genY = mcParticle.y();
      histos.fill(HIST("MCGenQA/Kinematics/hEtaVsPtBeforeCuts"), mcParticle.eta(), genPt);
      histos.fill(HIST("MCGenQA/Kinematics/hPhiVsRapidityBeforeCuts"), mcParticle.phi(), genY);

      if (genY < trackCuts.cfgMinY || genY > trackCuts.cfgMaxY) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 4.f);

      const float maxTrackEta = trackCuts.cfgCutEta;
      if (mcCfg.applyGeneratedDaughterEtaCut && (std::abs(genProton.eta()) >= maxTrackEta || std::abs(genKaon.eta()) >= maxTrackEta)) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 5.f);

      if (mcCfg.applyGeneratedEtaCut && std::abs(mcParticle.eta()) > mcCfg.generatedEtaMax) {
        continue;
      }
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 6.f);
      histos.fill(HIST("CutFlow/MCGen/hGeneratedLambda1520Selection"), 7.f);

      histos.fill(HIST("MCGenQA/Kinematics/hEtaVsPtAfterCuts"), mcParticle.eta(), genPt);
      histos.fill(HIST("MCGenQA/Kinematics/hPhiVsRapidityAfterCuts"), mcParticle.phi(), genY);
      histos.fill(HIST("MCGenQA/Kinematics/hRapidityVsPtAccepted"), genPt, genY);

      const float genMass = computeGeneratedMass(mcParticle);
      const float ptPr = genProton.pt();
      const float ptKa = genKaon.pt();

      if (pdg > 0) {
        histos.fill(HIST("InvMass/MCGen/Generated/hLambda1520"), genMass, genPt, genCentrality, genY);
        histos.fill(HIST("MCGenQA/DaughterPt/hProtonFromLambda1520"), ptPr);
        histos.fill(HIST("MCGenQA/DaughterPt/hKaonFromLambda1520"), ptKa);
      } else {
        histos.fill(HIST("InvMass/MCGen/Generated/hAntiLambda1520"), genMass, genPt, genCentrality, genY);
        histos.fill(HIST("MCGenQA/DaughterPt/hProtonFromAntiLambda1520"), ptPr);
        histos.fill(HIST("MCGenQA/DaughterPt/hKaonFromAntiLambda1520"), ptKa);
      }

      if (hasAcceptedReco) {
        if (pdg > 0) {
          histos.fill(HIST("InvMass/MCGen/GeneratedInAcceptedEvents/hLambda1520"), genMass, genPt, genCentrality, genY);
        } else {
          histos.fill(HIST("InvMass/MCGen/GeneratedInAcceptedEvents/hAntiLambda1520"), genMass, genPt, genCentrality, genY);
        }
      }
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMCGen, "Process MC truth: generated-level Lambda(1520) spectra and QA", false);

  void processMCGenSpinAlignment(McCollisionsCent::iterator const& mcCollision,
                                 aod::McParticles const& mcParticles,
                                 soa::SmallGroups<EventCandidatesMC> const& collisions)
  {
    if (std::abs(mcCollision.posZ()) > evSel.cfgCutVertex) {
      return;
    }
    if (evSel.cfgApplyMCTruthINELgt0 && !isTruthInelGt0(mcParticles)) {
      return;
    }
    float genCentrality = getTruthCentrality(mcCollision, mcParticles);
    if (genCentrality == InvalidCentrality) {
      return;
    }

    const AcceptedRecoCollision acceptedReco = findFirstAcceptedRecoCollision(collisions);
    const bool hasAcceptedReco = acceptedReco.found;
    if (hasAcceptedReco && !evSel.cfgUseMCTruthCentrality) {
      genCentrality = acceptedReco.centrality;
    }

    for (auto const& mcParticle : mcParticles) {
      if (!mcParticle.producedByGenerator()) {
        continue;
      }
      const int pdg = mcParticle.pdgCode();
      if (std::abs(pdg) != mLambda1520AbsPdgCode) {
        continue;
      }
      auto genProton = mcParticle;
      auto genKaon = mcParticle;
      if (!findProtonKaonDaughters(mcParticle, genProton, genKaon)) {
        continue;
      }

      const float genPt = mcParticle.pt();
      const float genY = mcParticle.y();
      if (genY < trackCuts.cfgMinY || genY > trackCuts.cfgMaxY) {
        continue;
      }

      const float maxTrackEta = trackCuts.cfgCutEta;
      if (mcCfg.applyGeneratedDaughterEtaCut && (std::abs(genProton.eta()) >= maxTrackEta || std::abs(genKaon.eta()) >= maxTrackEta)) {
        continue;
      }
      if (mcCfg.applyGeneratedEtaCut && std::abs(mcParticle.eta()) > mcCfg.generatedEtaMax) {
        continue;
      }

      const float genMass = computeGeneratedMass(mcParticle);
      const float cosTheta = computeCosThetaStar(makeGeneratedDecayKinematics(mcParticle, genProton, genKaon));

      if (pdg > 0) {
        histos.fill(HIST("CosThetaStar/MCGen/Generated/hLambda1520"), genMass, genPt, genCentrality, cosTheta);
      } else {
        histos.fill(HIST("CosThetaStar/MCGen/Generated/hAntiLambda1520"), genMass, genPt, genCentrality, cosTheta);
      }

      if (hasAcceptedReco) {
        if (pdg > 0) {
          histos.fill(HIST("CosThetaStar/MCGen/GeneratedInAcceptedEvents/hLambda1520"), genMass, genPt, genCentrality, cosTheta);
        } else {
          histos.fill(HIST("CosThetaStar/MCGen/GeneratedInAcceptedEvents/hAntiLambda1520"), genMass, genPt, genCentrality, cosTheta);
        }
      }
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processMCGenSpinAlignment, "Process MC truth: generated-level Lambda(1520) spin-alignment cos(theta*) with same cuts as processMCGen", false);

  void processEventFactor(EventCandidatesMC const& collisions,
                          McCollisionsCent const& mcCollisions,
                          aod::McParticles const& mcParticles)
  {
    std::unordered_map<int64_t, int> nRecoCollisionsPerMc;
    std::unordered_map<int64_t, float> acceptedMcCollisionRecoCent;
    for (auto const& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      ++nRecoCollisionsPerMc[collision.mcCollisionId()];
      if (!passesEventSelection(collision)) {
        continue;
      }
      acceptedMcCollisionRecoCent.emplace(collision.mcCollisionId(), getCentrality(collision));

      const auto& mcCollision = collision.mcCollision_as<McCollisionsCent>();
      auto mcParticlesOfRecoCollision = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());
      if (evSel.cfgApplyMCTruthINELgt0 && !isTruthInelGt0(mcParticlesOfRecoCollision)) {
        continue;
      }
      const float truthCentrality = getTruthCentrality(mcCollision, mcParticlesOfRecoCollision);
      if (truthCentrality == InvalidCentrality) {
        continue;
      }
      histos.fill(HIST("EventLoss/hAcceptedRecoEventsVsTruthCentrality"), truthCentrality);
    }

    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("CutFlow/EventLoss/hGeneratedEventSelection"), 0.f);

      auto mcParticlesOfGenCollision = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());
      const bool isInsideVertexRange = std::abs(mcCollision.posZ()) <= evSel.cfgCutVertex;
      const bool truthInelGt0 = isTruthInelGt0(mcParticlesOfGenCollision);
      const float truthCentrality = getTruthCentrality(mcCollision, mcParticlesOfGenCollision);

      if (!isInsideVertexRange) {
        continue;
      }
      histos.fill(HIST("CutFlow/EventLoss/hGeneratedEventSelection"), 1.f);

      if (evSel.cfgApplyMCTruthINELgt0 && !truthInelGt0) {
        continue;
      }
      histos.fill(HIST("CutFlow/EventLoss/hGeneratedEventSelection"), 2.f);

      if (truthCentrality == InvalidCentrality) {
        continue;
      }
      histos.fill(HIST("EventLoss/hGeneratedEventsVsTruthCentrality"), truthCentrality);

      const auto recoCollisionCountEntry = nRecoCollisionsPerMc.find(mcCollision.globalIndex());
      const int nRecoTotal = (recoCollisionCountEntry != nRecoCollisionsPerMc.end()) ? recoCollisionCountEntry->second : 0;
      histos.fill(HIST("EventLoss/hRecoCollisionMultiplicity"), nRecoTotal);

      if (nRecoTotal > 0) {
        histos.fill(HIST("CutFlow/EventLoss/hGeneratedEventSelection"), 3.f);
      }

      const auto acceptedRecoCentralityEntry = acceptedMcCollisionRecoCent.find(mcCollision.globalIndex());
      if (acceptedRecoCentralityEntry != acceptedMcCollisionRecoCent.end()) {
        histos.fill(HIST("CutFlow/EventLoss/hGeneratedEventSelection"), 4.f);
        histos.fill(HIST("CutFlow/EventLoss/hGeneratedEventSelection"), 5.f);
        histos.fill(HIST("EventLoss/hGeneratedAcceptedEventsVsTruthCentrality"), truthCentrality);
        histos.fill(HIST("EventLoss/hGeneratedAcceptedEventsVsRecoCentrality"), acceptedRecoCentralityEntry->second);
      }
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processEventFactor, "Process Event Loss / Event Splitting factors", false);

  void processSignalLoss(McCollisionsCent::iterator const& mcCollision,
                         aod::McParticles const& mcParticles,
                         soa::SmallGroups<EventCandidatesMC> const& collisions)
  {
    if (std::abs(mcCollision.posZ()) > evSel.cfgCutVertex) {
      return;
    }
    if (evSel.cfgApplyMCTruthINELgt0 && !isTruthInelGt0(mcParticles)) {
      return;
    }
    const float centrality = getTruthCentrality(mcCollision, mcParticles);
    if (centrality == InvalidCentrality) {
      return;
    }

    const bool hasAcceptedReco = findFirstAcceptedRecoCollision(collisions).found;

    const float minRapidity = trackCuts.cfgMinY;
    const float maxRapidity = trackCuts.cfgMaxY;
    const float maxAbsGeneratedEta = mcCfg.generatedEtaMax;

    for (auto const& part : mcParticles) {
      const int absPdgCode = std::abs(part.pdgCode());
      const bool isLambda1520 = absPdgCode == mLambda1520AbsPdgCode;
      const bool isProxySpecies = absPdgCode == lambda1520analysispo_analysis::PdgKaon ||
                                  absPdgCode == lambda1520analysispo_analysis::PdgProton ||
                                  absPdgCode == lambda1520analysispo_analysis::PdgLambda0 ||
                                  absPdgCode == lambda1520analysispo_analysis::PdgXiMinus ||
                                  absPdgCode == lambda1520analysispo_analysis::PdgOmegaMinus;
      if (!isLambda1520 && !isProxySpecies) {
        continue;
      }

      if (part.y() < minRapidity || part.y() > maxRapidity) {
        continue;
      }
      if (mcCfg.applyGeneratedEtaCut && std::abs(part.eta()) > maxAbsGeneratedEta) {
        continue;
      }

      if (isLambda1520) {
        if (!part.producedByGenerator()) {
          continue;
        }
        auto genProton = part;
        auto genKaon = part;
        if (!findProtonKaonDaughters(part, genProton, genKaon)) {
          continue;
        }
        if (part.pdgCode() > 0) {
          histos.fill(HIST("SignalLoss/Lambda1520/hGeneratedPtVsCentrality"), part.pt(), centrality);
          if (hasAcceptedReco) {
            histos.fill(HIST("SignalLoss/Lambda1520/hAcceptedEventPtVsCentrality"), part.pt(), centrality);
          }
        } else {
          histos.fill(HIST("SignalLoss/AntiLambda1520/hGeneratedPtVsCentrality"), part.pt(), centrality);
          if (hasAcceptedReco) {
            histos.fill(HIST("SignalLoss/AntiLambda1520/hAcceptedEventPtVsCentrality"), part.pt(), centrality);
          }
        }
        continue;
      }

      if (!part.isPhysicalPrimary()) {
        continue;
      }
      fillSignalLossProxySpecies<false>(part, centrality);
      if (hasAcceptedReco) {
        fillSignalLossProxySpecies<true>(part, centrality);
      }
    }
  }
  PROCESS_SWITCH(Lambda1520Analysispo, processSignalLoss, "Process signal loss (K, p, Lambda, Xi, Omega proxies and Lambda(1520))", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<Lambda1520Analysispo>(context)};
}
