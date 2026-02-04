/// \file lambdaInvMassTest.cxx
/// \brief Lambda analysis task using derived data as a test
///
/// \author Cicero Domenico Muncinelli <cicero.domenico.muncinelli@cern.ch>, Campinas State University, Brazil
//
// Lambda Invariant Mass test task
// ================
//
// This code loops over a V0Cores table and produces some
// standard analysis output. It is meant to be run over
// derived data.
// This is NOT meant to be in the ALICE O2 repository, 
// ever! It is just a small test for me to get up to 
// speed with O2 analyses!
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    cicero.domenico.muncinelli@cern.ch
//
// This code is heavily based on the derivedlambdakzeroanalysis.cxx code!

#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h" // This actually needs ONNX to be installed in the system! Or, at least, you should link the libraries properly in the CMakeLists.txt

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using namespace o2::aod::rctsel;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using DauMCTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackMCIds, aod::DauTrackTPCPIDs>;
// using V0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0K0ShortMLScores>;
using V0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores>;
// using V0McCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0MCCores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0MCCollRefs>;
// using V0McCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0K0ShortMLScores>;
using V0McCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores>;

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

// bool doPlainTopoQA = true;

struct lambdaInvMassTest{
    HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

    bool isRun3 = true;

    // master analysis switches
    Configurable<bool> analyseK0Short{"analyseK0Short", false, "process K0Short-like candidates"};
    Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
    Configurable<bool> analyseAntiLambda{"analyseAntiLambda", false, "process AntiLambda-like candidates"};
    Configurable<bool> calculateFeeddownMatrix{"calculateFeeddownMatrix", true, "fill feeddown matrix if MC"};

    Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};
    Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
    Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0C, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFT0CVariant1, 3:CentMFT, 4:CentNGlobal, 5:CentFV0A)"};

    Configurable<bool> doEventQA{"doEventQA", false, "do event QA histograms"};
    Configurable<bool> doCompleteTopoQA{"doCompleteTopoQA", false, "do topological variable QA histograms"};
    Configurable<bool> doTPCQA{"doTPCQA", false, "do TPC QA histograms"};
    Configurable<bool> doTOFQA{"doTOFQA", false, "do TOF QA histograms"};
    Configurable<int> doDetectPropQA{"doDetectPropQA", 0, "do Detector/ITS map QA: 0: no, 1: 4D, 2: 5D with mass; 3: plain in 3D"};
    Configurable<bool> doEtaPhiQA{"doEtaPhiQA", false, "do Eta/Phi QA histograms"};

    Configurable<bool> doPlainTopoQA{"doPlainTopoQA", true, "do simple 1D QA of candidates"};
    Configurable<float> qaMinPt{"qaMinPt", 0.0f, "minimum pT for QA plots"};
    Configurable<float> qaMaxPt{"qaMaxPt", 1000.0f, "maximum pT for QA plots"};
    Configurable<bool> qaCentrality{"qaCentrality", false, "qa centrality flag: check base raw values"};

    // for MC
    Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"}; // Will not do MC, so can keep this on regardless
    Configurable<bool> doTreatPiToMuon{"doTreatPiToMuon", false, "Take pi decay into muon into account in MC"};
    Configurable<bool> doCollisionAssociationQA{"doCollisionAssociationQA", true, "check collision association"};

    //   // Defining a configurable axis for easier manipulation later on:
    // ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f,
    // 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f,
    // 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f,
    // 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f,
    // 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
    // // ConfigurableAxis axisInvMassLambda{"axisInvMassLambda", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f,
    // // 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f,
    // // 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f,
    // // 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f,
    // // 30.0f, 35.0f, 40.0f, 50.0f}, "Invariant mass axis for Lambda"};

    // My own version of these axes:
    // Configurable<int> nBinsInvMass{"nBinsInvMass", 100, "Number of bins of invariant mass axis"};
    // Configurable<float> minInvMass{"minInvMass", 0.7f, "Lower bound of invariant mass axis"};
    // Configurable<float> maxInvMass{"maxInvMass", 1.3f, "Upper bound of invariant mass axis"};
        // From David's code:
    // ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};

    // Defining a configurable object group for the event selections:
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

        // Run 2 specific event selections
        Configurable<bool> requireSel7{"requireSel7", true, "require sel7 event selection (Run 2 only: event selection decision based on V0A & V0C)"};
        Configurable<bool> requireINT7{"requireINT7", true, "require INT7 trigger selection (Run 2 only)"};
        Configurable<bool> rejectIncompleteDAQ{"rejectIncompleteDAQ", true, "reject events with incomplete DAQ (Run 2 only)"};
        Configurable<bool> requireConsistentSPDAndTrackVtx{"requireConsistentSPDAndTrackVtx", true, "reject events with inconsistent in SPD and Track vertices (Run 2 only)"};
        Configurable<bool> rejectPileupFromSPD{"rejectPileupFromSPD", true, "reject events with pileup according to SPD vertexer (Run 2 only)"};
        Configurable<bool> rejectV0PFPileup{"rejectV0PFPileup", false, "reject events tagged as OOB pileup according to V0 past-future info (Run 2 only)"};
        Configurable<bool> rejectPileupInMultBins{"rejectPileupInMultBins", true, "reject events tagged as pileup according to multiplicity-differential pileup checks (Run 2 only)"};
        Configurable<bool> rejectPileupMV{"rejectPileupMV", true, "reject events tagged as pileup according to according to multi-vertexer (Run 2 only)"};
        Configurable<bool> rejectTPCPileup{"rejectTPCPileup", false, "reject events tagged as pileup according to pileup in TPC (Run 2 only)"};
        Configurable<bool> requireNoV0MOnVsOffPileup{"requireNoV0MOnVsOffPileup", false, "reject events tagged as OOB pileup according to online-vs-offline VOM correlation (Run 2 only)"};
        Configurable<bool> requireNoSPDOnVsOffPileup{"requireNoSPDOnVsOffPileup", false, "reject events tagged as pileup according to online-vs-offline SPD correlation (Run 2 only)"};
        Configurable<bool> requireNoSPDClsVsTklBG{"requireNoSPDClsVsTklBG", true, "reject events tagged as beam-gas and pileup according to cluster-vs-tracklet correlation (Run 2 only)"};

        Configurable<bool> useSPDTrackletsCent{"useSPDTrackletsCent", false, "Use SPD tracklets for estimating centrality? If not, use V0M-based centrality (Run 2 only)"};
    } eventSelections;

    static constexpr float DefaultLifetimeCuts[1][2] = {{30., 20.}};

    // Defining the configurable object that is going to be used: v0Selections, for selections later on
    struct : ConfigurableGroup {
        std::string prefix = "v0Selections"; // JSON group name
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
        Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {DefaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

        // // Additional selection on the AP plot (exclusive for K0Short)
        // // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
        // Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

        // Track quality
        Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
        Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
        Configurable<float> minTPCrowsOverFindableClusters{"minTPCrowsOverFindableClusters", -1, "minimum nbr of TPC crossed rows over findable clusters"};
        Configurable<float> minTPCfoundOverFindableClusters{"minTPCfoundOverFindableClusters", -1, "minimum nbr of found over findable TPC clusters"};
        Configurable<float> maxFractionTPCSharedClusters{"maxFractionTPCSharedClusters", 1e+09, "maximum fraction of TPC shared clusters"};
        Configurable<float> maxITSchi2PerNcls{"maxITSchi2PerNcls", 1e+09, "maximum ITS chi2 per clusters"};
        Configurable<float> maxTPCchi2PerNcls{"maxTPCchi2PerNcls", 1e+09, "maximum TPC chi2 per clusters"};
        Configurable<bool> skipTPConly{"skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
        Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
        Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};
        Configurable<bool> rejectPosITSafterburner{"rejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
        Configurable<bool> rejectNegITSafterburner{"rejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};
        Configurable<bool> requirePosITSafterburnerOnly{"requirePosITSafterburnerOnly", false, "require positive track formed out of afterburner ITS tracks"};
        Configurable<bool> requireNegITSafterburnerOnly{"requireNegITSafterburnerOnly", false, "require negative track formed out of afterburner ITS tracks"};
        Configurable<bool> rejectTPCsectorBoundary{"rejectTPCsectorBoundary", false, "reject tracks close to the TPC sector boundaries"};
        Configurable<std::string> phiLowCut{"phiLowCut", "0.06/x+pi/18.0-0.06", "Low azimuth cut parametrisation"};
        Configurable<std::string> phiHighCut{"phiHighCut", "0.1/x+pi/18.0+0.06", "High azimuth cut parametrisation"};

        // PID (TPC/TOF)
        Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
        Configurable<float> tofPidNsigmaCutLaPr{"tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
        Configurable<float> tofPidNsigmaCutLaPi{"tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
        Configurable<float> tofPidNsigmaCutK0Pi{"tofPidNsigmaCutK0Pi", 1e+6, "tofPidNsigmaCutK0Pi"};

        // PID (TOF)
        Configurable<float> maxDeltaTimeProton{"maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
        Configurable<float> maxDeltaTimePion{"maxDeltaTimePion", 1e+9, "check maximum allowed time"};
    } v0Selections;

    TF1* fPhiCutLow = new TF1("fPhiCutLow", v0Selections.phiLowCut.value.data(), 0, 100);
    TF1* fPhiCutHigh = new TF1("fPhiCutHigh", v0Selections.phiHighCut.value.data(), 0, 100);

    struct : ConfigurableGroup {
        std::string prefix = "rctConfigurations"; // JSON group name
        Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "", "Which detector condition requirements? (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
        Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "Include ZDC flags in the bit selection (for Pb-Pb only)"};
        Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
    } rctConfigurations;

    RCTFlagsChecker rctFlagsChecker{rctConfigurations.cfgRCTLabel.value};

    // Machine learning evaluation for pre-selection and corresponding information generation
    o2::ml::OnnxModel mlCustomModelK0Short;
    o2::ml::OnnxModel mlCustomModelLambda;
    o2::ml::OnnxModel mlCustomModelAntiLambda;
    o2::ml::OnnxModel mlCustomModelGamma;

    struct : ConfigurableGroup { // Kept the original configurable scores for K0Short and all else due to the line "if (lambdaScore > mlConfigurations.thresholdK0Short.value) (...)"
        std::string prefix = "mlConfigurations"; // JSON group name
        // ML classifiers: master flags to control whether we should use custom ML classifiers or the scores in the derived data
        Configurable<bool> useK0ShortScores{"useK0ShortScores", false, "use ML scores to select K0Short"};
        Configurable<bool> useLambdaScores{"useLambdaScores", false, "use ML scores to select Lambda"};
        Configurable<bool> useAntiLambdaScores{"useAntiLambdaScores", false, "use ML scores to select AntiLambda"};

        Configurable<bool> calculateK0ShortScores{"calculateK0ShortScores", false, "calculate K0Short ML scores"};
        Configurable<bool> calculateLambdaScores{"calculateLambdaScores", false, "calculate Lambda ML scores"};
        Configurable<bool> calculateAntiLambdaScores{"calculateAntiLambdaScores", false, "calculate AntiLambda ML scores"};

        // ML input for ML calculation
        Configurable<std::string> customModelPathCCDB{"customModelPathCCDB", "", "Custom ML Model path in CCDB"};
        Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
        Configurable<bool> loadCustomModelsFromCCDB{"loadCustomModelsFromCCDB", false, "Flag to enable or disable the loading of custom models from CCDB"};
        Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

        // Local paths for test purposes
        Configurable<std::string> localModelPathLambda{"localModelPathLambda", "Lambda_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
        Configurable<std::string> localModelPathAntiLambda{"localModelPathAntiLambda", "AntiLambda_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
        Configurable<std::string> localModelPathK0Short{"localModelPathK0Short", "KZeroShort_BDTModel.onnx", "(std::string) Path to the local .onnx file."};

        // Thresholds for choosing to populate V0Cores tables with pre-selections
        Configurable<float> thresholdLambda{"thresholdLambda", -1.0f, "Threshold to keep Lambda candidates"};
        Configurable<float> thresholdAntiLambda{"thresholdAntiLambda", -1.0f, "Threshold to keep AntiLambda candidates"};
        Configurable<float> thresholdK0Short{"thresholdK0Short", -1.0f, "Threshold to keep K0Short candidates"};
    } mlConfigurations;

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
        ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
        ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
        ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};
        ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles"};
        ConfigurableAxis axisIRBinning{"axisIRBinning", {500, 0, 50}, "Binning for the interaction rate (kHz)"};
        ConfigurableAxis axisMultFT0M{"axisMultFT0M", {500, 0.0f, +100000.0f}, "Multiplicity FT0M"};
        ConfigurableAxis axisMultFT0C{"axisMultFT0C", {500, 0.0f, +10000.0f}, "Multiplicity FT0C"};
        ConfigurableAxis axisMultFV0A{"axisMultFV0A", {500, 0.0f, +100000.0f}, "Multiplicity FV0A"};

        ConfigurableAxis axisRawCentrality{"axisRawCentrality", {VARIABLE_WIDTH, 0.000f, 52.320f, 75.400f, 95.719f, 115.364f, 135.211f, 155.791f, 177.504f, 200.686f, 225.641f, 252.645f, 281.906f, 313.850f, 348.302f, 385.732f, 426.307f, 470.146f, 517.555f, 568.899f, 624.177f, 684.021f, 748.734f, 818.078f, 892.577f, 973.087f, 1058.789f, 1150.915f, 1249.319f, 1354.279f, 1465.979f, 1584.790f, 1710.778f, 1844.863f, 1985.746f, 2134.643f, 2291.610f, 2456.943f, 2630.653f, 2813.959f, 3006.631f, 3207.229f, 3417.641f, 3637.318f, 3865.785f, 4104.997f, 4354.938f, 4615.786f, 4885.335f, 5166.555f, 5458.021f, 5762.584f, 6077.881f, 6406.834f, 6746.435f, 7097.958f, 7462.579f, 7839.165f, 8231.629f, 8635.640f, 9052.000f, 9484.268f, 9929.111f, 10389.350f, 10862.059f, 11352.185f, 11856.823f, 12380.371f, 12920.401f, 13476.971f, 14053.087f, 14646.190f, 15258.426f, 15890.617f, 16544.433f, 17218.024f, 17913.465f, 18631.374f, 19374.983f, 20136.700f, 20927.783f, 21746.796f, 22590.880f, 23465.734f, 24372.274f, 25314.351f, 26290.488f, 27300.899f, 28347.512f, 29436.133f, 30567.840f, 31746.818f, 32982.664f, 34276.329f, 35624.859f, 37042.588f, 38546.609f, 40139.742f, 41837.980f, 43679.429f, 45892.130f, 400000.000f}, "raw centrality signal"}; // for QA

        ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};

        // topological variable QA axes
        ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {20, 0.0f, 1.0f}, "DCA (cm)"};
        ConfigurableAxis axisDCAdau{"axisDCAdau", {20, 0.0f, 2.0f}, "DCA (cm)"};
        ConfigurableAxis axisPointingAngle{"axisPointingAngle", {20, 0.0f, 2.0f}, "pointing angle (rad)"};
        ConfigurableAxis axisV0Radius{"axisV0Radius", {20, 0.0f, 60.0f}, "V0 2D radius (cm)"};
        ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
        ConfigurableAxis axisTPCsignal{"axisTPCsignal", {200, 0.0f, 200.0f}, "TPC signal"};
        ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {200, -10.0f, 10.0f}, "N sigma TOF"};
        ConfigurableAxis axisTOFdeltaT{"axisTOFdeltaT", {200, -5000.0f, 5000.0f}, "TOF Delta T (ps)"};
        ConfigurableAxis axisPhi{"axisPhi", {18, 0.0f, constants::math::TwoPI}, "Azimuth angle (rad)"};
        ConfigurableAxis axisPhiMod{"axisPhiMod", {100, 0.0f, constants::math::TwoPI / 18}, "Azimuth angle wrt TPC sector (rad.)"};
        ConfigurableAxis axisEta{"axisEta", {10, -1.0f, 1.0f}, "#eta"};
        ConfigurableAxis axisITSchi2{"axisITSchi2", {100, 0.0f, 100.0f}, "#chi^{2} per ITS clusters"};
        ConfigurableAxis axisTPCchi2{"axisTPCchi2", {100, 0.0f, 100.0f}, "#chi^{2} per TPC clusters"};
        ConfigurableAxis axisTPCrowsOverFindable{"axisTPCrowsOverFindable", {120, 0.0f, 1.2f}, "Fraction of TPC crossed rows over findable clusters"};
        ConfigurableAxis axisTPCfoundOverFindable{"axisTPCfoundOverFindable", {120, 0.0f, 1.2f}, "Fraction of TPC found over findable clusters"};
        ConfigurableAxis axisTPCsharedClusters{"axisTPCsharedClusters", {101, -0.005f, 1.005f}, "Fraction of TPC shared clusters"};

        // UPC axes
        ConfigurableAxis axisSelGap{"axisSelGap", {4, -1.5, 2.5}, "Gap side"};

        // AP plot axes
        ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
        ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

        // Track quality axes
        ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
        ConfigurableAxis axisITSclus{"axisITSclus", {7, 0.0f, 7.0f}, "N ITS Clusters"};
        ConfigurableAxis axisITScluMap{"axisITScluMap", {128, -0.5f, 127.5f}, "ITS Cluster map"};
        ConfigurableAxis axisDetMap{"axisDetMap", {16, -0.5f, 15.5f}, "Detector use map"};
        ConfigurableAxis axisITScluMapCoarse{"axisITScluMapCoarse", {16, -3.5f, 12.5f}, "ITS Coarse cluster map"};
        ConfigurableAxis axisDetMapCoarse{"axisDetMapCoarse", {5, -0.5f, 4.5f}, "Detector Coarse user map"};

        // MC coll assoc QA axis
        ConfigurableAxis axisMonteCarloNch{"axisMonteCarloNch", {300, 0.0f, 3000.0f}, "N_{ch} MC"};
    } axisConfigurations;

    // These will not be used, but the variables were needed for some cut David did:
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


    // For manual sliceBy
    // Preslice<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;
    PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;
    PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCentsRun2, aod::StraEvSelsRun2, aod::StraCollLabels>> perMcCollisionRun2 = aod::v0data::straMCCollisionId;

    enum Selection : uint64_t { selCosPA = 0,
                                selRadius,
                                selRadiusMax,
                                selDCANegToPV,
                                selDCAPosToPV,
                                selDCAV0Dau,
                                selK0ShortRapidity,
                                selLambdaRapidity,
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
    uint64_t maskTrackProperties;

    uint64_t maskK0ShortSpecific;
    uint64_t maskLambdaSpecific;
    uint64_t maskAntiLambdaSpecific;

    uint64_t maskSelectionK0Short;
    uint64_t maskSelectionLambda;
    uint64_t maskSelectionAntiLambda;

    uint64_t secondaryMaskSelectionLambda;
    uint64_t secondaryMaskSelectionAntiLambda;


    void init(InitContext const&){
        // const AxisSpec axisCounter{1, 0, +1, ""};
        // // invMassAx = AxisSpec{nBinsInvMass, minInvMass, maxInvMass};

        // histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
        // // histos.add("LambdaInvMass1D", "Test LambdaInvMass 1D", kTH1F, {axisLambdaMass});
        histos.add("hMassLambda", "hMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
        // histos.add("ptQAHist", "ptQAHist", kTH1F, {axisPtQA});


        ///////////////////////////////////////////////////////////
            // Event Counters
        histos.add("hEventSelection", "hEventSelection", kTH1D, {{21, -0.5f, +20.5f}});
        if (isRun3) {
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
        } else {
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "sel7 cut");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kINT7");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kIsTriggerTVX");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "kNoIncompleteDAQ");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(7, "posZ cut");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(8, "kNoInconsistentVtx");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(9, "kNoPileupFromSPD");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(10, "kNoV0PFPileup");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(11, "kNoPileupInMultBins");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(12, "kNoPileupMV");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(13, "kNoPileupTPC");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(14, "kNoV0MOnVsOfPileup");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(15, "kNoSPDOnVsOfPileup");
        histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(16, "kNoSPDClsVsTklBG");
        if (doPPAnalysis) {
            histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(17, "INEL>0");
            histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(18, "INEL>1");
        }
        }

        histos.add("hEventCentrality", "hEventCentrality", kTH1D, {{101, 0.0f, 101.0f}});
        histos.add("hCentralityVsNch", "hCentralityVsNch", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisNch});
        if (doEventQA) {
        if (isRun3) {
            histos.add("hEventSelectionVsCentrality", "hEventSelectionVsCentrality", kTH2D, {{21, -0.5f, +20.5f}, {101, 0.0f, 101.0f}});
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(1, "All collisions");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(6, "posZ cut");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(7, "kIsVertexITSTPC");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(9, "kIsVertexTOFmatched");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(10, "kIsVertexTRDmatched");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(11, "kNoSameBunchPileup");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeStd");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeStrict");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(14, "kNoCollInTimeRangeNarrow");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(15, "kNoCollInRofStd");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(16, "kNoCollInRofStrict");
            if (doPPAnalysis) {
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(17, "INEL>0");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(18, "INEL>1");
            } else {
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(17, "Below min occup.");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(18, "Above max occup.");
            }
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(19, "Below min IR");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(20, "Above max IR");
            histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(21, "RCT flags");

            histos.add("hCentralityVsNGlobal", "hCentralityVsNGlobal", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisNch});
            histos.add("hEventCentVsMultFT0M", "hEventCentVsMultFT0M", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisMultFT0M});
            histos.add("hEventCentVsMultFT0C", "hEventCentVsMultFT0C", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisMultFT0C});
            histos.add("hEventCentVsMultNGlobal", "hEventCentVsMultNGlobal", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisNch});
            histos.add("hEventCentVsMultFV0A", "hEventCentVsMultFV0A", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisMultFV0A});
            histos.add("hEventMultFT0MvsMultNGlobal", "hEventMultFT0MvsMultNGlobal", kTH2D, {axisConfigurations.axisMultFT0M, axisConfigurations.axisNch});
            histos.add("hEventMultFT0CvsMultNGlobal", "hEventMultFT0CvsMultNGlobal", kTH2D, {axisConfigurations.axisMultFT0C, axisConfigurations.axisNch});
            histos.add("hEventMultFV0AvsMultNGlobal", "hEventMultFV0AvsMultNGlobal", kTH2D, {axisConfigurations.axisMultFV0A, axisConfigurations.axisNch});
            histos.add("hEventMultPVvsMultNGlobal", "hEventMultPVvsMultNGlobal", kTH2D, {axisConfigurations.axisNch, axisConfigurations.axisNch});
            histos.add("hEventMultFT0CvsMultFV0A", "hEventMultFT0CvsMultFV0A", kTH2D, {axisConfigurations.axisMultFT0C, axisConfigurations.axisMultFV0A});
        }
        }

        histos.add("hEventPVz", "hEventPVz", kTH1D, {{100, -20.0f, +20.0f}});
        histos.add("hCentralityVsPVz", "hCentralityVsPVz", kTH2D, {{101, 0.0f, 101.0f}, {100, -20.0f, +20.0f}});
        if (isRun3) {
        histos.add("hEventPVzMC", "hEventPVzMC", kTH1D, {{100, -20.0f, +20.0f}});
        histos.add("hCentralityVsPVzMC", "hCentralityVsPVzMC", kTH2D, {{101, 0.0f, 101.0f}, {100, -20.0f, +20.0f}});
        }

        histos.add("hEventOccupancy", "hEventOccupancy", kTH1D, {axisConfigurations.axisOccupancy});
        histos.add("hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisOccupancy});

        histos.add("hGapSide", "Gap side; Entries", kTH1D, {{5, -0.5, 4.5}});
        histos.add("hSelGapSide", "Selected gap side; Entries", kTH1D, {axisConfigurations.axisSelGap});
        histos.add("hEventCentralityVsSelGapSide", ";Centrality (%); Selected gap side", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisSelGap});

        histos.add("hInteractionRate", "hInteractionRate", kTH1D, {axisConfigurations.axisIRBinning});
        histos.add("hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisIRBinning});

        histos.add("hInteractionRateVsOccupancy", "hInteractionRateVsOccupancy", kTH2D, {axisConfigurations.axisIRBinning, axisConfigurations.axisOccupancy});

        // for QA and test purposes
        auto hRawCentrality = histos.add<TH1>("hRawCentrality", "hRawCentrality", kTH1D, {axisConfigurations.axisRawCentrality});

        for (int ii = 1; ii < 101; ii++) {
        float value = 100.5f - static_cast<float>(ii);
        hRawCentrality->SetBinContent(ii, value);
        }

        auto hSelectionV0s = histos.add<TH1>("GeneralQA/hSelectionV0s", "hSelectionV0s", kTH1D, {{static_cast<int>(selPhysPrimAntiLambda) + 3, -0.5f, static_cast<double>(selPhysPrimAntiLambda) + 2.5f}});
        hSelectionV0s->GetXaxis()->SetBinLabel(1, "All");
        hSelectionV0s->GetXaxis()->SetBinLabel(selCosPA + 2, "cosPA");
        hSelectionV0s->GetXaxis()->SetBinLabel(selRadius + 2, "Radius min.");
        hSelectionV0s->GetXaxis()->SetBinLabel(selRadiusMax + 2, "Radius max.");
        hSelectionV0s->GetXaxis()->SetBinLabel(selDCANegToPV + 2, "DCA neg. to PV");
        hSelectionV0s->GetXaxis()->SetBinLabel(selDCAPosToPV + 2, "DCA pos. to PV");
        hSelectionV0s->GetXaxis()->SetBinLabel(selDCAV0Dau + 2, "DCA V0 dau.");
        hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortRapidity + 2, "K^{0}_{S} rapidity");
        hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaRapidity + 2, "#Lambda rapidity");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTPCPIDPositivePion + 2, "TPC PID #pi^{+}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTPCPIDNegativePion + 2, "TPC PID #pi^{-}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTPCPIDPositiveProton + 2, "TPC PID p");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTPCPIDNegativeProton + 2, "TPC PID #bar{p}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFDeltaTPositiveProtonLambda + 2, "TOF #Delta t p from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFDeltaTPositivePionLambda + 2, "TOF #Delta t #pi^{+} from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFDeltaTPositivePionK0Short + 2, "TOF #Delta t #pi^{+} from K^{0}_{S}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFDeltaTNegativeProtonLambda + 2, "TOF #Delta t #bar{p} from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFDeltaTNegativePionLambda + 2, "TOF #Delta t #pi^{-} from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFDeltaTNegativePionK0Short + 2, "TOF #Delta t #pi^{-} from K^{0}_{S}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFNSigmaPositiveProtonLambda + 2, "TOF PID p from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFNSigmaPositivePionLambda + 2, "TOF PID #pi^{+} from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFNSigmaPositivePionK0Short + 2, "TOF PID #pi^{+} from K^{0}_{S}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFNSigmaNegativeProtonLambda + 2, "TOF PID #bar{p} from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFNSigmaNegativePionLambda + 2, "TOF PID #pi^{-} from #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selTOFNSigmaNegativePionK0Short + 2, "TOF PID #pi^{-} from K^{0}_{S}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortCTau + 2, "K^{0}_{S} lifetime");
        hSelectionV0s->GetXaxis()->SetBinLabel(selLambdaCTau + 2, "#Lambda lifetime");
        hSelectionV0s->GetXaxis()->SetBinLabel(selK0ShortArmenteros + 2, "Arm. pod. cut");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPosGoodTPCTrack + 2, "Pos. good TPC track");
        hSelectionV0s->GetXaxis()->SetBinLabel(selNegGoodTPCTrack + 2, "Neg. good TPC track");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPosGoodITSTrack + 2, "Pos. good ITS track");
        hSelectionV0s->GetXaxis()->SetBinLabel(selNegGoodITSTrack + 2, "Neg. good ITS track");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPosItsOnly + 2, "Pos. ITS-only");
        hSelectionV0s->GetXaxis()->SetBinLabel(selNegItsOnly + 2, "Neg. ITS-only");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPosNotTPCOnly + 2, "Pos. not TPC-only");
        hSelectionV0s->GetXaxis()->SetBinLabel(selNegNotTPCOnly + 2, "Neg. not TPC-only");
        hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderK0Short + 2, "True K^{0}_{S}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderLambda + 2, "True #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selConsiderAntiLambda + 2, "True #bar{#Lambda}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimK0Short + 2, "Phys. prim. K^{0}_{S}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimLambda + 2, "Phys. prim. #Lambda");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimAntiLambda + 2, "Phys. prim. #bar{#Lambda}");
        hSelectionV0s->GetXaxis()->SetBinLabel(selPhysPrimAntiLambda + 3, "Cand. selected");
        ///////////////////////////////////////////////////////////

        // From the analyseLambda flag:
        histos.add("h2dNbrOfLambdaVsCentrality", "h2dNbrOfLambdaVsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
        histos.add("h3dMassLambda", "h3dMassLambda", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        // Non-UPC info
        histos.add("h3dMassLambdaHadronic", "h3dMassLambdaHadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        // Not doing ultra-peripheral!
        // // UPC info
        // histos.add("h3dMassLambdaSGA", "h3dMassLambdaSGA", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        // histos.add("h3dMassLambdaSGC", "h3dMassLambdaSGC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        // histos.add("h3dMassLambdaDG", "h3dMassLambdaDG", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        
        // // For the doCompleteTopoQA analysis:
        // histos.add("Lambda/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAtoPV});
        // histos.add("Lambda/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAtoPV});
        // histos.add("Lambda/h4dDCADaughters", "h4dDCADaughters", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAdau});
        // histos.add("Lambda/h4dPointingAngle", "h4dPointingAngle", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPointingAngle});
        // histos.add("Lambda/h4dV0Radius", "h4dV0Radius", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisV0Radius});
        
        // For the doPlainTopoQA analysis:
            // All candidates received
        histos.add("hPosDCAToPV", "hPosDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("hNegDCAToPV", "hNegDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("hPointingAngle", "hPointingAngle", kTH1D, {axisConfigurations.axisPointingAngle});
        histos.add("hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("h2dPositivePtVsPhi", "h2dPositivePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
        histos.add("h2dNegativePtVsPhi", "h2dNegativePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});

            // Specifically for Lambda:
        histos.add("Lambda/hPosDCAToPV", "hPosDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("Lambda/hNegDCAToPV", "hNegDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("Lambda/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("Lambda/hPointingAngle", "hPointingAngle", kTH1D, {axisConfigurations.axisPointingAngle});
        histos.add("Lambda/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("Lambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("Lambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("Lambda/h2dPositivePtVsPhi", "h2dPositivePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
        histos.add("Lambda/h2dNegativePtVsPhi", "h2dNegativePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});

            // Check if doing the right thing in AP space please
        histos.add("GeneralQA/h2dArmenterosAll", "h2dArmenterosAll", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
        histos.add("GeneralQA/h2dArmenterosSelected", "h2dArmenterosSelected", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    }

    template <typename TCollision>
    void initCCDB(TCollision collision)
    {
        if (mRunNumber == collision.runNumber()) {
        return;
        }

        mRunNumber = collision.runNumber();

        // machine learning initialization if requested
        if (mlConfigurations.calculateLambdaScores) {
        int64_t timeStampML = collision.timestamp();
        if (mlConfigurations.timestampCCDB.value != -1)
            timeStampML = mlConfigurations.timestampCCDB.value;
        loadMachines(timeStampML);
        }
        // Fetching magnetic field if requested
        if (v0Selections.rejectTPCsectorBoundary) {
        // In case override, don't proceed, please - no CCDB access required
        if (ccdbConfigurations.useCustomMagField) {
            magField = ccdbConfigurations.customMagField;
        } else {
            grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, mRunNumber);
            if (!grpmag) {
            LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField and " << ccdbConfigurations.grpPath << " of object GRPObject for run " << mRunNumber;
            }
            // Fetch magnetic field from ccdb for current collision
            magField = std::lround(5.f * grpmag->getL3Current() / 30000.f);
            LOG(info) << "Retrieved GRP for run " << mRunNumber << " with magnetic field of " << magField << " kZG";
        }
        }
    }

    double computePhiMod(double phi, int sign)
    // Compute phi wrt to a TPC sector
    // Calculation taken from CF: https://github.com/AliceO2Group/O2Physics/blob/376392cb87349886a300c75fa2492b50b7f46725/PWGCF/Flow/Tasks/flowAnalysisGF.cxx#L470
    {
        if (magField < 0) // for negative polarity field
        phi = o2::constants::math::TwoPI - phi;
        if (sign < 0) // for negative charge
        phi = o2::constants::math::TwoPI - phi;
        if (phi < 0)
        LOGF(warning, "phi < 0: %g", phi);

        phi += o2::constants::math::PI / 18.0; // to center gap in the middle
        return fmod(phi, o2::constants::math::PI / 9.0);
    }

    bool isTrackFarFromTPCBoundary(double trackPt, double trackPhi, int sign)
    // check whether the track passes close to a TPC sector boundary
    {
        double phiModn = computePhiMod(trackPhi, sign);
        if (phiModn > fPhiCutHigh->Eval(trackPt))
        return true; // keep track
        if (phiModn < fPhiCutLow->Eval(trackPt))
        return true; // keep track
        return false;  // reject track
    }

    template <typename TV0, typename TCollision>
    // uint64_t computeReconstructionBitmap(TV0 v0, TCollision collision, float rapidityLambda, float rapidityK0Short, float /*pT*/)
    uint64_t computeReconstructionBitmap(TV0 v0, TCollision collision, float rapidityLambda, float /*pT*/)
    // precalculate this information so that a check is one mask operation, not many
    {
        uint64_t bitMap = 0;
        // Base topological variables
        if (v0.v0radius() > v0Selections.v0radius)
            BITSET(bitMap, selRadius);
        if (v0.v0radius() < v0Selections.v0radiusMax)
            BITSET(bitMap, selRadiusMax);
        if (std::abs(v0.dcapostopv()) > v0Selections.dcapostopv)
            BITSET(bitMap, selDCAPosToPV);
        if (std::abs(v0.dcanegtopv()) > v0Selections.dcanegtopv)
            BITSET(bitMap, selDCANegToPV);
        if (v0.v0cosPA() > v0Selections.v0cospa)
            BITSET(bitMap, selCosPA);
        if (v0.dcaV0daughters() < v0Selections.dcav0dau)
            BITSET(bitMap, selDCAV0Dau);

        // rapidity
        if (std::abs(rapidityLambda) < v0Selections.rapidityCut)
            BITSET(bitMap, selLambdaRapidity);

        auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
        auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

        // ITS quality flags
        bool posIsFromAfterburner = posTrackExtra.hasITSAfterburner();
        bool negIsFromAfterburner = negTrackExtra.hasITSAfterburner();

        // check minimum number of ITS clusters + maximum ITS chi2 per clusters + reject or select ITS afterburner tracks if requested
        if (posTrackExtra.itsNCls() >= v0Selections.minITSclusters &&             // check minium ITS clusters
            posTrackExtra.itsChi2NCl() < v0Selections.maxITSchi2PerNcls &&        // check maximum ITS chi2 per clusters
            (!v0Selections.rejectPosITSafterburner || !posIsFromAfterburner) &&   // reject afterburner track or not
            (!v0Selections.requirePosITSafterburnerOnly || posIsFromAfterburner)) // keep afterburner track or not
            BITSET(bitMap, selPosGoodITSTrack);
        if (negTrackExtra.itsNCls() >= v0Selections.minITSclusters &&             // check minium ITS clusters
            negTrackExtra.itsChi2NCl() < v0Selections.maxITSchi2PerNcls &&        // check maximum ITS chi2 per clusters
            (!v0Selections.rejectNegITSafterburner || !negIsFromAfterburner) &&   // reject afterburner track or not
            (!v0Selections.requireNegITSafterburnerOnly || negIsFromAfterburner)) // select only afterburner track or not
            BITSET(bitMap, selNegGoodITSTrack);

        // TPC quality flags
        if (posTrackExtra.tpcCrossedRows() >= v0Selections.minTPCrows &&                                                // check minimum TPC crossed rows
            posTrackExtra.tpcChi2NCl() < v0Selections.maxTPCchi2PerNcls &&                                              // check maximum TPC chi2 per clusters
            posTrackExtra.tpcCrossedRowsOverFindableCls() >= v0Selections.minTPCrowsOverFindableClusters &&             // check minimum fraction of TPC rows over findable
            posTrackExtra.tpcFoundOverFindableCls() >= v0Selections.minTPCfoundOverFindableClusters &&                  // check minimum fraction of found over findable TPC clusters
            posTrackExtra.tpcFractionSharedCls() < v0Selections.maxFractionTPCSharedClusters &&                         // check the maximum fraction of allowed shared TPC clusters
            (!v0Selections.rejectTPCsectorBoundary || isTrackFarFromTPCBoundary(v0.positivept(), v0.positivephi(), 1))) // reject track far from TPC sector boundary or not
            BITSET(bitMap, selPosGoodTPCTrack);
        if (negTrackExtra.tpcCrossedRows() >= v0Selections.minTPCrows &&                                                 // check minimum TPC crossed rows
            negTrackExtra.tpcChi2NCl() < v0Selections.maxTPCchi2PerNcls &&                                               // check maximum TPC chi2 per clusters
            negTrackExtra.tpcCrossedRowsOverFindableCls() >= v0Selections.minTPCrowsOverFindableClusters &&              // check minimum fraction of TPC rows over findable
            negTrackExtra.tpcFoundOverFindableCls() >= v0Selections.minTPCfoundOverFindableClusters &&                   // check minimum fraction of found over findable TPC clusters
            negTrackExtra.tpcFractionSharedCls() < v0Selections.maxFractionTPCSharedClusters &&                          // check the maximum fraction of allowed shared TPC clusters
            (!v0Selections.rejectTPCsectorBoundary || isTrackFarFromTPCBoundary(v0.negativept(), v0.negativephi(), -1))) // reject track far from TPC sector boundary or not
            BITSET(bitMap, selNegGoodTPCTrack);

        // TPC PID
        if (std::fabs(posTrackExtra.tpcNSigmaPi()) < v0Selections.tpcPidNsigmaCut)
            BITSET(bitMap, selTPCPIDPositivePion);
        if (std::fabs(posTrackExtra.tpcNSigmaPr()) < v0Selections.tpcPidNsigmaCut)
            BITSET(bitMap, selTPCPIDPositiveProton);
        if (std::fabs(negTrackExtra.tpcNSigmaPi()) < v0Selections.tpcPidNsigmaCut)
            BITSET(bitMap, selTPCPIDNegativePion);
        if (std::fabs(negTrackExtra.tpcNSigmaPr()) < v0Selections.tpcPidNsigmaCut)
            BITSET(bitMap, selTPCPIDNegativeProton);

        // TOF PID in DeltaT
        // Positive track
        if (!posTrackExtra.hasTOF() || std::fabs(v0.posTOFDeltaTLaPr()) < v0Selections.maxDeltaTimeProton)
            BITSET(bitMap, selTOFDeltaTPositiveProtonLambda);
        if (!posTrackExtra.hasTOF() || std::fabs(v0.posTOFDeltaTLaPi()) < v0Selections.maxDeltaTimePion)
            BITSET(bitMap, selTOFDeltaTPositivePionLambda);
        // Negative track
        if (!negTrackExtra.hasTOF() || std::fabs(v0.negTOFDeltaTLaPr()) < v0Selections.maxDeltaTimeProton)
            BITSET(bitMap, selTOFDeltaTNegativeProtonLambda);
        if (!negTrackExtra.hasTOF() || std::fabs(v0.negTOFDeltaTLaPi()) < v0Selections.maxDeltaTimePion)
            BITSET(bitMap, selTOFDeltaTNegativePionLambda);

        // TOF PID in NSigma
        // Positive track
        if (!posTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaLaPr()) < v0Selections.tofPidNsigmaCutLaPr)
            BITSET(bitMap, selTOFNSigmaPositiveProtonLambda);
        if (!posTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaALaPi()) < v0Selections.tofPidNsigmaCutLaPi)
            BITSET(bitMap, selTOFNSigmaPositivePionLambda);
        // Negative track
        if (!negTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaALaPr()) < v0Selections.tofPidNsigmaCutLaPr)
            BITSET(bitMap, selTOFNSigmaNegativeProtonLambda);
        if (!negTrackExtra.hasTOF() || std::fabs(v0.tofNSigmaLaPi()) < v0Selections.tofPidNsigmaCutLaPi)
            BITSET(bitMap, selTOFNSigmaNegativePionLambda);

        // ITS only tag
        if (posTrackExtra.tpcCrossedRows() < 1)
            BITSET(bitMap, selPosItsOnly);
        if (negTrackExtra.tpcCrossedRows() < 1)
            BITSET(bitMap, selNegItsOnly);

        // TPC only tag
        if (posTrackExtra.detectorMap() != o2::aod::track::TPC)
            BITSET(bitMap, selPosNotTPCOnly);
        if (negTrackExtra.detectorMap() != o2::aod::track::TPC)
            BITSET(bitMap, selNegNotTPCOnly);

        // proper lifetime
        if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < v0Selections.lifetimecut->get("lifetimecutLambda"))
            BITSET(bitMap, selLambdaCTau);

        return bitMap;
    }

    bool verifyMask(uint64_t bitmap, uint64_t mask)
    {
        return (bitmap & mask) == mask;
    }

    int computeITSclusBitmap(uint8_t itsClusMap, bool fromAfterburner)
    // Focus on the 12 dominant ITS cluster configurations
    {
        int bitMap = 0;

        if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS : x  x  x  x  x  x  x
        bitMap = 12;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :    x  x  x  x  x  x
        bitMap = 11;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :       x  x  x  x  x
        bitMap = 10;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :          x  x  x  x
        bitMap = 9;
        if (fromAfterburner)
            bitMap = -3;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 4) | (uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :             x  x  x
        bitMap = 8;
        if (fromAfterburner)
            bitMap = -2;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 5) | (uint8_t(1) << 6)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :                x  x
        bitMap = 7;
        if (fromAfterburner)
            bitMap = -1;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS : x  x  x  x  x  x
        bitMap = 6;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :    x  x  x  x  x
        bitMap = 5;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4) | (uint8_t(1) << 5)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :       x  x  x  x
        bitMap = 4;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS : x  x  x  x  x
        bitMap = 3;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3) | (uint8_t(1) << 4)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS :    x  x  x  x
        bitMap = 2;
        } else if (verifyMask(itsClusMap, ((uint8_t(1) << 0) | (uint8_t(1) << 1) | (uint8_t(1) << 2) | (uint8_t(1) << 3)))) {
        // ITS :    IB         OB
        // ITS : L0 L1 L2 L3 L4 L5 L6
        // ITS : x  x  x  x
        bitMap = 1;
        } else {
        // ITS : other configurations
        bitMap = 0;
        }

        return bitMap;
    }

    uint computeDetBitmap(uint8_t detMap)
    // Focus on the 4 dominant track configurations :
    //  Others
    //  ITS-TPC
    //  ITS-TPC-TRD
    //  ITS-TPC-TOF
    //  ITS-TPC-TRD-TOF
    {
        uint bitMap = 0;

        if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD | o2::aod::track::TOF))) {
        // ITS-TPC-TRD-TOF
        bitMap = 4;
        } else if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TOF))) {
        // ITS-TPC-TOF
        bitMap = 3;
        } else if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD))) {
        // ITS-TPC-TRD
        bitMap = 2;
        } else if (verifyMask(detMap, (o2::aod::track::ITS | o2::aod::track::TPC))) {
        // ITS-TPC
        bitMap = 1;
        }

        return bitMap;
    }

    // function to load models for ML-based classifiers
    void loadMachines(int64_t timeStampML)
    {
        if (mlConfigurations.loadCustomModelsFromCCDB) {
        ccdbApi.init(ccdbConfigurations.ccdbUrl);
        LOG(info) << "Fetching models for timestamp: " << timeStampML;

        if (mlConfigurations.calculateLambdaScores) {
            bool retrieveSuccessLambda = ccdbApi.retrieveBlob(mlConfigurations.customModelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathLambda.value);
            if (retrieveSuccessLambda) {
            mlCustomModelLambda.initModel(mlConfigurations.localModelPathLambda.value, mlConfigurations.enableOptimizations.value);
            } else {
            LOG(fatal) << "Error encountered while fetching/loading the Lambda model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
            }
        }
        }
        if (mlConfigurations.calculateLambdaScores)
            mlCustomModelLambda.initModel(mlConfigurations.localModelPathLambda.value, mlConfigurations.enableOptimizations.value);
        LOG(info) << "ML Models loaded.";
    }

    template <typename TV0>
    // void analyseCandidate(TV0 v0, float pt, float centrality, uint64_t selMap, uint8_t gapSide, int& nK0Shorts, int& nLambdas, int& nAntiLambdas)
    void analyseCandidate(TV0 v0, float pt, float centrality, uint64_t selMap, uint8_t gapSide, int& nLambdas)
    // precalculate this information so that a check is one mask operation, not many
    {
        bool passLambdaSelections = false;

        // machine learning is on, go for calculation of thresholds
        // FIXME THIS NEEDS ADJUSTING
        std::vector<float> inputFeatures{pt, 0.0f, 0.0f, v0.v0radius(), v0.v0cosPA(), v0.dcaV0daughters(), v0.dcapostopv(), v0.dcanegtopv()};

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

        auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
        auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

        bool posIsFromAfterburner = posTrackExtra.itsChi2PerNcl() < 0;
        bool negIsFromAfterburner = negTrackExtra.itsChi2PerNcl() < 0;

        uint posDetMap = computeDetBitmap(posTrackExtra.detectorMap());
        int posITSclusMap = computeITSclusBitmap(posTrackExtra.itsClusterMap(), posIsFromAfterburner);
        uint negDetMap = computeDetBitmap(negTrackExtra.detectorMap());
        int negITSclusMap = computeITSclusBitmap(negTrackExtra.itsClusterMap(), negIsFromAfterburner);

        // __________________________________________
        // fill with no selection if plain QA requested
        if (doPlainTopoQA) {
        histos.fill(HIST("hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("hPointingAngle"), std::acos(v0.v0cosPA()));
        histos.fill(HIST("hV0Radius"), v0.v0radius());
        histos.fill(HIST("h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
        histos.fill(HIST("h2dPositivePtVsPhi"), v0.positivept(), computePhiMod(v0.positivephi(), 1));
        histos.fill(HIST("h2dNegativePtVsPhi"), v0.negativept(), computePhiMod(v0.negativephi(), -1));
        }

        // Fill first bin: all candidates
        histos.fill(HIST("GeneralQA/hSelectionV0s"), 0);
        // Loop over all bits in the enum and fill if passed
        for (uint64_t i = 0; i <= selPhysPrimAntiLambda; i++) {
        if (BITCHECK(selMap, i)) {
            histos.fill(HIST("GeneralQA/hSelectionV0s"), i + 1); // +1 because bin 0 = "All"
        }
        }

        // __________________________________________
        // main analysis
        if (passLambdaSelections && analyseLambda) {
        histos.fill(HIST("GeneralQA/hSelectionV0s"), selPhysPrimAntiLambda + 2); //
        histos.fill(HIST("h3dMassLambda"), centrality, pt, v0.mLambda());
        if (gapSide == 0)
            histos.fill(HIST("h3dMassLambdaSGA"), centrality, pt, v0.mLambda());
        else if (gapSide == 1)
            histos.fill(HIST("h3dMassLambdaSGC"), centrality, pt, v0.mLambda());
        else if (gapSide == 2)
            histos.fill(HIST("h3dMassLambdaDG"), centrality, pt, v0.mLambda());
        else
            histos.fill(HIST("h3dMassLambdaHadronic"), centrality, pt, v0.mLambda());
        histos.fill(HIST("hMassLambda"), v0.mLambda());
        if (doPlainTopoQA) {
            histos.fill(HIST("Lambda/hPosDCAToPV"), v0.dcapostopv());
            histos.fill(HIST("Lambda/hNegDCAToPV"), v0.dcanegtopv());
            histos.fill(HIST("Lambda/hDCADaughters"), v0.dcaV0daughters());
            histos.fill(HIST("Lambda/hPointingAngle"), std::acos(v0.v0cosPA()));
            histos.fill(HIST("Lambda/hV0Radius"), v0.v0radius());
            histos.fill(HIST("Lambda/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
            histos.fill(HIST("Lambda/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
            histos.fill(HIST("Lambda/h2dPositivePtVsPhi"), v0.positivept(), computePhiMod(v0.positivephi(), 1));
            histos.fill(HIST("Lambda/h2dNegativePtVsPhi"), v0.negativept(), computePhiMod(v0.negativephi(), -1));
        }
        if (doDetectPropQA == 1) {
            histos.fill(HIST("Lambda/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt);
            histos.fill(HIST("Lambda/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt);
            histos.fill(HIST("Lambda/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt);
        }
        if (doDetectPropQA == 2) {
            histos.fill(HIST("Lambda/h7dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pt, v0.mLambda());
            histos.fill(HIST("Lambda/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pt, v0.mLambda());
            histos.fill(HIST("Lambda/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pt, v0.mLambda());
        }
        if (doDetectPropQA == 3) {
            histos.fill(HIST("Lambda/h3dITSchi2"), centrality, pt, std::max(posTrackExtra.itsChi2NCl(), negTrackExtra.itsChi2NCl()));
            histos.fill(HIST("Lambda/h3dTPCchi2"), centrality, pt, std::max(posTrackExtra.tpcChi2NCl(), negTrackExtra.tpcChi2NCl()));
            histos.fill(HIST("Lambda/h3dTPCFoundOverFindable"), centrality, pt, std::min(posTrackExtra.tpcFoundOverFindableCls(), negTrackExtra.tpcFoundOverFindableCls()));
            histos.fill(HIST("Lambda/h3dTPCrowsOverFindable"), centrality, pt, std::min(posTrackExtra.tpcCrossedRowsOverFindableCls(), negTrackExtra.tpcCrossedRowsOverFindableCls()));
            histos.fill(HIST("Lambda/h3dTPCsharedCls"), centrality, pt, std::max(posTrackExtra.tpcFractionSharedCls(), negTrackExtra.tpcFractionSharedCls()));
            histos.fill(HIST("Lambda/h3dPositiveITSchi2"), centrality, pt, posTrackExtra.itsChi2NCl());
            histos.fill(HIST("Lambda/h3dNegativeITSchi2"), centrality, pt, negTrackExtra.itsChi2NCl());
            histos.fill(HIST("Lambda/h3dPositiveTPCchi2"), centrality, pt, posTrackExtra.tpcChi2NCl());
            histos.fill(HIST("Lambda/h3dNegativeTPCchi2"), centrality, pt, negTrackExtra.tpcChi2NCl());
            histos.fill(HIST("Lambda/h3dPositiveITSclusters"), centrality, pt, posTrackExtra.itsNCls());
            histos.fill(HIST("Lambda/h3dNegativeITSclusters"), centrality, pt, negTrackExtra.itsNCls());
            histos.fill(HIST("Lambda/h3dPositiveTPCcrossedRows"), centrality, pt, posTrackExtra.tpcCrossedRows());
            histos.fill(HIST("Lambda/h3dNegativeTPCcrossedRows"), centrality, pt, negTrackExtra.tpcCrossedRows());
        }
        if (doTPCQA) {
            histos.fill(HIST("Lambda/h3dPosNsigmaTPC"), centrality, pt, posTrackExtra.tpcNSigmaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTPC"), centrality, pt, negTrackExtra.tpcNSigmaPi());
            histos.fill(HIST("Lambda/h3dPosTPCsignal"), centrality, pt, posTrackExtra.tpcSignal());
            histos.fill(HIST("Lambda/h3dNegTPCsignal"), centrality, pt, negTrackExtra.tpcSignal());
            histos.fill(HIST("Lambda/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.pfracpos() * v0.p(), posTrackExtra.tpcNSigmaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.pfracneg() * v0.p(), negTrackExtra.tpcNSigmaPi());
            histos.fill(HIST("Lambda/h3dPosTPCsignalVsTrackPtot"), centrality, v0.pfracpos() * v0.p(), posTrackExtra.tpcSignal());
            histos.fill(HIST("Lambda/h3dNegTPCsignalVsTrackPtot"), centrality, v0.pfracneg() * v0.p(), negTrackExtra.tpcSignal());
            histos.fill(HIST("Lambda/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPi());
            histos.fill(HIST("Lambda/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
            histos.fill(HIST("Lambda/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
        }
        if (doTOFQA) {
            histos.fill(HIST("Lambda/h3dPosNsigmaTOF"), centrality, pt, v0.tofNSigmaLaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTOF"), centrality, pt, v0.tofNSigmaLaPi());
            histos.fill(HIST("Lambda/h3dPosTOFdeltaT"), centrality, pt, v0.posTOFDeltaTLaPr());
            histos.fill(HIST("Lambda/h3dNegTOFdeltaT"), centrality, pt, v0.negTOFDeltaTLaPi());
            histos.fill(HIST("Lambda/h3dPosNsigmaTOFvsTrackPtot"), centrality, v0.pfracpos() * v0.p(), v0.tofNSigmaLaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTOFvsTrackPtot"), centrality, v0.pfracneg() * v0.p(), v0.tofNSigmaLaPi());
            histos.fill(HIST("Lambda/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.pfracpos() * v0.p(), v0.posTOFDeltaTLaPr());
            histos.fill(HIST("Lambda/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.pfracneg() * v0.p(), v0.negTOFDeltaTLaPi());
            histos.fill(HIST("Lambda/h3dPosNsigmaTOFvsTrackPt"), centrality, v0.positivept(), v0.tofNSigmaLaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTOFvsTrackPt"), centrality, v0.negativept(), v0.tofNSigmaLaPi());
            histos.fill(HIST("Lambda/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTLaPr());
            histos.fill(HIST("Lambda/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTLaPi());
        }
        if (doEtaPhiQA) {
            histos.fill(HIST("Lambda/h5dV0PhiVsEta"), centrality, pt, v0.mLambda(), v0.phi(), v0.eta());
            histos.fill(HIST("Lambda/h5dPosPhiVsEta"), centrality, v0.positivept(), v0.mLambda(), v0.positivephi(), v0.positiveeta());
            histos.fill(HIST("Lambda/h5dNegPhiVsEta"), centrality, v0.negativept(), v0.mLambda(), v0.negativephi(), v0.negativeeta());
        }
        nLambdas++;
        }

        // // __________________________________________
        // // do systematics / qa plots
        // if (doCompleteTopoQA) {
        // if (analyseLambda) {
        //     if (verifyMask(selMap, maskTopoNoV0Radius | maskLambdaSpecific))
        //     histos.fill(HIST("Lambda/h4dV0Radius"), centrality, pt, v0.mLambda(), v0.v0radius());
        //     if (verifyMask(selMap, maskTopoNoDCAPosToPV | maskLambdaSpecific))
        //     histos.fill(HIST("Lambda/h4dPosDCAToPV"), centrality, pt, v0.mLambda(), std::abs(v0.dcapostopv()));
        //     if (verifyMask(selMap, maskTopoNoDCANegToPV | maskLambdaSpecific))
        //     histos.fill(HIST("Lambda/h4dNegDCAToPV"), centrality, pt, v0.mLambda(), std::abs(v0.dcanegtopv()));
        //     if (verifyMask(selMap, maskTopoNoCosPA | maskLambdaSpecific))
        //     histos.fill(HIST("Lambda/h4dPointingAngle"), centrality, pt, v0.mLambda(), std::acos(v0.v0cosPA()));
        //     if (verifyMask(selMap, maskTopoNoDCAV0Dau | maskLambdaSpecific))
        //     histos.fill(HIST("Lambda/h4dDCADaughters"), centrality, pt, v0.mLambda(), v0.dcaV0daughters());
        // } // end systematics / qa
        // }
    }

    // ______________________________________________________
    // Return slicing output
    template <typename TCollision>
    auto getCentralityRun3(TCollision const& collision)
    {
        if (centralityEstimator == kCentFT0C)
        return collision.centFT0C();
        else if (centralityEstimator == kCentFT0M)
        return collision.centFT0M();
        else if (centralityEstimator == kCentFT0CVariant1)
        return collision.centFT0CVariant1();
        else if (centralityEstimator == kCentMFT)
        return collision.centMFT();
        else if (centralityEstimator == kCentNGlobal)
        return collision.centNGlobal();
        else if (centralityEstimator == kCentFV0A)
        return collision.centFV0A();

        return -1.f;
    }


    // ______________________________________________________
    // Reconstructed data processing
    // Fill reconstructed event information
    // Return centrality, occupancy, interaction rate, gap side and selGapside via reference-passing in arguments
    template <typename TCollision>
    void fillReconstructedEventProperties(TCollision const& collision, float& centrality, float& collisionOccupancy, double& interactionRate, int& gapSide, int& selGapSide)
    {
        if constexpr (requires { collision.centFT0C(); }) { // check if we are in Run 3
        centrality = getCentralityRun3(collision);
        collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
        // Fetch interaction rate only if required (in order to limit ccdb calls)
        interactionRate = !irSource.value.empty() ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource) * 1.e-3 : -1;

        if (qaCentrality) {
            auto hRawCentrality = histos.get<TH1>(HIST("hRawCentrality"));
            centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(doPPAnalysis ? collision.multFT0A() + collision.multFT0C() : collision.multFT0C()));
        }

        // gap side
        gapSide = collision.gapSide();
        // -1 --> Hadronic
        // 0 --> Single Gap - A side
        // 1 --> Single Gap - C side
        // 2 --> Double Gap - both A & C sides
        selGapSide = sgSelector.trueGap(collision, upcCuts.fv0Cut, upcCuts.ft0Acut, upcCuts.ft0Ccut, upcCuts.zdcCut);
        } else { // no, we are in Run 2
        centrality = eventSelections.useSPDTrackletsCent ? collision.centRun2SPDTracklets() : collision.centRun2V0M();
        }

        histos.fill(HIST("hGapSide"), gapSide);
        histos.fill(HIST("hSelGapSide"), selGapSide);
        histos.fill(HIST("hEventCentralityVsSelGapSide"), centrality, selGapSide <= 2 ? selGapSide : -1);

        histos.fill(HIST("hEventCentrality"), centrality);

        histos.fill(HIST("hCentralityVsNch"), centrality, collision.multNTracksPVeta1());
        if (doEventQA) {
        if constexpr (requires { collision.centFT0C(); }) { // check if we are in Run 3
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

    template <typename TCollision>
    bool isEventAccepted(TCollision collision, bool fillHists)
    // check whether the collision passes our collision selections
    {
        float centrality = -1.0f;
        if (fillHists) {
        histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
        if (doEventQA) {
            if constexpr (requires { collision.centFT0C(); }) { // check if we are in Run 3
            centrality = getCentralityRun3(collision);
            }
            histos.fill(HIST("hEventSelectionVsCentrality"), 0. /* all collisions */, centrality);
        }
        }

        if constexpr (requires { collision.centFT0C(); }) { // check if we are in Run 3
        if (eventSelections.requireSel8 && !collision.sel8()) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 1 /* sel8 collisions */, centrality);
            }
        }

        if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */, centrality);
            }
        }

        if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 3 /* Not at ITS ROF border */, centrality);
            }
        }

        if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 4 /* Not at TF border */, centrality);
            }
        }

        if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 5 /* vertex-Z selected */, centrality);
            }
        }

        if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 6 /* Contains at least one ITS-TPC track */, centrality);
            }
        }

        if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 7 /* PV position consistency check */, centrality);
            }
        }

        if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 8 /* PV with at least one contributor matched with TOF */, centrality);
            }
        }

        if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 9 /* PV with at least one contributor matched with TRD */, centrality);
            }
        }

        if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 10 /* Not at same bunch pile-up */, centrality);
            }
        }

        if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/, centrality);
            }
        }

        if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 12 /* No other collision within +/- 10 microseconds */, centrality);
            }
        }

        if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 13 /* No other collision within +/- 2 microseconds */, centrality);
            }
        }

        if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */, centrality);
            }
        }

        if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 15 /* No other collision within the same ITS ROF */, centrality);
            }
        }

        if (doPPAnalysis) { // we are in pp
            if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) {
            return false;
            }
            if (fillHists) {
            histos.fill(HIST("hEventSelection"), 16 /* INEL > 0 */);
            if (doEventQA) {
                histos.fill(HIST("hEventSelectionVsCentrality"), 16 /* INEL > 0 */, centrality);
            }
            }

            if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) {
            return false;
            }
            if (fillHists) {
            histos.fill(HIST("hEventSelection"), 17 /* INEL > 1 */);
            if (doEventQA) {
                histos.fill(HIST("hEventSelectionVsCentrality"), 17 /* INEL > 1 */, centrality);
            }
            }

        } else { // we are in Pb-Pb
            float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
            if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
            return false;
            }
            if (fillHists) {
            histos.fill(HIST("hEventSelection"), 16 /* Below min occupancy */);
            if (doEventQA) {
                histos.fill(HIST("hEventSelectionVsCentrality"), 16 /* Below min occupancy */, centrality);
            }
            }

            if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
            return false;
            }
            if (fillHists) {
            histos.fill(HIST("hEventSelection"), 17 /* Above max occupancy */);
            if (doEventQA) {
                histos.fill(HIST("hEventSelectionVsCentrality"), 17 /* Above max occupancy */, centrality);
            }
            }
        }

        // Fetch interaction rate only if required (in order to limit ccdb calls)
        double interactionRate = (eventSelections.minIR >= 0 || eventSelections.maxIR >= 0) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource) * 1.e-3 : -1;
        if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 18 /* Below min IR */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 18 /* Below min IR */, centrality);
            }
        }

        if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 19 /* Above max IR */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 19 /* Above max IR */, centrality);
            }
        }

        if (!rctConfigurations.cfgRCTLabel.value.empty() && !rctFlagsChecker(collision)) {
            return false;
        }
        if (fillHists) {
            histos.fill(HIST("hEventSelection"), 20 /* Pass CBT condition */);
            if (doEventQA) {
            histos.fill(HIST("hEventSelectionVsCentrality"), 20 /* Pass CBT condition */, centrality);
            }
        }

        } else { // we are in Run 2
        if (eventSelections.requireSel8 && !collision.sel8()) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);

        if (eventSelections.requireSel7 && !collision.sel7()) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 2 /* sel7 collisions */);

        if (eventSelections.requireINT7 && !collision.alias_bit(kINT7)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 3 /* INT7-triggered collisions */);

        if (eventSelections.requireTriggerTVX && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 4 /* FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level */);

        if (eventSelections.rejectIncompleteDAQ && !collision.selection_bit(o2::aod::evsel::kNoIncompleteDAQ)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 5 /* Complete events according to DAQ flags */);

        if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 6 /* vertex-Z selected */);

        if (eventSelections.requireConsistentSPDAndTrackVtx && !collision.selection_bit(o2::aod::evsel::kNoInconsistentVtx)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 7 /* No inconsistency in SPD and Track vertices */);

        if (eventSelections.rejectPileupFromSPD && !collision.selection_bit(o2::aod::evsel::kNoPileupFromSPD)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 8 /* No pileup according to SPD vertexer */);

        if (eventSelections.rejectV0PFPileup && !collision.selection_bit(o2::aod::evsel::kNoV0PFPileup)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 9 /* No out-of-bunch pileup according to V0 past-future info */);

        if (eventSelections.rejectPileupInMultBins && !collision.selection_bit(o2::aod::evsel::kNoPileupInMultBins)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 10 /* No pileup according to multiplicity-differential pileup checks */);

        if (eventSelections.rejectPileupMV && !collision.selection_bit(o2::aod::evsel::kNoPileupMV)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 11 /* No pileup according to multi-vertexer */);

        if (eventSelections.rejectTPCPileup && !collision.selection_bit(o2::aod::evsel::kNoPileupTPC)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 12 /* No pileup in TPC */);

        if (eventSelections.requireNoV0MOnVsOffPileup && !collision.selection_bit(o2::aod::evsel::kNoV0MOnVsOfPileup)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 13 /* No out-of-bunch pileup according to online-vs-offline VOM correlation */);

        if (eventSelections.requireNoSPDOnVsOffPileup && !collision.selection_bit(o2::aod::evsel::kNoSPDOnVsOfPileup)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 14 /* No out-of-bunch pileup according to online-vs-offline SPD correlation */);

        if (eventSelections.requireNoSPDClsVsTklBG && !collision.selection_bit(o2::aod::evsel::kNoSPDClsVsTklBG)) {
            return false;
        }
        if (fillHists)
            histos.fill(HIST("hEventSelection"), 15 /* No beam-gas according to cluster-vs-tracklet correlation */);

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
        }
        }

        return true;
    }

    // ______________________________________________________
    // Real data processing - no MC subscription
    template <typename TCollision, typename TV0s>
    void analyzeRecoedV0sInRealData(TCollision const& collision, TV0s const& fullV0s)
    {
        // Fire up CCDB
        if ((mlConfigurations.useLambdaScores && mlConfigurations.calculateLambdaScores) ||
            v0Selections.rejectTPCsectorBoundary) {
            initCCDB(collision);
        }

        if (!isEventAccepted(collision, true)) {
        return;
        }

        float centrality = -1;
        float collisionOccupancy = -2; // -1 already taken for the case where occupancy cannot be evaluated
        double interactionRate = -1;
        // gap side
        int gapSide = -1;
        int selGapSide = -1; // -1 --> Hadronic ; 0 --> Single Gap - A side ; 1 --> Single Gap - C side ; 2 --> Double Gap - both A & C sides
        // Fill recoed event properties
        fillReconstructedEventProperties(collision, centrality, collisionOccupancy, interactionRate, gapSide, selGapSide);

        histos.fill(HIST("hInteractionRateVsOccupancy"), interactionRate, collisionOccupancy);

        // __________________________________________
        // perform main analysis
        int nLambdas = 0;
        for (auto const& v0 : fullV0s) {
            if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut || std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
                continue; // remove acceptance that's badly reproduced by MC / superfluous in future

            if (v0.v0Type() != v0Selections.v0TypeSelection && v0Selections.v0TypeSelection > -1)
                continue; // skip V0s that are not standard

            // fill AP plot for all V0s
            histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.alpha(), v0.qtarm());

            // uint64_t selMap = computeReconstructionBitmap(v0, collision, v0.yLambda(), v0.yK0Short(), v0.pt());
            uint64_t selMap = computeReconstructionBitmap(v0, collision, v0.yLambda(), v0.pt()); // Removed unneeded K0Short info

            // consider for histograms for all species
            BITSET(selMap, selConsiderLambda);
            BITSET(selMap, selPhysPrimLambda);

            // analyseCandidate(v0, v0.pt(), centrality, selMap, selGapSide, nK0Shorts, nLambdas, nAntiLambdas);
            analyseCandidate(v0, v0.pt(), centrality, selMap, selGapSide, nLambdas);
        } // end v0 loop

        // fill the histograms with the number of reconstructed Lambda per collision
        if (analyseLambda) {
            histos.fill(HIST("h2dNbrOfLambdaVsCentrality"), centrality, nLambdas);
        }
    }
    
    // Subscriping to the appropriate tables and running the code:
    // ______________________________________________________
    // Real data processing in Run 3 - no MC subscription
    void processRealDataRun3(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator const& collision, V0Candidates const& fullV0s, DauTracks const&)
    {
        analyzeRecoedV0sInRealData(collision, fullV0s);
    }

    // Kept only the process switch that is relevant for this particular work:
    PROCESS_SWITCH(lambdaInvMassTest, processRealDataRun3, "process as if real data in Run 3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
    return WorkflowSpec{adaptAnalysisTask<lambdaInvMassTest>(cfgc)};
}