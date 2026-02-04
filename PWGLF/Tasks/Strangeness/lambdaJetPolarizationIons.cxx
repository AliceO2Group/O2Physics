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
/// \file lambdaJetPolarizationIons.cxx
/// \brief Lambda and antiLambda polarization analysis task using raw data
///
/// \author Cicero Domenico Muncinelli <cicero.domenico.muncinelli@cern.ch>, Campinas State University
//
// Jet Polarization Ions task
// ================
//
// This code loops over a V0Cores table and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    cicero.domenico.muncinelli@cern.ch
//

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "Common/DataModel/Centrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
// #include "Common/DataModel/Multiplicity.h"
// #include "Common/DataModel/PIDResponseTOF.h" // Not using right now

#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/PIDResponseTPC.h"
// #include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

// Jets:
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>


////// Others not included from strangederivedbuilder.cxx:
// #include "CommonConstants/PhysicsConstants.h"
// #include "DCAFitter/DCAFitterN.h"
// #include "DataFormatsParameters/GRPMagField.h"
// #include "DataFormatsParameters/GRPObject.h"
// #include "DetectorsBase/GeometryManager.h"
// #include "DetectorsBase/Propagator.h"
// #include "Framework/O2DatabasePDGPlugin.h"
// #include "Framework/RunningWorkflowInfo.h"
// #include "Framework/StaticFor.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

//////////////////////////////////////////////
// From Youpeng's:
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Framework/ASoA.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
// #include "TProfile2D.h"
// #include <TFile.h>
// // #include <TLorentzVector.h>
// #include <TMatrixD.h>
// #include <TTree.h>
//////////////////////////////////////////////

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::aod::rctsel;

// Aliases for joined tables:
using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>; // Estimates centrality from FT0M. Will add other variations for cross-check
                                                                               // Cross-check with derivedlambdakzeroanalysis' CentEstimator enumerable. It has a lot of centrality estimators! (TODO)
// using SelV0Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentNGlobals>; // Is this better for event selection? Maybe the pre-ordered centralities are useful?
using DauTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,
                                 aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
                                // , aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>; // Not using TOF right now due to some possible mismatches
// using SimCollisions = soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>;
// using DauTracksMC = soa::Join<DaughterTracks, aod::McTrackLabels>;

// (TODO: check what is gapSide and how to implement related selections)

//////////////////////////////////////////////
struct lambdajetpolarizationions {

    // HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
    // Define histogram registries
    HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
    HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

    // master analysis switches
    Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
    Configurable<bool> analyseAntiLambda{"analyseAntiLambda", false, "process AntiLambda-like candidates"}; // Will work only with Lambdas, in a first analysis

    Configurable<bool> doPPAnalysis{"doPPAnalysis", false, "if in pp, set to true. Default is HI"};
    Configurable<std::string> irSource{"irSource", "ZNChadronic", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNChadronic)"}; // Renamed David's "ZNC hadronic" to the proper code "ZNChadronic"
    // Configurable<int> centralityEstimator{"centralityEstimator", kCentFT0C, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFT0CVariant1, 3:CentMFT, 4:CentNGlobal, 5:CentFV0A)"}; // Default is FT0M

    ///////////////////////////////////////////////
    // QA block -- not implemented! (TODO)
    // Configurable<bool> doEventQA{"doEventQA", false, "do event QA histograms"};
    // Configurable<bool> doCompleteTopoQA{"doCompleteTopoQA", false, "do topological variable QA histograms"};
    // Configurable<bool> doTPCQA{"doTPCQA", false, "do TPC QA histograms"};
    // Configurable<bool> doTOFQA{"doTOFQA", false, "do TOF QA histograms"};
    // Configurable<int> doDetectPropQA{"doDetectPropQA", 0, "do Detector/ITS map QA: 0: no, 1: 4D, 2: 5D with mass; 3: plain in 3D"};
    // Configurable<bool> doEtaPhiQA{"doEtaPhiQA", false, "do Eta/Phi QA histograms"};

    // Configurable<bool> doPlainTopoQA{"doPlainTopoQA", true, "do simple 1D QA of candidates"};
    // Configurable<float> qaMinPt{"qaMinPt", 0.0f, "minimum pT for QA plots"};
    // Configurable<float> qaMaxPt{"qaMaxPt", 1000.0f, "maximum pT for QA plots"};
    // Configurable<bool> qaCentrality{"qaCentrality", false, "qa centrality flag: check base raw values"};
    ///////////////////////////////////////////////

    ///////////////////////////////////////////////
    // MC block -- not implemented! (TODO)
    // Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};
    // Configurable<bool> doTreatPiToMuon{"doTreatPiToMuon", false, "Take pi decay into muon into account in MC"};
    // Configurable<bool> doCollisionAssociationQA{"doCollisionAssociationQA", true, "check collision association"};
    ///////////////////////////////////////////////

    //////////////////////////////////////////////
    // Manual slice by: (TODO)
    // SliceCache cache;
    // Preslice<aod::V0Datas> V0perCollision = o2::aod::v0data::collisionId;
    // Preslice<aod::CascDatas> CascperCollision = o2::aod::cascdata::collisionId;
    // Preslice<aod::KFCascDatas> KFCascperCollision = o2::aod::cascdata::collisionId;
    // Preslice<aod::TraCascDatas> TraCascperCollision = o2::aod::cascdata::collisionId;
    // Preslice<aod::McParticles> mcParticlePerMcCollision = o2::aod::mcparticle::mcCollisionId;
    // Preslice<UDCollisionsFull> udCollisionsPerCollision = o2::aod::udcollision::collisionId;
    /////////////////////////////////////////////
    
    
    // Configurable groups:
    struct : ConfigurableGroup {
        std::string prefix = "eventSelections"; // JSON group name
        Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
        Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"}; // part of sel8, actually
        Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"}; // part of sel8, actually
        Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"}; // part of sel8, actually
        Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track (Run 3 only)"};
        Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"}; // o2::aod::evsel::kIsGoodZvtxFT0vsPV. Recommended for OO
        Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF (Run 3 only)"};
        Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD (Run 3 only)"};
        Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC (Run 3 only)"}; // o2::aod::evsel::kNoSameBunchPileup. Recommended for OO
        Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds (Run 3 only)"};
        Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds (Run 3 only)"};
        Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds (Run 3 only)"};
        Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold (Run 3 only)"};
        Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF (Run 3 only)"};
        Configurable<bool> requireINEL0{"requireINEL0", true, "require INEL>0 event selection"}; // Only truly useful in pp
        Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};

        Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};

        Configurable<bool> useEvtSelInDenomEff{"useEvtSelInDenomEff", false, "Consider event selections in the recoed <-> gen collision association for the denominator (or numerator) of the acc. x eff. (or signal loss)?"};
        Configurable<bool> applyZVtxSelOnMCPV{"applyZVtxSelOnMCPV", true, "Apply Z-vtx cut on the PV of the generated collision?"}; // I see no reason as to not do this by default
        Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
        // fast check on occupancy
        Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
        Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
        // fast check on interaction rate
        Configurable<float> minIR{"minIR", -1, "minimum IR collisions"};
        Configurable<float> maxIR{"maxIR", -1, "maximum IR collisions"};
    } eventSelections;

    struct : ConfigurableGroup {
        std::string prefix = "v0Selections"; // JSON group name
        Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};

        // Selection criteria: acceptance
        Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
        Configurable<float> daughterEtaCut{"daughterEtaCut", 0.9, "max eta for daughters"}; // Default is 0.8. Changed to 0.9 to agree with jet selection. TODO: test the impact/biasing of this!

        // Standard 5 topological criteria -- Closed a bit more for the Lambda analysis
        Configurable<float> v0cospa{"v0cospa", 0.995, "min V0 CosPA"}; // Default is 0.97
        Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"}; // Default is 1.0
        // Configurable<float> dcanegtopv{"dcanegtopv", .2, "min DCA Neg To PV (cm)"}; // Default is .05
        // Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"}; // Default is .05
        // Renamed for better consistency of candidate selection (the cut is not determined by charge, but by mass and how deflected the daughter is):
        Configurable<float> dcaPionToPV{"dcaPionToPV", .2, "min DCA pion-like daughter To PV (cm)"}; // Default is .05. Suppresses pion background.
        Configurable<float> dcaProtonToPV{"dcaProtonToPV", .05, "min DCA proton-like daughter To PV (cm)"}; // Default is .05
        Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"}; // Default is  1.2
        Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
        Configurable<float> lambdaLifetimeCut{"lambdaLifetimeCut", 30., "lifetime cut (c*tau) for Lambda (cm)"}

        // invariant mass selection
        Configurable<float> compMassRejection{"compMassRejection", -1, "Competing mass rejection (GeV/#it{c}^{2})"}; // This was creating some bumps in Youpeng's inv mass spectra. Turned off for now.

        // Track quality
        Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
        Configurable<int> minITSclusters{"minITSclusters", 3, "minimum ITS clusters"}; // Default is off
        Configurable<float> minTPCrowsOverFindableClusters{"minTPCrowsOverFindableClusters", -1, "minimum nbr of TPC crossed rows over findable clusters"};
        Configurable<float> minTPCfoundOverFindableClusters{"minTPCfoundOverFindableClusters", -1, "minimum nbr of found over findable TPC clusters"};
        Configurable<float> maxFractionTPCSharedClusters{"maxFractionTPCSharedClusters", 1e+09, "maximum fraction of TPC shared clusters"};
        Configurable<float> maxITSchi2PerNcls{"maxITSchi2PerNcls", 36.0f, "maximum ITS chi2 per clusters"}; // Default is 1e+09. New values from StraInJets recommendations
        Configurable<float> maxTPCchi2PerNcls{"maxTPCchi2PerNcls", 4.0f, "maximum TPC chi2 per clusters"}; // Default is 1e+09
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
        Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 3, "tpcPidNsigmaCut"}; // Default is 5. Reduced to agree with strangenessInJetsIons
        Configurable<float> tofPidNsigmaCutLaPr{"tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
        Configurable<float> tofPidNsigmaCutLaPi{"tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
        Configurable<float> tofPidNsigmaCutK0Pi{"tofPidNsigmaCutK0Pi", 1e+6, "tofPidNsigmaCutK0Pi"};

        // PID (TOF)
        Configurable<float> maxDeltaTimeProton{"maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
        Configurable<float> maxDeltaTimePion{"maxDeltaTimePion", 1e+9, "check maximum allowed time"};
    } v0Selections;

    // Helpers for the "isTrackFarFromTPCBoundary" function:
    TF1* fPhiCutLow = new TF1("fPhiCutLow", v0Selections.phiLowCut.value.data(), 0, 100);
    TF1* fPhiCutHigh = new TF1("fPhiCutHigh", v0Selections.phiHighCut.value.data(), 0, 100);

    // Run Condition Table (RCT) configurables
    struct : ConfigurableGroup {
        std::string prefix = "rctConfigurations"; // JSON group name
        Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "", "Which detector condition requirements? (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
        Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "Include ZDC flags in the bit selection (for Pb-Pb only)"};
        Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
    } rctConfigurations;
    RCTFlagsChecker rctFlagsChecker{rctConfigurations.cfgRCTLabel.value};

    // ML SELECTIONS BLOCK -- NOT IMPLEMENTED! (TODO)

    // CCDB options
    struct : ConfigurableGroup {
        std::string prefix = "ccdbConfigurations"; // JSON group name
        Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
        Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
        Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
        Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
        Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
        Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

        // manual magnetic field:
        Configurable<bool> useCustomMagField{"useCustomMagField", false, "Use custom magnetic field value"};
        Configurable<float> customMagField{"customMagField", 5.0f, "Manually set magnetic field"};
    } ccdbConfigurations;

    // Instantiating CCDB:
    o2::ccdb::CcdbApi ccdbApi;
    Service<o2::ccdb::BasicCCDBManager> ccdb;

    // Other useful variables:
    ctpRateFetcher rateFetcher;
    int mRunNumber;
    float magField;
    std::map<std::string, std::string> metadata;
    o2::parameters::GRPMagField* grpmag = nullptr;
    
    // Histogram axes configuration:
    struct : ConfigurableGroup {
        std::string prefix = "axisConfigurations"; // JSON group name
        ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
        ConfigurableAxis axisPtXi{"axisPtXi", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for feeddown from Xi"};
        ConfigurableAxis axisPtCoarse{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
        ConfigurableAxis axisLambdaMass{"axisLambdaMass", {450, 1.08f, 1.15f}, ""}; // Default is {200, 1.101f, 1.131f}
        ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};
        ConfigurableAxis axisNch{"axisNch", {500, 0.0f, +5000.0f}, "Number of charged particles"};
        ConfigurableAxis axisIRBinning{"axisIRBinning", {500, 0, 50}, "Binning for the interaction rate (kHz)"};
        ConfigurableAxis axisMultFT0M{"axisMultFT0M", {500, 0.0f, +100000.0f}, "Multiplicity FT0M"};
        // ConfigurableAxis axisMultFT0C{"axisMultFT0C", {500, 0.0f, +10000.0f}, "Multiplicity FT0C"}; // (TODO)
        // ConfigurableAxis axisMultFV0A{"axisMultFV0A", {500, 0.0f, +100000.0f}, "Multiplicity FV0A"};

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

    // Jet selection configuration:
    // (TODO: Add a configurable to select charged jets, neutral jets, full jets, photon-tagged jets and so on)
    struct : ConfigurableGroup {
        std::string prefix = "jetConfigurations"; // JSON group name
        Configurable<double> minJetPt{"minJetPt", 30.0f, "Minimum reconstructed pt of the jet (GeV/c)"}; // Something in between pp and PbPb minima
        Configurable<double> radiusJet{"rJet", 0.4f, "Jet resolution parameter (R)"}; // (TODO: check if the JE people don't define this as a rescaled int to not lose precision for stricter selections)
        // Notice that the maximum Eta of the jet will then be 0.9 - R to keep the jet contained within the ITS+TPC barrel.
        
        Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
        Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "Jet recombination scheme: E_scheme (0), pT-scheme (1), pt2-scheme (2), WTA_pt_scheme (7)"} // See PWGJE/JetFinders/jetFinder.h for more info.
        Configurable<int> bkgSubtractionStrategy{"bkgSubtractionStrategy", 0, "Jet background subtraction: Area (0), Constituent (1)"}; // Selection bool for background subtraction strategy
        Configurable<int> jetType{"jetType", 0, "Jet type: Charged Jet (0), Full Jet (1), Photon-tagged (2), Z-tagged (3)"} // (TODO: create a reasonable track selection for photon/Z-tagged jet tracks, including detector angular acceptance parameters)

        // // Configurables from JE PWG:
        // // (TODO: check the maximum pT of jets used in my analyses! If it is way too hard, it might not be the best jet to use!)
        // Configurable<float> jetEWSPtMin{"jetEWSPtMin", 0.0, "minimum event-wise subtracted jet pT"};
        // Configurable<float> jetEWSPtMax{"jetEWSPtMax", 1000.0, "maximum event-wise subtracted jet pT"};
        // Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
        // Configurable<int> ghostRepeat{"ghostRepeat", 0, "set to 0 to gain speed if you dont need area calculation"};
        // Configurable<bool> DoTriggering{"DoTriggering", false, "used for the charged jet trigger to remove the eta constraint on the jet axis"};
        // Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

        // (TODO: Check which of these configurables might be useful for the photon-tagged and regular analyses)
        // // event level configurables
        // Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
        // Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
        // Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};
        // Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};
        // // track level configurables
        // Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
        // Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};
        // // cluster level configurables
        // Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
        // Configurable<float> clusterEtaMin{"clusterEtaMin", -0.71, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
        // Configurable<float> clusterEtaMax{"clusterEtaMax", 0.71, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
        // Configurable<float> clusterPhiMin{"clusterPhiMin", 1.39, "minimum cluster phi"};
        // Configurable<float> clusterPhiMax{"clusterPhiMax", 3.27, "maximum cluster phi"};
        // Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
        // Configurable<float> clusterTimeMin{"clusterTimeMin", -25., "minimum Cluster time (ns)"};
        // Configurable<float> clusterTimeMax{"clusterTimeMax", 25., "maximum Cluster time (ns)"};
        // Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};
        // Configurable<int> hadronicCorrectionType{"hadronicCorrectionType", 0, "0 = no correction, 1 = CorrectedOneTrack1, 2 = CorrectedOneTrack2, 3 = CorrectedAllTracks1, 4 = CorrectedAllTracks2"};
        // Configurable<bool> doEMCALEventSelection{"doEMCALEventSelection", true, "apply the selection to the event alias_bit for full and neutral jets"};
        // Configurable<bool> doEMCALEventSelectionChargedJets{"doEMCALEventSelectionChargedJets", false, "apply the selection to the event alias_bit for charged jets"};
    } jetConfigurations;

    // Track analysis parameters -- A specific group that is different from the v0Selections. In jet analyses we need to control our PseudoJet candidates!
    // (TODO: include minimal selection criteria for electrons, muons and photons)
    // Notice you do NOT need any PID for the PseudoJet candidates! Only need is to know the 4-momentum appropriately. Thus removed nsigma checks on PID
    struct : ConfigurableGroup {
        std::string prefix = "pseudoJetCandidateTrackSelections"; // JSON group name
        Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "Minimum number of TPC crossed rows"};
        Configurable<int> minITSnCls{"minITSnCls", -1, "Minimum number of ITS clusters"};
        Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "Maximum chi2 per cluster TPC"};
        Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "Maximum chi2 per cluster ITS"};
        Configurable<float> etaCut{"etaCut", 0.9f, "Maximum eta absolute value"}; // (TODO: same test as the previous 0.8 eta cut)

        Configurable<float> minCandidatePt{"minCandidatePt", 0.15f,"Minimum track pT for pseudojet candidate (GeV/c)"}; // Reduces number of pseudojet candidates from IR radiation
        // (TODO: test these minimal ratios to suppress split tracks in high occupancy PbPb or OO)
        // Configurable<float> minTPCrowsOverFindableClusters{"minTPCrowsOverFindableClusters", -1, "minimum nbr of TPC crossed rows over findable clusters"};
        // Configurable<float> minTPCfoundOverFindableClusters{"minTPCfoundOverFindableClusters", 0.8f, "minimum nbr of found over findable TPC clusters"};

        // Jets typical cuts (suppress non-primary candidates):
        Configurable<bool> doDCAcuts{"doDCAcuts", false, "Apply DCA cuts to jet candidates (use with care in this analysis!)"}
        Configurable<float> maxDCAz{"maxDCAz", 3.2f, "Max DCAz to primary vertex [cm]"};
        
        // Configurable<float> maxDCAxy{"maxDCAxy", 2.4f,"Max DCAxy to primary vertex [cm]"};
        // Using same cuts as the StrangenessInJets analysis, with a pt dependence (which may bias high pt, so use with care):
        Configurable<float> dcaxyMaxTrackPar0{"dcaxyMaxTrackPar0", 0.0105f, ""};
        Configurable<float> dcaxyMaxTrackPar1{"dcaxyMaxTrackPar1", 0.035f, ""};
        Configurable<float> dcaxyMaxTrackPar2{"dcaxyMaxTrackPar2", 1.1f, ""};
    } pseudoJetCandidateTrackSelections;

    // struct : ConfigurableGroup {
    //     std::string prefix = "jetQAConfigurations"; // JSON group name
    //     Configurable
    // } jetQAConfigurations; // (TODO)

    // Instantiate utility class for jet background subtraction
    JetBkgSubUtils backgroundSub;

    // Lambda Ring Polarization axes configurable group:
    struct : ConfigurableGroup {
        std::string prefix = "ringPolConfigurations"; // JSON group name
        
    } ringPolConfigurations; // (TODO)


    // For manual sliceBy
    Preslice<aod::V0Datas> V0perCollision = o2::aod::v0data::collisionId;
    Preslice<aod::McParticles> mcParticlePerMcCollision = o2::aod::mcparticle::mcCollisionId;

    Service<o2::framework::O2DatabasePDG> pdg;

    // std::vector<uint32_t> genLambda;
    // std::vector<uint32_t> genAntiLambda;
    // std::vector<uint32_t> genXiMinus;
    // std::vector<uint32_t> genXiPlus;
    // std::vector<uint32_t> genOmegaMinus;
    // std::vector<uint32_t> genOmegaPlus;

    void init(InitContext const&){ // (TODO: add all useful histograms here! Add flags for QA plots and the such too)
    
    }

    /////////////////////////////////////////////
    // Computation helper functions:
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
        if (phiModn > fPhiCutHigh->Eval(trackPt)) return true; // keep track
        if (phiModn < fPhiCutLow->Eval(trackPt)) return true; // keep track

        return false;  // reject track
    }

    /////////////////////////////////////////////
    // Helper functions for event and candidate selection:
    template <typename TCollision> // (TODO: add fillHists and doEventQA capabilities from derivedlambdakzeroanalysis)
    bool isEventAccepted(TCollision const& collision){ // check whether the collision passes our collision selections
        if (eventSelections.requireSel8 && !collision.sel8()) return false;
        if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) return false;
        if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) return false;
        if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) return false;

        if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) return false;

        if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) return false;
        if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) return false;
        if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) return false;
        if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) return false;
        if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) return false;
        if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) return false;
        if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) return false;
        if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) return false;
        if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) return false;
        if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) return false;
        
        if (doPPAnalysis) { // we are in pp
            if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1) return false;
            if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2) return false;
        }
        else { // Performing selections as if in Pb-Pb:
            float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
            if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) return false;
            if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) return false;

            // Fetch interaction rate only if required (in order to limit ccdb calls)
            double interactionRate = (eventSelections.minIR >= 0 || eventSelections.maxIR >= 0) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource) * 1.e-3 : -1;
            if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR) return false;
            if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR) return false;
            if (!rctConfigurations.cfgRCTLabel.value.empty() && !rctFlagsChecker(collision)) return false;
        }
        return true;
    }


    template <typename JetCandidate>
    inline bool isCandidateForPseudojetAccepted(JetCandidate const& track){
        // ITS/TPC cuts:
        if (!track.hasTPC()) return false; // Always demand TPC
        if (pseudoJetCandidateTrackSelections.minITSnCls >= 0){
            if (track.itsNCls() < pseudoJetCandidateTrackSelections.minITSnCls) return false;
        }

        if (track.tpcNClsCrossedRows() < pseudoJetCandidateTrackSelections.minNCrossedRowsTPC) return false;

        if (track.tpcChi2NCl() > pseudoJetCandidateTrackSelections.maxChi2TPC) return false;
        if (track.itsChi2NCl() > pseudoJetCandidateTrackSelections.maxChi2ITS) return false;

        // Kinematics:
        const float pt = track.pt();
        if (pt < pseudoJetCandidateTrackSelections.minCandidatePt) return false;
        if (std::fabs(track.eta()) > pseudoJetCandidateTrackSelections.etaCut) return false;

        // DCA pseudojet candidate selections:
        if (pseudoJetCandidateTrackSelections.doDCAcuts){
            // if (std::fabs(track.dcaXY()) > pseudoJetCandidateTrackSelections.maxDCAxy) return false;
            if (std::fabs(track.dcaZ()) > pseudoJetCandidateTrackSelections.maxDCAz) return false;
                // Slightly more physics-motivated cut (parametrizes the DCA resolution as function of pt)
            if (std::fabs(track.dcaXY()) > (pseudoJetCandidateTrackSelections.dcaxyMaxTrackPar0 + 
                pseudoJetCandidateTrackSelections.dcaxyMaxTrackPar1 / std::pow(pt, pseudoJetCandidateTrackSelections.dcaxyMaxTrackPar2))) return false;
        }
        return true;
    }

    // Lambda selections
    // Tests the hypothesis of the V0 being a Lambda or of it being an antiLambda.
    // (TODO: implement Armenteros cuts for extra quality control?)
    template <typename TV0, typename TCollision>
    bool passesLambdaLambdaBarHypothesis(TV0 const& v0, TCollision const& collision, bool Lambda_hypothesis){
        // Base topological variables (high rejection, low cost checks)
        if (v0.v0radius() < v0Selections.v0radius) return false;
        if (v0.v0radius() > v0Selections.v0radiusMax) return false;

        const float dcaProtonToPV = isLambda ? std::abs(v0.dcapostopv()) : std::abs(v0.dcanegtopv());
        if (std::fabs(dcaProtonToPV) < v0Selections.dcaProtonToPV) return false;
        const float dcaPionToPV = isLambda ? std::abs(v0.dcanegtopv()) : std::abs(v0.dcapostopv());
        if (std::fabs(dcaPionToPV) < v0Selections.dcaPionToPV) return false;

        if (v0.v0cosPA() < v0Selections.v0cospa) return false;
        if (v0.dcaV0daughters() > v0Selections.dcav0dau) return false;

        // rapidity
        if (std::fabs(v0.yLambda()) > v0Selections.rapidityCut) return false;
        
        // competing mass rejection
        if (std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0Selections.compMassRejection) return false;

        const auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
        const auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

        // For the cuts dcapostopv, dcanegtopv and PID cuts to be properly applied
        // while also keeping this function general enough for Lambdas and AntiLambdas,
        // we identify the roles of proton-like and pion-like for the pos and neg
        // tracks accordingly:
        auto const& protonTrack = isLambda ? posTrackExtra : negTrackExtra;
        auto const& pionTrack = isLambda ? negTrackExtra : posTrackExtra;

        // ITS quality cuts
        bool posIsFromAfterburner = posTrackExtra.hasITSAfterburner();
        bool negIsFromAfterburner = negTrackExtra.hasITSAfterburner();

        // check minimum number of ITS clusters + maximum ITS chi2 per clusters + reject or select ITS afterburner tracks if requested
        if (posTrackExtra.itsNCls() < v0Selections.minITSclusters) return false; // check minimum ITS clusters
        if (posTrackExtra.itsChi2NCl() >= v0Selections.maxITSchi2PerNcls) return false; // check maximum ITS chi2 per clusters
        if (v0Selections.rejectPosITSafterburner && posIsFromAfterburner) return false; // reject afterburner track or not
        if (v0Selections.requirePosITSafterburnerOnly && !posIsFromAfterburner) return false; // keep afterburner track or not

        if (negTrackExtra.itsNCls() < v0Selections.minITSclusters) return false; // check minimum ITS clusters
        if (negTrackExtra.itsChi2NCl() >= v0Selections.maxITSchi2PerNcls) return false; // check maximum ITS chi2 per clusters
        if (v0Selections.rejectNegITSafterburner && negIsFromAfterburner) return false; // reject afterburner track or not
        if (v0Selections.requireNegITSafterburnerOnly && !negIsFromAfterburner) return false; // keep afterburner track or not

        // TPC quality cuts
        if (posTrackExtra.tpcCrossedRows() < v0Selections.minTPCrows) return false; // check minimum TPC crossed rows
        if (posTrackExtra.tpcChi2NCl() >= v0Selections.maxTPCchi2PerNcls) return false; // check maximum TPC chi2 per clusters
        if (posTrackExtra.tpcCrossedRowsOverFindableCls() < v0Selections.minTPCrowsOverFindableClusters) return false; // check minimum fraction of TPC rows over findable
        if (posTrackExtra.tpcFoundOverFindableCls() < v0Selections.minTPCfoundOverFindableClusters) return false; // check minimum fraction of found over findable TPC clusters
        if (posTrackExtra.tpcFractionSharedCls() >= v0Selections.maxFractionTPCSharedClusters) return false; // check the maximum fraction of allowed shared TPC clusters
        if (v0Selections.rejectTPCsectorBoundary && !isTrackFarFromTPCBoundary(v0.positivept(), v0.positivephi(), +1)) return false; // reject track far from TPC sector boundary or not

        if (negTrackExtra.tpcCrossedRows() < v0Selections.minTPCrows) return false; // check minimum TPC crossed rows
        if (negTrackExtra.tpcChi2NCl() >= v0Selections.maxTPCchi2PerNcls) return false; // check maximum TPC chi2 per clusters
        if (negTrackExtra.tpcCrossedRowsOverFindableCls() < v0Selections.minTPCrowsOverFindableClusters) return false; // check minimum fraction of TPC rows over findable
        if (negTrackExtra.tpcFoundOverFindableCls() < v0Selections.minTPCfoundOverFindableClusters) return false; // check minimum fraction of found over findable TPC clusters
        if (negTrackExtra.tpcFractionSharedCls() >= v0Selections.maxFractionTPCSharedClusters) return false; // check the maximum fraction of allowed shared TPC clusters
        if (v0Selections.rejectTPCsectorBoundary && !isTrackFarFromTPCBoundary(v0.negativept(), v0.negativephi(), -1)) return false; // reject track far from TPC sector boundary or not

        // ITS only tag
        if (v0Selections.requirePosITSonly && posTrackExtra.tpcCrossedRows() > 1) return false;
        if (v0Selections.requireNegITSonly && negTrackExtra.tpcCrossedRows() > 1) return false;

        // TPC only tag
        if (v0Selections.skipTPConly && posTrackExtra.detectorMap() == o2::aod::track::TPC) return false;
        if (v0Selections.skipTPConly && negTrackExtra.detectorMap() == o2::aod::track::TPC) return false;

        ///// Expensive PID checks come last:
        // TPC PID
        if (std::fabs(protonTrack.tpcNSigmaPr()) > v0Selections.tpcPidNsigmaCut) return false;
        if (std::fabs(pionTrack.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut) return false;

        // TOF PID in DeltaT (if TOF is not available, then uses the track. If is available, uses it. In this sense, TOF is optional)
        // const bool posHasTOF = posTrackExtra.hasTOF(); // For the older version, which worked only for Lambdas
        const bool protonHasTOF = protonTrack.hasTOF();
        const bool pionHasTOF = pionTrack.hasTOF();

        // Positive track
        if (protonHasTOF && std::abs(isLambda ? v0.posTOFDeltaTLaPr()
                      : v0.negTOFDeltaTLaPr()) > v0Selections.maxDeltaTimeProton) return false;
        // Negative track
        // if (pionHasTOF && std::fabs(v0.negTOFDeltaTLaPi()) > v0Selections.maxDeltaTimePion) return false; // Older version, for Lambda only

        // TOF PID in NSigma
        // Positive track
        if (protonHasTOF && std::fabs(v0.tofNSigmaLaPr()) > v0Selections.tofPidNsigmaCutLaPr) return false;
        // Negative track
        if (pionHasTOF && std::fabs(v0.tofNSigmaLaPi()) > v0Selections.tofPidNsigmaCutLaPi) return false;

        // proper lifetime
        if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 > v0Selections.lambdaLifetimeCut) return false;
    }

    template <typename TV0>
    bool isLambdaCandidate(TV0 const& v0){

    }

    template <typename Jet>
    bool isJetAccepted(const Jet& jet){
        
    }



    


    void processDataRing(){
        if (!isEventAccepted(collision)) return;



    }





    PROCESS_SWITCH(lambdajetpolarizationions, processDataRing, "Process Ring polarization in Run 3 Data", true);
    PROCESS_SWITCH(lambdajetpolarizationions, processMCRing, "Process Ring polarization in Run 3 MC", false);
    PROCESS_SWITCH(lambdajetpolarizationions, processJetQAData, "Process Jet QA in Run 3 Data", false);
};

    






