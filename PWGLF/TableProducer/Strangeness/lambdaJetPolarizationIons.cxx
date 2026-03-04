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
// This code loops over a V0Cores table and produces standard derived
// data as output. In the post-processing stage, this analysis aims
// to measure the formation of vorticity rings in HI collisions.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    cicero.domenico.muncinelli@cern.ch
//

// O2 Framework
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

// O2 CCDB / Conditions
#include "DataFormatsParameters/GRPMagField.h"
#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>

// O2 Reconstruction Data Formats
#include <ReconstructionDataFormats/Track.h>

// O2 Common Core
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"

// O2 Common DataModel
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h" // for pp
#include "Common/DataModel/PIDResponseTPC.h"
// For PID in raw data:
// #include "Common/DataModel/PIDResponseTOF.h" // Maybe switch this around with LFStrangenessPIDTables?
// #include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

// PWGJE
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

// PWGLF
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
// For V0TOFPIDs and NSigmas getters. Better for considering the daughters as coming from V0s instead of from PV?
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/lambdaJetPolarizationIons.h"
#include "PWGLF/DataModel/mcCentrality.h"

// External Libraries (FastJet)
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

// ROOT Math
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

// Standard Library
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::aod::rctsel;

///// Aliases for joined tables
/// Collisions:
using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As,
                                aod::PVMults, aod::FT0Mults, aod::FV0Mults>; // Added PVMults to get MultNTracksPVeta1 as centrality estimator
using SelCollisionsSimple = soa::Join<aod::Collisions, aod::EvSels>;         // Simpler, for jets

/// V0s and Daughter tracks:
// using V0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
// using V0CandidatesSimple = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>; // No TOF
/// To run in RAW data:
// using V0Candidates = aod::V0Datas; // TODO: possible quicker subscription for analysis that do not require TOF.
using V0CandidatesWithTOF = soa::Join<aod::V0Datas, aod::V0TOFPIDs, aod::V0TOFNSigmas>; // Tables created by o2-analysis-lf-strangenesstofpid
// using DauTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
// Actually used subscriptions (smaller memory usage):
using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullPr>;

/// Jets:
using PseudoJetTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksDCA>; // Simpler tracks access. (Not using TracksIU and TracksCovIU. Did not use their info for now)
                                                                                                 // , aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>; // Not using TOF right now due to some possible mismatches

/// MC:
// using SimCollisions = soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>;
// using DauTracksMC = soa::Join<DaughterTracks, aod::McTrackLabels>;

enum CentEstimator {
  kCentFT0C = 0,
  kCentFT0M,
  kCentFV0A
};

enum JetAlgorithm {
  kKt = 0,
  kCambridgeAachen,
  kAntiKt
};

enum JetRecombScheme {
  kEScheme = 0,
  kPtScheme = 1,
  kPt2Scheme = 2,
  kWTAScheme = 7
};

enum JetType {
  kChargedJet = 0,
  kFullJet,
  kPhotonJet,
  kZJet
};

enum BkgSubtraction {
  kNoSubtraction = 0,
  kAreaBased,
  kConstituentBased
};

//////////////////////////////////////////////
struct lambdajetpolarizationions {

  // struct : ProducesGroup {
  // } products;
  Produces<o2::aod::RingLaV0s> tableV0s;
  Produces<o2::aod::RingJets> tableJets;
  Produces<o2::aod::RingLeadP> tableLeadParticles;
  Produces<o2::aod::RingCollisions> tableCollisions;

  // Define histogram registries:
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", false, "process AntiLambda-like candidates"}; // Will work only with Lambdas, in a first analysis

  Configurable<bool> doPPAnalysis{"doPPAnalysis", false, "if in pp, set to true. Default is HI"};
  Configurable<std::string> irSource{"irSource", "ZNChadronic", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNChadronic)"}; // Renamed David's "ZNC hadronic" to the proper code "ZNChadronic"
  Configurable<int> centralityEstimatorForQA{"centralityEstimatorForQA", kCentFT0M, "Run 3 centrality estimator (0:CentFT0C, 1:CentFT0M, 2:CentFV0A)"};  // Default is FT0M
  // (Now saving all centralities at the derived data level -- Makes them all available for consumer)
  // (But still using this variable for QA histograms)

  /////////////////////////////////////////////
  Configurable<bool> doEventQA{"doEventQA", false, "do event QA histograms"};
  // Configurable<bool> qaCentrality{"qaCentrality", false, "qa centrality flag: check base raw values"};
  Configurable<bool> doCompleteTopoQA{"doCompleteTopoQA", false, "do topological variables QA histograms"}; // Includes doPlainTopoQA from derivedlambdakzeroanalysis
  Configurable<bool> doV0KinematicQA{"doV0KinematicQA", false, "do kinematic variables QA histograms"};
  Configurable<bool> doArmenterosQA{"doArmenterosQA", false, "do Armenteros QA histograms"};
  Configurable<bool> doTPCQA{"doTPCQA", false, "do TPC QA histograms"};
  Configurable<bool> doTOFQA{"doTOFQA", false, "do TOF QA histograms"};
  Configurable<bool> doEtaPhiQA{"doEtaPhiQA", false, "do Eta/Phi QA histograms for V0s and daughters"};
  Configurable<bool> doJetKinematicsQA{"doJetKinematicsQA", false, "do pT,Eta,Phi QA histograms for jets"};
  /////////////////////////////////////////////

  // /////////////////////////////////////////////
  // MC block -- not implemented! (TODO)
  // Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};
  // Configurable<bool> doTreatPiToMuon{"doTreatPiToMuon", false, "Take pi decay into muon into account in MC"};
  // Configurable<bool> doCollisionAssociationQA{"doCollisionAssociationQA", true, "check collision association"};
  // /////////////////////////////////////////////

  // TODO: COMPLEMENTARY ANALYSES TO STUDY SPURIOUS POLARIZATION SOURCES!
  // TODO: add an event plane selection procedure to get an angle between the global polarization axis and the jet axis to uncouple polarizations?
  // TODO: (related to previous comment) if we already have event plane, also estimate v_2-caused polarization. Hydro papers indicate observable is unsensitive to this spurious polarization, but this is a perfect consistency check.
  // TODO: add a longitudinal polarization block of code to estimate other sources of polarization (and possibly study their differential dependence on the angle wrlt the jets and their rings)?
  // TODO: add a block of code that calculates polarization from Lambda fragmentation to estimate the contamination of this third source of polarization

  // Configurable groups:
  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"}; // part of sel8, actually
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"};                             // part of sel8, actually
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"};                                          // part of sel8, actually
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
    Configurable<float> rapidityCut{"rapidityCut", 1.0f, "rapidity"};
    Configurable<float> v0EtaCut{"v0EtaCut", 0.9f, "eta cut for v0"};
    Configurable<float> daughterEtaCut{"daughterEtaCut", 0.9f, "max eta for daughters"}; // Default is 0.8. Changed to 0.9 to agree with jet selection. TODO: test the impact/biasing of this!

    // Standard 5 topological criteria -- Closed a bit more for the Lambda analysis
    Configurable<float> v0cospa{"v0cospa", 0.995, "min V0 CosPA"};              // Default is 0.97
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"}; // Default is 1.0
    // Configurable<float> dcanegtopv{"dcanegtopv", .2, "min DCA Neg To PV (cm)"}; // Default is .05
    // Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"}; // Default is .05
    // Renamed for better consistency of candidate selection (the cut is not determined by charge, but by mass and how deflected the daughter is):
    Configurable<float> dcaPionToPV{"dcaPionToPV", .2, "min DCA pion-like daughter To PV (cm)"};        // Default is .05. Suppresses pion background.
    Configurable<float> dcaProtonToPV{"dcaProtonToPV", .05, "min DCA proton-like daughter To PV (cm)"}; // Default is .05
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};                            // Default is  1.2
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
    Configurable<float> lambdaLifetimeCut{"lambdaLifetimeCut", 30., "lifetime cut (c*tau) for Lambda (cm)"};

    // invariant mass selection
    Configurable<float> compMassRejection{"compMassRejection", -1, "Competing mass rejection (GeV/#it{c}^{2})"}; // This was creating bumps in the pp analysis code's invariant mass. Turned off for now.

    // Track quality
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"minITSclusters", 3, "minimum ITS clusters"}; // Default is off
    Configurable<float> minTPCrowsOverFindableClusters{"minTPCrowsOverFindableClusters", -1, "minimum nbr of TPC crossed rows over findable clusters"};
    Configurable<float> minTPCfoundOverFindableClusters{"minTPCfoundOverFindableClusters", -1, "minimum nbr of found over findable TPC clusters"};
    Configurable<float> maxFractionTPCSharedClusters{"maxFractionTPCSharedClusters", 1e+09, "maximum fraction of TPC shared clusters"};
    Configurable<float> maxITSchi2PerNcls{"maxITSchi2PerNcls", 36.0f, "maximum ITS chi2 per clusters"}; // Default is 1e+09. New values from StraInJets recommendations
    Configurable<float> maxTPCchi2PerNcls{"maxTPCchi2PerNcls", 4.0f, "maximum TPC chi2 per clusters"};  // Default is 1e+09
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
    ConfigurableAxis axisPhi{"axisPhi", {50, 0.0f, constants::math::TwoPI}, "Azimuth angle (rad)"};
    ConfigurableAxis axisPhiMod{"axisPhiMod", {100, 0.0f, constants::math::TwoPI / 18}, "Azimuth angle wrt TPC sector (rad.)"};
    ConfigurableAxis axisEta{"axisEta", {50, -1.0f, 1.0f}, "#eta"};
    ConfigurableAxis axisRapidity{"axisRapidity", {50, -1.0f, 1.0f}, "y"};
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

    // Jet QA axes:
    ConfigurableAxis JetsPerEvent{"JetsPerEvent", {20, 0, 20}, "Jets per event"};

    ConfigurableAxis axisLeadingParticlePt{"axisLeadingParticlePt", {200, 0.f, 200.f}, "Leading particle p_{T} (GeV/c)"}; // Simpler version!
    ConfigurableAxis axisJetPt{"axisJetPt", {200, 0.f, 200.f}, "Jet p_{t} (GeV)"};
    ConfigurableAxis axisCosTheta{"axisCosTheta", {50, -1.f, 1.f}, "cos(#Delta #theta_{jet})"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {50, -constants::math::PI, constants::math::PI}, "#Delta #phi"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {50, -1.5f, 1.5f}, "#Delta #phi"};                     // Calculated as twice the subtraction "eta_max=0.9 - R_min=0.2", with a margin
    ConfigurableAxis axisDeltaR{"axisDeltaR", {50, 0, 3.5f}, "#Delta R"};                                // From 0 to about the maximum Delta R possible with R = 0.2
    ConfigurableAxis axisEnergy{"axisEnergy", {200, 0.f, 200.f}, "E_{jet} (GeV) (#pi mass hypothesis)"}; // Jet energy is not that well defined here, due to track mass hypothesis being of pions! This is just to include px,py,pz in full!
  } axisConfigurations;

  // Jet selection configuration:
  // (TODO: create a reasonable track selection for full, photon, and Z-tagged jet tracks, including detector angular acceptance parameters for EMCal)
  struct : ConfigurableGroup {
    std::string prefix = "jetConfigurations";                                                        // JSON group name
    Configurable<double> minJetPt{"minJetPt", 30.0f, "Minimum reconstructed pt of the jet (GeV/c)"}; // Something in between pp and PbPb minima. Change for bkgSubtraction true or false!
    Configurable<double> radiusJet{"radiusJet", 0.4f, "Jet resolution parameter (R)"};               // (TODO: check if the JE people don't define this as a rescaled int to not lose precision for stricter selections)
    // Notice that the maximum Eta of the jet will then be 0.9 - R to keep the jet contained within the ITS+TPC barrel.

    Configurable<int> jetAlgorithm{"jetAlgorithm", kAntiKt, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
    Configurable<int> jetRecombScheme{"jetRecombScheme", kEScheme, "Jet recombination scheme: 0: E_scheme, 1: pT-scheme, 2: pt2-scheme, 7: WTA_pt_scheme"};    // See PWGJE/JetFinders/jetFinder.h for more info.
    Configurable<int> bkgSubtraction{"bkgSubtraction", kNoSubtraction, "Jet background subtraction: No subtraction (false), Area (true), Constituent (TODO)"}; // Selection bool for background subtraction strategy
    Configurable<float> GhostedAreaSpecRapidity{"GhostedAreaSpecRapidity", 1.1, "Max ghost particle rapidity for jet area estimates"};                         // At least 1.0 for tracks and jets within the |eta| < 0.9 window of ITS+TPC
                                                                                                                                                               // Using an enum for readability:
    Configurable<int> jetType{"jetType", kChargedJet, "Jet type: 0: Charged Jet, 1: Full Jet, 2: Photon-tagged, 3: Z-tagged"};                                 // TODO: implement full, photon and Z jets
    // (TODO: check the maximum pT of jets used in my analyses! If it is way too hard, it might not be the best jet to use!)

    // (TODO: Check which of these configurables might be useful for the photon-tagged and regular analyses)
    // // Configurables from JE PWG:
    // Configurable<float> jetEWSPtMin{"jetEWSPtMin", 0.0, "minimum event-wise subtracted jet pT"};
    // Configurable<float> jetEWSPtMax{"jetEWSPtMax", 1000.0, "maximum event-wise subtracted jet pT"};
    // Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
    // Configurable<int> ghostRepeat{"ghostRepeat", 0, "set to 0 to gain speed if you dont need area calculation"};
    // Configurable<bool> DoTriggering{"DoTriggering", false, "used for the charged jet trigger to remove the eta constraint on the jet axis"};
    // Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
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

    Configurable<float> minLeadParticlePt{"minLeadParticlePt", 2.0f, "Minimum Pt for a lead track to be considered a valid proxy for a jet"}; // For OO, about 2 or 3 should be enough (z~0.3 of jet), and for PbPb maybe 8 GeV
  } jetConfigurations;

  // Creating a short map to make sure the proper FastJet enums are used (safeguard against possible updates in FastJet indices):
  fastjet::JetAlgorithm mapFJAlgorithm(int algoIdx)
  {
    switch (algoIdx) {
      case 0:
        return fastjet::kt_algorithm;
      case 1:
        return fastjet::cambridge_algorithm;
      case 2:
        return fastjet::antikt_algorithm;
      default:
        throw std::invalid_argument("Unknown jet algorithm");
    }
  }
  fastjet::RecombinationScheme mapFJRecombScheme(int schemeIdx)
  {
    switch (schemeIdx) {
      case 0:
        return fastjet::E_scheme;
      case 1:
        return fastjet::pt_scheme;
      case 2:
        return fastjet::pt2_scheme;
      case 7:
        return fastjet::WTA_pt_scheme;
      default:
        throw std::invalid_argument("Unknown recombination scheme");
    }
  }

  // Track analysis parameters -- A specific group that is different from the v0Selections. In jet analyses we need to control our PseudoJet candidates!
  // (TODO: include minimal selection criteria for electrons, muons and photons)
  // Notice you do NOT need any PID for the PseudoJet candidates! Only need is to know the 4-momentum appropriately. Thus removed nsigma checks on PID
  struct : ConfigurableGroup {
    std::string prefix = "pseudoJetCandidateTrackSelections"; // JSON group name
    Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "Minimum number of TPC crossed rows"};
    Configurable<int> minITSnCls{"minITSnCls", -1, "Minimum number of ITS clusters"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 5.0f, "Maximum chi2 per cluster TPC"}; // Loose cuts for pseudojet candidate selection
    Configurable<float> maxChi2ITS{"maxChi2ITS", 40.0f, "Maximum chi2 per cluster ITS"};
    Configurable<float> etaCut{"etaCut", 0.9f, "Maximum eta absolute value"}; // (TODO: same test as the previous 0.8 eta cut)

    Configurable<float> minCandidatePt{"minCandidatePt", 0.15f, "Minimum track pt for pseudojet candidate (GeV/c)"}; // Reduces number of pseudojet candidates from IR radiation
    // (TODO: test these minimal ratios to suppress split tracks in high occupancy PbPb or OO)
    // Configurable<float> minTPCrowsOverFindableClusters{"minTPCrowsOverFindableClusters", -1, "minimum nbr of TPC crossed rows over findable clusters"};
    // Configurable<float> minTPCfoundOverFindableClusters{"minTPCfoundOverFindableClusters", 0.8f, "minimum nbr of found over findable TPC clusters"};

    // Jets typical cuts (suppress non-primary candidates):
    Configurable<bool> doDCAcuts{"doDCAcuts", false, "Apply DCA cuts to jet candidates (biases towards primary-vertex/prompt hadron jets)"};
    Configurable<float> maxDCAz{"maxDCAz", 3.2f, "Max DCAz to primary vertex [cm] (remove pileup influence)"};

    // Configurable<float> maxDCAxy{"maxDCAxy", 2.4f,"Max DCAxy to primary vertex [cm]"};
    // Using same cuts as the StrangenessInJets analysis, with a pt dependence (which may bias high pt, so use with care):
    Configurable<float> dcaxyMaxTrackPar0{"dcaxyMaxTrackPar0", 0.0105f, "Asymptotic DCA resolution at high pt [cm]"};
    Configurable<float> dcaxyMaxTrackPar1{"dcaxyMaxTrackPar1", 0.035f, "Low pt multiple-scattering term for DCA resolution [cm*(GeV/c)^Par2]"};
    Configurable<float> dcaxyMaxTrackPar2{"dcaxyMaxTrackPar2", 1.1f, "Exponent of pt dependence of DCA resolution"};
  } pseudoJetCandidateTrackSelections;

  JetBkgSubUtils backgroundSub;

  void init(InitContext const&)
  {
    // setting CCDB service
    ccdb->setURL(ccdbConfigurations.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    // Initialise the RCTFlagsChecker
    rctFlagsChecker.init(rctConfigurations.cfgRCTLabel.value, rctConfigurations.cfgCheckZDC, rctConfigurations.cfgTreatLimitedAcceptanceAsBad);

    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1D, {{23, -0.5f, +20.5f}});
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
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(22, "hasRingJet");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(23, "hasRingV0");
    // (notice we lack a hasRingJet AND hasRingV0 bin because the tasks run separately on all events!)
    // (this QA number can be obtained at derived data level with ease)

    histos.add("Centrality/hEventCentrality", "hEventCentrality", kTH1D, {{101, 0.0f, 101.0f}});
    histos.add("Centrality/hCentralityVsNch", "hCentralityVsNch", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisNch});
    if (doEventQA) {
      histos.add("hEventSelectionVsCentrality", "hEventSelectionVsCentrality", kTH2D, {{23, -0.5f, +20.5f}, {101, 0.0f, 101.0f}});
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
      histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(22, "hasRingJet");
      histos.get<TH2>(HIST("hEventSelectionVsCentrality"))->GetXaxis()->SetBinLabel(23, "hasRingV0");

      // Centrality:
      histos.add("Centrality/hEventCentVsMultFT0M", "hEventCentVsMultFT0M", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisMultFT0M});
      histos.add("Centrality/hEventCentVsMultFT0C", "hEventCentVsMultFT0C", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisMultFT0C});
      histos.add("Centrality/hEventCentVsMultFV0A", "hEventCentVsMultFV0A", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisMultFV0A});
      histos.add("Centrality/hEventMultFT0CvsMultFV0A", "hEventMultFT0CvsMultFV0A", kTH2D, {axisConfigurations.axisMultFT0C, axisConfigurations.axisMultFV0A});
    }

    histos.add("hEventPVz", "hEventPVz", kTH1D, {{100, -20.0f, +20.0f}});
    histos.add("hCentralityVsPVz", "hCentralityVsPVz", kTH2D, {{101, 0.0f, 101.0f}, {100, -20.0f, +20.0f}});

    // (TODO: add MC centrality vs PVz histos)

    histos.add("hEventOccupancy", "hEventOccupancy", kTH1D, {axisConfigurations.axisOccupancy});
    histos.add("hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisOccupancy});
    histos.add("hInteractionRate", "hInteractionRate", kTH1D, {axisConfigurations.axisIRBinning});
    histos.add("hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2D, {{101, 0.0f, 101.0f}, axisConfigurations.axisIRBinning});
    histos.add("hInteractionRateVsOccupancy", "hInteractionRateVsOccupancy", kTH2D, {axisConfigurations.axisIRBinning, axisConfigurations.axisOccupancy});

    // for QA and test purposes
    // auto hRawCentrality = histos.add<TH1>("Centrality/hRawCentrality", "hRawCentrality", kTH1D, {axisConfigurations.axisRawCentrality});

    // for (int ii = 1; ii < 101; ii++) {
    //     float value = 100.5f - static_cast<float>(ii);
    //     hRawCentrality->SetBinContent(ii, value);
    // }

    //////////////////////////////////////////////////////////////
    /// Lambda / AntiLambda V0 selection QA
    //////////////////////////////////////////////////////////////
    struct CutLabel {
      std::string label;
      bool enabled;
    }; // A method of hiding labels of selections which were not used!
    std::vector<CutLabel> v0LambdaSelectionLabels = {
      {"All V0 candidates", true},
      {"V0 radius (min)", true},
      {"V0 radius (max)", true},
      {"V0 cosPA", true},
      {"DCA_{V0 daughters}", true},
      {"|y_{#Lambda}|", v0Selections.rapidityCut > 0.f},
      {"K^{0}_{S} mass rejection", v0Selections.compMassRejection >= 0.f},
      {"ITS clusters (pos)", v0Selections.minITSclusters > 0},
      {"ITS #chi^{2}/N_{cls} (pos)", v0Selections.maxITSchi2PerNcls < 1e8},
      {"Reject ITS afterburner (pos)", v0Selections.rejectPosITSafterburner},
      {"Require ITS afterburner (pos)", v0Selections.requirePosITSafterburnerOnly},
      {"ITS clusters (neg)", v0Selections.minITSclusters > 0},
      {"ITS #chi^{2}/N_{cls} (neg)", v0Selections.maxITSchi2PerNcls < 1e8},
      {"Reject ITS afterburner (neg)", v0Selections.rejectNegITSafterburner},
      {"Require ITS afterburner (neg)", v0Selections.requireNegITSafterburnerOnly},
      {"TPC crossed rows (pos)", v0Selections.minTPCrows > 0},
      {"TPC #chi^{2}/N_{cls} (pos)", v0Selections.maxTPCchi2PerNcls < 1e8},
      {"TPC rows / findable (pos)", v0Selections.minTPCrowsOverFindableClusters >= 0},
      {"TPC found / findable (pos)", v0Selections.minTPCfoundOverFindableClusters >= 0},
      {"TPC shared clusters (pos)", v0Selections.maxFractionTPCSharedClusters < 1e8},
      {"TPC sector boundary (pos)", v0Selections.rejectTPCsectorBoundary},
      {"TPC crossed rows (neg)", v0Selections.minTPCrows > 0},
      {"TPC #chi^{2}/N_{cls} (neg)", v0Selections.maxTPCchi2PerNcls < 1e8},
      {"TPC rows / findable (neg)", v0Selections.minTPCrowsOverFindableClusters >= 0},
      {"TPC found / findable (neg)", v0Selections.minTPCfoundOverFindableClusters >= 0},
      {"TPC shared clusters (neg)", v0Selections.maxFractionTPCSharedClusters < 1e8},
      {"TPC sector boundary (neg)", v0Selections.rejectTPCsectorBoundary},
      {"Require ITS-only (pos)", v0Selections.requirePosITSonly},
      {"Require ITS-only (neg)", v0Selections.requireNegITSonly},
      {"Reject TPC-only (pos)", v0Selections.skipTPConly},
      {"Reject TPC-only (neg)", v0Selections.skipTPConly},
    }; // First, the labels that are hypothesis-agnostic
    // Adding the Lambda or AntiLambda hypothesis labels as needed:
    auto addHypothesis = [&](bool isLambda, bool analysisEnabled) {
      if (!analysisEnabled)
        return; // i.e., don't add these labels if not analyzing said particle type
      std::string p = isLambda ? "#Lambda: " : "#bar{#Lambda}: ";
      v0LambdaSelectionLabels.insert(v0LambdaSelectionLabels.end(), {{p + "DCA_{p} to PV", true},
                                                                     {p + "DCA_{#pi} to PV", true},
                                                                     {p + "TPC PID p", v0Selections.tpcPidNsigmaCut > 0},
                                                                     {p + "TPC PID #pi", v0Selections.tpcPidNsigmaCut > 0},
                                                                     {p + "TOF #Delta t p", v0Selections.maxDeltaTimeProton < 1e+9},
                                                                     {p + "TOF #Delta t #pi", v0Selections.maxDeltaTimePion < 1e+9},
                                                                     {p + "TOF PID p", v0Selections.tofPidNsigmaCutLaPr < 1e+6},
                                                                     {p + "TOF PID #pi", v0Selections.tofPidNsigmaCutLaPi < 1e+6},
                                                                     {p + "c#tau", v0Selections.lambdaLifetimeCut > 0}});
    };
    constexpr bool Lambda = true;      // Some constexpr to make it more readable (works at compile level)
    constexpr bool AntiLambda = false; // "false" is just a flag for this addHypothesis function! It just means fill "AntiLambda" labels
    addHypothesis(Lambda, analyseLambda);
    addHypothesis(AntiLambda, analyseAntiLambda);

    auto hSelectionV0s = histos.add<TH1>("GeneralQA/hSelectionV0s", "V0 #rightarrow #Lambda / #bar{#Lambda} selection flow", kTH1D,
                                         {{(int)v0LambdaSelectionLabels.size(), -0.5, (double)v0LambdaSelectionLabels.size() - 0.5}});
    for (size_t i = 0; i < v0LambdaSelectionLabels.size(); ++i) {
      auto lbl = v0LambdaSelectionLabels[i].label;
      if (!v0LambdaSelectionLabels[i].enabled)
        lbl = "#color[16]{(off) " + lbl + "}";
      hSelectionV0s->GetXaxis()->SetBinLabel(i + 1, lbl.c_str()); // First non-underflow bin is bin 1
    }
    ////////////////////////////////////////////////
    // Jet track candidate selection flow (analogous to hSelectionV0s):
    // Each label's "enabled" flag reflects whether the corresponding configurable
    // makes that cut active, so disabled stages are shown in grey in the output.
    std::vector<CutLabel> jetTrackSelectionLabels = {
      {"All track candidates", true},
      {"ITS clusters (min)", pseudoJetCandidateTrackSelections.minITSnCls >= 0},
      {"TPC crossed rows (min)", pseudoJetCandidateTrackSelections.minNCrossedRowsTPC > 0},
      {"TPC #chi^{2}/N_{cls} (max)", pseudoJetCandidateTrackSelections.maxChi2TPC < 1.e8f},
      {"ITS #chi^{2}/N_{cls} (max)", pseudoJetCandidateTrackSelections.maxChi2ITS < 1.e8f},
      {"p_{T} min", pseudoJetCandidateTrackSelections.minCandidatePt > 0.f},
      {"|#eta| cut", pseudoJetCandidateTrackSelections.etaCut < 1.5f},
      {"DCA_{z} to PV", pseudoJetCandidateTrackSelections.doDCAcuts.value},
      {"DCA_{xy} to PV (parametric)", pseudoJetCandidateTrackSelections.doDCAcuts.value},
    };
    auto hSelectionJetTracks = histos.add<TH1>("GeneralQA/hSelectionJetTracks", "Charged pseudojet candidate selection flow", kTH1D,
                                               {{(int)jetTrackSelectionLabels.size(), -0.5, (double)jetTrackSelectionLabels.size() - 0.5}});
    for (size_t i = 0; i < jetTrackSelectionLabels.size(); ++i) {
      auto lbl = jetTrackSelectionLabels[i].label;
      if (!jetTrackSelectionLabels[i].enabled)
        lbl = "#color[16]{(off) " + lbl + "}";
      hSelectionJetTracks->GetXaxis()->SetBinLabel(i + 1, lbl.c_str());
    }
    ////////////////////////////////////////////////

    // Histograms versus mass:
    if (analyseLambda) {
      histos.add("Lambda/h2dNbrOfLambdaVsCentrality", "h2dNbrOfLambdaVsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("Lambda/h3dMassLambda", "h3dMassLambda", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      // Non-UPC info
      histos.add("Lambda/h3dMassLambdaHadronic", "h3dMassLambdaHadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      if (doTPCQA) {
        histos.add("Lambda/h3dPosNsigmaTPC", "h3dPosNsigmaTPC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("Lambda/h3dNegNsigmaTPC", "h3dNegNsigmaTPC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("Lambda/h3dPosTPCsignal", "h3dPosTPCsignal", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("Lambda/h3dNegTPCsignal", "h3dNegTPCsignal", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("Lambda/h3dPosNsigmaTPCvsTrackPtot", "h3dPosNsigmaTPCvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("Lambda/h3dNegNsigmaTPCvsTrackPtot", "h3dNegNsigmaTPCvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("Lambda/h3dPosTPCsignalVsTrackPtot", "h3dPosTPCsignalVsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("Lambda/h3dNegTPCsignalVsTrackPtot", "h3dNegTPCsignalVsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("Lambda/h3dPosNsigmaTPCvsTrackPt", "h3dPosNsigmaTPCvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("Lambda/h3dNegNsigmaTPCvsTrackPt", "h3dNegNsigmaTPCvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("Lambda/h3dPosTPCsignalVsTrackPt", "h3dPosTPCsignalVsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("Lambda/h3dNegTPCsignalVsTrackPt", "h3dNegTPCsignalVsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
      }
      if (doTOFQA) {
        histos.add("Lambda/h3dPosNsigmaTOF", "h3dPosNsigmaTOF", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("Lambda/h3dNegNsigmaTOF", "h3dNegNsigmaTOF", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("Lambda/h3dPosTOFdeltaT", "h3dPosTOFdeltaT", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("Lambda/h3dNegTOFdeltaT", "h3dNegTOFdeltaT", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("Lambda/h3dPosNsigmaTOFvsTrackPtot", "h3dPosNsigmaTOFvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("Lambda/h3dNegNsigmaTOFvsTrackPtot", "h3dNegNsigmaTOFvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("Lambda/h3dPosTOFdeltaTvsTrackPtot", "h3dPosTOFdeltaTvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("Lambda/h3dNegTOFdeltaTvsTrackPtot", "h3dNegTOFdeltaTvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("Lambda/h3dPosNsigmaTOFvsTrackPt", "h3dPosNsigmaTOFvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("Lambda/h3dNegNsigmaTOFvsTrackPt", "h3dNegNsigmaTOFvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("Lambda/h3dPosTOFdeltaTvsTrackPt", "h3dPosTOFdeltaTvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("Lambda/h3dNegTOFdeltaTvsTrackPt", "h3dNegTOFdeltaTvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
      }
      // (TODO: add collision association capabilities in MC)
      if (doEtaPhiQA) {
        histos.add("Lambda/h5dV0PhiVsEta", "h5dV0PhiVsEta", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPhi, axisConfigurations.axisEta});
        histos.add("Lambda/h5dPosPhiVsEta", "h5dPosPhiVsEta", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPhi, axisConfigurations.axisEta});
        histos.add("Lambda/h5dNegPhiVsEta", "h5dNegPhiVsEta", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPhi, axisConfigurations.axisEta});
      }
    }
    if (analyseAntiLambda) {
      histos.add("AntiLambda/h2dNbrOfAntiLambdaVsCentrality", "h2dNbrOfAntiLambdaVsCentrality", kTH2D, {axisConfigurations.axisCentrality, {10, -0.5f, 9.5f}});
      histos.add("AntiLambda/h3dMassAntiLambda", "h3dMassAntiLambda", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      // Non-UPC info
      histos.add("AntiLambda/h3dMassAntiLambdaHadronic", "h3dMassAntiLambdaHadronic", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
      if (doTPCQA) {
        histos.add("AntiLambda/h3dPosNsigmaTPC", "h3dPosNsigmaTPC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("AntiLambda/h3dNegNsigmaTPC", "h3dNegNsigmaTPC", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("AntiLambda/h3dPosTPCsignal", "h3dPosTPCsignal", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("AntiLambda/h3dNegTPCsignal", "h3dNegTPCsignal", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("AntiLambda/h3dPosNsigmaTPCvsTrackPtot", "h3dPosNsigmaTPCvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("AntiLambda/h3dNegNsigmaTPCvsTrackPtot", "h3dNegNsigmaTPCvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("AntiLambda/h3dPosTPCsignalVsTrackPtot", "h3dPosTPCsignalVsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("AntiLambda/h3dNegTPCsignalVsTrackPtot", "h3dNegTPCsignalVsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("AntiLambda/h3dPosNsigmaTPCvsTrackPt", "h3dPosNsigmaTPCvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("AntiLambda/h3dNegNsigmaTPCvsTrackPt", "h3dNegNsigmaTPCvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTPC});
        histos.add("AntiLambda/h3dPosTPCsignalVsTrackPt", "h3dPosTPCsignalVsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
        histos.add("AntiLambda/h3dNegTPCsignalVsTrackPt", "h3dNegTPCsignalVsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTPCsignal});
      }
      if (doTOFQA) {
        histos.add("AntiLambda/h3dPosNsigmaTOF", "h3dPosNsigmaTOF", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("AntiLambda/h3dNegNsigmaTOF", "h3dNegNsigmaTOF", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("AntiLambda/h3dPosTOFdeltaT", "h3dPosTOFdeltaT", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("AntiLambda/h3dNegTOFdeltaT", "h3dNegTOFdeltaT", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("AntiLambda/h3dPosNsigmaTOFvsTrackPtot", "h3dPosNsigmaTOFvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("AntiLambda/h3dNegNsigmaTOFvsTrackPtot", "h3dNegNsigmaTOFvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("AntiLambda/h3dPosTOFdeltaTvsTrackPtot", "h3dPosTOFdeltaTvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("AntiLambda/h3dNegTOFdeltaTvsTrackPtot", "h3dNegTOFdeltaTvsTrackPtot", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("AntiLambda/h3dPosNsigmaTOFvsTrackPt", "h3dPosNsigmaTOFvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("AntiLambda/h3dNegNsigmaTOFvsTrackPt", "h3dNegNsigmaTOFvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisNsigmaTOF});
        histos.add("AntiLambda/h3dPosTOFdeltaTvsTrackPt", "h3dPosTOFdeltaTvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
        histos.add("AntiLambda/h3dNegTOFdeltaTvsTrackPt", "h3dNegTOFdeltaTvsTrackPt", kTH3D, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisTOFdeltaT});
      }
      if (doEtaPhiQA) {
        histos.add("AntiLambda/h5dV0PhiVsEta", "h5dV0PhiVsEta", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPhi, axisConfigurations.axisEta});
        histos.add("AntiLambda/h5dPosPhiVsEta", "h5dPosPhiVsEta", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPhi, axisConfigurations.axisEta});
        histos.add("AntiLambda/h5dNegPhiVsEta", "h5dNegPhiVsEta", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPhi, axisConfigurations.axisEta});
      }
    }

    if (analyseLambda) {
      histos.add("hMassLambda", "hMassLambda", kTH1D, {axisConfigurations.axisLambdaMass});
      histos.add("Lambda/hLambdasPerEvent", "hLambdasPerEvent", kTH1D, {{15, 0, 15}});
    }
    if (analyseAntiLambda) {
      histos.add("hMassAntiLambda", "hMassAntiLambda", kTH1D, {axisConfigurations.axisLambdaMass});
      histos.add("AntiLambda/hAntiLambdasPerEvent", "hAntiLambdasPerEvent", kTH1D, {{15, 0, 15}});
    };
    if (analyseLambda && analyseAntiLambda) {
      histos.add("hAmbiguousLambdaCandidates", "hAmbiguousLambdaCandidates", kTH1D, {{1, 0, 1}});
      histos.add("hAmbiguousPerEvent", "hAmbiguousPerEvent", kTH1D, {{15, 0, 15}});
    }

    // QA histograms if requested
    if (doV0KinematicQA) {
      if (analyseLambda) {
        // --- Basic kinematics ---
        histos.add("V0KinematicQA/Lambda/hPt", "Lambda p_{T}", kTH1D, {axisConfigurations.axisPt});
        histos.add("V0KinematicQA/Lambda/hY", "Lambda rapidity", kTH1D, {axisConfigurations.axisRapidity});
        histos.add("V0KinematicQA/Lambda/hPhi", "Lambda #varphi", kTH1D, {axisConfigurations.axisPhi});
        // --- Mass correlations ---
        histos.add("V0KinematicQA/Lambda/hMassVsPt", "Lambda mass vs p_{T}", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        histos.add("V0KinematicQA/Lambda/hMassVsY", "Lambda mass vs y", kTH2D, {axisConfigurations.axisRapidity, axisConfigurations.axisLambdaMass});
        histos.add("V0KinematicQA/Lambda/hMassVsPhi", "Lambda mass vs #varphi", kTH2D, {axisConfigurations.axisPhi, axisConfigurations.axisLambdaMass});
        // --- Kinematic correlations ---
        histos.add("V0KinematicQA/Lambda/hYVsPt", "Lambda y vs p_{T}", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisRapidity});
        histos.add("V0KinematicQA/Lambda/hPhiVsPt", "Lambda #varphi vs p_{T}", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPhi});
      }
      if (analyseAntiLambda) {
        // --- Basic kinematics ---
        histos.add("V0KinematicQA/AntiLambda/hPt", "AntiLambda p_{T}", kTH1D, {axisConfigurations.axisPt});
        histos.add("V0KinematicQA/AntiLambda/hY", "AntiLambda rapidity", kTH1D, {axisConfigurations.axisRapidity});
        histos.add("V0KinematicQA/AntiLambda/hPhi", "AntiLambda #varphi", kTH1D, {axisConfigurations.axisPhi});
        // --- Mass correlations ---
        histos.add("V0KinematicQA/AntiLambda/hMassVsPt", "AntiLambda mass vs p_{T}", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisLambdaMass});
        histos.add("V0KinematicQA/AntiLambda/hMassVsY", "AntiLambda mass vs y", kTH2D, {axisConfigurations.axisRapidity, axisConfigurations.axisLambdaMass});
        histos.add("V0KinematicQA/AntiLambda/hMassVsPhi", "AntiLambda mass vs #varphi", kTH2D, {axisConfigurations.axisPhi, axisConfigurations.axisLambdaMass});
        // --- Kinematic correlations ---
        histos.add("V0KinematicQA/AntiLambda/hYVsPt", "AntiLambda y vs p_{T}", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisRapidity});
        histos.add("V0KinematicQA/AntiLambda/hPhiVsPt", "AntiLambda #varphi vs p_{T}", kTH2D, {axisConfigurations.axisPt, axisConfigurations.axisPhi});
      }
    }

    if (doCompleteTopoQA) {
      if (analyseLambda) {
        histos.add("Lambda/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAtoPV});
        histos.add("Lambda/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAtoPV});
        histos.add("Lambda/h4dDCADaughters", "h4dDCADaughters", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAdau});
        histos.add("Lambda/h4dPointingAngle", "h4dPointingAngle", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPointingAngle});
        histos.add("Lambda/h4dV0Radius", "h4dV0Radius", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisV0Radius});
      }
      if (analyseAntiLambda) {
        histos.add("AntiLambda/h4dPosDCAToPV", "h4dPosDCAToPV", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAtoPV});
        histos.add("AntiLambda/h4dNegDCAToPV", "h4dNegDCAToPV", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAtoPV});
        histos.add("AntiLambda/h4dDCADaughters", "h4dDCADaughters", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisDCAdau});
        histos.add("AntiLambda/h4dPointingAngle", "h4dPointingAngle", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisPointingAngle});
        histos.add("AntiLambda/h4dV0Radius", "h4dV0Radius", kTHnD, {axisConfigurations.axisCentrality, axisConfigurations.axisPtCoarse, axisConfigurations.axisLambdaMass, axisConfigurations.axisV0Radius});
      }

      // For all received candidates:
      histos.add("V0KinematicQA/hPosDCAToPV", "hPosDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
      histos.add("V0KinematicQA/hNegDCAToPV", "hNegDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
      histos.add("V0KinematicQA/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
      histos.add("V0KinematicQA/hPointingAngle", "hPointingAngle", kTH1D, {axisConfigurations.axisPointingAngle});
      histos.add("V0KinematicQA/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
      histos.add("V0KinematicQA/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
      histos.add("V0KinematicQA/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
      histos.add("V0KinematicQA/h2dPositivePtVsPhi", "h2dPositivePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
      histos.add("V0KinematicQA/h2dNegativePtVsPhi", "h2dNegativePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
      if (analyseLambda) {
        histos.add("Lambda/hPosDCAToPV", "hPosDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("Lambda/hNegDCAToPV", "hNegDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("Lambda/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("Lambda/hPointingAngle", "hPointingAngle", kTH1D, {axisConfigurations.axisPointingAngle});
        histos.add("Lambda/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("Lambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("Lambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("Lambda/h2dPositivePtVsPhi", "h2dPositivePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
        histos.add("Lambda/h2dNegativePtVsPhi", "h2dNegativePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
      }
      if (analyseAntiLambda) {
        histos.add("AntiLambda/hPosDCAToPV", "hPosDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("AntiLambda/hNegDCAToPV", "hNegDCAToPV", kTH1D, {axisConfigurations.axisDCAtoPV});
        histos.add("AntiLambda/hDCADaughters", "hDCADaughters", kTH1D, {axisConfigurations.axisDCAdau});
        histos.add("AntiLambda/hPointingAngle", "hPointingAngle", kTH1D, {axisConfigurations.axisPointingAngle});
        histos.add("AntiLambda/hV0Radius", "hV0Radius", kTH1D, {axisConfigurations.axisV0Radius});
        histos.add("AntiLambda/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("AntiLambda/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2D, {axisConfigurations.axisTPCrows, axisConfigurations.axisITSclus});
        histos.add("AntiLambda/h2dPositivePtVsPhi", "h2dPositivePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
        histos.add("AntiLambda/h2dNegativePtVsPhi", "h2dNegativePtVsPhi", kTH2D, {axisConfigurations.axisPtCoarse, axisConfigurations.axisPhiMod});
      }
    }

    // Check ambiguous candidates in AP space:
    histos.add("GeneralQA/h2dArmenterosAll", "h2dArmenterosAll", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosKinematicSelected", "h2dArmenterosKinematicSelected", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosFullSelected", "h2dArmenterosFullSelected", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosFullSelectedLambda", "h2dArmenterosFullSelectedLambda", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosFullSelectedAntiLambda", "h2dArmenterosFullSelectedAntiLambda", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});
    histos.add("GeneralQA/h2dArmenterosFullSelectedAmbiguous", "h2dArmenterosFullSelectedAmbiguous", kTH2D, {axisConfigurations.axisAPAlpha, axisConfigurations.axisAPQt});

    // Jets histograms:
    // Histogram that needs to be present even out of QA:
    histos.add("hEventsWithJet", "hEventsWithJet", kTH1D, {{1, 0, 1}});
    histos.add("hJetsPerEvent", "hJetsPerEvent", kTH1D, {axisConfigurations.JetsPerEvent});
    // counter of events with jet (could be interesting to compare with the minimum pT cut or between the background subtraction vs no background subtraction cases)
    // number of jets per event
    if (doJetKinematicsQA) {
      histos.add("JetKinematicsQA/hJetPt", "hJetPt", kTH1D, {axisConfigurations.axisJetPt});
      histos.add("JetKinematicsQA/hJetEta", "hJetEta", kTH1D, {axisConfigurations.axisEta});
      histos.add("JetKinematicsQA/hJetPhi", "hJetPhi", kTH1D, {axisConfigurations.axisPhi});

      histos.add("JetKinematicsQA/hCosThetaToLeadingJet", "hCosThetaToLeadingJet", kTH1D, {axisConfigurations.axisCosTheta});
      histos.add("JetKinematicsQA/hDeltaPhiToLeadingJet", "hDeltaPhiToLeadingJet", kTH1D, {axisConfigurations.axisDeltaPhi});
      histos.add("JetKinematicsQA/hDeltaEtaToLeadingJet", "hDeltaEtaToLeadingJet", kTH1D, {axisConfigurations.axisDeltaEta});
      histos.add("JetKinematicsQA/hDeltaRToLeadingJet", "hDeltaRToLeadingJet", kTH1D, {axisConfigurations.axisDeltaR});

      histos.add("JetKinematicsQA/hLeadingJetPt", "hLeadingJetPt", kTH1D, {axisConfigurations.axisJetPt});
      histos.add("JetKinematicsQA/hLeadingJetEta", "hLeadingJetEta", kTH1D, {axisConfigurations.axisEta});
      histos.add("JetKinematicsQA/hLeadingJetPhi", "hLeadingJetPhi", kTH1D, {axisConfigurations.axisPhi});

      // 2D correlations:
      histos.add("JetKinematicsQA/h2dJetsPerEventvsLeadJetPt", "h2dJetsPerEventvsLeadJetPt", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisJetPt});
      histos.add("JetKinematicsQA/h2dJetsPerEventvsJetPt", "h2dJetsPerEventvsJetPt", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisJetPt});
      histos.add("JetKinematicsQA/h2dCosThetaToLeadvsDeltaPhiToLead", "h2dCosThetaToLeadvsDeltaPhiToLead", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaPhi});
      histos.add("JetKinematicsQA/h2dCosThetaToLeadvsDeltaEtaToLead", "h2dCosThetaToLeadvsDeltaEtaToLead", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaEta});
      histos.add("JetKinematicsQA/h2dCosThetaToLeadvsDeltaRToLead", "h2dCosThetaToLeadvsDeltaRToLead", kTH2D, {axisConfigurations.axisCosTheta, axisConfigurations.axisDeltaR});
      histos.add("JetKinematicsQA/h2dDeltaPhiToLeadvsDeltaEtaToLead", "h2dDeltaPhiToLeadvsDeltaEtaToLead", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDeltaEta}); // to see existence of back-to-back jets, and in which window

      // Comparisons to jet energy:
      histos.add("JetKinematicsQA/h2dJetPtvsDeltaPhiToLead", "h2dJetPtvsDeltaPhiToLead", kTH2D, {axisConfigurations.axisJetPt, axisConfigurations.axisDeltaPhi});
      histos.add("JetKinematicsQA/h2dJetEnergyvsDeltaPhiToLead", "h2dJetEnergyvsDeltaPhiToLead", kTH2D, {axisConfigurations.axisEnergy, axisConfigurations.axisDeltaPhi});
      histos.add("JetKinematicsQA/h2dJetEnergyvsCosThetaToLead", "h2dJetEnergyvsCosThetaToLead", kTH2D, {axisConfigurations.axisEnergy, axisConfigurations.axisCosTheta});

      // Jets per event vs correlation to lead jet
      histos.add("JetKinematicsQA/h2dJetsPerEventvsDeltaPhiToLead", "h2dJetsPerEventvsDeltaPhiToLead", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisDeltaPhi});
      histos.add("JetKinematicsQA/h2dJetsPerEventvsDeltaEtaToLead", "h2dJetsPerEventvsDeltaEtaToLead", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisDeltaEta});
      histos.add("JetKinematicsQA/h2dJetsPerEventvsCosThetaToLead", "h2dJetsPerEventvsCosThetaToLead", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisCosTheta});

      ////////////////////////////
      // Leading particle 1D QA:
      histos.add("JetVsLeadingParticleQA/hLeadingParticlePt", "hLeadingParticlePt", kTH1D, {axisConfigurations.axisLeadingParticlePt});
      histos.add("JetVsLeadingParticleQA/hLeadingParticleEta", "hLeadingParticleEta", kTH1D, {axisConfigurations.axisEta});
      histos.add("JetVsLeadingParticleQA/hLeadingParticlePhi", "hLeadingParticlePhi", kTH1D, {axisConfigurations.axisPhi});

      // 1D correlations to lead jet:
      histos.add("JetVsLeadingParticleQA/hCosThetaLeadParticleToJet", "hCosThetaLeadParticleToJet", kTH1D, {axisConfigurations.axisCosTheta});
      histos.add("JetVsLeadingParticleQA/hDeltaPhiLeadParticleToJet", "hDeltaPhiLeadParticleToJet", kTH1D, {axisConfigurations.axisDeltaPhi});
      histos.add("JetVsLeadingParticleQA/hDeltaEtaToLeadParticleToJet", "hDeltaEtaToLeadParticleToJet", kTH1D, {axisConfigurations.axisDeltaEta});

      // Leading particle correlations:
      histos.add("JetVsLeadingParticleQA/h2dDeltaPhiParticleToLeadvsDeltaEtaParticleToLead", "h2dDeltaPhiParticleToLeadvsDeltaEtaParticleToLead", kTH2D, {axisConfigurations.axisDeltaPhi, axisConfigurations.axisDeltaEta});

      // Jets-per-event vs particle-to-lead correlations:
      histos.add("JetVsLeadingParticleQA/h2dJetsPerEventvsDeltaPhiParticleToLead", "h2dJetsPerEventvsDeltaPhiParticleToLead", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisDeltaPhi});
      histos.add("JetVsLeadingParticleQA/h2dJetsPerEventvsDeltaEtaParticleToLead", "h2dJetsPerEventvsDeltaEtaParticleToLead", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisDeltaEta});
      histos.add("JetVsLeadingParticleQA/h2dJetsPerEventvsCosThetaParticleToLead", "h2dJetsPerEventvsCosThetaParticleToLead", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisCosTheta});

      // Main "Leading jet vs leading particle" correlations:
      histos.add("JetVsLeadingParticleQA/h2dJetsPerEventvsLeadParticlePt", "h2dJetsPerEventvsLeadParticlePt", kTH2D, {axisConfigurations.JetsPerEvent, axisConfigurations.axisLeadingParticlePt});
      histos.add("JetVsLeadingParticleQA/h2dLeadJetPtvsLeadParticlePt", "h2dLeadJetPtvsLeadParticlePt", kTH2D, {axisConfigurations.axisJetPt, axisConfigurations.axisLeadingParticlePt});
      histos.add("JetVsLeadingParticleQA/h2dLeadJetPtvsCosThetaParticleToLead", "h2dLeadJetPtvsCosThetaParticleToLead", kTH2D, {axisConfigurations.axisJetPt, axisConfigurations.axisCosTheta});
      histos.add("JetVsLeadingParticleQA/h2dLeadParticlePtvsCosThetaParticleToLead", "h2dLeadParticlePtvsCosThetaParticleToLead", kTH2D, {axisConfigurations.axisLeadingParticlePt, axisConfigurations.axisCosTheta});
      histos.add("JetVsLeadingParticleQA/h2dLeadJetPtvsDeltaPhiParticleToLead", "h2dLeadJetPtvsDeltaPhiParticleToLead", kTH2D, {axisConfigurations.axisJetPt, axisConfigurations.axisDeltaPhi});
      histos.add("JetVsLeadingParticleQA/h2dLeadParticlePtvsDeltaPhiParticleToLead", "h2dLeadParticlePtvsDeltaPhiParticleToLead", kTH2D, {axisConfigurations.axisLeadingParticlePt, axisConfigurations.axisDeltaPhi});
    }

    // inspect histogram sizes, please
    histos.print();
  }

  template <typename TCollision>
  auto getCentrality(TCollision const& collision)
  {
    if (centralityEstimatorForQA == kCentFT0M)
      return collision.centFT0M();
    else if (centralityEstimatorForQA == kCentFT0C)
      return collision.centFT0C();
    else if (centralityEstimatorForQA == kCentFV0A)
      return collision.centFV0A();
    return -1.f;
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();

    // Fetching magnetic field as requested
    // In case override, don't proceed, please - no CCDB access required
    if (ccdbConfigurations.useCustomMagField) {
      magField = ccdbConfigurations.customMagField;
    } else {
      grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, mRunNumber);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField and " << ccdbConfigurations.grpPath << " of object GRPObject for run " << mRunNumber;
      }
      // Fetch magnetic field from ccdb for current bc
      magField = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for run " << mRunNumber << " with magnetic field of " << magField << " kZG";
    }
  }

  // Minimal helper to fill the hSelectionV0s histogram without having to deal with bins by myself
  // (CAUTION! If you change selection order, change this too!)
  struct V0SelectionFlowCounter {        // Using struct to keep internal bin counter over different functions
    int binValue = -1;                   // Starts at x=-1, which will go to bin 0 (underflow) in the definition of hSelectionV0s
                                         // Made it like this because we use ++binValue when filling, so the first filled
                                         // bin will always be x=0 due to operator precedence.
    HistogramRegistry* histos = nullptr; // Had to pass the histos group to this struct, as it was not visible to the members of this struct

    void resetForNewV0() { binValue = -1; }
    void fill() { histos->fill(HIST("GeneralQA/hSelectionV0s"), ++binValue); } // Hardcoded hSelectionV0s histogram, as it will not change. Increments before filling, by default
  };
  V0SelectionFlowCounter V0SelCounter{0, &histos};

  // Minimal helper to fill hSelectionJetTracks, mirroring V0SelectionFlowCounter.
  // Reset once per track candidate, fill once per passed cut stage.
  struct JetTrackSelectionFlowCounter {
    int binValue = -1; // Same convention as V0: starts at -1, first fill goes to bin x=0
    HistogramRegistry* histos = nullptr;
    void resetForNewTrack() { binValue = -1; }
    void fill() { histos->fill(HIST("GeneralQA/hSelectionJetTracks"), ++binValue); }
  };
  JetTrackSelectionFlowCounter JetTrackSelCounter{0, &histos};

  // Short inlined helper to simplify QA
  inline void fillEventSelectionQA(int bin, float centrality)
  {
    histos.fill(HIST("hEventSelection"), bin);
    histos.fill(HIST("hEventSelectionVsCentrality"), bin, centrality);
  }

  // Fill reconstructed event centrality information
  // Based off fillReconstructedEventProperties, but optimized to avoid re-accessing information already present on isEventAccepted!
  template <typename TCollision>
  void fillCentralityProperties(TCollision const& collision, float centrality)
  {
    // if (qaCentrality) {
    //     auto hRawCentrality = histos.get<TH1>(HIST("Centrality/hRawCentrality"));
    //     centrality = hRawCentrality->GetBinContent(hRawCentrality->FindBin(doPPAnalysis ? collision.multFT0A() + collision.multFT0C() : collision.multFT0C()));
    // }
    histos.fill(HIST("Centrality/hEventCentrality"), centrality);
    histos.fill(HIST("Centrality/hCentralityVsNch"), centrality, collision.multNTracksPVeta1());
    if (doEventQA) {
      histos.fill(HIST("Centrality/hEventCentVsMultFT0M"), collision.centFT0M(), collision.multFT0A() + collision.multFT0C());
      histos.fill(HIST("Centrality/hEventCentVsMultFT0C"), collision.centFT0C(), collision.multFT0C());
      histos.fill(HIST("Centrality/hEventCentVsMultFV0A"), collision.centFV0A(), collision.multFV0A());
      histos.fill(HIST("Centrality/hEventMultFT0CvsMultFV0A"), collision.multFT0C(), collision.multFV0A());
    }
    return;
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
    if (phiModn > fPhiCutHigh->Eval(trackPt))
      return true; // keep track
    if (phiModn < fPhiCutLow->Eval(trackPt))
      return true; // keep track

    return false; // reject track
  }

  inline float cosThetaJets(const fastjet::PseudoJet& a, const fastjet::PseudoJet& b)
  {
    const double dot = a.px() * b.px() + a.py() * b.py() + a.pz() * b.pz();
    const double magA = std::sqrt(a.px() * a.px() + a.py() * a.py() + a.pz() * a.pz());
    const double magB = std::sqrt(b.px() * b.px() + b.py() * b.py() + b.pz() * b.pz());
    return dot / (magA * magB);
  }

  /////////////////////////////////////////////
  // Helper functions for event and candidate selection:
  template <typename TCollision, typename TBC>
  bool isEventAccepted(TCollision const& collision, TBC const& bc, float centrality, bool fillHists)
  {                       // check whether the collision passes our collision selections
    int selectionIdx = 0; // To loop over QA histograms. First bin is already filled: first call will already increment this index (not actually the bin index, but a value in the X axis).
    if (eventSelections.requireSel8 && !collision.sel8())
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);

    const float collisionPVz = collision.posZ();
    if (std::abs(collisionPVz) > eventSelections.maxZVtxPosition)
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);

    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);
    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict))
      return false;
    if (fillHists)
      fillEventSelectionQA(++selectionIdx, centrality);

    if (doPPAnalysis) { // we are in pp
      if constexpr (requires { collision.multNTracksPVeta1(); }) {
        // Only considers compiling this block when the collision type actually
        // has multNTracksPVeta1(). This is done to reduce collision-table
        // subscriptions in the jet processing function.
        // This is a compile-time check: since the function is templated, it
        // is instantiated separately for Jets and V0s, and this block will be
        // properly compiled for each use case and table subscription automatically!
        if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < 1)
          return false;
        if (fillHists)
          fillEventSelectionQA(++selectionIdx, centrality);
        if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < 2)
          return false;
        if (fillHists)
          fillEventSelectionQA(++selectionIdx, centrality);
      }
    } else { // Performing selections as if in Pb-Pb:
      const float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy)
        return false;
      if (fillHists)
        fillEventSelectionQA(++selectionIdx, centrality);
      if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy)
        return false;
      if (fillHists)
        fillEventSelectionQA(++selectionIdx, centrality);

      // Fetch interaction rate only if required (in order to limit ccdb calls)
      const double interactionRate = (eventSelections.minIR >= 0 || eventSelections.maxIR >= 0) ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3 : -1;
      if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR)
        return false;
      if (fillHists)
        fillEventSelectionQA(++selectionIdx, centrality);
      if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR)
        return false;
      if (fillHists)
        fillEventSelectionQA(++selectionIdx, centrality);
      if (!rctConfigurations.cfgRCTLabel.value.empty() && !rctFlagsChecker(collision))
        return false;
      if (fillHists)
        fillEventSelectionQA(++selectionIdx, centrality);

      // Filling histograms previously filled in fillReconstructedEventProperties here, to avoid re-accessing data:
      if (fillHists) {
        histos.fill(HIST("hEventOccupancy"), collisionOccupancy);
        histos.fill(HIST("hCentralityVsOccupancy"), centrality, collisionOccupancy);
        histos.fill(HIST("hInteractionRate"), interactionRate);
        histos.fill(HIST("hCentralityVsInteractionRate"), centrality, interactionRate);
        histos.fill(HIST("hInteractionRateVsOccupancy"), interactionRate, collisionOccupancy);
      }
    }

    if (fillHists) {
      histos.fill(HIST("hCentralityVsPVz"), centrality, collisionPVz);
      histos.fill(HIST("hEventPVz"), collisionPVz);
    }
    return true;
  }

  template <typename JetCandidate>
  bool isCandidateForChargedPseudojetAccepted(JetCandidate const& track)
  { // (TODO: add an equivalent for photon jets and Z jets, which don't consider charged particles)
    // if (track.sign() == 0) return false; // Tracks are always either positive or negative, at least in TPC and ITS, which are the ones used (not looking at photon-jets right now)
    // ITS/TPC cuts:
    if (pseudoJetCandidateTrackSelections.minITSnCls >= 0) {
      if (track.itsNCls() < pseudoJetCandidateTrackSelections.minITSnCls)
        return false;
    }
    JetTrackSelCounter.fill(); // bin: ITS clusters (min)

    if (track.tpcNClsCrossedRows() < pseudoJetCandidateTrackSelections.minNCrossedRowsTPC)
      return false;
    JetTrackSelCounter.fill();

    if (track.tpcChi2NCl() > pseudoJetCandidateTrackSelections.maxChi2TPC)
      return false;
    JetTrackSelCounter.fill();
    if (track.itsChi2NCl() > pseudoJetCandidateTrackSelections.maxChi2ITS)
      return false;
    JetTrackSelCounter.fill();

    // Kinematics:
    const float pt = track.pt();
    if (pt < pseudoJetCandidateTrackSelections.minCandidatePt)
      return false;
    JetTrackSelCounter.fill();
    if (std::fabs(track.eta()) > pseudoJetCandidateTrackSelections.etaCut)
      return false;
    JetTrackSelCounter.fill();

    // DCA pseudojet candidate selections -- These select primary vertex particles for the jet:
    if (pseudoJetCandidateTrackSelections.doDCAcuts) {
      // if (std::fabs(track.dcaXY()) > pseudoJetCandidateTrackSelections.maxDCAxy) return false;
      if (std::fabs(track.dcaZ()) > pseudoJetCandidateTrackSelections.maxDCAz)
        return false;
      JetTrackSelCounter.fill();
      // Slightly more physics-motivated cut (parametrizes the DCA resolution as function of pt)
      if (std::fabs(track.dcaXY()) > (pseudoJetCandidateTrackSelections.dcaxyMaxTrackPar0 +
                                      pseudoJetCandidateTrackSelections.dcaxyMaxTrackPar1 / std::pow(pt, pseudoJetCandidateTrackSelections.dcaxyMaxTrackPar2)))
        return false;
      JetTrackSelCounter.fill();
    }
    return true;
  }

  // Lambda selections:
  template <typename TV0>
  bool passesGenericV0Cuts(TV0 const& v0)
  {
    // Base topological variables (high rejection, low cost checks)
    if (v0.v0radius() < v0Selections.v0radius)
      return false;
    V0SelCounter.fill();
    if (v0.v0radius() > v0Selections.v0radiusMax)
      return false;
    V0SelCounter.fill();
    if (v0.v0cosPA() < v0Selections.v0cospa)
      return false;
    V0SelCounter.fill();
    if (v0.dcaV0daughters() > v0Selections.dcav0dau)
      return false;
    V0SelCounter.fill();

    // pseudorapidity cuts:
    if (std::fabs(v0.yLambda()) > v0Selections.rapidityCut)
      return false;
    // if (std::fabs(v0.eta()) > v0Selections.v0EtaCut) return false;
    V0SelCounter.fill();
    // if (std::fabs(v0.eta()) > v0Selections.daughterEtaCut) return false; // (TODO: properly consider this in daughter selection!)

    // competing mass rejection (if compMassRejection < 0, this cut does nothing)
    if (std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0Selections.compMassRejection)
      return false;
    V0SelCounter.fill();

    const auto posTrackExtra = v0.template posTrack_as<DauTracks>(); // (TODO: is it worth it to cache these transformations outside of the function? They are reused in the Lambda hypothesis checks)
    const auto negTrackExtra = v0.template negTrack_as<DauTracks>();

    // ITS quality cuts
    bool posIsFromAfterburner = posTrackExtra.isITSAfterburner();
    bool negIsFromAfterburner = negTrackExtra.isITSAfterburner();

    // check minimum number of ITS clusters + maximum ITS chi2 per clusters + reject or select ITS afterburner tracks if requested
    if (posTrackExtra.itsNCls() < v0Selections.minITSclusters)
      return false; // check minimum ITS clusters
    V0SelCounter.fill();
    if (posTrackExtra.itsChi2NCl() >= v0Selections.maxITSchi2PerNcls)
      return false; // check maximum ITS chi2 per clusters
    V0SelCounter.fill();
    if (v0Selections.rejectPosITSafterburner && posIsFromAfterburner)
      return false; // reject afterburner track or not
    V0SelCounter.fill();
    if (v0Selections.requirePosITSafterburnerOnly && !posIsFromAfterburner)
      return false; // keep afterburner track or not
    V0SelCounter.fill();

    if (negTrackExtra.itsNCls() < v0Selections.minITSclusters)
      return false; // check minimum ITS clusters
    V0SelCounter.fill();
    if (negTrackExtra.itsChi2NCl() >= v0Selections.maxITSchi2PerNcls)
      return false; // check maximum ITS chi2 per clusters
    V0SelCounter.fill();
    if (v0Selections.rejectNegITSafterburner && negIsFromAfterburner)
      return false; // reject afterburner track or not
    V0SelCounter.fill();
    if (v0Selections.requireNegITSafterburnerOnly && !negIsFromAfterburner)
      return false; // keep afterburner track or not
    V0SelCounter.fill();

    // TPC quality cuts
    if (posTrackExtra.tpcNClsCrossedRows() < v0Selections.minTPCrows)
      return false; // check minimum TPC crossed rows
    V0SelCounter.fill();
    if (posTrackExtra.tpcChi2NCl() >= v0Selections.maxTPCchi2PerNcls)
      return false; // check maximum TPC chi2 per clusters
    V0SelCounter.fill();
    if (posTrackExtra.tpcCrossedRowsOverFindableCls() < v0Selections.minTPCrowsOverFindableClusters)
      return false; // check minimum fraction of TPC rows over findable
    V0SelCounter.fill();
    if (posTrackExtra.tpcFoundOverFindableCls() < v0Selections.minTPCfoundOverFindableClusters)
      return false; // check minimum fraction of found over findable TPC clusters
    V0SelCounter.fill();
    if (posTrackExtra.tpcFractionSharedCls() >= v0Selections.maxFractionTPCSharedClusters)
      return false; // check the maximum fraction of allowed shared TPC clusters
    V0SelCounter.fill();
    if (v0Selections.rejectTPCsectorBoundary && !isTrackFarFromTPCBoundary(v0.positivept(), v0.positivephi(), +1))
      return false; // reject track far from TPC sector boundary or not
    V0SelCounter.fill();

    if (negTrackExtra.tpcNClsCrossedRows() < v0Selections.minTPCrows)
      return false; // check minimum TPC crossed rows
    V0SelCounter.fill();
    if (negTrackExtra.tpcChi2NCl() >= v0Selections.maxTPCchi2PerNcls)
      return false; // check maximum TPC chi2 per clusters
    V0SelCounter.fill();
    if (negTrackExtra.tpcCrossedRowsOverFindableCls() < v0Selections.minTPCrowsOverFindableClusters)
      return false; // check minimum fraction of TPC rows over findable
    V0SelCounter.fill();
    if (negTrackExtra.tpcFoundOverFindableCls() < v0Selections.minTPCfoundOverFindableClusters)
      return false; // check minimum fraction of found over findable TPC clusters
    V0SelCounter.fill();
    if (negTrackExtra.tpcFractionSharedCls() >= v0Selections.maxFractionTPCSharedClusters)
      return false; // check the maximum fraction of allowed shared TPC clusters
    V0SelCounter.fill();
    if (v0Selections.rejectTPCsectorBoundary && !isTrackFarFromTPCBoundary(v0.negativept(), v0.negativephi(), -1))
      return false; // reject track far from TPC sector boundary or not
    V0SelCounter.fill();

    // ITS only tag
    if (v0Selections.requirePosITSonly && posTrackExtra.tpcNClsCrossedRows() > 1)
      return false;
    V0SelCounter.fill();
    if (v0Selections.requireNegITSonly && negTrackExtra.tpcNClsCrossedRows() > 1)
      return false;
    V0SelCounter.fill();

    // TPC only tag
    if (v0Selections.skipTPConly && posTrackExtra.detectorMap() == o2::aod::track::TPC)
      return false;
    V0SelCounter.fill();
    if (v0Selections.skipTPConly && negTrackExtra.detectorMap() == o2::aod::track::TPC)
      return false;
    V0SelCounter.fill();

    return true;
  }

  // Tests the hypothesis of the V0 being a Lambda or of it being an antiLambda.
  template <typename TV0, typename TCollision>
  bool passesLambdaLambdaBarHypothesis(TV0 const& v0, TCollision const& collision, bool Lambda_hypothesis)
  {
    // Remaining topological cuts that were charge-dependent:
    // (there is no real gain in doing a looser version of these in the passesGenericV0Cuts function.
    //  The DCA check will be done anyways and is very unexpensive)
    // (even though they are high rejection, they demand a Lambda vs AntiLambda hypothesis, so they
    //  only appear here...)
    const float dcaProtonToPV = Lambda_hypothesis ? std::abs(v0.dcapostopv()) : std::abs(v0.dcanegtopv());
    if (dcaProtonToPV < v0Selections.dcaProtonToPV)
      return false;
    V0SelCounter.fill();
    const float dcaPionToPV = Lambda_hypothesis ? std::abs(v0.dcanegtopv()) : std::abs(v0.dcapostopv()); // Checks Lambda_hypothesis twice, but compiler can handle it cleanly.
    if (dcaPionToPV < v0Selections.dcaPionToPV)
      return false;
    V0SelCounter.fill();

    const auto posTrackExtra = v0.template posTrack_as<DauTracks>();
    const auto negTrackExtra = v0.template negTrack_as<DauTracks>();

    // For the PID cuts to be properly applied while also keeping this function
    // general enough for Lambdas and AntiLambdas, we identify the roles of
    // proton-like and pion-like for the pos and neg tracks accordingly:
    auto const& protonTrack = Lambda_hypothesis ? posTrackExtra : negTrackExtra;
    auto const& pionTrack = Lambda_hypothesis ? negTrackExtra : posTrackExtra;

    ///// Expensive PID checks come last:
    // TPC PID
    if (std::fabs(protonTrack.tpcNSigmaPr()) > v0Selections.tpcPidNsigmaCut)
      return false;
    V0SelCounter.fill();
    if (std::fabs(pionTrack.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCut)
      return false;
    V0SelCounter.fill();

    // TOF PID in DeltaT (if TOF is not available, then uses the track. If is available, uses it. In this sense, TOF is optional)
    // const bool posHasTOF = posTrackExtra.hasTOF(); // For the older version, which worked only for Lambdas
    const bool protonHasTOF = protonTrack.hasTOF(); // Should work even without PIDResponseTOF.h, as it is a TracksExtra property
    const bool pionHasTOF = pionTrack.hasTOF();

    // Proton-like track
    if (protonHasTOF && std::abs(Lambda_hypothesis ? v0.posTOFDeltaTLaPr() : v0.negTOFDeltaTLaPr()) > v0Selections.maxDeltaTimeProton)
      return false;
    V0SelCounter.fill();
    // Pion-like track
    if (pionHasTOF && std::abs(Lambda_hypothesis ? v0.negTOFDeltaTLaPi() : v0.posTOFDeltaTLaPi()) > v0Selections.maxDeltaTimePion)
      return false;
    V0SelCounter.fill();

    // TOF PID in NSigma (TODO: add asymmetric NSigma windows for purity tuning?)
    // Proton-like track
    if (protonHasTOF && std::fabs(v0.tofNSigmaLaPr()) > v0Selections.tofPidNsigmaCutLaPr)
      return false; // (No need to select which candidate is which with the Lambda_hypothesis. Automatically done already!)
    V0SelCounter.fill();
    // Pion-like track
    if (pionHasTOF && std::fabs(v0.tofNSigmaLaPi()) > v0Selections.tofPidNsigmaCutLaPi)
      return false;
    V0SelCounter.fill();

    // (CAUTION!) You cannot use the getter for raw data's PIDResponseTOF.h instead of LFStrangenessPIDTables.h (as below)
    // If you do use, TOF will just try to identify that track as a proton, instead of using the correct path length from the
    // V0s PV-DCA and the such! In other words, it is a naive estimator of TOF PID, because it does not correct for the V0
    // mother's travel time and considers all tracks as if they came from the PV!
    // if (protonHasTOF && std::fabs(protonTrack.tofNSigmaPr()) > v0Selections.tofPidNsigmaCutLaPr) return false;
    // To properly use the LFStrangenessPIDTables version, you need to call o2-analysis-lf-strangenesstofpid too.

    // proper lifetime
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 > v0Selections.lambdaLifetimeCut)
      return false;
    V0SelCounter.fill();

    return true;
  }

  // Function to help distinguish ambiguous candidates (via Armenteros) that pass both
  // the Lambda_hypothesis true (i.e., a Lambda) or false (i.e., an AntiLambda) checks
  // (This function is only called in about 1-3% of the Lambda-Like V0s which remain ambiguous after all other cuts)
  // int isCandidateArmenterosLambda(const float alpha, const float qt){
  //     // Remove K0s band
  //     if (std::abs(alpha) < v0Selections.armK0AlphaThreshold && qt < v0Selections.armK0QtThreshold) return kIsArmenterosK0;
  //     // std::abs(alpha) < 0.2 && qt < 0.1
  //     if (std::abs(alpha) < v0Selections.armMinAlpha) return kArmenterosAmbiguous;
  //     // std::abs(alpha) < 0.01f
  //     // Lambda selection
  //     if (alpha > 0) return kIsArmenterosLambda;
  //     else return kIsArmenterosAntiLambda;
  // }

  // TODO: another possible check that could be done (if not implemented already inside mLambda() getters)
  // template <typename TV0>
  // int isCandidateMassLambda(TV0 const& v0) {
  //     float m1 = v0.mLambda();      // proton=positive
  //     float m2 = v0.mAntiLambda();  // proton=negative
  //     float d = std::abs(m1 - mLambdaTrue) - std::abs(m2 - mLambdaTrue);
  //     if (d < 0.f) return +1;   // Lambda
  //     else return -1;   // AntiLambda
  // }

  void processJetsData(SelCollisionsSimple::iterator const& collision, PseudoJetTracks const& tracks, aod::BCsWithTimestamps const& bcs)
  {                           // Uses BCsWithTimestamps to get timestamps for rejectTPCsectorBoundary
    float centrality = -1.0f; // Just a placeholder

    // For event QA the last two indices never change for NEv_withJets and NEv_withV0s
    // (Not the best way to initialize this: runs once per collision! TODO: think of a better way to do it)
    int lastBinEvSel = histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->GetNbins();
    bool validJetAlreadyFound = false; // Do not fill Event QA more than once

    auto bc = bcs.iteratorAt(collision.bcId()); // Got the iteratorAt() idea from O2Physics/PWGUD/Core/UDHelpers.h
    if (!isEventAccepted(collision, bc, centrality, false))
      return; // Uses return instead of continue, as there is no explicit loop here
    const uint64_t collIdx = collision.globalIndex();

    // Loop over reconstructed tracks:
    std::vector<fastjet::PseudoJet> fjParticles;
    int leadingParticleIdx = -1; // Initialized as -1, but could leave it unitialized as well. We reject any invalid events where this could pose a problem (e.g., pT<=0)
    float leadingParticlePt = 0;
    for (auto const& track : tracks) {
      JetTrackSelCounter.resetForNewTrack(); // reset bin counter for this candidate
      JetTrackSelCounter.fill();             // bin: "All track candidates"

      // Require that tracks pass selection criteria
      if (!isCandidateForChargedPseudojetAccepted(track))
        continue;

      // Constructing pseudojet candidates vector:
      // Using pion mass as hypothesis for track energy estimate (before PID, all particles treated as if with the same invariant mass)
      // (TODO: study the possibility of using identified PseudoJet candidates for this estimate)
      fastjet::PseudoJet candidate(track.px(), track.py(), track.pz(), track.energy(o2::constants::physics::MassPionCharged));
      fjParticles.emplace_back(candidate);

      // Calculating leading particle
      float pt = candidate.pt();
      if (pt > leadingParticlePt) {
        leadingParticlePt = pt;
        leadingParticleIdx = fjParticles.size() - 1;
      }
    }
    // Reject empty events
    if (fjParticles.size() < 1)
      return;

    auto const& leadingParticle = fjParticles[leadingParticleIdx];
    if (leadingParticle.pt() > jetConfigurations.minLeadParticlePt) { // If not, leading particle is probably a bad proxy
      tableLeadParticles(collIdx, leadingParticle.pt(), leadingParticle.eta(), leadingParticle.phi());
    }

    // Start jet clusterization:
    // Cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(mapFJAlgorithm(jetConfigurations.jetAlgorithm), jetConfigurations.radiusJet, mapFJRecombScheme(jetConfigurations.jetRecombScheme));
    // std::vector<float> jets_pt, jets_eta, jets_phi; // Not worth it to store 4-vectors: the tracks assume pion mass hypothesis, so energy and rapidity are not right.
    if (jetConfigurations.bkgSubtraction == kAreaBased) {
      fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(jetConfigurations.GhostedAreaSpecRapidity));
      fastjet::ClusterSequenceArea clustSeq(fjParticles, jetDef, areaDef);                     // Attributes an area for each pseudojet in the list
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(clustSeq.inclusive_jets()); // No minimum pt before background subtraction
      if (jets.empty())
        return;
      // Perpendicular cone area subtraction, not the traditional subtraction (TODO: include an option for traditional area subtraction)
      auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], jetConfigurations.radiusJet); // This uses a geometric, pi*R^2 area, not exactly a ghost-based area!

      // Loop over clustered jets:
      int selectedJets = 0;

      fastjet::PseudoJet leadingJetSub;
      // bool hasLeadingJet = false; // Not needed: if the event has any jet, that is the leading jet. Check is superseded by the selectedJets information
      float leadingJetPt = -1.f;
      for (const auto& jet : jets) {
        // Jet must be fully contained in the acceptance (0.9 for ITS+TPC barrel)
        const float jet_eta = jet.eta();
        if (std::fabs(jet_eta) > (0.9f - jetConfigurations.radiusJet))
          continue;

        auto jetForSub = jet;
        // Subtracts same background estimated for highest pt jet, but every jet might have a slightly different area
        // (TODO: check possible problems with OO and physics impacts of this particular cone method and choice of single background estimator based on leading jet)
        // (TODO: improve for Pb-Pb, specially central!)
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
        // Jet pt must be larger than threshold:
        if (jetMinusBkg.pt() < jetConfigurations.minJetPt)
          continue;
        selectedJets++;

        // Store jet:
        tableJets(collIdx,
                  jetMinusBkg.pt(),
                  jetMinusBkg.eta(), // Using eta instead of rapidity
                  jetMinusBkg.phi(),
                  jetMinusBkg.constituents().size());

        // Finding the leading jet after subtraction (leading jet is NOT known a priori!):
        if (jetMinusBkg.pt() > leadingJetPt) {
          leadingJetPt = jetMinusBkg.pt();
          leadingJetSub = jetMinusBkg;
        }
      }
      histos.fill(HIST("hJetsPerEvent"), selectedJets);
      if (selectedJets == 0)
        return;
      histos.fill(HIST("hEventsWithJet"), 0.5);
      // Another version of this counter, which is already integrated in the Event Selection flow:
      if (doEventQA && !validJetAlreadyFound)
        fillEventSelectionQA(lastBinEvSel - 1, centrality); // hasRingJet passes
      validJetAlreadyFound = true;

      if (doJetKinematicsQA) {
        histos.fill(HIST("JetKinematicsQA/hLeadingJetPt"), leadingJetSub.pt());
        histos.fill(HIST("JetKinematicsQA/hLeadingJetEta"), leadingJetSub.eta());
        histos.fill(HIST("JetKinematicsQA/hLeadingJetPhi"), leadingJetSub.phi());

        // Now looping through jets again to calculate the correlations:
        for (const auto& jet : jets) {
          // Will recalculated background subtraction during QA to avoid storing jets in memory when running in non-QA cases:
          auto jetForSub = jet;
          fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);

          if (jetMinusBkg.pt() < jetConfigurations.minJetPt)
            continue;

          float cosTheta = cosThetaJets(leadingJetSub, jetMinusBkg);
          float deltaPhi = RecoDecay::constrainAngle(leadingJetSub.phi() - jetMinusBkg.phi(), -o2::constants::math::PI);
          float deltaEta = leadingJetSub.eta() - jetMinusBkg.eta();
          float deltaR = std::sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);

          histos.fill(HIST("JetKinematicsQA/hCosThetaToLeadingJet"), cosTheta);
          histos.fill(HIST("JetKinematicsQA/hDeltaPhiToLeadingJet"), deltaPhi);
          histos.fill(HIST("JetKinematicsQA/hDeltaEtaToLeadingJet"), deltaEta);
          histos.fill(HIST("JetKinematicsQA/hDeltaRToLeadingJet"), deltaR);

          // 2D correlations:
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsLeadJetPt"), selectedJets, leadingJetSub.pt());
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsJetPt"), selectedJets, jetMinusBkg.pt());
          histos.fill(HIST("JetKinematicsQA/h2dCosThetaToLeadvsDeltaPhiToLead"), cosTheta, deltaPhi);
          histos.fill(HIST("JetKinematicsQA/h2dCosThetaToLeadvsDeltaEtaToLead"), cosTheta, deltaEta);
          histos.fill(HIST("JetKinematicsQA/h2dCosThetaToLeadvsDeltaRToLead"), cosTheta, deltaR);
          histos.fill(HIST("JetKinematicsQA/h2dDeltaPhiToLeadvsDeltaEtaToLead"), deltaPhi, deltaEta);

          histos.fill(HIST("JetKinematicsQA/h2dJetPtvsDeltaPhiToLead"), jetMinusBkg.pt(), deltaPhi);    // Can't really get the energy of the jet, just the pt to make this comparison
          histos.fill(HIST("JetKinematicsQA/h2dJetEnergyvsDeltaPhiToLead"), jetMinusBkg.E(), deltaPhi); // Just a different scale
          histos.fill(HIST("JetKinematicsQA/h2dJetEnergyvsCosThetaToLead"), jetMinusBkg.E(), cosTheta);

          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsDeltaPhiToLead"), selectedJets, deltaPhi);
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsDeltaEtaToLead"), selectedJets, deltaEta);
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsCosThetaToLead"), selectedJets, cosTheta);
        }
        // Leading particle comparisons:
        histos.fill(HIST("JetVsLeadingParticleQA/hLeadingParticlePt"), leadingParticle.pt());
        histos.fill(HIST("JetVsLeadingParticleQA/hLeadingParticleEta"), leadingParticle.eta());
        histos.fill(HIST("JetVsLeadingParticleQA/hLeadingParticlePhi"), leadingParticle.phi());

        float deltaPhiParticleToJet = RecoDecay::constrainAngle(leadingJetSub.phi() - leadingParticle.phi(), -o2::constants::math::PI);
        float deltaEtaParticleToJet = leadingJetSub.eta() - leadingParticle.eta();
        float cosThetaParticleToJet = cosThetaJets(leadingJetSub, leadingParticle); // Takes advantage of the fact that this leading particle is a PseudoJet object

        histos.fill(HIST("JetVsLeadingParticleQA/hCosThetaLeadParticleToJet"), cosThetaParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/hDeltaPhiLeadParticleToJet"), deltaPhiParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/hDeltaEtaToLeadParticleToJet"), deltaEtaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dDeltaPhiParticleToLeadvsDeltaEtaParticleToLead"), deltaPhiParticleToJet, deltaEtaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsDeltaPhiParticleToLead"), selectedJets, deltaPhiParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsDeltaEtaParticleToLead"), selectedJets, deltaEtaParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsCosThetaParticleToLead"), selectedJets, cosThetaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsLeadParticlePt"), selectedJets, leadingParticle.pt());
        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadJetPtvsLeadParticlePt"), leadingJetSub.pt(), leadingParticle.pt());

        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadJetPtvsCosThetaParticleToLead"), leadingJetSub.pt(), cosThetaParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadParticlePtvsCosThetaParticleToLead"), leadingParticle.pt(), cosThetaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadJetPtvsDeltaPhiParticleToLead"), leadingJetSub.pt(), deltaPhiParticleToJet); // To see if there is any backgound in phi due to soft jets (or soft particles below)
        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadParticlePtvsDeltaPhiParticleToLead"), leadingParticle.pt(), deltaPhiParticleToJet);
      }
    } else { // Otherwise, simple jet clustering (TODO: this is the fall back for kConstituentBased while not implemented)
      fastjet::ClusterSequence clustSeq(fjParticles, jetDef);
      // Jet pt must be larger than threshold:
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(clustSeq.inclusive_jets(jetConfigurations.minJetPt));

      const int jetsInEvent = jets.size();
      histos.fill(HIST("hJetsPerEvent"), jetsInEvent); // Fills even in empty events, as this is a useful number to know!

      if (jetsInEvent == 0)
        return;
      histos.fill(HIST("hEventsWithJet"), 0.5);
      // Another version of this counter, which is already integrated in the Event Selection flow:
      if (doEventQA && !validJetAlreadyFound)
        fillEventSelectionQA(lastBinEvSel - 1, centrality); // hasRingJet passes
      validJetAlreadyFound = true;

      const auto& leadingJet = jets[0];
      for (const auto& jet : jets) {
        // Jet must be fully contained in the acceptance (0.9 for ITS+TPC barrel)
        const float jet_eta = jet.eta();
        if (std::fabs(jet_eta) > (0.9f - jetConfigurations.radiusJet))
          continue;

        tableJets(collIdx,
                  jet.pt(),
                  jet_eta, // Using eta instead of rapidity
                  jet.phi(),
                  jet.constituents().size());

        if (doJetKinematicsQA) {
          histos.fill(HIST("JetKinematicsQA/hJetPt"), jet.pt());
          histos.fill(HIST("JetKinematicsQA/hJetEta"), jet_eta);
          histos.fill(HIST("JetKinematicsQA/hJetPhi"), jet.phi());

          // Calculate angle to leading jet:
          float cosTheta = cosThetaJets(leadingJet, jet);

          // Calculate angular separation in projected angles:
          float deltaPhi = RecoDecay::constrainAngle(leadingJet.phi() - jet.phi(), -o2::constants::math::PI);
          float deltaEta = leadingJet.eta() - jet_eta;
          float deltaR = std::sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta); // 2D angular distance in the eta-phi plane

          histos.fill(HIST("JetKinematicsQA/hCosThetaToLeadingJet"), cosTheta); // Measuring the cosine, not angle, because it is faster!
          histos.fill(HIST("JetKinematicsQA/hDeltaPhiToLeadingJet"), deltaPhi);
          histos.fill(HIST("JetKinematicsQA/hDeltaEtaToLeadingJet"), deltaEta);
          histos.fill(HIST("JetKinematicsQA/hDeltaRToLeadingJet"), deltaR);

          // 2D correlations:
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsLeadJetPt"), jetsInEvent, leadingJet.pt());
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsJetPt"), jetsInEvent, jet.pt());
          histos.fill(HIST("JetKinematicsQA/h2dCosThetaToLeadvsDeltaPhiToLead"), cosTheta, deltaPhi);
          histos.fill(HIST("JetKinematicsQA/h2dCosThetaToLeadvsDeltaEtaToLead"), cosTheta, deltaEta);
          histos.fill(HIST("JetKinematicsQA/h2dCosThetaToLeadvsDeltaRToLead"), cosTheta, deltaR);
          histos.fill(HIST("JetKinematicsQA/h2dDeltaPhiToLeadvsDeltaEtaToLead"), deltaPhi, deltaEta);

          histos.fill(HIST("JetKinematicsQA/h2dJetPtvsDeltaPhiToLead"), jet.pt(), deltaPhi);    // Can't really get the energy of the jet, just the pt to make this comparison
          histos.fill(HIST("JetKinematicsQA/h2dJetEnergyvsDeltaPhiToLead"), jet.E(), deltaPhi); // Just a different scale
          histos.fill(HIST("JetKinematicsQA/h2dJetEnergyvsCosThetaToLead"), jet.E(), cosTheta);

          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsDeltaPhiToLead"), jetsInEvent, deltaPhi);
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsDeltaEtaToLead"), jetsInEvent, deltaEta);
          histos.fill(HIST("JetKinematicsQA/h2dJetsPerEventvsCosThetaToLead"), jetsInEvent, cosTheta);
        }
      }
      if (doJetKinematicsQA) {
        histos.fill(HIST("JetKinematicsQA/hLeadingJetPt"), leadingJet.pt());
        histos.fill(HIST("JetKinematicsQA/hLeadingJetEta"), leadingJet.eta());
        histos.fill(HIST("JetKinematicsQA/hLeadingJetPhi"), leadingJet.phi());

        // Leading particle comparisons:
        histos.fill(HIST("JetVsLeadingParticleQA/hLeadingParticlePt"), leadingParticle.pt());
        histos.fill(HIST("JetVsLeadingParticleQA/hLeadingParticleEta"), leadingParticle.eta());
        histos.fill(HIST("JetVsLeadingParticleQA/hLeadingParticlePhi"), leadingParticle.phi());

        double deltaPhiParticleToJet = RecoDecay::constrainAngle(leadingJet.phi() - leadingParticle.phi(), -o2::constants::math::PI);
        double deltaEtaParticleToJet = leadingJet.eta() - leadingParticle.eta();
        double cosThetaParticleToJet = cosThetaJets(leadingJet, leadingParticle); // Takes advantage of the fact that this leading particle is a PseudoJet object

        histos.fill(HIST("JetVsLeadingParticleQA/hCosThetaLeadParticleToJet"), cosThetaParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/hDeltaPhiLeadParticleToJet"), deltaPhiParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/hDeltaEtaToLeadParticleToJet"), deltaEtaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dDeltaPhiParticleToLeadvsDeltaEtaParticleToLead"), deltaPhiParticleToJet, deltaEtaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsDeltaPhiParticleToLead"), jetsInEvent, deltaPhiParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsDeltaEtaParticleToLead"), jetsInEvent, deltaEtaParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsCosThetaParticleToLead"), jetsInEvent, cosThetaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dJetsPerEventvsLeadParticlePt"), jetsInEvent, leadingParticle.pt());
        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadJetPtvsLeadParticlePt"), leadingJet.pt(), leadingParticle.pt());

        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadJetPtvsCosThetaParticleToLead"), leadingJet.pt(), cosThetaParticleToJet);
        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadParticlePtvsCosThetaParticleToLead"), leadingParticle.pt(), cosThetaParticleToJet);

        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadJetPtvsDeltaPhiParticleToLead"), leadingJet.pt(), deltaPhiParticleToJet); // To see if there is any backgound in phi due to soft jets (or soft particles below)
        histos.fill(HIST("JetVsLeadingParticleQA/h2dLeadParticlePtvsDeltaPhiParticleToLead"), leadingParticle.pt(), deltaPhiParticleToJet);
      }
    }
  }

  // Had to include DauTracks in subscription, even though I don't loop in it, for the indices
  // to resolve, avoiding " Exception while running: Index pointing to Tracks is not bound!"
  // Added the compiler option [[maybe_unused]] to avoid triggering any warnings because of this
  void processV0sData(SelCollisions::iterator const& collision, V0CandidatesWithTOF const& fullV0s, aod::BCsWithTimestamps const& bcs, [[maybe_unused]] DauTracks const& V0DauTracks)
  {
    float centrality = getCentrality(collision); // Strictly for QA. We save other types of centrality estimators in the derived data!

    // For event QA the last two indices never change for NEv_withJets and NEv_withV0s
    // (Not the best way to initialize this: runs once per collision! TODO: think of a better way to do it)
    int lastBinEvSel = histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->GetNbins();
    bool validV0AlreadyFound = false;

    histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    histos.fill(HIST("hEventSelectionVsCentrality"), 0. /* all collisions */, centrality);

    auto bc = bcs.iteratorAt(collision.bcId());
    if (!isEventAccepted(collision, bc, centrality, doEventQA))
      return; // Uses return instead of continue, as there is no explicit loop here

    if (doEventQA)
      fillCentralityProperties(collision, centrality);
    const uint64_t collIdx = collision.globalIndex();
    if (v0Selections.rejectTPCsectorBoundary)
      initCCDB(bc); // Substituted call from collision to bc for raw data

    // Fill event table:
    tableCollisions(collIdx,
                    collision.centFT0M(),
                    collision.centFT0C(),
                    collision.centFV0A()); // (TODO: add InteractionRate info and other useful cuts for later on in the analysis?)

    uint NLambdas = 0; // Counting particles per event
    uint NAntiLambdas = 0;
    uint NAmbiguous = 0;
    for (auto const& v0 : fullV0s) {
      V0SelCounter.resetForNewV0();
      V0SelCounter.fill(); // Fill for all v0 candidates
      if (doArmenterosQA)
        histos.fill(HIST("GeneralQA/h2dArmenterosAll"), v0.alpha(), v0.qtarm()); // fill AP plot for all V0s
      if (!passesGenericV0Cuts(v0))
        continue;

      if (doArmenterosQA)
        histos.fill(HIST("GeneralQA/h2dArmenterosKinematicSelected"), v0.alpha(), v0.qtarm());

      // Else, just continue the loop:
      bool isLambda = false;
      bool isAntiLambda = false;
      if (analyseLambda)
        isLambda = passesLambdaLambdaBarHypothesis(v0, collision, true);
      if (analyseAntiLambda)
        isAntiLambda = passesLambdaLambdaBarHypothesis(v0, collision, false);

      if (!isLambda && !isAntiLambda)
        continue; // Candidate is not considered to be a Lambda

      if (isLambda)
        NLambdas++;
      if (isAntiLambda)
        NAntiLambdas++;

      if (doArmenterosQA)
        histos.fill(HIST("GeneralQA/h2dArmenterosFullSelected"), v0.alpha(), v0.qtarm()); // cross-check
      if (isLambda && !isAntiLambda)
        histos.fill(HIST("GeneralQA/h2dArmenterosFullSelectedLambda"), v0.alpha(), v0.qtarm());
      if (!isLambda && isAntiLambda)
        histos.fill(HIST("GeneralQA/h2dArmenterosFullSelectedAntiLambda"), v0.alpha(), v0.qtarm());

      // int lambdaIdx = -1; // No need to pass armenteros
      if (isLambda && isAntiLambda) {
        NAmbiguous++;
        histos.fill(HIST("hAmbiguousLambdaCandidates"), 0);
        if (doArmenterosQA)
          histos.fill(HIST("GeneralQA/h2dArmenterosFullSelectedAmbiguous"), v0.alpha(), v0.qtarm()); // To know the discerning power of Armenteros in an Ambiguous Lambda vs AntiLambda case

        // Armenteros cut is not worth it! From QA histograms, only about 0.05% of ambiguous candidates are in the regions probable to be Lamda/AntiLambdas!
        // The statistics gain is not worth it.
        // // Third and final check to distinguish between Lambda and AntiLambda ambiguous v0s:
        // // (This check is only performed to recycle AMBIGUOUS candidates! Not a hard cut on all candidates!)
        // lambdaIdx = isCandidateArmenterosLambda(v0.alpha(), v0.qtarm());
      }
      // if (lambdaIdx == kIsArmenterosK0) continue; // Should just skip this step then!

      if (doEventQA)
        fillEventSelectionQA(lastBinEvSel, centrality); // hasRingV0 passes

      // // Extra competing mass rejection of Lambdas // (TODO: test competing mass cuts)
      // v0.mLambda()

      // Saving the Lambdas into a derived data column:
      auto const v0pt = v0.pt();
      const auto posTrackExtra = v0.template posTrack_as<DauTracks>();
      const auto negTrackExtra = v0.template negTrack_as<DauTracks>();
      tableV0s(collIdx,
               v0pt, v0.eta(), v0.phi(), // Using eta instead of rapidity
               isLambda, isAntiLambda,
               v0.mLambda(), v0.mAntiLambda(),
               v0.positivept(), v0.positiveeta(), v0.positivephi(),
               v0.negativept(), v0.negativeeta(), v0.negativephi(),
               posTrackExtra.tpcNSigmaPr(), posTrackExtra.tpcNSigmaPi(),
               negTrackExtra.tpcNSigmaPr(), negTrackExtra.tpcNSigmaPi(),
               v0.v0cosPA(), v0.v0radius(), v0.dcaV0daughters(), v0.dcapostopv(), v0.dcanegtopv());
      if (doEventQA && !validV0AlreadyFound)
        fillEventSelectionQA(lastBinEvSel, centrality); // hasRingV0 passes
      validV0AlreadyFound = true;

      if (doV0KinematicQA) {
        // Cache kinematics once
        const float v0y = v0.yLambda();
        const float v0phi = v0.phi();
        const float mLambda = v0.mLambda();
        const float mAntiLambda = v0.mAntiLambda();
        if (analyseLambda && isLambda) {
          // --- Basic kinematics ---
          histos.fill(HIST("V0KinematicQA/Lambda/hPt"), v0pt);
          histos.fill(HIST("V0KinematicQA/Lambda/hY"), v0y);
          histos.fill(HIST("V0KinematicQA/Lambda/hPhi"), v0phi);
          // --- Mass correlations ---
          histos.fill(HIST("V0KinematicQA/Lambda/hMassVsPt"), v0pt, mLambda);
          histos.fill(HIST("V0KinematicQA/Lambda/hMassVsY"), v0y, mLambda);
          histos.fill(HIST("V0KinematicQA/Lambda/hMassVsPhi"), v0phi, mLambda);
          // --- Kinematic correlations ---
          histos.fill(HIST("V0KinematicQA/Lambda/hYVsPt"), v0pt, v0y);
          histos.fill(HIST("V0KinematicQA/Lambda/hPhiVsPt"), v0pt, v0phi);
        }
        if (analyseAntiLambda && isAntiLambda) {
          // --- Basic kinematics ---
          histos.fill(HIST("V0KinematicQA/AntiLambda/hPt"), v0pt);
          histos.fill(HIST("V0KinematicQA/AntiLambda/hY"), v0y);
          histos.fill(HIST("V0KinematicQA/AntiLambda/hPhi"), v0phi);
          // --- Mass correlations ---
          histos.fill(HIST("V0KinematicQA/AntiLambda/hMassVsPt"), v0pt, mAntiLambda);
          histos.fill(HIST("V0KinematicQA/AntiLambda/hMassVsY"), v0y, mAntiLambda);
          histos.fill(HIST("V0KinematicQA/AntiLambda/hMassVsPhi"), v0phi, mAntiLambda);
          // --- Kinematic correlations ---
          histos.fill(HIST("V0KinematicQA/AntiLambda/hYVsPt"), v0pt, v0y);
          histos.fill(HIST("V0KinematicQA/AntiLambda/hPhiVsPt"), v0pt, v0phi);
        }
      }

      if (doCompleteTopoQA) {
        // Remaking these variables outside of the passesLambdaLambdaBarHypothesis. Loses performance, but that should be OK for QA
        histos.fill(HIST("V0KinematicQA/hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("V0KinematicQA/hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("V0KinematicQA/hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("V0KinematicQA/hPointingAngle"), std::acos(v0.v0cosPA()));
        histos.fill(HIST("V0KinematicQA/hV0Radius"), v0.v0radius());
        histos.fill(HIST("V0KinematicQA/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcNClsCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("V0KinematicQA/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcNClsCrossedRows(), negTrackExtra.itsNCls());
        histos.fill(HIST("V0KinematicQA/h2dPositivePtVsPhi"), v0.positivept(), computePhiMod(v0.positivephi(), 1));
        histos.fill(HIST("V0KinematicQA/h2dNegativePtVsPhi"), v0.negativept(), computePhiMod(v0.negativephi(), -1));
        if (isLambda && analyseLambda) {
          histos.fill(HIST("hMassLambda"), v0.mLambda());
          histos.fill(HIST("Lambda/h3dMassLambda"), centrality, v0pt, v0.mLambda());
          histos.fill(HIST("Lambda/hPosDCAToPV"), v0.dcapostopv());
          histos.fill(HIST("Lambda/hNegDCAToPV"), v0.dcanegtopv());
          histos.fill(HIST("Lambda/hDCADaughters"), v0.dcaV0daughters());
          histos.fill(HIST("Lambda/hPointingAngle"), std::acos(v0.v0cosPA()));
          histos.fill(HIST("Lambda/hV0Radius"), v0.v0radius());
          histos.fill(HIST("Lambda/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcNClsCrossedRows(), posTrackExtra.itsNCls());
          histos.fill(HIST("Lambda/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcNClsCrossedRows(), negTrackExtra.itsNCls());
          histos.fill(HIST("Lambda/h2dPositivePtVsPhi"), v0.positivept(), computePhiMod(v0.positivephi(), 1));
          histos.fill(HIST("Lambda/h2dNegativePtVsPhi"), v0.negativept(), computePhiMod(v0.negativephi(), -1));
          if (doTPCQA) {
            histos.fill(HIST("Lambda/h3dPosNsigmaTPC"), centrality, v0pt, posTrackExtra.tpcNSigmaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTPC"), centrality, v0pt, negTrackExtra.tpcNSigmaPi());
            histos.fill(HIST("Lambda/h3dPosTPCsignal"), centrality, v0pt, posTrackExtra.tpcSignal());
            histos.fill(HIST("Lambda/h3dNegTPCsignal"), centrality, v0pt, negTrackExtra.tpcSignal());
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
            histos.fill(HIST("Lambda/h3dPosNsigmaTOF"), centrality, v0pt, v0.tofNSigmaLaPr());
            histos.fill(HIST("Lambda/h3dNegNsigmaTOF"), centrality, v0pt, v0.tofNSigmaLaPi());
            histos.fill(HIST("Lambda/h3dPosTOFdeltaT"), centrality, v0pt, v0.posTOFDeltaTLaPr());
            histos.fill(HIST("Lambda/h3dNegTOFdeltaT"), centrality, v0pt, v0.negTOFDeltaTLaPi());
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
            histos.fill(HIST("Lambda/h5dV0PhiVsEta"), centrality, v0pt, v0.mLambda(), v0.phi(), v0.eta());
            histos.fill(HIST("Lambda/h5dPosPhiVsEta"), centrality, v0.positivept(), v0.mLambda(), v0.positivephi(), v0.positiveeta());
            histos.fill(HIST("Lambda/h5dNegPhiVsEta"), centrality, v0.negativept(), v0.mLambda(), v0.negativephi(), v0.negativeeta());
          }
        }
        if (isAntiLambda && analyseAntiLambda) {
          histos.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
          histos.fill(HIST("AntiLambda/h3dMassAntiLambda"), centrality, v0pt, v0.mAntiLambda());
          histos.fill(HIST("AntiLambda/hPosDCAToPV"), v0.dcapostopv());
          histos.fill(HIST("AntiLambda/hNegDCAToPV"), v0.dcanegtopv());
          histos.fill(HIST("AntiLambda/hDCADaughters"), v0.dcaV0daughters());
          histos.fill(HIST("AntiLambda/hPointingAngle"), std::acos(v0.v0cosPA()));
          histos.fill(HIST("AntiLambda/hV0Radius"), v0.v0radius());
          histos.fill(HIST("AntiLambda/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcNClsCrossedRows(), posTrackExtra.itsNCls());
          histos.fill(HIST("AntiLambda/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcNClsCrossedRows(), negTrackExtra.itsNCls());
          histos.fill(HIST("AntiLambda/h2dPositivePtVsPhi"), v0.positivept(), computePhiMod(v0.positivephi(), 1));
          histos.fill(HIST("AntiLambda/h2dNegativePtVsPhi"), v0.negativept(), computePhiMod(v0.negativephi(), -1));
          if (doTPCQA) {
            histos.fill(HIST("AntiLambda/h3dPosNsigmaTPC"), centrality, v0pt, posTrackExtra.tpcNSigmaPi());
            histos.fill(HIST("AntiLambda/h3dNegNsigmaTPC"), centrality, v0pt, negTrackExtra.tpcNSigmaPr());
            histos.fill(HIST("AntiLambda/h3dPosTPCsignal"), centrality, v0pt, posTrackExtra.tpcSignal());
            histos.fill(HIST("AntiLambda/h3dNegTPCsignal"), centrality, v0pt, negTrackExtra.tpcSignal());
            histos.fill(HIST("AntiLambda/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.pfracpos() * v0.p(), posTrackExtra.tpcNSigmaPi());
            histos.fill(HIST("AntiLambda/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.pfracneg() * v0.p(), negTrackExtra.tpcNSigmaPr());
            histos.fill(HIST("AntiLambda/h3dPosTPCsignalVsTrackPtot"), centrality, v0.pfracpos() * v0.p(), posTrackExtra.tpcSignal());
            histos.fill(HIST("AntiLambda/h3dNegTPCsignalVsTrackPtot"), centrality, v0.pfracneg() * v0.p(), negTrackExtra.tpcSignal());
            histos.fill(HIST("AntiLambda/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPi());
            histos.fill(HIST("AntiLambda/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPr());
            histos.fill(HIST("AntiLambda/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
            histos.fill(HIST("AntiLambda/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
          }
          if (doTOFQA) {
            histos.fill(HIST("AntiLambda/h3dPosNsigmaTOF"), centrality, v0pt, v0.tofNSigmaALaPi());
            histos.fill(HIST("AntiLambda/h3dNegNsigmaTOF"), centrality, v0pt, v0.tofNSigmaALaPr());
            histos.fill(HIST("AntiLambda/h3dPosTOFdeltaT"), centrality, v0pt, v0.posTOFDeltaTLaPi());
            histos.fill(HIST("AntiLambda/h3dNegTOFdeltaT"), centrality, v0pt, v0.negTOFDeltaTLaPr());
            histos.fill(HIST("AntiLambda/h3dPosNsigmaTOFvsTrackPtot"), centrality, v0.pfracpos() * v0.p(), v0.tofNSigmaALaPi());
            histos.fill(HIST("AntiLambda/h3dNegNsigmaTOFvsTrackPtot"), centrality, v0.pfracneg() * v0.p(), v0.tofNSigmaALaPr());
            histos.fill(HIST("AntiLambda/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.pfracpos() * v0.p(), v0.posTOFDeltaTLaPi());
            histos.fill(HIST("AntiLambda/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.pfracneg() * v0.p(), v0.negTOFDeltaTLaPr());
            histos.fill(HIST("AntiLambda/h3dPosNsigmaTOFvsTrackPt"), centrality, v0.positivept(), v0.tofNSigmaALaPi());
            histos.fill(HIST("AntiLambda/h3dNegNsigmaTOFvsTrackPt"), centrality, v0.negativept(), v0.tofNSigmaALaPr());
            histos.fill(HIST("AntiLambda/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTLaPi());
            histos.fill(HIST("AntiLambda/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTLaPr());
          }
          if (doEtaPhiQA) {
            histos.fill(HIST("AntiLambda/h5dV0PhiVsEta"), centrality, v0pt, v0.mAntiLambda(), v0.phi(), v0.eta());
            histos.fill(HIST("AntiLambda/h5dPosPhiVsEta"), centrality, v0.positivept(), v0.mAntiLambda(), v0.positivephi(), v0.positiveeta());
            histos.fill(HIST("AntiLambda/h5dNegPhiVsEta"), centrality, v0.negativept(), v0.mAntiLambda(), v0.negativephi(), v0.negativeeta());
          }
        }
      } // end CompleteTopoQA
    } // end V0s loop

    // Fill histograms on a per-event level:
    histos.fill(HIST("Lambda/hLambdasPerEvent"), NLambdas);
    histos.fill(HIST("AntiLambda/hAntiLambdasPerEvent"), NAntiLambdas);
    histos.fill(HIST("hAmbiguousPerEvent"), NAmbiguous);
    histos.fill(HIST("Lambda/h2dNbrOfLambdaVsCentrality"), centrality, NLambdas);
    histos.fill(HIST("AntiLambda/h2dNbrOfAntiLambdaVsCentrality"), centrality, NAntiLambdas);
  }

  PROCESS_SWITCH(lambdajetpolarizationions, processJetsData, "Process jets and produce derived data in Run 3 Data", true);
  PROCESS_SWITCH(lambdajetpolarizationions, processV0sData, "Process V0s and produce derived data in Run 3 Data", true);
  // PROCESS_SWITCH(lambdajetpolarizationions, processJetsMC, "Process jets and produced derived data in Run 3 MC", true);
  // PROCESS_SWITCH(lambdajetpolarizationions, processV0sMC, "Process V0s and produce derived data in Run 3 MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdajetpolarizationions>(cfgc)};
}
 