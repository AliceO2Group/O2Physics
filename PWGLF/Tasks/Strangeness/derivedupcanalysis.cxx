// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \file derivedupcanalysis.cxx
/// \brief Analysis of strangeness production in UPC collisions
/// \author Roman Nepeivoda (roman.nepeivoda@cern.ch)

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/strangenessMasks.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TFile.h>
#include <TH2F.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

using std::array;

using DauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using DauMCTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackMCIds, aod::DauTrackTPCPIDs>;

using V0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
using V0CandidatesMC = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0CoreMCLabels, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers>;

using CascadeCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFPIDs, aod::CascTOFNSigmas>;
using CascadeCandidatesMC = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascCoreMCLabels>;

using NeutronsMC = soa::Join<aod::ZDCNMCCollRefs, aod::ZDCNeutrons>;

using CascMCCoresFull = soa::Join<aod::CascMCCores, aod::CascMCCollRefs>;

using StraCollisonsFull = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>;
using StraCollisonFull = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps>::iterator;

using StraCollisonsFullMC = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels>;
using StraCollisonFullMC = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels>::iterator;

using StraMCCollisionsFull = soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>;
using V0MCCoresFull = soa::Join<aod::V0MCCores, aod::V0MCCollRefs>;

struct Derivedupcanalysis {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> analyseK0Short{"analyseK0Short", true, "process K0Short-like candidates"};
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", true, "process AntiLambda-like candidates"};
  Configurable<bool> analyseXi{"analyseXi", true, "process Xi-like candidates"};
  Configurable<bool> analyseAntiXi{"analyseAntiXi", true, "process AntiXi-like candidates"};
  Configurable<bool> analyseOmega{"analyseOmega", true, "process Omega-like candidates"};
  Configurable<bool> analyseAntiOmega{"analyseAntiOmega", true, "process AntiOmega-like candidates"};

  Configurable<std::vector<int>> generatorIds{"generatorIds", std::vector<int>{-1}, "MC generatorIds to process"};

  // Event selections
  struct : ConfigurableGroup {
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireIsTriggerTVX{"requireIsTriggerTVX", false, "require coincidence in FT0A and FT0C"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", false, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", true, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> studyUPConly{"studyUPConly", true, "is UPC-only analysis"};
    Configurable<bool> useUPCflag{"useUPCflag", false, "select UPC flagged events"};

    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", true, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerFV0Check{"cfgEvtRCTFlagCheckerFV0Check", true, "Evt sel: RCT flag checker FV0 check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } evSels;

  RCTFlagsChecker rctChecker;

  // Custom grouping
  std::vector<std::vector<int>> v0sGrouped;
  std::vector<std::vector<int>> cascadesGrouped;

  Configurable<bool> verbose{"verbose", false, "additional printouts"};

  // Acceptance selections
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
  Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

  Configurable<bool> doDaughterDCA{"doDaughterDCA", true, "dcaXY cut for daughter tracks"};

  // Standard V0 topological criteria
  struct : ConfigurableGroup {
    Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.5, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>std::fabs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};
    Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
  } v0cuts;
  static constexpr float kNCtauCutsV0[1][2] = {{6, 6.}};
  Configurable<LabeledArray<float>> nCtauCutV0{"nCtauCutV0", {kNCtauCutsV0[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "nCtauCutV0"};

  // Standard cascade topological criteria
  struct : ConfigurableGroup {
    Configurable<double> casccospa{"casccospa", 0.97, "Casc CosPA"};
    Configurable<float> dcacascdau{"dcacascdau", 1.2, "DCA Casc Daughters"};
    Configurable<float> cascradius{"cascradius", 0.6, "minimum cascade radius (cm)"};
    Configurable<float> cascradiusMax{"cascradiusMax", 1E5, "maximum cascade radius (cm)"};
    Configurable<float> bachbaryoncospa{"bachbaryoncospa", 2, "Bachelor baryon CosPA"};
    Configurable<float> bachbaryondcaxytopv{"bachbaryondcaxytopv", -1, "DCA bachelor baryon to PV"};
    Configurable<float> dcamesontopv{"dcamesontopv", 0.1, "DCA of meson doughter track To PV"};
    Configurable<float> dcabaryontopv{"dcabaryontopv", 0.05, "DCA of baryon doughter track To PV"};
    Configurable<float> dcabachtopv{"dcabachtopv", 0.04, "DCA Bach To PV"};
    Configurable<float> dcav0topv{"dcav0topv", 0.06, "DCA V0 To PV"};
    // Cascade specific selections
    Configurable<float> masswin{"masswin", 0.05, "mass window limit"};
    Configurable<float> lambdamasswin{"lambdamasswin", 0.005, "V0 Mass window limit"};
    Configurable<float> rejcomp{"rejcomp", 0.008, "competing Cascade rejection"};
  } casccuts;
  Configurable<bool> doBachelorBaryonCut{"doBachelorBaryonCut", false, "Enable Bachelor-Baryon cut "};
  static constexpr float kNCtauCutsCasc[1][2] = {{6., 6.}};
  Configurable<LabeledArray<float>> nCtauCutCasc{"nCtauCutCasc", {kNCtauCutsCasc[0], 2, {"lifetimecutXi", "lifetimecutOmega"}}, "nCtauCutCasc"};

  // UPC selections
  SGSelector sgSelector;
  struct : ConfigurableGroup {
    Configurable<float> fv0a{"fv0a", 50., "FV0A threshold"};
    Configurable<float> ft0a{"ft0a", 100., "FT0A threshold"};
    Configurable<float> ft0c{"ft0c", 50., "FT0C threshold"};
    Configurable<float> zdc{"zdc", 1., "ZDC threshold"};
    Configurable<int> genGapSide{"genGapSide", 0, "0 -- A, 1 -- C, 2 -- double"};
  } upcCuts;

  // Track quality
  struct : ConfigurableGroup {
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
    Configurable<bool> skipTPConly{"skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
    Configurable<bool> requireBachITSonly{"requireBachITSonly", false, "require that bachelor track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};
  } TrackConfigurations;

  // PID (TPC/TOF)
  struct : ConfigurableGroup {
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 1e+6, "tpcPidNsigmaCut"};
    Configurable<float> tofPidNsigmaCutLaPr{"tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
    Configurable<float> tofPidNsigmaCutLaPi{"tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
    Configurable<float> tofPidNsigmaCutK0Pi{"tofPidNsigmaCutK0Pi", 1e+6, "tofPidNsigmaCutK0Pi"};

    Configurable<float> tofPidNsigmaCutXiPi{"tofPidNsigmaCutXiPi", 1e+6, "tofPidNsigmaCutXiPi"};
    Configurable<float> tofPidNsigmaCutOmegaKaon{"tofPidNsigmaCutOmegaKaon", 1e+6, "tofPidNsigmaCutOmegaKaon"};

    Configurable<bool> doTPCQA{"doTPCQA", false, "do TPC QA histograms"};
    Configurable<bool> doTOFQA{"doTOFQA", false, "do TOF QA histograms"};

    // PID (TOF)
    Configurable<float> maxDeltaTimeProton{"maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimePion{"maxDeltaTimePion", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimeKaon{"maxDeltaTimeKaon", 1e+9, "check maximum allowed time"};
  } PIDConfigurations;

  Configurable<bool> doKienmaticQA{"doKienmaticQA", true, "do Kinematic QA histograms"};
  Configurable<int> doDetectPropQA{"doDetectPropQA", 0, "do Detector/ITS map QA: 0: no, 1: 4D, 2: 5D with mass"};
  Configurable<bool> doPlainTopoQA{"doPlainTopoQA", true, "do simple 1D QA of candidates"};

  struct : ConfigurableGroup {
    ConfigurableAxis axisFT0Aampl{"axisFT0Aampl", {100, 0.0f, 2000.0f}, "FT0Aamplitude"};
    ConfigurableAxis axisFT0Campl{"axisFT0Campl", {100, 0.0f, 2000.0f}, "FT0Camplitude"};
    ConfigurableAxis axisFT0ampl{"axisFT0ampl", {2002, -1.5f, 2000.5f}, "axisFT0ampl"};
    ConfigurableAxis axisFV0Aampl{"axisFV0Aampl", {100, 0.0f, 2000.0f}, "FV0Aamplitude"};
    ConfigurableAxis axisFDDAampl{"axisFDDAampl", {100, 0.0f, 2000.0f}, "FDDAamplitude"};
    ConfigurableAxis axisFDDCampl{"axisFDDCampl", {100, 0.0f, 2000.0f}, "FDDCamplitude"};
    ConfigurableAxis axisZNAampl{"axisZNAampl", {100, 0.0f, 250.0f}, "ZNAamplitude"};
    ConfigurableAxis axisZNCampl{"axisZNCampl", {100, 0.0f, 250.0f}, "ZNCamplitude"};
  } axisDetectors;

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};
  Configurable<bool> doTreatPiToMuon{"doTreatPiToMuon", false, "Take pi decay into muon into account in MC"};
  Configurable<bool> calculateFeeddownMatrix{"calculateFeeddownMatrix", true, "fill feeddown matrix if MC"};
  ConfigurableAxis axisGeneratorIds{"axisGeneratorIds", {256, -0.5f, 255.5f}, "axis for generatorIds"};
  Configurable<bool> checkNeutronsInMC{"checkNeutronsInMC", true, "require no neutrons for single-gap in MC"};
  Configurable<float> neutronEtaCut{"neutronEtaCut", 8.8, "ZN acceptance"};

  // Occupancy cut
  Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
  Configurable<float> maxOccupancy{"maxOccupancy", 1000, "maximum occupancy from neighbouring collisions"};

  // z vertex cut
  Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10.0f, "max Z vtx position"};

  // Kinematic axes
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for v0 analysis"};
  ConfigurableAxis axisPtXi{"axisPtXi", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for cascade analysis"};
  ConfigurableAxis axisPtCoarse{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};
  ConfigurableAxis axisEta{"axisEta", {100, -2.0f, 2.0f}, "#eta"};
  ConfigurableAxis axisRap{"axisRap", {100, -2.0f, 2.0f}, "y"};

  // Invariant mass axes
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.28f, 1.36f}, ""};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.59f, 1.75f}, ""};
  std::vector<ConfigurableAxis> axisInvMass = {axisK0Mass,
                                               axisLambdaMass,
                                               axisLambdaMass,
                                               axisXiMass,
                                               axisXiMass,
                                               axisOmegaMass,
                                               axisOmegaMass};

  ConfigurableAxis axisNTracksGlobal{"axisNTracksGlobal", {101, -1.5f, 99.5f}, "Number of global tracks"};
  ConfigurableAxis axisNTracksPVeta1{"axisNTracksPVeta1", {100, -0.5f, 99.5f}, "Number of PV contributors in |eta| < 1"};
  ConfigurableAxis axisNTracksPVeta05{"axisNTracksPVeta05", {100, -0.5f, 99.5f}, "Number of PV contributors in |eta| < 0.5"};
  ConfigurableAxis axisNAssocColl{"axisNAssocColl", {10, -0.5f, 9.5f}, "Number of assoc. rec. collisions"};
  ConfigurableAxis axisNTracksTotalExceptITSonly{"axisNTracksTotalExceptITSonly", {100, -0.5f, 99.5f}, "Number of ITS-TPC and TPC only tracks"};
  ConfigurableAxis axisNchInvMass{"axisNchInvMass", {201, -1.5f, 199.5f}, "Number of charged particles for kTHnSparseF"};

  ConfigurableAxis axisFT0Cqa{"axisFT0Cqa",
                              {VARIABLE_WIDTH, -1.5, -0.5, 0., 1., 5, 10, 20, 30, 40, 50, 60, 70, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 105.5},
                              "FT0C (%)"};

  ConfigurableAxis axisFT0C{"axisFT0C",
                            {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 105.5},
                            "FT0C (%)"};

  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};

  // UPC axes
  ConfigurableAxis axisSelGap{"axisSelGap", {4, -1.5, 2.5}, "Gap side"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // MC coll assoc QA axis
  ConfigurableAxis axisMonteCarloNch{"axisMonteCarloNch", {300, 0.0f, 3000.0f}, "N_{ch} MC"};

  // Track quality axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, -0.5f, 159.5f}, "N TPC rows"};
  ConfigurableAxis axisITSclus{"axisITSclus", {7, -0.5f, 6.5f}, "N ITS Clusters"};
  ConfigurableAxis axisITScluMap{"axisITScluMap", {128, -0.5f, 127.5f}, "ITS Cluster map"};
  ConfigurableAxis axisDetMap{"axisDetMap", {16, -0.5f, 15.5f}, "Detector use map"};
  ConfigurableAxis axisITScluMapCoarse{"axisITScluMapCoarse", {16, -3.5f, 12.5f}, "ITS Coarse cluster map"};
  ConfigurableAxis axisDetMapCoarse{"axisDetMapCoarse", {5, -0.5f, 4.5f}, "Detector Coarse user map"};

  // Topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {80, -4.0f, 4.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {48, 0.0f, 1.2f}, "DCA (cm)"};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {100, 0.0f, 0.5f}, "pointing angle (rad)"};
  ConfigurableAxis axisCosPA{"axisCosPA", {300, 0.97f, 1.0f}, "cosPA"};
  ConfigurableAxis axisV0Radius{"axisV0Radius", {100, 0.0f, 10.0f}, "V0 2D radius (cm)"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {200, 0.0f, 200.0f}, "TPC signal"};
  ConfigurableAxis axisTOFdeltaT{"axisTOFdeltaT", {200, -5000.0f, 5000.0f}, "TOF Delta T (ps)"};
  ConfigurableAxis axisCtau{"axisCtau", {200, 0.0f, 20.0f}, "c x tau (cm)"};

  static constexpr std::string_view kParticlenames[] = {"K0Short", "Lambda", "AntiLambda", "Xi", "AntiXi", "Omega", "AntiOmega"};

  void setBits(std::bitset<kSelNum>& mask, std::initializer_list<int> selections)
  {
    for (const int& sel : selections) {
      mask.set(sel);
    }
  }

  template <int partID>
  void addTopoHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    if (isCascade) {
      histos.add(Form("%s/hCascCosPA", kParticlenames[partID].data()), "hCascCosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
      histos.add(Form("%s/hDCACascDaughters", kParticlenames[partID].data()), "hDCACascDaughters", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
      histos.add(Form("%s/hCascRadius", kParticlenames[partID].data()), "hCascRadius", kTH2D, {axisPtCoarse, {500, 0.0f, 50.0f}});
      histos.add(Form("%s/hMesonDCAToPV", kParticlenames[partID].data()), "hMesonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
      histos.add(Form("%s/hBaryonDCAToPV", kParticlenames[partID].data()), "hBaryonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
      histos.add(Form("%s/hBachDCAToPV", kParticlenames[partID].data()), "hBachDCAToPV", kTH2F, {axisPtCoarse, {200, -1.0f, 1.0f}});
      histos.add(Form("%s/hV0CosPA", kParticlenames[partID].data()), "hV0CosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
      histos.add(Form("%s/hV0Radius", kParticlenames[partID].data()), "hV0Radius", kTH2D, {axisPtCoarse, axisV0Radius});
      histos.add(Form("%s/hDCAV0Daughters", kParticlenames[partID].data()), "hDCAV0Daughters", kTH2F, {axisPtCoarse, axisDCAdau});
      histos.add(Form("%s/hDCAV0ToPV", kParticlenames[partID].data()), "hDCAV0ToPV", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
      histos.add(Form("%s/hMassLambdaDau", kParticlenames[partID].data()), "hMassLambdaDau", kTH2F, {axisPtCoarse, axisLambdaMass});
      histos.add(Form("%s/hCtau", kParticlenames[partID].data()), "hCtau", kTH2F, {axisPtCoarse, axisCtau});
      if (doBachelorBaryonCut) {
        histos.add(Form("%s/hBachBaryonCosPA", kParticlenames[partID].data()), "hBachBaryonCosPA", kTH2F, {axisPtCoarse, {100, 0.0f, 1.0f}});
        histos.add(Form("%s/hBachBaryonDCAxyToPV", kParticlenames[partID].data()), "hBachBaryonDCAxyToPV", kTH2F, {axisPtCoarse, {300, -3.0f, 3.0f}});
      }
    } else {
      histos.add(Form("%s/hPosDCAToPV", kParticlenames[partID].data()), "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add(Form("%s/hNegDCAToPV", kParticlenames[partID].data()), "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add(Form("%s/hDCADaughters", kParticlenames[partID].data()), "hDCADaughters", kTH1F, {axisDCAdau});
      histos.add(Form("%s/hPointingAngle", kParticlenames[partID].data()), "hPointingAngle", kTH1F, {axisPointingAngle});
      histos.add(Form("%s/hCosPA", kParticlenames[partID].data()), "hCosPA", kTH1F, {axisCosPA});
      histos.add(Form("%s/hV0Radius", kParticlenames[partID].data()), "hV0Radius", kTH1F, {axisV0Radius});
      histos.add(Form("%s/h2dPositiveITSvsTPCpts", kParticlenames[partID].data()), "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add(Form("%s/h2dNegativeITSvsTPCpts", kParticlenames[partID].data()), "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add(Form("%s/hCtau", kParticlenames[partID].data()), "hCtau", kTH2F, {axisPtCoarse, axisCtau});
    }
  }

  template <int partID>
  void addTPCQAHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    histos.add(Form("%s/h3dPosNsigmaTPC", kParticlenames[partID].data()), "h3dPosNsigmaTPC", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dNegNsigmaTPC", kParticlenames[partID].data()), "h3dNegNsigmaTPC", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dPosTPCsignal", kParticlenames[partID].data()), "h3dPosTPCsignal", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});
    histos.add(Form("%s/h3dNegTPCsignal", kParticlenames[partID].data()), "h3dNegTPCsignal", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});

    histos.add(Form("%s/h3dPosNsigmaTPCvsTrackPtot", kParticlenames[partID].data()), "h3dPosNsigmaTPCvsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dNegNsigmaTPCvsTrackPtot", kParticlenames[partID].data()), "h3dNegNsigmaTPCvsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});

    histos.add(Form("%s/h3dPosTPCsignalVsTrackPtot", kParticlenames[partID].data()), "h3dPosTPCsignalVsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});
    histos.add(Form("%s/h3dNegTPCsignalVsTrackPtot", kParticlenames[partID].data()), "h3dNegTPCsignalVsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});

    histos.add(Form("%s/h3dPosNsigmaTPCvsTrackPt", kParticlenames[partID].data()), "h3dPosNsigmaTPCvsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dNegNsigmaTPCvsTrackPt", kParticlenames[partID].data()), "h3dNegNsigmaTPCvsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});

    histos.add(Form("%s/h3dPosTPCsignalVsTrackPt", kParticlenames[partID].data()), "h3dPosTPCsignalVsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});
    histos.add(Form("%s/h3dNegTPCsignalVsTrackPt", kParticlenames[partID].data()), "h3dNegTPCsignalVsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});

    if (isCascade) {
      histos.add(Form("%s/h3dBachTPCsignal", kParticlenames[partID].data()), "h3dBachTPCsignal", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});
      histos.add(Form("%s/h3dBachNsigmaTPC", kParticlenames[partID].data()), "h3dBachNsigmaTPC", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});
      histos.add(Form("%s/h3dBachNsigmaTPCvsTrackPtot", kParticlenames[partID].data()), "h3dBachNsigmaTPCvsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});
      histos.add(Form("%s/h3dBachTPCsignalVsTrackPtot", kParticlenames[partID].data()), "h3dBachTPCsignalVsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});
      histos.add(Form("%s/h3dBachNsigmaTPCvsTrackPt", kParticlenames[partID].data()), "h3dBachNsigmaTPCvsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisNsigmaTPC});
      histos.add(Form("%s/h3dBachTPCsignalVsTrackPt", kParticlenames[partID].data()), "h3dBachTPCsignalVsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTPCsignal});
    }
  }

  template <int partID>
  void addTOFQAHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    histos.add(Form("%s/h3dPosTOFdeltaT", kParticlenames[partID].data()), "h3dPosTOFdeltaT", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dNegTOFdeltaT", kParticlenames[partID].data()), "h3dNegTOFdeltaT", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dPosTOFdeltaTvsTrackPtot", kParticlenames[partID].data()), "h3dPosTOFdeltaTvsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dNegTOFdeltaTvsTrackPtot", kParticlenames[partID].data()), "h3dNegTOFdeltaTvsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dPosTOFdeltaTvsTrackPt", kParticlenames[partID].data()), "h3dPosTOFdeltaTvsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dNegTOFdeltaTvsTrackPt", kParticlenames[partID].data()), "h3dNegTOFdeltaTvsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
    if (isCascade) {
      histos.add(Form("%s/h3dBachTOFdeltaT", kParticlenames[partID].data()), "h3dBachTOFdeltaT", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
      histos.add(Form("%s/h3dBachTOFdeltaTvsTrackPtot", kParticlenames[partID].data()), "h3dBachTOFdeltaTvsTrackPtot", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
      histos.add(Form("%s/h3dBachTOFdeltaTvsTrackPt", kParticlenames[partID].data()), "h3dBachTOFdeltaTvsTrackPt", kTH3F, {axisDetectors.axisFT0ampl, axisPtCoarse, axisTOFdeltaT});
    }
  }

  template <int partID>
  void addKinematicQAHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    histos.add(Form("%s/h3dPosEtaPt", kParticlenames[partID].data()), "h3dPosEtaPt", kTH3F, {axisPtCoarse, axisEta, axisSelGap});
    histos.add(Form("%s/h3dNegEtaPt", kParticlenames[partID].data()), "h3dNegEtaPt", kTH3F, {axisPtCoarse, axisEta, axisSelGap});
    histos.add(Form("%s/h3dRapPt", kParticlenames[partID].data()), "h3dRapPt", kTH3F, {axisPtCoarse, axisRap, axisSelGap});
    if (isCascade) {
      histos.add(Form("%s/h3dBachEtaPt", kParticlenames[partID].data()), "h3dBachEtaPt", kTH3F, {axisPtCoarse, axisEta, axisSelGap});
    }
  }

  template <int partID>
  void addDetectorPropHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    if (doDetectPropQA == 1) {
      if (isCascade) {
        histos.add(Form("%s/h8dDetectPropVsCentrality", kParticlenames[partID].data()), "h8dDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
      } else {
        histos.add(Form("%s/h6dDetectPropVsCentrality", kParticlenames[partID].data()), "h6dDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
      }
      histos.add(Form("%s/h4dPosDetectPropVsCentrality", kParticlenames[partID].data()), "h4dPosDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMap, axisITScluMap, axisPtCoarse});
      histos.add(Form("%s/h4dNegDetectPropVsCentrality", kParticlenames[partID].data()), "h4dNegDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMap, axisITScluMap, axisPtCoarse});
      histos.add(Form("%s/h4dBachDetectPropVsCentrality", kParticlenames[partID].data()), "h4dBachDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMap, axisITScluMap, axisPtCoarse});
    }
    if (doDetectPropQA == 2) {
      if (isCascade) {
        histos.add(Form("%s/h9dDetectPropVsCentrality", kParticlenames[partID].data()), "h9dDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisInvMass.at(partID)});
      } else {
        histos.add(Form("%s/h7dDetectPropVsCentrality", kParticlenames[partID].data()), "h7dDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisInvMass.at(partID)});
      }
      histos.add(Form("%s/h5dPosDetectPropVsCentrality", kParticlenames[partID].data()), "h5dPosDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMap, axisITScluMap, axisPtCoarse, axisInvMass.at(partID)});
      histos.add(Form("%s/h5dNegDetectPropVsCentrality", kParticlenames[partID].data()), "h5dNegDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMap, axisITScluMap, axisPtCoarse, axisInvMass.at(partID)});
      histos.add(Form("%s/h5dBachDetectPropVsCentrality", kParticlenames[partID].data()), "h5dBachDetectPropVsCentrality", kTHnSparseF, {axisDetectors.axisFT0ampl, axisDetMap, axisITScluMap, axisPtCoarse, axisInvMass.at(partID)});
    }
  }

  template <int partID>
  void addHistograms(HistogramRegistry& histos)
  {
    histos.add(Form("%s/h7dMass", kParticlenames[partID].data()), "h7dMass", kTHnSparseF, {axisDetectors.axisFT0ampl, axisPt, axisInvMass.at(partID), axisSelGap, axisNchInvMass, axisRap, axisEta});
    histos.add(Form("%s/h2dMass", kParticlenames[partID].data()), "h2dMass", kTH2F, {axisInvMass.at(partID), axisSelGap});
    if (doPlainTopoQA) {
      addTopoHistograms<partID>(histos);
    }
    if (PIDConfigurations.doTPCQA) {
      addTPCQAHistograms<partID>(histos);
    }
    if (PIDConfigurations.doTOFQA) {
      addTOFQAHistograms<partID>(histos);
    }
    if (doKienmaticQA) {
      addKinematicQAHistograms<partID>(histos);
    }
    addDetectorPropHistograms<partID>(histos);
  }

  template <int partID, typename TCand, typename TCollision>
  void fillHistogramsV0(TCand cand, TCollision coll, int gap)
  {
    float invMass = 0;
    float ft0ampl = -1.f;
    if (gap == 0) {
      ft0ampl = coll.totalFT0AmplitudeC();
    } else if (gap == 1) {
      ft0ampl = coll.totalFT0AmplitudeA();
    }
    float pT = cand.pt();
    float rapidity = 1e6;

    // c x tau
    float ctau = 0;

    float tpcNsigmaPos = 0;
    float tpcNsigmaNeg = 0;
    float tofDeltaTPos = 0;
    float tofDeltaTNeg = 0;

    auto posTrackExtra = cand.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = cand.template negTrackExtra_as<DauTracks>();

    bool posIsFromAfterburner = posTrackExtra.itsChi2PerNcl() < 0;
    bool negIsFromAfterburner = negTrackExtra.itsChi2PerNcl() < 0;

    uint posDetMap = computeDetBitmap(posTrackExtra.detectorMap());
    int posITSclusMap = computeITSclusBitmap(posTrackExtra.itsClusterMap(), posIsFromAfterburner);
    uint negDetMap = computeDetBitmap(negTrackExtra.detectorMap());
    int negITSclusMap = computeITSclusBitmap(negTrackExtra.itsClusterMap(), negIsFromAfterburner);

    if (partID == 0) {
      histos.fill(HIST("generalQA/h2dArmenterosSelected"), cand.alpha(), cand.qtarm());
      invMass = cand.mK0Short();
      rapidity = cand.yK0Short();
      ctau = cand.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * o2::constants::physics::MassK0Short;
      if (PIDConfigurations.doTOFQA) {
        tofDeltaTPos = cand.posTOFDeltaTK0Pi();
        tofDeltaTNeg = cand.negTOFDeltaTK0Pi();
      }
      if (PIDConfigurations.doTPCQA) {
        tpcNsigmaPos = posTrackExtra.tpcNSigmaPi();
        tpcNsigmaNeg = negTrackExtra.tpcNSigmaPi();
      }
    } else if (partID == 1) {
      invMass = cand.mLambda();
      rapidity = cand.yLambda();
      ctau = cand.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * o2::constants::physics::MassLambda0;
      if (PIDConfigurations.doTOFQA) {
        tofDeltaTPos = cand.posTOFDeltaTLaPr();
        tofDeltaTNeg = cand.negTOFDeltaTLaPi();
      }
      if (PIDConfigurations.doTPCQA) {
        tpcNsigmaPos = posTrackExtra.tpcNSigmaPr();
        tpcNsigmaNeg = negTrackExtra.tpcNSigmaPi();
      }
    } else if (partID == 2) {
      invMass = cand.mAntiLambda();
      rapidity = cand.yLambda();
      ctau = cand.distovertotmom(coll.posX(), coll.posY(), coll.posZ()) * o2::constants::physics::MassLambda0Bar;
      if (PIDConfigurations.doTOFQA) {
        tofDeltaTPos = cand.posTOFDeltaTLaPi();
        tofDeltaTNeg = cand.negTOFDeltaTLaPr();
      }
      if (PIDConfigurations.doTPCQA) {
        tpcNsigmaPos = posTrackExtra.tpcNSigmaPi();
        tpcNsigmaNeg = negTrackExtra.tpcNSigmaPr();
      }
    } else {
      LOG(fatal) << "Particle is unknown!";
    }

    histos.fill(HIST(kParticlenames[partID]) + HIST("/h2dMass"), invMass, gap);
    histos.fill(HIST(kParticlenames[partID]) + HIST("/h7dMass"), ft0ampl, pT, invMass, gap, coll.multNTracksGlobal(), rapidity, cand.eta());
    if (doKienmaticQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosEtaPt"), pT, cand.positiveeta(), gap);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegEtaPt"), pT, cand.negativeeta(), gap);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dRapPt"), pT, rapidity, gap);
    }
    if (doPlainTopoQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hPosDCAToPV"), cand.dcapostopv());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hNegDCAToPV"), cand.dcanegtopv());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hDCADaughters"), cand.dcaV0daughters());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hPointingAngle"), std::acos(cand.v0cosPA()));
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hCosPA"), cand.v0cosPA());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hV0Radius"), cand.v0radius());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hCtau"), pT, ctau);
    }
    if (doDetectPropQA == 1) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h6dDetectPropVsCentrality"), ft0ampl, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pT);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h4dPosDetectPropVsCentrality"), ft0ampl, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h4dNegDetectPropVsCentrality"), ft0ampl, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT);
    }
    if (doDetectPropQA == 2) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h7dPosDetectPropVsCentrality"), ft0ampl, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pT, invMass);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h5dPosDetectPropVsCentrality"), ft0ampl, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT, invMass);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h5dNegDetectPropVsCentrality"), ft0ampl, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT, invMass);
    }
    if (PIDConfigurations.doTPCQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTPCsignal"), ft0ampl, pT, posTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTPCsignal"), ft0ampl, pT, negTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPtot"), ft0ampl, cand.positivept() * std::cosh(cand.positiveeta()), posTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPtot"), ft0ampl, cand.negativept() * std::cosh(cand.negativeeta()), negTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPt"), ft0ampl, cand.positivept(), posTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPt"), ft0ampl, cand.negativept(), negTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPt"), ft0ampl, cand.positivept(), posTrackExtra.tpcNSigmaPi());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPt"), ft0ampl, cand.negativept(), negTrackExtra.tpcNSigmaPi());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPtot"), ft0ampl, cand.positivept() * std::cosh(cand.positiveeta()), tpcNsigmaPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPtot"), ft0ampl, cand.negativept() * std::cosh(cand.negativeeta()), tpcNsigmaNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosNsigmaTPC"), ft0ampl, pT, tpcNsigmaPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegNsigmaTPC"), ft0ampl, pT, tpcNsigmaNeg);
    }
    if (PIDConfigurations.doTOFQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPt"), ft0ampl, cand.positivept(), tofDeltaTPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPt"), ft0ampl, cand.negativept(), tofDeltaTNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTOFdeltaT"), ft0ampl, pT, tofDeltaTPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTOFdeltaT"), ft0ampl, pT, tofDeltaTNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPtot"), ft0ampl, cand.positivept() * std::cosh(cand.positiveeta()), tofDeltaTPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPtot"), ft0ampl, cand.negativept() * std::cosh(cand.negativeeta()), tofDeltaTNeg);
    }
  }

  template <int partID, typename TCand, typename TCollision>
  void fillHistogramsCasc(TCand cand, TCollision coll, const int gap)
  {
    float invMass = 0;
    float centrality = -1.f;
    if (gap == 0) {
      centrality = coll.totalFT0AmplitudeC();
    } else if (gap == 1) {
      centrality = coll.totalFT0AmplitudeA();
    }
    float pT = cand.pt();
    float rapidity = 1e6;

    // Access daughter tracks
    auto posTrackExtra = cand.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = cand.template negTrackExtra_as<DauTracks>();
    auto bachTrackExtra = cand.template bachTrackExtra_as<DauTracks>();

    bool posIsFromAfterburner = posTrackExtra.itsChi2PerNcl() < 0;
    bool negIsFromAfterburner = negTrackExtra.itsChi2PerNcl() < 0;
    bool bachIsFromAfterburner = bachTrackExtra.itsChi2PerNcl() < 0;

    uint posDetMap = computeDetBitmap(posTrackExtra.detectorMap());
    int posITSclusMap = computeITSclusBitmap(posTrackExtra.itsClusterMap(), posIsFromAfterburner);
    uint negDetMap = computeDetBitmap(negTrackExtra.detectorMap());
    int negITSclusMap = computeITSclusBitmap(negTrackExtra.itsClusterMap(), negIsFromAfterburner);
    uint bachDetMap = computeDetBitmap(bachTrackExtra.detectorMap());
    int bachITSclusMap = computeITSclusBitmap(bachTrackExtra.itsClusterMap(), bachIsFromAfterburner);

    // c x tau
    float decayPos = std::hypot(cand.x() - coll.posX(), cand.y() - coll.posY(), cand.z() - coll.posZ());
    float totalMom = std::hypot(cand.px(), cand.py(), cand.pz());

    float ctau = 0;

    float tpcNsigmaPos = 0;
    float tpcNsigmaNeg = 0;
    float tpcNsigmaBach = 0;
    float tofDeltaTPos = 0;
    float tofDeltaTNeg = 0;
    float tofDeltaTBach = 0;

    if (partID == 3) {
      invMass = cand.mXi();
      ctau = totalMom != 0 ? o2::constants::physics::MassXiMinus * decayPos / (totalMom * ctauxiPDG) : 1e6;
      rapidity = cand.yXi();

      if (PIDConfigurations.doTPCQA) {
        tpcNsigmaPos = posTrackExtra.tpcNSigmaPr();
        tpcNsigmaNeg = negTrackExtra.tpcNSigmaPi();
        tpcNsigmaBach = bachTrackExtra.tpcNSigmaPi();
      }
      if (PIDConfigurations.doTOFQA) {
        tofDeltaTPos = cand.posTOFDeltaTXiPr();
        tofDeltaTNeg = cand.negTOFDeltaTXiPi();
        tofDeltaTBach = cand.bachTOFDeltaTXiPi();
      }
    } else if (partID == 4) {
      invMass = cand.mXi();
      ctau = totalMom != 0 ? o2::constants::physics::MassXiPlusBar * decayPos / (totalMom * ctauxiPDG) : 1e6;
      rapidity = cand.yXi();

      if (PIDConfigurations.doTPCQA) {
        tpcNsigmaPos = posTrackExtra.tpcNSigmaPi();
        tpcNsigmaNeg = negTrackExtra.tpcNSigmaPr();
        tpcNsigmaBach = bachTrackExtra.tpcNSigmaPi();
      }
      if (PIDConfigurations.doTOFQA) {
        tofDeltaTPos = cand.posTOFDeltaTXiPi();
        tofDeltaTNeg = cand.negTOFDeltaTXiPr();
        tofDeltaTBach = cand.bachTOFDeltaTXiPi();
      }

    } else if (partID == 5) {
      invMass = cand.mOmega();
      ctau = totalMom != 0 ? o2::constants::physics::MassOmegaMinus * decayPos / (totalMom * ctauomegaPDG) : 1e6;
      rapidity = cand.yOmega();

      if (PIDConfigurations.doTPCQA) {
        tpcNsigmaPos = posTrackExtra.tpcNSigmaPr();
        tpcNsigmaNeg = negTrackExtra.tpcNSigmaPi();
        tpcNsigmaBach = bachTrackExtra.tpcNSigmaKa();
      }
      if (PIDConfigurations.doTOFQA) {
        tofDeltaTPos = cand.posTOFDeltaTOmPi();
        tofDeltaTNeg = cand.posTOFDeltaTOmPr();
        tofDeltaTBach = cand.bachTOFDeltaTOmKa();
      }

    } else if (partID == 6) {
      invMass = cand.mOmega();
      ctau = totalMom != 0 ? o2::constants::physics::MassOmegaPlusBar * decayPos / (totalMom * ctauomegaPDG) : 1e6;
      rapidity = cand.yOmega();

      if (PIDConfigurations.doTPCQA) {
        tpcNsigmaPos = posTrackExtra.tpcNSigmaPi();
        tpcNsigmaNeg = negTrackExtra.tpcNSigmaPr();
        tpcNsigmaBach = bachTrackExtra.tpcNSigmaKa();
      }
      if (PIDConfigurations.doTOFQA) {
        tofDeltaTPos = cand.posTOFDeltaTOmPr();
        tofDeltaTNeg = cand.posTOFDeltaTOmPi();
        tofDeltaTBach = cand.bachTOFDeltaTOmKa();
      }
    }
    histos.fill(HIST(kParticlenames[partID]) + HIST("/h2dMass"), invMass, gap);
    histos.fill(HIST(kParticlenames[partID]) + HIST("/h7dMass"), centrality, pT, invMass, gap, coll.multNTracksGlobal(), rapidity, cand.eta());
    if (doKienmaticQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosEtaPt"), pT, cand.positiveeta(), gap);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegEtaPt"), pT, cand.negativeeta(), gap);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachEtaPt"), pT, cand.bacheloreta(), gap);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dRapPt"), pT, rapidity, gap);
    }
    if (doPlainTopoQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hCascCosPA"), pT, cand.casccosPA(coll.posX(), coll.posY(), coll.posZ()));
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hDCACascDaughters"), pT, cand.dcacascdaughters());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hCascRadius"), pT, cand.cascradius());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hMesonDCAToPV"), pT, cand.dcanegtopv());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hBaryonDCAToPV"), pT, cand.dcapostopv());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hBachDCAToPV"), pT, cand.dcabachtopv());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hV0CosPA"), pT, cand.v0cosPA(coll.posX(), coll.posY(), coll.posZ()));
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hV0Radius"), pT, cand.v0radius());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hDCAV0Daughters"), pT, cand.dcaV0daughters());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hDCAV0ToPV"), pT, std::fabs(cand.dcav0topv(coll.posX(), coll.posY(), coll.posZ())));
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hMassLambdaDau"), pT, cand.mLambda());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/hCtau"), pT, ctau);
    }
    if (PIDConfigurations.doTPCQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosNsigmaTPC"), centrality, pT, tpcNsigmaPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegNsigmaTPC"), centrality, pT, tpcNsigmaNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachNsigmaTPC"), centrality, pT, tpcNsigmaBach);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTPCsignal"), centrality, pT, posTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTPCsignal"), centrality, pT, negTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachTPCsignal"), centrality, pT, bachTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPtot"), centrality, cand.positivept() * std::cosh(cand.positiveeta()), tpcNsigmaPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPtot"), centrality, cand.negativept() * std::cosh(cand.negativeeta()), tpcNsigmaNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachNsigmaTPCvsTrackPtot"), centrality, cand.bachelorpt() * std::cosh(cand.bacheloreta()), tpcNsigmaBach);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPtot"), centrality, cand.positivept() * std::cosh(cand.positiveeta()), posTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPtot"), centrality, cand.negativept() * std::cosh(cand.negativeeta()), negTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachTPCsignalVsTrackPtot"), centrality, cand.bachelorpt() * std::cosh(cand.bacheloreta()), bachTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPt"), centrality, cand.positivept(), tpcNsigmaPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPt"), centrality, cand.negativept(), tpcNsigmaNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachNsigmaTPCvsTrackPt"), centrality, cand.bachelorpt(), tpcNsigmaBach);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPt"), centrality, cand.positivept(), posTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPt"), centrality, cand.negativept(), negTrackExtra.tpcSignal());
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachTPCsignalVsTrackPt"), centrality, cand.bachelorpt(), bachTrackExtra.tpcSignal());
    }
    if (PIDConfigurations.doTOFQA) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTOFdeltaT"), centrality, pT, tofDeltaTPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTOFdeltaT"), centrality, pT, tofDeltaTNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachTOFdeltaT"), centrality, pT, tofDeltaTBach);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPtot"), centrality, cand.positivept() * std::cosh(cand.positiveeta()), tofDeltaTPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPtot"), centrality, cand.negativept() * std::cosh(cand.negativeeta()), tofDeltaTNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachTOFdeltaTvsTrackPtot"), centrality, cand.bachelorpt() * std::cosh(cand.bacheloreta()), tofDeltaTBach);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPt"), centrality, cand.positivept(), tofDeltaTPos);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPt"), centrality, cand.negativept(), tofDeltaTNeg);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h3dBachTOFdeltaTvsTrackPt"), centrality, cand.bachelorpt(), tofDeltaTBach);
    }
    if (doDetectPropQA == 1) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h8dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, bachDetMap, bachITSclusMap, pT);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h4dBachDetectPropVsCentrality"), centrality, bachTrackExtra.detectorMap(), bachTrackExtra.itsClusterMap(), pT);
    }
    if (doDetectPropQA == 2) {
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h9dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, bachDetMap, bachITSclusMap, pT, invMass);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT, invMass);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT, invMass);
      histos.fill(HIST(kParticlenames[partID]) + HIST("/h5dBachDetectPropVsCentrality"), centrality, bachTrackExtra.detectorMap(), bachTrackExtra.itsClusterMap(), pT, invMass);
    }
  }

  void init(InitContext const&)
  {
    if (doprocessV0s && doprocessCascades) {
      LOG(fatal) << "Unable to analyze both v0s and cascades simultaneously. Please enable only one process at a time";
    }

    if ((doprocessV0sMC || doprocessCascadesMC || doprocessGenerated) && (doprocessV0s || doprocessCascades)) {
      LOG(fatal) << "Cannot analyze both data and MC simultaneously. Please select one of them.";
    }

    rctChecker.init(evSels.cfgEvtRCTFlagCheckerLabel, evSels.cfgEvtRCTFlagCheckerZDCCheck, evSels.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    if (evSels.cfgEvtRCTFlagCheckerFV0Check) {
      rctChecker.set(o2::aod::rctsel::kFV0Bad);
    }

    // initialise bit masks
    setBits(maskTopologicalV0, {selV0CosPA, selDCANegToPV, selDCAPosToPV, selDCAV0Dau, selV0Radius, selV0RadiusMax});
    setBits(maskTopologicalCasc, {selCascCosPA, selDCACascDau, selCascRadius, selCascRadiusMax, selBachToPV, selMesonToPV, selBaryonToPV,
                                  selDCAV0ToPV, selV0CosPA, selDCAV0Dau, selV0Radius, selV0RadiusMax, selLambdaMassWin});

    if (doBachelorBaryonCut)
      maskTopologicalCasc.set(selBachBaryon);

    setBits(maskKinematicV0, {selPosEta, selNegEta});
    setBits(maskKinematicCasc, {selPosEta, selNegEta, selBachEta});

    if (doDaughterDCA) {
      maskKinematicV0.set(selDauDCA);
      maskKinematicCasc.set(selDauDCA);
    }

    // Specific masks
    setBits(maskK0ShortSpecific, {selK0ShortRapidity, selK0ShortCTau, selK0ShortArmenteros, selConsiderK0Short});
    setBits(maskLambdaSpecific, {selLambdaRapidity, selLambdaCTau, selConsiderLambda});
    setBits(maskAntiLambdaSpecific, {selLambdaRapidity, selLambdaCTau, selConsiderAntiLambda});
    setBits(maskXiSpecific, {selXiRapidity, selXiCTau, selRejCompXi, selMassWinXi, selConsiderXi});
    setBits(maskAntiXiSpecific, {selXiRapidity, selXiCTau, selRejCompXi, selMassWinXi, selConsiderAntiXi});
    setBits(maskOmegaSpecific, {selOmegaRapidity, selOmegaCTau, selRejCompOmega, selMassWinOmega, selConsiderOmega});
    setBits(maskAntiOmegaSpecific, {selOmegaRapidity, selOmegaCTau, selRejCompOmega, selMassWinOmega, selConsiderAntiOmega});

    // ask for specific TPC/TOF PID selections
    // positive track
    if (TrackConfigurations.requirePosITSonly) {
      setBits(maskTrackPropertiesV0, {selPosItsOnly, selPosGoodITSTrack});
    } else {
      setBits(maskTrackPropertiesV0, {selPosGoodTPCTrack, selPosGoodITSTrack});
      // TPC signal is available: ask for positive track PID
      if (PIDConfigurations.tpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific.set(selTPCPIDPositivePion);
        maskLambdaSpecific.set(selTPCPIDPositiveProton);
        maskAntiLambdaSpecific.set(selTPCPIDPositivePion);

        maskXiSpecific.set(selTPCPIDPositiveProton);
        maskAntiXiSpecific.set(selTPCPIDPositivePion);
        maskOmegaSpecific.set(selTPCPIDPositiveProton);
        maskAntiOmegaSpecific.set(selTPCPIDPositivePion);
      }
      // TOF PID
      if (PIDConfigurations.tofPidNsigmaCutK0Pi < 1e+5) { // safeguard for no cut
        setBits(maskK0ShortSpecific, {selTOFNSigmaPositivePionK0Short, selTOFDeltaTPositivePionK0Short});
      }
      if (PIDConfigurations.tofPidNsigmaCutLaPr < 1e+5) { // safeguard for no cut
        setBits(maskLambdaSpecific, {selTOFNSigmaPositiveProtonLambda, selTOFDeltaTPositiveProtonLambda});
        setBits(maskXiSpecific, {selTOFNSigmaPositiveProtonLambdaXi, selTOFDeltaTPositiveProtonLambdaXi});
        setBits(maskOmegaSpecific, {selTOFNSigmaPositiveProtonLambdaOmega, selTOFDeltaTPositiveProtonLambdaOmega});
      }
      if (PIDConfigurations.tofPidNsigmaCutLaPi < 1e+5) { // safeguard for no cut
        setBits(maskAntiLambdaSpecific, {selTOFNSigmaPositivePionLambda, selTOFDeltaTPositivePionLambda});
        setBits(maskAntiXiSpecific, {selTOFNSigmaPositivePionLambdaXi, selTOFDeltaTPositivePionLambdaXi});
        setBits(maskAntiOmegaSpecific, {selTOFNSigmaPositivePionLambdaOmega, selTOFDeltaTPositivePionLambdaOmega});
      }
    }
    // negative track
    if (TrackConfigurations.requireNegITSonly) {
      setBits(maskTrackPropertiesV0, {selNegItsOnly, selNegGoodITSTrack});
    } else {
      setBits(maskTrackPropertiesV0, {selNegGoodTPCTrack, selNegGoodITSTrack});
      // TPC signal is available: ask for negative track PID
      if (PIDConfigurations.tpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific.set(selTPCPIDNegativePion);
        maskLambdaSpecific.set(selTPCPIDNegativePion);
        maskAntiLambdaSpecific.set(selTPCPIDNegativeProton);

        maskXiSpecific.set(selTPCPIDNegativePion);
        maskAntiXiSpecific.set(selTPCPIDPositiveProton);
        maskOmegaSpecific.set(selTPCPIDNegativePion);
        maskAntiOmegaSpecific.set(selTPCPIDPositiveProton);
      }
      // TOF PID
      if (PIDConfigurations.tofPidNsigmaCutK0Pi < 1e+5) { // safeguard for no cut
        setBits(maskK0ShortSpecific, {selTOFNSigmaNegativePionK0Short, selTOFDeltaTNegativePionK0Short});
      }
      if (PIDConfigurations.tofPidNsigmaCutLaPr < 1e+5) { // safeguard for no cut
        setBits(maskAntiLambdaSpecific, {selTOFNSigmaNegativeProtonLambda, selTOFDeltaTNegativeProtonLambda});
        setBits(maskAntiXiSpecific, {selTOFNSigmaNegativeProtonLambdaXi, selTOFDeltaTNegativeProtonLambdaXi});
        setBits(maskAntiOmegaSpecific, {selTOFNSigmaNegativeProtonLambdaOmega, selTOFDeltaTNegativeProtonLambdaOmega});
      }
      if (PIDConfigurations.tofPidNsigmaCutLaPi < 1e+5) { // safeguard for no cut
        setBits(maskLambdaSpecific, {selTOFNSigmaNegativePionLambda, selTOFDeltaTNegativePionLambda});
        setBits(maskXiSpecific, {selTOFNSigmaNegativePionLambdaXi, selTOFDeltaTNegativePionLambdaXi});
        setBits(maskOmegaSpecific, {selTOFNSigmaNegativePionLambdaOmega, selTOFDeltaTNegativePionLambdaOmega});
      }
    }
    // bachelor track
    maskTrackPropertiesCasc = maskTrackPropertiesV0;
    if (TrackConfigurations.requireBachITSonly) {
      setBits(maskTrackPropertiesCasc, {selBachItsOnly, selBachGoodITSTrack});
    } else {
      setBits(maskTrackPropertiesCasc, {selBachGoodTPCTrack, selBachGoodITSTrack});
      // TPC signal is available: ask for positive track PID
      if (PIDConfigurations.tpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskXiSpecific.set(selTPCPIDBachPion);
        maskAntiXiSpecific.set(selTPCPIDBachPion);
        maskOmegaSpecific.set(selTPCPIDBachKaon);
        maskAntiOmegaSpecific.set(selTPCPIDBachKaon);
      }
      // TOF PID
      if (PIDConfigurations.tofPidNsigmaCutXiPi < 1e+5) { // safeguard for no cut
        setBits(maskXiSpecific, {selTOFNSigmaBachPionXi, selTOFDeltaTBachPionXi});
        setBits(maskAntiXiSpecific, {selTOFNSigmaBachPionXi, selTOFDeltaTBachPionXi});
      }
      if (PIDConfigurations.tofPidNsigmaCutOmegaKaon < 1e+5) { // safeguard for no cut
        setBits(maskOmegaSpecific, {selTOFNSigmaBachKaonOmega, selTOFDeltaTBachKaonOmega});
        setBits(maskAntiOmegaSpecific, {selTOFNSigmaBachKaonOmega, selTOFDeltaTBachKaonOmega});
      }
    }

    if (TrackConfigurations.skipTPConly) {
      setBits(maskK0ShortSpecific, {selPosNotTPCOnly, selNegNotTPCOnly});
      setBits(maskLambdaSpecific, {selPosNotTPCOnly, selNegNotTPCOnly});
      setBits(maskAntiLambdaSpecific, {selPosNotTPCOnly, selNegNotTPCOnly});
      setBits(maskXiSpecific, {selPosNotTPCOnly, selNegNotTPCOnly, selBachNotTPCOnly});
      setBits(maskOmegaSpecific, {selPosNotTPCOnly, selNegNotTPCOnly, selBachNotTPCOnly});
      setBits(maskAntiXiSpecific, {selPosNotTPCOnly, selNegNotTPCOnly, selBachNotTPCOnly});
      setBits(maskAntiOmegaSpecific, {selPosNotTPCOnly, selNegNotTPCOnly, selBachNotTPCOnly});
    }

    // Primary particle selection, central to analysis
    maskSelectionK0Short = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskK0ShortSpecific | (std::bitset<kSelNum>(1) << selPhysPrimK0Short);
    maskSelectionLambda = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskLambdaSpecific | (std::bitset<kSelNum>(1) << selPhysPrimLambda);
    maskSelectionAntiLambda = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskAntiLambdaSpecific | (std::bitset<kSelNum>(1) << selPhysPrimAntiLambda);
    maskSelectionXi = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskXiSpecific | (std::bitset<kSelNum>(1) << selPhysPrimXi);
    maskSelectionAntiXi = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskAntiXiSpecific | (std::bitset<kSelNum>(1) << selPhysPrimAntiXi);
    maskSelectionOmega = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskOmegaSpecific | (std::bitset<kSelNum>(1) << selPhysPrimOmega);
    maskSelectionAntiOmega = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskAntiOmegaSpecific | (std::bitset<kSelNum>(1) << selPhysPrimAntiOmega);

    // No primary requirement for feeddown matrix
    secondaryMaskSelectionLambda = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskLambdaSpecific;
    secondaryMaskSelectionAntiLambda = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskAntiLambdaSpecific;

    // Event Counter
    histos.add("eventQA/hEventSelection", "hEventSelection", kTH1D, {{17, -0.5f, 16.5f}});
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(2, "kIsTriggerTVX");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(11, "kNoCollInTimeRangeStd");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeNarrow");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(13, "Below min occup.");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(14, "Above max occup.");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(15, "RCTFlagsChecker");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(16, "isUPC");
    histos.get<TH1>(HIST("eventQA/hEventSelection"))->GetXaxis()->SetBinLabel(17, "has UPC flag");

    // Event QA
    histos.add("eventQA/hCentralityVsFT0ampl", "hCentralityVsFT0ampl", kTH3D, {axisFT0Cqa, axisDetectors.axisFT0Aampl, axisSelGap});
    histos.add("eventQA/hCentrality", "hCentrality", kTH2D, {axisFT0Cqa, axisSelGap});
    histos.add("eventQA/hFT0ampl", "hFT0ampl", kTH2D, {axisDetectors.axisFT0Aampl, axisSelGap});
    histos.add("eventQA/hCentralityVsTracksPVeta1", "hCentralityVsTracksPVeta1", kTH3D, {axisFT0Cqa, axisNTracksPVeta1, axisSelGap});
    histos.add("eventQA/hCentralityVsTracksTotalExceptITSonly", "hCentralityVsTracksTotalExceptITSonly", kTH3D, {axisFT0Cqa, axisNTracksTotalExceptITSonly, axisSelGap});
    histos.add("eventQA/hOccupancy", "hOccupancy", kTH2D, {axisOccupancy, axisSelGap});
    histos.add("eventQA/hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH3D, {axisFT0Cqa, axisOccupancy, axisSelGap});
    histos.add("eventQA/hTracksPVeta1VsTracksGlobal", "hTracksPVeta1VsTracksGlobal", kTH3D, {axisNTracksPVeta1, axisNTracksGlobal, axisSelGap});
    histos.add("eventQA/hCentralityVsTracksGlobal", "hCentralityVsTracksGlobal", kTH3D, {axisFT0Cqa, axisNTracksGlobal, axisSelGap});
    histos.add("eventQA/hFT0amplVsTracksGlobal", "hFT0amplVsTracksGlobal", kTH3D, {axisDetectors.axisFT0Aampl, axisNTracksGlobal, axisSelGap});
    histos.add("eventQA/hRawGapSide", "Raw Gap side; Entries", kTH1D, {{6, -1.5, 4.5}});
    histos.add("eventQA/hSelGapSide", "Selected gap side (with n); Entries", kTH1D, {axisSelGap});
    histos.add("eventQA/hPosX", "Vertex position in x", kTH2D, {{100, -0.1, 0.1}, axisSelGap});
    histos.add("eventQA/hPosY", "Vertex position in y", kTH2D, {{100, -0.1, 0.1}, axisSelGap});
    histos.add("eventQA/hPosZ", "Vertex position in z", kTH2D, {{100, -20., 20.}, axisSelGap});
    histos.add("eventQA/hFT0", "hFT0", kTH3D, {axisDetectors.axisFT0Aampl, axisDetectors.axisFT0Campl, axisSelGap});
    histos.add("eventQA/hFDD", "hFDD", kTH3D, {axisDetectors.axisFDDAampl, axisDetectors.axisFDDCampl, axisSelGap});
    histos.add("eventQA/hZN", "hZN", kTH3D, {axisDetectors.axisZNAampl, axisDetectors.axisZNCampl, axisSelGap});

    if (doprocessGenerated) {
      histos.add("eventQA/mc/hEventSelectionMC", "hEventSelectionMC", kTH3D, {{3, -0.5, 2.5}, axisNTracksPVeta1, axisGeneratorIds});
      histos.get<TH3>(HIST("eventQA/mc/hEventSelectionMC"))->GetXaxis()->SetBinLabel(1, "All collisions");
      histos.get<TH3>(HIST("eventQA/mc/hEventSelectionMC"))->GetXaxis()->SetBinLabel(2, "posZ cut");
      histos.get<TH3>(HIST("eventQA/mc/hEventSelectionMC"))->GetXaxis()->SetBinLabel(3, "rec. at least once");
      histos.add("eventQA/mc/hTracksGlobalvsMCNParticlesEta08gen", "hTracksGlobalvsMCNParticlesEta08gen", kTH2D, {axisNTracksGlobal, axisNTracksGlobal});
      histos.add("eventQA/mc/hTracksGlobalVsNcoll_beforeEvSel", "hTracksGlobalVsNcoll_beforeEvSel", kTH2D, {axisNTracksGlobal, axisNAssocColl});
      histos.add("eventQA/mc/hTracksGlobalVsNcoll_afterEvSel", "hTracksGlobalVsNcoll_afterEvSel", kTH2D, {axisNTracksGlobal, axisNAssocColl});
      histos.add("eventQA/mc/hTracksGlobalVsPVzMC", "hTracksGlobalVsPVzMC", kTH2D, {axisNTracksGlobal, {100, -20., 20.}});
      histos.add("eventQA/mc/hEventPVzMC", "hEventPVzMC", kTH1D, {{100, -20., 20.}});
      histos.add("eventQA/mc/hGenEventFT0ampl", "hGenEventFT0ampl", kTH1D, {axisDetectors.axisFT0Aampl});
      histos.add("eventQA/mc/hGenEventCentrality", "hGenEventCentrality", kTH1D, {axisFT0Cqa});
      histos.add("eventQA/mc/hGeneratorsId", "hGeneratorsId", kTH1D, {axisGeneratorIds});
      histos.add("eventQA/mc/hSelGeneratorsId", "hSelGeneratorsId", kTH1D, {axisGeneratorIds});
    }

    if (doprocessV0sMC || doprocessCascadesMC) {
      // Event QA
      histos.add("eventQA/mc/hFakeEvents", "hFakeEvents", {kTH1D, {{1, -0.5f, 0.5f}}});
      histos.add("eventQA/mc/hNTracksGlobalvsMCNParticlesEta08rec", "hNTracksGlobalvsMCNParticlesEta08rec", kTH2D, {axisNTracksGlobal, axisNTracksGlobal});
      histos.add("eventQA/mc/hNTracksPVeta1vsMCNParticlesEta10rec", "hNTracksPVeta1vsMCNParticlesEta10rec", kTH2D, {axisNTracksPVeta1, axisNTracksPVeta1});
      histos.add("eventQA/mc/hNTracksGlobalvstotalMultMCParticles", "hNTracksGlobalvstotalMultMCParticles", kTH2D, {axisNTracksGlobal, axisNchInvMass});
      histos.add("eventQA/mc/hNTracksPVeta1vstotalMultMCParticles", "hNTracksPVeta1vstotalMultMCParticles", kTH2D, {axisNTracksPVeta1, axisNchInvMass});
      histos.add("eventQA/hSelGapSideNoNeutrons", "Selected gap side (no n); Entries", kTH1D, {{5, -0.5, 4.5}});
    }

    if (doprocessV0sMC) {
      if (analyseLambda && calculateFeeddownMatrix)
        histos.add(Form("%s/h3dLambdaFeeddown", kParticlenames[1].data()), "h3dLambdaFeeddown", kTH3F, {axisNTracksGlobal, axisPt, axisPt});
      if (analyseAntiLambda && calculateFeeddownMatrix)
        histos.add(Form("%s/h3dAntiLambdaFeeddown", kParticlenames[2].data()), "h3dAntiLambdaFeeddown", kTH3F, {axisNTracksGlobal, axisPt, axisPt});
    }

    if (doprocessGenerated) {
      for (int partID = 0; partID <= 6; partID++) {
        histos.add(Form("%s/mc/h7dGen", kParticlenames[partID].data()), "h7dGen", kTHnSparseF, {axisDetectors.axisFT0ampl, axisNchInvMass, axisNchInvMass, axisPt, axisSelGap, axisRap, axisGeneratorIds});
      }
    }

    if (doprocessV0s || doprocessV0sMC) {
      // For all candidates
      if (doPlainTopoQA) {
        histos.add("generalQA/hPt", "hPt", kTH1F, {axisPtCoarse});
        histos.add("generalQA/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("generalQA/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("generalQA/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
        histos.add("generalQA/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
        histos.add("generalQA/hCosPA", "hCosPA", kTH1F, {axisCosPA});
        histos.add("generalQA/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
        histos.add("generalQA/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
        histos.add("generalQA/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
        histos.add("generalQA/h2dArmenterosAll", "h2dArmenterosAll", kTH2F, {axisAPAlpha, axisAPQt});
        histos.add("generalQA/h2dArmenterosSelected", "h2dArmenterosSelected", kTH2F, {axisAPAlpha, axisAPQt});
      }

      // K0s
      if (analyseK0Short) {
        addHistograms<0>(histos);
      }

      // Lambda
      if (analyseLambda) {
        addHistograms<1>(histos);
      }

      // Anti-Lambda
      if (analyseAntiLambda) {
        addHistograms<2>(histos);
      }
    }

    if (doprocessCascades || doprocessCascadesMC) {
      // For all candidates
      if (doPlainTopoQA) {
        histos.add("generalQA/hPt", "hPt", kTH1F, {axisPtCoarse});
        histos.add("generalQA/hCascCosPA", "hCascCosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
        histos.add("generalQA/hDCACascDaughters", "hDCACascDaughters", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
        histos.add("generalQA/hCascRadius", "hCascRadius", kTH2D, {axisPtCoarse, {500, 0.0f, 50.0f}});
        histos.add("generalQA/hMesonDCAToPV", "hMesonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
        histos.add("generalQA/hBaryonDCAToPV", "hBaryonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
        histos.add("generalQA/hBachDCAToPV", "hBachDCAToPV", kTH2F, {axisPtCoarse, {200, -1.0f, 1.0f}});
        histos.add("generalQA/hV0CosPA", "hV0CosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
        histos.add("generalQA/hV0Radius", "hV0Radius", kTH2D, {axisPtCoarse, axisV0Radius});
        histos.add("generalQA/hDCAV0Daughters", "hDCAV0Daughters", kTH2F, {axisPtCoarse, axisDCAdau});
        histos.add("generalQA/hDCAV0ToPV", "hDCAV0ToPV", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
        histos.add("generalQA/hMassLambdaDau", "hMassLambdaDau", kTH2F, {axisPtCoarse, axisLambdaMass});
        histos.add("generalQA/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
        histos.add("generalQA/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
        histos.add("generalQA/h2dBachITSvsTPCpts", "h2dBachITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      }

      // Xi
      if (analyseXi) {
        addHistograms<3>(histos);
      }

      // Anti-Xi
      if (analyseAntiXi) {
        addHistograms<4>(histos);
      }

      // Omega
      if (analyseOmega) {
        addHistograms<5>(histos);
      }

      // Anti-Omega
      if (analyseAntiOmega) {
        addHistograms<6>(histos);
      }
    }

    if (verbose) {
      histos.print();
    }
  }

  template <typename TCollision>
  int getGapSide(TCollision const& collision)
  {
    int selGapSide = sgSelector.trueGap(collision, upcCuts.fv0a, upcCuts.ft0a, upcCuts.ft0c, upcCuts.zdc);
    return selGapSide;
  }

  template <typename TCollision>
  void fillHistogramsQA(TCollision const& collision, int const& gap)
  {
    // QA histograms
    float centrality = collision.centFT0C();
    float ft0ampl = -1.f;

    if (gap == 0) {
      ft0ampl = collision.totalFT0AmplitudeC();
    } else if (gap == 1) {
      ft0ampl = collision.totalFT0AmplitudeA();
    }

    histos.fill(HIST("eventQA/hCentralityVsFT0ampl"), centrality, ft0ampl, gap);
    histos.fill(HIST("eventQA/hCentrality"), centrality, gap);
    histos.fill(HIST("eventQA/hFT0ampl"), ft0ampl, gap);
    histos.fill(HIST("eventQA/hCentralityVsTracksTotalExceptITSonly"), centrality, collision.multAllTracksTPCOnly() + collision.multAllTracksITSTPC(), gap);
    histos.fill(HIST("eventQA/hCentralityVsTracksPVeta1"), centrality, collision.multNTracksPVeta1(), gap);
    histos.fill(HIST("eventQA/hOccupancy"), collision.trackOccupancyInTimeRange(), gap);
    histos.fill(HIST("eventQA/hCentralityVsOccupancy"), centrality, collision.trackOccupancyInTimeRange(), gap);
    histos.fill(HIST("eventQA/hTracksPVeta1VsTracksGlobal"), collision.multNTracksPVeta1(), collision.multNTracksGlobal(), gap);
    histos.fill(HIST("eventQA/hCentralityVsTracksGlobal"), centrality, collision.multNTracksGlobal(), gap);
    histos.fill(HIST("eventQA/hFT0amplVsTracksGlobal"), ft0ampl, collision.multNTracksGlobal(), gap);
    histos.fill(HIST("eventQA/hPosX"), collision.posX(), gap);
    histos.fill(HIST("eventQA/hPosY"), collision.posY(), gap);
    histos.fill(HIST("eventQA/hPosZ"), collision.posZ(), gap);

    histos.fill(HIST("eventQA/hSelGapSide"), gap);
    histos.fill(HIST("eventQA/hFT0"), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), gap);
    histos.fill(HIST("eventQA/hFDD"), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(), gap);

    auto zna = collision.energyCommonZNA();
    auto znc = collision.energyCommonZNC();
    constexpr float inf_f = std::numeric_limits<float>::infinity();

    if (zna == -inf_f)
      histos.fill(HIST("eventQA/hZN"), -1, znc, gap);
    else if (znc == -inf_f)
      histos.fill(HIST("eventQA/hZN"), zna, -1, gap);
    else if (zna == -999 && znc == -999)
      histos.fill(HIST("eventQA/hZN"), -2, -2, gap);
    else if (zna == -999 || znc == -999)
      LOG(warning) << "Only one ZDC signal is -999";
  }

  template <typename TCollision>
  bool acceptEvent(TCollision const& collision, bool fillQA)
  {
    struct SelectionCheck {
      bool selection;
      bool condition;
      float qaBin;
    };

    const std::array<SelectionCheck, 15> selections = {{
      {true, true, 0.0},                                                                                                // All collisions
      {evSels.requireIsTriggerTVX, collision.selection_bit(aod::evsel::kIsTriggerTVX), 1.0},                            // Triggered by FT0M
      {true, std::fabs(collision.posZ()) < maxZVtxPosition, 2.0},                                                       // Vertex-Z selected
      {evSels.rejectITSROFBorder, collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder), 3.0},                   // Not at ITS ROF border
      {evSels.rejectTFBorder, collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder), 4.0},                        // Not at TF border
      {evSels.requireIsVertexITSTPC, collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC), 5.0},                    // At least one ITS-TPC track
      {evSels.requireIsGoodZvtxFT0VsPV, collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV), 6.0},              // PV position consistency
      {evSels.requireIsVertexTOFmatched, collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched), 7.0},            // PV with TOF match
      {evSels.requireIsVertexTRDmatched, collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched), 8.0},            // PV with TRD match
      {evSels.rejectSameBunchPileup, collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup), 9.0},                 // No same-bunch pileup
      {evSels.requireNoCollInTimeRangeStd, collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard), 10.0},  // No collision within +-10 µs
      {evSels.requireNoCollInTimeRangeNarrow, collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow), 11.0}, // No collision within +-4 µs
      {minOccupancy >= 0, collision.trackOccupancyInTimeRange() >= minOccupancy, 12.0},                                 // Above min occupancy
      {maxOccupancy > 0, collision.trackOccupancyInTimeRange() < maxOccupancy, 13.0},                                   // Below max occupancy
      {evSels.requireRCTFlagChecker, rctChecker(collision), 14.0},                                                      // Verify collision in RCT
    }};

    for (const auto& sel : selections) {
      if (sel.selection && !sel.condition) {
        return false;
      }
      if (fillQA && sel.selection) {
        histos.fill(HIST("eventQA/hEventSelection"), sel.qaBin);
      }
    }

    if (evSels.studyUPConly && !collision.isUPC()) {
      return false;
    } else if (collision.isUPC() && fillQA) {
      histos.fill(HIST("eventQA/hEventSelection"), 15.0); // is UPC compatible
    }

    // Additional check for UPC collision flag
    if (evSels.useUPCflag && collision.flags() < 1) {
      return false;
    }
    if (collision.flags() >= 1 && fillQA) {
      histos.fill(HIST("eventQA/hEventSelection"), 16.0); // UPC event
    }

    return true;
  }

  bool verifyMask(std::bitset<kSelNum> bitmap, std::bitset<kSelNum> mask)
  {
    return (bitmap & mask) == mask;
  }

  int computeITSclusBitmap(uint8_t itsClusMap, bool fromAfterburner)
  {
    int bitMap = 0;

    struct MaskBitmapPair {
      uint8_t mask;
      int bitmap;
      int afterburnerBitmap;
    };

    constexpr MaskBitmapPair kConfigs[] = {
      // L6 <-- L0
      {0x7F, 12, 12}, // 01111 111 (L0 to L6)
      {0x7E, 11, 11}, // 01111 110 (L1 to L6)
      {0x7C, 10, 10}, // 01111 100 (L2 to L6)
      {0x78, 9, -3},  // 01111 000 (L3 to L6)
      {0x70, 8, -2},  // 01110 000 (L4 to L6)
      {0x60, 7, -1},  // 01100 000 (L5 to L6)
      {0x3F, 6, 6},   // 00111 111 (L0 to L5)
      {0x3E, 5, 5},   // 00111 110 (L1 to L5)
      {0x3C, 4, 4},   // 00111 100 (L2 to L5)
      {0x1F, 3, 3},   // 00011 111 (L0 to L4)
      {0x1E, 2, 2},   // 00011 110 (L1 to L4)
      {0x0F, 1, 1},   // 00001 111 (L0 to L3)
    };

    for (const auto& config : kConfigs) {
      if (verifyMask(itsClusMap, config.mask)) {
        bitMap = fromAfterburner ? config.afterburnerBitmap : config.bitmap;
        break;
      }
    }

    return bitMap;
  }

  uint computeDetBitmap(uint8_t detMap)
  {
    uint bitMap = 0;

    struct MaskBitmapPair {
      uint8_t mask;
      int bitmap;
    };

    constexpr MaskBitmapPair kConfigs[] = {
      {o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD | o2::aod::track::TOF, 4}, // ITS-TPC-TRD-TOF
      {o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TOF, 3},                       // ITS-TPC-TOF
      {o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD, 2},                       // ITS-TPC-TRD
      {o2::aod::track::ITS | o2::aod::track::TPC, 1}                                              // ITS-TPC
    };

    for (const auto& config : kConfigs) {
      if (verifyMask(detMap, config.mask)) {
        bitMap = config.bitmap;
        break;
      }
    }

    return bitMap;
  }

  template <typename TCasc, typename TCollision>
  std::bitset<kSelNum> computeBitmapCascade(TCasc const& casc, TCollision const& coll)
  {
    float rapidityXi = casc.yXi();
    float rapidityOmega = casc.yOmega();

    // Access daughter tracks
    auto posTrackExtra = casc.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = casc.template negTrackExtra_as<DauTracks>();
    auto bachTrackExtra = casc.template bachTrackExtra_as<DauTracks>();

    // c x tau
    float decayPos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
    float totalMom = std::hypot(casc.px(), casc.py(), casc.pz());
    float ctauXi = totalMom != 0 ? o2::constants::physics::MassXiMinus * decayPos / totalMom : 1e6;
    float ctauOmega = totalMom != 0 ? o2::constants::physics::MassOmegaMinus * decayPos / totalMom : 1e6;

    std::bitset<kSelNum> bitMap = 0;

    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) > casccuts.casccospa)
      bitMap.set(selCascCosPA);
    if (casc.dcacascdaughters() < casccuts.dcacascdau)
      bitMap.set(selDCACascDau);
    if (casc.cascradius() > casccuts.cascradius)
      bitMap.set(selCascRadius);
    if (casc.cascradius() < casccuts.cascradiusMax)
      bitMap.set(selCascRadiusMax);
    if (doBachelorBaryonCut && (casc.bachBaryonCosPA() < casccuts.bachbaryoncospa) && (std::fabs(casc.bachBaryonDCAxyToPV()) > casccuts.bachbaryondcaxytopv))
      bitMap.set(selBachBaryon);
    if (std::fabs(casc.dcabachtopv()) > casccuts.dcabachtopv)
      bitMap.set(selBachToPV);

    if (casc.sign() > 0) {
      if (std::fabs(casc.dcanegtopv()) > casccuts.dcabaryontopv)
        bitMap.set(selBaryonToPV);
      if (std::fabs(casc.dcapostopv()) > casccuts.dcamesontopv)
        bitMap.set(selMesonToPV);
    } else { // no sign == 0, in principle
      if (std::fabs(casc.dcapostopv()) > casccuts.dcabaryontopv)
        bitMap.set(selBaryonToPV);
      if (std::fabs(casc.dcanegtopv()) > casccuts.dcamesontopv)
        bitMap.set(selMesonToPV);
    }

    if (std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) < casccuts.masswin)
      bitMap.set(selMassWinXi);
    if (std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < casccuts.masswin)
      bitMap.set(selMassWinOmega);
    if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) < casccuts.lambdamasswin)
      bitMap.set(selLambdaMassWin);

    if (std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) > casccuts.dcav0topv)
      bitMap.set(selDCAV0ToPV);
    if (casc.v0radius() > v0cuts.v0radius)
      bitMap.set(selV0Radius);
    if (casc.v0radius() < v0cuts.v0radiusMax)
      bitMap.set(selV0RadiusMax);
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) > v0cuts.v0cospa)
      bitMap.set(selV0CosPA);
    if (casc.dcaV0daughters() < v0cuts.dcav0dau)
      bitMap.set(selDCAV0Dau);

    // proper lifetime
    if (ctauXi < nCtauCutCasc->get("lifetimecutXi") * ctauxiPDG)
      bitMap.set(selXiCTau);
    if (ctauOmega < nCtauCutCasc->get("lifetimecutOmega") * ctauomegaPDG)
      bitMap.set(selOmegaCTau);

    auto poseta = RecoDecay::eta(std::array{casc.pxpos(), casc.pypos(), casc.pzpos()});
    auto negeta = RecoDecay::eta(std::array{casc.pxneg(), casc.pyneg(), casc.pzneg()});
    auto bacheta = RecoDecay::eta(std::array{casc.pxbach(), casc.pybach(), casc.pzbach()});

    // kinematic
    if (std::fabs(rapidityXi) < rapidityCut)
      bitMap.set(selXiRapidity);
    if (std::fabs(rapidityOmega) < rapidityCut)
      bitMap.set(selOmegaRapidity);
    if (std::fabs(poseta) < daughterEtaCut)
      bitMap.set(selNegEta);
    if (std::fabs(negeta) < daughterEtaCut)
      bitMap.set(selPosEta);
    if (std::fabs(bacheta) < daughterEtaCut)
      bitMap.set(selBachEta);

    // DCA cuts
    auto pospt = std::sqrt(std::pow(casc.pxpos(), 2) + std::pow(casc.pypos(), 2));
    auto negpt = std::sqrt(std::pow(casc.pxneg(), 2) + std::pow(casc.pyneg(), 2));
    auto bachpt = std::sqrt(std::pow(casc.pxbach(), 2) + std::pow(casc.pybach(), 2));

    double posDcaXYLimit = 0.0105f + 0.035f / std::pow(pospt, 1.1f);
    double negDcaXYLimit = 0.0105f + 0.035f / std::pow(negpt, 1.1f);
    double bachDcaXYLimit = 0.0105f + 0.035f / std::pow(bachpt, 1.1f);

    // TODO: separate xy and z //
    if ((std::abs(casc.dcapostopv()) > posDcaXYLimit) &&
        (std::abs(casc.dcanegtopv()) > negDcaXYLimit) &&
        (std::abs(casc.dcabachtopv()) > bachDcaXYLimit)) {
      bitMap.set(selDauDCA);
    }

    // ITS quality flags
    if (posTrackExtra.itsNCls() >= TrackConfigurations.minITSclusters)
      bitMap.set(selPosGoodITSTrack);
    if (negTrackExtra.itsNCls() >= TrackConfigurations.minITSclusters)
      bitMap.set(selNegGoodITSTrack);
    if (bachTrackExtra.itsNCls() >= TrackConfigurations.minITSclusters)
      bitMap.set(selBachGoodITSTrack);

    // TPC quality flags
    if (posTrackExtra.tpcCrossedRows() >= TrackConfigurations.minTPCrows)
      bitMap.set(selPosGoodTPCTrack);
    if (negTrackExtra.tpcCrossedRows() >= TrackConfigurations.minTPCrows)
      bitMap.set(selNegGoodTPCTrack);
    if (bachTrackExtra.tpcCrossedRows() >= TrackConfigurations.minTPCrows)
      bitMap.set(selBachGoodTPCTrack);

    // TPC PID
    // positive track
    if (std::fabs(posTrackExtra.tpcNSigmaPi()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositivePion);
    if (std::fabs(posTrackExtra.tpcNSigmaPr()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositiveProton);
    // negative track
    if (std::fabs(negTrackExtra.tpcNSigmaPi()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativePion);
    if (std::fabs(negTrackExtra.tpcNSigmaPr()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativeProton);
    // bachelor track
    if (std::fabs(bachTrackExtra.tpcNSigmaPi()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDBachPion);
    if (std::fabs(bachTrackExtra.tpcNSigmaKa()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDBachKaon);

    // TOF PID in DeltaT
    // positive track
    if (std::fabs(casc.posTOFDeltaTXiPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTPositiveProtonLambdaXi);
    if (std::fabs(casc.posTOFDeltaTXiPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionLambdaXi);
    if (std::fabs(casc.posTOFDeltaTOmPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTPositiveProtonLambdaOmega);
    if (std::fabs(casc.posTOFDeltaTOmPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionLambdaOmega);
    // negative track
    if (std::fabs(casc.negTOFDeltaTXiPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTNegativeProtonLambdaXi);
    if (std::fabs(casc.negTOFDeltaTXiPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionLambdaXi);
    if (std::fabs(casc.negTOFDeltaTOmPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTNegativeProtonLambdaOmega);
    if (std::fabs(casc.negTOFDeltaTOmPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionLambdaOmega);
    // bachelor track
    if (std::fabs(casc.bachTOFDeltaTOmKa()) < PIDConfigurations.maxDeltaTimeKaon)
      bitMap.set(selTOFDeltaTBachKaonOmega);
    if (std::fabs(casc.bachTOFDeltaTXiPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTBachPionXi);

    // TOF PID in NSigma
    // meson track
    if (std::fabs(casc.tofNSigmaXiLaPi()) < PIDConfigurations.tofPidNsigmaCutLaPi) {
      bitMap.set(selTOFNSigmaPositivePionLambdaXi);
      bitMap.set(selTOFNSigmaNegativePionLambdaXi);
    }
    if (std::fabs(casc.tofNSigmaOmLaPi()) < PIDConfigurations.tofPidNsigmaCutLaPi) {
      bitMap.set(selTOFNSigmaPositivePionLambdaOmega);
      bitMap.set(selTOFNSigmaNegativePionLambdaOmega);
    }
    // baryon track
    if (std::fabs(casc.tofNSigmaXiLaPr()) < PIDConfigurations.tofPidNsigmaCutLaPr) {
      bitMap.set(selTOFNSigmaNegativeProtonLambdaXi);
      bitMap.set(selTOFNSigmaPositiveProtonLambdaXi);
    }
    if (std::fabs(casc.tofNSigmaOmLaPr()) < PIDConfigurations.tofPidNsigmaCutLaPr) {
      bitMap.set(selTOFNSigmaNegativePionLambdaOmega);
      bitMap.set(selTOFNSigmaPositivePionLambdaOmega);
    }
    // bachelor track
    if (std::fabs(casc.tofNSigmaXiPi()) < PIDConfigurations.tofPidNsigmaCutXiPi) {
      bitMap.set(selTOFNSigmaBachPionXi);
    }
    if (std::fabs(casc.tofNSigmaOmKa()) < PIDConfigurations.tofPidNsigmaCutOmegaKaon) {
      bitMap.set(selTOFNSigmaBachKaonOmega);
    }

    // ITS only tag
    if (posTrackExtra.tpcCrossedRows() < 1)
      bitMap.set(selPosItsOnly);
    if (negTrackExtra.tpcCrossedRows() < 1)
      bitMap.set(selNegItsOnly);
    if (bachTrackExtra.tpcCrossedRows() < 1)
      bitMap.set(selBachItsOnly);

    // rej. comp.
    if (std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) > casccuts.rejcomp)
      bitMap.set(selRejCompXi);
    if (std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) > casccuts.rejcomp)
      bitMap.set(selRejCompOmega);

    // TPC only tag
    if (posTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitMap.set(selPosNotTPCOnly);
    if (negTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitMap.set(selNegNotTPCOnly);
    if (bachTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitMap.set(selBachNotTPCOnly);

    return bitMap;
  }

  template <typename TV0, typename TCollision>
  std::bitset<kSelNum> computeBitmapV0(TV0 const& v0, TCollision const& collision)
  {
    float rapidityLambda = v0.yLambda();
    float rapidityK0Short = v0.yK0Short();

    std::bitset<kSelNum> bitMap = 0;

    // base topological variables
    if (v0.v0radius() > v0cuts.v0radius)
      bitMap.set(selV0Radius);
    if (v0.v0radius() < v0cuts.v0radiusMax)
      bitMap.set(selV0RadiusMax);
    if (std::fabs(v0.dcapostopv()) > v0cuts.dcapostopv)
      bitMap.set(selDCAPosToPV);
    if (std::fabs(v0.dcanegtopv()) > v0cuts.dcanegtopv)
      bitMap.set(selDCANegToPV);
    if (v0.v0cosPA() > v0cuts.v0cospa)
      bitMap.set(selV0CosPA);
    if (v0.dcaV0daughters() < v0cuts.dcav0dau)
      bitMap.set(selDCAV0Dau);

    // kinematic
    if (std::fabs(rapidityLambda) < rapidityCut)
      bitMap.set(selLambdaRapidity);
    if (std::fabs(rapidityK0Short) < rapidityCut)
      bitMap.set(selK0ShortRapidity);
    if (std::fabs(v0.negativeeta()) < daughterEtaCut)
      bitMap.set(selNegEta);
    if (std::fabs(v0.positiveeta()) < daughterEtaCut)
      bitMap.set(selPosEta);

    // c x tau
    float ctauK0short = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;

    auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

    // DCA cuts
    auto pospt = std::sqrt(std::pow(v0.pxpos(), 2) + std::pow(v0.pypos(), 2));
    auto negpt = std::sqrt(std::pow(v0.pxneg(), 2) + std::pow(v0.pyneg(), 2));

    double posDcaXYLimit = 0.0105f + 0.035f / std::pow(pospt, 1.1f);
    double negDcaXYLimit = 0.0105f + 0.035f / std::pow(negpt, 1.1f);

    // TODO: separate xy and z //
    if ((std::abs(v0.dcapostopv()) > posDcaXYLimit) &&
        (std::abs(v0.dcanegtopv()) > negDcaXYLimit)) {
      bitMap.set(selDauDCA);
    }

    // ITS quality flags
    if (posTrackExtra.itsNCls() >= TrackConfigurations.minITSclusters)
      bitMap.set(selPosGoodITSTrack);
    if (negTrackExtra.itsNCls() >= TrackConfigurations.minITSclusters)
      bitMap.set(selNegGoodITSTrack);

    // TPC quality flags
    if (posTrackExtra.tpcCrossedRows() >= TrackConfigurations.minTPCrows)
      bitMap.set(selPosGoodTPCTrack);
    if (negTrackExtra.tpcCrossedRows() >= TrackConfigurations.minTPCrows)
      bitMap.set(selNegGoodTPCTrack);

    // TPC PID
    if (std::fabs(posTrackExtra.tpcNSigmaPi()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositivePion);
    if (std::fabs(posTrackExtra.tpcNSigmaPr()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositiveProton);
    if (std::fabs(negTrackExtra.tpcNSigmaPi()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativePion);
    if (std::fabs(negTrackExtra.tpcNSigmaPr()) < PIDConfigurations.tpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativeProton);

    // TOF PID in DeltaT
    // positive track
    if (std::fabs(v0.posTOFDeltaTLaPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTPositiveProtonLambda);
    if (std::fabs(v0.posTOFDeltaTLaPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionLambda);
    if (std::fabs(v0.posTOFDeltaTK0Pi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionK0Short);
    // negative track
    if (std::fabs(v0.negTOFDeltaTLaPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTNegativeProtonLambda);
    if (std::fabs(v0.negTOFDeltaTLaPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionLambda);
    if (std::fabs(v0.negTOFDeltaTK0Pi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionK0Short);

    // TOF PID in NSigma
    // positive track
    if (std::fabs(v0.tofNSigmaLaPr()) < PIDConfigurations.tofPidNsigmaCutLaPr)
      bitMap.set(selTOFNSigmaPositiveProtonLambda);
    if (std::fabs(v0.tofNSigmaALaPi()) < PIDConfigurations.tofPidNsigmaCutLaPi)
      bitMap.set(selTOFNSigmaPositivePionLambda);
    if (std::fabs(v0.tofNSigmaK0PiPlus()) < PIDConfigurations.tofPidNsigmaCutK0Pi)
      bitMap.set(selTOFNSigmaPositivePionK0Short);
    // negative track
    if (std::fabs(v0.tofNSigmaALaPr()) < PIDConfigurations.tofPidNsigmaCutLaPr)
      bitMap.set(selTOFNSigmaNegativeProtonLambda);
    if (std::fabs(v0.tofNSigmaLaPi()) < PIDConfigurations.tofPidNsigmaCutLaPi)
      bitMap.set(selTOFNSigmaNegativePionLambda);
    if (std::fabs(v0.tofNSigmaK0PiMinus()) < PIDConfigurations.tofPidNsigmaCutK0Pi)
      bitMap.set(selTOFNSigmaNegativePionK0Short);

    // ITS only tag
    if (posTrackExtra.tpcCrossedRows() < 1)
      bitMap.set(selPosItsOnly);
    if (negTrackExtra.tpcCrossedRows() < 1)
      bitMap.set(selNegItsOnly);

    // TPC only tag
    if (posTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitMap.set(selPosNotTPCOnly);
    if (negTrackExtra.detectorMap() != o2::aod::track::TPC)
      bitMap.set(selNegNotTPCOnly);

    // proper lifetime
    if (ctauLambda < nCtauCutV0->get("lifetimecutLambda") * ctaulambdaPDG)
      bitMap.set(selLambdaCTau);
    if (ctauK0short < nCtauCutV0->get("lifetimecutK0S") * ctauk0shortPDG)
      bitMap.set(selK0ShortCTau);

    // armenteros
    if (v0.qtarm() * v0cuts.armPodCut > std::fabs(v0.alpha()) || v0cuts.armPodCut < 1e-4)
      bitMap.set(selK0ShortArmenteros);

    return bitMap;
  }

  template <typename TCasc, typename TCollision>
  void analyseCascCandidate(TCasc const& casc, TCollision const& coll, int const& gap, std::bitset<kSelNum> const& selMap)
  {
    // Access daughter tracks
    auto posTrackExtra = casc.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = casc.template negTrackExtra_as<DauTracks>();
    auto bachTrackExtra = casc.template bachTrackExtra_as<DauTracks>();

    if (doPlainTopoQA) {
      histos.fill(HIST("generalQA/hPt"), casc.pt());
      histos.fill(HIST("generalQA/hCascCosPA"), casc.pt(), casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()));
      histos.fill(HIST("generalQA/hDCACascDaughters"), casc.pt(), casc.dcacascdaughters());
      histos.fill(HIST("generalQA/hCascRadius"), casc.pt(), casc.cascradius());
      histos.fill(HIST("generalQA/hMesonDCAToPV"), casc.pt(), casc.dcanegtopv());
      histos.fill(HIST("generalQA/hBaryonDCAToPV"), casc.pt(), casc.dcapostopv());
      histos.fill(HIST("generalQA/hBachDCAToPV"), casc.pt(), casc.dcabachtopv());
      histos.fill(HIST("generalQA/hV0CosPA"), casc.pt(), casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()));
      histos.fill(HIST("generalQA/hV0Radius"), casc.pt(), casc.v0radius());
      histos.fill(HIST("generalQA/hDCAV0Daughters"), casc.pt(), casc.dcaV0daughters());
      histos.fill(HIST("generalQA/hDCAV0ToPV"), casc.pt(), std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())));
      histos.fill(HIST("generalQA/hMassLambdaDau"), casc.pt(), casc.mLambda());
      histos.fill(HIST("generalQA/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
      histos.fill(HIST("generalQA/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      histos.fill(HIST("generalQA/h2dBachITSvsTPCpts"), bachTrackExtra.tpcCrossedRows(), bachTrackExtra.itsNCls());
    }

    // Xi
    if (verifyMask(selMap, maskSelectionXi) && analyseXi) {
      fillHistogramsCasc<3>(casc, coll, gap);
    }

    // Anti-Xi
    if (verifyMask(selMap, maskSelectionAntiXi) && analyseAntiXi) {
      fillHistogramsCasc<4>(casc, coll, gap);
    }

    // Omega
    if (verifyMask(selMap, maskSelectionOmega) && analyseOmega) {
      fillHistogramsCasc<5>(casc, coll, gap);
    }

    // Anti-Omega
    if (verifyMask(selMap, maskSelectionAntiOmega) && analyseAntiOmega) {
      fillHistogramsCasc<6>(casc, coll, gap);
    }
  }

  template <typename TV0>
  void computeV0MCAssociation(const TV0& v0, std::bitset<kSelNum>& bitMap)
  {
    const int pdgPos = v0.pdgCodePositive();
    const int pdgNeg = v0.pdgCodeNegative();
    const int pdgV0 = v0.pdgCode();
    const bool isPhysPrim = v0.isPhysicalPrimary();

    const bool isPositiveProton = (pdgPos == PDG_t::kProton);
    const bool isPositivePion = (pdgPos == PDG_t::kPiPlus) || (doTreatPiToMuon && pdgPos == PDG_t::kMuonPlus);
    const bool isNegativeProton = (pdgNeg == kProtonBar);
    const bool isNegativePion = (pdgNeg == PDG_t::kPiMinus) || (doTreatPiToMuon && pdgNeg == PDG_t::kMuonMinus);

    switch (pdgV0) {
      case PDG_t::kK0Short: // K0Short
        if (isPositivePion && isNegativePion) {
          bitMap.set(selConsiderK0Short);
          if (isPhysPrim)
            bitMap.set(selPhysPrimK0Short);
        }
        break;
      case PDG_t::kLambda0: // Lambda
        if (isPositiveProton && isNegativePion) {
          bitMap.set(selConsiderLambda);
          if (isPhysPrim)
            bitMap.set(selPhysPrimLambda);
        }
        break;
      case PDG_t::kLambda0Bar: // AntiLambda
        if (isPositivePion && isNegativeProton) {
          bitMap.set(selConsiderAntiLambda);
          if (isPhysPrim)
            bitMap.set(selPhysPrimAntiLambda);
        }
        break;
    }
  }

  template <typename TCasc>
  void computeCascadeMCAssociation(const TCasc& casc, std::bitset<kSelNum>& bitMap)
  {
    const int pdgPos = casc.pdgCodePositive();
    const int pdgNeg = casc.pdgCodeNegative();
    const int pdgBach = casc.pdgCodeBachelor();
    const int pdgCasc = casc.pdgCode();
    const bool isPhysPrim = casc.isPhysicalPrimary();

    const bool isPositiveProton = (pdgPos == PDG_t::kProton);

    const bool isPositivePion = (pdgPos == PDG_t::kPiPlus);
    const bool isBachelorPositivePion = (pdgBach == PDG_t::kPiPlus);

    const bool isNegativeProton = (pdgNeg == kProtonBar);

    const bool isNegativePion = (pdgNeg == PDG_t::kPiMinus);
    const bool isBachelorNegativePion = (pdgBach == PDG_t::kPiMinus);

    const bool isBachelorPositiveKaon = (pdgBach == PDG_t::kKPlus);
    const bool isBachelorNegativeKaon = (pdgBach == PDG_t::kKMinus);

    switch (pdgCasc) {
      case PDG_t::kXiMinus: // Xi
        if (isPositiveProton && isNegativePion && isBachelorNegativePion) {
          bitMap.set(selConsiderXi);
          if (isPhysPrim)
            bitMap.set(selPhysPrimXi);
        }
        break;
      case PDG_t::kXiPlusBar: // Anti-Xi
        if (isNegativeProton && isPositivePion && isBachelorPositivePion) {
          bitMap.set(selConsiderAntiXi);
          if (isPhysPrim)
            bitMap.set(selPhysPrimAntiXi);
        }
        break;
      case PDG_t::kOmegaMinus: // Omega
        if (isPositiveProton && isNegativePion && isBachelorNegativeKaon) {
          bitMap.set(selConsiderOmega);
          if (isPhysPrim)
            bitMap.set(selPhysPrimOmega);
        }
        break;
      case PDG_t::kOmegaPlusBar: // Anti-Omega
        if (isNegativeProton && isPositivePion && isBachelorPositiveKaon) {
          bitMap.set(selConsiderAntiOmega);
          if (isPhysPrim)
            bitMap.set(selPhysPrimAntiOmega);
        }
        break;
    }
  }

  template <typename TV0, typename TCollision>
  void analyseV0Candidate(TV0 const& v0, TCollision const& coll, int const& gap, std::bitset<kSelNum> const& selMap)
  {
    auto posTrackExtra = v0.template posTrackExtra_as<DauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<DauTracks>();

    // QA plots
    if (doPlainTopoQA) {
      histos.fill(HIST("generalQA/hPt"), v0.pt());
      histos.fill(HIST("generalQA/hPosDCAToPV"), v0.dcapostopv());
      histos.fill(HIST("generalQA/hNegDCAToPV"), v0.dcanegtopv());
      histos.fill(HIST("generalQA/hDCADaughters"), v0.dcaV0daughters());
      histos.fill(HIST("generalQA/hPointingAngle"), std::acos(v0.v0cosPA()));
      histos.fill(HIST("generalQA/hCosPA"), v0.v0cosPA());
      histos.fill(HIST("generalQA/hV0Radius"), v0.v0radius());
      histos.fill(HIST("generalQA/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
      histos.fill(HIST("generalQA/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
    }

    histos.fill(HIST("generalQA/h2dArmenterosAll"), v0.alpha(), v0.qtarm());

    // K0s
    if (verifyMask(selMap, maskSelectionK0Short) && analyseK0Short) {
      fillHistogramsV0<0>(v0, coll, gap);
    }

    // Lambda
    if (verifyMask(selMap, maskSelectionLambda) && analyseLambda) {
      fillHistogramsV0<1>(v0, coll, gap);
    }

    // Anti-Lambda
    if (verifyMask(selMap, maskSelectionAntiLambda) && analyseAntiLambda) {
      fillHistogramsV0<2>(v0, coll, gap);
    }
  }

  PresliceUnsorted<StraCollisonsFullMC> perMcCollision = aod::v0data::straMCCollisionId;
  PresliceUnsorted<NeutronsMC> neutronsPerMcCollision = aod::zdcneutrons::straMCCollisionId;

  std::vector<int> getListOfRecoCollIds(StraMCCollisionsFull const& mcCollisions,
                                        StraCollisonsFullMC const& collisions,
                                        NeutronsMC const& neutrons)
  {
    std::vector<int> listBestCollisionIds(mcCollisions.size(), -1);

    for (auto const& mcCollision : mcCollisions) {
      if (std::find(generatorIds->begin(), generatorIds->end(), mcCollision.generatorsID()) == generatorIds->end()) {
        continue;
      }

      // Group collisions and neutrons by MC collision index
      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto groupedNeutrons = neutrons.sliceBy(neutronsPerMcCollision, mcCollision.globalIndex());
      // Find the collision with the biggest nbr of PV contributors
      // Follows what was done here: https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/mcCollsExtra.cxx#L93
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      for (auto const& collision : groupedCollisions) {
        if (!acceptEvent(collision, false)) {
          continue;
        }

        int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
        if (checkNeutronsInMC) {
          for (const auto& neutron : groupedNeutrons) {
            if (selGapSide < -0.5)
              break;

            const float eta = neutron.eta();
            switch (selGapSide) {
              case 0: // SGA
                if (eta > neutronEtaCut)
                  selGapSide = -1;
                break;
              case 1: // SGC
                if (eta < -neutronEtaCut)
                  selGapSide = -1;
                break;
              case 2: // DG
                if (eta > neutronEtaCut)
                  selGapSide = 1;
                else if (eta < -neutronEtaCut)
                  selGapSide = 0;
                break;
            }
          }
        }

        if (evSels.studyUPConly && (selGapSide != static_cast<int>(upcCuts.genGapSide)))
          continue;

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          bestCollisionIndex = collision.globalIndex();
        }
      }
      listBestCollisionIds[mcCollision.globalIndex()] = bestCollisionIndex;
    }

    return listBestCollisionIds;
  }

  void fillGenMCHistogramsQA(StraMCCollisionsFull const& mcCollisions,
                             StraCollisonsFullMC const& collisions,
                             NeutronsMC const& neutrons)
  {
    for (auto const& mcCollision : mcCollisions) {
      // LOGF(info, "Generator ID is %i", mcCollision.generatorsID());
      histos.fill(HIST("eventQA/mc/hGeneratorsId"), mcCollision.generatorsID());

      if (std::find(generatorIds->begin(), generatorIds->end(), mcCollision.generatorsID()) == generatorIds->end()) {
        continue;
      }

      histos.fill(HIST("eventQA/mc/hSelGeneratorsId"), mcCollision.generatorsID());

      histos.fill(HIST("eventQA/mc/hEventSelectionMC"), 0.0, mcCollision.multMCNParticlesEta08(), mcCollision.generatorsID());

      if (std::abs(mcCollision.posZ()) > maxZVtxPosition)
        continue;

      histos.fill(HIST("eventQA/mc/hEventSelectionMC"), 1.0, mcCollision.multMCNParticlesEta08(), mcCollision.generatorsID());

      // Group collisions and neutrons by MC collision index
      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto groupedNeutrons = neutrons.sliceBy(neutronsPerMcCollision, mcCollision.globalIndex());

      bool atLeastOne = false;
      float centrality = -1.f;
      float ft0ampl = -1.f;
      int nCollisions = 0;
      int biggestNContribs = -1;
      int nTracksGlobal = -1;

      // Find the max contributors and count accepted collisions
      for (auto const& collision : groupedCollisions) {
        if (!acceptEvent(collision, false)) {
          continue;
        }

        int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
        if (checkNeutronsInMC) {
          for (const auto& neutron : groupedNeutrons) {
            if (selGapSide < -0.5)
              break;

            const float eta = neutron.eta();
            switch (selGapSide) {
              case 0: // SGA
                if (eta > neutronEtaCut)
                  selGapSide = -1;
                break;
              case 1: // SGC
                if (eta < -neutronEtaCut)
                  selGapSide = -1;
                break;
              case 2: // DG
                if (eta > neutronEtaCut)
                  selGapSide = 1;
                else if (eta < -neutronEtaCut)
                  selGapSide = 0;
                break;
            }
          }
        }

        if (evSels.studyUPConly && (selGapSide != static_cast<int>(upcCuts.genGapSide)))
          continue;

        ++nCollisions;
        atLeastOne = true;

        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          if (static_cast<int>(upcCuts.genGapSide) == 0) {
            ft0ampl = collision.totalFT0AmplitudeC();
            centrality = collision.centFT0C();
          } else if (static_cast<int>(upcCuts.genGapSide) == 1) {
            ft0ampl = collision.totalFT0AmplitudeA();
            centrality = collision.centFT0A();
          }
          nTracksGlobal = collision.multNTracksGlobal();
        }
      }

      // Fill histograms
      histos.fill(HIST("eventQA/mc/hTracksGlobalVsNcoll_beforeEvSel"), nTracksGlobal, groupedCollisions.size());
      histos.fill(HIST("eventQA/mc/hTracksGlobalVsNcoll_afterEvSel"), nTracksGlobal, nCollisions);
      histos.fill(HIST("eventQA/mc/hTracksGlobalvsMCNParticlesEta08gen"), nTracksGlobal, mcCollision.multMCNParticlesEta08());
      histos.fill(HIST("eventQA/mc/hTracksGlobalVsPVzMC"), nTracksGlobal, mcCollision.posZ());
      histos.fill(HIST("eventQA/mc/hEventPVzMC"), mcCollision.posZ());

      if (atLeastOne) {
        histos.fill(HIST("eventQA/mc/hEventSelectionMC"), 2.0, mcCollision.multMCNParticlesEta08(), mcCollision.generatorsID());
        histos.fill(HIST("eventQA/mc/hGenEventCentrality"), centrality);
        histos.fill(HIST("eventQA/mc/hGenEventFT0ampl"), ft0ampl);
      }
    }
  }

  template <typename TCollision, typename TV0>
  void fillFeeddownMatrix(TCollision const& collision, TV0 const& v0, std::bitset<kSelNum> const& selMap)
  {
    if (!v0.has_motherMCPart()) {
      return;
    }

    const auto v0mother = v0.template motherMCPart_as<aod::MotherMCParts>();
    if (v0mother.size() < 1) {
      return;
    }

    const float rapidityXi = RecoDecay::y(std::array{v0mother.px(), v0mother.py(), v0mother.pz()}, o2::constants::physics::MassXiMinus);
    if (std::fabs(rapidityXi) > 0.5f) {
      return;
    }

    const float mult = collision.multNTracksGlobal();
    const float v0pt = v0.pt();
    const float motherPt = std::hypot(v0mother.px(), v0mother.py());

    if (analyseLambda && verifyMask(selMap, secondaryMaskSelectionLambda) &&
        (v0mother.pdgCode() == PDG_t::kXiMinus) && v0mother.isPhysicalPrimary()) {
      histos.fill(HIST(kParticlenames[1]) + HIST("/h3dLambdaFeeddown"), mult, v0pt, motherPt);
    }

    if (analyseAntiLambda && verifyMask(selMap, secondaryMaskSelectionAntiLambda) &&
        (v0mother.pdgCode() == PDG_t::kXiPlusBar) && v0mother.isPhysicalPrimary()) {
      histos.fill(HIST(kParticlenames[2]) + HIST("/h3dAntiLambdaFeeddown"), mult, v0pt, motherPt);
    }
  }

  void processV0s(StraCollisonsFull const& collisions, V0Candidates const& fullV0s, DauTracks const&)
  {
    v0sGrouped.clear();
    v0sGrouped.resize(collisions.size());
    for (const auto& v0 : fullV0s) {
      v0sGrouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }
    for (const auto& collision : collisions) {
      if (!acceptEvent(collision, true)) {
        continue;
      } // event is accepted

      histos.fill(HIST("eventQA/hRawGapSide"), collision.gapSide());

      int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
      if (evSels.studyUPConly && (selGapSide < -0.5))
        continue;

      fillHistogramsQA(collision, selGapSide);

      std::size_t nV0sThisColl = v0sGrouped[collision.globalIndex()].size();

      for (std::size_t i = 0; i < nV0sThisColl; i++) {
        auto v0 = fullV0s.rawIteratorAt(v0sGrouped[collision.globalIndex()][i]);
        if ((v0.v0Type() != v0cuts.v0TypeSelection) && (v0cuts.v0TypeSelection > 0))
          continue; // skip V0s that are not standard

        std::bitset<kSelNum> selMap = computeBitmapV0(v0, collision);

        // consider all species for the candidate
        setBits(selMap, {selConsiderK0Short, selConsiderLambda, selConsiderAntiLambda,
                         selPhysPrimK0Short, selPhysPrimLambda, selPhysPrimAntiLambda});

        analyseV0Candidate(v0, collision, selGapSide, selMap);
      } // end v0 loop
    }
  }

  void processV0sMC(StraCollisonsFullMC const& collisions,
                    V0CandidatesMC const& fullV0s,
                    DauTracks const&,
                    aod::MotherMCParts const&,
                    StraMCCollisionsFull const&,
                    V0MCCoresFull const&,
                    NeutronsMC const& neutrons)
  {
    v0sGrouped.clear();
    v0sGrouped.resize(collisions.size());
    for (const auto& v0 : fullV0s) {
      v0sGrouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }

    for (const auto& collision : collisions) {
      if (!collision.has_straMCCollision()) {
        histos.fill(HIST("eventQA/mc/hFakeEvents"), 0); // no assoc. MC collisions
        continue;
      }

      const auto& mcCollision = collision.straMCCollision_as<StraMCCollisionsFull>(); // take gen. collision associated to the rec. collision

      if (std::find(generatorIds->begin(), generatorIds->end(), mcCollision.generatorsID()) == generatorIds->end()) {
        continue;
      }

      if (!acceptEvent(collision, true)) {
        continue;
      } // event is accepted

      histos.fill(HIST("eventQA/hRawGapSide"), collision.gapSide());

      int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
      int selGapSideNoNeutrons = selGapSide;

      auto groupedNeutrons = neutrons.sliceBy(neutronsPerMcCollision, mcCollision.globalIndex());
      if (checkNeutronsInMC) {
        for (const auto& neutron : groupedNeutrons) {
          if (selGapSide < -0.5)
            break;

          const float eta = neutron.eta();
          switch (selGapSide) {
            case 0: // SGA
              if (eta > neutronEtaCut)
                selGapSide = -1;
              break;
            case 1: // SGC
              if (eta < -neutronEtaCut)
                selGapSide = -1;
              break;
            case 2: // DG
              if (eta > neutronEtaCut)
                selGapSide = 1;
              else if (eta < -neutronEtaCut)
                selGapSide = 0;
              break;
          }
        }
      }

      if (evSels.studyUPConly && (selGapSide < -0.5))
        continue;

      histos.fill(HIST("eventQA/hSelGapSideNoNeutrons"), selGapSideNoNeutrons);
      fillHistogramsQA(collision, selGapSide);

      histos.fill(HIST("eventQA/mc/hNTracksGlobalvsMCNParticlesEta08rec"), collision.multNTracksGlobal(), mcCollision.multMCNParticlesEta08());
      histos.fill(HIST("eventQA/mc/hNTracksPVeta1vsMCNParticlesEta10rec"), collision.multNTracksPVeta1(), mcCollision.multMCNParticlesEta10());
      histos.fill(HIST("eventQA/mc/hNTracksGlobalvstotalMultMCParticles"), collision.multNTracksGlobal(), mcCollision.totalMultMCParticles());
      histos.fill(HIST("eventQA/mc/hNTracksPVeta1vstotalMultMCParticles"), collision.multNTracksPVeta1(), mcCollision.totalMultMCParticles());

      std::size_t nV0sThisColl = v0sGrouped[collision.globalIndex()].size();

      for (std::size_t i = 0; i < nV0sThisColl; i++) {
        auto v0 = fullV0s.rawIteratorAt(v0sGrouped[collision.globalIndex()][i]);
        if ((v0.v0Type() != v0cuts.v0TypeSelection) && (v0cuts.v0TypeSelection > 0))
          continue; // skip V0s that are not standard

        std::bitset<kSelNum> selMap = computeBitmapV0(v0, collision);

        if (doMCAssociation) {
          if (v0.has_v0MCCore()) {
            const auto& v0MC = v0.v0MCCore_as<V0MCCoresFull>();
            computeV0MCAssociation(v0MC, selMap);
            if (calculateFeeddownMatrix) {
              fillFeeddownMatrix(collision, v0, selMap);
            }
          }
        } else {
          // consider all species for the candidate
          setBits(selMap, {selConsiderK0Short, selConsiderLambda, selConsiderAntiLambda,
                           selPhysPrimK0Short, selPhysPrimLambda, selPhysPrimAntiLambda});
        }

        analyseV0Candidate(v0, collision, selGapSide, selMap);
      } // end v0 loop
    }
  }

  void processCascades(StraCollisonsFull const& collisions,
                       CascadeCandidates const& fullCascades,
                       DauTracks const&)
  {
    cascadesGrouped.clear();
    cascadesGrouped.resize(collisions.size());
    for (const auto& cascade : fullCascades) {
      cascadesGrouped[cascade.straCollisionId()].push_back(cascade.globalIndex());
    }

    for (const auto& collision : collisions) {
      if (!acceptEvent(collision, true)) {
        continue;
      } // event is accepted

      histos.fill(HIST("eventQA/hRawGapSide"), collision.gapSide());

      int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
      if (evSels.studyUPConly && (selGapSide < -0.5))
        continue;

      fillHistogramsQA(collision, selGapSide);

      std::size_t nCascadesThisColl = cascadesGrouped[collision.globalIndex()].size();

      for (std::size_t i = 0; i < nCascadesThisColl; i++) {
        auto casc = fullCascades.rawIteratorAt(cascadesGrouped[collision.globalIndex()][i]);
        std::bitset<kSelNum> selMap = computeBitmapCascade(casc, collision);
        // the candidate may belong to any particle species
        setBits(selMap, {selConsiderXi, selConsiderAntiXi, selConsiderOmega, selConsiderAntiOmega,
                         selPhysPrimXi, selPhysPrimAntiXi, selPhysPrimOmega, selPhysPrimAntiOmega});

        analyseCascCandidate(casc, collision, selGapSide, selMap);
      } // end casc loop
    }
  }

  void processCascadesMC(StraCollisonsFullMC const& collisions,
                         CascadeCandidatesMC const& fullCascades,
                         DauTracks const&,
                         aod::MotherMCParts const&,
                         StraMCCollisionsFull const&,
                         CascMCCoresFull const&,
                         NeutronsMC const& neutrons)
  {
    cascadesGrouped.clear();
    cascadesGrouped.resize(collisions.size());
    for (const auto& cascade : fullCascades) {
      cascadesGrouped[cascade.straCollisionId()].push_back(cascade.globalIndex());
    }

    for (const auto& collision : collisions) {
      if (!collision.has_straMCCollision()) {
        histos.fill(HIST("eventQA/mc/hFakeEvents"), 0); // no assoc. MC collisions
        continue;
      }

      const auto& mcCollision = collision.straMCCollision_as<StraMCCollisionsFull>(); // take gen. collision associated to the rec. collision

      if (std::find(generatorIds->begin(), generatorIds->end(), mcCollision.generatorsID()) == generatorIds->end()) {
        continue;
      }

      if (!acceptEvent(collision, true)) {
        continue;
      } // event is accepted

      histos.fill(HIST("eventQA/hRawGapSide"), collision.gapSide());

      int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
      int selGapSideNoNeutrons = selGapSide;

      auto groupedNeutrons = neutrons.sliceBy(neutronsPerMcCollision, mcCollision.globalIndex());
      if (checkNeutronsInMC) {
        for (const auto& neutron : groupedNeutrons) {
          if (selGapSide < -0.5)
            break;

          const float eta = neutron.eta();
          switch (selGapSide) {
            case 0: // SGA
              if (eta > neutronEtaCut)
                selGapSide = -1;
              break;
            case 1: // SGC
              if (eta < -neutronEtaCut)
                selGapSide = -1;
              break;
            case 2: // DG
              if (eta > neutronEtaCut)
                selGapSide = 1;
              else if (eta < -neutronEtaCut)
                selGapSide = 0;
              break;
          }
        }
      }

      if (evSels.studyUPConly && (selGapSide < -0.5))
        continue;

      histos.fill(HIST("eventQA/hSelGapSideNoNeutrons"), selGapSideNoNeutrons);
      fillHistogramsQA(collision, selGapSide);

      histos.fill(HIST("eventQA/mc/hNTracksGlobalvsMCNParticlesEta08rec"), collision.multNTracksGlobal(), mcCollision.multMCNParticlesEta08());
      histos.fill(HIST("eventQA/mc/hNTracksPVeta1vsMCNParticlesEta10rec"), collision.multNTracksPVeta1(), mcCollision.multMCNParticlesEta10());
      histos.fill(HIST("eventQA/mc/hNTracksGlobalvstotalMultMCParticles"), collision.multNTracksGlobal(), mcCollision.totalMultMCParticles());
      histos.fill(HIST("eventQA/mc/hNTracksPVeta1vstotalMultMCParticles"), collision.multNTracksPVeta1(), mcCollision.totalMultMCParticles());

      std::size_t nCascadesThisColl = cascadesGrouped[collision.globalIndex()].size();

      for (std::size_t i = 0; i < nCascadesThisColl; i++) {
        auto casc = fullCascades.rawIteratorAt(cascadesGrouped[collision.globalIndex()][i]);
        std::bitset<kSelNum> selMap = computeBitmapCascade(casc, collision);

        if (doMCAssociation) {
          if (casc.has_cascMCCore()) {
            const auto& cascMC = casc.cascMCCore_as<CascMCCoresFull>();
            computeCascadeMCAssociation(cascMC, selMap);
          }
        } else {
          // the candidate may belong to any particle species
          setBits(selMap, {selConsiderXi, selConsiderAntiXi, selConsiderOmega, selConsiderAntiOmega,
                           selPhysPrimXi, selPhysPrimAntiXi, selPhysPrimOmega, selPhysPrimAntiOmega});
        }

        analyseCascCandidate(casc, collision, selGapSide, selMap);
      } // end casc loop
    }
  }

  void processGenerated(StraMCCollisionsFull const& mcCollisions,
                        V0MCCoresFull const& V0MCCores,
                        CascMCCoresFull const& CascMCCores,
                        StraCollisonsFullMC const& collisions,
                        NeutronsMC const& neutrons)
  {
    fillGenMCHistogramsQA(mcCollisions, collisions, neutrons);
    std::vector<int> listBestCollisionIds = getListOfRecoCollIds(mcCollisions, collisions, neutrons);
    // V0 start
    for (auto const& v0MC : V0MCCores) {
      // Consider only primaries
      if (!v0MC.has_straMCCollision() || !v0MC.isPhysicalPrimary())
        continue;

      // Kinematics (|y| < rapidityCut)
      float pTmc = v0MC.ptMC();
      float ymc = 1e3;
      if (v0MC.pdgCode() == PDG_t::kK0Short)
        ymc = v0MC.rapidityMC(0);
      else if ((v0MC.pdgCode() == PDG_t::kLambda0) || (v0MC.pdgCode() == PDG_t::kLambda0Bar))
        ymc = v0MC.rapidityMC(1);
      if (std::abs(ymc) > rapidityCut)
        continue;

      const auto& mcCollision = v0MC.straMCCollision_as<StraMCCollisionsFull>(); // take gen. collision

      if (std::abs(mcCollision.posZ()) > maxZVtxPosition)
        continue;

      // Collision is of the proccess of interest
      if (std::find(generatorIds->begin(), generatorIds->end(), mcCollision.generatorsID()) == generatorIds->end()) {
        continue;
      }

      // float centrality = -1.f;
      float ft0ampl = -1.f;
      int nTracksGlobal = -1;

      if (listBestCollisionIds[mcCollision.globalIndex()] > -1) {
        auto collision = collisions.iteratorAt(listBestCollisionIds[mcCollision.globalIndex()]);
        // centrality = collision.centFT0C();
        if (static_cast<int>(upcCuts.genGapSide) == 0) {
          ft0ampl = collision.totalFT0AmplitudeC();
        } else if (static_cast<int>(upcCuts.genGapSide) == 1) {
          ft0ampl = collision.totalFT0AmplitudeA();
        }
        nTracksGlobal = collision.multNTracksGlobal();
      }

      const int pdgPos = v0MC.pdgCodePositive();
      const int pdgNeg = v0MC.pdgCodeNegative();
      const int pdgV0 = v0MC.pdgCode();

      const bool isPositiveProton = (pdgPos == PDG_t::kProton);
      const bool isPositivePion = (pdgPos == PDG_t::kPiPlus) || (doTreatPiToMuon && pdgPos == PDG_t::kMuonPlus);
      const bool isNegativeProton = (pdgNeg == kProtonBar);
      const bool isNegativePion = (pdgNeg == PDG_t::kPiMinus) || (doTreatPiToMuon && pdgNeg == PDG_t::kMuonMinus);

      // Fill histograms
      if ((pdgV0 == PDG_t::kK0Short) && isPositivePion && isNegativePion) {
        histos.fill(HIST(kParticlenames[0]) + HIST("/mc/h7dGen"), ft0ampl, nTracksGlobal, mcCollision.multMCNParticlesEta08(), pTmc, static_cast<int>(upcCuts.genGapSide), ymc, mcCollision.generatorsID());
      }
      if ((pdgV0 == PDG_t::kLambda0) && isPositiveProton && isNegativePion) {
        histos.fill(HIST(kParticlenames[1]) + HIST("/mc/h7dGen"), ft0ampl, nTracksGlobal, mcCollision.multMCNParticlesEta08(), pTmc, static_cast<int>(upcCuts.genGapSide), ymc, mcCollision.generatorsID());
      }
      if ((pdgV0 == PDG_t::kLambda0Bar) && isPositivePion && isNegativeProton) {
        histos.fill(HIST(kParticlenames[2]) + HIST("/mc/h7dGen"), ft0ampl, nTracksGlobal, mcCollision.multMCNParticlesEta08(), pTmc, static_cast<int>(upcCuts.genGapSide), ymc, mcCollision.generatorsID());
      }
    } // V0 end

    // Cascade start
    for (auto const& cascMC : CascMCCores) {
      // Consider only primaries
      if (!cascMC.has_straMCCollision() || !cascMC.isPhysicalPrimary())
        continue;
      // Kinematics (|y| < rapidityCut)
      float pTmc = cascMC.ptMC();
      float ymc = 1e3;
      if ((cascMC.pdgCode() == PDG_t::kXiMinus) || (cascMC.pdgCode() == PDG_t::kXiPlusBar)) {
        ymc = cascMC.rapidityMC(0);
      } else if ((cascMC.pdgCode() == PDG_t::kOmegaMinus) || (cascMC.pdgCode() == PDG_t::kOmegaPlusBar)) {
        ymc = cascMC.rapidityMC(2);
      }
      if (std::abs(ymc) > rapidityCut)
        continue;

      const auto& mcCollision = cascMC.straMCCollision_as<StraMCCollisionsFull>(); // take gen. collision
      if (std::abs(mcCollision.posZ()) > maxZVtxPosition)
        continue;

      // float centrality = -1.f;
      float ft0ampl = -1.f;
      int nTracksGlobal = -1;

      if (listBestCollisionIds[mcCollision.globalIndex()] > -1) {
        auto collision = collisions.iteratorAt(listBestCollisionIds[mcCollision.globalIndex()]);
        // centrality = collision.centFT0C();
        if (static_cast<int>(upcCuts.genGapSide) == 0) {
          ft0ampl = collision.totalFT0AmplitudeC();
        } else if (static_cast<int>(upcCuts.genGapSide) == 1) {
          ft0ampl = collision.totalFT0AmplitudeA();
        }
      }

      const int pdgPos = cascMC.pdgCodePositive();
      const int pdgNeg = cascMC.pdgCodeNegative();
      const int pdgBach = cascMC.pdgCodeBachelor();
      const int pdgCasc = cascMC.pdgCode();

      const bool isPositiveProton = (pdgPos == PDG_t::kProton);

      const bool isPositivePion = (pdgPos == PDG_t::kPiPlus);
      const bool isBachelorPositivePion = (pdgBach == PDG_t::kPiPlus);

      const bool isNegativeProton = (pdgNeg == kProtonBar);

      const bool isNegativePion = (pdgNeg == PDG_t::kPiMinus);
      const bool isBachelorNegativePion = (pdgBach == PDG_t::kPiMinus);

      const bool isBachelorPositiveKaon = (pdgBach == PDG_t::kKPlus);
      const bool isBachelorNegativeKaon = (pdgBach == PDG_t::kKMinus);

      // Fill histograms
      if ((pdgCasc == PDG_t::kXiMinus) && isPositiveProton && isNegativePion && isBachelorNegativePion) {
        histos.fill(HIST(kParticlenames[3]) + HIST("/mc/h7dGen"), ft0ampl, nTracksGlobal, mcCollision.multMCNParticlesEta08(), pTmc, static_cast<int>(upcCuts.genGapSide), ymc, mcCollision.generatorsID());
      }
      if ((pdgCasc == PDG_t::kXiPlusBar) && isNegativeProton && isPositivePion && isBachelorPositivePion) {
        histos.fill(HIST(kParticlenames[4]) + HIST("/mc/h7dGen"), ft0ampl, nTracksGlobal, mcCollision.multMCNParticlesEta08(), pTmc, static_cast<int>(upcCuts.genGapSide), ymc, mcCollision.generatorsID());
      }
      if ((pdgCasc == PDG_t::kOmegaMinus) && isPositiveProton && isNegativePion && isBachelorNegativeKaon) {
        histos.fill(HIST(kParticlenames[5]) + HIST("/mc/h7dGen"), ft0ampl, nTracksGlobal, mcCollision.multMCNParticlesEta08(), pTmc, static_cast<int>(upcCuts.genGapSide), ymc, mcCollision.generatorsID());
      }
      if ((pdgCasc == PDG_t::kOmegaPlusBar) && isNegativeProton && isPositivePion && isBachelorPositiveKaon) {
        histos.fill(HIST(kParticlenames[6]) + HIST("/mc/h7dGen"), ft0ampl, nTracksGlobal, mcCollision.multMCNParticlesEta08(), pTmc, static_cast<int>(upcCuts.genGapSide), ymc, mcCollision.generatorsID());
      }
    } // Cascade end
  }

  PROCESS_SWITCH(Derivedupcanalysis, processV0s, "Process V0s", true);
  PROCESS_SWITCH(Derivedupcanalysis, processV0sMC, "Process V0s MC", false);
  PROCESS_SWITCH(Derivedupcanalysis, processCascades, "Process Cascades", false);
  PROCESS_SWITCH(Derivedupcanalysis, processCascadesMC, "Process Cascades MC", false);
  PROCESS_SWITCH(Derivedupcanalysis, processGenerated, "Process Generated Level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Derivedupcanalysis>(cfgc)};
}
