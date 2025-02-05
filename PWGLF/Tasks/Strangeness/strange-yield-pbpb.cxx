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
// Strangeness in UPC analysis task
// ================
// This code is meant to be run over derived data.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    roman.nepeivoda@cern.ch
//

#include <bitset>
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/StaticFor.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGLF/Utils/strangenessMasks.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using dauMCTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackMCIds, aod::DauTrackTPCPIDs>;

using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
using v0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0MCCores, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0MCCollRefs>;

using cascadeCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascTOFPIDs, aod::CascTOFNSigmas>;

using straCollisonFull = soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator;

struct strangeYieldPbPb {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master analysis switches
  Configurable<bool> analyseK0Short{"analyseK0Short", true, "process K0Short-like candidates"};
  Configurable<bool> analyseLambda{"analyseLambda", true, "process Lambda-like candidates"};
  Configurable<bool> analyseAntiLambda{"analyseAntiLambda", true, "process AntiLambda-like candidates"};
  Configurable<bool> analyseXi{"analyseXi", true, "process Xi-like candidates"};
  Configurable<bool> analyseAntiXi{"analyseAntiXi", true, "process AntiXi-like candidates"};
  Configurable<bool> analyseOmega{"analyseOmega", true, "process Omega-like candidates"};
  Configurable<bool> analyseAntiOmega{"analyseAntiOmega", true, "process AntiOmega-like candidates"};

  // Event selections
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
  Configurable<bool> requireIsTriggerTVX{"requireIsTriggerTVX", true, "require coincidence in FT0A and FT0C"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
  Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
  Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
  Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
  Configurable<bool> studyUPConly{"studyUPConly", false, "is UPC-only analysis"};
  Configurable<bool> useUPCflag{"useUPCflag", false, "select UPC flagged events"};

  Configurable<bool> verbose{"verbose", false, "additional printouts"};

  // Acceptance selections
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
  Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

  // Standard V0 topological criteria
  struct : ConfigurableGroup {
    Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.5, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>fabs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};
    Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};
  } v0cuts;
  static constexpr float lifetimeCutsV0[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecutV0{"lifetimecutV0", {lifetimeCutsV0[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecutV0"};

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
  static constexpr float nCtauCutsCasc[1][2] = {{6., 6.}};
  Configurable<LabeledArray<float>> nCtauCutCasc{"nCtauCutCasc", {nCtauCutsCasc[0], 2, {"lifetimecutXi", "lifetimecutOmega"}}, "nCtauCutCasc"};

  // UPC selections
  SGSelector sgSelector;
  struct : ConfigurableGroup {
    Configurable<float> FV0cut{"FV0cut", 100., "FV0A threshold"};
    Configurable<float> FT0Acut{"FT0Acut", 200., "FT0A threshold"};
    Configurable<float> FT0Ccut{"FT0Ccut", 100., "FT0C threshold"};
    Configurable<float> ZDCcut{"ZDCcut", 10., "ZDC threshold"};
    // Configurable<float> gapSel{"gapSel", 2, "Gap selection"};
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
    Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 1e+6, "TpcPidNsigmaCut"};
    Configurable<float> TofPidNsigmaCutLaPr{"TofPidNsigmaCutLaPr", 1e+6, "TofPidNsigmaCutLaPr"};
    Configurable<float> TofPidNsigmaCutLaPi{"TofPidNsigmaCutLaPi", 1e+6, "TofPidNsigmaCutLaPi"};
    Configurable<float> TofPidNsigmaCutK0Pi{"TofPidNsigmaCutK0Pi", 1e+6, "TofPidNsigmaCutK0Pi"};

    Configurable<float> TofPidNsigmaCutXiPi{"TofPidNsigmaCutXiPi", 1e+6, "TofPidNsigmaCutXiPi"};
    Configurable<float> TofPidNsigmaCutOmegaKaon{"TofPidNsigmaCutOmegaKaon", 1e+6, "TofPidNsigmaCutOmegaKaon"};

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
    ConfigurableAxis axisFV0Aampl{"axisFV0Aampl", {100, 0.0f, 2000.0f}, "FV0Aamplitude"};
    ConfigurableAxis axisFDDAampl{"axisFDDAampl", {100, 0.0f, 2000.0f}, "FDDAamplitude"};
    ConfigurableAxis axisFDDCampl{"axisFDDCampl", {100, 0.0f, 2000.0f}, "FDDCamplitude"};
    ConfigurableAxis axisZNAampl{"axisZNAampl", {100, 0.0f, 250.0f}, "ZNAamplitude"};
    ConfigurableAxis axisZNCampl{"axisZNCampl", {100, 0.0f, 250.0f}, "ZNCamplitude"};
  } axisDetectors;

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", false, "if MC, do MC association"};
  Configurable<bool> doCollisionAssociationQA{"doCollisionAssociationQA", false, "check collision association"};

  // fast check on occupancy
  Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
  Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};

  // Kinematic axes
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for v0 analysis"};
  ConfigurableAxis axisPtXi{"axisPtCasc", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for cascade analysis"};
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

  ConfigurableAxis axisNTracksGlobal{"axisNTracksGlobal", {100, -0.5f, 99.5f}, "Number of global tracks"};
  ConfigurableAxis axisNTracksPVeta1{"axisNTracksPVeta1", {100, -0.5f, 99.5f}, "Number of PV contributors in |eta| < 1"};
  ConfigurableAxis axisNTracksTotalExceptITSonly{"axisNTracksTotalExceptITSonly", {100, -0.5f, 99.5f}, "Number of ITS-TPC and TPC only tracks"};
  ConfigurableAxis axisNchInvMass{"axisNchInvMass", {200, -0.5f, 199.5f}, "Number of charged particles for kTHnSparseF"};

  ConfigurableAxis axisFT0C_QA{"axisFT0C_QA",
                               {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 105.5},
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
  ConfigurableAxis axisITScluMap{"axisITSMap", {128, -0.5f, 127.5f}, "ITS Cluster map"};
  ConfigurableAxis axisDetMap{"axisDetMap", {16, -0.5f, 15.5f}, "Detector use map"};
  ConfigurableAxis axisITScluMapCoarse{"axisITScluMapCoarse", {16, -3.5f, 12.5f}, "ITS Coarse cluster map"};
  ConfigurableAxis axisDetMapCoarse{"axisDetMapCoarse", {5, -0.5f, 4.5f}, "Detector Coarse user map"};

  // Topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {80, -4.0f, 4.0f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {24, 0.0f, 1.2f}, "DCA (cm)"};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {100, 0.0f, 0.5f}, "pointing angle (rad)"};
  ConfigurableAxis axisV0Radius{"axisV0Radius", {60, 0.0f, 60.0f}, "V0 2D radius (cm)"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {200, -10.0f, 10.0f}, "N sigma TPC"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {200, 0.0f, 200.0f}, "TPC signal"};
  ConfigurableAxis axisTOFdeltaT{"axisTOFdeltaT", {200, -5000.0f, 5000.0f}, "TOF Delta T (ps)"};
  ConfigurableAxis axisNctau{"axisNctau", {100, 0.0f, 10.0f}, "n c x tau"};

  // PDG database
  Service<o2::framework::O2DatabasePDG> pdgDB;

  static constexpr std::string_view particlenames[] = {"K0Short", "Lambda", "AntiLambda", "Xi", "AntiXi", "Omega", "AntiOmega"};

  void setBits(std::bitset<selNum>& mask, std::initializer_list<int> selections)
  {
    for (int sel : selections) {
      mask.set(sel);
    }
  }

  template <int partID>
  void addTopoHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    if (isCascade) {
      histos.add(Form("%s/hCascCosPA", particlenames[partID].data()), "hCascCosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
      histos.add(Form("%s/hDCACascDaughters", particlenames[partID].data()), "hDCACascDaughters", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
      histos.add(Form("%s/hCascRadius", particlenames[partID].data()), "hCascRadius", kTH2D, {axisPtCoarse, {500, 0.0f, 50.0f}});
      histos.add(Form("%s/hMesonDCAToPV", particlenames[partID].data()), "hMesonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
      histos.add(Form("%s/hBaryonDCAToPV", particlenames[partID].data()), "hBaryonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
      histos.add(Form("%s/hBachDCAToPV", particlenames[partID].data()), "hBachDCAToPV", kTH2F, {axisPtCoarse, {200, -1.0f, 1.0f}});
      histos.add(Form("%s/hV0CosPA", particlenames[partID].data()), "hV0CosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
      histos.add(Form("%s/hV0Radius", particlenames[partID].data()), "hV0Radius", kTH2D, {axisPtCoarse, axisV0Radius});
      histos.add(Form("%s/hDCAV0Daughters", particlenames[partID].data()), "hDCAV0Daughters", kTH2F, {axisPtCoarse, axisDCAdau});
      histos.add(Form("%s/hDCAV0ToPV", particlenames[partID].data()), "hDCAV0ToPV", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
      histos.add(Form("%s/hMassLambdaDau", particlenames[partID].data()), "hMassLambdaDau", kTH2F, {axisPtCoarse, axisLambdaMass});
      histos.add(Form("%s/hNctau", particlenames[partID].data()), "hNctau", kTH2F, {axisPtCoarse, axisNctau});
      if (doBachelorBaryonCut) {
        histos.add(Form("%s/hBachBaryonCosPA", particlenames[partID].data()), "hBachBaryonCosPA", kTH2F, {axisPtCoarse, {100, 0.0f, 1.0f}});
        histos.add(Form("%s/hBachBaryonDCAxyToPV", particlenames[partID].data()), "hBachBaryonDCAxyToPV", kTH2F, {axisPtCoarse, {300, -3.0f, 3.0f}});
      }
    } else {
      histos.add(Form("%s/hPosDCAToPV", particlenames[partID].data()), "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add(Form("%s/hNegDCAToPV", particlenames[partID].data()), "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add(Form("%s/hDCADaughters", particlenames[partID].data()), "hDCADaughters", kTH1F, {axisDCAdau});
      histos.add(Form("%s/hPointingAngle", particlenames[partID].data()), "hPointingAngle", kTH1F, {axisPointingAngle});
      histos.add(Form("%s/hV0Radius", particlenames[partID].data()), "hV0Radius", kTH1F, {axisV0Radius});
      histos.add(Form("%s/h2dPositiveITSvsTPCpts", particlenames[partID].data()), "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add(Form("%s/h2dNegativeITSvsTPCpts", particlenames[partID].data()), "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    }
  }

  template <int partID>
  void addTPCQAHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    histos.add(Form("%s/h3dPosNsigmaTPC", particlenames[partID].data()), "h3dPosNsigmaTPC", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dNegNsigmaTPC", particlenames[partID].data()), "h3dNegNsigmaTPC", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dPosTPCsignal", particlenames[partID].data()), "h3dPosTPCsignal", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    histos.add(Form("%s/h3dNegTPCsignal", particlenames[partID].data()), "h3dNegTPCsignal", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});

    histos.add(Form("%s/h3dPosNsigmaTPCvsTrackPtot", particlenames[partID].data()), "h3dPosNsigmaTPCvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dNegNsigmaTPCvsTrackPtot", particlenames[partID].data()), "h3dNegNsigmaTPCvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});

    histos.add(Form("%s/h3dPosTPCsignalVsTrackPtot", particlenames[partID].data()), "h3dPosTPCsignalVsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    histos.add(Form("%s/h3dNegTPCsignalVsTrackPtot", particlenames[partID].data()), "h3dNegTPCsignalVsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});

    histos.add(Form("%s/h3dPosNsigmaTPCvsTrackPt", particlenames[partID].data()), "h3dPosNsigmaTPCvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(Form("%s/h3dNegNsigmaTPCvsTrackPt", particlenames[partID].data()), "h3dNegNsigmaTPCvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});

    histos.add(Form("%s/h3dPosTPCsignalVsTrackPt", particlenames[partID].data()), "h3dPosTPCsignalVsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    histos.add(Form("%s/h3dNegTPCsignalVsTrackPt", particlenames[partID].data()), "h3dNegTPCsignalVsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});

    if (isCascade) {
      histos.add(Form("%s/h3dBachTPCsignal", particlenames[partID].data()), "h3dBachTPCsignal", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
      histos.add(Form("%s/h3dBachNsigmaTPC", particlenames[partID].data()), "h3dBachNsigmaTPC", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
      histos.add(Form("%s/h3dBachNsigmaTPCvsTrackPtot", particlenames[partID].data()), "h3dBachNsigmaTPCvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
      histos.add(Form("%s/h3dBachTPCsignalVsTrackPtot", particlenames[partID].data()), "h3dBachTPCsignalVsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
      histos.add(Form("%s/h3dBachNsigmaTPCvsTrackPt", particlenames[partID].data()), "h3dBachNsigmaTPCvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
      histos.add(Form("%s/h3dBachTPCsignalVsTrackPt", particlenames[partID].data()), "h3dBachTPCsignalVsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    }
  }

  template <int partID>
  void addTOFQAHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    histos.add(Form("%s/h3dPosTOFdeltaT", particlenames[partID].data()), "h3dPosTOFdeltaT", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dNegTOFdeltaT", particlenames[partID].data()), "h3dNegTOFdeltaT", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dPosTOFdeltaTvsTrackPtot", particlenames[partID].data()), "h3dPosTOFdeltaTvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dNegTOFdeltaTvsTrackPtot", particlenames[partID].data()), "h3dNegTOFdeltaTvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dPosTOFdeltaTvsTrackPt", particlenames[partID].data()), "h3dPosTOFdeltaTvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(Form("%s/h3dNegTOFdeltaTvsTrackPt", particlenames[partID].data()), "h3dNegTOFdeltaTvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    if (isCascade) {
      histos.add(Form("%s/h3dBachTOFdeltaT", particlenames[partID].data()), "h3dBachTOFdeltaT", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
      histos.add(Form("%s/h3dBachTOFdeltaTvsTrackPtot", particlenames[partID].data()), "h3dBachTOFdeltaTvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
      histos.add(Form("%s/h3dBachTOFdeltaTvsTrackPt", particlenames[partID].data()), "h3dBachTOFdeltaTvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    }
  }

  template <int partID>
  void addKinematicQAHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    histos.add(Form("%s/h3dPosEtaPt", particlenames[partID].data()), "h3dPosEtaPt", kTH3F, {axisPtCoarse, axisEta, axisSelGap});
    histos.add(Form("%s/h3dNegEtaPt", particlenames[partID].data()), "h3dNegEtaPt", kTH3F, {axisPtCoarse, axisEta, axisSelGap});
    histos.add(Form("%s/h3dRapPt", particlenames[partID].data()), "h3dRapPt", kTH3F, {axisPtCoarse, axisRap, axisSelGap});
    if (isCascade) {
      histos.add(Form("%s/h3dBachEtaPt", particlenames[partID].data()), "h3dBachEtaPt", kTH3F, {axisPtCoarse, axisEta, axisSelGap});
    }
  }

  template <int partID>
  void addDetectorPropHistograms(HistogramRegistry& histos)
  {
    const bool isCascade = (partID > 2.5) ? true : false;
    if (doDetectPropQA == 1) {
      if (isCascade) {
        histos.add(Form("%s/h8dDetectPropVsCentrality", particlenames[partID].data()), "h8dDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
      } else {
        histos.add(Form("%s/h6dDetectPropVsCentrality", particlenames[partID].data()), "h6dDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
      }
      histos.add(Form("%s/h4dPosDetectPropVsCentrality", particlenames[partID].data()), "h4dPosDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse});
      histos.add(Form("%s/h4dNegDetectPropVsCentrality", particlenames[partID].data()), "h4dNegDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse});
      histos.add(Form("%s/h4dBachDetectPropVsCentrality", particlenames[partID].data()), "h4dBachDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse});
    }
    if (doDetectPropQA == 2) {
      if (isCascade) {
        histos.add(Form("%s/h9dDetectPropVsCentrality", particlenames[partID].data()), "h9dDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisInvMass.at(partID)});
      } else {
        histos.add(Form("%s/h7dDetectPropVsCentrality", particlenames[partID].data()), "h7dDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisInvMass.at(partID)});
      }
      histos.add(Form("%s/h5dPosDetectPropVsCentrality", particlenames[partID].data()), "h5dPosDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse, axisInvMass.at(partID)});
      histos.add(Form("%s/h5dNegDetectPropVsCentrality", particlenames[partID].data()), "h5dNegDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse, axisInvMass.at(partID)});
      histos.add(Form("%s/h5dBachDetectPropVsCentrality", particlenames[partID].data()), "h5dBachDetectPropVsCentrality", kTHnSparseF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse, axisInvMass.at(partID)});
    }
  }

  template <int partID>
  void addHistograms(HistogramRegistry& histos)
  {
    histos.add(Form("%s/h7dMass", particlenames[partID].data()), "h7dMass", kTHnSparseF, {axisFT0C, axisPt, axisInvMass.at(partID), axisSelGap, axisNchInvMass, axisRap, axisEta});
    histos.add(Form("%s/h2dMass", particlenames[partID].data()), "h2dMass", kTH2F, {axisInvMass.at(partID), axisSelGap});
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

  template <int partID>
  void addCollisionAssocHistograms(HistogramRegistry& histos)
  {
    histos.add(Form("%s/h2dPtVsNch", particlenames[partID].data()), "h2dPtVsNch", kTH2F, {axisMonteCarloNch, axisPt});
    histos.add(Form("%s/h2dPtVsNch_BadCollAssig", particlenames[partID].data()), "h2dPtVsNch_BadCollAssig", kTH2F, {axisMonteCarloNch, axisPt});
  }

  template <int partID, typename TCand, typename TCollision>
  void fillHistogramsV0(TCand cand, TCollision coll, int gap)
  {
    float invMass = 0;
    float centrality = coll.centFT0C();
    float pT = cand.pt();
    float rapidity = 1e6;

    float tpcNsigmaPos = 0;
    float tpcNsigmaNeg = 0;
    float tofDeltaTPos = 0;
    float tofDeltaTNeg = 0;

    auto posTrackExtra = cand.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = cand.template negTrackExtra_as<dauTracks>();

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

    histos.fill(HIST(particlenames[partID]) + HIST("/h2dMass"), invMass, gap);
    histos.fill(HIST(particlenames[partID]) + HIST("/h7dMass"), centrality, pT, invMass, gap, coll.multNTracksGlobal(), rapidity, cand.eta());
    if (doKienmaticQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosEtaPt"), pT, cand.positiveeta(), gap);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegEtaPt"), pT, cand.negativeeta(), gap);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dRapPt"), pT, rapidity, gap);
    }
    if (doPlainTopoQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/hPosDCAToPV"), cand.dcapostopv());
      histos.fill(HIST(particlenames[partID]) + HIST("/hNegDCAToPV"), cand.dcanegtopv());
      histos.fill(HIST(particlenames[partID]) + HIST("/hDCADaughters"), cand.dcaV0daughters());
      histos.fill(HIST(particlenames[partID]) + HIST("/hPointingAngle"), TMath::ACos(cand.v0cosPA()));
      histos.fill(HIST(particlenames[partID]) + HIST("/hV0Radius"), cand.v0radius());
      histos.fill(HIST(particlenames[partID]) + HIST("/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
      histos.fill(HIST(particlenames[partID]) + HIST("/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
    }
    if (doDetectPropQA == 1) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pT);
      histos.fill(HIST(particlenames[partID]) + HIST("/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT);
      histos.fill(HIST(particlenames[partID]) + HIST("/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT);
    }
    if (doDetectPropQA == 2) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h7dPosDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, pT, invMass);
      histos.fill(HIST(particlenames[partID]) + HIST("/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT, invMass);
      histos.fill(HIST(particlenames[partID]) + HIST("/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT, invMass);
    }
    if (PIDConfigurations.doTPCQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTPCsignal"), centrality, pT, posTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTPCsignal"), centrality, pT, negTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPtot"), centrality, cand.positivept() * TMath::CosH(cand.positiveeta()), posTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPtot"), centrality, cand.negativept() * TMath::CosH(cand.negativeeta()), negTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPt"), centrality, cand.positivept(), posTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPt"), centrality, cand.negativept(), negTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPt"), centrality, cand.positivept(), posTrackExtra.tpcNSigmaPi());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPt"), centrality, cand.negativept(), negTrackExtra.tpcNSigmaPi());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPtot"), centrality, cand.positivept() * TMath::CosH(cand.positiveeta()), tpcNsigmaPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPtot"), centrality, cand.negativept() * TMath::CosH(cand.negativeeta()), tpcNsigmaNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosNsigmaTPC"), centrality, pT, tpcNsigmaPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegNsigmaTPC"), centrality, pT, tpcNsigmaNeg);
    }
    if (PIDConfigurations.doTOFQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPt"), centrality, cand.positivept(), tofDeltaTPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPt"), centrality, cand.negativept(), tofDeltaTNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTOFdeltaT"), centrality, pT, tofDeltaTPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTOFdeltaT"), centrality, pT, tofDeltaTNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPtot"), centrality, cand.positivept() * TMath::CosH(cand.positiveeta()), tofDeltaTPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPtot"), centrality, cand.negativept() * TMath::CosH(cand.negativeeta()), tofDeltaTNeg);
    }
  }

  template <int partID, typename TCand, typename TCollision>
  void fillHistogramsCasc(TCand cand, TCollision coll, int gap)
  {
    float invMass = 0;
    float centrality = coll.centFT0C();
    float pT = cand.pt();
    float rapidity = 1e6;

    // Access daughter tracks
    auto posTrackExtra = cand.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = cand.template negTrackExtra_as<dauTracks>();
    auto bachTrackExtra = cand.template bachTrackExtra_as<dauTracks>();

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
      ctau = totalMom != 0 ? pdgDB->Mass(3312) * decayPos / (totalMom * ctauxiPDG) : 1e6;
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
      ctau = totalMom != 0 ? pdgDB->Mass(3312) * decayPos / (totalMom * ctauxiPDG) : 1e6;
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
      ctau = totalMom != 0 ? pdgDB->Mass(3334) * decayPos / (totalMom * ctauomegaPDG) : 1e6;
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
      ctau = totalMom != 0 ? pdgDB->Mass(3334) * decayPos / (totalMom * ctauomegaPDG) : 1e6;
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
    histos.fill(HIST(particlenames[partID]) + HIST("/h2dMass"), invMass, gap);
    histos.fill(HIST(particlenames[partID]) + HIST("/h7dMass"), centrality, pT, invMass, gap, coll.multNTracksGlobal(), rapidity, cand.eta());
    if (doKienmaticQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosEtaPt"), pT, cand.positiveeta(), gap);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegEtaPt"), pT, cand.negativeeta(), gap);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachEtaPt"), pT, cand.bacheloreta(), gap);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dRapPt"), pT, rapidity, gap);
    }
    if (doPlainTopoQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/hCascCosPA"), pT, cand.casccosPA(coll.posX(), coll.posY(), coll.posZ()));
      histos.fill(HIST(particlenames[partID]) + HIST("/hDCACascDaughters"), pT, cand.dcacascdaughters());
      histos.fill(HIST(particlenames[partID]) + HIST("/hCascRadius"), pT, cand.cascradius());
      histos.fill(HIST(particlenames[partID]) + HIST("/hMesonDCAToPV"), pT, cand.dcanegtopv());
      histos.fill(HIST(particlenames[partID]) + HIST("/hBaryonDCAToPV"), pT, cand.dcapostopv());
      histos.fill(HIST(particlenames[partID]) + HIST("/hBachDCAToPV"), pT, cand.dcabachtopv());
      histos.fill(HIST(particlenames[partID]) + HIST("/hV0CosPA"), pT, cand.v0cosPA(coll.posX(), coll.posY(), coll.posZ()));
      histos.fill(HIST(particlenames[partID]) + HIST("/hV0Radius"), pT, cand.v0radius());
      histos.fill(HIST(particlenames[partID]) + HIST("/hDCAV0Daughters"), pT, cand.dcaV0daughters());
      histos.fill(HIST(particlenames[partID]) + HIST("/hDCAV0ToPV"), pT, fabs(cand.dcav0topv(coll.posX(), coll.posY(), coll.posZ())));
      histos.fill(HIST(particlenames[partID]) + HIST("/hMassLambdaDau"), pT, cand.mLambda());
      histos.fill(HIST(particlenames[partID]) + HIST("/hNctau"), pT, ctau);
    }
    if (PIDConfigurations.doTPCQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosNsigmaTPC"), centrality, pT, tpcNsigmaPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegNsigmaTPC"), centrality, pT, tpcNsigmaNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachNsigmaTPC"), centrality, pT, tpcNsigmaBach);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTPCsignal"), centrality, pT, posTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTPCsignal"), centrality, pT, negTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachTPCsignal"), centrality, pT, bachTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPtot"), centrality, cand.positivept() * TMath::CosH(cand.positiveeta()), tpcNsigmaPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPtot"), centrality, cand.negativept() * TMath::CosH(cand.negativeeta()), tpcNsigmaNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachNsigmaTPCvsTrackPtot"), centrality, cand.bachelorpt() * TMath::CosH(cand.bacheloreta()), tpcNsigmaBach);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPtot"), centrality, cand.positivept() * TMath::CosH(cand.positiveeta()), posTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPtot"), centrality, cand.negativept() * TMath::CosH(cand.negativeeta()), negTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachTPCsignalVsTrackPtot"), centrality, cand.bachelorpt() * TMath::CosH(cand.bacheloreta()), bachTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosNsigmaTPCvsTrackPt"), centrality, cand.positivept(), tpcNsigmaPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegNsigmaTPCvsTrackPt"), centrality, cand.negativept(), tpcNsigmaNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachNsigmaTPCvsTrackPt"), centrality, cand.bachelorpt(), tpcNsigmaBach);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTPCsignalVsTrackPt"), centrality, cand.positivept(), posTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTPCsignalVsTrackPt"), centrality, cand.negativept(), negTrackExtra.tpcSignal());
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachTPCsignalVsTrackPt"), centrality, cand.bachelorpt(), bachTrackExtra.tpcSignal());
    }
    if (PIDConfigurations.doTOFQA) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTOFdeltaT"), centrality, pT, tofDeltaTPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTOFdeltaT"), centrality, pT, tofDeltaTNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachTOFdeltaT"), centrality, pT, tofDeltaTBach);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPtot"), centrality, cand.positivept() * TMath::CosH(cand.positiveeta()), tofDeltaTPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPtot"), centrality, cand.negativept() * TMath::CosH(cand.negativeeta()), tofDeltaTNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachTOFdeltaTvsTrackPtot"), centrality, cand.bachelorpt() * TMath::CosH(cand.bacheloreta()), tofDeltaTBach);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dPosTOFdeltaTvsTrackPt"), centrality, cand.positivept(), tofDeltaTPos);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dNegTOFdeltaTvsTrackPt"), centrality, cand.negativept(), tofDeltaTNeg);
      histos.fill(HIST(particlenames[partID]) + HIST("/h3dBachTOFdeltaTvsTrackPt"), centrality, cand.bachelorpt(), tofDeltaTBach);
    }
    if (doDetectPropQA == 1) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h8dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, bachDetMap, bachITSclusMap, pT);
      histos.fill(HIST(particlenames[partID]) + HIST("/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT);
      histos.fill(HIST(particlenames[partID]) + HIST("/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT);
      histos.fill(HIST(particlenames[partID]) + HIST("/h4dBachDetectPropVsCentrality"), centrality, bachTrackExtra.detectorMap(), bachTrackExtra.itsClusterMap(), pT);
    }
    if (doDetectPropQA == 2) {
      histos.fill(HIST(particlenames[partID]) + HIST("/h9dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, bachDetMap, bachITSclusMap, pT, invMass);
      histos.fill(HIST(particlenames[partID]) + HIST("/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), pT, invMass);
      histos.fill(HIST(particlenames[partID]) + HIST("/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), pT, invMass);
      histos.fill(HIST(particlenames[partID]) + HIST("/h5dBachDetectPropVsCentrality"), centrality, bachTrackExtra.detectorMap(), bachTrackExtra.itsClusterMap(), pT, invMass);
    }
  }

  void init(InitContext const&)
  {
    if ((doprocessV0s == true) && (doprocessCascades == true)) {
      LOG(fatal) << "Unable to analyze both v0s and cascades simultaneously. Please enable only one process at a time";
    }

    // initialise bit masks
    setBits(maskTopologicalV0, {selV0CosPA, selDCANegToPV, selDCAPosToPV, selDCAV0Dau, selV0Radius, selV0RadiusMax});
    setBits(maskTopologicalCasc, {selCascCosPA, selDCACascDau, selCascRadius, selCascRadiusMax, selBachToPV, selMesonToPV, selBaryonToPV,
                                  selDCAV0ToPV, selV0CosPA, selDCAV0Dau, selV0Radius, selV0RadiusMax, selLambdaMassWin});

    if (doBachelorBaryonCut)
      maskTopologicalCasc.set(selBachBaryon);

    setBits(maskKinematicV0, {selPosEta, selNegEta});
    setBits(maskKinematicCasc, {selPosEta, selNegEta, selBachEta});

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
      if (PIDConfigurations.TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific.set(selTPCPIDPositivePion);
        maskLambdaSpecific.set(selTPCPIDPositiveProton);
        maskAntiLambdaSpecific.set(selTPCPIDPositivePion);

        maskXiSpecific.set(selTPCPIDPositiveProton);
        maskAntiXiSpecific.set(selTPCPIDPositivePion);
        maskOmegaSpecific.set(selTPCPIDPositiveProton);
        maskAntiOmegaSpecific.set(selTPCPIDPositivePion);
      }
      // TOF PID
      if (PIDConfigurations.TofPidNsigmaCutK0Pi < 1e+5) { // safeguard for no cut
        setBits(maskK0ShortSpecific, {selTOFNSigmaPositivePionK0Short, selTOFDeltaTPositivePionK0Short});
      }
      if (PIDConfigurations.TofPidNsigmaCutLaPr < 1e+5) { // safeguard for no cut
        setBits(maskLambdaSpecific, {selTOFNSigmaPositiveProtonLambda, selTOFDeltaTPositiveProtonLambda});
        setBits(maskXiSpecific, {selTOFNSigmaPositiveProtonLambdaXi, selTOFDeltaTPositiveProtonLambdaXi});
        setBits(maskOmegaSpecific, {selTOFNSigmaPositiveProtonLambdaOmega, selTOFDeltaTPositiveProtonLambdaOmega});
      }
      if (PIDConfigurations.TofPidNsigmaCutLaPi < 1e+5) { // safeguard for no cut
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
      if (PIDConfigurations.TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskK0ShortSpecific.set(selTPCPIDNegativePion);
        maskLambdaSpecific.set(selTPCPIDNegativePion);
        maskAntiLambdaSpecific.set(selTPCPIDNegativeProton);

        maskXiSpecific.set(selTPCPIDNegativePion);
        maskAntiXiSpecific.set(selTPCPIDPositiveProton);
        maskOmegaSpecific.set(selTPCPIDNegativePion);
        maskAntiOmegaSpecific.set(selTPCPIDPositiveProton);
      }
      // TOF PID
      if (PIDConfigurations.TofPidNsigmaCutK0Pi < 1e+5) { // safeguard for no cut
        setBits(maskK0ShortSpecific, {selTOFNSigmaNegativePionK0Short, selTOFDeltaTNegativePionK0Short});
      }
      if (PIDConfigurations.TofPidNsigmaCutLaPr < 1e+5) { // safeguard for no cut
        setBits(maskAntiLambdaSpecific, {selTOFNSigmaNegativeProtonLambda, selTOFDeltaTNegativeProtonLambda});
        setBits(maskAntiXiSpecific, {selTOFNSigmaNegativeProtonLambdaXi, selTOFDeltaTNegativeProtonLambdaXi});
        setBits(maskAntiOmegaSpecific, {selTOFNSigmaNegativeProtonLambdaOmega, selTOFDeltaTNegativeProtonLambdaOmega});
      }
      if (PIDConfigurations.TofPidNsigmaCutLaPi < 1e+5) { // safeguard for no cut
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
      if (PIDConfigurations.TpcPidNsigmaCut < 1e+5) { // safeguard for no cut
        maskXiSpecific.set(selTPCPIDBachPion);
        maskAntiXiSpecific.set(selTPCPIDBachPion);
        maskOmegaSpecific.set(selTPCPIDBachKaon);
        maskAntiOmegaSpecific.set(selTPCPIDBachKaon);
      }
      // TOF PID
      if (PIDConfigurations.TofPidNsigmaCutXiPi < 1e+5) { // safeguard for no cut
        setBits(maskXiSpecific, {selTOFNSigmaBachPionXi, selTOFDeltaTBachPionXi});
        setBits(maskAntiXiSpecific, {selTOFNSigmaBachPionXi, selTOFDeltaTBachPionXi});
      }
      if (PIDConfigurations.TofPidNsigmaCutOmegaKaon < 1e+5) { // safeguard for no cut
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
    maskSelectionK0Short = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskK0ShortSpecific | (std::bitset<selNum>(1) << selPhysPrimK0Short);
    maskSelectionLambda = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskLambdaSpecific | (std::bitset<selNum>(1) << selPhysPrimLambda);
    maskSelectionAntiLambda = maskTopologicalV0 | maskKinematicV0 | maskTrackPropertiesV0 | maskAntiLambdaSpecific | (std::bitset<selNum>(1) << selPhysPrimAntiLambda);
    maskSelectionXi = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskXiSpecific | (std::bitset<selNum>(1) << selPhysPrimXi);
    maskSelectionAntiXi = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskAntiXiSpecific | (std::bitset<selNum>(1) << selPhysPrimAntiXi);
    maskSelectionOmega = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskOmegaSpecific | (std::bitset<selNum>(1) << selPhysPrimOmega);
    maskSelectionAntiOmega = maskTopologicalCasc | maskKinematicCasc | maskTrackPropertiesCasc | maskAntiOmegaSpecific | (std::bitset<selNum>(1) << selPhysPrimAntiOmega);

    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{16, -0.5f, +15.5f}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "kIsTriggerTVX");
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
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(15, "isUPC");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(16, "has UPC flag");

    // Event QA
    histos.add("eventQA/hCentrality", "hCentrality", kTH2F, {axisFT0C_QA, axisSelGap});
    histos.add("eventQA/hCentralityVsTracksPVeta1", "hCentralityVsTracksPVeta1", kTH3F, {axisFT0C_QA, axisNTracksPVeta1, axisSelGap});
    histos.add("eventQA/hCentralityVsTracksTotalExceptITSonly", "hCentralityVsTracksTotalExceptITSonly", kTH3F, {axisFT0C_QA, axisNTracksTotalExceptITSonly, axisSelGap});
    histos.add("eventQA/hOccupancy", "hOccupancy", kTH2F, {axisOccupancy, axisSelGap});
    histos.add("eventQA/hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH3F, {axisFT0C_QA, axisOccupancy, axisSelGap});
    histos.add("eventQA/hTracksPVeta1VsTracksGlobal", "hTracksPVeta1VsTracksGlobal", kTH3F, {axisNTracksPVeta1, axisNTracksGlobal, axisSelGap});
    histos.add("eventQA/hCentralityVsTracksGlobal", "hCentralityVsTracksGlobal", kTH3F, {axisFT0C_QA, axisNTracksGlobal, axisSelGap});
    histos.add("eventQA/hGapSide", "Gap side; Entries", kTH1F, {{5, -0.5, 4.5}});
    histos.add("eventQA/hSelGapSide", "Selected gap side; Entries", kTH1F, {axisSelGap});
    histos.add("eventQA/hPosX", "Vertex position in x", kTH2F, {{100, -0.1, 0.1}, axisSelGap});
    histos.add("eventQA/hPosY", "Vertex position in y", kTH2F, {{100, -0.1, 0.1}, axisSelGap});
    histos.add("eventQA/hPosZ", "Vertex position in z", kTH2F, {{100, -20., 20.}, axisSelGap});
    histos.add("eventQA/hFT0", "hFT0", kTH3F, {axisDetectors.axisFT0Aampl, axisDetectors.axisFT0Campl, axisSelGap});
    histos.add("eventQA/hFDD", "hFDD", kTH3F, {axisDetectors.axisFDDAampl, axisDetectors.axisFDDCampl, axisSelGap});
    histos.add("eventQA/hZN", "hZN", kTH3F, {axisDetectors.axisZNAampl, axisDetectors.axisZNCampl, axisSelGap});

    if (doprocessV0s) {
      // For all candidates
      if (doPlainTopoQA) {
        histos.add("generalQA/hPt", "hPt", kTH1F, {axisPtCoarse});
        histos.add("generalQA/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("generalQA/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
        histos.add("generalQA/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
        histos.add("generalQA/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
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
        if (doCollisionAssociationQA) {
          addCollisionAssocHistograms<1>(histos);
        }
      }

      // Anti-Lambda
      if (analyseAntiLambda) {
        addHistograms<2>(histos);
        if (doCollisionAssociationQA) {
          addCollisionAssocHistograms<2>(histos);
        }
      }
    }

    if (doprocessCascades) {
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
    int selGapSide = sgSelector.trueGap(collision, upcCuts.FV0cut, upcCuts.FT0Acut, upcCuts.FT0Ccut, upcCuts.ZDCcut);
    histos.fill(HIST("eventQA/hGapSide"), collision.gapSide());
    histos.fill(HIST("eventQA/hSelGapSide"), selGapSide);
    histos.fill(HIST("eventQA/hFT0"), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), selGapSide);
    histos.fill(HIST("eventQA/hFDD"), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC(), selGapSide);
    histos.fill(HIST("eventQA/hZN"), collision.energyCommonZNA(), collision.energyCommonZNC(), selGapSide);
    return selGapSide;
  }

  template <typename TCollision>
  void fillHistogramsQA(TCollision const& collision, int const& gap)
  {
    // QA histograms
    float centrality = collision.centFT0C();
    histos.fill(HIST("eventQA/hCentrality"), centrality, gap);
    histos.fill(HIST("eventQA/hCentralityVsTracksTotalExceptITSonly"), centrality, collision.multAllTracksTPCOnly() + collision.multAllTracksITSTPC(), gap);
    histos.fill(HIST("eventQA/hCentralityVsTracksPVeta1"), centrality, collision.multNTracksPVeta1(), gap);
    histos.fill(HIST("eventQA/hOccupancy"), collision.trackOccupancyInTimeRange(), gap);
    histos.fill(HIST("eventQA/hCentralityVsOccupancy"), centrality, collision.trackOccupancyInTimeRange(), gap);
    histos.fill(HIST("eventQA/hTracksPVeta1VsTracksGlobal"), collision.multNTracksPVeta1(), collision.multNTracksGlobal(), gap);
    histos.fill(HIST("eventQA/hCentralityVsTracksGlobal"), centrality, collision.multNTracksGlobal(), gap);
    histos.fill(HIST("eventQA/hPosX"), collision.posX(), gap);
    histos.fill(HIST("eventQA/hPosY"), collision.posY(), gap);
    histos.fill(HIST("eventQA/hPosZ"), collision.posZ(), gap);
  }

  template <typename TCollision>
  bool acceptEvent(TCollision const& collision)
  {
    histos.fill(HIST("hEventSelection"), 0. /* all collisions */);

    if (requireIsTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 1 /* triggered by FT0M */);

    if (fabs(collision.posZ()) > 10.f) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 2 /* vertex-Z selected */);

    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);

    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);

    if (requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 5 /* Contains at least one ITS-TPC track */);

    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 6 /* PV position consistency check */);

    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 7 /* PV with at least one contributor matched with TOF */);

    if (requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TRD */);

    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 9 /* Not at same bunch pile-up */);

    if (requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 10 /* No other collision within +/- 10 microseconds */);

    if (requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 4 microseconds */);

    if (minOccupancy > 0 && collision.trackOccupancyInTimeRange() < minOccupancy) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 12 /* Above min occupancy */);

    if (maxOccupancy > 0 && collision.trackOccupancyInTimeRange() > maxOccupancy) {
      return false;
    }
    histos.fill(HIST("hEventSelection"), 13 /* Below max occupancy */);

    if (studyUPConly && !collision.isUPC()) {
      return false;
    } else if (collision.isUPC()) {
      histos.fill(HIST("hEventSelection"), 14 /* is UPC compatible */);
    }

    if (useUPCflag && (collision.flags() < 1)) {
      return false;
    } else if (collision.flags() >= 1) {
      histos.fill(HIST("hEventSelection"), 15 /* UPC event */);
    }

    return true;
  }

  bool verifyMask(std::bitset<selNum> bitmap, std::bitset<selNum> mask)
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

    constexpr MaskBitmapPair configs[] = {
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

    for (const auto& config : configs) {
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

    constexpr MaskBitmapPair configs[] = {
      {o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD | o2::aod::track::TOF, 4}, // ITS-TPC-TRD-TOF
      {o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TOF, 3},                       // ITS-TPC-TOF
      {o2::aod::track::ITS | o2::aod::track::TPC | o2::aod::track::TRD, 2},                       // ITS-TPC-TRD
      {o2::aod::track::ITS | o2::aod::track::TPC, 1}                                              // ITS-TPC
    };

    for (const auto& config : configs) {
      if (verifyMask(detMap, config.mask)) {
        bitMap = config.bitmap;
        break;
      }
    }

    return bitMap;
  }

  template <typename TCasc, typename TCollision>
  std::bitset<selNum> computeBitmapCascade(TCasc const& casc, TCollision const& coll)
  {
    float rapidityXi = casc.yXi();
    float rapidityOmega = casc.yOmega();

    // Access daughter tracks
    auto posTrackExtra = casc.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = casc.template negTrackExtra_as<dauTracks>();
    auto bachTrackExtra = casc.template bachTrackExtra_as<dauTracks>();

    // c x tau
    float decayPos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
    float totalMom = std::hypot(casc.px(), casc.py(), casc.pz());
    float ctauXi = totalMom != 0 ? pdgDB->Mass(3312) * decayPos / totalMom : 1e6;
    float ctauOmega = totalMom != 0 ? pdgDB->Mass(3334) * decayPos / totalMom : 1e6;

    std::bitset<selNum> bitMap = 0;

    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) > casccuts.casccospa)
      bitMap.set(selCascCosPA);
    if (casc.dcacascdaughters() < casccuts.dcacascdau)
      bitMap.set(selDCACascDau);
    if (casc.cascradius() > casccuts.cascradius)
      bitMap.set(selCascRadius);
    if (casc.cascradius() < casccuts.cascradiusMax)
      bitMap.set(selCascRadiusMax);
    if (doBachelorBaryonCut && (casc.bachBaryonCosPA() < casccuts.bachbaryoncospa) && (fabs(casc.bachBaryonDCAxyToPV()) > casccuts.bachbaryondcaxytopv))
      bitMap.set(selBachBaryon);
    if (fabs(casc.dcabachtopv()) > casccuts.dcabachtopv)
      bitMap.set(selBachToPV);

    if (casc.sign() > 0) {
      if (fabs(casc.dcanegtopv()) > casccuts.dcabaryontopv)
        bitMap.set(selBaryonToPV);
      if (fabs(casc.dcapostopv()) > casccuts.dcamesontopv)
        bitMap.set(selMesonToPV);
    } else { // no sign == 0, in principle
      if (fabs(casc.dcapostopv()) > casccuts.dcabaryontopv)
        bitMap.set(selBaryonToPV);
      if (fabs(casc.dcanegtopv()) > casccuts.dcamesontopv)
        bitMap.set(selMesonToPV);
    }

    if (fabs(casc.mXi() - pdgDB->Mass(3312)) < casccuts.masswin)
      bitMap.set(selMassWinXi);
    if (fabs(casc.mOmega() - pdgDB->Mass(3334)) < casccuts.masswin)
      bitMap.set(selMassWinOmega);
    if (fabs(casc.mLambda() - pdgDB->Mass(3122)) < casccuts.lambdamasswin)
      bitMap.set(selLambdaMassWin);

    if (fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) > casccuts.dcav0topv)
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
    if (fabs(rapidityXi) < rapidityCut)
      bitMap.set(selXiRapidity);
    if (fabs(rapidityOmega) < rapidityCut)
      bitMap.set(selOmegaRapidity);
    if (fabs(poseta) < daughterEtaCut)
      bitMap.set(selNegEta);
    if (fabs(negeta) < daughterEtaCut)
      bitMap.set(selPosEta);
    if (fabs(bacheta) < daughterEtaCut)
      bitMap.set(selBachEta);

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
    if (fabs(posTrackExtra.tpcNSigmaPi()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositivePion);
    if (fabs(posTrackExtra.tpcNSigmaPr()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositiveProton);
    // negative track
    if (fabs(negTrackExtra.tpcNSigmaPi()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativePion);
    if (fabs(negTrackExtra.tpcNSigmaPr()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativeProton);
    // bachelor track
    if (fabs(bachTrackExtra.tpcNSigmaPi()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDBachPion);
    if (fabs(bachTrackExtra.tpcNSigmaKa()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDBachKaon);

    // TOF PID in DeltaT
    // positive track
    if (fabs(casc.posTOFDeltaTXiPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTPositiveProtonLambdaXi);
    if (fabs(casc.posTOFDeltaTXiPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionLambdaXi);
    if (fabs(casc.posTOFDeltaTOmPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTPositiveProtonLambdaOmega);
    if (fabs(casc.posTOFDeltaTOmPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionLambdaOmega);
    // negative track
    if (fabs(casc.negTOFDeltaTXiPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTNegativeProtonLambdaXi);
    if (fabs(casc.negTOFDeltaTXiPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionLambdaXi);
    if (fabs(casc.negTOFDeltaTOmPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTNegativeProtonLambdaOmega);
    if (fabs(casc.negTOFDeltaTOmPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionLambdaOmega);
    // bachelor track
    if (fabs(casc.bachTOFDeltaTOmKa()) < PIDConfigurations.maxDeltaTimeKaon)
      bitMap.set(selTOFDeltaTBachKaonOmega);
    if (fabs(casc.bachTOFDeltaTXiPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTBachPionXi);

    // TOF PID in NSigma
    // meson track
    if (fabs(casc.tofNSigmaXiLaPi()) < PIDConfigurations.TofPidNsigmaCutLaPi) {
      bitMap.set(selTOFNSigmaPositivePionLambdaXi);
      bitMap.set(selTOFNSigmaNegativePionLambdaXi);
    }
    if (fabs(casc.tofNSigmaOmLaPi()) < PIDConfigurations.TofPidNsigmaCutLaPi) {
      bitMap.set(selTOFNSigmaPositivePionLambdaOmega);
      bitMap.set(selTOFNSigmaNegativePionLambdaOmega);
    }
    // baryon track
    if (fabs(casc.tofNSigmaXiLaPr()) < PIDConfigurations.TofPidNsigmaCutLaPr) {
      bitMap.set(selTOFNSigmaNegativeProtonLambdaXi);
      bitMap.set(selTOFNSigmaPositiveProtonLambdaXi);
    }
    if (fabs(casc.tofNSigmaOmLaPr()) < PIDConfigurations.TofPidNsigmaCutLaPr) {
      bitMap.set(selTOFNSigmaNegativePionLambdaOmega);
      bitMap.set(selTOFNSigmaPositivePionLambdaOmega);
    }
    // bachelor track
    if (fabs(casc.tofNSigmaXiPi()) < PIDConfigurations.TofPidNsigmaCutXiPi) {
      bitMap.set(selTOFNSigmaBachPionXi);
    }
    if (fabs(casc.tofNSigmaOmKa()) < PIDConfigurations.TofPidNsigmaCutOmegaKaon) {
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
    if (fabs(casc.mOmega() - pdgDB->Mass(3334)) > casccuts.rejcomp)
      bitMap.set(selRejCompXi);
    if (fabs(casc.mXi() - pdgDB->Mass(3312)) > casccuts.rejcomp)
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
  std::bitset<selNum> computeBitmapV0(TV0 const& v0, TCollision const& collision)
  {
    float rapidityLambda = v0.yLambda();
    float rapidityK0Short = v0.yK0Short();

    std::bitset<selNum> bitMap = 0;

    // base topological variables
    if (v0.v0radius() > v0cuts.v0radius)
      bitMap.set(selV0Radius);
    if (v0.v0radius() < v0cuts.v0radiusMax)
      bitMap.set(selV0RadiusMax);
    if (fabs(v0.dcapostopv()) > v0cuts.dcapostopv)
      bitMap.set(selDCAPosToPV);
    if (fabs(v0.dcanegtopv()) > v0cuts.dcanegtopv)
      bitMap.set(selDCANegToPV);
    if (v0.v0cosPA() > v0cuts.v0cospa)
      bitMap.set(selV0CosPA);
    if (v0.dcaV0daughters() < v0cuts.dcav0dau)
      bitMap.set(selDCAV0Dau);

    // kinematic
    if (fabs(rapidityLambda) < rapidityCut)
      bitMap.set(selLambdaRapidity);
    if (fabs(rapidityK0Short) < rapidityCut)
      bitMap.set(selK0ShortRapidity);
    if (fabs(v0.negativeeta()) < daughterEtaCut)
      bitMap.set(selNegEta);
    if (fabs(v0.positiveeta()) < daughterEtaCut)
      bitMap.set(selPosEta);

    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

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
    if (fabs(posTrackExtra.tpcNSigmaPi()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositivePion);
    if (fabs(posTrackExtra.tpcNSigmaPr()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDPositiveProton);
    if (fabs(negTrackExtra.tpcNSigmaPi()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativePion);
    if (fabs(negTrackExtra.tpcNSigmaPr()) < PIDConfigurations.TpcPidNsigmaCut)
      bitMap.set(selTPCPIDNegativeProton);

    // TOF PID in DeltaT
    // positive track
    if (fabs(v0.posTOFDeltaTLaPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTPositiveProtonLambda);
    if (fabs(v0.posTOFDeltaTLaPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionLambda);
    if (fabs(v0.posTOFDeltaTK0Pi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTPositivePionK0Short);
    // negative track
    if (fabs(v0.negTOFDeltaTLaPr()) < PIDConfigurations.maxDeltaTimeProton)
      bitMap.set(selTOFDeltaTNegativeProtonLambda);
    if (fabs(v0.negTOFDeltaTLaPi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionLambda);
    if (fabs(v0.negTOFDeltaTK0Pi()) < PIDConfigurations.maxDeltaTimePion)
      bitMap.set(selTOFDeltaTNegativePionK0Short);

    // TOF PID in NSigma
    // positive track
    if (fabs(v0.tofNSigmaLaPr()) < PIDConfigurations.TofPidNsigmaCutLaPr)
      bitMap.set(selTOFNSigmaPositiveProtonLambda);
    if (fabs(v0.tofNSigmaALaPi()) < PIDConfigurations.TofPidNsigmaCutLaPi)
      bitMap.set(selTOFNSigmaPositivePionLambda);
    if (fabs(v0.tofNSigmaK0PiPlus()) < PIDConfigurations.TofPidNsigmaCutK0Pi)
      bitMap.set(selTOFNSigmaPositivePionK0Short);
    // negative track
    if (fabs(v0.tofNSigmaALaPr()) < PIDConfigurations.TofPidNsigmaCutLaPr)
      bitMap.set(selTOFNSigmaNegativeProtonLambda);
    if (fabs(v0.tofNSigmaLaPi()) < PIDConfigurations.TofPidNsigmaCutLaPi)
      bitMap.set(selTOFNSigmaNegativePionLambda);
    if (fabs(v0.tofNSigmaK0PiMinus()) < PIDConfigurations.TofPidNsigmaCutK0Pi)
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
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecutV0->get("lifetimecutLambda"))
      bitMap.set(selLambdaCTau);
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecutV0->get("lifetimecutK0S"))
      bitMap.set(selK0ShortCTau);

    // armenteros
    if (v0.qtarm() * v0cuts.armPodCut > fabs(v0.alpha()) || v0cuts.armPodCut < 1e-4)
      bitMap.set(selK0ShortArmenteros);

    return bitMap;
  }

  template <typename TCasc, typename TCollision>
  void analyseCascCandidate(TCasc const& casc, TCollision const& coll, int const& gap, std::bitset<selNum> const& selMap)
  {
    // Access daughter tracks
    auto posTrackExtra = casc.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = casc.template negTrackExtra_as<dauTracks>();
    auto bachTrackExtra = casc.template bachTrackExtra_as<dauTracks>();

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
      histos.fill(HIST("generalQA/hDCAV0ToPV"), casc.pt(), fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())));
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

  template <typename TV0, typename TCollision>
  void analyseV0Candidate(TV0 const& v0, TCollision const& coll, int const& gap, std::bitset<selNum> const& selMap)
  {
    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

    // QA plots
    if (doPlainTopoQA) {
      histos.fill(HIST("generalQA/hPt"), v0.pt());
      histos.fill(HIST("generalQA/hPosDCAToPV"), v0.dcapostopv());
      histos.fill(HIST("generalQA/hNegDCAToPV"), v0.dcanegtopv());
      histos.fill(HIST("generalQA/hDCADaughters"), v0.dcaV0daughters());
      histos.fill(HIST("generalQA/hPointingAngle"), TMath::ACos(v0.v0cosPA()));
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

  void processV0s(straCollisonFull const& collision, v0Candidates const& fullV0s, dauTracks const&)
  {
    if (!acceptEvent(collision)) {
      return;
    } // event is accepted

    int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
    if (studyUPConly && (selGapSide < -0.5))
      return;

    fillHistogramsQA(collision, selGapSide);

    for (auto& v0 : fullV0s) {
      if ((v0.v0Type() != v0cuts.v0TypeSelection) && (v0cuts.v0TypeSelection > 0))
        continue; // skip V0s that are not standard

      std::bitset<selNum> selMap = computeBitmapV0(v0, collision);

      // consider for histograms for all species
      setBits(selMap, {selConsiderK0Short, selConsiderLambda, selConsiderAntiLambda,
                       selPhysPrimK0Short, selPhysPrimLambda, selPhysPrimAntiLambda});

      analyseV0Candidate(v0, collision, selGapSide, selMap);
    } // end v0 loop
  }

  void processCascades(straCollisonFull const& collision, cascadeCandidates const& fullCascades, dauTracks const&)
  {
    if (!acceptEvent(collision)) {
      return;
    } // event is accepted

    int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
    if (studyUPConly && (selGapSide < -0.5))
      return;

    fillHistogramsQA(collision, selGapSide);

    for (auto& casc : fullCascades) {
      std::bitset<selNum> selMap = computeBitmapCascade(casc, collision);
      // consider for histograms for all species
      setBits(selMap, {selConsiderXi, selConsiderAntiXi, selConsiderOmega, selConsiderAntiOmega,
                       selPhysPrimXi, selPhysPrimAntiXi, selPhysPrimOmega, selPhysPrimAntiOmega});

      analyseCascCandidate(casc, collision, selGapSide, selMap);
    } // end casc loop
  }

  PROCESS_SWITCH(strangeYieldPbPb, processV0s, "Process V0s", true);
  PROCESS_SWITCH(strangeYieldPbPb, processCascades, "Process Cascades", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeYieldPbPb>(cfgc)};
}
