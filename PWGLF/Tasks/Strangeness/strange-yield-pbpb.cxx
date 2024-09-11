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
  Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", true, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
  Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
  Configurable<bool> studyUPConly{"studyUPConly", false, "is UPC-only analysis"};

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
    Configurable<int> doDetectPropQA{"doDetectPropQA", 0, "do Detector/ITS map QA: 0: no, 1: 4D, 2: 5D with mass"};

    Configurable<bool> doPlainTopoQA{"doPlainTopoQA", true, "do simple 1D QA of candidates"};
    Configurable<float> qaMinPt{"qaMinPt", 0.0f, "minimum pT for QA plots"};
    Configurable<float> qaMaxPt{"qaMaxPt", 1000.0f, "maximum pT for QA plots"};

    // PID (TOF)
    Configurable<float> maxDeltaTimeProton{"maxDeltaTimeProton", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimePion{"maxDeltaTimePion", 1e+9, "check maximum allowed time"};
    Configurable<float> maxDeltaTimeKaon{"maxDeltaTimeKaon", 1e+9, "check maximum allowed time"};
  } PIDConfigurations;

  struct : ConfigurableGroup {
    ConfigurableAxis axisFT0Aampl{"FT0Aamplitude", {100, 0.0f, 2000.0f}, "FT0Aamplitude"};
    ConfigurableAxis axisFT0Campl{"FT0Camplitude", {100, 0.0f, 2000.0f}, "FT0Camplitude"};
    ConfigurableAxis axisFV0Aampl{"FV0Aamplitude", {100, 0.0f, 2000.0f}, "FV0Aamplitude"};
    ConfigurableAxis axisFDDAampl{"FDDAamplitude", {100, 0.0f, 2000.0f}, "FDDAamplitude"};
    ConfigurableAxis axisFDDCampl{"FDDCamplitude", {100, 0.0f, 2000.0f}, "FDDCamplitude"};
    ConfigurableAxis axisZNAampl{"ZNAamplitude", {100, 0.0f, 250.0f}, "ZNAamplitude"};
    ConfigurableAxis axisZNCampl{"ZNCamplitude", {100, 0.0f, 250.0f}, "ZNCamplitude"};
  } axisDetectors;

  // for MC
  Configurable<bool> doMCAssociation{"doMCAssociation", true, "if MC, do MC association"};
  Configurable<bool> doCollisionAssociationQA{"doCollisionAssociationQA", true, "check collision association"};

  // fast check on occupancy
  Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
  Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for v0 analysis"};
  ConfigurableAxis axisPtXi{"axisPtCasc", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for cascade analysis"};
  ConfigurableAxis axisPtCoarse{"axisPtCoarse", {VARIABLE_WIDTH, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 7.0f, 10.0f, 15.0f}, "pt axis for QA"};

  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.28f, 1.36f}, ""};
  ConfigurableAxis axisOmegaMass{"axisOmegaMass", {200, 1.59f, 1.75f}, ""};

  ConfigurableAxis axisNch{"axisNch", {1000, -0.5f, 999.5f}, "Number of charged particles"};
  ConfigurableAxis axisFT0C{"FT0C",
                            {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 105.5},
                            "FT0C (%)"};

  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};

  // UPC axes
  ConfigurableAxis axisSelGap{"axisSelGap", {4, -1.5, 2.5}, "Gap side"};

  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // Track quality axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, -0.5f, 159.5f}, "N TPC rows"};
  ConfigurableAxis axisITSclus{"axisITSclus", {7, -0.5f, 6.5f}, "N ITS Clusters"};
  ConfigurableAxis axisITScluMap{"axisITSMap", {128, -0.5f, 127.5f}, "ITS Cluster map"};
  ConfigurableAxis axisDetMap{"axisDetMap", {16, -0.5f, 15.5f}, "Detector use map"};
  ConfigurableAxis axisITScluMapCoarse{"axisITScluMapCoarse", {16, -3.5f, 12.5f}, "ITS Coarse cluster map"};
  ConfigurableAxis axisDetMapCoarse{"axisDetMapCoarse", {5, -0.5f, 4.5f}, "Detector Coarse user map"};

  // MC coll assoc QA axis
  ConfigurableAxis axisMonteCarloNch{"axisMonteCarloNch", {300, 0.0f, 3000.0f}, "N_{ch} MC"};
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

  void setBits(std::bitset<selNum>& mask, std::initializer_list<int> selections)
  {
    for (int sel : selections) {
      mask.set(sel);
    }
  }

  void addTopoHistograms(HistogramRegistry& histos, const std::string& particleName, const bool isCascade)
  {
    if (isCascade) {
      histos.add(particleName + "/hCascCosPA", "hCascCosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
      histos.add(particleName + "/hDCACascDaughters", "hDCACascDaughters", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
      histos.add(particleName + "/hCascRadius", "hCascRadius", kTH2D, {axisPtCoarse, {500, 0.0f, 50.0f}});
      histos.add(particleName + "/hMesonDCAToPV", "hMesonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
      histos.add(particleName + "/hBaryonDCAToPV", "hBaryonDCAToPV", kTH2F, {axisPtCoarse, axisDCAtoPV});
      histos.add(particleName + "/hBachDCAToPV", "hBachDCAToPV", kTH2F, {axisPtCoarse, {200, -1.0f, 1.0f}});
      histos.add(particleName + "/hV0CosPA", "hV0CosPA", kTH2F, {axisPtCoarse, {100, 0.9f, 1.0f}});
      histos.add(particleName + "/hV0Radius", "hV0Radius", kTH2D, {axisPtCoarse, axisV0Radius});
      histos.add(particleName + "/hDCAV0Daughters", "hDCAV0Daughters", kTH2F, {axisPtCoarse, axisDCAdau});
      histos.add(particleName + "/hDCAV0ToPV", "hDCAV0ToPV", kTH2F, {axisPtCoarse, {44, 0.0f, 2.2f}});
      histos.add(particleName + "/hMassLambdaDau", "hMassLambdaDau", kTH2F, {axisPtCoarse, axisLambdaMass});
      histos.add(particleName + "/hNctau", "hNctau", kTH2F, {axisPtCoarse, axisNctau});
      if (doBachelorBaryonCut) {
        histos.add(particleName + "/hBachBaryonCosPA", "hBachBaryonCosPA", kTH2F, {axisPtCoarse, {100, 0.0f, 1.0f}});
        histos.add(particleName + "/hBachBaryonDCAxyToPV", "hBachBaryonDCAxyToPV", kTH2F, {axisPtCoarse, {300, -3.0f, 3.0f}});
      }
    } else {
      histos.add(particleName + "/hPosDCAToPV", "hPosDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add(particleName + "/hNegDCAToPV", "hNegDCAToPV", kTH1F, {axisDCAtoPV});
      histos.add(particleName + "/hDCADaughters", "hDCADaughters", kTH1F, {axisDCAdau});
      histos.add(particleName + "/hPointingAngle", "hPointingAngle", kTH1F, {axisPointingAngle});
      histos.add(particleName + "/hV0Radius", "hV0Radius", kTH1F, {axisV0Radius});
      histos.add(particleName + "/h2dPositiveITSvsTPCpts", "h2dPositiveITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
      histos.add(particleName + "/h2dNegativeITSvsTPCpts", "h2dNegativeITSvsTPCpts", kTH2F, {axisTPCrows, axisITSclus});
    }
  }

  void addTPCQAHistograms(HistogramRegistry& histos, const std::string& particleName, const bool isCascade)
  {
    histos.add(particleName + "/h3dPosNsigmaTPC", "h3dPosNsigmaTPC", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(particleName + "/h3dNegNsigmaTPC", "h3dNegNsigmaTPC", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(particleName + "/h3dPosTPCsignal", "h3dPosTPCsignal", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    histos.add(particleName + "/h3dNegTPCsignal", "h3dNegTPCsignal", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});

    histos.add(particleName + "/h3dPosNsigmaTPCvsTrackPtot", "h3dPosNsigmaTPCvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(particleName + "/h3dNegNsigmaTPCvsTrackPtot", "h3dNegNsigmaTPCvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});

    histos.add(particleName + "/h3dPosTPCsignalVsTrackPtot", "h3dPosTPCsignalVsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    histos.add(particleName + "/h3dNegTPCsignalVsTrackPtot", "h3dNegTPCsignalVsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});

    histos.add(particleName + "/h3dPosNsigmaTPCvsTrackPt", "h3dPosNsigmaTPCvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
    histos.add(particleName + "/h3dNegNsigmaTPCvsTrackPt", "h3dNegNsigmaTPCvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});

    histos.add(particleName + "/h3dPosTPCsignalVsTrackPt", "h3dPosTPCsignalVsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    histos.add(particleName + "/h3dNegTPCsignalVsTrackPt", "h3dNegTPCsignalVsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});

    if (isCascade) {
      histos.add(particleName + "/h3dBachTPCsignal", "h3dBachTPCsignal", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
      histos.add(particleName + "/h3dBachNsigmaTPC", "h3dBachNsigmaTPC", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
      histos.add(particleName + "/h3dBachNsigmaTPCvsTrackPtot", "h3dBachNsigmaTPCvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
      histos.add(particleName + "/h3dBachTPCsignalVsTrackPtot", "h3dBachTPCsignalVsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
      histos.add(particleName + "/h3dBachNsigmaTPCvsTrackPt", "h3dBachNsigmaTPCvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisNsigmaTPC});
      histos.add(particleName + "/h3dBachTPCsignalVsTrackPt", "h3dBachTPCsignalVsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTPCsignal});
    }
  }

  void addTOFQAHistograms(HistogramRegistry& histos, const std::string& particleName, const bool isCascade)
  {
    histos.add(particleName + "/h3dPosTOFdeltaT", "h3dPosTOFdeltaT", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(particleName + "/h3dNegTOFdeltaT", "h3dNegTOFdeltaT", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(particleName + "/h3dPosTOFdeltaTvsTrackPtot", "h3dPosTOFdeltaTvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(particleName + "/h3dNegTOFdeltaTvsTrackPtot", "h3dNegTOFdeltaTvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(particleName + "/h3dPosTOFdeltaTvsTrackPt", "h3dPosTOFdeltaTvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    histos.add(particleName + "/h3dNegTOFdeltaTvsTrackPt", "h3dNegTOFdeltaTvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    if (isCascade) {
      histos.add(particleName + "/h3dBachTOFdeltaT", "h3dBachTOFdeltaT", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
      histos.add(particleName + "/h3dBachTOFdeltaTvsTrackPtot", "h3dBachTOFdeltaTvsTrackPtot", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
      histos.add(particleName + "/h3dBachTOFdeltaTvsTrackPt", "h3dBachTOFdeltaTvsTrackPt", kTH3F, {axisFT0C, axisPtCoarse, axisTOFdeltaT});
    }
  }

  void addDetectorPropHistograms(HistogramRegistry& histos, const std::string& particleName, const bool isCascade)
  {
    if (PIDConfigurations.doDetectPropQA == 1) {
      if (isCascade) {
        histos.add(particleName + "/h8dDetectPropVsCentrality", "h8dDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
      } else {
        histos.add(particleName + "/h6dDetectPropVsCentrality", "h6dDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse});
      }
      histos.add(particleName + "/h4dPosDetectPropVsCentrality", "h4dPosDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse});
      histos.add(particleName + "/h4dNegDetectPropVsCentrality", "h4dNegDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse});
      histos.add(particleName + "/h4dBachDetectPropVsCentrality", "h4dBachDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse});
    }
    if (PIDConfigurations.doDetectPropQA == 2) {
      if (isCascade) {
        histos.add(particleName + "/h9dDetectPropVsCentrality", "h9dDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisXiMass});
      } else {
        histos.add(particleName + "/h7dDetectPropVsCentrality", "h7dDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMapCoarse, axisITScluMapCoarse, axisDetMapCoarse, axisITScluMapCoarse, axisPtCoarse, axisLambdaMass});
      }
      histos.add(particleName + "/h5dPosDetectPropVsCentrality", "h5dPosDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse, axisXiMass});
      histos.add(particleName + "/h5dNegDetectPropVsCentrality", "h5dNegDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse, axisXiMass});
      histos.add(particleName + "/h5dBachDetectPropVsCentrality", "h5dBachDetectPropVsCentrality", kTHnF, {axisFT0C, axisDetMap, axisITScluMap, axisPtCoarse, axisXiMass});
    }
  }

  void addHistograms(HistogramRegistry& histos, const std::string& particleName, const bool isCascade)
  {
    if (PIDConfigurations.doPlainTopoQA) {
      addTopoHistograms(histos, particleName, isCascade);
    }
    if (PIDConfigurations.doTPCQA) {
      addTPCQAHistograms(histos, particleName, isCascade);
    }
    if (PIDConfigurations.doTOFQA) {
      addTOFQAHistograms(histos, particleName, isCascade);
    }
    addDetectorPropHistograms(histos, particleName, isCascade);
  }

  void addCollisionAssocHistograms(HistogramRegistry& histos, const std::string& particleName)
  {
    histos.add(particleName + "/h2dPtVsNch", "h2dPtVsNch", kTH2F, {axisMonteCarloNch, axisPt});
    histos.add(particleName + "/h2dPtVsNch_BadCollAssig", "h2dPtVsNch_BadCollAssig", kTH2F, {axisMonteCarloNch, axisPt});
  }

  // template <typename TCasc, typename TCollision>
  // void fillTopoHistograms(TCasc casc, TCollision coll, int gap, std::bitset<selNum> selMap) {
  // }

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
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{15, -0.5f, +14.5f}});
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

    // Event QA
    histos.add("eventQA/hEventCentrality", "hEventCentrality", kTH1F, {axisFT0C});
    histos.add("eventQA/hCentralityVsNch", "hCentralityVsNch", kTH2F, {axisFT0C, axisNch});
    histos.add("eventQA/hEventOccupancy", "hEventOccupancy", kTH1F, {axisOccupancy});
    histos.add("eventQA/hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH2F, {axisFT0C, axisOccupancy});
    histos.add("eventQA/hGapSide", "Gap side; Entries", kTH1F, {{5, -0.5, 4.5}});
    histos.add("eventQA/hSelGapSide", "Selected gap side; Entries", kTH1F, {axisSelGap});
    histos.add("eventQA/hPosX", "Vertex position in x", kTH1F, {{100, -0.1, 0.1}});
    histos.add("eventQA/hPosY", "Vertex position in y", kTH1F, {{100, -0.1, 0.1}});
    histos.add("eventQA/hPosZ", "Vertex position in z", kTH1F, {{100, -20., 20.}});
    histos.add("eventQA/hFT0", "hFT0", kTH3F, {axisDetectors.axisFT0Aampl, axisDetectors.axisFT0Campl, axisSelGap});
    histos.add("eventQA/hFDD", "hFDD", kTH3F, {axisDetectors.axisFDDAampl, axisDetectors.axisFDDCampl, axisSelGap});
    histos.add("eventQA/hZN", "hZN", kTH3F, {axisDetectors.axisZNAampl, axisDetectors.axisZNCampl, axisSelGap});

    if (doprocessV0s) {
      // For all candidates
      if (PIDConfigurations.doPlainTopoQA) {
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
        histos.add("K0Short/h4dMass", "h4dMass", kTHnF, {axisFT0C, axisPt, axisK0Mass, axisSelGap});
        histos.add("K0Short/h2dMass", "h2dMass", kTH2F, {axisK0Mass, axisSelGap});
        addHistograms(histos, "K0Short", false);
      }

      // Lambda
      if (analyseLambda) {
        histos.add("Lambda/h4dMass", "h4dMass", kTHnF, {axisFT0C, axisPt, axisLambdaMass, axisSelGap});
        histos.add("Lambda/h2dMass", "h2dMass", kTH2F, {axisLambdaMass, axisSelGap});
        addHistograms(histos, "Lambda", false);
        if (doCollisionAssociationQA) {
          addCollisionAssocHistograms(histos, "Lambda");
        }
      }

      // Anti-Lambda
      if (analyseAntiLambda) {
        histos.add("AntiLambda/h4dMass", "h4dMass", kTHnF, {axisFT0C, axisPt, axisLambdaMass, axisSelGap});
        histos.add("AntiLambda/h2dMass", "h2dMass", kTH2F, {axisLambdaMass, axisSelGap});
        addHistograms(histos, "AntiLambda", false);
        if (doCollisionAssociationQA) {
          addCollisionAssocHistograms(histos, "AntiLambda");
        }
      }
    }

    if (doprocessCascades) {
      // For all candidates
      if (PIDConfigurations.doPlainTopoQA) {
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
        histos.add("Xi/h4dMass", "h4dMass", kTHnF, {axisFT0C, axisPt, axisXiMass, axisSelGap});
        histos.add("Xi/h2dMass", "h2dMass", kTH2F, {axisXiMass, axisSelGap});
        addHistograms(histos, "Xi", true);
      }

      // Anti-Xi
      if (analyseAntiXi) {
        histos.add("AntiXi/h2dMass", "h2dMass", kTH2F, {axisXiMass, axisSelGap});
        addHistograms(histos, "AntiXi", true);
      }

      // Omega
      if (analyseOmega) {
        histos.add("Omega/h2dMass", "h2dMass", kTH2F, {axisOmegaMass, axisSelGap});
        addHistograms(histos, "Omega", true);
      }

      // Anti-Omega
      if (analyseAntiOmega) {
        histos.add("AntiOmega/h2dMass", "h2dMass", kTH2F, {axisOmegaMass, axisSelGap});
        addHistograms(histos, "AntiOmega", true);
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

    // QA histograms
    float centrality = collision.centFT0C();
    histos.fill(HIST("eventQA/hEventCentrality"), centrality);
    histos.fill(HIST("eventQA/hCentralityVsNch"), centrality, collision.multNTracksPVeta1());
    histos.fill(HIST("eventQA/hEventOccupancy"), collision.trackOccupancyInTimeRange());
    histos.fill(HIST("eventQA/hCentralityVsOccupancy"), centrality, collision.trackOccupancyInTimeRange());
    histos.fill(HIST("eventQA/hPosX"), collision.posX());
    histos.fill(HIST("eventQA/hPosY"), collision.posY());
    histos.fill(HIST("eventQA/hPosZ"), collision.posZ());

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
    if (doBachelorBaryonCut && casc.bachBaryonCosPA() < casccuts.bachbaryoncospa && fabs(casc.bachBaryonDCAxyToPV()) > casccuts.bachbaryondcaxytopv)
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
    float decayPos = std::hypot(casc.x() - coll.posX(), casc.y() - coll.posY(), casc.z() - coll.posZ());
    float totalMom = std::hypot(casc.px(), casc.py(), casc.pz());
    float ctauXi = totalMom != 0 ? pdgDB->Mass(3312) * decayPos / totalMom : 1e6;
    float ctauOmega = totalMom != 0 ? pdgDB->Mass(3334) * decayPos / totalMom : 1e6;

    float centrality = coll.centFT0C();

    if (PIDConfigurations.doPlainTopoQA) {
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
      histos.fill(HIST("Xi/h2dMass"), casc.mXi(), gap);
      histos.fill(HIST("Xi/h4dMass"), centrality, casc.pt(), casc.mXi(), gap);
      if (PIDConfigurations.doPlainTopoQA) {
        histos.fill(HIST("Xi/hCascCosPA"), casc.pt(), casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()));
        histos.fill(HIST("Xi/hDCACascDaughters"), casc.pt(), casc.dcacascdaughters());
        histos.fill(HIST("Xi/hCascRadius"), casc.pt(), casc.cascradius());
        histos.fill(HIST("Xi/hMesonDCAToPV"), casc.pt(), casc.dcanegtopv());
        histos.fill(HIST("Xi/hBaryonDCAToPV"), casc.pt(), casc.dcapostopv());
        histos.fill(HIST("Xi/hBachDCAToPV"), casc.pt(), casc.dcabachtopv());
        histos.fill(HIST("Xi/hV0CosPA"), casc.pt(), casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()));
        histos.fill(HIST("Xi/hV0Radius"), casc.pt(), casc.v0radius());
        histos.fill(HIST("Xi/hDCAV0Daughters"), casc.pt(), casc.dcaV0daughters());
        histos.fill(HIST("Xi/hDCAV0ToPV"), casc.pt(), fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())));
        histos.fill(HIST("Xi/hMassLambdaDau"), casc.pt(), casc.mLambda());
        histos.fill(HIST("Xi/hNctau"), casc.pt(), ctauXi / ctauxiPDG);
      }
      if (PIDConfigurations.doTPCQA) {
        histos.fill(HIST("Xi/h3dPosNsigmaTPC"), centrality, casc.pt(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dNegNsigmaTPC"), centrality, casc.pt(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dBachNsigmaTPC"), centrality, casc.pt(), bachTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dPosTPCsignal"), centrality, casc.pt(), posTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dNegTPCsignal"), centrality, casc.pt(), negTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dBachTPCsignal"), centrality, casc.pt(), bachTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dPosNsigmaTPCvsTrackPtot"), centrality, casc.positivept() * TMath::CosH(casc.positiveeta()), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dNegNsigmaTPCvsTrackPtot"), centrality, casc.negativept() * TMath::CosH(casc.negativeeta()), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dBachNsigmaTPCvsTrackPtot"), centrality, casc.bachelorpt() * TMath::CosH(casc.bacheloreta()), bachTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dPosTPCsignalVsTrackPtot"), centrality, casc.positivept() * TMath::CosH(casc.positiveeta()), posTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dNegTPCsignalVsTrackPtot"), centrality, casc.negativept() * TMath::CosH(casc.negativeeta()), negTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dBachTPCsignalVsTrackPtot"), centrality, casc.bachelorpt() * TMath::CosH(casc.bacheloreta()), bachTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dPosNsigmaTPCvsTrackPt"), centrality, casc.positivept(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dNegNsigmaTPCvsTrackPt"), centrality, casc.negativept(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dBachNsigmaTPCvsTrackPt"), centrality, casc.bachelorpt(), bachTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Xi/h3dPosTPCsignalVsTrackPt"), centrality, casc.positivept(), posTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dNegTPCsignalVsTrackPt"), centrality, casc.negativept(), negTrackExtra.tpcSignal());
        histos.fill(HIST("Xi/h3dBachTPCsignalVsTrackPt"), centrality, casc.bachelorpt(), bachTrackExtra.tpcSignal());
      }
      if (PIDConfigurations.doTOFQA) {
        histos.fill(HIST("Xi/h3dPosTOFdeltaT"), centrality, casc.pt(), casc.posTOFDeltaTXiPr());
        histos.fill(HIST("Xi/h3dNegTOFdeltaT"), centrality, casc.pt(), casc.negTOFDeltaTXiPi());
        histos.fill(HIST("Xi/h3dBachTOFdeltaT"), centrality, casc.pt(), casc.bachTOFDeltaTXiPi());
        histos.fill(HIST("Xi/h3dPosTOFdeltaTvsTrackPtot"), centrality, casc.positivept() * TMath::CosH(casc.positiveeta()), casc.posTOFDeltaTXiPr());
        histos.fill(HIST("Xi/h3dNegTOFdeltaTvsTrackPtot"), centrality, casc.negativept() * TMath::CosH(casc.negativeeta()), casc.negTOFDeltaTXiPi());
        histos.fill(HIST("Xi/h3dBachTOFdeltaTvsTrackPtot"), centrality, casc.bachelorpt() * TMath::CosH(casc.bacheloreta()), casc.bachTOFDeltaTXiPi());
        histos.fill(HIST("Xi/h3dPosTOFdeltaTvsTrackPt"), centrality, casc.positivept(), casc.posTOFDeltaTXiPr());
        histos.fill(HIST("Xi/h3dNegTOFdeltaTvsTrackPt"), centrality, casc.negativept(), casc.negTOFDeltaTXiPi());
        histos.fill(HIST("Xi/h3dBachTOFdeltaTvsTrackPt"), centrality, casc.bachelorpt(), casc.bachTOFDeltaTXiPi());
      }
      if (PIDConfigurations.doDetectPropQA == 1) {
        histos.fill(HIST("Xi/h8dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, bachDetMap, bachITSclusMap, casc.pt());
        histos.fill(HIST("Xi/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), casc.pt());
        histos.fill(HIST("Xi/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), casc.pt());
        histos.fill(HIST("Xi/h4dBachDetectPropVsCentrality"), centrality, bachTrackExtra.detectorMap(), bachTrackExtra.itsClusterMap(), casc.pt());
      }
      if (PIDConfigurations.doDetectPropQA == 2) {
        histos.fill(HIST("Xi/h9dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, bachDetMap, bachITSclusMap, casc.pt(), casc.mXi());
        histos.fill(HIST("Xi/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), casc.pt(), casc.mXi());
        histos.fill(HIST("Xi/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), casc.pt(), casc.mXi());
        histos.fill(HIST("Xi/h5dBachDetectPropVsCentrality"), centrality, bachTrackExtra.detectorMap(), bachTrackExtra.itsClusterMap(), casc.pt(), casc.mXi());
      }
    }

    // Anti-Xi
    if (verifyMask(selMap, maskSelectionAntiXi) && analyseAntiXi) {
      histos.fill(HIST("AntiXi/h2dMass"), casc.mXi(), gap);
    }

    // Omega
    if (verifyMask(selMap, maskSelectionOmega) && analyseOmega) {
      histos.fill(HIST("Omega/h2dMass"), casc.mOmega(), gap);
    }

    // Anti-Omega
    if (verifyMask(selMap, maskSelectionAntiOmega) && analyseAntiOmega) {
      histos.fill(HIST("AntiOmega/h2dMass"), casc.mOmega(), gap);
    }
  }

  template <typename TV0>
  void analyseV0Candidate(TV0 v0, float centrality, int gap, std::bitset<selNum> selMap)
  {
    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

    bool posIsFromAfterburner = posTrackExtra.itsChi2PerNcl() < 0;
    bool negIsFromAfterburner = negTrackExtra.itsChi2PerNcl() < 0;

    uint posDetMap = computeDetBitmap(posTrackExtra.detectorMap());
    int posITSclusMap = computeITSclusBitmap(posTrackExtra.itsClusterMap(), posIsFromAfterburner);
    uint negDetMap = computeDetBitmap(negTrackExtra.detectorMap());
    int negITSclusMap = computeITSclusBitmap(negTrackExtra.itsClusterMap(), negIsFromAfterburner);

    // QA plots
    if (PIDConfigurations.doPlainTopoQA) {
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
      histos.fill(HIST("K0Short/h2dMass"), v0.mK0Short(), gap);
      histos.fill(HIST("K0Short/h4dMass"), centrality, v0.pt(), v0.mK0Short(), gap);
      histos.fill(HIST("generalQA/h2dArmenterosSelected"), v0.alpha(), v0.qtarm());
      if (PIDConfigurations.doPlainTopoQA) {
        histos.fill(HIST("K0Short/hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("K0Short/hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("K0Short/hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("K0Short/hPointingAngle"), TMath::ACos(v0.v0cosPA()));
        histos.fill(HIST("K0Short/hV0Radius"), v0.v0radius());
        histos.fill(HIST("K0Short/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("K0Short/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      }
      if (PIDConfigurations.doDetectPropQA == 1) {
        histos.fill(HIST("K0Short/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, v0.pt());
        histos.fill(HIST("K0Short/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), v0.pt());
        histos.fill(HIST("K0Short/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), v0.pt());
      }
      if (PIDConfigurations.doDetectPropQA == 2) {
        histos.fill(HIST("K0Short/h7dPosDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, v0.pt(), v0.mK0Short());
        histos.fill(HIST("K0Short/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), v0.pt(), v0.mK0Short());
        histos.fill(HIST("K0Short/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), v0.pt(), v0.mK0Short());
      }
      if (PIDConfigurations.doTPCQA) {
        histos.fill(HIST("K0Short/h3dPosNsigmaTPC"), centrality, v0.pt(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dNegNsigmaTPC"), centrality, v0.pt(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dPosTPCsignal"), centrality, v0.pt(), posTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dNegTPCsignal"), centrality, v0.pt(), negTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dPosTPCsignalVsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dNegTPCsignalVsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("K0Short/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
        histos.fill(HIST("K0Short/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
      }
      if (PIDConfigurations.doTOFQA) {
        histos.fill(HIST("K0Short/h3dPosTOFdeltaT"), centrality, v0.pt(), v0.posTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dNegTOFdeltaT"), centrality, v0.pt(), v0.negTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), v0.posTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), v0.negTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTK0Pi());
        histos.fill(HIST("K0Short/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTK0Pi());
      }
    }

    // Lambda
    if (verifyMask(selMap, maskSelectionLambda) && analyseLambda) {
      histos.fill(HIST("Lambda/h2dMass"), v0.mLambda(), gap);
      histos.fill(HIST("Lambda/h4dMass"), centrality, v0.pt(), v0.mLambda(), gap);
      if (PIDConfigurations.doPlainTopoQA) {
        histos.fill(HIST("Lambda/hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("Lambda/hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("Lambda/hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("Lambda/hPointingAngle"), TMath::ACos(v0.v0cosPA()));
        histos.fill(HIST("Lambda/hV0Radius"), v0.v0radius());
        histos.fill(HIST("Lambda/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("Lambda/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      }
      if (PIDConfigurations.doDetectPropQA == 1) {
        histos.fill(HIST("Lambda/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, v0.pt());
        histos.fill(HIST("Lambda/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), v0.pt());
        histos.fill(HIST("Lambda/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), v0.pt());
      }
      if (PIDConfigurations.doDetectPropQA == 2) {
        histos.fill(HIST("Lambda/h7dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, v0.pt(), v0.mLambda());
        histos.fill(HIST("Lambda/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), v0.pt(), v0.mLambda());
        histos.fill(HIST("Lambda/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), v0.pt(), v0.mLambda());
      }
      if (PIDConfigurations.doTPCQA) {
        histos.fill(HIST("Lambda/h3dPosNsigmaTPC"), centrality, v0.pt(), posTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("Lambda/h3dNegNsigmaTPC"), centrality, v0.pt(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Lambda/h3dPosTPCsignal"), centrality, v0.pt(), posTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dNegTPCsignal"), centrality, v0.pt(), negTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("Lambda/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Lambda/h3dPosTPCsignalVsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dNegTPCsignalVsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("Lambda/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("Lambda/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
        histos.fill(HIST("Lambda/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
      }
      if (PIDConfigurations.doTOFQA) {
        histos.fill(HIST("Lambda/h3dPosTOFdeltaT"), centrality, v0.pt(), v0.posTOFDeltaTLaPr());
        histos.fill(HIST("Lambda/h3dNegTOFdeltaT"), centrality, v0.pt(), v0.negTOFDeltaTLaPi());
        histos.fill(HIST("Lambda/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), v0.posTOFDeltaTLaPr());
        histos.fill(HIST("Lambda/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), v0.negTOFDeltaTLaPi());
        histos.fill(HIST("Lambda/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTLaPr());
        histos.fill(HIST("Lambda/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTLaPi());
      }
    }

    // Anti-Lambda
    if (verifyMask(selMap, maskSelectionAntiLambda) && analyseAntiLambda) {
      histos.fill(HIST("AntiLambda/h2dMass"), v0.mAntiLambda(), gap);
      histos.fill(HIST("AntiLambda/h4dMass"), centrality, v0.pt(), v0.mAntiLambda(), gap);
      if (PIDConfigurations.doPlainTopoQA) {
        histos.fill(HIST("AntiLambda/hPosDCAToPV"), v0.dcapostopv());
        histos.fill(HIST("AntiLambda/hNegDCAToPV"), v0.dcanegtopv());
        histos.fill(HIST("AntiLambda/hDCADaughters"), v0.dcaV0daughters());
        histos.fill(HIST("AntiLambda/hPointingAngle"), TMath::ACos(v0.v0cosPA()));
        histos.fill(HIST("AntiLambda/hV0Radius"), v0.v0radius());
        histos.fill(HIST("AntiLambda/h2dPositiveITSvsTPCpts"), posTrackExtra.tpcCrossedRows(), posTrackExtra.itsNCls());
        histos.fill(HIST("AntiLambda/h2dNegativeITSvsTPCpts"), negTrackExtra.tpcCrossedRows(), negTrackExtra.itsNCls());
      }
      if (PIDConfigurations.doDetectPropQA == 1) {
        histos.fill(HIST("AntiLambda/h6dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, v0.pt());
        histos.fill(HIST("AntiLambda/h4dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), v0.pt());
        histos.fill(HIST("AntiLambda/h4dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), v0.pt());
      }
      if (PIDConfigurations.doDetectPropQA == 2) {
        histos.fill(HIST("AntiLambda/h7dDetectPropVsCentrality"), centrality, posDetMap, posITSclusMap, negDetMap, negITSclusMap, v0.pt(), v0.mAntiLambda());
        histos.fill(HIST("AntiLambda/h5dPosDetectPropVsCentrality"), centrality, posTrackExtra.detectorMap(), posTrackExtra.itsClusterMap(), v0.pt(), v0.mAntiLambda());
        histos.fill(HIST("AntiLambda/h5dNegDetectPropVsCentrality"), centrality, negTrackExtra.detectorMap(), negTrackExtra.itsClusterMap(), v0.pt(), v0.mAntiLambda());
      }
      if (PIDConfigurations.doTPCQA) {
        histos.fill(HIST("AntiLambda/h3dPosNsigmaTPC"), centrality, v0.pt(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("AntiLambda/h3dNegNsigmaTPC"), centrality, v0.pt(), negTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("AntiLambda/h3dPosTPCsignal"), centrality, v0.pt(), posTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dNegTPCsignal"), centrality, v0.pt(), negTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dPosNsigmaTPCvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("AntiLambda/h3dNegNsigmaTPCvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("AntiLambda/h3dPosTPCsignalVsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), posTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dNegTPCsignalVsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), negTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dPosNsigmaTPCvsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcNSigmaPi());
        histos.fill(HIST("AntiLambda/h3dNegNsigmaTPCvsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcNSigmaPr());
        histos.fill(HIST("AntiLambda/h3dPosTPCsignalVsTrackPt"), centrality, v0.positivept(), posTrackExtra.tpcSignal());
        histos.fill(HIST("AntiLambda/h3dNegTPCsignalVsTrackPt"), centrality, v0.negativept(), negTrackExtra.tpcSignal());
      }
      if (PIDConfigurations.doTOFQA) {
        histos.fill(HIST("AntiLambda/h3dPosTOFdeltaT"), centrality, v0.pt(), v0.posTOFDeltaTLaPi());
        histos.fill(HIST("AntiLambda/h3dNegTOFdeltaT"), centrality, v0.pt(), v0.negTOFDeltaTLaPr());
        histos.fill(HIST("AntiLambda/h3dPosTOFdeltaTvsTrackPtot"), centrality, v0.positivept() * TMath::CosH(v0.positiveeta()), v0.posTOFDeltaTLaPi());
        histos.fill(HIST("AntiLambda/h3dNegTOFdeltaTvsTrackPtot"), centrality, v0.negativept() * TMath::CosH(v0.negativeeta()), v0.negTOFDeltaTLaPr());
        histos.fill(HIST("AntiLambda/h3dPosTOFdeltaTvsTrackPt"), centrality, v0.positivept(), v0.posTOFDeltaTLaPi());
        histos.fill(HIST("AntiLambda/h3dNegTOFdeltaTvsTrackPt"), centrality, v0.negativept(), v0.negTOFDeltaTLaPr());
      }
    }
  }

  void processV0s(straCollisonFull const& collision, v0Candidates const& fullV0s, dauTracks const&)
  {
    if (!acceptEvent(collision)) {
      return;
    } // event is accepted

    int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
    if (studyUPConly && selGapSide < -0.5)
      return;

    for (auto& v0 : fullV0s) {
      if (v0.v0Type() != v0cuts.v0TypeSelection && v0cuts.v0TypeSelection > 0)
        continue; // skip V0s that are not standard

      std::bitset<selNum> selMap = computeBitmapV0(v0, collision);

      // consider for histograms for all species
      setBits(selMap, {selConsiderK0Short, selConsiderLambda, selConsiderAntiLambda,
                       selPhysPrimK0Short, selPhysPrimLambda, selPhysPrimAntiLambda});

      analyseV0Candidate(v0, collision.centFT0C(), selGapSide, selMap);
    } // end v0 loop
  }

  void processCascades(straCollisonFull const& collision, cascadeCandidates const& fullCascades, dauTracks const&)
  {
    if (!acceptEvent(collision)) {
      return;
    } // event is accepted

    int selGapSide = collision.isUPC() ? getGapSide(collision) : -1;
    if (studyUPConly && selGapSide < -0.5)
      return;

    for (auto& casc : fullCascades) {
      std::bitset<selNum> selMap = computeBitmapCascade(casc, collision);
      // consider for histograms for all species
      setBits(selMap, {selConsiderXi, selConsiderAntiXi, selConsiderOmega, selConsiderAntiOmega,
                       selPhysPrimXi, selPhysPrimAntiXi, selPhysPrimOmega, selPhysPrimAntiOmega});

      analyseCascCandidate(casc, collision, selGapSide, selMap);
    } // end casc loop
  }

  PROCESS_SWITCH(strangeYieldPbPb, processV0s, "Process V0s in UPC", true);
  PROCESS_SWITCH(strangeYieldPbPb, processCascades, "Process Cascades in UPC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeYieldPbPb>(cfgc)};
}
