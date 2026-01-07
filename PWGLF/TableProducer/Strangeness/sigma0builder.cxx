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
// This is a task that employs the standard V0 tables and attempts to combine
// two V0s into a Sigma0 -> Lambda + gamma candidate.
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  Sigma0 builder task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//

#include "PWGLF/DataModel/LFSigmaTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector3D.h"
#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0StandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0DerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0TOFStandardDerivedDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;
using V0TOFDerivedMCDatas = soa::Join<aod::V0Cores, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0MCMothers, aod::V0CoreMCLabels, aod::V0LambdaMLScores, aod::V0AntiLambdaMLScores, aod::V0GammaMLScores>;

static const std::vector<std::string> DirList = {"V0BeforeSel", "PhotonSel", "LambdaSel"};
;

struct sigma0builder {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  //___________________________________________________
  // KStar Specific

  Produces<aod::KStarCores> kstarcores;                // kstar candidates info for analysis
  Produces<aod::KShortExtras> kshortExtras; // lambdas from sigma0 candidates info
  Produces<aod::KStarPhotonExtras> kstarPhotonExtras; // photons from kstar candidates info
  Produces<aod::KStarCollRef> kstarCollRefs;      // references collisions from kstarcores
  Produces<aod::KStarMCCores> kstarmccores;          // Reco sigma0 MC properties
  Produces<aod::KStarGens> kstarGens;                // Generated sigma0s
  Produces<aod::KStarGenCollRef> kstarGenCollRefs;    // references collisions from sigma0Gens


  //__________________________________________________
  // Sigma0 specific
  Produces<aod::SigmaIndices> sigmaIndices;            // references V0Cores from sigma0Gens
  Produces<aod::Sigma0Cores> sigma0cores;              // sigma0 candidates info for analysis
  Produces<aod::Sigma0PhotonExtras> sigmaPhotonExtras; // photons from sigma0 candidates info
  Produces<aod::Sigma0LambdaExtras> sigmaLambdaExtras; // lambdas from sigma0 candidates info
  Produces<aod::SigmaCollRef> sigma0CollRefs;          // references collisions from Sigma0Cores
  Produces<aod::Sigma0MCCores> sigma0mccores;          // Reco sigma0 MC properties
  Produces<aod::SigmaMCLabels> sigma0mclabel;          // Link of reco sigma0 to mcparticles
  Produces<aod::Sigma0Gens> sigma0Gens;                // Generated sigma0s
  Produces<aod::SigmaGenCollRef> sigma0GenCollRefs;    // references collisions from sigma0Gens[

  //__________________________________________________
  // Pi0 specific
  Produces<aod::Pi0Cores> pi0cores;            // pi0 candidates info for analysis
  Produces<aod::Pi0CollRef> pi0coresRefs;      // references collisions from photonpair
  Produces<aod::Pi0CoresMC> pi0coresmc;        // Reco pi0 MC properties
  Produces<aod::Pi0Gens> pi0Gens;              // Generated pi0s
  Produces<aod::Pi0GenCollRef> pi0GenCollRefs; // references collisions from pi0Gens

  //__________________________________________________
  // pack track quality but separte also afterburner
  // dynamic range: 0-31
  enum selection : int { hasTPC = 0,
                         hasITSTracker,
                         hasITSAfterburner,
                         hasTRD,
                         hasTOF };

  // Histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> fFillV03DPositionHistos{"fFillV03DPositionHistos", false, "Fill XYZ histo for Photons and Lambdas."};
  Configurable<bool> fFillNoSelV0Histos{"fFillNoSelV0Histos", false, "Fill QA histos for input V0s."};
  Configurable<bool> fFillSelPhotonHistos{"fFillSelPhotonHistos", true, "Fill QA histos for sel photons."};
  Configurable<bool> fFillSelLambdaHistos{"fFillSelLambdaHistos", true, "Fill QA histos for sel lambdas."};

  Configurable<bool> doAssocStudy{"doAssocStudy", false, "Do v0 to collision association study."};
  Configurable<bool> doPPAnalysis{"doPPAnalysis", true, "if in pp, set to true"};

  Configurable<bool> fGetIR{"fGetIR", false, "Flag to retrieve the IR info."};
  Configurable<bool> fIRCrashOnNull{"fIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash."};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> fUseEventSelection{"fUseEventSelection", false, "Apply event selection cuts"};
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", true, "require events with at least one ITS-TPC track"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", false, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeVzDep{"requireNoCollInTimeRangeVzDep", false, "reject collisions corrupted by the cannibalism, with other collisions with pvZ of drifting TPC tracks from past/future collisions within 2.5 cm the current pvZ"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF"};
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

  // Tables to fill
  Configurable<bool> fillPi0Tables{"fillPi0Tables", false, "fill pi0 tables for QA"};
  Configurable<bool> fillSigma0Tables{"fillSigma0Tables", true, "fill sigma0 tables for analysis"};
  Configurable<bool> fillKStarTables{"fillKStarTables", true, "fill kstar tables for analysis"};

  // For ML Selection
  Configurable<bool> useMLScores{"useMLScores", false, "use ML scores to select candidates"};

  // For standard approach:

  // Lambda criteria:
  struct : ConfigurableGroup {
    std::string prefix = "lambdaSelections"; // JSON group name
    Configurable<float> Lambda_MLThreshold{"Lambda_MLThreshold", 0.1, "Decision Threshold value to select lambdas"};
    Configurable<float> AntiLambda_MLThreshold{"AntiLambda_MLThreshold", 0.1, "Decision Threshold value to select antilambdas"};
    Configurable<float> doMCAssociation{"doMCAssociation", false, "if MC, select true lambda/alambdas only"};
    Configurable<float> LambdaMinDCANegToPv{"LambdaMinDCANegToPv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> LambdaMinDCAPosToPv{"LambdaMinDCAPosToPv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> LambdaMaxDCAV0Dau{"LambdaMaxDCAV0Dau", 2.5, "Max DCA V0 Daughters (cm)"};
    Configurable<float> LambdaMinv0radius{"LambdaMinv0radius", 0.0, "Min V0 radius (cm)"};
    Configurable<float> LambdaMaxv0radius{"LambdaMaxv0radius", 40, "Max V0 radius (cm)"};
    Configurable<float> LambdaMinQt{"LambdaMinQt", 0.01, "Min lambda qt value (AP plot) (GeV/c)"};
    Configurable<float> LambdaMaxQt{"LambdaMaxQt", 0.17, "Max lambda qt value (AP plot) (GeV/c)"};
    Configurable<float> LambdaMinAlpha{"LambdaMinAlpha", 0.25, "Min lambda alpha absolute value (AP plot)"};
    Configurable<float> LambdaMaxAlpha{"LambdaMaxAlpha", 1.0, "Max lambda alpha absolute value (AP plot)"};
    Configurable<float> LambdaMinv0cospa{"LambdaMinv0cospa", 0.95, "Min V0 CosPA"};
    Configurable<float> LambdaMaxLifeTime{"LambdaMaxLifeTime", 30, "Max lifetime"};
    Configurable<float> LambdaWindow{"LambdaWindow", 0.015, "Mass window around expected (in GeV/c2). Leave negative to disable"};
    Configurable<float> LambdaMinRapidity{"LambdaMinRapidity", -0.5, "v0 min rapidity"};
    Configurable<float> LambdaMaxRapidity{"LambdaMaxRapidity", 0.5, "v0 max rapidity"};
    Configurable<float> LambdaDauEtaMin{"LambdaDauEtaMin", -0.8, "Min pseudorapidity of daughter tracks"};
    Configurable<float> LambdaDauEtaMax{"LambdaDauEtaMax", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<float> LambdaMinZ{"LambdaMinZ", -240, "Min lambda decay point z value (cm)"};
    Configurable<float> LambdaMaxZ{"LambdaMaxZ", 240, "Max lambda decay point z value (cm)"};
    Configurable<int> LambdaMinTPCCrossedRows{"LambdaMinTPCCrossedRows", 50, "Min daughter TPC Crossed Rows"};
    Configurable<int> LambdaMinITSclusters{"LambdaMinITSclusters", 1, "minimum ITS clusters"};
    Configurable<bool> LambdaRejectPosITSafterburner{"LambdaRejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
    Configurable<bool> LambdaRejectNegITSafterburner{"LambdaRejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};
  } lambdaSelections;

  //// Photon criteria:
  struct : ConfigurableGroup {
    std::string prefix = "photonSelections"; // JSON group name
    Configurable<float> Gamma_MLThreshold{"Gamma_MLThreshold", 0.1, "Decision Threshold value to select gammas"};
    Configurable<float> doMCAssociation{"doMCAssociation", false, "if MC, select true photons only"};
    Configurable<int> Photonv0TypeSel{"Photonv0TypeSel", 7, "select on a certain V0 type (leave negative if no selection desired)"};
    Configurable<float> PhotonMinDCADauToPv{"PhotonMinDCADauToPv", 0.0, "Min DCA daughter To PV (cm)"};
    Configurable<float> PhotonMaxDCAV0Dau{"PhotonMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
    Configurable<int> PhotonMinTPCCrossedRows{"PhotonMinTPCCrossedRows", 30, "Min daughter TPC Crossed Rows"};
    Configurable<float> PhotonMinTPCNSigmas{"PhotonMinTPCNSigmas", -7, "Min TPC NSigmas for daughters"};
    Configurable<float> PhotonMaxTPCNSigmas{"PhotonMaxTPCNSigmas", 7, "Max TPC NSigmas for daughters"};
    Configurable<float> PhotonMinRapidity{"PhotonMinRapidity", -0.5, "v0 min rapidity"};
    Configurable<float> PhotonMaxRapidity{"PhotonMaxRapidity", 0.5, "v0 max rapidity"};
    Configurable<float> PhotonDauEtaMin{"PhotonDauEtaMin", -0.8, "Min pseudorapidity of daughter tracks"};
    Configurable<float> PhotonDauEtaMax{"PhotonDauEtaMax", 0.8, "Max pseudorapidity of daughter tracks"};
    Configurable<float> PhotonMinRadius{"PhotonMinRadius", 3.0, "Min photon conversion radius (cm)"};
    Configurable<float> PhotonMaxRadius{"PhotonMaxRadius", 115, "Max photon conversion radius (cm)"};
    Configurable<float> PhotonMinZ{"PhotonMinZ", -240, "Min photon conversion point z value (cm)"};
    Configurable<float> PhotonMaxZ{"PhotonMaxZ", 240, "Max photon conversion point z value (cm)"};
    Configurable<float> PhotonMaxQt{"PhotonMaxQt", 0.08, "Max photon qt value (AP plot) (GeV/c)"};
    Configurable<float> PhotonMaxAlpha{"PhotonMaxAlpha", 1.0, "Max photon alpha absolute value (AP plot)"};
    Configurable<float> PhotonMinV0cospa{"PhotonMinV0cospa", 0.80, "Min V0 CosPA"};
    Configurable<float> PhotonMaxMass{"PhotonMaxMass", 0.10, "Max photon mass (GeV/c^{2})"};
    Configurable<float> PhotonPhiMin1{"PhotonPhiMin1", -1, "Phi min value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> PhotonPhiMax1{"PhotonPhiMax1", -1, "Phi max value to reject photons, region 1 (leave negative if no selection desired)"};
    Configurable<float> PhotonPhiMin2{"PhotonPhiMin2", -1, "Phi max value to reject photons, region 2 (leave negative if no selection desired)"};
    Configurable<float> PhotonPhiMax2{"PhotonPhiMax2", -1, "Phi min value to reject photons, region 2 (leave negative if no selection desired)"};
  } photonSelections;

  /// KShort criteria:
  Configurable<float> V0Rapidity{"V0Rapidity", 0.5, "v0 rapidity"};
  Configurable<float> K0ShortDauPseudoRap{"K0SDauPseudoRap", 1.5, "Max pseudorapidity of daughter tracks"};
  Configurable<float> K0ShortMinDCANegToPv{"K0SMinDCANegToPv", 0.0, "min DCA Neg To PV (cm)"};
  Configurable<float> K0ShortMinDCAPosToPv{"K0SMinDCAPosToPv", 0.0, "min DCA Pos To PV (cm)"};
  Configurable<float> K0ShortMaxDCAV0Dau{"K0SMaxDCAV0Dau", 3.5, "Max DCA V0 Daughters (cm)"};
  Configurable<float> K0ShortMinv0radius{"K0SMinv0radius", 0.0, "Min V0 radius (cm)"};
  Configurable<float> K0ShortMaxv0radius{"K0SMaxv0radius", 60, "Max V0 radius (cm)"};
  Configurable<float> K0ShortWindow{"K0SWindow", 0.1, "Mass window around expected (in GeV/c2)"};

  // KStar criteria: 
  Configurable<float> KStarWindow{"KStarWindow", 0.1, "Mass window around expected (in GeV/c2)"};
  Configurable<float> KStarMaxRap{"KStarMaxRap", 0.8, "Max kstar rapidity"};

  //// Sigma0 criteria:
  Configurable<float> Sigma0Window{"Sigma0Window", 0.1, "Mass window around expected (in GeV/c2)"};
  Configurable<float> SigmaMaxRap{"SigmaMaxRap", 0.8, "Max sigma0 rapidity"};

  //// Pi0 criteria::
  Configurable<float> Pi0MaxRap{"Pi0MaxRap", 0.8, "Max Pi0 Rapidity"};
  Configurable<float> Pi0MassWindow{"Pi0MassWindow", 0.115, "Mass window around expected (in GeV/c2)"};

  //// Generated particles criteria:
  struct : ConfigurableGroup {
    std::string prefix = "genSelections"; // JSON group name
    Configurable<bool> doQA{"doQA", true, "If True, fill QA histos"};
    Configurable<bool> mc_keepOnlyFromGenerator{"mc_keepOnlyFromGenerator", false, "Keep only mcparticles from the generator"};
    Configurable<bool> mc_keepOnlyFromTransport{"mc_keepOnlyFromTransport", false, "Keep only mcparticles from the transport code"};
    Configurable<int> mc_selectMCProcess{"mc_selectMCProcess", -1, "Keep only mcparticles produced in the selected MC process"};
    Configurable<float> mc_rapidityMin{"mc_rapidityMin", -0.5, "Min generated particle rapidity"};
    Configurable<float> mc_rapidityMax{"mc_rapidityMax", 0.5, "Max generated particle rapidity"};
  } genSelections;

  // Axis
  // base properties
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Centrality"};
  ConfigurableAxis axisNch{"axisNch", {300, 0.0f, 3000.0f}, "N_{ch}"};

  // Invariant Mass
  ConfigurableAxis axisSigmaMass{"axisSigmaMass", {500, 1.10f, 1.30f}, "M_{#Sigma^{0}} (GeV/c^{2})"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, "M_{#Lambda} (GeV/c^{2})"};
  ConfigurableAxis axisPhotonMass{"axisPhotonMass", {200, 0.0f, 0.3f}, "M_{#Gamma}"};
  ConfigurableAxis axisK0SMass{"axisK0SMass", {200, 0.4f, 0.6f}, "M_{K^{0}}"};
  ConfigurableAxis axisKStarMass{"axisKStarMass", {500, 0.6f, 1.6f}, "M_{K^{*}} (GeV/c^{2})"};
  // AP plot axes
  ConfigurableAxis axisAPAlpha{"axisAPAlpha", {220, -1.1f, 1.1f}, "V0 AP alpha"};
  ConfigurableAxis axisAPQt{"axisAPQt", {220, 0.0f, 0.5f}, "V0 AP alpha"};

  // topological variable QA axes
  ConfigurableAxis axisTPCrows{"axisTPCrows", {160, 0.0f, 160.0f}, "N TPC rows"};
  ConfigurableAxis axisNCls{"axisNCls", {8, -0.5, 7.5}, "NCls"};
  ConfigurableAxis axisTPCNSigma{"axisTPCNSigma", {40, -10, 10}, "TPC NSigma"};
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {500, 0.0f, 50.0f}, "DCA (cm)"};
  ConfigurableAxis axisXY{"axisXY", {120, -120.0f, 120.0f}, "XY axis"};
  ConfigurableAxis axisZ{"axisZ", {120, -120.0f, 120.0f}, "V0 Z position (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {50, 0.0f, 5.0f}, "DCA (cm)"};
  ConfigurableAxis axisCosPA{"axisCosPA", {200, 0.5f, 1.0f}, "Cosine of pointing angle"};
  ConfigurableAxis axisRadius{"axisRadius", {240, 0.0f, 120.0f}, "V0 radius (cm)"};
  ConfigurableAxis axisPhi{"axisPhi", {200, 0, 2 * o2::constants::math::PI}, "Phi for photons"};
  ConfigurableAxis axisPA{"axisPA", {100, 0.0f, 1}, "Pointing angle"};
  ConfigurableAxis axisRapidity{"axisRapidity", {100, -2.0f, 2.0f}, "Rapidity"};
  ConfigurableAxis axisCandSel{"axisCandSel", {15, 0.5f, +15.5f}, "Candidate Selection"};
  ConfigurableAxis axisIRBinning{"axisIRBinning", {151, -10, 1500}, "Binning for the interaction rate (kHz)"};

  // For manual sliceBy (necessary to calculate the correction factors)
  PresliceUnsorted<soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraCollLabels>> perMcCollision = aod::v0data::straMCCollisionId;

  void init(InitContext const&)
  {
    LOGF(info, "Initializing now: cross-checking correctness...");
    if (doprocessRealData +
          doprocessRealDataWithTOF +
          doprocessMonteCarlo +
          doprocessMonteCarloWithTOF +
          doprocessPhotonLambdaQA +
          doprocessPhotonLambdaMCQA >
        1) {
      LOGF(fatal, "You have enabled more than one process function. Please check your configuration! Aborting now.");
    }

    // setting CCDB service
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    histos.add("hEventCentrality", "hEventCentrality", kTH1D, {axisCentrality});

    if (eventSelections.fUseEventSelection) {
      histos.add("hEventSelection", "hEventSelection", kTH1D, {{21, -0.5f, +20.5f}});
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

      if (fGetIR) {
        histos.add("GeneralQA/hRunNumberNegativeIR", "", kTH1D, {{1, 0., 1.}});
        histos.add("GeneralQA/hInteractionRate", "hInteractionRate", kTH1D, {axisIRBinning});
        histos.add("GeneralQA/hCentralityVsInteractionRate", "hCentralityVsInteractionRate", kTH2D, {axisCentrality, axisIRBinning});
      }
    }

    histos.add("K0ShortSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
    histos.get<TH1>(HIST("K0ShortSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("K0ShortSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "K0Short Mass Cut");
    histos.get<TH1>(HIST("K0ShortSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "K0Short Eta/Y Cut");
    histos.get<TH1>(HIST("K0ShortSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(4, "K0Short DCAToPV Cut");
    histos.get<TH1>(HIST("K0ShortSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(5, "K0Short Radius Cut");
    histos.get<TH1>(HIST("K0ShortSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(6, "K0Short DCADau Cut");

    histos.add("K0ShortSel/hK0ShortMass", "hK0ShortMass", kTH1F, {axisK0SMass});
    histos.add("K0ShortSel/hK0ShortNegEta", "hK0ShortNegEta", kTH1F, {axisRapidity});
    histos.add("K0ShortSel/hK0ShortPosEta", "hK0ShortPosEta", kTH1F, {axisRapidity});
    histos.add("K0ShortSel/hK0ShortY", "hK0ShortY", kTH1F, {axisRapidity});
    histos.add("K0ShortSel/hK0ShortDCANegToPV", "hK0ShortDCANegToPV", kTH1F, {axisDCAtoPV});
    histos.add("K0ShortSel/hK0ShortDCAPosToPV", "hK0ShortDCAPosToPV", kTH1F, {axisDCAtoPV});
    histos.add("K0ShortSel/hK0ShortDCADau", "hK0ShortDCADau", kTH1F, {axisDCAdau});
    histos.add("K0ShortSel/hK0ShortRadius", "hK0ShortRadius", kTH1F, {axisRadius});
    histos.add("K0ShortSel/h3dK0ShortMass", "h3dK0ShortMass", kTH3D, {axisCentrality, axisPt, axisK0SMass});

    histos.add("KStarSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
    histos.get<TH1>(HIST("KStarSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("KStarSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "KStar Mass Window");
    histos.get<TH1>(HIST("KStarSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "KStar Y Window");
    
    histos.add("KStarSel/hKStarMassSelected", "hKStarMassSelected", kTH1F, {axisKStarMass});

    for (const auto& histodir : DirList) {
      if ((histodir == "V0BeforeSel" && !fFillNoSelV0Histos) ||
          (histodir == "PhotonSel" && !fFillSelPhotonHistos) ||
          (histodir == "LambdaSel" && !fFillSelLambdaHistos)) {
        continue;
      }

      histos.add(histodir + "/hpT", "hpT", kTH1D, {axisPt});
      histos.add(histodir + "/hV0Type", "hV0Type", kTH1D, {{8, 0.5f, 8.5f}});
      histos.add(histodir + "/hNegEta", "hNegEta", kTH1D, {axisRapidity});
      histos.add(histodir + "/hPosEta", "hPosEta", kTH1D, {axisRapidity});
      histos.add(histodir + "/hDCANegToPV", "hDCANegToPV", kTH1D, {axisDCAtoPV});
      histos.add(histodir + "/hDCAPosToPV", "hDCAPosToPV", kTH1D, {axisDCAtoPV});
      histos.add(histodir + "/hDCADau", "hnDCADau", kTH1D, {axisDCAdau});
      histos.add(histodir + "/hRadius", "hnRadius", kTH1D, {axisRadius});
      histos.add(histodir + "/hZ", "hZ", kTH1D, {axisZ});
      histos.add(histodir + "/hCosPA", "hCosPA", kTH1D, {axisCosPA});
      histos.add(histodir + "/hPhi", "hPhi", kTH1D, {axisPhi});
      histos.add(histodir + "/hPosTPCCR", "hPosTPCCR", kTH1D, {axisTPCrows});
      histos.add(histodir + "/hNegTPCCR", "hNegTPCCR", kTH1D, {axisTPCrows});
      histos.add(histodir + "/hPosITSNCls", "hPosITSNCls", kTH1D, {axisNCls});
      histos.add(histodir + "/hNegITSNCls", "hNegITSNCls", kTH1D, {axisNCls});
      histos.add(histodir + "/hPosTPCNSigmaEl", "hPosTPCNSigmaEl", kTH1D, {axisTPCNSigma});
      histos.add(histodir + "/hNegTPCNSigmaEl", "hNegTPCNSigmaEl", kTH1D, {axisTPCNSigma});
      histos.add(histodir + "/hPosTPCNSigmaPi", "hPosTPCNSigmaPi", kTH1D, {axisTPCNSigma});
      histos.add(histodir + "/hNegTPCNSigmaPi", "hNegTPCNSigmaPi", kTH1D, {axisTPCNSigma});
      histos.add(histodir + "/hPosTPCNSigmaPr", "hPosTPCNSigmaPr", kTH1D, {axisTPCNSigma});
      histos.add(histodir + "/hNegTPCNSigmaPr", "hNegTPCNSigmaPr", kTH1D, {axisTPCNSigma});
      histos.add(histodir + "/h2dArmenteros", "h2dArmenteros", kTH2D, {axisAPAlpha, axisAPQt});

      histos.add(histodir + "/hPhotonY", "hPhotonY", kTH1D, {axisRapidity});
      histos.add(histodir + "/hPhotonMass", "hPhotonMass", kTH1D, {axisPhotonMass});
      histos.add(histodir + "/h2dMassPhotonVsK0S", "h2dMassPhotonVsK0S", kTH2D, {axisPhotonMass, axisK0SMass});
      histos.add(histodir + "/h2dMassPhotonVsLambda", "h2dMassPhotonVsLambda", kTH2D, {axisPhotonMass, axisLambdaMass});
      histos.add(histodir + "/hLifeTime", "hLifeTime", kTH1D, {axisRapidity});
      histos.add(histodir + "/hLambdaY", "hLambdaY", kTH1D, {axisRapidity});
      histos.add(histodir + "/hLambdaMass", "hLambdaMass", kTH1D, {axisLambdaMass});
      histos.add(histodir + "/hALambdaMass", "hALambdaMass", kTH1D, {axisLambdaMass});
      histos.add(histodir + "/h2dMassLambdaVsK0S", "h2dMassLambdaVsK0S", kTH2D, {axisLambdaMass, axisK0SMass});
      histos.add(histodir + "/h2dMassLambdaVsGamma", "h2dMassLambdaVsGamma", kTH2D, {axisLambdaMass, axisPhotonMass});

      if (histodir != "V0BeforeSel" && fFillV03DPositionHistos) // We dont want this for all reco v0s!
        histos.add(histodir + "/h3dV0XYZ", "h3dV0XYZ", kTH3D, {axisXY, axisXY, axisZ});
    }

    histos.add("PhotonSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "Mass");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "Y");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(4, "Neg Eta");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(5, "Pos Eta");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(6, "DCAToPV");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(7, "DCADau");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(8, "Radius");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(9, "Z");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(10, "CosPA");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(11, "Phi");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(12, "TPCCR");
    histos.get<TH1>(HIST("PhotonSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(13, "TPC NSigma");

    histos.add("LambdaSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "Mass");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "Y");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(4, "Neg Eta");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(5, "Pos Eta");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(6, "DCAToPV");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(7, "Radius");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(8, "Z");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(9, "DCADau");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(10, "Armenteros");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(11, "CosPA");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(12, "TPCCR");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(13, "ITSNCls");
    histos.get<TH1>(HIST("LambdaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(14, "Lifetime");

    if (doprocessRealData || doprocessRealDataWithTOF || doprocessMonteCarlo || doprocessMonteCarloWithTOF) {
      histos.add("SigmaSel/hSigma0DauDeltaIndex", "hSigma0DauDeltaIndex", kTH1F, {{100, -49.5f, 50.5f}});
      histos.add("SigmaSel/hSelectionStatistics", "hSelectionStatistics", kTH1D, {axisCandSel});
      histos.get<TH1>(HIST("SigmaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(1, "No Sel");
      histos.get<TH1>(HIST("SigmaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(2, "Sigma Mass Window");
      histos.get<TH1>(HIST("SigmaSel/hSelectionStatistics"))->GetXaxis()->SetBinLabel(3, "Sigma Y Window");

      histos.add("SigmaSel/hSigmaMassSelected", "hSigmaMassSelected", kTH1F, {axisSigmaMass});
    }

    if (doAssocStudy && (doprocessMonteCarlo || doprocessMonteCarloWithTOF)) {
      histos.add("V0AssoQA/h2dIRVsPt_TrueGamma", "h2dIRVsPt_TrueGamma", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueGamma", "h3dPAVsIRVsPt_TrueGamma", kTH3F, {axisPA, axisIRBinning, axisPt});
      histos.add("V0AssoQA/h2dIRVsPt_TrueGamma_BadCollAssig", "h2dIRVsPt_TrueGamma_BadCollAssig", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueGamma_BadCollAssig", "h3dPAVsIRVsPt_TrueGamma_BadCollAssig", kTH3F, {axisPA, axisIRBinning, axisPt});

      histos.add("V0AssoQA/h2dIRVsPt_TrueLambda", "h2dIRVsPt_TrueLambda", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueLambda", "h3dPAVsIRVsPt_TrueLambda", kTH3F, {axisPA, axisIRBinning, axisPt});
      histos.add("V0AssoQA/h2dIRVsPt_TrueLambda_BadCollAssig", "h2dIRVsPt_TrueLambda_BadCollAssig", kTH2F, {axisIRBinning, axisPt});
      histos.add("V0AssoQA/h3dPAVsIRVsPt_TrueLambda_BadCollAssig", "h3dPAVsIRVsPt_TrueLambda_BadCollAssig", kTH3F, {axisPA, axisIRBinning, axisPt});
    }

    // MC
    if (doprocessMonteCarlo || doprocessMonteCarloWithTOF) {
      histos.add("MCQA/h2dPhotonNMothersVsPDG", "h2dPhotonNMothersVsPDG", kTHnSparseD, {{10, -0.5f, +9.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("MCQA/h2dPhotonNMothersVsMCProcess", "h2dPhotonNMothersVsMCProcess", kTH2D, {{10, -0.5f, +9.5f}, {50, -0.5f, 49.5f}});
      histos.add("MCQA/hPhotonMotherSize", "hPhotonMotherSize", kTH1D, {{10, -0.5f, +9.5f}});
      histos.add("MCQA/hPhotonMCProcess", "hPhotonMCProcess", kTH1D, {{50, -0.5f, 49.5f}});
      histos.add("MCQA/hPhotonMotherMCProcess", "hPhotonMotherMCProcess", kTH1D, {{50, -0.5f, 49.5f}});
      histos.add("MCQA/hLambdaMotherSize", "hLambdaMotherSize", kTH1D, {{10, -0.5f, +9.5f}});
      histos.add("MCQA/hLambdaMCProcess", "hLambdaMCProcess", kTH1D, {{50, -0.5f, 49.5f}});
      histos.add("MCQA/hLambdaMotherMCProcess", "hLambdaMotherMCProcess", kTH1D, {{50, -0.5f, 49.5f}});
      histos.add("MCQA/hSigma0MCCheck", "hSigma0MCCheck", kTH1D, {{4, -0.5f, +3.5f}});
      histos.add("MCQA/hNoV0MCCores", "hNoV0MCCores", kTH1D, {{4, -0.5f, +3.5f}});
    }

    if (doprocessGeneratedRun3 && genSelections.doQA) {

      // Pi0s
      histos.add("GenQA/hGenPi0", "hGenPi0", kTH1D, {axisPt});

      auto hPrimaryPi0s = histos.add<TH1>("GenQA/hPrimaryPi0s", "hPrimaryPi0s", kTH1D, {{2, -0.5f, 1.5f}});
      hPrimaryPi0s->GetXaxis()->SetBinLabel(1, "All Pi0s");
      hPrimaryPi0s->GetXaxis()->SetBinLabel(2, "Primary Pi0s");

      histos.add("GenQA/h2dPi0MCSourceVsPDGMother", "h2dPi0MCSourceVsPDGMother", kTHnSparseD, {{2, -0.5f, 1.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("GenQA/h2dPi0NDaughtersVsPDG", "h2dPi0NDaughtersVsPDG", kTHnSparseD, {{10, -0.5f, +9.5f}, {10001, -5000.5f, +5000.5f}});

      auto h2DGenPi0TypeVsProducedByGen = histos.add<TH2>("GenQA/h2DGenPi0TypeVsProducedByGen", "h2DGenPi0TypeVsProducedByGen", kTH2D, {{2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
      h2DGenPi0TypeVsProducedByGen->GetXaxis()->SetBinLabel(1, "Sterile");
      h2DGenPi0TypeVsProducedByGen->GetXaxis()->SetBinLabel(2, "Non-Sterile");
      h2DGenPi0TypeVsProducedByGen->GetYaxis()->SetBinLabel(1, "Generator");
      h2DGenPi0TypeVsProducedByGen->GetYaxis()->SetBinLabel(2, "Transport");

      // ______________________________________________________
      // Sigma0s
      histos.add("GenQA/hGenSigma0", "hGenSigma0", kTH1D, {axisPt});
      histos.add("GenQA/hGenAntiSigma0", "hGenAntiSigma0", kTH1D, {axisPt});

      histos.add("GenQA/h3dGenSigma0_pTMap", "h3dGenSigma0_pTMap", kTH3D, {axisPt, axisPt, axisPt});
      histos.add("GenQA/h3dGenASigma0_pTMap", "h3dGenASigma0_pTMap", kTH3D, {axisPt, axisPt, axisPt});

      histos.add("GenQA/h2dGenSigma0xy_Generator", "hGenSigma0xy_Generator", kTH2D, {axisXY, axisXY});
      histos.add("GenQA/h2dGenSigma0xy_Transport", "hGenSigma0xy_Transport", kTH2D, {axisXY, axisXY});
      histos.add("GenQA/hGenSigma0Radius_Generator", "hGenSigma0Radius_Generator", kTH1D, {axisRadius});
      histos.add("GenQA/hGenSigma0Radius_Transport", "hGenSigma0Radius_Transport", kTH1D, {axisRadius});

      histos.add("GenQA/h2dSigma0MCSourceVsPDGMother", "h2dSigma0MCSourceVsPDGMother", kTHnSparseD, {{2, -0.5f, 1.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("GenQA/h2dSigma0NDaughtersVsPDG", "h2dSigma0NDaughtersVsPDG", kTHnSparseD, {{10, -0.5f, +9.5f}, {10001, -5000.5f, +5000.5f}});

      auto hPrimarySigma0s = histos.add<TH1>("GenQA/hPrimarySigma0s", "hPrimarySigma0s", kTH1D, {{2, -0.5f, 1.5f}});
      hPrimarySigma0s->GetXaxis()->SetBinLabel(1, "All Sigma0s");
      hPrimarySigma0s->GetXaxis()->SetBinLabel(2, "Primary Sigma0s");

      auto hGenSpecies = histos.add<TH1>("GenQA/hGenSpecies", "hGenSpecies", kTH1D, {{4, -0.5f, 3.5f}});
      hGenSpecies->GetXaxis()->SetBinLabel(1, "All Prim. Lambda");
      hGenSpecies->GetXaxis()->SetBinLabel(2, "All Prim. ALambda");
      hGenSpecies->GetXaxis()->SetBinLabel(5, "All Sigma0s");
      hGenSpecies->GetXaxis()->SetBinLabel(6, "All ASigma0s");

      histos.add("GenQA/hSigma0NDau", "hSigma0NDau", kTH1D, {{10, -0.5f, +9.5f}});
      histos.add("GenQA/h2dSigma0NDauVsProcess", "h2dSigma0NDauVsProcess", kTH2D, {{10, -0.5f, +9.5f}, {50, -0.5f, 49.5f}});

      auto h2DGenSigma0TypeVsProducedByGen = histos.add<TH2>("GenQA/h2DGenSigma0TypeVsProducedByGen", "h2DGenSigma0TypeVsProducedByGen", kTH2D, {{2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
      h2DGenSigma0TypeVsProducedByGen->GetXaxis()->SetBinLabel(1, "Sterile");
      h2DGenSigma0TypeVsProducedByGen->GetXaxis()->SetBinLabel(2, "Non-Sterile");
      h2DGenSigma0TypeVsProducedByGen->GetYaxis()->SetBinLabel(1, "Generator");
      h2DGenSigma0TypeVsProducedByGen->GetYaxis()->SetBinLabel(2, "Transport");

      // ______________________________________________________
      // KStar
      histos.add("GenQA/hGenKStar", "hGenKStar", kTH1D, {axisPt});

      histos.add("GenQA/h2dGenKStarxy_Generator", "hGenKStarxy_Generator", kTH2D, {axisXY, axisXY});
      histos.add("GenQA/h2dGenKStarxy_Transport", "hGenKStarxy_Transport", kTH2D, {axisXY, axisXY});
      histos.add("GenQA/hGenKStarRadius_Generator", "hGenKStarRadius_Generator", kTH1D, {axisRadius});
      histos.add("GenQA/hGenKStarRadius_Transport", "hGenKStarRadius_Transport", kTH1D, {axisRadius});

      histos.add("GenQA/h2dKStarMCSourceVsPDGMother", "h2dKStarMCSourceVsPDGMother", kTHnSparseD, {{2, -0.5f, 1.5f}, {10001, -5000.5f, +5000.5f}});
      histos.add("GenQA/h2dKStarNDaughtersVsPDG", "h2dKStarNDaughtersVsPDG", kTHnSparseD, {{10, -0.5f, +9.5f}, {10001, -5000.5f, +5000.5f}});

      auto hPrimaryKStars = histos.add<TH1>("GenQA/hPrimaryKStars", "hPrimaryKStars", kTH1D, {{2, -0.5f, 1.5f}});
      hPrimaryKStars->GetXaxis()->SetBinLabel(1, "All KStars");
      hPrimaryKStars->GetXaxis()->SetBinLabel(2, "Primary KStars");

      auto hGenSpeciesKStar = histos.add<TH1>("GenQA/hGenSpeciesKStar", "hGenSpeciesKStar", kTH1D, {{4, -0.5f, 3.5f}});
      hGenSpeciesKStar->GetXaxis()->SetBinLabel(1, "All Prim. KShort");
      hGenSpeciesKStar->GetXaxis()->SetBinLabel(5, "All KStars");

      histos.add("GenQA/hKStarNDau", "hKStarNDau", kTH1D, {{10, -0.5f, +9.5f}});
      histos.add("GenQA/h2dKStarNDauVsProcess", "h2dKStarNDauVsProcess", kTH2D, {{10, -0.5f, +9.5f}, {50, -0.5f, 49.5f}});

      auto h2DGenKStarTypeVsProducedByGen = histos.add<TH2>("GenQA/h2DGenKStarTypeVsProducedByGen", "h2DGenKStarTypeVsProducedByGen", kTH2D, {{2, -0.5f, 1.5f}, {2, -0.5f, 1.5f}});
      h2DGenKStarTypeVsProducedByGen->GetXaxis()->SetBinLabel(1, "Sterile");
      h2DGenKStarTypeVsProducedByGen->GetXaxis()->SetBinLabel(2, "Non-Sterile");
      h2DGenKStarTypeVsProducedByGen->GetYaxis()->SetBinLabel(1, "Generator");
      h2DGenKStarTypeVsProducedByGen->GetYaxis()->SetBinLabel(2, "Transport");


    }

    if (doprocessPhotonLambdaQA || doprocessPhotonLambdaMCQA) {

      // Event selection:
      histos.add("PhotonLambdaQA/hEventCentrality", "hEventCentrality", kTH1D, {axisCentrality});

      // Photon part:
      histos.add("PhotonLambdaQA/h3dPhotonMass", "h3dPhotonMass", kTH3D, {axisCentrality, axisPt, axisPhotonMass});
      histos.add("PhotonLambdaQA/h3dYPhotonMass", "h3dYPhotonMass", kTH3D, {axisRapidity, axisPt, axisPhotonMass});
      histos.add("PhotonLambdaQA/h3dYPhotonRadius", "h3dYPhotonRadius", kTH3D, {axisRapidity, axisPt, axisRadius});

      histos.add("PhotonLambdaQA/h3dTruePhotonMass", "h3dTruePhotonMass", kTH3D, {axisCentrality, axisPt, axisPhotonMass});
      histos.add("PhotonLambdaQA/h2dTrueSigma0PhotonMass", "h2dTrueSigma0PhotonMass", kTH2D, {axisPt, axisPhotonMass});

      // Lambda part:
      histos.add("PhotonLambdaQA/h3dLambdaMass", "h3dLambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});
      histos.add("PhotonLambdaQA/h3dTrueLambdaMass", "h3dTrueLambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});
      histos.add("PhotonLambdaQA/h3dYLambdaMass", "h3dYLambdaMass", kTH3D, {axisRapidity, axisPt, axisLambdaMass});
      histos.add("PhotonLambdaQA/h3dYRLambdaMass", "h3dYRLambdaMass", kTH3D, {axisRapidity, axisRadius, axisLambdaMass});

      histos.add("PhotonLambdaQA/h2dTrueSigma0LambdaMass", "h2dTrueSigma0LambdaMass", kTH2D, {axisPt, axisLambdaMass});

      // AntiLambda part:
      histos.add("PhotonLambdaQA/h3dALambdaMass", "h3dALambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});
      histos.add("PhotonLambdaQA/h3dTrueALambdaMass", "h3dTrueALambdaMass", kTH3D, {axisCentrality, axisPt, axisLambdaMass});
      histos.add("PhotonLambdaQA/h3dYALambdaMass", "h3dYALambdaMass", kTH3D, {axisRapidity, axisPt, axisLambdaMass});
      histos.add("PhotonLambdaQA/h3dYRALambdaMass", "h3dYRALambdaMass", kTH3D, {axisRapidity, axisRadius, axisLambdaMass});

      histos.add("PhotonLambdaQA/h2dTrueASigma0ALambdaMass", "h2dTrueASigma0ALambdaMass", kTH2D, {axisPt, axisLambdaMass});
    }

    if (doprocessPhotonLambdaGenerated) {

      histos.add("PhotonLambdaQA/hGenEvents", "hGenEvents", kTH2D, {{axisNch}, {2, -0.5f, +1.5f}});
      histos.get<TH2>(HIST("PhotonLambdaQA/hGenEvents"))->GetYaxis()->SetBinLabel(1, "All gen. events");
      histos.get<TH2>(HIST("PhotonLambdaQA/hGenEvents"))->GetYaxis()->SetBinLabel(2, "Gen. with at least 1 rec. events");
      histos.add("PhotonLambdaQA/hGenEventCentrality", "hGenEventCentrality", kTH1D, {{101, 0.0f, 101.0f}});

      histos.add("PhotonLambdaQA/hCentralityVsNcoll_beforeEvSel", "hCentralityVsNcoll_beforeEvSel", kTH2D, {axisCentrality, {50, -0.5f, 49.5f}});
      histos.add("PhotonLambdaQA/hCentralityVsNcoll_afterEvSel", "hCentralityVsNcoll_afterEvSel", kTH2D, {axisCentrality, {50, -0.5f, 49.5f}});

      histos.add("PhotonLambdaQA/hCentralityVsMultMC", "hCentralityVsMultMC", kTH2D, {{101, 0.0f, 101.0f}, axisNch});
      histos.add("PhotonLambdaQA/hEventPVzMC", "hEventPVzMC", kTH1D, {{100, -20.0f, +20.0f}});
      histos.add("PhotonLambdaQA/hCentralityVsPVzMC", "hCentralityVsPVzMC", kTH2D, {{101, 0.0f, 101.0f}, {100, -20.0f, +20.0f}});

      histos.add("PhotonLambdaQA/h2dGenPhoton", "h2dGenPhoton", kTH2D, {axisCentrality, axisPt});
      histos.add("PhotonLambdaQA/h2dGenLambda", "h2dGenLambda", kTH2D, {axisCentrality, axisPt});
      histos.add("PhotonLambdaQA/h2dGenAntiLambda", "h2dGenAntiLambda", kTH2D, {axisCentrality, axisPt});

      histos.add("PhotonLambdaQA/h2dGenPhotonVsMultMC_RecoedEvt", "h2dGenPhotonVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
      histos.add("PhotonLambdaQA/h2dGenLambdaVsMultMC_RecoedEvt", "h2dGenLambdaVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});
      histos.add("PhotonLambdaQA/h2dGenAntiLambdaVsMultMC_RecoedEvt", "h2dGenAntiLambdaVsMultMC_RecoedEvt", kTH2D, {axisNch, axisPt});

      histos.add("PhotonLambdaQA/h2dGenPhotonVsMultMC", "h2dGenPhotonVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("PhotonLambdaQA/h2dGenLambdaVsMultMC", "h2dGenLambdaVsMultMC", kTH2D, {axisNch, axisPt});
      histos.add("PhotonLambdaQA/h2dGenAntiLambdaVsMultMC", "h2dGenAntiLambdaVsMultMC", kTH2D, {axisNch, axisPt});
    }

    // inspect histogram sizes, please
    histos.print();
  }

  // ______________________________________________________
  // Struct to store V0Pair properties
  struct V0PairTopoInfo {
    float X = -999.f;
    float Y = -999.f;
    float Z = -999.f;
    float DCADau = -999.f;
    float CosPA = -1.f;
  };

  // ______________________________________________________
  // Struct to store V0Pair MC properties
  struct V0PairMCInfo {
    bool fIsV01CorrectlyAssign = false;
    bool fIsV02CorrectlyAssign = false;
    bool fIsV01Primary = false;
    bool fIsV02Primary = false;
    bool fV0PairProducedByGenerator = false;
    int V01PDGCode = 0;
    int V02PDGCode = 0;
    int V01PDGCodeMother = 0;
    int V02PDGCodeMother = 0;
    int V0PairPDGCode = 0;
    int V0PairPDGCodeMother = 0;
    int V0PairMCProcess = -1;
    int V0PairMCParticleID = -1;
    float V01MCpx = -999.f;
    float V01MCpy = -999.f;
    float V01MCpz = -999.f;
    float V02MCpx = -999.f;
    float V02MCpy = -999.f;
    float V02MCpz = -999.f;
    float V0PairMCRadius = -999.f;
  };

  // ______________________________________________________
  // Struct to store V0Pair Generated properties
  struct V0PairGenInfo {
    bool IsPrimary = false;
    bool IsV0Lambda = false;
    bool IsV0AntiLambda = false;
    bool IsV0KShort = false;
    bool IsPi0 = false;
    bool IsSigma0 = false;
    bool IsAntiSigma0 = false;
    bool IsKStar = false;
    bool IsProducedByGenerator = false;
    bool IsSterile = false;
    int MCProcess = -1;
    int MCCollId = -1;
    int PDGCodeMother = 0;
    int NDaughters = -1;
    float MCPt = -999.f;
    float MCDau1Pt = -999.f;
    float MCDau2Pt = -999.f;
    float MCvx = 999.f;
    float MCvy = 999.f;
  };

  template <typename TV01, typename TV02>
  V0PairTopoInfo propagateV0PairToDCA(TV01 const& v01, TV02 const& v02)
  {
    V0PairTopoInfo info;

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

    // Calculate properties and fill struct
    info.DCADau = (cross.Mag2() > 0) ? std::abs(posdiff.Dot(cross)) / cross.R() : 999.f;
    info.CosPA = v01momentumNorm.Dot(v02momentumNorm);

    if (d < 1e-5f) {                  // Parallel or nearly parallel lines
      info.X = info.Y = info.Z = 0.f; // should we use another dummy value? Perhaps 999.f?
      return info;
    }

    info.X = PCA.X();
    info.Y = PCA.Y();
    info.Z = PCA.Z();

    return info;
  }

  template <typename TV01, typename TV02, typename TCollision, typename TMCParticles>
  V0PairMCInfo getV0PairMCInfo(TV01 const& v01, TV02 const& v02, TCollision const& collision, TMCParticles const& mcparticles)
  {
    V0PairMCInfo MCinfo;

    if (!v01.has_v0MCCore() || !v02.has_v0MCCore()) {
      histos.fill(HIST("MCQA/hNoV0MCCores"), 1);
      return MCinfo;
    }

    auto v01MC = v01.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
    auto v02MC = v02.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

    // Sanity check: Is V0Pair <-> Mother assignment correct?
    bool fIsSigma0 = false;
    if ((v01MC.pdgCode() == 22) && (v01MC.pdgCodeMother() == 3212) && (v02MC.pdgCode() == 3122) && (v02MC.pdgCodeMother() == 3212) && (v01.motherMCPartId() == v02.motherMCPartId()))
      fIsSigma0 = true;

    // Check collision assignment
    if (collision.has_straMCCollision()) {
      auto MCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
      MCinfo.fIsV01CorrectlyAssign = (v01MC.straMCCollisionId() == MCCollision.globalIndex());
      MCinfo.fIsV02CorrectlyAssign = (v02MC.straMCCollisionId() == MCCollision.globalIndex());
    }

    // Basic kinematic info
    MCinfo.V01MCpx = v01MC.pxMC();
    MCinfo.V01MCpy = v01MC.pyMC();
    MCinfo.V01MCpz = v01MC.pzMC();
    MCinfo.V02MCpx = v02MC.pxMC();
    MCinfo.V02MCpy = v02MC.pyMC();
    MCinfo.V02MCpz = v02MC.pzMC();

    // MC association info
    MCinfo.fIsV01Primary = v01MC.isPhysicalPrimary();
    MCinfo.fIsV02Primary = v02MC.isPhysicalPrimary();
    MCinfo.V01PDGCode = v01MC.pdgCode();
    MCinfo.V02PDGCode = v02MC.pdgCode();
    MCinfo.V01PDGCodeMother = v01MC.pdgCodeMother();
    MCinfo.V02PDGCodeMother = v02MC.pdgCodeMother();

    // Get corresponding entries in MCParticles table
    auto MCParticle_v01 = mcparticles.rawIteratorAt(v01MC.particleIdMC());
    auto MCParticle_v02 = mcparticles.rawIteratorAt(v02MC.particleIdMC());

    // Get MC Mothers
    auto const& MCMothersList_v01 = MCParticle_v01.template mothers_as<aod::McParticles>();
    auto const& MCMothersList_v02 = MCParticle_v02.template mothers_as<aod::McParticles>();

    if (!MCMothersList_v01.empty() && !MCMothersList_v02.empty()) { // Are there mothers?
      auto const& MCMother_v01 = MCMothersList_v01.front(); // First mother
      auto const& MCMother_v02 = MCMothersList_v02.front(); // First mother

      if (MCMother_v01.globalIndex() == MCMother_v02.globalIndex()) { // Is it the same mother?

        MCinfo.fV0PairProducedByGenerator = MCMother_v01.producedByGenerator();
        MCinfo.V0PairPDGCode = MCMother_v01.pdgCode();
        MCinfo.V0PairMCProcess = MCMother_v01.getProcess();
        MCinfo.V0PairMCParticleID = MCMother_v01.globalIndex();
        MCinfo.V0PairMCRadius = std::hypot(MCMother_v01.vx(), MCMother_v01.vy()); // production position radius

        auto const& v0pairmothers = MCMother_v01.template mothers_as<aod::McParticles>(); // Get mothers
        if (!v0pairmothers.empty()) {
          auto& v0PairMother = v0pairmothers.front(); // V0Pair mother, V0s grandmother
          MCinfo.V0PairPDGCodeMother = v0PairMother.pdgCode();
        }

        // MC QA histograms
        // Parenthood check for sigma0-like candidate
        if (MCParticle_v01.pdgCode() == 22 && TMath::Abs(MCParticle_v02.pdgCode()) == 3122) {
          for (const auto& mother1 : MCMothersList_v01) { // Photon mothers
            histos.fill(HIST("MCQA/h2dPhotonNMothersVsPDG"), MCMothersList_v01.size(), mother1.pdgCode());
            histos.fill(HIST("MCQA/h2dPhotonNMothersVsMCProcess"), MCMothersList_v01.size(), mother1.getProcess());

            for (const auto& mother2 : MCMothersList_v02) {         // Lambda mothers
              if (mother1.globalIndex() == mother2.globalIndex()) { // Match found: same physical mother

                if (mother1.globalIndex() == MCMother_v01.globalIndex()) {
                  histos.fill(HIST("MCQA/hPhotonMotherSize"), MCMothersList_v01.size());
                  histos.fill(HIST("MCQA/hPhotonMCProcess"), MCParticle_v01.getProcess());
                  histos.fill(HIST("MCQA/hPhotonMotherMCProcess"), mother1.getProcess());
                }

                if (mother2.globalIndex() == MCMother_v02.globalIndex()) {
                  histos.fill(HIST("MCQA/hLambdaMotherSize"), MCMothersList_v02.size());
                  histos.fill(HIST("MCQA/hLambdaMCProcess"), MCParticle_v02.getProcess());
                  histos.fill(HIST("MCQA/hLambdaMotherMCProcess"), mother2.getProcess());
                }
              }
            }
          }
        }

        // Check association correctness
        if (fIsSigma0 && (MCinfo.V0PairPDGCode == 3212))
          histos.fill(HIST("MCQA/hSigma0MCCheck"), 1); // match
        if (fIsSigma0 && !(MCinfo.V0PairPDGCode == 3212))
          histos.fill(HIST("MCQA/hSigma0MCCheck"), 2); // mismatch
        if (!fIsSigma0 && (MCinfo.V0PairPDGCode == 3212))
          histos.fill(HIST("MCQA/hSigma0MCCheck"), 3); // mismatch
      }
    }

    return MCinfo;
  }

  // ______________________________________________________
  // Check whether the collision passes our collision selections
  // Should work with collisions, mccollisions, stracollisions and stramccollisions tables!
  template <typename TCollision>
  bool IsEventAccepted(TCollision const& collision, bool fillHists)
  {
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 0. /* all collisions */);
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);
    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
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
    } else { // we are in Pb-Pb
      float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
      if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 16 /* Below min occupancy */);
      if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
        return false;
      }
      if (fillHists)
        histos.fill(HIST("hEventSelection"), 17 /* Above max occupancy */);
    }

    // Fetch interaction rate only if required (in order to limit ccdb calls)
    float interactionRate = (fGetIR) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource, fIRCrashOnNull) * 1.e-3 : -1;
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    if (fGetIR) {
      if (interactionRate < 0)
        histos.get<TH1>(HIST("GeneralQA/hRunNumberNegativeIR"))->Fill(Form("%d", collision.runNumber()), 1); // This lists all run numbers without IR info!

      histos.fill(HIST("GeneralQA/hInteractionRate"), interactionRate);
      histos.fill(HIST("GeneralQA/hCentralityVsInteractionRate"), centrality, interactionRate);
    }

    if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 18 /* Below min IR */);

    if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR) {
      return false;
    }
    if (fillHists)
      histos.fill(HIST("hEventSelection"), 19 /* Above max IR */);

    // Fill centrality histogram after event selection
    if (fillHists)
      histos.fill(HIST("hEventCentrality"), centrality);

    return true;
  }

  // ______________________________________________________
  // Simulated processing
  // Return the list of indices to the recoed collision associated to a given MC collision.
  template <typename TMCollisions, typename TCollisions>
  std::vector<int> getListOfRecoCollIndices(TMCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      int biggestNContribs = -1;
      int bestCollisionIndex = -1;
      for (auto const& collision : groupedCollisions) {
        // consider event selections in the recoed <-> gen collision association, for the denominator (or numerator) of the efficiency (or signal loss)?
        if (eventSelections.useEvtSelInDenomEff && eventSelections.fUseEventSelection) {
          if (!IsEventAccepted(collision, false))
            continue;
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
  // Simulated processing
  // Fill generated event information (for event loss/splitting estimation)
  template <typename TMCCollisions, typename TCollisions>
  void fillGeneratedEventProperties(TMCCollisions const& mcCollisions, TCollisions const& collisions)
  {
    std::vector<int> listBestCollisionIdx(mcCollisions.size());
    for (auto const& mcCollision : mcCollisions) {
      // Apply selections on MC collisions
      if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
        continue;
      }
      if (doPPAnalysis) { // we are in pp
        if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
          continue;
        }

        if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
          continue;
        }
      }

      histos.fill(HIST("PhotonLambdaQA/hGenEvents"), mcCollision.multMCNParticlesEta05(), 0 /* all gen. events*/);

      auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCollision.globalIndex());
      // Check if there is at least one of the reconstructed collisions associated to this MC collision
      // If so, we consider it
      bool atLeastOne = false;
      int biggestNContribs = -1;
      float centrality = 100.5f;
      int nCollisions = 0;
      for (auto const& collision : groupedCollisions) {
        if (eventSelections.fUseEventSelection) {
          if (!IsEventAccepted(collision, false))
            continue;
        }
        if (biggestNContribs < collision.multPVTotalContributors()) {
          biggestNContribs = collision.multPVTotalContributors();
          centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();
        }

        nCollisions++;
        atLeastOne = true;
      }

      histos.fill(HIST("PhotonLambdaQA/hCentralityVsNcoll_beforeEvSel"), centrality, groupedCollisions.size());
      histos.fill(HIST("PhotonLambdaQA/hCentralityVsNcoll_afterEvSel"), centrality, nCollisions);
      histos.fill(HIST("PhotonLambdaQA/hCentralityVsMultMC"), centrality, mcCollision.multMCNParticlesEta05());
      histos.fill(HIST("PhotonLambdaQA/hCentralityVsPVzMC"), centrality, mcCollision.posZ());
      histos.fill(HIST("PhotonLambdaQA/hEventPVzMC"), mcCollision.posZ());

      if (atLeastOne) {
        histos.fill(HIST("PhotonLambdaQA/hGenEvents"), mcCollision.multMCNParticlesEta05(), 1 /* at least 1 rec. event*/);
        histos.fill(HIST("PhotonLambdaQA/hGenEventCentrality"), centrality);
      }
    }
    return;
  }

  // ______________________________________________________
  // Auxiliary function to get generated photon and lambda info
  template <typename TMCCollisions, typename TV0MCs, typename TCollisions>
  void runGenPhotonLambdaQA(TMCCollisions const& mcCollisions, TV0MCs const& V0MCCores, TCollisions const& collisions)
  {
    fillGeneratedEventProperties(mcCollisions, collisions);
    std::vector<int> listBestCollisionIdx = getListOfRecoCollIndices(mcCollisions, collisions);
    for (auto const& v0MC : V0MCCores) {
      if (!v0MC.has_straMCCollision())
        continue;

      if (!v0MC.isPhysicalPrimary())
        continue;

      float ptmc = v0MC.ptMC();
      float ymc = 1e3;
      if (v0MC.pdgCode() == PDG_t::kGamma)
        ymc = RecoDecay::y(std::array{v0MC.pxMC(), v0MC.pyMC(), v0MC.pzMC()}, o2::constants::physics::MassGamma);
      else if (std::abs(v0MC.pdgCode()) == PDG_t::kLambda0)
        ymc = v0MC.rapidityMC(1);

      if ((ymc < genSelections.mc_rapidityMin) || (ymc > genSelections.mc_rapidityMax))
        continue;

      auto mcCollision = v0MC.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
      if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
        continue;
      }
      if (doPPAnalysis) { // we are in pp
        if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
          continue;
        }

        if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
          continue;
        }
      }

      float centrality = 100.5f;
      if (listBestCollisionIdx[mcCollision.globalIndex()] > -1) {
        auto collision = collisions.iteratorAt(listBestCollisionIdx[mcCollision.globalIndex()]);

        centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

        if (v0MC.pdgCode() == PDG_t::kGamma) {
          histos.fill(HIST("PhotonLambdaQA/h2dGenPhotonVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (v0MC.pdgCode() == PDG_t::kLambda0) {
          histos.fill(HIST("PhotonLambdaQA/h2dGenLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
        if (v0MC.pdgCode() == PDG_t::kLambda0Bar) {
          histos.fill(HIST("PhotonLambdaQA/h2dGenAntiLambdaVsMultMC_RecoedEvt"), mcCollision.multMCNParticlesEta05(), ptmc);
        }
      }

      if (v0MC.pdgCode() == PDG_t::kGamma) {
        histos.fill(HIST("PhotonLambdaQA/h2dGenPhoton"), centrality, ptmc);
        histos.fill(HIST("PhotonLambdaQA/h2dGenPhotonVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (v0MC.pdgCode() == PDG_t::kLambda0) {
        histos.fill(HIST("PhotonLambdaQA/h2dGenLambda"), centrality, ptmc);
        histos.fill(HIST("PhotonLambdaQA/h2dGenLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
      if (v0MC.pdgCode() == PDG_t::kLambda0Bar) {
        histos.fill(HIST("PhotonLambdaQA/h2dGenAntiLambda"), centrality, ptmc);
        histos.fill(HIST("PhotonLambdaQA/h2dGenAntiLambdaVsMultMC"), mcCollision.multMCNParticlesEta05(), ptmc);
      }
    }
  }

  // ______________________________________________________
  // MC-specific
  // Analyze v0-to-collision association
  template <typename TCollision, typename TV0Object>
  void analyzeV0CollAssoc(TCollision const& collision, TV0Object const& fullv0s, std::vector<int> selV0Indices, bool isPhotonAnalysis)
  {
    auto v0MCCollision = collision.template straMCCollision_as<soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>>();
    float IR = (fGetIR) ? rateFetcher.fetch(ccdb.service, collision.timestamp(), collision.runNumber(), irSource, fIRCrashOnNull) * 1.e-3 : -1;

    for (size_t i = 0; i < selV0Indices.size(); ++i) {
      auto v0 = fullv0s.rawIteratorAt(selV0Indices[i]);
      auto v0MC = v0.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

      float V0MCpT = RecoDecay::pt(array<float, 2>{v0MC.pxMC(), v0MC.pyMC()});
      float V0PA = TMath::ACos(v0.v0cosPA());
      bool fIsV0CorrectlyAssigned = (v0MC.straMCCollisionId() == v0MCCollision.globalIndex());
      bool isPrimary = v0MC.isPhysicalPrimary();

      if ((v0MC.pdgCode() == 22) && isPhotonAnalysis && isPrimary) { // True Gamma
        histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueGamma"), IR, V0MCpT);
        histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueGamma"), V0PA, IR, V0MCpT);

        if (!fIsV0CorrectlyAssigned) {
          histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueGamma_BadCollAssig"), IR, V0MCpT);
          histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueGamma_BadCollAssig"), V0PA, IR, V0MCpT);
        }
      }
      if ((v0MC.pdgCode() == 3122) && !isPhotonAnalysis && isPrimary) { // True Lambda
        histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueLambda"), IR, V0MCpT);
        histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueLambda"), V0PA, IR, V0MCpT);

        if (!fIsV0CorrectlyAssigned) {
          histos.fill(HIST("V0AssoQA/h2dIRVsPt_TrueLambda_BadCollAssig"), IR, V0MCpT);
          histos.fill(HIST("V0AssoQA/h3dPAVsIRVsPt_TrueLambda_BadCollAssig"), V0PA, IR, V0MCpT);
        }
      }
    }
  }

  template <typename TMCParticle>
  V0PairGenInfo getV0PairGenInfo(TMCParticle const& mcParticle)
  {
    V0PairGenInfo GenInfo; // auxiliary struct to store info

    // Fill with properties
    GenInfo.IsPrimary = mcParticle.isPhysicalPrimary();
    GenInfo.IsV0Lambda = mcParticle.pdgCode() == 3122;
    GenInfo.IsV0AntiLambda = mcParticle.pdgCode() == -3122;
    GenInfo.IsV0KShort = mcParticle.pdgCode() == 310;
    GenInfo.IsPi0 = mcParticle.pdgCode() == 111;
    GenInfo.IsSigma0 = mcParticle.pdgCode() == 3212;
    GenInfo.IsAntiSigma0 = mcParticle.pdgCode() == -3212;
    GenInfo.IsKStar = mcParticle.pdgCode() == 313;
    GenInfo.IsProducedByGenerator = mcParticle.producedByGenerator();
    GenInfo.MCProcess = mcParticle.getProcess();
    GenInfo.MCPt = mcParticle.pt();
    GenInfo.MCvx = mcParticle.vx(); // production position X
    GenInfo.MCvy = mcParticle.vy(); // production position Y

    if (mcParticle.has_mcCollision())
      GenInfo.MCCollId = mcParticle.mcCollisionId(); // save this reference, please

    // Checking decay mode if sigma0 or pi0 (it is easier here)
    if (GenInfo.IsSigma0 || GenInfo.IsAntiSigma0 || GenInfo.IsPi0 || GenInfo.IsKStar) {

      // This is a costly operation, so we do it only for pi0s and sigma0s
      auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
      GenInfo.NDaughters = daughters.size();
      GenInfo.IsSterile = daughters.size() == 0;

      auto const& GenMothersList = mcParticle.template mothers_as<aod::McParticles>();
      GenInfo.PDGCodeMother = (!GenMothersList.empty()) ? GenMothersList.front().pdgCode() : 0;

      if ((GenInfo.IsSigma0 || GenInfo.IsAntiSigma0) && genSelections.doQA) {
        histos.fill(HIST("GenQA/h2dSigma0MCSourceVsPDGMother"), GenInfo.IsProducedByGenerator, GenInfo.PDGCodeMother);

        // Checking decay modes and getting daughter pTs
        for (auto& daughter : daughters) {
          histos.fill(HIST("GenQA/h2dSigma0NDaughtersVsPDG"), daughters.size(), daughter.pdgCode());

          if (GenInfo.NDaughters == 2) {
            if (daughter.pdgCode() == 22)
              GenInfo.MCDau1Pt = daughter.pt();

            if (TMath::Abs(daughter.pdgCode()) == 3122)
              GenInfo.MCDau2Pt = daughter.pt();
          }
        }
      }

      if ((GenInfo.IsKStar) && genSelections.doQA) {
        histos.fill(HIST("GenQA/h2dKStarMCSourceVsPDGMother"), GenInfo.IsProducedByGenerator, GenInfo.PDGCodeMother);
        for (auto& daughter : daughters) // checking decay modes
          histos.fill(HIST("GenQA/h2dKStarNDaughtersVsPDG"), daughters.size(), daughter.pdgCode());
      }


      if (GenInfo.IsPi0 && genSelections.doQA) {
        histos.fill(HIST("GenQA/h2dPi0MCSourceVsPDGMother"), GenInfo.IsProducedByGenerator, GenInfo.PDGCodeMother);
        for (auto& daughter : daughters) // checking decay modes
          histos.fill(HIST("GenQA/h2dPi0NDaughtersVsPDG"), daughters.size(), daughter.pdgCode());
      }
    }
    return GenInfo;
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  void fillGenQAHistos(V0PairGenInfo const& GenInfo)
  {
    if (GenInfo.IsPi0) {
      histos.fill(HIST("GenQA/hGenPi0"), GenInfo.MCPt);
      histos.fill(HIST("GenQA/hPrimaryPi0s"), 0);
      if (GenInfo.IsPrimary)
        histos.fill(HIST("GenQA/hPrimaryPi0s"), 1);

      if (GenInfo.IsSterile) {
        if (GenInfo.IsProducedByGenerator)
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 0, 0);
        else
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 0, 1);
      } else {
        if (GenInfo.IsProducedByGenerator)
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 1, 0);
        else
          histos.fill(HIST("GenQA/h2DGenPi0TypeVsProducedByGen"), 1, 1);
      }
    }

    if (GenInfo.IsV0Lambda && GenInfo.IsPrimary)
      histos.fill(HIST("GenQA/hGenSpecies"), 0);
    if (GenInfo.IsV0AntiLambda && GenInfo.IsPrimary)
      histos.fill(HIST("GenQA/hGenSpecies"), 1);
    if (GenInfo.IsV0KShort && GenInfo.IsPrimary)
      histos.fill(HIST("GenQA/hGenSpeciesKStar"), 0);
      

    // Checking decay mode
    if (GenInfo.IsSigma0 || GenInfo.IsAntiSigma0) {
      histos.fill(HIST("GenQA/hSigma0NDau"), GenInfo.NDaughters);
      histos.fill(HIST("GenQA/h2dSigma0NDauVsProcess"), GenInfo.NDaughters, GenInfo.MCProcess);

      const auto radius = std::hypot(GenInfo.MCvx, GenInfo.MCvy);
      // Sigma0 XY and radius (separate histos for Gen/Transport)
      if (GenInfo.IsProducedByGenerator) {
        histos.fill(HIST("GenQA/h2dGenSigma0xy_Generator"), GenInfo.MCvx, GenInfo.MCvy);
        histos.fill(HIST("GenQA/hGenSigma0Radius_Generator"), radius);
      } else {
        histos.fill(HIST("GenQA/h2dGenSigma0xy_Transport"), GenInfo.MCvx, GenInfo.MCvy);
        histos.fill(HIST("GenQA/hGenSigma0Radius_Transport"), radius);
      }

      // Sigma0 type vs origin (single 2D histo)
      const int genIndex = GenInfo.IsProducedByGenerator ? 0 : 1; // 0 = Generator, 1 = Transport
      const int typeIndex = GenInfo.IsSterile ? 0 : 1;            // 0 = Sterile,   1 = Normal
      histos.fill(HIST("GenQA/h2DGenSigma0TypeVsProducedByGen"), typeIndex, genIndex);

      // Fill histograms
      if (GenInfo.IsSigma0) {
        histos.fill(HIST("GenQA/hGenSpecies"), 2);
        histos.fill(HIST("GenQA/hGenSigma0"), GenInfo.MCPt);
        histos.fill(HIST("GenQA/h3dGenSigma0_pTMap"), GenInfo.MCPt, GenInfo.MCDau1Pt, GenInfo.MCDau2Pt);

        histos.fill(HIST("GenQA/hPrimarySigma0s"), 0);
        if (GenInfo.IsPrimary)
          histos.fill(HIST("GenQA/hPrimarySigma0s"), 1);
      }
      if (GenInfo.IsAntiSigma0) {
        histos.fill(HIST("GenQA/hGenSpecies"), 3);
        histos.fill(HIST("GenQA/hGenAntiSigma0"), GenInfo.MCPt);
        histos.fill(HIST("GenQA/h3dGenASigma0_pTMap"), GenInfo.MCPt, GenInfo.MCDau1Pt, GenInfo.MCDau2Pt);
      }
    }

  
    // Checking decay mode
    if (GenInfo.IsKStar) {
      histos.fill(HIST("GenQA/hKStarNDau"), GenInfo.NDaughters);
      histos.fill(HIST("GenQA/h2dKStarNDauVsProcess"), GenInfo.NDaughters, GenInfo.MCProcess);

      const auto radius = std::hypot(GenInfo.MCvx, GenInfo.MCvy);
      // KStar XY and radius (separate histos for Gen/Transport)
      if (GenInfo.IsProducedByGenerator) {
        histos.fill(HIST("GenQA/h2dGenKStarxy_Generator"), GenInfo.MCvx, GenInfo.MCvy);
        histos.fill(HIST("GenQA/hGenKStarRadius_Generator"), radius);
      } else {
        histos.fill(HIST("GenQA/h2dGenKStarxy_Transport"), GenInfo.MCvx, GenInfo.MCvy);
        histos.fill(HIST("GenQA/hGenKStarRadius_Transport"), radius);
      }

      // kstar type vs origin (single 2D histo)
      const int genIndex = GenInfo.IsProducedByGenerator ? 0 : 1; // 0 = Generator, 1 = Transport
      const int typeIndex = GenInfo.IsSterile ? 0 : 1;            // 0 = Sterile,   1 = Normal
      histos.fill(HIST("GenQA/h2DGenKStarTypeVsProducedByGen"), typeIndex, genIndex);

      // Fill histograms
      if (GenInfo.IsKStar) {
        histos.fill(HIST("GenQA/hGenSpeciesKStar"), 2);
        histos.fill(HIST("GenQA/hGenKStar"), GenInfo.MCPt);
        histos.fill(HIST("GenQA/hPrimaryKStars"), 0);
        if (GenInfo.IsPrimary)
          histos.fill(HIST("GenQA/hPrimaryKStars"), 1);
      }
    }
  
  }

  // ______________________________________________________
  // Simulated processing (subscribes to MC information too)
  template <typename TMCParticles>
  void genProcess(TMCParticles const& mcParticles)
  {
    for (auto& mcParticle : mcParticles) {
      // Rapidity selection
      if ((mcParticle.y() < genSelections.mc_rapidityMin) || (mcParticle.y() > genSelections.mc_rapidityMax))
        continue;

      // Selection on the source (generator/transport)
      if (genSelections.mc_keepOnlyFromGenerator && !genSelections.mc_keepOnlyFromTransport) {
        if (!mcParticle.producedByGenerator())
          continue;
      }

      if (genSelections.mc_keepOnlyFromTransport && !genSelections.mc_keepOnlyFromGenerator) {
        if (mcParticle.producedByGenerator())
          continue;
      }

      // MC Process selection
      if ((genSelections.mc_selectMCProcess >= 0) && (genSelections.mc_selectMCProcess != mcParticle.getProcess()))
        continue;

      // Get generated particle info
      auto MCGenInfo = getV0PairGenInfo(mcParticle);

      // Fill QA histos
      if (genSelections.doQA)
        fillGenQAHistos(MCGenInfo);

      // Fill tables
      // Pi0
      if (fillPi0Tables && MCGenInfo.IsPi0) {
        pi0Gens(MCGenInfo.IsProducedByGenerator, MCGenInfo.MCPt, mcParticle.y()); // optional table to store generated pi0 candidates. Be careful, this is a large table!
        pi0GenCollRefs(MCGenInfo.MCCollId);                       // link to stramccollision table
      }

      // Sigma0/ASigma0
      if (fillSigma0Tables && (MCGenInfo.IsSigma0 || MCGenInfo.IsAntiSigma0)) {
        sigma0Gens(MCGenInfo.IsSigma0, MCGenInfo.IsProducedByGenerator, MCGenInfo.MCPt, mcParticle.y());
        sigma0GenCollRefs(MCGenInfo.MCCollId); // link to stramccollision table
      }

      // KStar
      if (fillKStarTables && MCGenInfo.IsKStar) {
        kstarGens(MCGenInfo.IsKStar, MCGenInfo.IsProducedByGenerator, MCGenInfo.MCPt);
        kstarGenCollRefs(MCGenInfo.MCCollId);                       // link to stramccollision table
      }
    }
  }

  // Function to fill QA histograms. mode = 0 (before selections, all v0s), 1 (photon candidates), 2 (lambda/alambda candidates)
  template <int mode, typename TV0Object, typename TCollision>
  void fillHistos(TV0Object const& v0, TCollision const& collision)
  {
    // Check whether it is before or after selections
    static constexpr std::string_view MainDir[] = {"V0BeforeSel", "PhotonSel", "LambdaSel"};

    auto posTrack = v0.template posTrackExtra_as<dauTracks>();
    auto negTrack = v0.template negTrackExtra_as<dauTracks>();

    float Phi = RecoDecay::phi(v0.px(), v0.py());
    float PhotonY = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);
    float fLambdaLifeTime = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;

    // Fill histos
    histos.fill(HIST(MainDir[mode]) + HIST("/hpT"), v0.pt());
    histos.fill(HIST(MainDir[mode]) + HIST("/hV0Type"), v0.v0Type());
    histos.fill(HIST(MainDir[mode]) + HIST("/hNegEta"), v0.negativeeta());
    histos.fill(HIST(MainDir[mode]) + HIST("/hPosEta"), v0.positiveeta());
    histos.fill(HIST(MainDir[mode]) + HIST("/hDCANegToPV"), v0.dcanegtopv());
    histos.fill(HIST(MainDir[mode]) + HIST("/hDCAPosToPV"), v0.dcapostopv());
    histos.fill(HIST(MainDir[mode]) + HIST("/hDCADau"), v0.dcaV0daughters());
    histos.fill(HIST(MainDir[mode]) + HIST("/hRadius"), v0.v0radius());
    histos.fill(HIST(MainDir[mode]) + HIST("/hZ"), v0.z());
    histos.fill(HIST(MainDir[mode]) + HIST("/hCosPA"), v0.v0cosPA());
    histos.fill(HIST(MainDir[mode]) + HIST("/hPhi"), Phi);
    histos.fill(HIST(MainDir[mode]) + HIST("/hPosTPCCR"), posTrack.tpcCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/hNegTPCCR"), negTrack.tpcCrossedRows());
    histos.fill(HIST(MainDir[mode]) + HIST("/hPosITSNCls"), posTrack.itsNCls());
    histos.fill(HIST(MainDir[mode]) + HIST("/hNegITSNCls"), negTrack.itsNCls());
    histos.fill(HIST(MainDir[mode]) + HIST("/hPosTPCNSigmaEl"), posTrack.tpcNSigmaEl());
    histos.fill(HIST(MainDir[mode]) + HIST("/hNegTPCNSigmaEl"), negTrack.tpcNSigmaEl());
    histos.fill(HIST(MainDir[mode]) + HIST("/hPosTPCNSigmaPi"), posTrack.tpcNSigmaPi());
    histos.fill(HIST(MainDir[mode]) + HIST("/hNegTPCNSigmaPi"), negTrack.tpcNSigmaPi());
    histos.fill(HIST(MainDir[mode]) + HIST("/hPosTPCNSigmaPr"), posTrack.tpcNSigmaPr());
    histos.fill(HIST(MainDir[mode]) + HIST("/hNegTPCNSigmaPr"), negTrack.tpcNSigmaPr());

    histos.fill(HIST(MainDir[mode]) + HIST("/h2dArmenteros"), v0.alpha(), v0.qtarm());

    if (fFillV03DPositionHistos && (MainDir[mode] != "V0BeforeSel")) // don't fill for all v0s before selection
      histos.fill(HIST(MainDir[mode]) + HIST("/h3dV0XYZ"), v0.x(), v0.y(), v0.z());

    // Photon-specific
    histos.fill(HIST(MainDir[mode]) + HIST("/hPhotonY"), PhotonY);
    histos.fill(HIST(MainDir[mode]) + HIST("/hPhotonMass"), v0.mGamma());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dMassPhotonVsK0S"), v0.mGamma(), v0.mK0Short());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dMassPhotonVsLambda"), v0.mGamma(), v0.mLambda());

    // Lambda/ALambda-specific
    histos.fill(HIST(MainDir[mode]) + HIST("/hLifeTime"), fLambdaLifeTime);
    histos.fill(HIST(MainDir[mode]) + HIST("/hLambdaY"), v0.yLambda());
    histos.fill(HIST(MainDir[mode]) + HIST("/hLambdaMass"), v0.mLambda());
    histos.fill(HIST(MainDir[mode]) + HIST("/hALambdaMass"), v0.mAntiLambda());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dMassLambdaVsK0S"), v0.mLambda(), v0.mK0Short());
    histos.fill(HIST(MainDir[mode]) + HIST("/h2dMassLambdaVsGamma"), v0.mLambda(), v0.mGamma());
  }

  //_______________________________________________
  // Process photon candidate
  template <typename TV0Object>
  bool processPhotonCandidate(TV0Object const& gamma)
  {
    // Optional MC selection
    if (photonSelections.doMCAssociation) {
      if constexpr (requires { gamma.motherMCPartId(); }) {
        if (gamma.has_v0MCCore()) {
          auto gammaMC = gamma.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
          if (gammaMC.pdgCode() != 22)
            return false;
        }
      }
    }

    // V0 type selection
    if (gamma.v0Type() != photonSelections.Photonv0TypeSel && photonSelections.Photonv0TypeSel > -1)
      return false;

    float PhotonY = RecoDecay::y(std::array{gamma.px(), gamma.py(), gamma.pz()}, o2::constants::physics::MassGamma);

    if (useMLScores) {
      if (gamma.gammaBDTScore() <= photonSelections.Gamma_MLThreshold)
        return false;

    } else {
      // Standard selection
      // Gamma basic selection criteria:
      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 1.);
      if ((gamma.mGamma() < 0) || (gamma.mGamma() > photonSelections.PhotonMaxMass))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 2.);
      if ((PhotonY < photonSelections.PhotonMinRapidity) || (PhotonY > photonSelections.PhotonMaxRapidity))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 3.);
      if (gamma.negativeeta() < photonSelections.PhotonDauEtaMin || gamma.negativeeta() > photonSelections.PhotonDauEtaMax)
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 4.);
      if (gamma.positiveeta() < photonSelections.PhotonDauEtaMin || gamma.positiveeta() > photonSelections.PhotonDauEtaMax)
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 5.);
      if ((TMath::Abs(gamma.dcapostopv()) < photonSelections.PhotonMinDCADauToPv) || (TMath::Abs(gamma.dcanegtopv()) < photonSelections.PhotonMinDCADauToPv))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 6.);
      if (TMath::Abs(gamma.dcaV0daughters()) > photonSelections.PhotonMaxDCAV0Dau)
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 7.);
      if ((gamma.v0radius() < photonSelections.PhotonMinRadius) || (gamma.v0radius() > photonSelections.PhotonMaxRadius))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 8.);
      if ((gamma.z() < photonSelections.PhotonMinZ) || (gamma.z() > photonSelections.PhotonMaxZ))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 9.);
      if (gamma.v0cosPA() < photonSelections.PhotonMinV0cospa)
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 10.);
      float PhotonPhi = RecoDecay::phi(gamma.px(), gamma.py());
      if ((((PhotonPhi > photonSelections.PhotonPhiMin1) && (PhotonPhi < photonSelections.PhotonPhiMax1)) || ((PhotonPhi > photonSelections.PhotonPhiMin2) && (PhotonPhi < photonSelections.PhotonPhiMax2))) && ((photonSelections.PhotonPhiMin1 != -1) && (photonSelections.PhotonPhiMax1 != -1) && (photonSelections.PhotonPhiMin2 != -1) && (photonSelections.PhotonPhiMax2 != -1)))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 11.);
      if (gamma.qtarm() > photonSelections.PhotonMaxQt)
        return false;

      if (TMath::Abs(gamma.alpha()) > photonSelections.PhotonMaxAlpha)
        return false;

      auto posTrackGamma = gamma.template posTrackExtra_as<dauTracks>();
      auto negTrackGamma = gamma.template negTrackExtra_as<dauTracks>();

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 12.);
      if ((posTrackGamma.tpcCrossedRows() < photonSelections.PhotonMinTPCCrossedRows) || (negTrackGamma.tpcCrossedRows() < photonSelections.PhotonMinTPCCrossedRows))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 13.);
      if (((posTrackGamma.tpcNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) || (posTrackGamma.tpcNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)))
        return false;

      if (((negTrackGamma.tpcNSigmaEl() < photonSelections.PhotonMinTPCNSigmas) || (negTrackGamma.tpcNSigmaEl() > photonSelections.PhotonMaxTPCNSigmas)))
        return false;

      histos.fill(HIST("PhotonSel/hSelectionStatistics"), 14.);
    }

    return true;
  }

  //_______________________________________________
  // Process lambda candidate
  template <typename TV0Object, typename TCollision>
  bool processLambdaCandidate(TV0Object const& lambda, TCollision const& collision)
  {
    // Optional MC selection
    if (lambdaSelections.doMCAssociation) {
      if constexpr (requires { lambda.motherMCPartId(); }) {
        if (lambda.has_v0MCCore()) {
          auto lambdaMC = lambda.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
          if (TMath::Abs(lambdaMC.pdgCode()) != 3122)
            return false;
        }
      }
    }

    // V0 type selection
    if (lambda.v0Type() != 1)
      return false;

    if (useMLScores) {
      if ((lambda.lambdaBDTScore() <= lambdaSelections.Lambda_MLThreshold) && (lambda.antiLambdaBDTScore() <= lambdaSelections.AntiLambda_MLThreshold))
        return false;

    } else {
      // Lambda basic selection criteria:
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 1.);
      if ((TMath::Abs(lambda.mLambda() - o2::constants::physics::MassLambda0) > lambdaSelections.LambdaWindow) && (TMath::Abs(lambda.mAntiLambda() - o2::constants::physics::MassLambda0) > lambdaSelections.LambdaWindow) && lambdaSelections.LambdaWindow > 0)
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 2.);
      if ((lambda.yLambda() < lambdaSelections.LambdaMinRapidity) || (lambda.yLambda() > lambdaSelections.LambdaMaxRapidity))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 3.);
      if ((lambda.negativeeta() < lambdaSelections.LambdaDauEtaMin) || (lambda.negativeeta() > lambdaSelections.LambdaDauEtaMax))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 4.);
      if ((lambda.positiveeta() < lambdaSelections.LambdaDauEtaMin) || (lambda.positiveeta() > lambdaSelections.LambdaDauEtaMax))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 5.);
      if ((TMath::Abs(lambda.dcapostopv()) < lambdaSelections.LambdaMinDCAPosToPv) || (TMath::Abs(lambda.dcanegtopv()) < lambdaSelections.LambdaMinDCANegToPv))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 6.);
      if ((lambda.v0radius() < lambdaSelections.LambdaMinv0radius) || (lambda.v0radius() > lambdaSelections.LambdaMaxv0radius))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 7.);
      if ((lambda.z() < lambdaSelections.LambdaMinZ) || (lambda.z() > lambdaSelections.LambdaMaxZ))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 8.);
      if (TMath::Abs(lambda.dcaV0daughters()) > lambdaSelections.LambdaMaxDCAV0Dau)
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 9.);
      if ((lambda.qtarm() < lambdaSelections.LambdaMinQt) || (lambda.qtarm() > lambdaSelections.LambdaMaxQt))
        return false;

      if ((TMath::Abs(lambda.alpha()) < lambdaSelections.LambdaMinAlpha) || (TMath::Abs(lambda.alpha()) > lambdaSelections.LambdaMaxAlpha))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 10.);
      if (lambda.v0cosPA() < lambdaSelections.LambdaMinv0cospa)
        return false;

      auto posTrackLambda = lambda.template posTrackExtra_as<dauTracks>();
      auto negTrackLambda = lambda.template negTrackExtra_as<dauTracks>();

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 11.);
      if ((posTrackLambda.tpcCrossedRows() < lambdaSelections.LambdaMinTPCCrossedRows) || (negTrackLambda.tpcCrossedRows() < lambdaSelections.LambdaMinTPCCrossedRows))
        return false;

      // MinITSCls
      bool posIsFromAfterburner = posTrackLambda.itsChi2PerNcl() < 0;
      bool negIsFromAfterburner = negTrackLambda.itsChi2PerNcl() < 0;
      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 12.);
      if (posTrackLambda.itsNCls() < lambdaSelections.LambdaMinITSclusters && (!lambdaSelections.LambdaRejectPosITSafterburner || posIsFromAfterburner))
        return false;
      if (negTrackLambda.itsNCls() < lambdaSelections.LambdaMinITSclusters && (!lambdaSelections.LambdaRejectNegITSafterburner || negIsFromAfterburner))
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 13.);
      float fLambdaLifeTime = lambda.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      if (fLambdaLifeTime > lambdaSelections.LambdaMaxLifeTime)
        return false;

      histos.fill(HIST("LambdaSel/hSelectionStatistics"), 14.);
    }

    return true;
  }


  template <typename TV0Object, typename TCollision>
  bool processKShortCandidate(TV0Object const& kshort, TCollision const& collision)
  {

    if (kshort.v0Type() != 1)
      return false;
    
    float centrality = doPPAnalysis ? collision.centFT0M() : collision.centFT0C();

    histos.fill(HIST("K0ShortSel/h3dK0ShortMass"), centrality, kshort.pt(), kshort.mK0Short());

    histos.fill(HIST("K0ShortSel/hSelectionStatistics"), 1.);
    histos.fill(HIST("K0ShortSel/hK0ShortMass"), kshort.mK0Short());

    if ((TMath::Abs(kshort.mK0Short() - o2::constants::physics::MassK0Short) > K0ShortWindow))
      return false;
    histos.fill(HIST("K0ShortSel/hK0ShortNegEta"), kshort.negativeeta());
    histos.fill(HIST("K0ShortSel/hK0ShortPosEta"), kshort.positiveeta());
    histos.fill(HIST("K0ShortSel/hK0ShortY"), kshort.yK0Short());
    histos.fill(HIST("K0ShortSel/hSelectionStatistics"), 2.);
    
    if ((TMath::Abs(kshort.yK0Short()) > V0Rapidity) || (TMath::Abs(kshort.negativeeta()) > K0ShortDauPseudoRap) || (TMath::Abs(kshort.positiveeta()) > K0ShortDauPseudoRap))
      return false;
    histos.fill(HIST("K0ShortSel/hK0ShortDCANegToPV"), kshort.dcanegtopv());
    histos.fill(HIST("K0ShortSel/hK0ShortDCAPosToPV"), kshort.dcapostopv());
    histos.fill(HIST("K0ShortSel/hSelectionStatistics"), 3.);

    if ((TMath::Abs(kshort.dcapostopv()) < K0ShortMinDCAPosToPv) || (TMath::Abs(kshort.dcanegtopv()) < K0ShortMinDCANegToPv))
      return false;
    histos.fill(HIST("K0ShortSel/hK0ShortRadius"), kshort.v0radius());
    histos.fill(HIST("K0ShortSel/hSelectionStatistics"), 4.);

    if ((kshort.v0radius() < K0ShortMinv0radius) || (kshort.v0radius() > K0ShortMaxv0radius))
      return false;
    histos.fill(HIST("K0ShortSel/hK0ShortDCADau"), kshort.dcaV0daughters());
    histos.fill(HIST("K0ShortSel/hSelectionStatistics"), 5.);

    if (TMath::Abs(kshort.dcaV0daughters()) > K0ShortMaxDCAV0Dau)
      return false;
    histos.fill(HIST("K0ShortSel/hSelectionStatistics"), 6.);
        // MC Processing (if available)
    if constexpr (requires { kshort.motherMCPartId(); }) {
      if (kshort.has_v0MCCore()) {
        auto kshortMC = kshort.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
      //   if (kshortMC.pdgCode() == 310) // Is K0Short
      //     histos.fill(HIST("MC/h2dPtVsCentrality_MCAssocKShort"), centrality, kshort.pt());
      }
    }

    return true;

  }


  //_______________________________________________
  // Build pi0 candidate for QA
  template <typename TV0Object, typename TCollision, typename TMCParticles>
  bool buildPi0(TV0Object const& gamma1, TV0Object const& gamma2, TCollision const& collision, TMCParticles const& mcparticles)
  {
    //_______________________________________________
    // Check if both V0s are made of the same tracks
    if (gamma1.posTrackExtraId() == gamma2.posTrackExtraId() ||
        gamma1.negTrackExtraId() == gamma2.negTrackExtraId()) {
      return false;
    }

    //_______________________________________________
    // Calculate pi0 properties
    std::array<float, 3> pVecGamma1{gamma1.px(), gamma1.py(), gamma1.pz()};
    std::array<float, 3> pVecGamma2{gamma2.px(), gamma2.py(), gamma2.pz()};
    std::array arrpi0{pVecGamma1, pVecGamma2};
    float pi0Mass = RecoDecay::m(arrpi0, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassPhoton});
    float pi0Y = RecoDecay::y(std::array{gamma1.px() + gamma2.px(), gamma1.py() + gamma2.py(), gamma1.pz() + gamma2.pz()}, o2::constants::physics::MassPi0);

    //_______________________________________________
    // Pi0-specific selections:
    if (TMath::Abs(pi0Y) > Pi0MaxRap)
      return false;

    if (TMath::Abs(pi0Mass - o2::constants::physics::MassPi0) > Pi0MassWindow)
      return false;

    // Fill optional tables for QA
    // Define the table!
    auto posTrackGamma1 = gamma1.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma1 = gamma1.template negTrackExtra_as<dauTracks>();
    auto posTrackGamma2 = gamma2.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma2 = gamma2.template negTrackExtra_as<dauTracks>();

    // Calculate Pi0 topological info
    auto pi0TopoInfo = propagateV0PairToDCA(gamma1, gamma2);

    // Check if MC data and populate corresponding table
    if constexpr (requires { gamma1.motherMCPartId(); gamma2.motherMCPartId(); }) {
      auto pi0MCInfo = getV0PairMCInfo(gamma1, gamma2, collision, mcparticles);

      pi0coresmc(pi0MCInfo.V0PairMCRadius, pi0MCInfo.V0PairPDGCode, pi0MCInfo.V0PairPDGCodeMother, pi0MCInfo.V0PairMCProcess, pi0MCInfo.fV0PairProducedByGenerator,
                 pi0MCInfo.V01MCpx, pi0MCInfo.V01MCpy, pi0MCInfo.V01MCpz,
                 pi0MCInfo.fIsV01Primary, pi0MCInfo.V01PDGCode, pi0MCInfo.V01PDGCodeMother, pi0MCInfo.fIsV01CorrectlyAssign,
                 pi0MCInfo.V02MCpx, pi0MCInfo.V02MCpy, pi0MCInfo.V02MCpz,
                 pi0MCInfo.fIsV02Primary, pi0MCInfo.V02PDGCode, pi0MCInfo.V02PDGCodeMother, pi0MCInfo.fIsV02CorrectlyAssign);
    }

    pi0cores(pi0TopoInfo.X, pi0TopoInfo.Y, pi0TopoInfo.Z, pi0TopoInfo.DCADau, pi0TopoInfo.CosPA,
             gamma1.px(), gamma1.py(), gamma1.pz(),
             gamma1.mGamma(), gamma1.qtarm(), gamma1.alpha(), gamma1.dcapostopv(), gamma1.dcanegtopv(), gamma1.dcaV0daughters(),
             gamma1.negativeeta(), gamma1.positiveeta(), gamma1.v0cosPA(), gamma1.v0radius(), gamma1.z(),
             posTrackGamma1.tpcCrossedRows(), negTrackGamma1.tpcCrossedRows(), posTrackGamma1.tpcNSigmaEl(), negTrackGamma1.tpcNSigmaEl(), gamma1.v0Type(),
             gamma2.px(), gamma2.py(), gamma2.pz(),
             gamma2.mGamma(), gamma2.qtarm(), gamma2.alpha(), gamma2.dcapostopv(), gamma2.dcanegtopv(), gamma2.dcaV0daughters(),
             gamma2.negativeeta(), gamma2.positiveeta(), gamma2.v0cosPA(), gamma2.v0radius(), gamma2.z(),
             posTrackGamma2.tpcCrossedRows(), negTrackGamma2.tpcCrossedRows(), posTrackGamma2.tpcNSigmaEl(), negTrackGamma2.tpcNSigmaEl(), gamma2.v0Type());

    pi0coresRefs(collision.globalIndex());

    return true;
  }

  //_______________________________________________
  // Build sigma0 candidate
  template <typename TV0Object, typename TCollision, typename TMCParticles>
  bool buildSigma0(TV0Object const& lambda, TV0Object const& gamma, TCollision const& collision, TMCParticles const& mcparticles)
  {
    //_______________________________________________
    // Checking if both V0s are made of the very same tracks
    if (gamma.posTrackExtraId() == lambda.posTrackExtraId() ||
        gamma.negTrackExtraId() == lambda.negTrackExtraId()) {
      return false;
    }

    //_______________________________________________
    // Sigma0 pre-selections
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecLambda{lambda.px(), lambda.py(), lambda.pz()};

    auto arrMom = std::array{pVecPhotons, pVecLambda};
    float sigmaMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassLambda0});
    float sigmaY = -999.f;

    if constexpr (requires { gamma.pxMC(); lambda.pxMC(); }) // If MC
      sigmaY = RecoDecay::y(std::array{gamma.pxMC() + lambda.pxMC(), gamma.pyMC() + lambda.pyMC(), gamma.pzMC() + lambda.pzMC()}, o2::constants::physics::MassSigma0);
    else // If DATA
      sigmaY = RecoDecay::y(std::array{gamma.px() + lambda.px(), gamma.py() + lambda.py(), gamma.pz() + lambda.pz()}, o2::constants::physics::MassSigma0);

    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 1.);
    if (TMath::Abs(sigmaMass - o2::constants::physics::MassSigma0) > Sigma0Window)
      return false;

    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 2.);
    if (TMath::Abs(sigmaY) > SigmaMaxRap)
      return false;

    histos.fill(HIST("SigmaSel/hSigmaMassSelected"), sigmaMass);
    histos.fill(HIST("SigmaSel/hSelectionStatistics"), 3.);
    //_______________________________________________
    // Calculate properties & Fill tables

    // Sigma0 topological info
    auto sigma0TopoInfo = propagateV0PairToDCA(gamma, lambda);

    sigma0cores(sigma0TopoInfo.X, sigma0TopoInfo.Y, sigma0TopoInfo.Z, sigma0TopoInfo.DCADau,
                gamma.px(), gamma.py(), gamma.pz(), gamma.mGamma(), lambda.px(), lambda.py(), lambda.pz(), lambda.mLambda(), lambda.mAntiLambda());

    // MC properties
    if constexpr (requires { gamma.motherMCPartId(); lambda.motherMCPartId(); }) {
      auto sigma0MCInfo = getV0PairMCInfo(gamma, lambda, collision, mcparticles);

      sigma0mccores(sigma0MCInfo.V0PairMCRadius, sigma0MCInfo.V0PairPDGCode, sigma0MCInfo.V0PairPDGCodeMother, sigma0MCInfo.V0PairMCProcess, sigma0MCInfo.fV0PairProducedByGenerator,
                    sigma0MCInfo.V01MCpx, sigma0MCInfo.V01MCpy, sigma0MCInfo.V01MCpz,
                    sigma0MCInfo.fIsV01Primary, sigma0MCInfo.V01PDGCode, sigma0MCInfo.V01PDGCodeMother, sigma0MCInfo.fIsV01CorrectlyAssign,
                    sigma0MCInfo.V02MCpx, sigma0MCInfo.V02MCpy, sigma0MCInfo.V02MCpz,
                    sigma0MCInfo.fIsV02Primary, sigma0MCInfo.V02PDGCode, sigma0MCInfo.V02PDGCodeMother, sigma0MCInfo.fIsV02CorrectlyAssign);

      sigma0mclabel(sigma0MCInfo.V0PairMCParticleID);
    }

    // Sigma0s -> stracollisions link
    histos.fill(HIST("SigmaSel/hSigma0DauDeltaIndex"), gamma.globalIndex() - lambda.globalIndex());
    sigma0CollRefs(collision.globalIndex());
    sigmaIndices(gamma.globalIndex(), lambda.globalIndex());

    //_______________________________________________
    // Photon extra properties
    auto posTrackGamma = gamma.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma = gamma.template negTrackExtra_as<dauTracks>();

    uint8_t fPhotonPosTrackCode = ((uint8_t(posTrackGamma.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrackGamma.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrackGamma.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrackGamma.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrackGamma.hasTOF()) << hasTOF));

    uint8_t fPhotonNegTrackCode = ((uint8_t(negTrackGamma.hasTPC()) << hasTPC) |
                                   (uint8_t(negTrackGamma.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(negTrackGamma.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(negTrackGamma.hasTRD()) << hasTRD) |
                                   (uint8_t(negTrackGamma.hasTOF()) << hasTOF));

    sigmaPhotonExtras(gamma.qtarm(), gamma.alpha(), gamma.v0cosPA(), gamma.dcaV0daughters(), gamma.dcanegtopv(), gamma.dcapostopv(), gamma.v0radius(), gamma.z(),
                      posTrackGamma.tpcNSigmaEl(), negTrackGamma.tpcNSigmaEl(), posTrackGamma.tpcCrossedRows(), negTrackGamma.tpcCrossedRows(),
                      gamma.positiveeta(), gamma.negativeeta(), gamma.psipair(), posTrackGamma.itsNCls(), negTrackGamma.itsNCls(), posTrackGamma.itsChi2PerNcl(), negTrackGamma.itsChi2PerNcl(),
                      fPhotonPosTrackCode, fPhotonNegTrackCode, gamma.v0Type());

    //_______________________________________________
    // Lambda extra properties
    auto posTrackLambda = lambda.template posTrackExtra_as<dauTracks>();
    auto negTrackLambda = lambda.template negTrackExtra_as<dauTracks>();

    uint8_t fLambdaPosTrackCode = ((uint8_t(posTrackLambda.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrackLambda.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrackLambda.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrackLambda.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrackLambda.hasTOF()) << hasTOF));

    uint8_t fLambdaNegTrackCode = ((uint8_t(negTrackLambda.hasTPC()) << hasTPC) |
                                   (uint8_t(negTrackLambda.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(negTrackLambda.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(negTrackLambda.hasTRD()) << hasTRD) |
                                   (uint8_t(negTrackLambda.hasTOF()) << hasTOF));

    float fLambdaLifeTime = lambda.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    float fLambdaPrTOFNSigma = -999.f;
    float fLambdaPiTOFNSigma = -999.f;
    float fALambdaPrTOFNSigma = -999.f;
    float fALambdaPiTOFNSigma = -999.f;

    if constexpr (requires { lambda.tofNSigmaLaPr(); }) { // If TOF info avaiable
      fLambdaPrTOFNSigma = lambda.tofNSigmaLaPr();
      fLambdaPiTOFNSigma = lambda.tofNSigmaLaPi();
      fALambdaPrTOFNSigma = lambda.tofNSigmaALaPr();
      fALambdaPiTOFNSigma = lambda.tofNSigmaALaPi();
    }

    sigmaLambdaExtras(lambda.qtarm(), lambda.alpha(), fLambdaLifeTime, lambda.v0radius(), lambda.v0cosPA(), lambda.dcaV0daughters(), lambda.dcanegtopv(), lambda.dcapostopv(),
                      posTrackLambda.tpcNSigmaPr(), posTrackLambda.tpcNSigmaPi(), negTrackLambda.tpcNSigmaPr(), negTrackLambda.tpcNSigmaPi(),
                      fLambdaPrTOFNSigma, fLambdaPiTOFNSigma, fALambdaPrTOFNSigma, fALambdaPiTOFNSigma,
                      posTrackLambda.tpcCrossedRows(), negTrackLambda.tpcCrossedRows(),
                      lambda.positiveeta(), lambda.negativeeta(),
                      posTrackLambda.itsNCls(), negTrackLambda.itsNCls(), posTrackLambda.itsChi2PerNcl(), negTrackLambda.itsChi2PerNcl(),
                      fLambdaPosTrackCode, fLambdaNegTrackCode, lambda.v0Type());

    return true;
  }


  //_______________________________________________
  // Build kstar candidate
  template <typename TV0Object, typename TCollision, typename TMCParticles> 
  bool buildKStar(TV0Object const& kshort, TV0Object const& gamma, TCollision const & collision, TMCParticles const& mcparticles)
  {

    //_______________________________________________
    // Checking if both V0s are made of the very same tracks
    if (gamma.posTrackExtraId() == kshort.posTrackExtraId() ||
        gamma.negTrackExtraId() == kshort.negTrackExtraId()) {
      return false;
    }

    //_______________________________________________
    // Kstar pre-selections
    std::array<float, 3> pVecPhotons{gamma.px(), gamma.py(), gamma.pz()};
    std::array<float, 3> pVecKShort{kshort.px(), kshort.py(), kshort.pz()};

    auto arrMom = std::array{pVecPhotons, pVecKShort};
    float kstarMass = RecoDecay::m(arrMom, std::array{o2::constants::physics::MassPhoton, o2::constants::physics::MassK0Short});
    float kstarY = -999.f;

    if constexpr (requires { gamma.pxMC(); kshort.pxMC(); }) // If MC
      kstarY = RecoDecay::y(std::array{gamma.pxMC() + kshort.pxMC(), gamma.pyMC() + kshort.pyMC(), gamma.pzMC() + kshort.pzMC()}, o2::constants::physics::MassK0Star892);
    else // If DATA
      kstarY = RecoDecay::y(std::array{gamma.px() + kshort.px(), gamma.py() + kshort.py(), gamma.pz() + kshort.pz()}, o2::constants::physics::MassK0Star892);

    histos.fill(HIST("KStarSel/hSelectionStatistics"), 1.);
    if (TMath::Abs(kstarMass - o2::constants::physics::MassK0Star892) > KStarWindow)
      return false;
    
    histos.fill(HIST("KStarSel/hSelectionStatistics"), 2.);
    if (TMath::Abs(kstarY) > KStarMaxRap)
      return false;

    histos.fill(HIST("KStarSel/hKStarMassSelected"), kstarMass);
    histos.fill(HIST("KStarSel/hSelectionStatistics"), 3.);

    auto kstarTopoInfo = propagateV0PairToDCA(gamma, kshort);

    kstarcores(kstarTopoInfo.X, kstarTopoInfo.Y, kstarTopoInfo.Z, kstarTopoInfo.DCADau,
                gamma.px(), gamma.py(), gamma.pz(), gamma.mGamma(), kshort.px(), kshort.py(), kshort.pz(), kshort.mK0Short());

    
    // MC properties
    if constexpr (requires { gamma.motherMCPartId(); kshort.motherMCPartId(); }) {
      auto kstarMCInfo = getV0PairMCInfo(gamma, kshort, collision, mcparticles);

      kstarmccores(kstarMCInfo.V0PairMCRadius, kstarMCInfo.V0PairPDGCode, kstarMCInfo.V0PairPDGCodeMother, kstarMCInfo.V0PairMCProcess, kstarMCInfo.fV0PairProducedByGenerator,
                    kstarMCInfo.V01MCpx, kstarMCInfo.V01MCpy, kstarMCInfo.V01MCpz,
                    kstarMCInfo.fIsV01Primary, kstarMCInfo.V01PDGCode, kstarMCInfo.V01PDGCodeMother, kstarMCInfo.fIsV01CorrectlyAssign,
                    kstarMCInfo.V02MCpx, kstarMCInfo.V02MCpy, kstarMCInfo.V02MCpz,
                    kstarMCInfo.fIsV02Primary, kstarMCInfo.V02PDGCode, kstarMCInfo.V02PDGCodeMother, kstarMCInfo.fIsV02CorrectlyAssign);
    }

    // KStar -> stracollisions link
    kstarCollRefs(collision.globalIndex());

    //_______________________________________________
    // Photon extra properties

    auto posTrackGamma = gamma.template posTrackExtra_as<dauTracks>();
    auto negTrackGamma = gamma.template negTrackExtra_as<dauTracks>();

    uint8_t fPhotonPosTrackCode = ((uint8_t(posTrackGamma.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrackGamma.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrackGamma.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrackGamma.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrackGamma.hasTOF()) << hasTOF));

    uint8_t fPhotonNegTrackCode = ((uint8_t(negTrackGamma.hasTPC()) << hasTPC) |
                                    (uint8_t(negTrackGamma.hasITSTracker()) << hasITSTracker) |
                                    (uint8_t(negTrackGamma.hasITSAfterburner()) << hasITSAfterburner) |
                                    (uint8_t(negTrackGamma.hasTRD()) << hasTRD) |
                                    (uint8_t(negTrackGamma.hasTOF()) << hasTOF));

    kstarPhotonExtras(gamma.qtarm(), gamma.alpha(), gamma.v0cosPA(), gamma.dcaV0daughters(), gamma.dcanegtopv(), gamma.dcapostopv(), gamma.v0radius(), gamma.z(),
                      posTrackGamma.tpcNSigmaEl(), negTrackGamma.tpcNSigmaEl(), posTrackGamma.tpcCrossedRows(), negTrackGamma.tpcCrossedRows(),
                      gamma.positiveeta(), gamma.negativeeta(), gamma.psipair(), posTrackGamma.itsNCls(), negTrackGamma.itsNCls(), posTrackGamma.itsChi2PerNcl(), negTrackGamma.itsChi2PerNcl(),
                      fPhotonPosTrackCode, fPhotonNegTrackCode, gamma.v0Type());


    //_______________________________________________
    // KShort extra properties

    auto posTrackKShort = kshort.template posTrackExtra_as<dauTracks>();
    auto negTrackKShort = kshort.template negTrackExtra_as<dauTracks>();

    uint8_t fKShortPosTrackCode = ((uint8_t(posTrackKShort.hasTPC()) << hasTPC) |
                                   (uint8_t(posTrackKShort.hasITSTracker()) << hasITSTracker) |
                                   (uint8_t(posTrackKShort.hasITSAfterburner()) << hasITSAfterburner) |
                                   (uint8_t(posTrackKShort.hasTRD()) << hasTRD) |
                                   (uint8_t(posTrackKShort.hasTOF()) << hasTOF));

    uint8_t fKShortNegTrackCode = ((uint8_t(negTrackKShort.hasTPC()) << hasTPC) |
                                    (uint8_t(negTrackKShort.hasITSTracker()) << hasITSTracker) |
                                    (uint8_t(negTrackKShort.hasITSAfterburner()) << hasITSAfterburner) |
                                    (uint8_t(negTrackKShort.hasTRD()) << hasTRD) |
                                    (uint8_t(negTrackKShort.hasTOF()) << hasTOF));

    float fKShortLifeTime = kshort.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    float fKShortPiTOFNSigma = -999.f;

    // float postrackTofNSigmaPi = -999.f;
    // float negtrackTofNSigmaPi = -999.f;

    if constexpr (requires { kshort.tofNSigmaK0Pi(); }) { // If TOF info avaiable

      // postrackTofNSigmaPi = posTrackKShort.tofNSigmaPi();
      // negtrackTofNSigmaPi = negTrackKShort.tofNSigmaPi();
      fKShortPiTOFNSigma = kshort.tofNSigmaLaPi();

    }
 
    kshortExtras(kshort.qtarm(), kshort.alpha(), fKShortLifeTime, kshort.v0radius(), kshort.v0cosPA(), kshort.dcaV0daughters(), kshort.dcanegtopv(), kshort.dcapostopv(),
                      posTrackKShort.tpcNSigmaPi(), negTrackKShort.tpcNSigmaPi(), fKShortPiTOFNSigma,
                      posTrackKShort.tpcCrossedRows(), negTrackKShort.tpcCrossedRows(),
                      kshort.positiveeta(), kshort.negativeeta(),
                      posTrackKShort.itsNCls(), negTrackKShort.itsNCls(), posTrackKShort.itsChi2PerNcl(), negTrackKShort.itsChi2PerNcl(),
                      fKShortPosTrackCode, fKShortNegTrackCode, kshort.v0Type());

    return true;
  }




  // Process photon and lambda candidates to build sigma0 candidates
  template <typename TCollision, typename TV0s, typename TMCParticles>
  void dataProcess(TCollision const& collisions, TV0s const& fullV0s, TMCParticles const& mcparticles)
  {
    //_______________________________________________
    // Initial setup
    // Auxiliary vectors to store best candidates
    std::vector<int> bestGammasArray;
    std::vector<int> bestLambdasArray;

    std::vector<int> bestKStarGammasArray;
    std::vector<int> bestKShortsArray;

    // Custom grouping
    std::vector<std::vector<int>> v0grouped(collisions.size());

    for (const auto& v0 : fullV0s) {
      v0grouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }

    //_______________________________________________
    // Collisions loop
    for (const auto& coll : collisions) {
      // Event selection
      if (eventSelections.fUseEventSelection) {
        if (!IsEventAccepted(coll, true))
          continue;
      }

      // Clear vectors
      bestGammasArray.clear();
      bestLambdasArray.clear();

      bestKStarGammasArray.clear();
      bestKShortsArray.clear();

      float centrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();
      histos.fill(HIST("hEventCentrality"), centrality);

      //_______________________________________________
      // V0s loop
      for (size_t i = 0; i < v0grouped[coll.globalIndex()].size(); i++) {
        auto v0 = fullV0s.rawIteratorAt(v0grouped[coll.globalIndex()][i]);

        if (fFillNoSelV0Histos)
          fillHistos<0>(v0, coll); // Filling "all V0s" histograms

        if (processPhotonCandidate(v0)) { // selecting photons
          if (fFillSelPhotonHistos)
            fillHistos<1>(v0, coll);                   // QA histos
          bestGammasArray.push_back(v0.globalIndex()); // Save indices of best gamma candidates
        }

        if (processLambdaCandidate(v0, coll)) { // selecting lambdas
          if (fFillSelLambdaHistos)
            fillHistos<2>(v0, coll);                    // QA histos
          bestLambdasArray.push_back(v0.globalIndex()); // Save indices of best lambda candidates
        }

        if (processKShortCandidate(v0, coll))           // selecting kshorts
          bestKShortsArray.push_back(v0.globalIndex()); // Save indices of best kshort candidates

      }

      //_______________________________________________
      // Wrongly collision association study (MC-specific)
      if constexpr (requires { coll.StraMCCollisionId(); }) {
        if (doAssocStudy) {
          analyzeV0CollAssoc(coll, fullV0s, bestGammasArray, true);   // Photon-analysis
          analyzeV0CollAssoc(coll, fullV0s, bestLambdasArray, false); // Lambda-analysis
        }
      }

      //_______________________________________________
      // V0 nested loop
      for (size_t i = 0; i < bestGammasArray.size(); ++i) {
        auto gamma1 = fullV0s.rawIteratorAt(bestGammasArray[i]);

        //_______________________________________________
        // Sigma0 loop
        if (fillSigma0Tables) {
          for (size_t j = 0; j < bestLambdasArray.size(); ++j) {
            auto lambda = fullV0s.rawIteratorAt(bestLambdasArray[j]);

            // Building sigma0 candidate & filling tables
            if (!buildSigma0(lambda, gamma1, coll, mcparticles))
              continue;
          }
        }

        //_______________________________________________
        // KStar loop 
        if (fillKStarTables) {
          for (size_t j = 0; j < bestKShortsArray.size(); ++j) {
            auto kshort = fullV0s.rawIteratorAt(bestKShortsArray[j]);

            // Building kstar candidate & filling tables
            if (!buildKStar(kshort, gamma1, coll, mcparticles))
              continue;
          }
        }


        //_______________________________________________
        // pi0 loop
        if (fillPi0Tables) {
          for (size_t j = i + 1; j < bestGammasArray.size(); ++j) {
            auto gamma2 = fullV0s.rawIteratorAt(bestGammasArray[j]);

            // Building pi0 candidate & filling tables
            if (!buildPi0(gamma1, gamma2, coll, mcparticles))
              continue;
          }
        }
      }
    }
  }

  // Process photon and lambda candidates
  template <typename TCollision, typename TV0s>
  void runPhotonLambdaQA(TCollision const& collisions, TV0s const& fullV0s)
  {
    // Custom grouping
    std::vector<std::vector<int>> v0grouped(collisions.size());

    for (const auto& v0 : fullV0s) {
      v0grouped[v0.straCollisionId()].push_back(v0.globalIndex());
    }

    //_______________________________________________
    // Collisions loop
    for (const auto& coll : collisions) {
      // Event selection
      if (eventSelections.fUseEventSelection) {
        if (!IsEventAccepted(coll, true))
          continue;
      }

      float centrality = doPPAnalysis ? coll.centFT0M() : coll.centFT0C();
      histos.fill(HIST("PhotonLambdaQA/hEventCentrality"), centrality);

      //_______________________________________________
      // V0s loop
      for (size_t i = 0; i < v0grouped[coll.globalIndex()].size(); i++) {
        bool fPassPhotonSel = false;
        bool fPassLambdaSel = false;

        // Get V0 object
        auto v0 = fullV0s.rawIteratorAt(v0grouped[coll.globalIndex()][i]);

        if (fFillNoSelV0Histos)
          fillHistos<0>(v0, coll); // Filling "all V0s" histograms

        // Selection process
        fPassPhotonSel = processPhotonCandidate(v0);
        fPassLambdaSel = processLambdaCandidate(v0, coll);

        // Reco part:
        if (fPassPhotonSel) {
          if (fFillSelPhotonHistos)
            fillHistos<1>(v0, coll);

          // Fill analysis Histos:
          float PhotonY = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);
          histos.fill(HIST("PhotonLambdaQA/h3dPhotonMass"), centrality, v0.pt(), v0.mGamma());
          histos.fill(HIST("PhotonLambdaQA/h3dYPhotonMass"), PhotonY, v0.pt(), v0.mGamma());
          histos.fill(HIST("PhotonLambdaQA/h3dYPhotonRadius"), PhotonY, v0.pt(), v0.v0radius());
        }
        if (fPassLambdaSel) {
          if (fFillSelLambdaHistos)
            fillHistos<2>(v0, coll);

          // Fill analysis Histos:
          histos.fill(HIST("PhotonLambdaQA/h3dLambdaMass"), centrality, v0.pt(), v0.mLambda());
          histos.fill(HIST("PhotonLambdaQA/h3dALambdaMass"), centrality, v0.pt(), v0.mAntiLambda());

          histos.fill(HIST("PhotonLambdaQA/h3dYLambdaMass"), v0.yLambda(), v0.pt(), v0.mLambda());
          histos.fill(HIST("PhotonLambdaQA/h3dYALambdaMass"), v0.yLambda(), v0.pt(), v0.mAntiLambda());

          histos.fill(HIST("PhotonLambdaQA/h3dYRLambdaMass"), v0.yLambda(), v0.v0radius(), v0.mLambda());
          histos.fill(HIST("PhotonLambdaQA/h3dYRALambdaMass"), v0.yLambda(), v0.v0radius(), v0.mAntiLambda());
        }

        // MC part:
        if constexpr (requires { v0.motherMCPartId(); }) {
          if (v0.has_v0MCCore()) {
            auto v0MC = v0.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

            // True Photon
            if (v0MC.pdgCode() == 22 && fPassPhotonSel) {
              histos.fill(HIST("PhotonLambdaQA/h3dTruePhotonMass"), centrality, v0MC.ptMC(), v0.mGamma());
              if (TMath::Abs(v0MC.pdgCodeMother()) == 3212) { // If from sigma0 or ASigma0
                histos.fill(HIST("PhotonLambdaQA/h2dTrueSigma0PhotonMass"), v0MC.ptMC(), v0.mGamma());
              }
            }

            // True Lambda
            if (v0MC.pdgCode() == 3122 && fPassLambdaSel) {
              histos.fill(HIST("PhotonLambdaQA/h3dTrueLambdaMass"), centrality, v0MC.ptMC(), v0.mLambda());
              if (v0MC.pdgCodeMother() == 3212) { // If from sigma0
                histos.fill(HIST("PhotonLambdaQA/h2dTrueSigma0LambdaMass"), v0MC.ptMC(), v0.mLambda());
              }
            }

            // True ALambda
            if (v0MC.pdgCode() == -3122 && fPassLambdaSel) {
              histos.fill(HIST("PhotonLambdaQA/h3dTrueALambdaMass"), centrality, v0MC.ptMC(), v0.mAntiLambda());
              if (v0MC.pdgCodeMother() == -3212) { // If from asigma0
                histos.fill(HIST("PhotonLambdaQA/h2dTrueASigma0ALambdaMass"), v0MC.ptMC(), v0.mAntiLambda());
              }
            }
          }
        }
      }
    }
  }

  // Sigma0 processing part
  void processRealData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0StandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    dataProcess(collisions, fullV0s, nullptr);
  }

  void processRealDataWithTOF(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0TOFStandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    dataProcess(collisions, fullV0s, nullptr);
  }

  void processMonteCarlo(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0DerivedMCDatas const& fullV0s, aod::McParticles const& mcParticles, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    dataProcess(collisions, fullV0s, mcParticles);
  }

  void processMonteCarloWithTOF(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0TOFDerivedMCDatas const& fullV0s, aod::McParticles const& mcParticles, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    dataProcess(collisions, fullV0s, mcParticles);
  }

  void processGeneratedRun3(aod::McParticles const& mcParticles)
  {
    genProcess(mcParticles);
  }

  // Photon and lambda-specific part (QA)
  void processPhotonLambdaQA(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps> const& collisions, V0StandardDerivedDatas const& fullV0s, dauTracks const&)
  {
    runPhotonLambdaQA(collisions, fullV0s);
  }

  void processPhotonLambdaMCQA(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions, V0DerivedMCDatas const& fullV0s, dauTracks const&, aod::MotherMCParts const&, soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    runPhotonLambdaQA(collisions, fullV0s);
  }

  void processPhotonLambdaGenerated(soa::Join<aod::StraMCCollisions, aod::StraMCCollMults> const& mcCollisions, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const& V0MCCores, soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraCollLabels> const& collisions)
  {
    runGenPhotonLambdaQA(mcCollisions, V0MCCores, collisions);
  }

  PROCESS_SWITCH(sigma0builder, processRealData, "process as if real data", true);
  PROCESS_SWITCH(sigma0builder, processRealDataWithTOF, "process as if real data", false);
  PROCESS_SWITCH(sigma0builder, processMonteCarlo, "process as if MC data", false);
  PROCESS_SWITCH(sigma0builder, processMonteCarloWithTOF, "process as if MC data, uses TOF PID info", false);
  PROCESS_SWITCH(sigma0builder, processGeneratedRun3, "process generated MC info", false);
  PROCESS_SWITCH(sigma0builder, processPhotonLambdaQA, "process QA of lambdas and photons", false);
  PROCESS_SWITCH(sigma0builder, processPhotonLambdaMCQA, "process QA of lambdas and photons", false);
  PROCESS_SWITCH(sigma0builder, processPhotonLambdaGenerated, "process QA of gen lambdas and photons", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<sigma0builder>(cfgc)};
}
