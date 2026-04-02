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

/// \file mftReassociationValidation.cxx
/// \brief validation task for MFT DCAxy and DCAxyz reassociation
/// \author Alexian Lejeune <alexian.lejeune@cern.ch >, Czech Technical University in Prague

#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <THn.h>
#include <TString.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// add a temporary change

using namespace o2;
using namespace o2::aod::rctsel;
using namespace o2::aod::track;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum Reassociation2dMftTracks {
  AllAmbiguousTracksAfterTrackSelectionsFor2d = 0,
  NotReassociated2dMftTracks,
  Reassociated2dMftTracks,
  NReassociation2dMftTracksSteps
};

enum Reassociation3dMftTracks {
  AllAmbiguousTracksAfterTrackSelectionsFor3d = 0,
  NotReassociated3dMftTracks,
  Reassociated3dMftTracks,
  NReassociation3dMftTracksSteps
};

enum MatchedToTrueCollisionStep {
  AllTracks = 0,
  IsNotMatchedToTrueCollision,
  IsMatchedToTrueCollision,
  NMatchedToTrueCollisionSteps
};

enum DataType {
  Data,
  Mc
};

enum SpecificEventSelectionStep {
  AllEventsPrecise = 0,
  HasMcCollision,
  IsNotSplitVertex,
  IsSel8,
  IsNoSameBunchPileup,
  IsGoodItsLayersAll,
  IsGoodZvtxFT0vsPV,
  IsNoCollInRofStandard,
  IsNoCollInRofStrict,
  IsNoCollInTimeRangeStandard,
  IsNoCollInTimeRangeStrict,
  IsNoHighMultCollInPrevRof,
  IsRctFlagChecked,
  IsWithinZvtxWindow,
  NSpecificEventSelectionSteps
};

enum SpecificTrackSelectionStep {
  AllTracksPrecise = 0,
  IsTrueTrack,
  AfterOrphanCut,
  AfterEtaCut,
  AfterClusterCut,
  AfterPtCut,
  AfterDcaXYCut,
  AfterDcaZCut,
  IsCATrack,
  IsLTFTrack,
  HasMcParticle,
  NSpecificTrackSelectionSteps
};

enum TrackAmbiguityCheckStep {
  AllTracksCheck = 0,
  IsAmbDegreeEqualToCompatibleCollIdsSize,
  NTrackAmbiguityCheckSteps
};

enum MonteCarloEventSelectionStep {
  AllMonteCarloEvents = 0,
  MonteCarloEventsAfterEventSelection,
  HasMonteCarloCollision,
  HasNotMonteCarloCollision,
  NMonteCarloEventSelectionSteps
};

enum MonteCarloTrackSelectionStep {
  AllMonteCarloTracks = 0,
  MonteCarloTracksAfterTrackSelection,
  HasMonteCarloParticle,
  HasNotMonteCarloParticle,
  NMonteCarloTrackSelectionSteps
};

enum MftTrackAmbiguityStep {
  AllMftTracks = 0,
  AfterTrackSelection,
  NumberOfAmbiguousTracks,
  NumberOfNonAmbiguousTracks,
  NMftAmbiguitySteps
};

enum MftAmbiguousAndMatchedToTrueCollisionStep {
  IsAmbiguous = 0,
  IsAmbiguousAndMatchedToTrueCollision,
  IsAmbiguousAndNotMatchedToTrueCollision,
  NMftAmbiguousAndMatchedToTrueCollisionSteps
};

enum MftNonAmbiguousAndMatchedToTrueCollisionStep {
  IsNonAmbiguous = 0,
  IsNonAmbiguousAndMatchedToTrueCollision,
  IsNonAmbiguousAndNotMatchedToTrueCollision,
  NMftNonAmbiguousAndMatchedToTrueCollisionSteps
};

enum Mft2dReassociatedAndMatchedToTrueCollisionStep {
  Is2dReassociated = 0,
  Is2dReassociatedAndMatchedToTrueCollision,
  Is2dReassociatedAndNotMatchedToTrueCollision,
  NMft2dReassociatedAndMatchedToTrueCollisionSteps
};

enum Mft3dReassociatedAndMatchedToTrueCollisionStep {
  Is3dReassociated = 0,
  Is3dReassociatedAndMatchedToTrueCollision,
  Is3dReassociatedAndNotMatchedToTrueCollision,
  NMft3dReassociatedAndMatchedToTrueCollisionSteps
};

enum MftNot2dReassociatedAndMatchedToTrueCollisionStep {
  IsNot2dReassociated = 0,
  IsNot2dReassociatedAndMatchedToTrueCollision,
  IsNot2dReassociatedAndNotMatchedToTrueCollision,
  NMftNot2dReassociatedAndMatchedToTrueCollisionSteps
};

enum MftNot3dReassociatedAndMatchedToTrueCollisionStep {
  IsNot3dReassociated = 0,
  IsNot3dReassociatedAndMatchedToTrueCollision,
  IsNot3dReassociatedAndNotMatchedToTrueCollision,
  NMftNot3dReassociatedAndMatchedToTrueCollisionSteps
};

enum MftIsTrueCollisionAmongCompatibleCollisionsStep {
  AllWronglyAssociatedTracks = 0,
  IsTrueCollisionAmongCompatibleCollisions,
  IsNotTrueCollisionAmongCompatibleCollisions,
  NMftIsTrueCollisionAmongCompatibleCollisionsSteps
};

enum MftTrackSelectionStep {
  NoSelection = 0,
  Eta,
  Cluster,
  Pt,
  IsCA,
  IsLTF,
  NMftTrackSelectionSteps
};

enum MultiplicityEstimators {
  MultNTracksPV = 0,
  MultNumContrib,
  MultFT0C,
  MultFT0M
};

static constexpr std::string_view WhatDataType[] = {"Data/", "MC/"};
static constexpr std::string_view WhatMultiplicityEstimator[] = {"multNTracksPV", "multNumContrib", "multFT0C", "multFT0M"};
std::unordered_map<int, float> recoVtxX;
std::unordered_map<int, float> recoVtxY;
std::unordered_map<int, float> recoVtxZ;
std::unordered_map<int, float> recoMcCollisionId;
std::unordered_map<int, float> recoMcCollBestCollisionIndex;

struct MftReassociationValidation {

  struct : ConfigurableGroup {
    std::string prefix = "ConfigCcdb_group";
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } configCcdb;

  //  configurables for processing options

  struct : ConfigurableGroup {
    std::string prefix = "ConfigTask_group";
    Configurable<bool> centralityBinsForMc{"centralityBinsForMc", false, "falsce = OFF, true = ON for data like multiplicity/centrality bins for MC steps"};
    Configurable<bool> keepOnlyPhysicalPrimary{"keepOnlyPhysicalPrimary", true, "Keep only physical primary particles"};
    Configurable<bool> keepOnlySecondaries{"keepOnlySecondaries", false, "Keep only secondary particles"};
    Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", true, "flag to apply z shift from CCDB"};
    Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
    Configurable<float> cfgManualZShift{"cfgManualZShift", 0.0f, "manual z-shift for propagation of global muon to PV"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } configTask;

  //   configurables for collisions
  struct : ConfigurableGroup {
    std::string prefix = "ConfigCollision_group";
    Configurable<bool> isApplyGoodItsLayersAll{"isApplyGoodItsLayersAll", false, "Enable GoodITSLayersAll"};
    Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", false, "Enable GoodZvtxFT0vsPV cut"};
    Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", false, "Enable SameBunchPileup cut"};
    Configurable<int> maxMultiplicity{"maxMultiplicity", 300, "maximum multiplicity selection for collision"};
    Configurable<int> minMultiplicity{"minMultiplicity", 0, "minimum multiplicity selection for collision"};
    Configurable<int> multiplicityEstimator{"multiplicityEstimator", 0, "0: multNTracksPV, 1: numContrib, 2: multFT0C, 3: multFT0M, 4: centFT0C, 5: centFT0CVariants1s, 6: centFT0M, 7: centFV0A, 8: centNTracksPV, 9: centNGlobal, 10: centMFT"};
    Configurable<bool> isApplyNoCollInTimeRangeStrict{"isApplyNoCollInTimeRangeStrict", false, ""};
    Configurable<bool> isApplyNoCollInTimeRangeStandard{"isApplyNoCollInTimeRangeStandard", false, ""};
    Configurable<bool> isApplyNoCollInRofStrict{"isApplyNoCollInRofStrict", false, ""};
    Configurable<bool> isApplyNoCollInRofStandard{"isApplyNoCollInRofStandard", false, ""};
    Configurable<bool> isApplyNoHighMultCollInPrevRof{"isApplyNoHighMultCollInPrevRof", false, ""};
    Configurable<float> zVertexMaxInFilter{"zVertexMaxInFilter", 30.0f, "Accepted z-vertex range"};
    Configurable<float> zVertexMax{"zVertexMax", 10.0f, "Accepted z-vertex range"};
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", false, "Check event quality in run condition table"};
    Configurable<bool> requireCorrelationAnalysisRCTFlagChecker{"requireCorrelationAnalysisRCTFlagChecker", false, "Check event quality in run condition table for correlation analysis"};
    Configurable<bool> requireMcCollBestCollIndex{"requireMcCollBestCollIndex", false, "check for split vertices"};
    Configurable<std::string> setRCTFlagCheckerLabel{"setRCTFlagCheckerLabel", "CBT_muon_global", "Evt sel: RCT flag checker label"};
    Configurable<bool> requireRCTFlagCheckerLimitAcceptanceAsBad{"requireRCTFlagCheckerLimitAcceptanceAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
    Configurable<bool> requireZDCCheck{"requireZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
  } configCollision;

  //   configurables for MFT tracks
  struct : ConfigurableGroup {
    std::string prefix = "ConfigMft_group";
    Configurable<int> cutBestCollisionId{"cutBestCollisionId", 0, "cut on the best collision Id used in a filter"};
    Configurable<bool> cutFakeTracks{"cutFakeTracks", false, "if true, cut fake tracks using McMask"};
    Configurable<bool> cutOrphanTracksExplicitly{"cutOrphanTracksExplicitly", false, "if true, cut orphan tracks explicitly"};
    Configurable<bool> cutOnDcaXY{"cutOnDcaXY", false, "if true, cut on DCA XY"};
    Configurable<bool> cutOnDcaZ{"cutOnDcaZ", false, "if true, cut on DCA Z"};
    Configurable<float> dcaXYMax{"dcaXYMax", 2.0f, "Maximum value for DCA XY"};
    Configurable<float> dcaZMax{"dcaZMax", 10.0f, "Maximum value for DCA Z"};
    Configurable<float> etaMftTrackMax{"etaMftTrackMax", -2.4f, "Maximum value for the eta of MFT tracks when used in cut function"};
    Configurable<float> etaMftTrackMin{"etaMftTrackMin", -3.36f, "Minimum value for the eta of MFT tracks when used in cut function"};
    Configurable<float> etaMftTrackMaxFilter{"etaMftTrackMaxFilter", -2.0f, "Maximum value for the eta of MFT tracks when used in filter"};
    Configurable<float> etaMftTrackMinFilter{"etaMftTrackMinFilter", -3.9f, "Minimum value for the eta of MFT tracks when used in filter"};
    Configurable<float> mftMaxDCAxy{"mftMaxDCAxy", 2.0f, "Cut on dcaXY for MFT tracks"};
    Configurable<float> mftMaxDCAz{"mftMaxDCAz", 2.0f, "Cut on dcaZ for MFT tracks"};
    Configurable<int> nClustersMftTrack{"nClustersMftTrack", 5, "Minimum number of clusters for the reconstruction of MFT tracks"};
    Configurable<float> ptMftTrackMax{"ptMftTrackMax", 10.0f, "max value of MFT tracks pT when used in cut function"};
    Configurable<float> ptMftTrackMin{"ptMftTrackMin", 0.f, "min value of MFT tracks pT when used in cut function"};
    Configurable<bool> useMftPtCut{"useMftPtCut", false, "if true, use the Mft pt function cut"};
    Configurable<bool> useOnlyCATracks{"useOnlyCATracks", false, "if true, use strictly MFT tracks reconstructed with CA algo."};
    Configurable<bool> useOnlyLTFTracks{"useOnlyLTFTracks", false, "if true, use strictly MFT tracks reconstructed with LTF algo."};
  } configMft;

  float mZShift = 0;                                     // z-vertex shift
  float bZ = 0;                                          // Magnetic field for MFT
  static constexpr double CcenterMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  o2::parameters::GRPMagField* grpmag = nullptr;
  RCTFlagsChecker rctChecker;
  RCTFlagsChecker correlationAnalysisRctChecker{kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kMFTBad};
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hZVtxDiffAmbiguousTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hZVtxDiffNonAmbiguousTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hZVtxDiff2dReassociatedTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hZVtxDiffNot2dReassociatedTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hZVtxDiff3dReassociatedTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hZVtxDiffNot3dReassociatedTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hDcaAmbiguousTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hDcaNonAmbiguousTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hDca2dReassociatedTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hDcaNot2dReassociatedTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hDca3dReassociatedTracks;
  std::array<std::shared_ptr<THnSparse>, MatchedToTrueCollisionStep::NMatchedToTrueCollisionSteps> hDcaNot3dReassociatedTracks;

  // =========================
  //      using declarations : DATA
  // =========================

  using FilteredCollisionsWSelMult = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using FilteredMftTracks = soa::Filtered<aod::MFTTracks>;
  using FilteredMftTracksWColls = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>>;

  // =========================
  //      using declarations : MONTE-CARLO
  // =========================

  using FilteredCollisionsWSelMultMcLabels = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>>;
  using FilteredMftTracksWCollsMcLabels = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>>;
  // using MftReasso3dTracksWCollsMcLabels = soa::Join<aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;
  // using MftReasso2dTracksWCollsMcLabels = soa::Join<aod::BestCollisionsFwd, aod::McMFTTrackLabels>;
  // using MftReasso2dTracksWCollsMcLabels = soa::Join<aod::MFTTracks, aod::BestCollisionsFwd, aod::McMFTTrackLabels>;
  // using MftReasso3dTracksWCollsMcLabels = soa::Join<aod::MFTTracks, aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;
  using FilteredMcParticles = soa::Filtered<aod::McParticles>;

  // =========================
  //      Filters & partitions : DATA
  // =========================

  //  Collision filters
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < configCollision.zVertexMaxInFilter;

  Filter mftTrackEtaFilter = ((aod::fwdtrack::eta < configMft.etaMftTrackMaxFilter) && (aod::fwdtrack::eta > configMft.etaMftTrackMinFilter));

  Filter mftTrackCollisionIdFilter = (aod::fwdtrack::bestCollisionId >= 0);
  // Filter mftTrackDcaXYFilter = (nabs(aod::fwdtrack::bestDCAXY) < configMft.mftMaxDCAxy);
  // Filter mftTrackDcaZFilter = (nabs(aod::fwdtrack::bestDCAZ) < configMft.mftMaxDCAz);

  // =========================
  //      Filters & partitions : MONTE-CARLO
  // =========================

  Filter mcParticleEtaFilter = (aod::mcparticle::eta < configMft.etaMftTrackMaxFilter) && (aod::mcparticle::eta > configMft.etaMftTrackMinFilter);

  Partition<FilteredMcParticles> mcParticlesSample = (aod::mcparticle::eta < configMft.etaMftTrackMaxFilter) && (aod::mcparticle::eta > configMft.etaMftTrackMinFilter);

  // =========================
  //      Preslice : DATA
  // =========================

  Preslice<FilteredMftTracks> perColMftTracks = o2::aod::fwdtrack::collisionId;

  struct : ConfigurableGroup {
    std::string prefix = "ConfigAxis_group";
    ConfigurableAxis axisEta{"axisEta", {48, -4, -2}, "eta axis for MFT histograms"};
    ConfigurableAxis axisDcaX{"axisDcaX", {800, -1., 1.}, "DCAx binning (cm)"};
    ConfigurableAxis axisDcaY{"axisDcaY", {800, -1., 1.}, "DCAy binning (cm)"};
    ConfigurableAxis axisDcaZ{"axisDcaZ", {800, -1., 1.}, "DCAz binning (cm)"};
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity axis for histograms"};
    ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "phi axis for histograms"};
    ConfigurableAxis axisPt{"axisPt", {72, 0, 36}, "pt axis for histograms"};
    ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
    ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {1, -10, 10}, "vertex axis for efficiency histograms"};
  } configAxis;

  HistogramRegistry registry{"registry"};

  template <DataType DataType>
  void addMftHistograms()
  {
    registry.add(Form("%shAmbiguityOfMftTracks", WhatDataType[DataType].data()), "hAmbiguityOfMftTracks", {HistType::kTH1D, {{MftTrackAmbiguityStep::NMftAmbiguitySteps, -0.5, +MftTrackAmbiguityStep::NMftAmbiguitySteps - 0.5}}});
    std::string labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NMftAmbiguitySteps];
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::AllMftTracks] = "all MFT tracks";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::AfterTrackSelection] = "MFT tracks after selection";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NumberOfAmbiguousTracks] = "how much tracks are ambigous";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks] = "how much tracks are non-ambiguous";
    registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hAmbiguityOfMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftTrackAmbiguityStep::NMftAmbiguitySteps; iBin++) {
      registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hAmbiguityOfMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsAmbiguityOfMftTracks[iBin].data());
    }

    registry.add(Form("%shMftTracksSelection", WhatDataType[DataType].data()), "hMftTracksSelection", {HistType::kTH1D, {{MftTrackSelectionStep::NMftTrackSelectionSteps, -0.5, +MftTrackSelectionStep::NMftTrackSelectionSteps - 0.5}}});
    std::string labelsMftTracksSelection[MftTrackSelectionStep::NMftTrackSelectionSteps];
    labelsMftTracksSelection[MftTrackSelectionStep::NoSelection] = "all MFT tracks";
    labelsMftTracksSelection[MftTrackSelectionStep::Eta] = "MFT tracks after eta selection";
    labelsMftTracksSelection[MftTrackSelectionStep::Cluster] = "MFT tracks after clusters selection";
    labelsMftTracksSelection[MftTrackSelectionStep::Pt] = "MFT tracks after pT selection";
    labelsMftTracksSelection[MftTrackSelectionStep::IsCA] = "MFT tracks reconstructed with CA";
    labelsMftTracksSelection[MftTrackSelectionStep::IsLTF] = "MFT tracks reconstructed with LTF";
    registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hMftTracksSelection"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftTrackSelectionStep::NMftTrackSelectionSteps; iBin++) {
      registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hMftTracksSelection"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftTracksSelection[iBin].data());
    }

    registry.add(Form("%shReassociation2dMftTracks", WhatDataType[DataType].data()), "hReassociation2dMftTracks", {HistType::kTH1D, {{Reassociation2dMftTracks::NReassociation2dMftTracksSteps, -0.5, +Reassociation2dMftTracks::NReassociation2dMftTracksSteps - 0.5}}});
    std::string labelsReassociation2dMftTracks[Reassociation2dMftTracks::NReassociation2dMftTracksSteps];
    labelsReassociation2dMftTracks[Reassociation2dMftTracks::AllAmbiguousTracksAfterTrackSelectionsFor2d] = "Ambiguous MFT tracks after track selection";
    labelsReassociation2dMftTracks[Reassociation2dMftTracks::NotReassociated2dMftTracks] = "Not reassociated MFT tracks by DCAxy method";
    labelsReassociation2dMftTracks[Reassociation2dMftTracks::Reassociated2dMftTracks] = "Reassociated MFT tracks by DCAxy method";
    registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hReassociation2dMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < Reassociation2dMftTracks::NReassociation2dMftTracksSteps; iBin++) {
      registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hReassociation2dMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsReassociation2dMftTracks[iBin].data());
    }

    registry.add(Form("%shReassociation3dMftTracks", WhatDataType[DataType].data()), "hReassociation3dMftTracks", {HistType::kTH1D, {{Reassociation3dMftTracks::NReassociation3dMftTracksSteps, -0.5, +Reassociation3dMftTracks::NReassociation3dMftTracksSteps - 0.5}}});
    std::string labelsReassociation3dMftTracks[Reassociation3dMftTracks::NReassociation3dMftTracksSteps];
    labelsReassociation3dMftTracks[Reassociation3dMftTracks::AllAmbiguousTracksAfterTrackSelectionsFor3d] = "Ambiguous MFT tracks after track selection";
    labelsReassociation3dMftTracks[Reassociation3dMftTracks::NotReassociated3dMftTracks] = "Not reassociated MFT tracks by DCAxyz method";
    labelsReassociation3dMftTracks[Reassociation3dMftTracks::Reassociated3dMftTracks] = "Reassociated MFT tracks by DCAxyz method";
    registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hReassociation3dMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < Reassociation3dMftTracks::NReassociation3dMftTracksSteps; iBin++) {
      registry.get<TH1>(HIST(WhatDataType[DataType]) + HIST("hReassociation3dMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsReassociation3dMftTracks[iBin].data());
    }
  }

  void addMftMonteCarloHistograms()
  {
    // AmbiguousTracks ZVtxDiff dist. (contains matched to true and not matched to true collision)
    hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hZVtxDiffAmbiguousTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    // AmbiguousTracks NOT matched to true collisions ZVtxDiff dist.
    hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffAmbiguousTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    // AmbiguousTracks matched to true collisions ZVtxDiff dist.
    hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffAmbiguousTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});

    hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hZVtxDiffNonAmbiguousTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffNonAmbiguousTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffNonAmbiguousTrackMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});

    hZVtxDiff2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hZVtxDiff2dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiff2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiff2dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiff2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiff2dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});

    hZVtxDiffNot2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hZVtxDiffNot2dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiffNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffNot2dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiffNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffNot2dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});

    hZVtxDiff3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hZVtxDiff3dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiff3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiff3dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiff3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiff3dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});

    hZVtxDiffNot3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hZVtxDiffNot3dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiffNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffNot3dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});
    hZVtxDiffNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hZVtxDiffNot3dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{X}_{vtx}^{reco}#minus#it{X}_{vtx}^{gen} (cm);#it{Y}_{vtx}^{reco}#minus#it{Y}_{vtx}^{gen} (cm);#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaY, configAxis.axisDcaZ});

    hDcaAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hDcaAmbiguousTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDcaAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaAmbiguousTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDcaAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaAmbiguousTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});

    hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hDcaNonAmbiguousTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaNonAmbiguousTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaNonAmbiguousTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});

    hDca2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hDca2dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{XY}^{gen} (cm);", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaX});
    hDca2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDca2dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{XY}^{gen} (cm);", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaX});
    hDca2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDca2dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{XY}^{gen} (cm);", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaX});

    hDcaNot2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hDcaNot2dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{XY}^{gen} (cm);", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaX});
    hDcaNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaNot2dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{XY}^{gen} (cm);", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaX});
    hDcaNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaNot2dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{XY}^{gen} (cm);", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaX});

    hDca3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hDca3dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDca3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDca3dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDca3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDca3dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});

    hDcaNot3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks] = registry.add<THnSparse>("MC/hDcaNot3dReassociatedTracks", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDcaNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaNot3dReassociatedTracksNotMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});
    hDcaNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision] = registry.add<THnSparse>("MC/hDcaNot3dReassociatedTracksMatchedToTrueCollision", ";#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};DCA_{XY}^{reco} (cm);  DCA_{Z}^{reco} (cm); DCA_{XY}^{gen} (cm);  DCA_{Z}^{gen} (cm)", HistType::kTHnSparseF, {configAxis.axisPt, configAxis.axisEta, configAxis.axisDcaX, configAxis.axisDcaZ, configAxis.axisDcaX, configAxis.axisDcaZ});

    registry.add("MC/hIsAmbiguousTrackMatchedToTrueCollision", "hIsAmbiguousTrackMatchedToTrueCollision", {HistType::kTH1D, {{MftAmbiguousAndMatchedToTrueCollisionStep::NMftAmbiguousAndMatchedToTrueCollisionSteps, -0.5, +MftAmbiguousAndMatchedToTrueCollisionStep::NMftAmbiguousAndMatchedToTrueCollisionSteps - 0.5}}});
    std::string labelsMftAmbiguousAndMatchedToTrueCollisionStep[MftAmbiguousAndMatchedToTrueCollisionStep::NMftAmbiguousAndMatchedToTrueCollisionSteps];
    labelsMftAmbiguousAndMatchedToTrueCollisionStep[MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguous] = "number of MFT ambiguous tracks";
    labelsMftAmbiguousAndMatchedToTrueCollisionStep[MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguousAndMatchedToTrueCollision] = "number of MFT ambiguous tracks matched to true collision";
    labelsMftAmbiguousAndMatchedToTrueCollisionStep[MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguousAndNotMatchedToTrueCollision] = "number of MFT ambiguous tracks NOT matched to true collision";
    registry.get<TH1>(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftAmbiguousAndMatchedToTrueCollisionStep::NMftAmbiguousAndMatchedToTrueCollisionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftAmbiguousAndMatchedToTrueCollisionStep[iBin].data());
    }

    registry.add("MC/hIsNonAmbiguousTrackMatchedToTrueCollision", "hIsNonAmbiguousTrackMatchedToTrue", {HistType::kTH1D, {{MftNonAmbiguousAndMatchedToTrueCollisionStep::NMftNonAmbiguousAndMatchedToTrueCollisionSteps, -0.5, +MftNonAmbiguousAndMatchedToTrueCollisionStep::NMftNonAmbiguousAndMatchedToTrueCollisionSteps - 0.5}}});
    std::string labelsMftNonAmbiguousAndMatchedToTrueCollisionStep[MftNonAmbiguousAndMatchedToTrueCollisionStep::NMftNonAmbiguousAndMatchedToTrueCollisionSteps];
    labelsMftNonAmbiguousAndMatchedToTrueCollisionStep[MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguous] = "number of MFT Non ambiguous tracks";
    labelsMftNonAmbiguousAndMatchedToTrueCollisionStep[MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguousAndMatchedToTrueCollision] = "number of MFT Non ambiguous tracks matched to true collision";
    labelsMftNonAmbiguousAndMatchedToTrueCollisionStep[MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguousAndNotMatchedToTrueCollision] = "number of MFT Non ambiguous tracks NOT matched to true collision";
    registry.get<TH1>(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftNonAmbiguousAndMatchedToTrueCollisionStep::NMftNonAmbiguousAndMatchedToTrueCollisionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftNonAmbiguousAndMatchedToTrueCollisionStep[iBin].data());
    }

    registry.add("MC/hIs2dReassociatedAndMatchedToTrueCollision", "Is2dReassociatedAndMatchedToTrueCollision", {HistType::kTH1D, {{Mft2dReassociatedAndMatchedToTrueCollisionStep::NMft2dReassociatedAndMatchedToTrueCollisionSteps, -0.5, +Mft2dReassociatedAndMatchedToTrueCollisionStep::NMft2dReassociatedAndMatchedToTrueCollisionSteps - 0.5}}});
    std::string labelsMft2dReassociatedAndMatchedToTrueCollisionStep[Mft2dReassociatedAndMatchedToTrueCollisionStep::NMft2dReassociatedAndMatchedToTrueCollisionSteps];
    labelsMft2dReassociatedAndMatchedToTrueCollisionStep[Mft2dReassociatedAndMatchedToTrueCollisionStep::Is2dReassociated] = "number of MFT 2d reassociated tracks";
    labelsMft2dReassociatedAndMatchedToTrueCollisionStep[Mft2dReassociatedAndMatchedToTrueCollisionStep::Is2dReassociatedAndMatchedToTrueCollision] = "number of MFT 2d reassociated tracks matched to true collision";
    labelsMft2dReassociatedAndMatchedToTrueCollisionStep[Mft2dReassociatedAndMatchedToTrueCollisionStep::Is2dReassociatedAndNotMatchedToTrueCollision] = "number of MFT 2d reassociated tracks NOT matched to true collision";
    registry.get<TH1>(HIST("MC/hIs2dReassociatedAndMatchedToTrueCollision"))->SetMinimum(0);

    for (int iBin = 0; iBin < Mft2dReassociatedAndMatchedToTrueCollisionStep::NMft2dReassociatedAndMatchedToTrueCollisionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hIs2dReassociatedAndMatchedToTrueCollision"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMft2dReassociatedAndMatchedToTrueCollisionStep[iBin].data());
    }

    registry.add("MC/hIsNot2dReassociatedAndMatchedToTrueCollision", "IsNot2dReassociatedAndMatchedToTrueCollision", {HistType::kTH1D, {{MftNot2dReassociatedAndMatchedToTrueCollisionStep::NMftNot2dReassociatedAndMatchedToTrueCollisionSteps, -0.5, +MftNot2dReassociatedAndMatchedToTrueCollisionStep::NMftNot2dReassociatedAndMatchedToTrueCollisionSteps - 0.5}}});
    std::string labelsMftNot2dReassociatedAndMatchedToTrueCollisionStep[MftNot2dReassociatedAndMatchedToTrueCollisionStep::NMftNot2dReassociatedAndMatchedToTrueCollisionSteps];
    labelsMftNot2dReassociatedAndMatchedToTrueCollisionStep[MftNot2dReassociatedAndMatchedToTrueCollisionStep::IsNot2dReassociated] = "number of MFT NOT 2d reassociated tracks";
    labelsMftNot2dReassociatedAndMatchedToTrueCollisionStep[MftNot2dReassociatedAndMatchedToTrueCollisionStep::IsNot2dReassociatedAndMatchedToTrueCollision] = "number of MFT NOT 2d reassociated tracks matched to true collision";
    labelsMftNot2dReassociatedAndMatchedToTrueCollisionStep[MftNot2dReassociatedAndMatchedToTrueCollisionStep::IsNot2dReassociatedAndNotMatchedToTrueCollision] = "number of MFT NOT 2d reassociated tracks NOT matched to true collision";
    registry.get<TH1>(HIST("MC/hIsNot2dReassociatedAndMatchedToTrueCollision"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftNot2dReassociatedAndMatchedToTrueCollisionStep::NMftNot2dReassociatedAndMatchedToTrueCollisionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hIsNot2dReassociatedAndMatchedToTrueCollision"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftNot2dReassociatedAndMatchedToTrueCollisionStep[iBin].data());
    }

    registry.add("MC/hIs3dReassociatedAndMatchedToTrueCollision", "Is3dReassociatedAndMatchedToTrueCollision", {HistType::kTH1D, {{Mft3dReassociatedAndMatchedToTrueCollisionStep::NMft3dReassociatedAndMatchedToTrueCollisionSteps, -0.5, +Mft3dReassociatedAndMatchedToTrueCollisionStep::NMft3dReassociatedAndMatchedToTrueCollisionSteps - 0.5}}});
    std::string labelsMft3dReassociatedAndMatchedToTrueCollisionStep[Mft3dReassociatedAndMatchedToTrueCollisionStep::NMft3dReassociatedAndMatchedToTrueCollisionSteps];
    labelsMft3dReassociatedAndMatchedToTrueCollisionStep[Mft3dReassociatedAndMatchedToTrueCollisionStep::Is3dReassociated] = "number of MFT 3d reassociated tracks";
    labelsMft3dReassociatedAndMatchedToTrueCollisionStep[Mft3dReassociatedAndMatchedToTrueCollisionStep::Is3dReassociatedAndMatchedToTrueCollision] = "number of MFT 3d reassociated tracks matched to true collision";
    labelsMft3dReassociatedAndMatchedToTrueCollisionStep[Mft3dReassociatedAndMatchedToTrueCollisionStep::Is3dReassociatedAndNotMatchedToTrueCollision] = "number of MFT 3d reassociated tracks NOT matched to true collision";
    registry.get<TH1>(HIST("MC/hIs3dReassociatedAndMatchedToTrueCollision"))->SetMinimum(0);

    for (int iBin = 0; iBin < Mft3dReassociatedAndMatchedToTrueCollisionStep::NMft3dReassociatedAndMatchedToTrueCollisionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hIs3dReassociatedAndMatchedToTrueCollision"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMft3dReassociatedAndMatchedToTrueCollisionStep[iBin].data());
    }

    registry.add("MC/hIsNot3dReassociatedAndMatchedToTrueCollision", "IsNot3dReassociatedAndMatchedToTrueCollision", {HistType::kTH1D, {{MftNot3dReassociatedAndMatchedToTrueCollisionStep::NMftNot3dReassociatedAndMatchedToTrueCollisionSteps, -0.5, +MftNot3dReassociatedAndMatchedToTrueCollisionStep::NMftNot3dReassociatedAndMatchedToTrueCollisionSteps - 0.5}}});
    std::string labelsMftNot3dReassociatedAndMatchedToTrueCollisionStep[MftNot3dReassociatedAndMatchedToTrueCollisionStep::NMftNot3dReassociatedAndMatchedToTrueCollisionSteps];
    labelsMftNot3dReassociatedAndMatchedToTrueCollisionStep[MftNot3dReassociatedAndMatchedToTrueCollisionStep::IsNot3dReassociated] = "number of MFT NOT 3d reassociated tracks";
    labelsMftNot3dReassociatedAndMatchedToTrueCollisionStep[MftNot3dReassociatedAndMatchedToTrueCollisionStep::IsNot3dReassociatedAndMatchedToTrueCollision] = "number of MFT NOT 3d reassociated tracks matched to true collision";
    labelsMftNot3dReassociatedAndMatchedToTrueCollisionStep[MftNot3dReassociatedAndMatchedToTrueCollisionStep::IsNot3dReassociatedAndNotMatchedToTrueCollision] = "number of MFT NOT 3d reassociated tracks NOT matched to true collision";
    registry.get<TH1>(HIST("MC/hIsNot3dReassociatedAndMatchedToTrueCollision"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftNot3dReassociatedAndMatchedToTrueCollisionStep::NMftNot3dReassociatedAndMatchedToTrueCollisionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hIsNot3dReassociatedAndMatchedToTrueCollision"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftNot3dReassociatedAndMatchedToTrueCollisionStep[iBin].data());
    }

    registry.add("MC/hIsTrueCollisionAmongCompatibleCollisions", "IsTrueCollisionAmongCompatibleCollisions", {HistType::kTH1D, {{MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps, -0.5, +MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps - 0.5}}});
    registry.add("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated", "hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated", {HistType::kTH1D, {{MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps, -0.5, +MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps - 0.5}}});
    registry.add("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated", "hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated", {HistType::kTH1D, {{MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps, -0.5, +MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps - 0.5}}});
    std::string labelsMftIsTrueCollisionAmongCompatibleCollisionsStep[MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps];
    labelsMftIsTrueCollisionAmongCompatibleCollisionsStep[MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks] = "number of wrongly associated tracks";
    labelsMftIsTrueCollisionAmongCompatibleCollisionsStep[MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions] = "number of MFT tracks with True Coll. among Compatible";
    labelsMftIsTrueCollisionAmongCompatibleCollisionsStep[MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions] = "number of MFT tracks WITHOUT True Coll. among Compatible";
    registry.get<TH1>(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"))->SetMinimum(0);
    registry.get<TH1>(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"))->SetMinimum(0);
    registry.get<TH1>(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftIsTrueCollisionAmongCompatibleCollisionsStep::NMftIsTrueCollisionAmongCompatibleCollisionsSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftIsTrueCollisionAmongCompatibleCollisionsStep[iBin].data());
      registry.get<TH1>(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftIsTrueCollisionAmongCompatibleCollisionsStep[iBin].data());
      registry.get<TH1>(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftIsTrueCollisionAmongCompatibleCollisionsStep[iBin].data());
    }
  }

  //  =========================
  //      init()
  //  =========================
  void init(InitContext&)
  {
    ccdb->setURL(configCcdb.ccdbUrl);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    rctChecker.init(configCollision.setRCTFlagCheckerLabel, configCollision.requireZDCCheck, configCollision.requireRCTFlagCheckerLimitAcceptanceAsBad, true);
    correlationAnalysisRctChecker.init({kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kMFTBad}, configCollision.requireZDCCheck, configCollision.requireRCTFlagCheckerLimitAcceptanceAsBad, true);

    //  =========================
    //      Event histograms
    //  =========================

    registry.add("Data/hVtxZ", "v_{z} (cm)", {HistType::kTH1D, {configAxis.axisVertex}});
    // registry.add("Data/hNTracks", "", {HistType::kTH1F, {configAxis.axisMultiplicity}});
    registry.add(Form("Data/hMultiplicity_%s", WhatMultiplicityEstimator[configCollision.multiplicityEstimator].data()), "", {HistType::kTH1D, {configAxis.axisMultiplicity}});

    registry.add("hPreciseEventCounter", "hPreciseEventCounter", {HistType::kTH1D, {{SpecificEventSelectionStep::NSpecificEventSelectionSteps, -0.5, +SpecificEventSelectionStep::NSpecificEventSelectionSteps - 0.5}}});
    std::string labels[SpecificEventSelectionStep::NSpecificEventSelectionSteps];
    labels[SpecificEventSelectionStep::AllEventsPrecise] = "all";
    labels[SpecificEventSelectionStep::HasMcCollision] = "has MC coll?";
    labels[SpecificEventSelectionStep::IsNotSplitVertex] = "Is not split vertex (BestCollisionIndex)";
    labels[SpecificEventSelectionStep::IsSel8] = "sel8";
    labels[SpecificEventSelectionStep::IsNoSameBunchPileup] = "IsNoSameBunchPileup";
    labels[SpecificEventSelectionStep::IsGoodItsLayersAll] = "IsGoodItsLayersAll";
    labels[SpecificEventSelectionStep::IsGoodZvtxFT0vsPV] = "IsGoodZvtxFT0vsPV";
    labels[SpecificEventSelectionStep::IsNoCollInRofStandard] = "IsNoCollInRofStandard";
    labels[SpecificEventSelectionStep::IsNoCollInRofStrict] = "IsNoCollInRofStrict";
    labels[SpecificEventSelectionStep::IsNoCollInTimeRangeStandard] = "IsNoCollInTimeRangeStandard";
    labels[SpecificEventSelectionStep::IsNoCollInTimeRangeStrict] = "IsNoCollInTimeRangeStrict";
    labels[SpecificEventSelectionStep::IsNoHighMultCollInPrevRof] = "IsNoHighMultCollInPrevRof";
    labels[SpecificEventSelectionStep::IsRctFlagChecked] = "IsRctFlagChecked";
    labels[SpecificEventSelectionStep::IsWithinZvtxWindow] = "IsWithinZvtxWindow";
    registry.get<TH1>(HIST("hPreciseEventCounter"))->SetMinimum(0);

    for (int iBin = 0; iBin < SpecificEventSelectionStep::NSpecificEventSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("hPreciseEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    registry.add("MC/hPreciseTrackSelectionCounter", "hPreciseTrackSelectionCounter", {HistType::kTH1D, {{SpecificTrackSelectionStep::NSpecificTrackSelectionSteps, -0.5, +SpecificTrackSelectionStep::NSpecificTrackSelectionSteps - 0.5}}});
    std::string labelsTrackSelection[SpecificTrackSelectionStep::NSpecificTrackSelectionSteps];
    labelsTrackSelection[SpecificTrackSelectionStep::AllTracksPrecise] = "all tracks";
    labelsTrackSelection[SpecificTrackSelectionStep::IsTrueTrack] = "true tracks";
    labelsTrackSelection[SpecificTrackSelectionStep::AfterOrphanCut] = "after orphan cut";
    labelsTrackSelection[SpecificTrackSelectionStep::AfterEtaCut] = "after eta cut";
    labelsTrackSelection[SpecificTrackSelectionStep::AfterClusterCut] = "after cluster cut";
    labelsTrackSelection[SpecificTrackSelectionStep::AfterPtCut] = "after pt cut";
    labelsTrackSelection[SpecificTrackSelectionStep::AfterDcaXYCut] = "after DCA XY cut";
    labelsTrackSelection[SpecificTrackSelectionStep::AfterDcaZCut] = "after DCA Z cut";
    labelsTrackSelection[SpecificTrackSelectionStep::IsCATrack] = "is CA track";
    labelsTrackSelection[SpecificTrackSelectionStep::IsLTFTrack] = "is LTF track";
    labelsTrackSelection[SpecificTrackSelectionStep::HasMcParticle] = "has MC particle";
    registry.get<TH1>(HIST("MC/hPreciseTrackSelectionCounter"))->SetMinimum(0);

    for (int iBin = 0; iBin < SpecificTrackSelectionStep::NSpecificTrackSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hPreciseTrackSelectionCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labelsTrackSelection[iBin].data());
    }

    registry.add("MC/hTrackAmbiguityCheck", "hTrackAmbiguityCheck", {HistType::kTH1D, {{TrackAmbiguityCheckStep::NTrackAmbiguityCheckSteps, -0.5, +TrackAmbiguityCheckStep::NTrackAmbiguityCheckSteps - 0.5}}});
    std::string labelsTrackAmbiguityCheck[TrackAmbiguityCheckStep::NTrackAmbiguityCheckSteps];
    labelsTrackAmbiguityCheck[TrackAmbiguityCheckStep::AllTracksCheck] = "all tracks";
    labelsTrackAmbiguityCheck[TrackAmbiguityCheckStep::IsAmbDegreeEqualToCompatibleCollIdsSize] = "ambDegree == compatibleCollIds.size";
    registry.get<TH1>(HIST("MC/hTrackAmbiguityCheck"))->SetMinimum(0);

    for (int iBin = 0; iBin < TrackAmbiguityCheckStep::NTrackAmbiguityCheckSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hTrackAmbiguityCheck"))->GetXaxis()->SetBinLabel(iBin + 1, labelsTrackAmbiguityCheck[iBin].data());
    }

    registry.add("MC/hMonteCarloEventCounter", "hMonteCarloEventCounter", {HistType::kTH1D, {{MonteCarloEventSelectionStep::NMonteCarloEventSelectionSteps, -0.5, +MonteCarloEventSelectionStep::NMonteCarloEventSelectionSteps - 0.5}}});
    std::string labelsMonteCarloEvents[MonteCarloEventSelectionStep::NMonteCarloEventSelectionSteps];
    labelsMonteCarloEvents[MonteCarloEventSelectionStep::AllMonteCarloEvents] = "all collisions";
    labelsMonteCarloEvents[MonteCarloEventSelectionStep::MonteCarloEventsAfterEventSelection] = "collisions after event selection";
    labelsMonteCarloEvents[MonteCarloEventSelectionStep::HasMonteCarloCollision] = "has MC collision";
    labelsMonteCarloEvents[MonteCarloEventSelectionStep::HasNotMonteCarloCollision] = "has not MC collision";
    registry.get<TH1>(HIST("MC/hMonteCarloEventCounter"))->SetMinimum(0);

    for (int iBin = 0; iBin < MonteCarloEventSelectionStep::NMonteCarloEventSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hMonteCarloEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMonteCarloEvents[iBin].data());
    }

    registry.add("MC/hMonteCarloTrackCounter", "hMonteCarloTrackCounter", {HistType::kTH1D, {{MonteCarloTrackSelectionStep::NMonteCarloTrackSelectionSteps, -0.5, +MonteCarloTrackSelectionStep::NMonteCarloTrackSelectionSteps - 0.5}}});
    std::string labelsMonteCarloTracks[MonteCarloTrackSelectionStep::NMonteCarloTrackSelectionSteps];
    labelsMonteCarloTracks[MonteCarloTrackSelectionStep::AllMonteCarloTracks] = "all tracks";
    labelsMonteCarloTracks[MonteCarloTrackSelectionStep::MonteCarloTracksAfterTrackSelection] = "tracks after track selection";
    labelsMonteCarloTracks[MonteCarloTrackSelectionStep::HasMonteCarloParticle] = "has MC particle";
    labelsMonteCarloTracks[MonteCarloTrackSelectionStep::HasNotMonteCarloParticle] = "has not MC particle";
    registry.get<TH1>(HIST("MC/hMonteCarloTrackCounter"))->SetMinimum(0);

    for (int iBin = 0; iBin < MonteCarloTrackSelectionStep::NMonteCarloTrackSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("MC/hMonteCarloTrackCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMonteCarloTracks[iBin].data());
    }

    //  =========================
    //      Process functions initialization
    //  =========================

    if (doprocessData) {
      addMftHistograms<Data>();
    }

    if (doprocessMcReassociated2d) {
      addMftHistograms<Mc>();
      addMftMonteCarloHistograms();
    }

    if (doprocessMcReassociated3d) {
      addMftHistograms<Mc>();
      addMftMonteCarloHistograms();
    }

  } // End of init() function

  // =========================
  //      Helper functions
  // =========================

  template <typename TCollision>
  float getMultiplicityEstimator(TCollision collision, bool isSameEvent)
  {
    switch (configCollision.multiplicityEstimator) {
      case MultiplicityEstimators::MultNTracksPV:
        if (isSameEvent) {
          registry.fill(HIST("Data/hMultiplicity_multNTracksPV"), collision.multNTracksPV());
        }
        return collision.multNTracksPV();
      case MultiplicityEstimators::MultNumContrib:
        if (isSameEvent) {
          registry.fill(HIST("Data/hMultiplicity_multNumContrib"), collision.numContrib());
        }
        return collision.numContrib();
      case MultiplicityEstimators::MultFT0C:
        if (isSameEvent) {
          registry.fill(HIST("Data/hMultiplicity_multFT0C"), collision.multFT0C());
        }
        return collision.multFT0C();
      case MultiplicityEstimators::MultFT0M:
        if (isSameEvent) {
          registry.fill(HIST("Data/hMultiplicity_multFT0M"), collision.multFT0M());
        }
        return collision.multFT0M();
      default:
        return collision.multNTracksPV();
    }
  }

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TTrack>
  bool isTrueCollisionAmongCompatibleCollisions(TTrack track)
  {

    auto const& compatibleIds = track.compatibleCollIds();

    // For this track, check if *any* reco collision that corresponds to its TRUE MC collision
    // is present among the compatible collisions.
    bool recoOfTrueInCompatible = false;
    if (track.ambDegree() != 0) {
      const int mcTrueCollisionId = track.mcParticle().mcCollisionId();
      // Fast membership test using compatibleCollIds() if available.
      // If compatibleIds is empty, it means the table is missing or no compatible collisions were stored.
      if (!compatibleIds.empty()) {
        for (auto const& id : compatibleIds) {

          auto iterator = recoMcCollisionId.find(id);

          if (iterator != recoMcCollisionId.end()) {
            if (iterator->second == mcTrueCollisionId) {

              recoOfTrueInCompatible = true;
              return recoOfTrueInCompatible;
            }
          } else {
            return recoOfTrueInCompatible;
          }
        }
      }
    }
    return recoOfTrueInCompatible;
  }

  template <typename TBc>
  void loadZVertexShiftCorrection(TBc const& bc)
  {

    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(configTask.grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    bZ = field->getBz(CcenterMFT);
    LOG(info) << "The field at the center of the MFT is bZ = " << bZ;

    if (configTask.cfgApplyZShiftFromCCDB) {
      auto* zShift = ccdb->getForTimeStamp<std::vector<float>>(configTask.cfgZShiftPath, bc.timestamp());
      if (zShift != nullptr && !zShift->empty()) {
        LOGF(info, "reading z shift %f from %s", (*zShift)[0], configTask.cfgZShiftPath.value);
        mZShift = (*zShift)[0];
      } else {
        LOGF(info, "z shift is not found in ccdb path %s. set to 0 cm", configTask.cfgZShiftPath.value);
        mZShift = 0;
      }
    } else {
      LOGF(info, "z shift is manually set to %f cm", configTask.cfgManualZShift.value);
      mZShift = configTask.cfgManualZShift;
    }
  }

  // =========================
  //      Cuts with functions
  // =========================

  //  FIXME: Some collisions are rejected here, what causes (part of) differences with the D0 task
  template <typename TCollision>
  bool isAcceptedCollision(TCollision const& collision, bool fillHistograms = false)
  {

    if (!collision.sel8()) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsSel8);
    }
    if (configCollision.isApplySameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNoSameBunchPileup);
    }
    if (configCollision.isApplyGoodItsLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsGoodItsLayersAll);
    }
    if (configCollision.isApplyGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsGoodZvtxFT0vsPV);
    }
    if (configCollision.isApplyNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNoCollInRofStandard);
    }
    if (configCollision.isApplyNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNoCollInRofStrict);
    }
    if (configCollision.isApplyNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNoCollInTimeRangeStandard);
    }
    if (configCollision.isApplyNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNoCollInTimeRangeStrict);
    }
    if (configCollision.isApplyNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNoHighMultCollInPrevRof);
    }
    if (configCollision.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if (configCollision.requireCorrelationAnalysisRCTFlagChecker && !correlationAnalysisRctChecker(collision)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsRctFlagChecked);
    }
    if (collision.posZ() > configCollision.zVertexMax || collision.posZ() < (-configCollision.zVertexMax)) {
      return false;
    }
    if (fillHistograms) {
      registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsWithinZvtxWindow);
    }

    registry.fill(HIST("Data/hVtxZ"), collision.posZ());

    return true;
  }

  //  TODO: Check how to put this into a Filter
  // I tried to put it as a filter, but filters for normal TPC tracks also apply to MFT tracks I think
  // and it seems that they are not compatible
  template <typename TTrack>
  bool isAcceptedMftTrack(TTrack const& mftTrack, bool fillHistograms, bool isData, float dcaXY, float dcaZ)
  {
    // cut on the eta of MFT tracks
    if (mftTrack.eta() > configMft.etaMftTrackMax || mftTrack.eta() < configMft.etaMftTrackMin) {
      return false;
    }

    if (fillHistograms) {
      if (isData) {
        registry.fill(HIST("Data/hMftTracksSelection"), MftTrackSelectionStep::Eta);
      } else {
        registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AfterEtaCut);
      }
    }

    // cut on the number of clusters of the reconstructed MFT track
    if (mftTrack.nClusters() < configMft.nClustersMftTrack) {
      return false;
    }

    if (fillHistograms) {
      if (isData) {
        registry.fill(HIST("Data/hMftTracksSelection"), MftTrackSelectionStep::Cluster);
      } else {
        registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AfterClusterCut);
      }
    }

    // cut on the pT of MFT tracks (for test purposes)
    if (configMft.useMftPtCut && (mftTrack.pt() > configMft.ptMftTrackMax || mftTrack.pt() < configMft.ptMftTrackMin)) {
      return false;
    }

    if (fillHistograms) {
      if (isData) {
        registry.fill(HIST("Data/hMftTracksSelection"), MftTrackSelectionStep::Pt);
      } else {
        registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AfterPtCut);
      }
    }

    if (configMft.cutOnDcaXY && std::abs(dcaXY) > configMft.dcaXYMax) {
      return false;
    }
    registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AfterDcaXYCut);
    if (configMft.cutOnDcaZ && std::abs(dcaZ) > configMft.dcaZMax) {
      return false;
    }
    registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AfterDcaZCut);

    // cut on the track algorithm of MFT tracks
    if (mftTrack.isCA()) {
      if (fillHistograms) {
        if (isData) {
          registry.fill(HIST("Data/hMftTracksSelection"), MftTrackSelectionStep::IsCA);
        } else {
          registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::IsCATrack);
        }
      }

      if (configMft.useOnlyLTFTracks) {
        return false;
      }

    } else {
      if (fillHistograms) {
        if (isData) {
          registry.fill(HIST("Data/hMftTracksSelection"), MftTrackSelectionStep::IsLTF);
        } else {
          registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::IsLTFTrack);
        }
      }

      if (configMft.useOnlyCATracks) {
        return false;
      }
    }

    return true;
  }

  // Cut on ambiguous MFT tracks
  template <typename TTrack>
  bool isAmbiguousMftTrack(TTrack const& mftTrack, bool fillHistograms)
  {
    if (mftTrack.ambDegree() > 1) {
      if (fillHistograms) {
        registry.fill(HIST("Data/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfAmbiguousTracks);
      }
      return true;
      if (mftTrack.ambDegree() > 1) {
        if (fillHistograms) {
          registry.fill(HIST("Data/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfAmbiguousTracks);
        }
        return true;
      }

      if (fillHistograms) {
        registry.fill(HIST("Data/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks);
      }
      return false;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks);
    }
    return false;
  }

  //============================================================================================
  //
  // PROCESS FUNCTIONS
  //
  //============================================================================================

  void processData(FilteredCollisionsWSelMult::iterator const& collision,
                   FilteredMftTracks const& /*mftTracks*/,
                   soa::SmallGroups<aod::BestCollisionsFwd> const& reassociated2dMftTracks,
                   aod::BCsWithTimestamps const&)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    loadZVertexShiftCorrection(bc);

    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    // const auto multiplicity = getMultiplicityEstimator(collision, true);

    for (const auto& reassociated2dMftTrack : reassociated2dMftTracks) {

      registry.fill(HIST("Data/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
      auto templatedMftTrack = reassociated2dMftTrack.template mfttrack_as<FilteredMftTracks>();

      o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(templatedMftTrack, mZShift);
      std::array<double, 3> dcaInfOrig;
      trackPar.propagateToDCAhelix(bZ, {collision.posX(), collision.posY(), collision.posZ()}, dcaInfOrig);
      auto dcaXYoriginal = 999.f;
      dcaXYoriginal = std::sqrt(dcaInfOrig[0] * dcaInfOrig[0] + dcaInfOrig[1] * dcaInfOrig[1]);

      if (!isAcceptedMftTrack(templatedMftTrack, false, true, dcaXYoriginal, dcaInfOrig[2])) {
        continue;
      }

      registry.fill(HIST("Data/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);

      if (isAmbiguousMftTrack(reassociated2dMftTrack, true)) {
        registry.fill(HIST("Data/hReassociation2dMftTracks"), Reassociation2dMftTracks::NotReassociated2dMftTracks);
      }

    } // end of loop over reassociated MFT tracks

    // TO-DO the same for reassociated3d (change the histograms)
  }
  PROCESS_SWITCH(MftReassociationValidation, processData, "Process MFT reassociation validation for DATA", false);

  void processCreateLookupTable(FilteredCollisionsWSelMultMcLabels const& collisions, soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions)
  {
    recoVtxX.clear();
    recoVtxY.clear();
    recoVtxZ.clear();
    recoMcCollisionId.clear();
    recoMcCollBestCollisionIndex.clear();

    recoVtxX.reserve(collisions.size());
    recoVtxY.reserve(collisions.size());
    recoVtxZ.reserve(collisions.size());
    recoMcCollisionId.reserve(collisions.size());
    recoMcCollBestCollisionIndex.reserve(mcCollisions.size());

    for (auto const& col : collisions) {
      recoVtxX.emplace(col.globalIndex(), col.posX());
      recoVtxY.emplace(col.globalIndex(), col.posY());
      recoVtxZ.emplace(col.globalIndex(), col.posZ());
      recoMcCollisionId.emplace(col.globalIndex(), col.mcCollisionId());
    }

    for (auto const& mcCol : mcCollisions) {
      recoMcCollBestCollisionIndex.emplace(mcCol.globalIndex(), mcCol.bestCollisionIndex());
    }
  }
  PROCESS_SWITCH(MftReassociationValidation, processCreateLookupTable, "Process look uptable creation", false);

  void processMcReassociated2d(FilteredCollisionsWSelMultMcLabels::iterator const& collision,
                               FilteredMftTracksWCollsMcLabels const& /*mftTracks*/,
                               soa::SmallGroups<soa::Join<aod::BestCollisionsFwd, aod::McMFTTrackLabels, aod::MFTTrkCompColls>> const& reassociated2dMftTracks,
                               aod::McCollisions const& /*mcCollisions*/,
                               aod::McParticles const& /*particles*/,
                               aod::BCsWithTimestamps const& /*bcs*/)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    loadZVertexShiftCorrection(bc);

    registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::AllMonteCarloEvents);
    registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::AllEventsPrecise);

    if (!collision.has_mcCollision()) {
      registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::HasNotMonteCarloCollision);
      return;
    }

    registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::HasMonteCarloCollision);
    registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::HasMcCollision);

    const int mcCollisionId = collision.mcCollisionId();
    auto iteratorMcCollisionBestCollIndex = recoMcCollBestCollisionIndex.find(mcCollisionId);
    if (iteratorMcCollisionBestCollIndex == recoMcCollBestCollisionIndex.end()) {
      return;
    }
    const float mcCollisionBestCollIndex = iteratorMcCollisionBestCollIndex->second;

    if (configCollision.requireMcCollBestCollIndex && (collision.globalIndex() != mcCollisionBestCollIndex)) {
      return;
    }

    registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNotSplitVertex);

    if (!isAcceptedCollision(collision, true)) {
      return;
    }

    registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::MonteCarloEventsAfterEventSelection);

    for (auto const& reassociated2dMftTrack : reassociated2dMftTracks) {

      if (reassociated2dMftTrack.has_mcParticle()) {
        if (configTask.keepOnlyPhysicalPrimary && !reassociated2dMftTrack.mcParticle().isPhysicalPrimary()) {
          continue;
        }
        if (configTask.keepOnlySecondaries && reassociated2dMftTrack.mcParticle().isPhysicalPrimary()) {
          continue;
        }
      }

      registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
      registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::AllMonteCarloTracks);
      registry.fill(HIST("MC/hTrackAmbiguityCheck"), TrackAmbiguityCheckStep::AllTracksCheck);
      registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AllTracksPrecise);

      auto templatedTrack = reassociated2dMftTrack.template mfttrack_as<FilteredMftTracksWCollsMcLabels>();

      o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(templatedTrack, mZShift);
      std::array<double, 3> dcaInfOrig;
      trackPar.propagateToDCAhelix(bZ, {collision.posX(), collision.posY(), collision.posZ()}, dcaInfOrig);
      auto dcaXYoriginal = 999.f;
      dcaXYoriginal = std::sqrt(dcaInfOrig[0] * dcaInfOrig[0] + dcaInfOrig[1] * dcaInfOrig[1]);

      if (configMft.cutFakeTracks && templatedTrack.mcMask() != 0) {
        continue;
      }
      registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::IsTrueTrack);
      if (configMft.cutOrphanTracksExplicitly && reassociated2dMftTrack.ambDegree() == 0) {
        continue;
      }
      registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AfterOrphanCut);

      if (!isAcceptedMftTrack(templatedTrack, true, false, dcaXYoriginal, dcaInfOrig[2])) {
        continue;
      }

      if (reassociated2dMftTrack.ambDegree() == static_cast<int>(reassociated2dMftTrack.compatibleCollIds().size())) {
        registry.fill(HIST("MC/hTrackAmbiguityCheck"), TrackAmbiguityCheckStep::IsAmbDegreeEqualToCompatibleCollIdsSize);
      }

      registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);
      registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::MonteCarloTracksAfterTrackSelection);

      if (templatedTrack.has_mcParticle()) {
        registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::HasMonteCarloParticle);
        registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::HasMcParticle);

        auto particle = templatedTrack.template mcParticle_as<aod::McParticles>();
        float deltaX = -999.f;
        float deltaY = -999.f;
        float deltaZ = -999.f;
        float reassociatedDeltaX = -999.f;
        float reassociatedDeltaY = -999.f;
        float reassociatedDeltaZ = -999.f;
        auto collision = templatedTrack.collision_as<FilteredCollisionsWSelMultMcLabels>();
        auto xPosTrue = reassociated2dMftTrack.mcParticle().mcCollision().posX();
        auto yPosTrue = reassociated2dMftTrack.mcParticle().mcCollision().posY();
        auto zPosTrue = reassociated2dMftTrack.mcParticle().mcCollision().posZ();

        const int bestRecoColl = reassociated2dMftTrack.bestCollisionId();
        const int originalRecoColl = templatedTrack.collisionId();

        auto iteratorOriginalCollVtxX = recoVtxX.find(originalRecoColl);
        auto iteratorOriginalCollVtxY = recoVtxY.find(originalRecoColl);
        auto iteratorOriginalCollVtxZ = recoVtxZ.find(originalRecoColl);
        auto iteratorBestCollVtxX = recoVtxX.find(bestRecoColl);
        auto iteratorBestCollVtxY = recoVtxY.find(bestRecoColl);
        auto iteratorBestCollVtxZ = recoVtxZ.find(bestRecoColl);
        auto iteratorRecoMcCollisionId = recoMcCollisionId.find(bestRecoColl);

        if (iteratorOriginalCollVtxX == recoVtxX.end()) {
          // bestRecoColl not found in reco collisions map -> skip or count separately
          continue;
        }
        if (iteratorOriginalCollVtxY == recoVtxY.end()) {
          continue;
        }
        if (iteratorOriginalCollVtxZ == recoVtxZ.end()) {
          continue;
        }
        if (iteratorBestCollVtxX == recoVtxX.end()) {
          continue;
        }
        if (iteratorBestCollVtxY == recoVtxY.end()) {
          continue;
        }
        if (iteratorBestCollVtxZ == recoVtxZ.end()) {
          continue;
        }
        if (iteratorRecoMcCollisionId == recoMcCollisionId.end()) {
          continue;
        }

        const float xPosOriginalColl = iteratorOriginalCollVtxX->second;
        const float yPosOriginalColl = iteratorOriginalCollVtxY->second;
        const float zPosOriginalColl = iteratorOriginalCollVtxZ->second;
        const float xPosBestColl = iteratorBestCollVtxX->second;
        const float yPosBestColl = iteratorBestCollVtxY->second;
        const float zPosBestColl = iteratorBestCollVtxZ->second;
        const int mcCollisionIdReco = iteratorRecoMcCollisionId->second;

        deltaX = xPosOriginalColl - xPosTrue;
        deltaY = yPosOriginalColl - yPosTrue;
        deltaZ = zPosOriginalColl - zPosTrue;

        reassociatedDeltaX = xPosBestColl - xPosTrue;
        reassociatedDeltaY = yPosBestColl - yPosTrue;
        reassociatedDeltaZ = zPosBestColl - zPosTrue;

        const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
        const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
        const auto dcaZtruth(particle.vz() - particle.mcCollision().posZ());
        auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth);

        if (reassociated2dMftTrack.ambDegree() > 1) { // AMBIGUOUS TRACKS
          registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfAmbiguousTracks);
          registry.fill(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"), MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguous);
          registry.fill(HIST("MC/hReassociation2dMftTracks"), Reassociation2dMftTracks::AllAmbiguousTracksAfterTrackSelectionsFor2d);
          hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
          hDcaAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);

          if (collision.mcCollisionId() == particle.mcCollisionId()) {
            registry.fill(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"), MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguousAndMatchedToTrueCollision);
            hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          } else {
            registry.fill(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"), MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguousAndNotMatchedToTrueCollision);
            hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          }

          if (templatedTrack.collisionId() == reassociated2dMftTrack.bestCollisionId()) { // IS NOT 2D REASSOCIATED

            registry.fill(HIST("MC/hReassociation2dMftTracks"), Reassociation2dMftTracks::NotReassociated2dMftTracks);
            registry.fill(HIST("MC/hIsNot2dReassociatedAndMatchedToTrueCollision"), MftNot2dReassociatedAndMatchedToTrueCollisionStep::IsNot2dReassociated);
            hZVtxDiffNot2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
            hDcaNot2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated2dMftTrack.bestDCAXY(), dcaXYtruth);

            if (mcCollisionIdReco == particle.mcCollisionId()) {
              registry.fill(HIST("MC/hIsNot2dReassociatedAndMatchedToTrueCollision"), MftNot2dReassociatedAndMatchedToTrueCollisionStep::IsNot2dReassociatedAndMatchedToTrueCollision);
              hZVtxDiffNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDcaNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated2dMftTrack.bestDCAXY(), dcaXYtruth);
            } else {
              registry.fill(HIST("MC/hIsNot2dReassociatedAndMatchedToTrueCollision"), MftNot2dReassociatedAndMatchedToTrueCollisionStep::IsNot2dReassociatedAndNotMatchedToTrueCollision);
              hZVtxDiffNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDcaNot2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated2dMftTrack.bestDCAXY(), dcaXYtruth);

              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              if (isTrueCollisionAmongCompatibleCollisions(reassociated2dMftTrack)) {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
              } else {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
              }
            }

          } else { // IS 2D REASSOCIATED

            registry.fill(HIST("MC/hReassociation2dMftTracks"), Reassociation2dMftTracks::Reassociated2dMftTracks);
            registry.fill(HIST("MC/hIs2dReassociatedAndMatchedToTrueCollision"), Mft2dReassociatedAndMatchedToTrueCollisionStep::Is2dReassociated);
            hZVtxDiff2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
            hDca2dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated2dMftTrack.bestDCAXY(), dcaXYtruth);

            // is collision.mcCollisionId() the reassociated collision vertex ? or the initial collision
            if (mcCollisionIdReco == particle.mcCollisionId()) {
              registry.fill(HIST("MC/hIs2dReassociatedAndMatchedToTrueCollision"), Mft2dReassociatedAndMatchedToTrueCollisionStep::Is2dReassociatedAndMatchedToTrueCollision);
              hZVtxDiff2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDca2dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated2dMftTrack.bestDCAXY(), dcaXYtruth);
            } else {
              registry.fill(HIST("MC/hIs2dReassociatedAndMatchedToTrueCollision"), Mft2dReassociatedAndMatchedToTrueCollisionStep::Is2dReassociatedAndNotMatchedToTrueCollision);
              hZVtxDiff2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDca2dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated2dMftTrack.bestDCAXY(), dcaXYtruth);

              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              if (isTrueCollisionAmongCompatibleCollisions(reassociated2dMftTrack)) {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
              } else {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
              }
            }
          }

        } else { // NON AMBI TRACKS

          registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks);
          registry.fill(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"), MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguous);
          hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
          hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);

          if (collision.mcCollisionId() == particle.mcCollisionId()) {
            registry.fill(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"), MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguousAndMatchedToTrueCollision);
            hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          } else {
            registry.fill(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"), MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguousAndNotMatchedToTrueCollision);
            hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          }

        } // end of if non ambi
      } else {
        registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::HasNotMonteCarloParticle);
      }
    } // end of loop over reassociated2dMftTracks
  }
  PROCESS_SWITCH(MftReassociationValidation, processMcReassociated2d, "Process MFT reassociation2d validation for MONTE-CARLO", false);

  void processMcReassociated3d(FilteredCollisionsWSelMultMcLabels::iterator const& collision,
                               FilteredMftTracksWCollsMcLabels const& /*mftTracks*/,
                               soa::SmallGroups<soa::Join<aod::BestCollisionsFwd3d, aod::McMFTTrackLabels, aod::MFTTrkCompColls>> const& reassociated3dMftTracks,
                               aod::McCollisions const& /*mcCollisions*/,
                               aod::McParticles const& /*particles*/,
                               aod::BCsWithTimestamps const& /*bcs*/)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    loadZVertexShiftCorrection(bc);

    registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::AllMonteCarloEvents);
    registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::AllEventsPrecise);

    if (!collision.has_mcCollision()) {
      registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::HasNotMonteCarloCollision);
      return;
    }

    registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::HasMonteCarloCollision);
    registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::HasMcCollision);

    const int mcCollisionId = collision.mcCollisionId();
    auto iteratorMcCollisionBestCollIndex = recoMcCollBestCollisionIndex.find(mcCollisionId);
    if (iteratorMcCollisionBestCollIndex == recoMcCollBestCollisionIndex.end()) {
      return;
    }
    const float mcCollisionBestCollIndex = iteratorMcCollisionBestCollIndex->second;

    if (configCollision.requireMcCollBestCollIndex && (collision.globalIndex() != mcCollisionBestCollIndex)) {
      return;
    }

    registry.fill(HIST("hPreciseEventCounter"), SpecificEventSelectionStep::IsNotSplitVertex);

    if (!isAcceptedCollision(collision, true)) {
      return;
    }

    registry.fill(HIST("MC/hMonteCarloEventCounter"), MonteCarloEventSelectionStep::MonteCarloEventsAfterEventSelection);

    for (auto const& reassociated3dMftTrack : reassociated3dMftTracks) {

      if (reassociated3dMftTrack.has_mcParticle()) {
        if (configTask.keepOnlyPhysicalPrimary && !reassociated3dMftTrack.mcParticle().isPhysicalPrimary()) {
          continue;
        }
        if (configTask.keepOnlySecondaries && reassociated3dMftTrack.mcParticle().isPhysicalPrimary()) {
          continue;
        }
      }

      registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
      registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::AllMonteCarloTracks);
      registry.fill(HIST("MC/hTrackAmbiguityCheck"), TrackAmbiguityCheckStep::AllTracksCheck);
      registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AllTracksPrecise);

      auto templatedTrack = reassociated3dMftTrack.template mfttrack_as<FilteredMftTracksWCollsMcLabels>();

      o2::track::TrackParCovFwd trackPar = o2::aod::fwdtrackutils::getTrackParCovFwdShift(templatedTrack, mZShift);
      std::array<double, 3> dcaInfOrig;
      trackPar.propagateToDCAhelix(bZ, {collision.posX(), collision.posY(), collision.posZ()}, dcaInfOrig);
      auto dcaXYoriginal = 999.f;
      dcaXYoriginal = std::sqrt(dcaInfOrig[0] * dcaInfOrig[0] + dcaInfOrig[1] * dcaInfOrig[1]);

      if (configMft.cutFakeTracks && templatedTrack.mcMask() != 0) {
        continue;
      }
      registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::IsTrueTrack);
      if (configMft.cutOrphanTracksExplicitly && reassociated3dMftTrack.ambDegree() == 0) {
        continue;
      }
      registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::AfterOrphanCut);

      if (!isAcceptedMftTrack(templatedTrack, true, false, dcaXYoriginal, dcaInfOrig[2])) {
        continue;
      }

      if (reassociated3dMftTrack.ambDegree() == static_cast<int>(reassociated3dMftTrack.compatibleCollIds().size())) {
        registry.fill(HIST("MC/hTrackAmbiguityCheck"), TrackAmbiguityCheckStep::IsAmbDegreeEqualToCompatibleCollIdsSize);
      }

      registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);
      registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::MonteCarloTracksAfterTrackSelection);

      if (templatedTrack.has_mcParticle()) {
        registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::HasMonteCarloParticle);
        registry.fill(HIST("MC/hPreciseTrackSelectionCounter"), SpecificTrackSelectionStep::HasMcParticle);

        auto particle = templatedTrack.template mcParticle_as<aod::McParticles>();
        float deltaX = -999.f;
        float deltaY = -999.f;
        float deltaZ = -999.f;
        float reassociatedDeltaX = -999.f;
        float reassociatedDeltaY = -999.f;
        float reassociatedDeltaZ = -999.f;
        auto collision = templatedTrack.collision_as<FilteredCollisionsWSelMultMcLabels>();
        auto xPosTrue = reassociated3dMftTrack.mcParticle().mcCollision().posX();
        auto yPosTrue = reassociated3dMftTrack.mcParticle().mcCollision().posY();
        auto zPosTrue = reassociated3dMftTrack.mcParticle().mcCollision().posZ();

        const int bestRecoColl = reassociated3dMftTrack.bestCollisionId();
        const int originalRecoColl = templatedTrack.collisionId();

        auto iteratorOriginalCollVtxX = recoVtxX.find(originalRecoColl);
        auto iteratorOriginalCollVtxY = recoVtxY.find(originalRecoColl);
        auto iteratorOriginalCollVtxZ = recoVtxZ.find(originalRecoColl);
        auto iteratorBestCollVtxX = recoVtxX.find(bestRecoColl);
        auto iteratorBestCollVtxY = recoVtxY.find(bestRecoColl);
        auto iteratorBestCollVtxZ = recoVtxZ.find(bestRecoColl);
        auto iteratorRecoMcCollisionId = recoMcCollisionId.find(bestRecoColl);

        if (iteratorOriginalCollVtxX == recoVtxX.end()) {
          // bestRecoColl not found in reco collisions map -> skip or count separately
          continue;
        }
        if (iteratorOriginalCollVtxY == recoVtxY.end()) {
          continue;
        }
        if (iteratorOriginalCollVtxZ == recoVtxZ.end()) {
          continue;
        }
        if (iteratorBestCollVtxX == recoVtxX.end()) {
          continue;
        }
        if (iteratorBestCollVtxY == recoVtxY.end()) {
          continue;
        }
        if (iteratorBestCollVtxZ == recoVtxZ.end()) {
          continue;
        }
        if (iteratorRecoMcCollisionId == recoMcCollisionId.end()) {
          continue;
        }

        const float xPosOriginalColl = iteratorOriginalCollVtxX->second;
        const float yPosOriginalColl = iteratorOriginalCollVtxY->second;
        const float zPosOriginalColl = iteratorOriginalCollVtxZ->second;
        const float xPosBestColl = iteratorBestCollVtxX->second;
        const float yPosBestColl = iteratorBestCollVtxY->second;
        const float zPosBestColl = iteratorBestCollVtxZ->second;
        const int mcCollisionIdReco = iteratorRecoMcCollisionId->second;

        deltaX = xPosOriginalColl - xPosTrue;
        deltaY = yPosOriginalColl - yPosTrue;
        deltaZ = zPosOriginalColl - zPosTrue;

        reassociatedDeltaX = xPosBestColl - xPosTrue;
        reassociatedDeltaY = yPosBestColl - yPosTrue;
        reassociatedDeltaZ = zPosBestColl - zPosTrue;

        const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
        const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
        const auto dcaZtruth(particle.vz() - particle.mcCollision().posZ());
        auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth);

        if (reassociated3dMftTrack.ambDegree() > 1) { // AMBIGUOUS TRACKS
          registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfAmbiguousTracks);
          registry.fill(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"), MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguous);
          registry.fill(HIST("MC/hReassociation3dMftTracks"), Reassociation3dMftTracks::AllAmbiguousTracksAfterTrackSelectionsFor3d);
          hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
          hDcaAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);

          if (collision.mcCollisionId() == particle.mcCollisionId()) {
            registry.fill(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"), MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguousAndMatchedToTrueCollision);
            hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          } else {
            registry.fill(HIST("MC/hIsAmbiguousTrackMatchedToTrueCollision"), MftAmbiguousAndMatchedToTrueCollisionStep::IsAmbiguousAndNotMatchedToTrueCollision);
            hZVtxDiffAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          }

          if (templatedTrack.collisionId() == reassociated3dMftTrack.bestCollisionId()) { // IS NOT 3D REASSOCIATED

            registry.fill(HIST("MC/hReassociation3dMftTracks"), Reassociation3dMftTracks::NotReassociated3dMftTracks);
            registry.fill(HIST("MC/hIsNot3dReassociatedAndMatchedToTrueCollision"), MftNot3dReassociatedAndMatchedToTrueCollisionStep::IsNot3dReassociated);
            hZVtxDiffNot3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
            hDcaNot3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated3dMftTrack.bestDCAXY(), reassociated3dMftTrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

            if (mcCollisionIdReco == particle.mcCollisionId()) {
              registry.fill(HIST("MC/hIsNot3dReassociatedAndMatchedToTrueCollision"), MftNot3dReassociatedAndMatchedToTrueCollisionStep::IsNot3dReassociatedAndMatchedToTrueCollision);
              hZVtxDiffNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDcaNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated3dMftTrack.bestDCAXY(), reassociated3dMftTrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
            } else {
              registry.fill(HIST("MC/hIsNot3dReassociatedAndMatchedToTrueCollision"), MftNot3dReassociatedAndMatchedToTrueCollisionStep::IsNot3dReassociatedAndNotMatchedToTrueCollision);
              hZVtxDiffNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDcaNot3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated3dMftTrack.bestDCAXY(), reassociated3dMftTrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              if (isTrueCollisionAmongCompatibleCollisions(reassociated3dMftTrack)) {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
              } else {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaNotReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
              }
            }

          } else { // IS 3D REASSOCIATED

            registry.fill(HIST("MC/hReassociation3dMftTracks"), Reassociation3dMftTracks::Reassociated3dMftTracks);
            registry.fill(HIST("MC/hIs3dReassociatedAndMatchedToTrueCollision"), Mft3dReassociatedAndMatchedToTrueCollisionStep::Is3dReassociated);
            hZVtxDiff3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
            hDca3dReassociatedTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated3dMftTrack.bestDCAXY(), reassociated3dMftTrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

            if (mcCollisionIdReco == particle.mcCollisionId()) {
              registry.fill(HIST("MC/hIs3dReassociatedAndMatchedToTrueCollision"), Mft3dReassociatedAndMatchedToTrueCollisionStep::Is3dReassociatedAndMatchedToTrueCollision);
              hZVtxDiff3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDca3dReassociatedTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated3dMftTrack.bestDCAXY(), reassociated3dMftTrack.bestDCAZ(), dcaXYtruth, dcaZtruth);
            } else {
              registry.fill(HIST("MC/hIs3dReassociatedAndMatchedToTrueCollision"), Mft3dReassociatedAndMatchedToTrueCollisionStep::Is3dReassociatedAndNotMatchedToTrueCollision);
              hZVtxDiff3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociatedDeltaX, reassociatedDeltaY, reassociatedDeltaZ);
              hDca3dReassociatedTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), reassociated3dMftTrack.bestDCAXY(), reassociated3dMftTrack.bestDCAZ(), dcaXYtruth, dcaZtruth);

              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::AllWronglyAssociatedTracks);
              if (isTrueCollisionAmongCompatibleCollisions(reassociated3dMftTrack)) {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsTrueCollisionAmongCompatibleCollisions);
              } else {
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisions"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
                registry.fill(HIST("MC/hIsTrueCollisionAmongCompatibleCollisionsDcaReassociated"), MftIsTrueCollisionAmongCompatibleCollisionsStep::IsNotTrueCollisionAmongCompatibleCollisions);
              }
            }
          }

        } else { // NON AMBI TRACKS

          registry.fill(HIST("MC/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks);
          registry.fill(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"), MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguous);
          hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
          hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::AllTracks]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);

          if (collision.mcCollisionId() == particle.mcCollisionId()) {
            registry.fill(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"), MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguousAndMatchedToTrueCollision);
            hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::IsMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          } else {
            registry.fill(HIST("MC/hIsNonAmbiguousTrackMatchedToTrueCollision"), MftNonAmbiguousAndMatchedToTrueCollisionStep::IsNonAmbiguousAndNotMatchedToTrueCollision);
            hZVtxDiffNonAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), deltaX, deltaY, deltaZ);
            hDcaNonAmbiguousTracks[MatchedToTrueCollisionStep::IsNotMatchedToTrueCollision]->Fill(templatedTrack.pt(), templatedTrack.eta(), dcaXYoriginal, dcaInfOrig[2], dcaXYtruth, dcaZtruth);
          }

        } // end of if non ambi
      } else {
        registry.fill(HIST("MC/hMonteCarloTrackCounter"), MonteCarloTrackSelectionStep::HasNotMonteCarloParticle);
      }
    } // end of loop over reassociated3dMftTracks
  }
  PROCESS_SWITCH(MftReassociationValidation, processMcReassociated3d, "Process MFT reassociation3d validation for MONTE-CARLO", false);

}; // End of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MftReassociationValidation>(cfgc)};
}
