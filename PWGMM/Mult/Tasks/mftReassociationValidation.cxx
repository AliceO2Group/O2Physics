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
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include <THn.h>
#include <TPDGCode.h>
#include <TString.h>

#include <sys/types.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::aod::track;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum ReassociationMftTracks {
  NotReassociatedMftTracks = 0,
  ReassociatedMftTracks,
  NReassociationMftTracksSteps
};

enum EventSelectionStep {
  AllEvents = 0,
  AfterEventSelection,
  NEventSelectionSteps
};

enum MftTrackAmbiguityStep {
  AllMftTracks = 0,
  AfterTrackSelection,
  NumberOfAmbiguousTracks,
  NumberOfNonAmbiguousTracks,
  NMftAmbiguitySteps
};

enum MftTrackSelectionStep {
  NoSelection = 0,
  Eta,
  Cluster,
  Pt,
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

static constexpr TrackSelectionFlags::flagtype TrackSelectionIts =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;
static constexpr TrackSelectionFlags::flagtype TrackSelectionTpc =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDca =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDcaxyOnly =
  TrackSelectionFlags::kDCAxy;

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
    Configurable<bool> doHeavyFlavor{"doHeavyFlavor", false, "Flag to know we in the heavy flavor case or not"};
    Configurable<bool> doReferenceFlow{"doReferenceFlow", false, "Flag to know if reference flow should be done"};
    Configurable<bool> isReadoutCenter{"isReadoutCenter", false, "Enable Readout Center"};
    Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per event"};
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
    Configurable<float> zVertexMax{"zVertexMax", 10.0f, "Accepted z-vertex range"};
  } configCollision;

  //  configurables for central barrel tracks
  struct : ConfigurableGroup {
    std::string prefix = "ConfigCentral_group";
    Configurable<float> dcaZCentralTrackMax{"dcaZCentralTrackMax", 0.2f, "max dcaZ of central tracks"};
    Configurable<float> etaCentralTrackMax{"etaCentralTrackMax", 0.8f, "max. eta of central tracks"};
    Configurable<bool> isApplyConversionCut{"isApplyConversionCut", false, "apply pair conversion cuts"};
    Configurable<bool> isApplyTwoTrackCut{"isApplyTwoTrackCut", false, "apply two track cut"};
    Configurable<bool> isApplyIndexOrdering{"isApplyIndexOrdering", false, "apply track1.index() <= track2.index() cut"};
    Configurable<bool> isApplyPtOrderingSameEvent{"isApplyPtOrderingSameEvent", false, "apply track1.pt() <= track2.pt() cut"};
    Configurable<bool> isApplyPtOrderingMixedEvent{"isApplyPtOrderingMixedEvent", false, "apply track1.pt() <= track2.pt() cut"};
    Configurable<bool> isApplySameTrackCut{"isApplySameTrackCut", false, "apply track1 == track2 cut"};
    Configurable<float> maxChi2ItsClusters{"maxChi2ItsClusters", 36.f, "max chi2 per ITS clusters"};
    Configurable<float> maxChi2TpcClusters{"maxChi2TpcClusters", 2.5f, "max chi2 per TPC clusters"};
    Configurable<float> maxMergingRadius{"maxMergingRadius", 2.5, "max radius for merging cut"};
    Configurable<float> mergingCut{"mergingCut", 0.02, "merging cut on track merge"};
    Configurable<float> minItsClusters{"minItsClusters", 5.0f, "cut for minimum ITS clusters"};
    Configurable<float> minMergingRadius{"minMergingRadius", 0.8, "max radius for merging cut"};
    Configurable<float> minTpcClusters{"minTpcClusters", 50.0f, "cut for minimum TPC clusters"};
    Configurable<float> minTpcCrossedRows{"minTpcCrossedRows", 70.0f, "cut for minimum TOC crossed rows"};
    Configurable<float> ptCentralTrackMin{"ptCentralTrackMin", 0.2f, "min. pT of central tracks"};
    Configurable<float> ptCentralTrackMax{"ptCentralTrackMax", 10.0f, "max. pT of central tracks"};
    Configurable<int> trackSelectionType{"trackSelectionType", 1, "Track selection: 0 -> kGlobalTrack or isGlobalTrackSDD , 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> No globalTrack selection"};
  } configCentral;

  //   configurables for MFT tracks
  struct : ConfigurableGroup {
    std::string prefix = "ConfigMft_group";
    Configurable<int> cutBestCollisionId{"cutBestCollisionId", 0, "cut on the best collision Id used in a filter"};
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
  } configMft;

  TF1* fPtDepDCAxy = nullptr;

  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // =========================
  //      using declarations : DATA
  // =========================

  using FilteredCollisionsWSelMult = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using FilteredMftTracks = soa::Filtered<aod::MFTTracks>;
  using FilteredMftTracksWColls = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>>;
  using FilteredMftTracksWCollsMcLabels = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>>;
  using MftReasso2dTracksWCollsMcLabels = soa::Join<aod::MFTTracks, aod::BestCollisionsFwd, aod::McMFTTrackLabels>;
  using MftReasso3dTracksWCollsMcLabels = soa::Join<aod::MFTTracks, aod::BestCollisionsFwd3d, aod::McMFTTrackLabels>;

  // =========================
  //      Filters & partitions : DATA
  // =========================

  //  Collision filters
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < configCollision.zVertexMax;

  // Central tracks filter
  Filter trackFilter = (nabs(aod::track::eta) < configCentral.etaCentralTrackMax) &&
                       (aod::track::pt > configCentral.ptCentralTrackMin) &&
                       (aod::track::pt < configCentral.ptCentralTrackMax) &&
                       requireGlobalTrackInFilter();

  Filter centralTrackItsTpcMatchingFilter = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC), ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true) &&
                                            ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                            ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);

  Filter centralTrackDcaFilter = (ifnode(configCentral.dcaZCentralTrackMax.node() > 0.f, nabs(aod::track::dcaZ) <= configCentral.dcaZCentralTrackMax && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly), ncheckbit(aod::track::trackCutFlag, TrackSelectionDca)));

  Filter centralTrackChi2TpcClusterFilter = (aod::track::tpcChi2NCl < configCentral.maxChi2TpcClusters);

  Filter centralTrackChi2ItsClusterFilter = (aod::track::itsChi2NCl < configCentral.maxChi2ItsClusters);

  Filter mftTrackEtaFilter = ((aod::fwdtrack::eta < configMft.etaMftTrackMaxFilter) && (aod::fwdtrack::eta > configMft.etaMftTrackMinFilter));

  // Filters below will be used for uncertainties
  Filter mftTrackCollisionIdFilter = (aod::fwdtrack::bestCollisionId >= 0);
  Filter mftTrackDcaXYFilter = (nabs(aod::fwdtrack::bestDCAXY) < configMft.mftMaxDCAxy);
  // Filter mftTrackDcaZFilter = (nabs(aod::fwdtrack::bestDCAZ) < configMft.mftMaxDCAz);

  // =========================
  //      Preslice : DATA
  // =========================

  Preslice<FilteredMftTracks> perColMftTracks = o2::aod::fwdtrack::collisionId;

  struct : ConfigurableGroup {
    std::string prefix = "ConfigAxis_group";
    ConfigurableAxis axisMass{"axisMass", {1, 1.5848, 2.1848}, "axis of invariant mass of candidates"};
    ConfigurableAxis binsMixingMultiplicity{"binsMixingMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity bins for event mixing"};
    ConfigurableAxis binsMixingVertex{"binsMixingVertex", {20, -10, 10}, "vertex bins for event mixing"};
    ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {1, -1.0, 1.0}, "eta axis for efficiency histograms"};
    ConfigurableAxis axisEtaAssociated{"axisEtaAssociated", {48, -4, -2}, "eta axis for MFT histograms"};
    ConfigurableAxis axisEtaTrigger{"axisEtaTrigger", {48, -1, 1}, "eta axis for TPC histograms"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity axis for histograms"};
    ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "phi axis for histograms"};
    ConfigurableAxis axisPt{"axisPt", {72, 0, 36}, "pt axis for histograms"};
    ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt associated axis for histograms"};
    ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.2, 10}, "pt axis for efficiency histograms"};
    ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
    ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
    ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {1, -10, 10}, "vertex axis for efficiency histograms"};
  } configAxis;

  HistogramRegistry registry{"registry"};

  // template <DataType DataType, CorrelationCase CorrelationCase, CorrelatedParticles CorrelatedParticles>
  // void addHistograms()
  // {
  //   registry.add(Form("%s%s%shEtaTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {configAxis.axisEtaTrigger}});
  //   registry.add(Form("%s%s%shPhiTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {configAxis.axisPhi}});
  //   registry.add(Form("%s%s%shPtTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {configAxis.axisPt}});
  //   registry.add(Form("%s%s%shYieldsTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH3F, {configAxis.axisMultiplicity, configAxis.axisPt, configAxis.axisEtaTrigger}});
  //   registry.add(Form("%s%s%shEtaPhiTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH3F, {configAxis.axisMultiplicity, configAxis.axisEtaTrigger, configAxis.axisPhi}});
  //   registry.add(Form("%s%s%shEtaAssociated", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {configAxis.axisEtaAssociated}});
  //   registry.add(Form("%s%s%shPhiAssociated", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {configAxis.axisPhi}});
  //   registry.add(Form("%s%s%shEtaPhiAssociated", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH3F, {configAxis.axisMultiplicity, configAxis.axisEtaAssociated, configAxis.axisPhi}});
  // }

  void addMftHistograms()
  {
    registry.add("Data/Mft/hAmbiguityOfMftTracks", "hAmbiguityOfMftTracks", {HistType::kTH1D, {{MftTrackAmbiguityStep::NMftAmbiguitySteps, -0.5, +MftTrackAmbiguityStep::NMftAmbiguitySteps - 0.5}}});
    std::string labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NMftAmbiguitySteps];
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::AllMftTracks] = "all MFT tracks";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::AfterTrackSelection] = "MFT tracks after selection";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NumberOfAmbiguousTracks] = "how much tracks are ambigous";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks] = "how much tracks are non-ambiguous";
    registry.get<TH1>(HIST("Data/Mft/hAmbiguityOfMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftTrackAmbiguityStep::NMftAmbiguitySteps; iBin++) {
      registry.get<TH1>(HIST("Data/Mft/hAmbiguityOfMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsAmbiguityOfMftTracks[iBin].data());
    }

    registry.add("Data/Mft/hMftTracksSelection", "hMftTracksSelection", {HistType::kTH1D, {{MftTrackSelectionStep::NMftTrackSelectionSteps, -0.5, +MftTrackSelectionStep::NMftTrackSelectionSteps - 0.5}}});
    std::string labelsMftTracksSelection[MftTrackSelectionStep::NMftTrackSelectionSteps];
    labelsMftTracksSelection[MftTrackSelectionStep::NoSelection] = "all MFT tracks";
    labelsMftTracksSelection[MftTrackSelectionStep::Eta] = "MFT tracks after eta selection";
    labelsMftTracksSelection[MftTrackSelectionStep::Cluster] = "MFT tracks after clusters selection";
    labelsMftTracksSelection[MftTrackSelectionStep::Pt] = "MFT tracks after pT selection";
    registry.get<TH1>(HIST("Data/Mft/hMftTracksSelection"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftTrackSelectionStep::NMftTrackSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("Data/Mft/hMftTracksSelection"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftTracksSelection[iBin].data());
    }

    registry.add("Data/Mft/hReassociationMftTracks", "hReassociationMftTracks", {HistType::kTH1D, {{ReassociationMftTracks::NReassociationMftTracksSteps, -0.5, +ReassociationMftTracks::NReassociationMftTracksSteps - 0.5}}});
    std::string labelsReassociationMftTracks[ReassociationMftTracks::NReassociationMftTracksSteps];
    labelsReassociationMftTracks[ReassociationMftTracks::NotReassociatedMftTracks] = "Ambiguous MFT tracks after track selection";
    labelsReassociationMftTracks[ReassociationMftTracks::ReassociatedMftTracks] = "Reassociated MFT tracks by DCAxy method";
    registry.get<TH1>(HIST("Data/Mft/hReassociationMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < ReassociationMftTracks::NReassociationMftTracksSteps; iBin++) {
      registry.get<TH1>(HIST("Data/Mft/hReassociationMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsReassociationMftTracks[iBin].data());
    }

    registry.add("Data/Mft/hPtMft", "", {HistType::kTH1D, {configAxis.axisPt}});
    registry.add("Data/Mft/hNMftTracks", "", {HistType::kTH1F, {configAxis.axisMultiplicity}});
    registry.add("Data/Mft/hNBestCollisionFwd", "", {HistType::kTH1F, {configAxis.axisMultiplicity}});
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

    //  =========================
    //      Event histograms
    //  =========================

    registry.add("Data/hVtxZ", "v_{z} (cm)", {HistType::kTH1D, {configAxis.axisVertex}});
    registry.add("Data/hNTracks", "", {HistType::kTH1F, {configAxis.axisMultiplicity}});
    registry.add(Form("Data/hMultiplicity_%s", WhatMultiplicityEstimator[configCollision.multiplicityEstimator].data()), "", {HistType::kTH1D, {configAxis.axisMultiplicity}});

    registry.add("Data/hEventCounter", "hEventCounter", {HistType::kTH1D, {{EventSelectionStep::NEventSelectionSteps, -0.5, +EventSelectionStep::NEventSelectionSteps - 0.5}}});
    std::string labels[EventSelectionStep::NEventSelectionSteps];
    labels[EventSelectionStep::AllEvents] = "all";
    labels[EventSelectionStep::AfterEventSelection] = "after Physics selection";
    registry.get<TH1>(HIST("Data/hEventCounter"))->SetMinimum(0);

    for (int iBin = 0; iBin < EventSelectionStep::NEventSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("Data/hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

  } // End of init() function

  // =========================
  //      Quality assessment functions
  // =========================

  // template <DataType DataType, CorrelationCase CorrelationCase, CorrelatedParticles CorrelatedParticles>
  // void fillTriggerQa(float multiplicity, float const& eta, float const& phi, float const& pt)
  // {
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hPtTrigger"), pt);
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaTrigger"), eta);
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hPhiTrigger"), phi);
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hYieldsTrigger"), multiplicity, pt, eta);
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaPhiTrigger"), multiplicity, eta, phi);
  // }

  // template <DataType DataType, CorrelationCase CorrelationCase, CorrelatedParticles CorrelatedParticles>
  // void fillAssociatedQa(float multiplicity, float const& eta, float const& phi)
  // {
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaAssociated"), eta);
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hPhiAssociated"), phi);
  //   registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaPhiAssociated"), multiplicity, eta, phi);
  // }

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

  // =========================
  //      Cuts with functions
  // =========================

  //  FIXME: Some collisions are rejected here, what causes (part of) differences with the D0 task
  template <typename TCollision>
  bool isAcceptedCollision(TCollision const& collision, bool fillHistograms = false)
  {
    if (fillHistograms) {
      registry.fill(HIST("Data/hEventCounter"), EventSelectionStep::AllEvents);
    }

    if (!collision.sel8()) {
      return false;
    }
    if (configCollision.isApplySameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (configCollision.isApplyGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (configCollision.isApplyNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (configCollision.isApplyGoodItsLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/hEventCounter"), EventSelectionStep::AfterEventSelection);
    }

    registry.fill(HIST("Data/hVtxZ"), collision.posZ());

    return true;
  }

  template <typename TTrack>
  bool isAcceptedCentralTrack(TTrack const& track)
  {
    if (track.tpcNClsFound() < configCentral.minTpcClusters) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < configCentral.minTpcCrossedRows) {
      return false;
    }
    if (track.itsNCls() < configCentral.minItsClusters) {
      return false;
    }
    return true;
  }

  //  TODO: Check how to put this into a Filter
  // I tried to put it as a filter, but filters for normal TPC tracks also apply to MFT tracks I think
  // and it seems that they are not compatible
  template <typename TTrack>
  bool isAcceptedMftTrack(TTrack const& mftTrack, bool fillHistograms)
  {
    // cut on the eta of MFT tracks
    if (mftTrack.eta() > configMft.etaMftTrackMax || mftTrack.eta() < configMft.etaMftTrackMin) {
      return false;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/Mft/hMftTracksSelection"), MftTrackSelectionStep::Eta);
    }

    // cut on the number of clusters of the reconstructed MFT track
    if (mftTrack.nClusters() < configMft.nClustersMftTrack) {
      return false;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/Mft/hMftTracksSelection"), MftTrackSelectionStep::Cluster);
    }

    // cut on the pT of MFT tracks (for test purposes)
    if (configMft.useMftPtCut && (mftTrack.pt() > configMft.ptMftTrackMax || mftTrack.pt() < configMft.ptMftTrackMin)) {
      return false;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/Mft/hMftTracksSelection"), MftTrackSelectionStep::Pt);
    }

    return true;
  }

  // Cut on ambiguous MFT tracks
  template <typename TTrack>
  bool isAmbiguousMftTrack(TTrack const& mftTrack, bool fillHistograms)
  {
    if (mftTrack.ambDegree() > 1) {
      if (fillHistograms) {
        registry.fill(HIST("Data/Mft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfAmbiguousTracks);
      }
      return true;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/Mft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks);
    }
    return false;
  }

  //============================================================================================
  //
  // PROCESS FUNCTIONS
  //
  //============================================================================================

  void processData(FilteredCollisionsWSelMult::iterator const& collision,
                   soa::SmallGroups<aod::BestCollisionsFwd> const& reassociated2dMftTracks,
                   aod::BCsWithTimestamps const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    // const auto multiplicity = getMultiplicityEstimator(collision, true);

    for (const auto& reassociated2dMftTrack : reassociated2dMftTracks) {

      registry.fill(HIST("Data/Mft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
      auto templatedMftTrack = reassociated2dMftTrack.template mfttrack_as<FilteredMftTracks>();

      if (!isAcceptedMftTrack(templatedMftTrack, false)) {
        continue;
      }

      registry.fill(HIST("Data/Mft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);

      if (isAmbiguousMftTrack(reassociated2dMftTrack, true)) {
        registry.fill(HIST("Data/Mft/hReassociationMftTracks"), ReassociationMftTracks::NotReassociatedMftTracks);
      }

    } // end of loop over reassociated MFT tracks
  }
  PROCESS_SWITCH(MftReassociationValidation, processData, "Process MFT reassociation validation for DATA", false);

  // void processMc(FilteredCollisionsWSelMult::iterator const& collision,
  //                FilteredTracksWDcaSel const& tracks,
  //                aod::BCsWithTimestamps const&)
  // {

  // }
  // PROCESS_SWITCH(MftReassociationValidation, processMc, "Process MFT reassociation validation for MONTE-CARLO", false);

}; // End of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MftReassociationValidation>(cfgc)};
}
