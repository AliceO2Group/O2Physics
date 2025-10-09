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

/// \file taskFlow.cxx
/// \brief HF-h correlations in TPC-TPC and TPC-MFT
/// \author Alexian Lejeune <alexian.lejeune@cern.ch >, Czech Technical University in Prague
/// \author Katarina Krizkova Gajdosova <katarina.gajdosova@cern.ch>, CERN
/// \author Maja Kabus <maja.kabus@cern.ch>, CERN

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsPid.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DetectorsCommonDataFormats/AlignParam.h>
#include <FT0Base/Geometry.h>
#include <FV0Base/Geometry.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include <THn.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
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
using namespace o2::analysis;
using namespace o2::aod::pid_tpc_tof_utils;
using namespace o2::aod::track;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum MftTrackSelectionStep {
  NoSelection = 0,
  Eta,
  Cluster,
  NMftTrackSelectionSteps
};

enum MftTrackAmbiguityStep {
  AllMftTracks = 0,
  AfterTrackSelection,
  NumberOfAmbiguousTracks,
  NumberOfNonAmbiguousTracks,
  NMftAmbiguitySteps
};

enum MultiplicityEstimators {
  MultNTracksPV = 0,
  MultNumContrib,
  MultFT0C,
  MultFT0M
};

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

enum DataType {
  Data,
  Mc
};

// enum EventType {
//   SameEvent,
//   MixedEvent
// };

enum CorrelationCase {
  TpcTpc,
  TpcMft,
  TpcFv0a,
  MftFv0a,
  TpcFt0a,
  MftFt0a
};

enum CorrelatedParticles {
  ChPartChPart,
  D0ChPart,
  LcChPart
};

// static constexpr std::string_view whatEventType[] = {"SameEvent/", "MixedEvent/"};
static constexpr std::string_view WhatDataType[] = {"Data/", "MC/"};
static constexpr std::string_view WhatCorrelationCase[] = {"TpcTpc/", "TpcMft/", "TpcFv0a/", "MftFv0a/", "TpcFt0a/", "MftFt0a/"};
static constexpr std::string_view WhatParticles[] = {"ChPartChPart/", "D0ChPart/", "LcChPart/"};
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

// static constexpr float kPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};

struct HfTaskFlow {

  struct : ConfigurableGroup {
    std::string prefix = "ConfigCcdb_group";
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } configCcdb;

  //  configurables for processing options

  struct : ConfigurableGroup {
    std::string prefix = "ConfigTask_group";
    Configurable<bool> centralityBinsForMc{"centralityBinsForMc", false, "false = OFF, true = ON for data like multiplicity/centrality bins for MC steps"};
    Configurable<bool> doHeavyFlavor{"doHeavyFlavor", false, "Flag to know we in the heavy flavor case or not"};
    Configurable<bool> doReferenceFlow{"doReferenceFlow", false, "Flag to know if reference flow should be done"};
    Configurable<bool> isReadoutCenter{"isReadoutCenter", false, "Enable Readout Center"};
    // Configurable<float> doTwoTrackCut{"doTwoTrackCut", -1, "Two track cut: -1 = off; >0 otherwise distance value (suggested: 0.02)"};
    Configurable<bool> processMc{"processMc", false, "Flag to run on MC"};
    Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per event"};
    // Configurable<float> twoTrackCutMinRadius{"twoTrackCutMinRadius", 0.8f, "Two track cut : radius in m from which two tracks cuts are applied"};
  } configTask;

  //   configurables for collisions
  struct : ConfigurableGroup {
    std::string prefix = "ConfigCollision_group";
    Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", false, "Enable GoodZvtxFT0vsPV cut"};
    Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", false, "Enable SameBunchPileup cut"};
    Configurable<int> multiplicityEstimator{"multiplicityEstimator", 0, "0: multNTracksPV, 1: numContrib, 2: multFT0C, 3: multFT0M, 4: centFT0C, 5: centFT0CVariants1s, 6: centFT0M, 7: centFV0A, 8: centNTracksPV, 9: centNGlobal, 10: centMFT"};
    Configurable<bool> isApplyNoCollInTimeRangeStrict{"isApplyNoCollInTimeRangeStrict", false, ""};
    Configurable<float> zVertexMax{"zVertexMax", 7.0f, "Accepted z-vertex range"};
  } configCollision;

  //  configurables for central barrel tracks
  struct : ConfigurableGroup {
    std::string prefix = "ConfigCentralTracks_group";
    Configurable<float> etaCentralTrackMax{"etaCentralTrackMax", 0.8f, "max. eta of central tracks"};
    Configurable<float> ptCentralTrackMin{"ptCentralTrackMin", 0.2f, "min. pT of central tracks"};
    Configurable<float> ptCentralTrackMax{"ptCentralTrackMax", 10.0f, "max. pT of central tracks"};
    Configurable<float> dcaZCentralTrackMax{"dcaZCentralTrackMax", 0.2f, "max dcaZ of central tracks"};
  } configCentral;

  //  configurables for HF candidates
  struct : ConfigurableGroup {
    std::string prefix = "ConfigCandidates_group";
    Configurable<float> etaCandidateMax{"etaCandidateMax", 0.8f, "max. eta of HF candidate"};
    Configurable<std::vector<int>> mcTriggerPdgs{"mcTriggerPdgs", {421, -421}, "MC PDG codes to use exclusively as trigger particles. D0= +-421, Lc = +-4122"};
    Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for Hf candidates"};
    Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
    Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  } configCandidates;

  //   configurables for MFT tracks
  struct : ConfigurableGroup {
    std::string prefix = "ConfigMft_group";
    Configurable<float> etaMftTrackMax{"etaMftTrackMax", -2.4f, "Maximum value for the eta of MFT tracks"};
    Configurable<float> etaMftTrackMin{"etaMftTrackMin", -3.36f, "Minimum value for the eta of MFT tracks"};
    Configurable<float> etaMftTrackMaxFilter{"etaMftTrackMaxFilter", -2.0f, "Maximum value for the eta of MFT tracks"};
    Configurable<float> etaMftTrackMinFilter{"etaMftTrackMinFilter", -3.9f, "Minimum value for the eta of MFT tracks"};
    Configurable<float> mftMaxDCAxy{"mftMaxDCAxy", 2.0f, "Cut on dcaXY for MFT tracks"};
    Configurable<float> mftMaxDCAz{"mftMaxDCAz", 2.0f, "Cut on dcaZ for MFT tracks"};
    Configurable<int> nClustersMftTrack{"nClustersMftTrack", 5, "Minimum number of clusters for the reconstruction of MFT tracks"};
  } configMft;

  HfHelper hfHelper;
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::vector<o2::detectors::AlignParam>* offsetFT0{};
  std::vector<o2::detectors::AlignParam>* offsetFV0{};
  o2::ccdb::CcdbApi ccdbApi;
  o2::ft0::Geometry ft0Det;
  o2::fv0::Geometry* fv0Det{};
  std::vector<int> hfIndexCache;

  // =========================
  //      using declarations : DATA
  // =========================

  using FilteredCollisionsWSelMult = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using HfCandidatesSelD0 = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using HfCandidatesSelLc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using FilteredTracksWDcaSel = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;

  using FilteredMftTracks = soa::Filtered<aod::MFTTracks>;
  //  using FilteredMftTracksWColls = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>>;
  // using FilteredAndReassociatedMftTracks = soa::Filtered<soa::Join<aod::BestCollisionsFwd, aod::MFTTracks>>;

  // =========================
  //      using declarations : MONTE CARLO
  // =========================

  using FilteredCollisionsWSelMultMcLabels = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults>>;
  using FilteredMcCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::MultMCExtras, aod::McCollsExtra>>;
  using HfCandidatesSelD0McRec = soa::Join<HfCandidatesSelD0, aod::HfCand2ProngMcRec>;
  using HfCandidatesSelLcMcRec = soa::Join<HfCandidatesSelLc, aod::HfCand3ProngMcRec>;
  using McParticles = aod::McParticles;
  using McParticles2ProngMatched = soa::Join<McParticles, aod::HfCand2ProngMcGen>;
  using McParticles3ProngMatched = soa::Join<McParticles, aod::HfCand3ProngMcGen>;
  using MftTracksMcLabels = soa::Join<FilteredMftTracks, aod::McMFTTrackLabels>;
  using FilteredTracksWDcaSelMC = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::McTrackLabels>>;
  // using FilteredMftTracksWCollsMcLabels = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>>;

  // =========================
  //      Filters & partitions : DATA
  // =========================

  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilterD0 = (aod::hf_sel_candidate_d0::isSelD0 >= configCandidates.selectionFlagHf) ||
                             (aod::hf_sel_candidate_d0::isSelD0bar >= configCandidates.selectionFlagHf);

  Filter candidateFilterLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= configCandidates.selectionFlagHf) ||
                             (aod::hf_sel_candidate_lc::isSelLcToPiKP >= configCandidates.selectionFlagHf);

  //  Collision filters
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < configCollision.zVertexMax;

  // Central tracks filter
  Filter centralTrackEtaPtFilter = (nabs(aod::track::eta) < configCentral.etaCentralTrackMax) && (aod::track::pt > configCentral.ptCentralTrackMin);
  Filter centralTrackTpcFilter = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                        ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true);
  Filter centralTrackItsFilter = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                 ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);
  Filter centralTrackDcaFilter = ifnode(configCentral.dcaZCentralTrackMax.node() > 0.f, nabs(aod::track::dcaZ) <= configCentral.dcaZCentralTrackMax && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                        ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));

  Filter mftTrackFilter = (aod::fwdtrack::eta < configMft.etaMftTrackMaxFilter) &&
                          (aod::fwdtrack::eta > configMft.etaMftTrackMinFilter);

  // Filters below will be used for uncertainties
  Filter mftTrackCollisionIdFilter = (aod::fwdtrack::bestCollisionId >= 0);
  Filter mftTrackDcaXYFilter = (nabs(aod::fwdtrack::bestDCAXY) < configMft.mftMaxDCAxy);
  // Filter mftTrackDcaZFilter = (nabs(aod::fwdtrack::bestDCAZ) < configMft.mftMaxDCAz);

  // =========================
  //      Filters & partitions : MC
  // =========================

  Filter candidateFilterD0Mc = (aod::hf_sel_candidate_d0::isRecoHfFlag >= configCandidates.selectionFlagHf) ||
                               (aod::hf_sel_candidate_d0::isRecoHfFlag >= configCandidates.selectionFlagHf);

  Filter candidateFilterLcMc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= configCandidates.selectionFlagHf) ||
                               (aod::hf_sel_candidate_lc::isSelLcToPiKP >= configCandidates.selectionFlagHf);

  // From Katarina's code, but not sure if I use it
  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < configCollision.zVertexMax;

  // Filter mcParticlesFilter = (nabs(aod::mcparticle::eta) < configCentral.etaCentralTrackMax) &&
  //                            (aod::mcparticle::pt > configCentral.ptCentralTrackMin);

  // I didn't manage to make partitions work with my mixed event, as I am pair my tracks BEFORE looping over collisions
  // I am thus not able to group tracks with sliceBy and can't use this method
  // For now I am fine as I am doing only TPC-MFT correlations and using only McParticles with MFT acceptance
  // However at some point I will have to use tracks from the other side (FV0, FT0-A) and I will have to do something about it
  // TO-DO : either change how I do mixed event, or implement isAcceptedTpcMcParticle, isAcceptedMftMcParticle
  // Partition<aod::McParticles> mcParticlesMft = (aod::mcparticle::eta > configMft.etaMftTrackMin) && (aod::mcparticle::eta < configMft.etaMftTrackMax);
  // Partition<aod::McParticles> mcParticlesTpc = (nabs(aod::mcparticle::eta) < configCentral.etaCentralTrackMax) &&
  //                                             (aod::mcparticle::pt > configCentral.ptCentralTrackMin);

  // =========================
  //      Preslice : DATA
  // =========================

  Preslice<HfCandidatesSelD0> perColD0s = aod::track::collisionId;
  Preslice<HfCandidatesSelLc> perColLcs = aod::track::collisionId;
  Preslice<FilteredMftTracks> perColMftTracks = o2::aod::fwdtrack::collisionId;
  Preslice<FilteredTracksWDcaSel> perColTracks = aod::track::collisionId;

  // =========================
  //      Preslice : MC
  // =========================

  Preslice<MftTracksMcLabels> mftTracksPerCollision = aod::fwdtrack::collisionId;

  //  configurables for containers
  //  TODO: flow of HF will need to be done vs. invariant mass, in the signal and side-band regions
  //        either 1) add invariant mass axis or 2) define several containers for different inv. mass regions
  //        Note: don't forget to check inv. mass separately for D0 and D0bar candidate
  ConfigurableAxis axisMass{"axisMass", {120, 1.5848, 2.1848}, "axis of invariant mass of candidates"};
  ConfigurableAxis binsMixingMultiplicity{"binsMixingMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity bins for event mixing"};
  ConfigurableAxis binsMixingVertex{"binsMixingVertex", {14, -7, 7}, "vertex bins for event mixing"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisEtaAssociated{"axisEtaAssociated", {48, -4, -2}, "eta axis for MFT histograms"};
  ConfigurableAxis axisEtaTrigger{"axisEtaTrigger", {48, -1, 1}, "eta axis for TPC histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0, TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {72, 0, 36}, "pt axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisVertex{"axisVertex", {14, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};

  HistogramRegistry registry{"registry"};

  // Correlation containers used for data
  OutputObj<CorrelationContainer> sameEvent{"sameEvent"};
  OutputObj<CorrelationContainer> mixedEvent{"mixedEvent"};
  OutputObj<CorrelationContainer> sameEventHf{"sameEventHf"};
  OutputObj<CorrelationContainer> mixedEventHf{"mixedEventHf"};

  // Correlation containers used for Monte-Carlo
  OutputObj<CorrelationContainer> sameEventHfMc{"sameEventHfMc"};
  OutputObj<CorrelationContainer> mixedEventHfMc{"mixedEventHfMc"};

  template <DataType DataType, CorrelationCase CorrelationCase, CorrelatedParticles CorrelatedParticles>
  void addHistograms()
  {
    registry.add(Form("%s%s%shEtaTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {axisEtaTrigger}});
    registry.add(Form("%s%s%shPhiTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {axisPhi}});
    registry.add(Form("%s%s%shPtTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {axisPt}});
    registry.add(Form("%s%s%shYieldsTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTrigger}});
    registry.add(Form("%s%s%shEtaPhiTrigger", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH3F, {axisMultiplicity, axisEtaTrigger, axisPhi}});
    registry.add(Form("%s%s%shEtaAssociated", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {axisEtaAssociated}});
    registry.add(Form("%s%s%shPhiAssociated", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH1D, {axisPhi}});
    registry.add(Form("%s%s%shEtaPhiAssociated", WhatDataType[DataType].data(), WhatCorrelationCase[CorrelationCase].data(), WhatParticles[CorrelatedParticles].data()), "", {HistType::kTH3F, {axisMultiplicity, axisEtaAssociated, axisPhi}});
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
    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", configCcdb.noLaterThan.value);
    offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", configCcdb.noLaterThan.value);
    LOGF(info, "Offset for FT0A: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY(), (*offsetFT0)[0].getZ());
    LOGF(info, "Offset for FT0C: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY(), (*offsetFT0)[1].getZ());
    LOGF(info, "Offset for FV0-left: x = %.3f y = %.3f z = %.3f\n", (*offsetFV0)[0].getX(), (*offsetFV0)[0].getY(), (*offsetFV0)[0].getZ());
    LOGF(info, "Offset for FV0-right: x = %.3f y = %.3f z = %.3f\n", (*offsetFV0)[1].getX(), (*offsetFV0)[1].getY(), (*offsetFV0)[1].getZ());

    //  =========================
    //      Event histograms
    //  =========================

    registry.add("Data/hVtxZ", "v_{z} (cm)", {HistType::kTH1D, {axisVertex}});
    registry.add(Form("Data/hMultiplicity_%s", WhatMultiplicityEstimator[configCollision.multiplicityEstimator].data()), "", {HistType::kTH1D, {axisMultiplicity}});

    registry.add("Data/hEventCounter", "hEventCounter", {HistType::kTH1D, {{EventSelectionStep::NEventSelectionSteps, -0.5, +EventSelectionStep::NEventSelectionSteps - 0.5}}});
    std::string labels[EventSelectionStep::NEventSelectionSteps];
    labels[EventSelectionStep::AllEvents] = "all";
    labels[EventSelectionStep::AfterEventSelection] = "after Physics selection";
    registry.get<TH1>(HIST("Data/hEventCounter"))->SetMinimum(0);

    for (int iBin = 0; iBin < EventSelectionStep::NEventSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("Data/hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

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
    registry.get<TH1>(HIST("Data/Mft/hMftTracksSelection"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftTrackSelectionStep::NMftTrackSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("Data/Mft/hMftTracksSelection"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftTracksSelection[iBin].data());
    }

    registry.add("Data/Mft/hReassociationMftTracks", "hReassociationMftTracks", {HistType::kTH1D, {{ReassociationMftTracks::NReassociationMftTracksSteps, -0.5, +ReassociationMftTracks::NReassociationMftTracksSteps - 0.5}}});
    std::string labelsReassociationMftTracks[ReassociationMftTracks::NReassociationMftTracksSteps];
    labelsReassociationMftTracks[ReassociationMftTracks::NotReassociatedMftTracks] = "MFT tracks after track selection";
    labelsReassociationMftTracks[ReassociationMftTracks::ReassociatedMftTracks] = "Reassociated MFT tracks by DCAxy method";
    registry.get<TH1>(HIST("Data/Mft/hReassociationMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < ReassociationMftTracks::NReassociationMftTracksSteps; iBin++) {
      registry.get<TH1>(HIST("Data/Mft/hReassociationMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsReassociationMftTracks[iBin].data());
    }

    registry.add("Data/Mft/hNTracks", "", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/Mft/hNMftTracks", "", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/Mft/hNBestCollisionFwd", "", {HistType::kTH1F, {axisMultiplicity}});

    //  =========================
    //      Declaration of correlation containers and their respective axis
    //  =========================

    std::vector<AxisSpec> const corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                            {axisPtAssoc, "p_{T} (GeV/c)"},
                                            {axisPtTrigger, "p_{T} (GeV/c)"},
                                            {axisMultiplicity, "multiplicity"},
                                            {axisDeltaPhi, "#Delta#varphi (rad)"},
                                            {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> const effAxis = {{axisEtaEfficiency, "#eta"},
                                           {axisPtEfficiency, "p_{T} (GeV/c)"},
                                           {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> const userAxis = {{axisMass, "m_{inv} (GeV/c^{2})"}};

    fv0Det = o2::fv0::Geometry::instance(o2::fv0::Geometry::eUninitialized);

    //  =========================
    //  Initialization of histograms and CorrelationContainers for TpcTpc cases
    //  =========================

    if (doprocessSameTpcTpcChCh) {
      addHistograms<Data, TpcTpc, ChPartChPart>();
      sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
      mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    }

    if (doprocessSameTpcTpcD0Ch) {
      addHistograms<Data, TpcTpc, D0ChPart>();
      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    if (doprocessSameTpcTpcLcCh) {
      addHistograms<Data, TpcTpc, LcChPart>();
      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    //  =========================
    //  Initialization of histograms and CorrelationContainers for TpcMft cases
    //  =========================

    // if (doprocessSameTpcMftChCh || doprocessSameTpcMftChChReassociated || doprocessSameTpcMftChChReassociated3d || doprocessSameTpcMftChChNonAmbiguous) {
    if (doprocessSameTpcMftChCh || doprocessSameTpcMftChChReassociated || doprocessSameTpcMftChChNonAmbiguous) {
      addHistograms<Data, TpcMft, ChPartChPart>();

      // All MFT tracks
      registry.add("Data/Mft/kCFStepAll/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/Mft/kCFStepAll/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      // Only non-ambiguous MFT tracks
      registry.add("Data/Mft/kCFStepTracked/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/Mft/kCFStepTracked/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
      mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    }

    if (doprocessSameTpcMftD0Ch || doprocessSameTpcMftD0ChReassociated) {
      addHistograms<Data, TpcMft, D0ChPart>();

      // All MFT tracks
      registry.add("Data/Mft/kCFStepAll/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/Mft/kCFStepAll/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      // Only non-ambiguous MFT tracks
      registry.add("Data/Mft/kCFStepTracked/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/Mft/kCFStepTracked/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    if (doprocessSameTpcMftLcCh || doprocessSameTpcMftLcChReassociated) {
      addHistograms<Data, TpcMft, LcChPart>();

      // All MFT tracks
      registry.add("Data/Mft/kCFStepAll/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/Mft/kCFStepAll/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      // Only non-ambiguous MFT tracks
      registry.add("Data/Mft/kCFStepTracked/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/Mft/kCFStepTracked/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    //  =========================
    //  Initialization of histograms and CorrelationContainers for TpcFv0a cases
    //  =========================

    if (doprocessSameTpcFv0aChCh) {
      addHistograms<Data, TpcFv0a, ChPartChPart>();

      sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
      mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    }

    if (doprocessSameTpcFv0aD0Ch) {
      addHistograms<Data, TpcFv0a, D0ChPart>();

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    if (doprocessSameTpcFv0aLcCh) {
      addHistograms<Data, TpcFv0a, LcChPart>();

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    //  =========================
    //  Initialization of histograms and CorrelationContainers for MftFv0a cases
    //  =========================

    // if (doprocessSameMftFv0aChCh || doprocessSameMftFv0aChChReassociated || doprocessSameMftFv0aReassociated3d || doprocessSameMftFv0aChChNonAmbiguous) {
    if (doprocessSameMftFv0aChCh || doprocessSameMftFv0aChChReassociated || doprocessSameMftFv0aChChNonAmbiguous) {
      addHistograms<Data, MftFv0a, ChPartChPart>();

      sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
      mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    }

    //  =========================
    //  Initialization of histograms and CorrelationContainers for TpcFt0a cases
    //  =========================

    if (doprocessSameTpcFt0aChCh) {
      addHistograms<Data, TpcFt0a, ChPartChPart>();

      sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
      mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    }

    if (doprocessSameTpcFt0aD0Ch) {
      addHistograms<Data, TpcFt0a, D0ChPart>();

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    if (doprocessSameTpcFt0aLcCh) {
      addHistograms<Data, TpcFt0a, LcChPart>();

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    //  =========================
    //  Initialization of histograms and CorrelationContainers for MftFt0a cases
    //  =========================

    if (doprocessSameMftFt0aChCh || doprocessSameMftFt0aChChReassociated || doprocessSameMftFt0aChChNonAmbiguous) {
      addHistograms<Data, MftFt0a, ChPartChPart>();

      sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
      mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    }

    //  =========================
    //  Initialization of histograms and CorrelationContainers for TpcMft MONTE-CARLO cases
    //  =========================

    if (doprocessSameTpcMftD0ChMcGen) {
      addHistograms<Mc, TpcMft, D0ChPart>();
      sameEventHfMc.setObject(new CorrelationContainer("sameEventHfMc", "sameEventHfMc", corrAxis, effAxis, userAxis));
      mixedEventHfMc.setObject(new CorrelationContainer("mixedEventHfMc", "mixedEventHfMc", corrAxis, effAxis, userAxis));
    }

    if (doprocessSameTpcMftLcChMcGen) {
      addHistograms<Mc, TpcMft, LcChPart>();
      sameEventHfMc.setObject(new CorrelationContainer("sameEventHfMc", "sameEventHfMc", corrAxis, effAxis, userAxis));
      mixedEventHfMc.setObject(new CorrelationContainer("mixedEventHfMc", "mixedEventHfMc", corrAxis, effAxis, userAxis));
    }

  } // End of init() function

  // =========================
  //      Quality assessment functions
  // =========================

  template <DataType DataType, CorrelationCase CorrelationCase, CorrelatedParticles CorrelatedParticles>
  void fillTriggerQa(float multiplicity, float const& eta, float const& phi, float const& pt)
  {
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hPtTrigger"), pt);
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaTrigger"), eta);
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hPhiTrigger"), phi);
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hYieldsTrigger"), multiplicity, pt, eta);
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaPhiTrigger"), multiplicity, eta, phi);
  }

  template <DataType DataType, CorrelationCase CorrelationCase, CorrelatedParticles CorrelatedParticles>
  void fillAssociatedQa(float multiplicity, float const& eta, float const& phi)
  {
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaAssociated"), eta);
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hPhiAssociated"), phi);
    registry.fill(HIST(WhatDataType[DataType]) + HIST(WhatCorrelationCase[CorrelationCase]) + HIST(WhatParticles[CorrelatedParticles]) + HIST("hEtaPhiAssociated"), multiplicity, eta, phi);
  }

  // =========================
  //      Helper functions
  // =========================

  HfProngSpecies getSpecies(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case PDG_t::kPiPlus: // positive or negative pion
        return HfProngSpecies::Pion;
      case PDG_t::kKPlus: // positive or negative kaon
        return HfProngSpecies::Kaon;
      case PDG_t::kProton: // proton or proton bar
        return HfProngSpecies::Proton;
      default: // NOTE. The efficiency histogram is hardcoded to contain 4 species. Anything special will have the last slot.
        return HfProngSpecies::NHfProngSpecies;
    }
  }

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

  double getPhiFT0(uint chno, int i)
  {
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + (*offsetFT0)[i].getX(), chPos.Y() + (*offsetFT0)[i].getY());
  }

  double getPhiFV0(unsigned int chno) const
  {
    int cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool const isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);
    float offsetX, offsetY;
    if (isChnoInLeft) {
      offsetX = (*offsetFV0)[0].getX();
      offsetY = (*offsetFV0)[0].getY();
    } else {
      offsetX = (*offsetFV0)[1].getX();
      offsetY = (*offsetFV0)[1].getY();
    }

    o2::fv0::Point3Dsimple chPos{};
    chPos = fv0Det->getReadoutCenter(chno);

    // if (configTask.isReadoutCenter)
    //   chPos = fv0Det->getReadoutCenter(chno);
    // else
    //   chPos = fv0Det->getCellCenter(chno);

    return RecoDecay::phi(chPos.x + offsetX, chPos.y + offsetY);
  }

  double getEtaFT0(uint chno, int i)
  {
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  double getEtaFV0(unsigned int chno) const
  {
    int cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool const isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);
    float offsetX, offsetY, offsetZ;
    if (isChnoInLeft) {
      offsetX = (*offsetFV0)[0].getX();
      offsetY = (*offsetFV0)[0].getY();
      offsetZ = (*offsetFV0)[0].getZ();
    } else {
      offsetX = (*offsetFV0)[1].getX();
      offsetY = (*offsetFV0)[1].getY();
      offsetZ = (*offsetFV0)[1].getZ();
    }

    o2::fv0::Point3Dsimple chPos{};
    chPos = fv0Det->getReadoutCenter(chno);
    // if (configTask.isReadoutCenter)
    //   chPos = fv0Det->getReadoutCenter(chno);
    // else
    //   chPos = fv0Det->getCellCenter(chno);

    auto x = chPos.x + offsetX;
    auto y = chPos.y + offsetY;
    auto z = chPos.z + offsetZ;
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
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

    if (!configTask.processMc) {
      if (!collision.sel8()) {
        return false;
      }
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

    if (fillHistograms) {
      registry.fill(HIST("Data/hEventCounter"), EventSelectionStep::AfterEventSelection);
    }

    registry.fill(HIST("Data/hVtxZ"), collision.posZ());

    return true;
  }

  //  TODO: Check how to put this into a Filter
  template <typename TTrack>
  bool isAcceptedCandidate(TTrack const& candidate)
  {
    auto etaCandidate = candidate.eta();

    if constexpr (std::is_same_v<HfCandidatesSelLc, TTrack>) { // For now, that means we do LambdaC
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        return false;
      }
      if (configCandidates.etaCandidateMax >= 0. && std::abs(etaCandidate) > configCandidates.etaCandidateMax) {
        return false;
      }
      if (configCandidates.yCandRecoMax >= 0. && std::abs(hfHelper.yLc(candidate)) > configCandidates.yCandRecoMax) {
        return false;
      }
      return true;
    } else { // For now, that means we do D0
      // Doesn't this exclude D0bar ?
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        return false;
      }
      if (configCandidates.etaCandidateMax >= 0. && std::abs(etaCandidate) > configCandidates.etaCandidateMax) {
        return false;
      }
      if (configCandidates.yCandRecoMax >= 0. && std::abs(hfHelper.yD0(candidate)) > configCandidates.yCandRecoMax) {
        return false;
      }
      return true;
    }
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

  // I am not sure if to template McParticles is useful, I'll address this when doing the MC Gen case of HF-h correlations
  template <typename TMcTrack>
  bool isAcceptedMcCandidate(TMcTrack& mcCandidate)
  {
    auto etaCandidate = mcCandidate.eta();

    if constexpr (std::is_same_v<McParticles2ProngMatched, TMcTrack>) { // For now, that means we do D0
      if (std::abs(mcCandidate.flagMcMatchGen()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {

        if (configCandidates.etaCandidateMax >= 0. && std::abs(etaCandidate) > configCandidates.etaCandidateMax) {
          return false;
        }

        if (configCandidates.yCandGenMax >= 0. && std::abs(RecoDecay::y(mcCandidate.pVector(), o2::constants::physics::MassD0)) > configCandidates.yCandGenMax) {
          return false;
        }

        // Later on, if I want to add prompt/non-prompt selection, below is how to select prompt only
        // if (!(particle.originMcGen() == RecoDecay::OriginType::Prompt)){
        //   return false;
        // }
      }
    } else { // For now, that means we do LambdaC
      if (std::abs(mcCandidate.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {

        if (configCandidates.etaCandidateMax >= 0. && std::abs(etaCandidate) > configCandidates.etaCandidateMax) {
          return false;
        }

        if (configCandidates.yCandGenMax >= 0. && std::abs(RecoDecay::y(mcCandidate.pVector(), o2::constants::physics::MassLambdaCPlus)) > configCandidates.yCandGenMax) {
          return false;
        }

        // Later on, if I want to add prompt/non-prompt selection, below is how to select prompt only
        // if (!(particle.originMcGen() == RecoDecay::OriginType::Prompt)){
        //   return false;
        // }
      }
    }

    return true;
  }

  // I am not sure if to template McParticles is useful, I'll address this when doing the MC Gen case of HF-h correlations
  template <typename TMcParticle>
  bool isAcceptedMftMcParticle(TMcParticle& mcParticle)
  {
    //  remove MC particles with charge = 0
    TParticlePDG* pdgparticle = pdg->GetParticle(mcParticle.pdgCode());
    if (pdgparticle != nullptr) {
      if (pdgparticle->Charge() == 0) {
        return false;
      }
    }

    /*
    //  MC particle has to be primary
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
      return mcParticle.isPhysicalPrimary();
    }
    */

    if (mcParticle.eta() > configMft.etaMftTrackMax || mcParticle.eta() < configMft.etaMftTrackMin) {
      return false;
    }

    // return true;
    return mcParticle.isPhysicalPrimary();
  }

  // ===============================================================================================================================================================================
  // ===============================================================================================================================================================================
  //      Correlation functions
  // ===============================================================================================================================================================================
  // ===============================================================================================================================================================================

  // ===============================================================================================================================================================================
  //      fillCorrelations
  // ===============================================================================================================================================================================

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelations(TTarget target, CorrelationContainer::CFStep step,
                        TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                        float multiplicity, float posZ, bool sameEvent)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    // To avoid filling associated tracks QA many times
    //  I fill it only for the first trigger track of the collision
    auto loopCounter = 0;

    // TRIGGER PARTICLE
    for (const auto& track1 : tracks1) {

      loopCounter++;

      float const eta1 = track1.eta();
      float const pt1 = track1.pt();
      float const phi1 = track1.phi();

      //  TODO: add getter for NUE trigger efficiency here

      //  calculating inv. mass to be filled into the container below
      bool fillingHFcontainer = false;
      double invmass = 0;
      if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
        //  TODO: Check how to put this into a Filter -> Pretty sure it cannot be a filter
        if (!isAcceptedCandidate(track1)) {
          continue;
        }
        fillingHFcontainer = true;
        if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // If D0
          invmass = hfHelper.invMassD0ToPiK(track1);
        } else { // If Lc
          invmass = hfHelper.invMassLcToPKPi(track1);
        }
      }

      // Selections for MC GENERATED
      if constexpr (std::is_same_v<McParticles2ProngMatched, TTracksTrig> || std::is_same_v<McParticles3ProngMatched, TTracksTrig>) {
        //  TODO: Check how to put this into a Filter -> Pretty sure it cannot be a filter
        if (!isAcceptedMcCandidate(track1)) {
          continue;
        }
        fillingHFcontainer = true;
        if constexpr (std::is_same_v<McParticles2ProngMatched, TTracksTrig>) { // If D0
          invmass = o2::constants::physics::MassD0;
        } else { // If Lc
          invmass = o2::constants::physics::MassLambdaCPlus;
        }
      }

      //  fill single-track distributions
      if (!fillingHFcontainer) { // if not HF-h case
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      // FILL QA PLOTS for trigger particle
      if (sameEvent && (step == CorrelationContainer::kCFStepReconstructed)) {
        if (!configTask.processMc) {                                        // If DATA
          if constexpr (!std::is_same_v<FilteredMftTracks, TTracksAssoc>) { // IF TPC-TPC case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-TPC D0-h
              fillTriggerQa<Data, TpcTpc, D0ChPart>(multiplicity, eta1, phi1, pt1);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-TPC Lc-h
              fillTriggerQa<Data, TpcTpc, LcChPart>(multiplicity, eta1, phi1, pt1);
            } else { // IF NEITHER D0 NOR LC -> TPC-TPC h-h
              fillTriggerQa<Data, TpcTpc, ChPartChPart>(multiplicity, eta1, phi1, pt1);
            }
          } else {                                                          // IF TPC-MFT case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-MFT D0-h
              fillTriggerQa<Data, TpcMft, D0ChPart>(multiplicity, eta1, phi1, pt1);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-MFT Lc-h
              fillTriggerQa<Data, TpcMft, LcChPart>(multiplicity, eta1, phi1, pt1);
            } else { // IF NEITHER D0 NOR LC -> TPC-MFT h-h
              fillTriggerQa<Data, TpcMft, ChPartChPart>(multiplicity, eta1, phi1, pt1);
            } // end of if condition for TPC-TPC or TPC-MFT case
          }
          // Maybe I won't need it for MC (first files are way lighter in MC, but also I need to loop over all tracks in MC GEN)
        }
      }

      // ASSOCIATED PARTICLE
      for (const auto& track2 : tracks2) {

        // apply cuts for MFT tracks
        if constexpr (std::is_same_v<FilteredMftTracks, TTracksAssoc>) {

          if (sameEvent && loopCounter == 1) { // To avoid double counting, we fill the plots only the first time
            registry.fill(HIST("Data/Mft/hMftTracksSelection"), MftTrackSelectionStep::NoSelection);

            if (!isAcceptedMftTrack(track2, true)) {
              continue;
            }
          } else { // After the first loop, we don't fill the plots anymore but still do the selection
            if (!isAcceptedMftTrack(track2, false)) {
              continue;
            }
          }
        }

        //  case of h-h correlations where the two types of tracks are the same
        //  this avoids autocorrelations and double counting of particle pairs
        if constexpr (std::is_same_v<TTracksAssoc, TTracksTrig>) {
          if (track1.index() <= track2.index()) {
            continue;
          }
        }

        //  in case of HF-h correlations, remove candidate daughters from the pool of associated hadrons
        //  with which the candidate is being correlated (will not have to do it for TPC-MFT case)
        if constexpr (!std::is_same_v<FilteredMftTracks, TTracksAssoc>) { // if NOT TPC-MFT case -> TPC-TPC case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // Remove the 2 prong daughters
            if ((track1.prong0Id() == track2.globalIndex()) || (track1.prong1Id() == track2.globalIndex())) {
              continue;
            }
          }
          if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // Remove the 3 prong daughters
            if ((track1.prong0Id() == track2.globalIndex()) || (track1.prong1Id() == track2.globalIndex()) || (track1.prong2Id() == track2.globalIndex())) {
              continue;
            }
          }
        }

        //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
        // if (configTask.processMc) {
        if constexpr (std::is_same_v<McParticles, TTracksTrig> || std::is_same_v<McParticles, TTracksAssoc>) {
          if (!isAcceptedMftMcParticle(track2)) {
            continue;
          }
        }

        float const eta2 = track2.eta();
        float const pt2 = track2.pt();
        float phi2 = track2.phi();
        o2::math_utils::bringTo02Pi(phi2);

        //  TODO: add getter for NUE associated efficiency here

        //  TODO: add pair cuts on phi*

        float deltaPhi = phi1 - phi2;
        //  set range of delta phi in (-pi/2 , 3/2*pi)
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

        if (!fillingHFcontainer) {
          //  fill pair correlations
          target->getPairHist()->Fill(step, eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);
        } else {
          target->getPairHist()->Fill(step, eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ, invmass,
                                      triggerWeight * associatedWeight);
        }

        // FILL QA PLOTS for associated particle
        if (sameEvent && (loopCounter == 1) && (step == CorrelationContainer::kCFStepReconstructed)) {
          if constexpr (!std::is_same_v<FilteredMftTracks, TTracksAssoc>) { // IF TPC-TPC case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-TPC D0-h
              fillAssociatedQa<Data, TpcTpc, D0ChPart>(multiplicity, eta2, phi2);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-TPC Lc-h
              fillAssociatedQa<Data, TpcTpc, LcChPart>(multiplicity, eta2, phi2);
            }
            // No if condition if it is h-h, because it would be the same plots than for the trigger particle
          } else {                                                          // IF TPC-MFT case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-MFT D0-h
              fillAssociatedQa<Data, TpcMft, D0ChPart>(multiplicity, eta2, phi2);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-MFT Lc-h
              fillAssociatedQa<Data, TpcMft, LcChPart>(multiplicity, eta2, phi2);
            } else { // IF NEITHER D0 NOR LC -> TPC-MFT h-h
              fillAssociatedQa<Data, TpcMft, ChPartChPart>(multiplicity, eta2, phi2);
            } // end of if condition for TPC-TPC or TPC-MFT case
          }
        }

        if (sameEvent && (loopCounter == 1) && std::is_same_v<FilteredMftTracks, TTracksAssoc>) {
          // FILL USUAL MFT DISTRIBUTIONS
          registry.fill(HIST("Data/Mft/kCFStepAll/hEta"), eta2);
          registry.fill(HIST("Data/Mft/kCFStepAll/hPhi"), phi2);
        }

      } // end of loop over tracks2
    } // end of loop over tracks 1
  }

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelationsReassociatedMftTracks(TTarget target, CorrelationContainer::CFStep step,
                                             TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                                             float multiplicity, float posZ, bool sameEvent, bool cutAmbiguousTracks)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    // To avoid filling associated tracks QA many times
    //  I fill it only for the first trigger track of the collision
    auto loopCounter = 0;

    // TRIGGER PARTICLE
    for (const auto& track1 : tracks1) {

      loopCounter++;

      float const eta1 = track1.eta();
      float const pt1 = track1.pt();
      float const phi1 = track1.phi();

      bool fillingHFcontainer = false;
      double invmass = 0;
      if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
        if (!isAcceptedCandidate(track1)) {
          continue;
        }
        fillingHFcontainer = true;
        if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // If D0
          invmass = hfHelper.invMassD0ToPiK(track1);
        } else { // If Lc
          invmass = hfHelper.invMassLcToPKPi(track1);
        }
      }

      //  fill single-track distributions
      if (!fillingHFcontainer) { // if not HF-h case
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      // FILL QA PLOTS for trigger particle
      if (sameEvent && (!cutAmbiguousTracks)) {
        if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) {
          fillTriggerQa<Data, TpcMft, D0ChPart>(multiplicity, eta1, phi1, pt1);
        } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
          fillTriggerQa<Data, TpcMft, LcChPart>(multiplicity, eta1, phi1, pt1);
        } else {
          fillTriggerQa<Data, TpcMft, ChPartChPart>(multiplicity, eta1, phi1, pt1);
        }
      }

      // ASSOCIATED PARTICLE
      for (const auto& track2 : tracks2) {

        // Fill QA plot for all MFT tracks () (only if cutAmbiguousTracks is false to avoid double counting)
        if (!cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
          registry.fill(HIST("Data/Mft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
        }

        auto reassociatedMftTrack = track2.template mfttrack_as<FilteredMftTracks>();

        if (!isAcceptedMftTrack(reassociatedMftTrack, false)) {
          continue;
        }

        // Fill QA plot for MFT tracks after physical selection (eta + clusters)
        if (!cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
          registry.fill(HIST("Data/Mft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);
          registry.fill(HIST("Data/Mft/hReassociationMftTracks"), ReassociationMftTracks::NotReassociatedMftTracks);
        }

        // We check if the track is ambiguous or non-ambiguous (QA plots are filled in isAmbiguousMftTrack)
        // Fill plots only if cutAmbiguousTracks is false (to avoid double counting)
        if (isAmbiguousMftTrack(track2, (!cutAmbiguousTracks && sameEvent && (loopCounter == 1)))) {
          // If the MFT track is ambiguous we may cut or not on the ambiguous track
          if (cutAmbiguousTracks) {
            continue;
          }
        }

        if (reassociatedMftTrack.collisionId() != track2.bestCollisionId()) {
          if (sameEvent && (loopCounter == 1)) {
            registry.fill(HIST("Data/Mft/hReassociationMftTracks"), ReassociationMftTracks::ReassociatedMftTracks);
          }
        }

        //  in case of HF-h correlations, remove candidate daughters from the pool of associated hadrons
        //  with which the candidate is being correlated (will not have to do it for TPC-MFT case)
        if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // Remove the 2 prong daughters
          if ((track1.prong0Id() == reassociatedMftTrack.globalIndex()) || (track1.prong1Id() == reassociatedMftTrack.globalIndex())) {
            continue;
          }
        }
        if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // Remove the 3 prong daughters
          if ((track1.prong0Id() == reassociatedMftTrack.globalIndex()) || (track1.prong1Id() == reassociatedMftTrack.globalIndex()) || (track1.prong2Id() == reassociatedMftTrack.globalIndex())) {
            continue;
          }
        }

        float const eta2 = reassociatedMftTrack.eta();
        float const pt2 = reassociatedMftTrack.pt();
        float phi2 = reassociatedMftTrack.phi();
        o2::math_utils::bringTo02Pi(phi2);

        //  TODO: add getter for NUE associated efficiency here

        //  TODO: add pair cuts on phi*

        float deltaPhi = phi1 - phi2;
        //  set range of delta phi in (-pi/2 , 3/2*pi)
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

        if (!fillingHFcontainer) {
          target->getPairHist()->Fill(step, eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);
        } else {
          target->getPairHist()->Fill(step, eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ, invmass,
                                      triggerWeight * associatedWeight);
        }

        // FILL QA PLOTS for associated particle
        if (sameEvent && (loopCounter == 1) && (!cutAmbiguousTracks)) {
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) {
            fillAssociatedQa<Data, TpcMft, D0ChPart>(multiplicity, eta2, phi2);
          } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
            fillAssociatedQa<Data, TpcMft, LcChPart>(multiplicity, eta2, phi2);
          } else {
            fillAssociatedQa<Data, TpcMft, ChPartChPart>(multiplicity, eta2, phi2);
          }
        }

        // QA plots for basic MFT distributions for non-ambiguous tracks only (kCFStepTracked)
        if (cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
          registry.fill(HIST("Data/Mft/kCFStepTracked/hEta"), eta2);
          registry.fill(HIST("Data/Mft/kCFStepTracked/hPhi"), phi2);
        }

      } // end of loop over tracks2
    } // end of loop over tracks 1
  }

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc, typename TFits>
  void fillCorrelationsFIT(TTarget target, CorrelationContainer::CFStep step,
                           TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TFits const&,
                           float multiplicity, float posZ, bool sameEvent)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    // To avoid filling associated tracks QA many times
    //  I fill it only for the first trigger track of the collision
    auto loopCounter = 0;

    // TRIGGER PARTICLE
    for (auto const& track1 : tracks1) {

      loopCounter++;

      float const eta1 = track1.eta();
      float const pt1 = track1.pt();
      float phi1 = track1.phi();
      if constexpr (std::is_same_v<FilteredMftTracks, TTracksTrig>) {
        o2::math_utils::bringTo02Pi(phi1);
      }

      bool fillingHFcontainer = false;
      double invmass = 0;
      if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
        if (!isAcceptedCandidate(track1)) {
          continue;
        }
        fillingHFcontainer = true;
        if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // If D0
          invmass = hfHelper.invMassD0ToPiK(track1);
        } else { // If Lc
          invmass = hfHelper.invMassLcToPKPi(track1);
        }
      }

      //  fill single-track distributions
      if (!fillingHFcontainer) { // if not HF-h case
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      // FILL QA PLOTS for trigger particle
      if (sameEvent && (step == CorrelationContainer::kCFStepReconstructed)) {
        if (!configTask.processMc) {                                        // If DATA
          if constexpr (!std::is_same_v<FilteredMftTracks, TTracksTrig>) {  // If not FilteredMftTracks as trigger -> TPC-FV0a correlations
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-FV0a D0-h
              if constexpr (std::is_same_v<aod::FV0As, TFits>) {            // IF NEITHER D0 NOR LC ->
                fillTriggerQa<Data, TpcFv0a, D0ChPart>(multiplicity, eta1, phi1, pt1);
              } else if constexpr (std::is_same_v<aod::FT0s, TFits>) {
                fillTriggerQa<Data, TpcFt0a, D0ChPart>(multiplicity, eta1, phi1, pt1);
              }
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-FV0a Lc-h
              if constexpr (std::is_same_v<aod::FV0As, TFits>) {                   // IF NEITHER D0 NOR LC ->
                fillTriggerQa<Data, TpcFv0a, LcChPart>(multiplicity, eta1, phi1, pt1);
              } else if constexpr (std::is_same_v<aod::FT0s, TFits>) {
                fillTriggerQa<Data, TpcFt0a, LcChPart>(multiplicity, eta1, phi1, pt1);
              }
            } else if constexpr (std::is_same_v<aod::FV0As, TFits>) { // IF NEITHER D0 NOR LC -
              fillTriggerQa<Data, TpcFv0a, ChPartChPart>(multiplicity, eta1, phi1, pt1);
            } else if constexpr (std::is_same_v<aod::FT0s, TFits>) {
              fillTriggerQa<Data, TpcFt0a, ChPartChPart>(multiplicity, eta1, phi1, pt1);
            }
          } else { // If FilteredMftTracks as trigger
            if constexpr (std::is_same_v<aod::FV0As, TFits>) {
              fillTriggerQa<Data, MftFv0a, ChPartChPart>(multiplicity, eta1, phi1, pt1);
            } else if constexpr (std::is_same_v<aod::FT0s, TFits>) {
              fillTriggerQa<Data, MftFt0a, ChPartChPart>(multiplicity, eta1, phi1, pt1);
            }
          }
        }
      } // end of if condition to fill QA plots for trigger particle

      // ASSOCIATED PARTICLE IF USING FV0
      if constexpr (std::is_same_v<aod::FV0As, TFits>) {
        for (std::size_t indexChannel = 0; indexChannel < tracks2.channel().size(); indexChannel++) {

          auto channelId = tracks2.channel()[indexChannel];
          // float fv0Amplitude = tracks2.amplitude()[indexChannel];
          // if (fv0Amplitude <= 0) {
          //   continue;
          // }

          auto phi2 = getPhiFV0(channelId);
          auto eta2 = getEtaFV0(channelId);

          float deltaPhi = phi1 - phi2;
          //  set range of delta phi in (-pi/2 , 3/2*pi)
          deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

          if (!fillingHFcontainer) {
            target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ,
                                        triggerWeight * associatedWeight);
          } else {
            target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ, invmass,
                                        triggerWeight * associatedWeight);
          }

          // FILL QA PLOTS for associated particle
          if (sameEvent && (loopCounter == 1) && (step == CorrelationContainer::kCFStepReconstructed)) {
            if constexpr (!std::is_same_v<FilteredMftTracks, TTracksTrig>) {  // If not FilteredMftTracks as trigger -> TPC-FV0a correlations
              if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-FV0a D0-h
                fillAssociatedQa<Data, TpcFv0a, D0ChPart>(multiplicity, eta2, phi2);
              } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-FV0a Lc-h
                fillAssociatedQa<Data, TpcFv0a, LcChPart>(multiplicity, eta2, phi2);
              } else if constexpr (std::is_same_v<aod::FV0As, TFits>) { // IF NEITHER D0 NOR LC -> ch. part. - ch. part
                fillAssociatedQa<Data, TpcFv0a, ChPartChPart>(multiplicity, eta2, phi2);
              }
            } else { // If FilteredMftTracks as trigger -> MFT-FV0a (non reassoc/ambiguous) correlations
              fillAssociatedQa<Data, MftFv0a, ChPartChPart>(multiplicity, eta2, phi2);
            }
          } // end of if condition to fill QA plots for associated particle
        } // end of loop over FV0 channel indices
      } // end of if condition for FV0s

      // ASSOCIATED PARTICLE IF USING FT0
      if constexpr (std::is_same_v<aod::FT0s, TFits>) {
        for (std::size_t indexChannel = 0; indexChannel < tracks2.channelA().size(); indexChannel++) {

          auto channelId = tracks2.channelA()[indexChannel];
          // float fv0Amplitude = tracks2.amplitudeA()[indexChannel];
          // if (fv0Amplitude <= 0) {
          //   continue;
          // }

          auto phi2 = getPhiFT0(channelId, 0);
          auto eta2 = getEtaFT0(channelId, 0);

          float deltaPhi = phi1 - phi2;
          //  set range of delta phi in (-pi/2 , 3/2*pi)
          deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

          if (!fillingHFcontainer) {
            target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ,
                                        triggerWeight * associatedWeight);
          } else {
            target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ, invmass,
                                        triggerWeight * associatedWeight);
          }

          // FILL QA PLOTS for associated particle
          if (sameEvent && (loopCounter == 1) && (step == CorrelationContainer::kCFStepReconstructed)) {
            if constexpr (!std::is_same_v<FilteredMftTracks, TTracksTrig>) {  // If not FilteredMftTracks as trigger -> TPC-Ft0a correlations
              if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-FV0a D0-h
                fillAssociatedQa<Data, TpcFt0a, D0ChPart>(multiplicity, eta2, phi2);
              } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-FV0a Lc-h
                fillAssociatedQa<Data, TpcFt0a, LcChPart>(multiplicity, eta2, phi2);
              } else if constexpr (std::is_same_v<aod::FT0s, TFits>) { // IF NEITHER D0 NOR LC -> ch. part. - ch. part
                fillAssociatedQa<Data, TpcFt0a, ChPartChPart>(multiplicity, eta2, phi2);
              }
            } else if constexpr (std::is_same_v<aod::FT0s, TFits>) { // If FilteredMftTracks as trigger -> MFT-Ft0a (non reassoc/ambiguous) correlations
              fillAssociatedQa<Data, MftFt0a, ChPartChPart>(multiplicity, eta2, phi2);
            }
          } // end of if condition to fill QA plots for associated particle
        } // end of loop over FT0 channel indices
      } // end of if condition for FT0s
    } // end of loop over tracks 1
  }

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc, typename TFits>
  void fillCorrelationsFITReassociatedMftTracks(TTarget target, CorrelationContainer::CFStep step,
                                                TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TFits const&,
                                                float multiplicity, float posZ, bool sameEvent, bool cutAmbiguousTracks)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    // To avoid filling associated tracks QA many times
    //  I fill it only for the first trigger track of the collision
    auto loopCounter = 0;

    // TRIGGER PARTICLE
    for (auto const& track1 : tracks1) {
      loopCounter++;

      auto reassociatedMftTrack = track1.template mfttrack_as<FilteredMftTracks>();

      if (!isAcceptedMftTrack(reassociatedMftTrack, false)) {
        continue;
      }

      // We check if the track is ambiguous or non-ambiguous (QA plots are filled in isAmbiguousMftTrack)
      // Fill plots only if cutAmbiguousTracks is false (to avoid double counting)
      if (isAmbiguousMftTrack(track1, false)) {
        // If the MFT track is ambiguous we may cut or not on the ambiguous track
        if (cutAmbiguousTracks) {
          continue;
        }
      }

      float const eta1 = reassociatedMftTrack.eta();
      float const pt1 = reassociatedMftTrack.pt();
      float phi1 = reassociatedMftTrack.phi();
      o2::math_utils::bringTo02Pi(phi1);

      target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);

      // FILL QA PLOTS for trigger particle
      if (sameEvent && (step == CorrelationContainer::kCFStepReconstructed)) {
        if constexpr (std::is_same_v<aod::FV0As, TFits>) {
          fillTriggerQa<Data, MftFv0a, ChPartChPart>(multiplicity, eta1, phi1, pt1);
        } else if constexpr (std::is_same_v<aod::FT0s, TFits>) {
          fillTriggerQa<Data, MftFt0a, ChPartChPart>(multiplicity, eta1, phi1, pt1);
        }
      } // end of if condition to fill QA plots for trigger particle

      // ASSOCIATED PARTICLE FOR FV0s
      if constexpr (std::is_same_v<aod::FV0As, TFits>) {
        for (std::size_t indexChannel = 0; indexChannel < tracks2.channel().size(); indexChannel++) {

          auto channelId = tracks2.channel()[indexChannel];
          // float fv0Amplitude = tracks2.amplitude()[indexChannel];
          // if (fv0Amplitude <= 0) {
          //   continue;
          // }

          auto phi2 = getPhiFV0(channelId);
          auto eta2 = getEtaFV0(channelId);

          float deltaPhi = phi1 - phi2;
          //  set range of delta phi in (-pi/2 , 3/2*pi)
          deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

          target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);

          // FILL QA PLOTS for associated particle
          if (sameEvent && (loopCounter == 1) && (step == CorrelationContainer::kCFStepReconstructed)) {
            fillAssociatedQa<Data, MftFv0a, ChPartChPart>(multiplicity, eta2, phi2);
          } // end of if condition to fill QA plots for associated particle
        } // end of loop over FV0 channel indices
      } // end of if condition for FV0s

      // ASSOCIATED PARTICLE FOR FT0s
      if constexpr (std::is_same_v<aod::FT0s, TFits>) {
        for (std::size_t indexChannel = 0; indexChannel < tracks2.channelA().size(); indexChannel++) {

          auto channelId = tracks2.channelA()[indexChannel];
          // float ft0Amplitude = tracks2.amplitudeA()[indexChannel];
          // if (ft0Amplitude <= 0) {
          //   continue;
          // }

          auto phi2 = getPhiFT0(channelId, 0);
          auto eta2 = getEtaFT0(channelId, 0);

          float deltaPhi = phi1 - phi2;
          //  set range of delta phi in (-pi/2 , 3/2*pi)
          deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

          target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);

          // FILL QA PLOTS for associated particle
          if (sameEvent && (loopCounter == 1) && (step == CorrelationContainer::kCFStepReconstructed)) {
            fillAssociatedQa<Data, MftFt0a, ChPartChPart>(multiplicity, eta2, phi2);
          } // end of if condition to fill QA plots for associated particle
        } // end of loop over FT0 channel indices
      } // end of if condition for FT0s
    } // end of loop over tracks 1
  }

  // ===============================================================================================================================================================================
  //      mixCollisions for RECONSTRUCTED events
  // ===============================================================================================================================================================================

  template <typename TCollisions, typename TTracksTrig, typename TTracksAssoc>
  void mixCollisions(TCollisions const& collisions, CorrelationContainer::CFStep step,
                     TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                     OutputObj<CorrelationContainer>& corrContainer)
  {
    auto getMultiplicity = [this](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = getMultiplicityEstimator(collision, false);
      return multiplicity;
    };

    // The first one that I call "Data" should work for data and mc rec
    using BinningTypeData = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;

    BinningTypeData const binningWithTracksSize{{getMultiplicity}, {binsMixingVertex, binsMixingMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<TCollisions, TTracksTrig, TTracksAssoc, BinningTypeData> const pair{binningWithTracksSize, configTask.nMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if constexpr (!std::is_same_v<FilteredMcCollisions, TCollisions>) { // if NOT MC -> do collision cut
        if (!(isAcceptedCollision(collision1, false))) {
          continue;
        }
        if (!(isAcceptedCollision(collision2, false))) {
          continue;
        }
      }

      const auto multiplicity = getMultiplicityEstimator(collision1, false);

      corrContainer->fillEvent(multiplicity, step);
      fillCorrelations(corrContainer, step, tracks1, tracks2, multiplicity, collision1.posZ(), false);
    }
  }

  template <typename TCollisions, typename TTracksTrig, typename TTracksAssoc, typename TPreslice>
  void mixCollisionsFIT(TCollisions const& collisions, CorrelationContainer::CFStep step,
                        TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TPreslice const& preslice,
                        OutputObj<CorrelationContainer>& corrContainer)
  {
    auto getMultiplicity = [this](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = getMultiplicityEstimator(collision, false);
      return multiplicity;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    MixedBinning const binningOnVtxAndMult{{getMultiplicity}, {binsMixingVertex, binsMixingMultiplicity}, true};

    for (auto const& [collision1, collision2] : soa::selfCombinations(binningOnVtxAndMult, configTask.nMixedEvents, -1, collisions, collisions)) {

      if (!isAcceptedCollision(collision1) || !isAcceptedCollision(collision2)) {
        continue;
      }

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      if constexpr (std::is_same_v<aod::FV0As, TTracksAssoc>) { // IF ASSOCIATED PARTICLE FROM FV0A
        if (collision1.has_foundFV0() && collision2.has_foundFV0()) {

          const auto multiplicity = getMultiplicityEstimator(collision1, false);

          auto slicedTriggerTracks = tracks1.sliceBy(preslice, collision1.globalIndex());
          const auto& fv0 = collision2.foundFV0();

          corrContainer->fillEvent(multiplicity, step);
          fillCorrelationsFIT(corrContainer, step, slicedTriggerTracks, fv0, tracks2, multiplicity, collision1.posZ(), false);
        }
      } // end of if condition for FV0s

      if constexpr (std::is_same_v<aod::FT0s, TTracksAssoc>) {
        if (collision1.has_foundFT0() && collision2.has_foundFT0()) {

          const auto multiplicity = getMultiplicityEstimator(collision1, false);

          auto slicedTriggerTracks = tracks1.sliceBy(preslice, collision1.globalIndex());
          const auto& ft0 = collision2.foundFT0();

          corrContainer->fillEvent(multiplicity, step);
          fillCorrelationsFIT(corrContainer, step, slicedTriggerTracks, ft0, tracks2, multiplicity, collision1.posZ(), false);
        }
      } // end of if condition for FT0s
    } // end of for loop
  }

  /*
  template <typename TCollisions, typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisionsReassociatedMftTracks(TCollisions const& collisions, int step,
                     TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                     TLambda getPartsSize,
                     OutputObj<CorrelationContainer>& corrContainer,
                     bool cutAmbiguousTracks)
  {

    // The first one that I call "Data" should work for data and mc rec
    using BinningTypeData = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::collision::PosZ, decltype(getPartsSize)>;

    BinningTypeData binningWithTracksSize{{getPartsSize}, {binsMixingVertex, binsMixingMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<TCollisions, TTracksTrig, TTracksAssoc, BinningTypeData> pair{binningWithTracksSize, configTask.nMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if constexpr (!std::is_same_v<FilteredMcCollisions, TCollisions>) { // if NOT MC -> do collision cut
        if (!(isAcceptedCollision(collision1, false))) {
          continue;
        }
        if (!(isAcceptedCollision(collision2, false))) {
          continue;
        }
      }

      const auto multiplicity = collision1.multNTracksPV();

      corrContainer->fillEvent(multiplicity, step);
      fillCorrelationsReassociatedMftTracks(corrContainer, step, tracks1, tracks2, true, multiplicity, collision1.posZ(), false, cutAmbiguousTracks, field );
    }
  }
  */

  // ===============================================================================================================================================================================
  //      mixCollisions for GENERATED events
  // ===============================================================================================================================================================================

  template <typename TMcCollisions, typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisionsMcTruth(TMcCollisions const& mcCollisions,
                            TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                            TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  {
    using BinningTypeMcTruth = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::mccollision::PosZ, decltype(getPartsSize)>;

    BinningTypeMcTruth const binningWithTracksSize{{getPartsSize}, {binsMixingVertex, binsMixingMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<TMcCollisions, TTracksTrig, TTracksAssoc, BinningTypeMcTruth> const pair{binningWithTracksSize, configTask.nMixedEvents, -1, mcCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      const auto multiplicity = collision1.multMCPVz();

      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
      fillCorrelations(corrContainer, CorrelationContainer::CFStep::kCFStepAll, tracks1, tracks2, multiplicity, collision1.posZ(), false);
    }
  }

  // ===================================================================================================================================================================================================================================================================
  // ===================================================================================================================================================================================================================================================================
  //    SAME EVENT PROCESS FUNCTIONS
  // ===================================================================================================================================================================================================================================================================
  // ===================================================================================================================================================================================================================================================================

  // ===================================================================================================================================================================================================================================================================
  //    DATA
  // ===================================================================================================================================================================================================================================================================

  // =====================================
  //    DATA : process same event correlations: TPC-TPC h-h case
  // =====================================

  void processSameTpcTpcChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             FilteredTracksWDcaSel const& tracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    //  the event histograms below are only filled for h-h case
    //  because there is a possibility of double-filling if more correlation
    //  options are ran at the same time
    //  temporary solution, since other correlation options always have to be ran with h-h, too
    //  TODO: rewrite it in a more intelligent way
    const auto multiplicity = getMultiplicityEstimator(collision, true);
    // registry.fill(HIST(Form("Data/hMultiplicity_%s", WhatMultiplicityEstimator[HfTaskFlow::configCollision.multiplicityEstimator].data())), multiplicity);

    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, tracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcChCh, "DATA : Process same-event correlations for TPC-TPC h-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-TPC HF-h case for D0
  // =====================================

  void processSameTpcTpcD0Ch(FilteredCollisionsWSelMult::iterator const& collision,
                             FilteredTracksWDcaSel const& tracks,
                             HfCandidatesSelD0 const& candidates)
  {
    auto fillEventSelectionPlots = true;

    // When doing reference flow, two cases are used (HF-h, h-h) and thus eventSelectionPlots was filled twice
    if (configTask.doReferenceFlow) {
      fillEventSelectionPlots = false;
    }

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, tracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcD0Ch, "DATA : Process same-event correlations for TPC-TPC D0-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-TPC HF-h case for Lc
  // =====================================

  void processSameTpcTpcLcCh(FilteredCollisionsWSelMult::iterator const& collision,
                             FilteredTracksWDcaSel const& tracks,
                             HfCandidatesSelLc const& candidates)
  {
    auto fillEventSelectionPlots = true;

    // When doing reference flow, two cases are used (HF-h, h-h) and thus eventSelectionPlots was filled twice
    if (configTask.doReferenceFlow) {
      fillEventSelectionPlots = false;
    }

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, tracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcLcCh, "DATA : Process same-event correlations for TPC-TPC Lc-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT h-h case
  // =====================================

  void processSameTpcMftChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             FilteredTracksWDcaSel const& tracks,
                             FilteredMftTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    // I use kCFStepAll for running my code with all MFTTracks were the reassociation process was not applied
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChCh, "DATA : Process same-event correlations for TPC-MFT h-h case", false);

  void processSameTpcMftChChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         FilteredTracksWDcaSel const& tracks,
                                         FilteredMftTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    registry.fill(HIST("Data/Mft/hNTracks"), tracks.size());
    registry.fill(HIST("Data/Mft/hNMftTracks"), mftTracks.size());
    registry.fill(HIST("Data/Mft/hNBestCollisionFwd"), reassociatedMftTracks.size());

    // const auto multiplicity = collision.multNTracksPV();
    const auto multiplicity = getMultiplicityEstimator(collision, true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChChReassociated, "DATA : Process same-event correlations for TPC-MFT h-h case reassociated", false);

  /*
  void processSameTpcMftChChReassociated3d(FilteredCollisionsWSelMult::iterator const& collision,
                                           soa::SmallGroups<aod::BestCollisionsFwd3d> const& reassociatedMftTracks,
                                           FilteredTracksWDcaSel const& tracks,
                                           FilteredMftTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    registry.fill(HIST("Data/Mft/hNTracks"), tracks.size());
    registry.fill(HIST("Data/Mft/hNMftTracks"), mftTracks.size());
    registry.fill(HIST("Data/Mft/hNBestCollisionFwd"), reassociatedMftTracks.size());

    // const auto multiplicity = collision.multNTracksPV();
    const auto multiplicity = getMultiplicityEstimator(collision, true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChChReassociated3d, "DATA : Process same-event correlations for TPC-MFT h-h case 3d reassociated", false);
  */

  void processSameTpcMftChChNonAmbiguous(FilteredCollisionsWSelMult::iterator const& collision,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         FilteredTracksWDcaSel const& tracks,
                                         FilteredMftTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    registry.fill(HIST("Data/Mft/hNTracks"), tracks.size());
    registry.fill(HIST("Data/Mft/hNMftTracks"), mftTracks.size());
    registry.fill(HIST("Data/Mft/hNBestCollisionFwd"), reassociatedMftTracks.size());

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    // I use kCFStepTracked for running my code with only non-ambiguous MFTTracks
    // This is the same as running with reassociatedMftTracks, but applying one more cut in the fillCorrelations function
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, reassociatedMftTracks, multiplicity, collision.posZ(), true, true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChChNonAmbiguous, "DATA : Process same-event correlations for TPC-MFT h-h case with non-ambiguous tracks", false);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT HF-h case for D0
  // =====================================

  void processSameTpcMftD0Ch(FilteredCollisionsWSelMult::iterator const& collision,
                             HfCandidatesSelD0 const& candidates,
                             FilteredTracksWDcaSel const& /*tracks*/,
                             FilteredMftTracks const& mftTracks)
  {
    auto fillEventSelectionPlots = true;

    // When doing reference flow, two cases are used (HF-h, h-h) and thus eventSelectionPlots was filled twice
    if (configTask.doReferenceFlow) {
      fillEventSelectionPlots = false;
    }

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftD0Ch, "DATA : Process same-event correlations for TPC-MFT D0-h case", false);

  void processSameTpcMftD0ChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         HfCandidatesSelD0 const& candidates,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         FilteredMftTracks const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftD0ChReassociated, "DATA : Process same-event correlations for TPC-MFT D0-h case reassociated", false);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT HF-h case for Lc
  // =====================================

  void processSameTpcMftLcCh(FilteredCollisionsWSelMult::iterator const& collision,
                             HfCandidatesSelLc const& candidates,
                             FilteredTracksWDcaSel const& /*tracks*/,
                             FilteredMftTracks const& mftTracks)
  {
    auto fillEventSelectionPlots = true;

    // When doing reference flow, two cases are used (HF-h, h-h) and thus eventSelectionPlots was filled twice
    if (configTask.doReferenceFlow) {
      fillEventSelectionPlots = false;
    }

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftLcCh, "DATA : Process same-event correlations for TPC-MFT Lc-h case", false);

  void processSameTpcMftLcChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         HfCandidatesSelLc const& candidates,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         FilteredMftTracks const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftLcChReassociated, "DATA : Process same-event correlations for TPC-MFT Lc-h case reassociated", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FV0A Ch. Part. - Ch. Part
  // =====================================

  void processSameTpcFv0aChCh(FilteredCollisionsWSelMult::iterator const& collision,
                              FilteredTracksWDcaSel const& tracks,
                              aod::FV0As const& fv0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, fv0, fv0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFv0aChCh, "DATA : Process same-event correlations for TPC-FV0-A h-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FV0A D0 - Ch. Part
  // =====================================

  void processSameTpcFv0aD0Ch(FilteredCollisionsWSelMult::iterator const& collision,
                              HfCandidatesSelD0 const& candidates,
                              aod::FV0As const& fv0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, fv0, fv0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFv0aD0Ch, "DATA : Process same-event correlations for TPC-FV0-A D0-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FV0A Lc - Ch. Part
  // =====================================

  void processSameTpcFv0aLcCh(FilteredCollisionsWSelMult::iterator const& collision,
                              HfCandidatesSelLc const& candidates,
                              aod::FV0As const& fv0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, fv0, fv0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFv0aLcCh, "DATA : Process same-event correlations for TPC-FV0-A Lc-h case", false);

  // =====================================
  //    DATA : process same event correlations: MFT-FV0A Ch. Part. - Ch. Part
  // =====================================

  void processSameMftFv0aChCh(FilteredCollisionsWSelMult::iterator const& collision,
                              FilteredMftTracks const& mftTracks,
                              aod::FV0As const& fv0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, mftTracks, fv0, fv0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFv0aChCh, "DATA : Process same-event correlations for MFT-FV0-A h-h case", false);

  void processSameMftFv0aChChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                          soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                          FilteredMftTracks const&,
                                          aod::FV0As const& fv0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFITReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, reassociatedMftTracks, fv0, fv0as, multiplicity, collision.posZ(), true, false);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFv0aChChReassociated, "DATA : Process same-event correlations for MFT-FV0a h-h case reassociated", false);

  /*
  void processSameMftFv0aChChReassociated3d(FilteredCollisionsWSelMult::iterator const& collision,
                                            soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                            FilteredMftTracks const&,
                                            aod::FV0As const& fv0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFITReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, reassociatedMftTracks, fv0, fv0as, multiplicity, collision.posZ(), true, false);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFv0aChChReassociated3d, "DATA : Process same-event correlations for MFT-FV0a h-h case 3d reassociated", false);
  */

  void processSameMftFv0aChChNonAmbiguous(FilteredCollisionsWSelMult::iterator const& collision,
                                          soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                          FilteredMftTracks const&,
                                          aod::FV0As const& fv0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFITReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, reassociatedMftTracks, fv0, fv0as, multiplicity, collision.posZ(), true, true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFv0aChChNonAmbiguous, "DATA : Process same-event correlations for MFT-FV0a h-h non-ambiguous case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FT0A Ch. Part. - Ch. Part
  // =====================================

  void processSameTpcFt0aChCh(FilteredCollisionsWSelMult::iterator const& collision,
                              FilteredTracksWDcaSel const& tracks,
                              aod::FT0s const& ft0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, ft0, ft0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFt0aChCh, "DATA : Process same-event correlations for TPC-FT0-A h-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FT0A Ch. Part. - Ch. Part
  // =====================================

  void processSameTpcFt0aD0Ch(FilteredCollisionsWSelMult::iterator const& collision,
                              HfCandidatesSelD0 const& candidates,
                              aod::FT0s const& ft0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();

      sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, ft0, ft0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFt0aD0Ch, "DATA : Process same-event correlations for TPC-FT0-A D0-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FT0A Ch. Part. - Ch. Part
  // =====================================

  void processSameTpcFt0aLcCh(FilteredCollisionsWSelMult::iterator const& collision,
                              HfCandidatesSelLc const& candidates,
                              aod::FT0s const& ft0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();

      sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, ft0, ft0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFt0aLcCh, "DATA : Process same-event correlations for TPC-FT0-A Lc-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FT0A Ch. Part. - Ch. Part
  // =====================================

  void processSameMftFt0aChCh(FilteredCollisionsWSelMult::iterator const& collision,
                              FilteredMftTracks const& mftTracks,
                              aod::FT0s const& ft0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFIT(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, mftTracks, ft0, ft0as, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFt0aChCh, "DATA : Process same-event correlations for MFT-FT0-A h-h case", false);

  void processSameMftFt0aChChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                          soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                          FilteredMftTracks const&,
                                          aod::FT0s const& ft0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFITReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, reassociatedMftTracks, ft0, ft0as, multiplicity, collision.posZ(), true, false);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFt0aChChReassociated, "DATA : Process same-event correlations for MFT-FT0-A h-h case reassociated", false);

  void processSameMftFt0aChChNonAmbiguous(FilteredCollisionsWSelMult::iterator const& collision,
                                          soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                          FilteredMftTracks const&,
                                          aod::FT0s const& ft0as)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = getMultiplicityEstimator(collision, true);

    if (collision.has_foundFT0()) {
      const auto& ft0 = collision.foundFT0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFITReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, reassociatedMftTracks, ft0, ft0as, multiplicity, collision.posZ(), true, true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFt0aChChNonAmbiguous, "DATA : Process same-event correlations for MFT-FT0-A h-h case non ambiguous", false);

  // ===================================================================================================================================================================================================================================================================
  //    MONTE-CARLO
  // ===================================================================================================================================================================================================================================================================

  // =====================================
  //    MONTE-CARLO GENERATED : process same event correlations : TPC-MFT D0-ch. part. case
  // =====================================

  void processSameTpcMftD0ChMcGen(FilteredMcCollisions::iterator const& mcCollision,
                                  McParticles2ProngMatched const& mcParticles2Prong,
                                  McParticles const& mcParticles)
  {
    const auto multiplicity = mcCollision.multMCPVz();

    BinningPolicyBase<2> const baseBinning{{axisVertex, axisMultiplicity}, true};

    sameEventHfMc->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations(sameEventHfMc, CorrelationContainer::CFStep::kCFStepAll, mcParticles2Prong, mcParticles, multiplicity, mcCollision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftD0ChMcGen, "MONTE-CARLO : Process same-event correlations for TPC-MFT D0-h case", false);

  // =====================================
  //    MONTE-CARLO GENERATED : process same event correlations : TPC-MFT Lc-ch. part. case
  // =====================================

  void processSameTpcMftLcChMcGen(FilteredMcCollisions::iterator const& mcCollision,
                                  McParticles3ProngMatched const& mcParticles3Prong,
                                  McParticles const& mcParticles)
  {
    const auto multiplicity = mcCollision.multMCPVz();

    BinningPolicyBase<2> const baseBinning{{axisVertex, axisMultiplicity}, true};

    sameEventHfMc->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations(sameEventHfMc, CorrelationContainer::CFStep::kCFStepAll, mcParticles3Prong, mcParticles, multiplicity, mcCollision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftLcChMcGen, "MONTE-CARLO : Process same-event correlations for TPC-MFT Lc-h case", false);

  // ===================================================================================================================================================================================================================================================================
  // ===================================================================================================================================================================================================================================================================
  //    MIXED EVENT PROCESS FUNCTIONS
  // ===================================================================================================================================================================================================================================================================
  // ===================================================================================================================================================================================================================================================================

  // ===================================================================================================================================================================================================================================================================
  //    DATA
  // ===================================================================================================================================================================================================================================================================

  // =====================================
  //    DATA : process mixed event correlations:TPC-TPC h-h case
  // =====================================

  void processMixedTpcTpcChCh(FilteredCollisionsWSelMult const& collisions,
                              FilteredTracksWDcaSel const& tracks)
  {
    mixCollisions(collisions, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, tracks, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChCh, "DATA : Process mixed-event correlations for TPC-TPC h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-TPC HF-h case for D0
  // =====================================

  void processMixedTpcTpcD0Ch(FilteredCollisionsWSelMult const& collisions,
                              FilteredTracksWDcaSel const& tracks,
                              HfCandidatesSelD0 const& candidates)
  {
    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, tracks, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcD0Ch, "DATA : Process mixed-event correlations for TPC-TPC D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-TPC HF-h case for Lc
  // =====================================

  void processMixedTpcTpcLcCh(FilteredCollisionsWSelMult const& collisions,
                              FilteredTracksWDcaSel const& tracks,
                              HfCandidatesSelLc const& candidates)
  {
    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, tracks, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcLcCh, "DATA : Process mixed-event correlations for TPC-TPC Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT h-h case
  // =====================================

  void processMixedTpcMftChCh(FilteredCollisionsWSelMult const& collisions,
                              FilteredTracksWDcaSel const& tracks,
                              FilteredMftTracks const& mftTracks)
  {
    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, tracks, mftTracks, mixedEvent);
    // mixCollisions(collisions, CorrelationContainer::kCFStepAll, tracks, mftTracks, getMultiplicity, mixedEvent);

    // The next following two lines were supposed to be used to do mixed event with the reassociated MFT tracks
    // However it seems the O2physics framework cannot handle how these combinations requests grouping according to Anton Alkin
    // So I leave them commented for now until it is solved, and put the "normal" mixCollisions back with kCFStepReconstructed

    // mixCollisionsReassociatedMftTracks(collisions, CorrelationContainer::kCFStepReconstructed, tracks, reassociatedMftTracks, getMultiplicity, mixedEvent, false);

    // mixCollisionsReassociatedMftTracks(collisions, CorrelationContainer::kCFStepTracked, tracks, reassociatedMftTracks, getMultiplicity, mixedEvent, true);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftChCh, "DATA : Process mixed-event correlations for TPC-MFT h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT HF-h case for D0
  // =====================================

  void processMixedTpcMftD0Ch(FilteredCollisionsWSelMult const& collisions,
                              HfCandidatesSelD0 const& candidates,
                              FilteredMftTracks const& mftTracks,
                              FilteredTracksWDcaSel const& /*tracks*/)
  {
    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, mftTracks, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftD0Ch, "DATA : Process mixed-event correlations for TPC-MFT D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT HF-h case
  // =====================================

  void processMixedTpcMftLcCh(FilteredCollisionsWSelMult const& collisions,
                              HfCandidatesSelLc const& candidates,
                              FilteredMftTracks const& mftTracks)
  {
    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, mftTracks, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftLcCh, "DATA : Process mixed-event correlations for TPC-MFT Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A ch part. - ch. part. case
  // =====================================

  void processMixedTpcFv0aChCh(FilteredCollisionsWSelMult const& collisions,
                               FilteredTracksWDcaSel const& tracks,
                               aod::FV0As const& fv0as)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, tracks, fv0as, perColTracks, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFv0aChCh, "DATA : Process mixed-event correlations for TPC-FV0-A h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A D0 - ch. part. case
  // =====================================

  void processMixedTpcFv0aD0Ch(FilteredCollisionsWSelMult const& collisions,
                               HfCandidatesSelD0 const& candidates,
                               aod::FV0As const& fv0as)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, candidates, fv0as, perColD0s, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFv0aD0Ch, "DATA : Process mixed-event correlations for TPC-FV0-A D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A Lc - ch. part. case
  // =====================================

  void processMixedTpcFv0aLcCh(FilteredCollisionsWSelMult const& collisions,
                               HfCandidatesSelLc const& candidates,
                               aod::FV0As const& fv0as)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, candidates, fv0as, perColLcs, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFv0aLcCh, "DATA : Process mixed-event correlations for TPC-FV0-A Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A ch part. - ch. part. case
  // =====================================

  void processMixedMftFv0aChCh(FilteredCollisionsWSelMult const& collisions,
                               FilteredMftTracks const& mftTracks,
                               aod::FV0As const& fv0as)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, mftTracks, fv0as, perColMftTracks, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedMftFv0aChCh, "DATA : Process mixed-event correlations for Mft-FV0-A h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FT0-A ch part. - ch. part. case
  // =====================================

  void processMixedTpcFt0aChCh(FilteredCollisionsWSelMult const& collisions,
                               FilteredTracksWDcaSel const& tracks,
                               aod::FT0s const& ft0s)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, tracks, ft0s, perColTracks, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFt0aChCh, "DATA : Process mixed-event correlations for TPC-FT0-A h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FT0-A D0 - ch. part. case
  // =====================================

  void processMixedTpcFt0aD0Ch(FilteredCollisionsWSelMult const& collisions,
                               HfCandidatesSelD0 const& candidates,
                               aod::FT0s const& ft0s)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, candidates, ft0s, perColD0s, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFt0aD0Ch, "DATA : Process mixed-event correlations for TPC-FT0-A D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FT0-A Lc - ch. part. case
  // =====================================

  void processMixedTpcFt0aLcCh(FilteredCollisionsWSelMult const& collisions,
                               HfCandidatesSelLc const& candidates,
                               aod::FT0s const& ft0s)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, candidates, ft0s, perColLcs, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFt0aLcCh, "DATA : Process mixed-event correlations for TPC-FT0-A Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FT0-A ch part. - ch. part. case
  // =====================================

  void processMixedMftFt0aChCh(FilteredCollisionsWSelMult const& collisions,
                               FilteredMftTracks const& mftTracks,
                               aod::FT0s const& ft0s)
  {
    mixCollisionsFIT(collisions, CorrelationContainer::kCFStepReconstructed, mftTracks, ft0s, perColMftTracks, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedMftFt0aChCh, "DATA : Process mixed-event correlations for MFT-FT0-A h-h case", false);

  // ===================================================================================================================================================================================================================================================================
  //    MONTE-CARLO
  // ===================================================================================================================================================================================================================================================================

  // =====================================
  //    MONTE-CARLO GENERATED : process mixed event correlations: TPC-MFT D0-ch. part. case
  // =====================================

  void processMixedTpcMftD0ChMcGen(FilteredMcCollisions const& mcCollisions,
                                   McParticles2ProngMatched const& mcParticles2Prong,
                                   McParticles const& mcParticles)
  {
    auto getMultiplicity = [](FilteredMcCollisions::iterator const& mcCollision) {
      auto multiplicity = mcCollision.multMCPVz();
      return multiplicity;
    };

    mixCollisionsMcTruth(mcCollisions, mcParticles2Prong, mcParticles, getMultiplicity, mixedEventHfMc);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftD0ChMcGen, "MONTE-CARLO : Process mixed-event correlations for TPC-MFT D0-h case", false);

  // =====================================
  //    MONTE-CARLO GENERATED : process mixed event correlations: TPC-MFT Lc-ch. part. case
  // =====================================

  void processMixedTpcMftLcChMcGen(FilteredMcCollisions const& mcCollisions,
                                   McParticles3ProngMatched const& mcParticles3Prong,
                                   McParticles const& mcParticles)
  {
    auto getMultiplicity = [](FilteredMcCollisions::iterator const& mcCollision) {
      auto multiplicity = mcCollision.multMCPVz();
      return multiplicity;
    };

    mixCollisionsMcTruth(mcCollisions, mcParticles3Prong, mcParticles, getMultiplicity, mixedEventHfMc);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftLcChMcGen, "MONTE-CARLO : Process mixed-event correlations for TPC-MFT D0-h case", false);

  // ===================================================================================================================================================================================================================================================================
  // ===================================================================================================================================================================================================================================================================
  //    EFFICIENCIES PROCESS FUNCTIONS
  // ===================================================================================================================================================================================================================================================================
  // ===================================================================================================================================================================================================================================================================

  // NOTE SmallGroups includes soa::Filtered always -> in the smallGroups there is the equivalent of FilteredCollisionsWSelMultMcLabels
  void processMcEfficiencyMft(FilteredMcCollisions::iterator const& mcCollision,
                              McParticles const& mcParticles,
                              soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults>> const& collisionsMcLabels,
                              MftTracksMcLabels const& mftTTracksMcLabels)
  {
    LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisionsMcLabels.size());

    auto multiplicity = mcCollision.multMCPVz();
    if (configTask.centralityBinsForMc) {
      if (collisionsMcLabels.size() == 0) {
        return;
      }
      for (const auto& collision : collisionsMcLabels) {
        multiplicity = collision.multNTracksPV();
      }
    }

    // Primaries
    for (const auto& mcParticle : mcParticles) {
      if (!isAcceptedMftMcParticle(mcParticle)) {
        sameEventHf->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
      }
    }
    for (const auto& collision : collisionsMcLabels) {
      auto groupedMftTTracksMcLabels = mftTTracksMcLabels.sliceBy(mftTracksPerCollision, collision.globalIndex());
      LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());
      LOGF(info, "  which has %d mft tracks", groupedMftTTracksMcLabels.size());

      for (const auto& mftTrack : groupedMftTTracksMcLabels) {
        if (mftTrack.has_mcParticle()) {
          const auto& mcParticle = mftTrack.mcParticle();
          if (!isAcceptedMftMcParticle(mcParticle)) {
            sameEventHf->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          }
          sameEventHf->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
        } else {
          // fake track
          // In the MFT the measurement of pT is not precise
          sameEventHf->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake, mftTrack.eta(), mftTrack.pt(), 0, multiplicity, mcCollision.posZ());
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processMcEfficiencyMft, "MONTE-CARLO : Extract efficiencies for MFT tracks", false);

}; // End of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskFlow>(cfgc)};
}
