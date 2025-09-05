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
#include "PWGCF/Core/PairCuts.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsPid.h"
#include "PWGMM/Mult/DataModel/Index.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FV0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include <TComplex.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <THn.h>
#include <TMath.h>

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::pid_tpc_tof_utils;
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
};

enum CorrelatedParticles {
  ChPartChPart,
  D0ChPart,
  LcChPart
};

static constexpr std::string_view WhatDataType[] = {"Data/", "MC/"};
// static constexpr std::string_view whatEventType[] = {"SameEvent/", "MixedEvent/"};
static constexpr std::string_view WhatCorrelationCase[] = {"TpcTpc/", "TpcMft/", "TpcFv0a/", "MftFv0a/"};
static constexpr std::string_view WhatParticles[] = {"ChPartChPart/", "D0ChPart/", "LcChPart/"};

// static constexpr float kPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};

struct HfTaskFlow {

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  //  configurables for processing options

  Configurable<bool> centralityBinsForMc{"centralityBinsForMc", false, "false = OFF, true = ON for data like multiplicity/centrality bins for MC steps"};
  Configurable<float> mftMaxDCAxy{"mftMaxDCAxy", 2.0f, "Cut on dcaXY for MFT tracks"};
  Configurable<bool> doHeavyFlavor{"doHeavyFlavor", false, "Flag to know we in the heavy flavor case or not"};
  Configurable<bool> doReferenceFlow{"doReferenceFlow", false, "Flag to know if reference flow should be done"};
  Configurable<bool> isReadoutCenter{"isReadoutCenter", false, "Enable Readout Center"};
  // Configurable<float> doTwoTrackCut{"doTwoTrackCut", -1, "Two track cut: -1 = off; >0 otherwise distance value (suggested: 0.02)"};
  Configurable<bool> processRun2{"processRun2", false, "Flag to run on Run 2 data"};
  Configurable<bool> processRun3{"processRun3", true, "Flag to run on Run 3 data"};
  Configurable<bool> processMc{"processMc", false, "Flag to run on MC"};
  Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per event"};
  // Configurable<float> twoTrackCutMinRadius{"twoTrackCutMinRadius", 0.8f, "Two track cut : radius in m from which two tracks cuts are applied"};
  //   configurables for collisions
  Configurable<float> zVertexMax{"zVertexMax", 7.0f, "Accepted z-vertex range"};
  //  configurables for TPC tracks
  Configurable<float> etaTpcTrackMax{"etaTpcTrackMax", 0.8f, "max. eta of TPC tracks"};
  Configurable<float> ptTpcTrackMin{"ptTpcTrackMin", 0.5f, "min. pT of TPC tracks"};
  //  configurables for HF candidates
  Configurable<float> etaCandidateMax{"etaCandidateMax", 0.8f, "max. eta of HF candidate"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for Hf candidates"};
  // Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  // Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  // Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for LambdaC"};
  // Configurable<int> selectionFlagLcToPiKP{"selectionFlagLcToPiKP", 1, "Selection Flag for LambdaC bar"};
  //   configurables for MFT tracks
  Configurable<float> etaMftTrackMax{"etaMftTrackMax", -2.4f, "Maximum value for the eta of MFT tracks"};
  Configurable<float> etaMftTrackMin{"etaMftTrackMin", -3.36f, "Minimum value for the eta of MFT tracks"};
  Configurable<std::vector<int>> mcTriggerPdgs{"mcTriggerPdgs", {421, -421}, "MC PDG codes to use exclusively as trigger particles. D0= +-421, Lc = +-4122"};
  Configurable<int> nClustersMftTrack{"nClustersMftTrack", 5, "Minimum number of clusters for the reconstruction of MFT tracks"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};

  HfHelper hfHelper;
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::vector<o2::detectors::AlignParam>* offsetFV0;
  o2::ccdb::CcdbApi ccdbApi;
  o2::fv0::Geometry* fv0Det;
  std::vector<int> hfIndexCache;

  // =========================
  //      using declarations : DATA
  // =========================

  using FilteredCollisionsWSelMult = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using HfCandidatesSelD0 = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using HfCandidatesSelLc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;

  // using FilteredMftTracks = soa::Filtered<aod::MFTTracks>;
  //  using FilteredMftTracksWColls = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>>;
  // using FilteredAndReassociatedMftTracks = soa::Filtered<soa::Join<aod::BestCollisionsFwd, aod::MFTTracks>>;

  using FilteredTracksWDcaSel = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection>>;

  // =========================
  //      using declarations : MONTE CARLO
  // =========================

  // Even add McCollisions in the join ?
  // Kata adds subscribes to it but do not add it in the join
  // using FilteredCollisionsWSelMultMcLabels = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::McCollisions>>;

  using FilteredCollisionsWSelMultMcLabels = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults>>;
  using FilteredMcCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::MultMCExtras, aod::McCollsExtra>>;
  using HfCandidatesSelD0McRec = soa::Join<HfCandidatesSelD0, aod::HfCand2ProngMcRec>;
  using HfCandidatesSelLcMcRec = soa::Join<HfCandidatesSelLc, aod::HfCand3ProngMcRec>;
  using McParticles = aod::McParticles;
  using McParticles2ProngMatched = soa::Join<McParticles, aod::HfCand2ProngMcGen>;
  using McParticles3ProngMatched = soa::Join<McParticles, aod::HfCand3ProngMcGen>;
  // using FilteredMftTracksWCollsMcLabels = soa::Filtered<soa::Join<aod::MFTTracks, aod::MFTTrkCompColls, aod::McMFTTrackLabels>>;
  using MftTracksMcLabels = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using FilteredTracksWDcaSelMC = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::McTrackLabels>>;

  // =========================
  //      Filters & partitions : DATA
  // =========================

  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilterD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagHf) ||
                             (aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagHf);

  Filter candidateFilterLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagHf) ||
                             (aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagHf);

  //  Collision filters
  //  FIXME: The filter is applied also on the candidates! Beware!
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < zVertexMax;

  Filter trackFilter = (nabs(aod::track::eta) < etaTpcTrackMax) &&
                       (aod::track::pt > ptTpcTrackMin) &&
                       requireGlobalTrackWoPtEtaInFilter();

  // Filter mftTrackEtaFilter = (aod::fwdtrack::eta < etaMftTrackMax) &&
  //                            (aod::fwdtrack::eta > etaMftTrackMin);

  // Filter mftTrackHasCollision = aod::fwdtrack::collisionId > 0;

  // Filters below will be used for uncertainties
  // Filter mftTrackCollisionIdFilter = (aod::fwdtrack::bestCollisionId >= 0);
  // Filter mftTrackDcaFilter = (nabs(aod::fwdtrack::bestDCAXY) < mftMaxDCAxy);

  // =========================
  //      Filters & partitions : MC
  // =========================

  Filter candidateFilterD0Mc = (aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf) ||
                               (aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf);

  Filter candidateFilterLcMc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagHf) ||
                               (aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagHf);

  // From Katarina's code, but not sure if I use it
  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < zVertexMax;

  // Filter mcParticlesFilter = (nabs(aod::mcparticle::eta) < etaTpcTrackMax) &&
  //                            (aod::mcparticle::pt > ptTpcTrackMin);

  // I didn't manage to make partitions work with my mixed event, as I am pair my tracks BEFORE looping over collisions
  // I am thus not able to group tracks with sliceBy and can't use this method
  // For now I am fine as I am doing only TPC-MFT correlations and using only McParticles with MFT acceptance
  // However at some point I will have to use tracks from the other side (FV0, FT0-A) and I will have to do something about it
  // TO-DO : either change how I do mixed event, or implement isAcceptedTpcMcParticle, isAcceptedMftMcParticle
  // Partition<aod::McParticles> mcParticlesMft = (aod::mcparticle::eta > etaMftTrackMin) && (aod::mcparticle::eta < etaMftTrackMax);
  // Partition<aod::McParticles> mcParticlesTpc = (nabs(aod::mcparticle::eta) < etaTpcTrackMax) &&
  //                                             (aod::mcparticle::pt > ptTpcTrackMin);

  // =========================
  //      Preslice : DATA
  // =========================

  Preslice<aod::MFTTracks> perColMftTracks = o2::aod::fwdtrack::collisionId;
  Preslice<FilteredTracksWDcaSel> perColTracks = aod::track::collisionId;
  Preslice<HfCandidatesSelD0> perColD0s = aod::track::collisionId;
  Preslice<HfCandidatesSelLc> perColLcs = aod::track::collisionId;

  // =========================
  //      Preslice : MC
  // =========================

  Preslice<MftTracksMcLabels> mftTracksPerCollision = aod::fwdtrack::collisionId;
  // Preslice<HfCandidatesSelD0McRec> d0CandidatesPerCollision = aod::hf_cand::collisionId;
  // Preslice<McParticles> mcPerCol = aod::mcparticle::mcCollisionId;
  // PresliceUnsorted<FilteredCollisionsWSelMultMcLabels> collisionsMcLabelPerMcCollision = aod::mccollisionlabel::mcCollisionId;

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

  template <DataType dataType, CorrelationCase correlationCase, CorrelatedParticles correlatedParticles>
  void addHistograms()
  {
    registry.add(Form("%s%s%shEtaTrigger", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH1D, {axisEtaTrigger}});
    registry.add(Form("%s%s%shPhiTrigger", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH1D, {axisPhi}});
    registry.add(Form("%s%s%shPtTrigger", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH1D, {axisPt}});
    registry.add(Form("%s%s%shYieldsTrigger", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTrigger}});
    registry.add(Form("%s%s%shEtaPhiTrigger", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH3F, {axisMultiplicity, axisEtaTrigger, axisPhi}});
    registry.add(Form("%s%s%shEtaAssociated", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH1D, {axisEtaAssociated}});
    registry.add(Form("%s%s%shPhiAssociated", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH1D, {axisPhi}});
    registry.add(Form("%s%s%shEtaPhiAssociated", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH3F, {axisMultiplicity, axisEtaAssociated, axisPhi}});
    registry.add(Form("%s%s%shMultiplicity", WhatDataType[dataType].data(), WhatCorrelationCase[correlationCase].data(), WhatParticles[correlatedParticles].data()), "", {HistType::kTH1D, {axisMultiplicity}});
  }

  //  =========================
  //      init()
  //  =========================
  void init(InitContext&)
  {
    // const int nBinsMix = axisMultiplicity->size() * axisVertex->size();
    ccdb->setURL(ccdbUrl);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", noLaterThan.value);
    LOGF(info, "Offset for FV0-left: x = %.3f y = %.3f z = %.3f\n", (*offsetFV0)[0].getX(), (*offsetFV0)[0].getY(), (*offsetFV0)[0].getZ());
    LOGF(info, "Offset for FV0-right: x = %.3f y = %.3f z = %.3f\n", (*offsetFV0)[1].getX(), (*offsetFV0)[1].getY(), (*offsetFV0)[1].getZ());

    //  =========================
    //      Event histograms
    //  =========================

    registry.add("Data/hEventCounter", "hEventCounter", {HistType::kTH1D, {{EventSelectionStep::NEventSelectionSteps, -0.5, +EventSelectionStep::NEventSelectionSteps - 0.5}}});
    std::string labels[EventSelectionStep::NEventSelectionSteps];
    labels[EventSelectionStep::AllEvents] = "all";
    labels[EventSelectionStep::AfterEventSelection] = "after Physics selection";
    registry.get<TH1>(HIST("Data/hEventCounter"))->SetMinimum(0);

    for (int iBin = 0; iBin < EventSelectionStep::NEventSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("Data/hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    registry.add("Data/TpcMft/hAmbiguityOfMftTracks", "hAmbiguityOfMftTracks", {HistType::kTH1D, {{MftTrackAmbiguityStep::NMftAmbiguitySteps, -0.5, +MftTrackAmbiguityStep::NMftAmbiguitySteps - 0.5}}});
    std::string labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NMftAmbiguitySteps];
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::AllMftTracks] = "all MFT tracks";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::AfterTrackSelection] = "MFT tracks after selection";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NumberOfAmbiguousTracks] = "how much tracks are ambigous";
    labelsAmbiguityOfMftTracks[MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks] = "how much tracks are non-ambiguous";
    registry.get<TH1>(HIST("Data/TpcMft/hAmbiguityOfMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftTrackAmbiguityStep::NMftAmbiguitySteps; iBin++) {
      registry.get<TH1>(HIST("Data/TpcMft/hAmbiguityOfMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsAmbiguityOfMftTracks[iBin].data());
    }

    registry.add("Data/TpcMft/hMftTracksSelection", "hMftTracksSelection", {HistType::kTH1D, {{MftTrackSelectionStep::NMftTrackSelectionSteps, -0.5, +MftTrackSelectionStep::NMftTrackSelectionSteps - 0.5}}});
    std::string labelsMftTracksSelection[MftTrackSelectionStep::NMftTrackSelectionSteps];
    labelsMftTracksSelection[MftTrackSelectionStep::NoSelection] = "all MFT tracks";
    labelsMftTracksSelection[MftTrackSelectionStep::Eta] = "MFT tracks after eta selection";
    labelsMftTracksSelection[MftTrackSelectionStep::Cluster] = "MFT tracks after clusters selection";
    registry.get<TH1>(HIST("Data/TpcMft/hMftTracksSelection"))->SetMinimum(0);

    for (int iBin = 0; iBin < MftTrackSelectionStep::NMftTrackSelectionSteps; iBin++) {
      registry.get<TH1>(HIST("Data/TpcMft/hMftTracksSelection"))->GetXaxis()->SetBinLabel(iBin + 1, labelsMftTracksSelection[iBin].data());
    }

    registry.add("Data/TpcMft/hReassociationMftTracks", "hReassociationMftTracks", {HistType::kTH1D, {{ReassociationMftTracks::NReassociationMftTracksSteps, -0.5, +ReassociationMftTracks::NReassociationMftTracksSteps - 0.5}}});
    std::string labelsReassociationMftTracks[ReassociationMftTracks::NReassociationMftTracksSteps];
    labelsReassociationMftTracks[ReassociationMftTracks::NotReassociatedMftTracks] = "MFT tracks after track selection";
    labelsReassociationMftTracks[ReassociationMftTracks::ReassociatedMftTracks] = "Reassociated MFT tracks by DCAxy method";
    registry.get<TH1>(HIST("Data/TpcMft/hReassociationMftTracks"))->SetMinimum(0);

    for (int iBin = 0; iBin < ReassociationMftTracks::NReassociationMftTracksSteps; iBin++) {
      registry.get<TH1>(HIST("Data/TpcMft/hReassociationMftTracks"))->GetXaxis()->SetBinLabel(iBin + 1, labelsReassociationMftTracks[iBin].data());
    }

    registry.add("Data/TpcMft/hNTracks", "", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/TpcMft/hNMftTracks", "", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/TpcMft/hNBestCollisionFwd", "", {HistType::kTH1F, {axisMultiplicity}});

    //  =========================
    //      Declaration of correlation containers and their respective axis
    //  =========================

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis = {{axisMass, "m_{inv} (GeV/c^{2})"}};

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

    if (doprocessSameTpcMftChCh || doprocessSameTpcMftChChReassociated || doprocessSameTpcMftChChNonAmbiguous) {
      addHistograms<Data, TpcMft, ChPartChPart>();

      // All MFT tracks
      registry.add("Data/TpcMft/kCFStepAll/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/TpcMft/kCFStepAll/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      // Only non-ambiguous MFT tracks
      registry.add("Data/TpcMft/kCFStepTracked/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/TpcMft/kCFStepTracked/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
      mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    }

    if (doprocessSameTpcMftD0Ch || doprocessSameTpcMftD0ChReassociated) {
      addHistograms<Data, TpcMft, D0ChPart>();

      // All MFT tracks
      registry.add("Data/TpcMft/kCFStepAll/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/TpcMft/kCFStepAll/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      // Only non-ambiguous MFT tracks
      registry.add("Data/TpcMft/kCFStepTracked/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/TpcMft/kCFStepTracked/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    if (doprocessSameTpcMftLcCh || doprocessSameTpcMftLcChReassociated) {
      addHistograms<Data, TpcMft, LcChPart>();

      // All MFT tracks
      registry.add("Data/TpcMft/kCFStepAll/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/TpcMft/kCFStepAll/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      // Only non-ambiguous MFT tracks
      registry.add("Data/TpcMft/kCFStepTracked/hEta", "eta", {HistType::kTH1D, {axisEtaAssociated}});
      registry.add("Data/TpcMft/kCFStepTracked/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

      sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
      mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));
    }

    //  =========================
    //  Initialization of histograms and CorrelationContainers for TpcFv0a cases
    //  =========================

    if (doprocessSameTpcFv0aChCh || doprocessSameTpcFv0aChCh) {
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

    if (doprocessSameMftFv0aChCh || doprocessSameMftFv0aChChReassociated || doprocessSameMftFv0aChChNonAmbiguous) {
      addHistograms<Data, MftFv0a, ChPartChPart>();

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
  //      Quality Assessment plots for TpcTpc cases
  // =========================

  // ---- DATA : TPC-TPC h-h Same Event QA ----
  template <typename TTrack>
  void fillTpcTpcChChSameEventQa(float multiplicity, TTrack const& track)
  {
    registry.fill(HIST("Data/TpcTpc/ChPartChPart/hPtTrigger"), track.pt());
    registry.fill(HIST("Data/TpcTpc/ChPartChPart/hEtaTrigger"), track.eta());
    registry.fill(HIST("Data/TpcTpc/ChPartChPart/hPhiTrigger"), track.phi());
    registry.fill(HIST("Data/TpcTpc/ChPartChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
    registry.fill(HIST("Data/TpcTpc/ChPartChPart/hEtaPhiTrigger"), multiplicity, track.eta(), track.phi());
    registry.fill(HIST("Data/TpcTpc/ChPartChPart/hMultiplicity"), multiplicity);
  }

  // ---- DATA : TPC-MFT HF-h Same Event candidates QA ----
  template <typename TTrack>
  void fillTpcTpcHfChSameEventCandidateQa(float multiplicity, TTrack const& track, bool isD0)
  {
    if (isD0) {
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hPtTrigger"), track.pt());
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hEtaTrigger"), track.eta());
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hPhiTrigger"), track.phi());
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hEtaPhiTrigger"), multiplicity, track.eta(), track.phi());
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hMultiplicity"), multiplicity);
    } else {
      registry.fill(HIST("Data/TpcTpc/LcChPart/hPtTrigger"), track.pt());
      registry.fill(HIST("Data/TpcTpc/LcChPart/hEtaTrigger"), track.eta());
      registry.fill(HIST("Data/TpcTpc/LcChPart/hPhiTrigger"), track.phi());
      registry.fill(HIST("Data/TpcTpc/LcChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcTpc/LcChPart/hEtaPhiTrigger"), multiplicity, track.eta(), track.phi());
      registry.fill(HIST("Data/TpcTpc/LcChPart/hMultiplicity"), multiplicity);
    }
  }

  // ---- DATA : TPC-TPC HF-h Same Event associated tracks QA ----
  template <typename TTrack>
  void fillTpcTpcHfChSameEventAssociatedQa(float multiplicity, TTrack const& track, bool isD0)
  {
    if (isD0) {
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hEtaAssociated"), track.eta());
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hPhiAssociated"), track.phi());
      registry.fill(HIST("Data/TpcTpc/D0ChPart/hEtaPhiAssociated"), multiplicity, track.eta(), track.phi());
    } else {
      registry.fill(HIST("Data/TpcTpc/LcChPart/hEtaAssociated"), track.eta());
      registry.fill(HIST("Data/TpcTpc/LcChPart/hPhiAssociated"), track.phi());
      registry.fill(HIST("Data/TpcTpc/LcChPart/hEtaPhiAssociated"), multiplicity, track.eta(), track.phi());
    }
  }

  // =========================
  //      Quality Assessment plots for TpcMft cases
  // =========================

  // ---- DATA : TPC-MFT h-h Same Event QA ----
  template <typename TTrack>
  void fillTpcMftChChSameEventQa(float multiplicity, TTrack const& track, bool isTPC)
  {
    float phi = track.phi();
    o2::math_utils::bringTo02Pi(phi);

    if (isTPC) { // trigger hadron from TPC
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hEtaTrigger"), track.eta());
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hPhiTrigger"), phi);
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hEtaPhiTrigger"), multiplicity, track.eta(), phi);
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hPtTrigger"), track.pt());
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hMultiplicity"), multiplicity);
    } else { // associated hadron from MFT
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hEtaAssociated"), track.eta());
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hPhiAssociated"), phi);
      registry.fill(HIST("Data/TpcMft/ChPartChPart/hEtaPhiAssociated"), multiplicity, track.eta(), phi);
    }
  }

  // ---- DATA : TPC-MFT HF-h Same Event associated tracks QA ----
  template <typename TTrack>
  void fillTpcMftHfChSameEventCandidateQa(float multiplicity, TTrack const& track, bool isD0)
  {
    if (isD0) {
      registry.fill(HIST("Data/TpcMft/D0ChPart/hPtTrigger"), track.pt());
      registry.fill(HIST("Data/TpcMft/D0ChPart/hEtaTrigger"), track.eta());
      registry.fill(HIST("Data/TpcMft/D0ChPart/hPhiTrigger"), track.phi());
      registry.fill(HIST("Data/TpcMft/D0ChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcMft/D0ChPart/hEtaPhiTrigger"), multiplicity, track.eta(), track.phi());
      registry.fill(HIST("Data/TpcMft/D0ChPart/hMultiplicity"), multiplicity);
    } else {
      registry.fill(HIST("Data/TpcMft/LcChPart/hPtTrigger"), track.pt());
      registry.fill(HIST("Data/TpcMft/LcChPart/hEtaTrigger"), track.eta());
      registry.fill(HIST("Data/TpcMft/LcChPart/hPhiTrigger"), track.phi());
      registry.fill(HIST("Data/TpcMft/LcChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcMft/LcChPart/hEtaPhiTrigger"), multiplicity, track.eta(), track.phi());
      registry.fill(HIST("Data/TpcMft/LcChPart/hMultiplicity"), multiplicity);
    }
  }

  // ---- DATA : TPC-MFT HF-h Same Event (Candidates) QA ----
  template <typename TTrack>
  void fillTpcMftHfChSameEventAssociatedQa(float multiplicity, TTrack const& track, bool isD0)
  {
    if (isD0) {
      registry.fill(HIST("Data/TpcMft/D0ChPart/hEtaAssociated"), track.eta());
      registry.fill(HIST("Data/TpcMft/D0ChPart/hPhiAssociated"), track.phi());
      registry.fill(HIST("Data/TpcMft/D0ChPart/hEtaPhiAssociated"), multiplicity, track.eta(), track.phi());
    } else {
      registry.fill(HIST("Data/TpcMft/LcChPart/hEtaAssociated"), track.eta());
      registry.fill(HIST("Data/TpcMft/LcChPart/hPhiAssociated"), track.phi());
      registry.fill(HIST("Data/TpcMft/LcChPart/hEtaPhiAssociated"), multiplicity, track.eta(), track.phi());
    }
  }

  // =========================
  //      Quality Assessment plots for TpcFv0 cases
  // =========================

  // ---- DATA : QA for FV0a ----
  void fillFv0aQa(float multiplicity, float const& eta, float const& phi, bool isTPC)
  {
    if (isTPC) {
      registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hEtaAssociated"), eta);
      registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hPhiAssociated"), phi);
      registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hEtaPhiAssociated"), multiplicity, eta, phi);
    } else {
      registry.fill(HIST("Data/MftFv0a/ChPartChPart/hEtaAssociated"), eta);
      registry.fill(HIST("Data/MftFv0a/ChPartChPart/hPhiAssociated"), phi);
      registry.fill(HIST("Data/MftFv0a/ChPartChPart/hEtaPhiAssociated"), multiplicity, eta, phi);
    }
  }

  // ---- DATA : TPC-MFT h-h Same Event trigger QA ----
  template <typename TTrack>
  void fillTpcFv0aChChSameEventTriggerQa(float multiplicity, TTrack const& track)
  {
    float phi = track.phi();
    o2::math_utils::bringTo02Pi(phi);
    registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hEtaTrigger"), track.eta());
    registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hPhiTrigger"), phi);
    registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hEtaPhiTrigger"), multiplicity, track.eta(), phi);
    registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hPtTrigger"), track.pt());
    registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
    registry.fill(HIST("Data/TpcFv0a/ChPartChPart/hMultiplicity"), multiplicity);
  }

  // ---- DATA : TPC-MFT HF-h Same Event associated tracks QA ----
  template <typename TTrack>
  void fillTpcFv0aHfChSameEventCandidateQa(float multiplicity, TTrack const& track, bool isD0)
  {
    if (isD0) {
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hPtTrigger"), track.pt());
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hEtaTrigger"), track.eta());
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hPhiTrigger"), track.phi());
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hEtaPhiTrigger"), multiplicity, track.eta(), track.phi());
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hMultiplicity"), multiplicity);
    } else {
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hPtTrigger"), track.pt());
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hEtaTrigger"), track.eta());
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hPhiTrigger"), track.phi());
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hEtaPhiTrigger"), multiplicity, track.eta(), track.phi());
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hMultiplicity"), multiplicity);
    }
  }

  // ---- DATA : TPC-MFT HF-h Same Event (Candidates) QA ----
  void fillTpcFv0aHfChSameEventAssociatedQa(float multiplicity, float const& eta, float const& phi, bool isD0)
  {
    if (isD0) {
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hEtaAssociated"), eta);
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hPhiAssociated"), phi);
      registry.fill(HIST("Data/TpcFv0a/D0ChPart/hEtaPhiAssociated"), multiplicity, eta, phi);
    } else {
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hEtaAssociated"), eta);
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hPhiAssociated"), phi);
      registry.fill(HIST("Data/TpcFv0a/LcChPart/hEtaPhiAssociated"), multiplicity, eta, phi);
    }
  }

  // =========================
  //      Quality Assessment plots for MftFv0 cases
  // =========================

  // ---- DATA : MFT-FV0a h-h Same Event QA for MFT ----
  template <typename TTrack>
  void fillMftFv0aChChSameEventQa(float multiplicity, TTrack const& track)
  {
    float phi = track.phi();
    o2::math_utils::bringTo02Pi(phi);

    registry.fill(HIST("Data/MftFv0a/ChPartChPart/hEtaTrigger"), track.eta());
    registry.fill(HIST("Data/MftFv0a/ChPartChPart/hPhiTrigger"), phi);
    registry.fill(HIST("Data/MftFv0a/ChPartChPart/hEtaPhiTrigger"), multiplicity, track.eta(), phi);
    registry.fill(HIST("Data/MftFv0a/ChPartChPart/hPtTrigger"), track.pt());
    registry.fill(HIST("Data/MftFv0a/ChPartChPart/hYieldsTrigger"), multiplicity, track.pt(), track.eta());
    registry.fill(HIST("Data/MftFv0a/ChPartChPart/hMultiplicity"), multiplicity);
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

  double getPhiFV0(unsigned int chno)
  {
    int cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);
    float offsetX, offsetY;
    if (isChnoInLeft) {
      offsetX = (*offsetFV0)[0].getX();
      offsetY = (*offsetFV0)[0].getY();
    } else {
      offsetX = (*offsetFV0)[1].getX();
      offsetY = (*offsetFV0)[1].getY();
    }

    o2::fv0::Point3Dsimple chPos;
    chPos = fv0Det->getReadoutCenter(chno);

    // if (isReadoutCenter)
    //   chPos = fv0Det->getReadoutCenter(chno);
    // else
    //   chPos = fv0Det->getCellCenter(chno);

    return RecoDecay::phi(chPos.x + offsetX, chPos.y + offsetY);
  }

  double getEtaFV0(unsigned int chno)
  {
    int cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);
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

    o2::fv0::Point3Dsimple chPos;
    chPos = fv0Det->getReadoutCenter(chno);
    // if (isReadoutCenter)
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

    if (processMc == false) {
      if (!collision.sel8()) {
        return false;
      }
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/hEventCounter"), EventSelectionStep::AfterEventSelection);
    }

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
      if (etaCandidateMax >= 0. && std::abs(etaCandidate) > etaCandidateMax) {
        return false;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandRecoMax) {
        return false;
      }
      return true;
    } else { // For now, that means we do D0
      // Doesn't this exclude D0bar ?
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        return false;
      }
      if (etaCandidateMax >= 0. && std::abs(etaCandidate) > etaCandidateMax) {
        return false;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yD0(candidate)) > yCandRecoMax) {
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
    if (mftTrack.eta() > etaMftTrackMax || mftTrack.eta() < etaMftTrackMin) {
      return false;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/TpcMft/hMftTracksSelection"), MftTrackSelectionStep::Eta);
    }

    // cut on the number of clusters of the reconstructed MFT track
    if (mftTrack.nClusters() < nClustersMftTrack) {
      return false;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/TpcMft/hMftTracksSelection"), MftTrackSelectionStep::Cluster);
    }

    return true;
  }

  // Cut on ambiguous MFT tracks
  template <typename TTrack>
  bool isAmbiguousMftTrack(TTrack const& mftTrack, bool fillHistograms)
  {
    if (mftTrack.ambDegree() > 1) {
      if (fillHistograms) {
        registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfAmbiguousTracks);
      }
      return true;
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::NumberOfNonAmbiguousTracks);
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

        if (etaCandidateMax >= 0. && std::abs(etaCandidate) > etaCandidateMax) {
          return false;
        }

        if (yCandGenMax >= 0. && std::abs(RecoDecay::y(mcCandidate.pVector(), o2::constants::physics::MassD0)) > yCandGenMax) {
          return false;
        }

        // Later on, if I want to add prompt/non-prompt selection, below is how to select prompt only
        // if (!(particle.originMcGen() == RecoDecay::OriginType::Prompt)){
        //   return false;
        // }
      }
    } else { // For now, that means we do LambdaC
      if (std::abs(mcCandidate.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {

        if (etaCandidateMax >= 0. && std::abs(etaCandidate) > etaCandidateMax) {
          return false;
        }

        if (yCandGenMax >= 0. && std::abs(RecoDecay::y(mcCandidate.pVector(), o2::constants::physics::MassLambdaCPlus)) > yCandGenMax) {
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

    if (mcParticle.eta() > etaMftTrackMax || mcParticle.eta() < etaMftTrackMin) {
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

    //
    // TRIGGER PARTICLE
    //
    for (const auto& track1 : tracks1) {

      loopCounter++;

      float eta1 = track1.eta();
      float pt1 = track1.pt();
      float phi1 = track1.phi();
      o2::math_utils::bringTo02Pi(phi1);

      //  TODO: add getter for NUE trigger efficiency here

      //  calculating inv. mass to be filled into the container below
      //  Note: this is needed only in case of HF-hadron correlations
      // TO DO ? Add one more if condition if its MC ?
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
        if (processMc == false) {                                           // If DATA
          if constexpr (!std::is_same_v<aod::MFTTracks, TTracksAssoc>) {    // IF TPC-TPC case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-TPC D0-h
              fillTpcTpcHfChSameEventCandidateQa(multiplicity, track1, true);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-TPC Lc-h
              fillTpcTpcHfChSameEventCandidateQa(multiplicity, track1, false);
            } else { // IF NEITHER D0 NOR LC -> TPC-TPC h-h
              fillTpcTpcChChSameEventQa(multiplicity, track1);
            }
          } else {                                                          // IF TPC-MFT case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-MFT D0-h
              fillTpcMftHfChSameEventCandidateQa(multiplicity, track1, true);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-MFT Lc-h
              fillTpcMftHfChSameEventCandidateQa(multiplicity, track1, false);
            } else { // IF NEITHER D0 NOR LC -> TPC-MFT h-h
              fillTpcMftChChSameEventQa(multiplicity, track1, true);
            } // end of if condition for TPC-TPC or TPC-MFT case
          }
          // Maybe I won't need it for MC (first files are way lighter in MC, but also I need to loop over all tracks in MC GEN)
        } else {                                                            // If MC (add cases later)
          if constexpr (!std::is_same_v<MftTracksMcLabels, TTracksAssoc>) { // IF TPC-TPC case
            // fillTpcTpcChChSameEventQaMc(multiplicity, track1);
          }
        }
      }

      //
      // ASSOCIATED PARTICLE
      //
      for (const auto& track2 : tracks2) {

        // apply cuts for MFT tracks
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {

          if (sameEvent && loopCounter == 1) { // To avoid double counting, we fill the plots only the first time
            registry.fill(HIST("Data/TpcMft/hMftTracksSelection"), MftTrackSelectionStep::NoSelection);

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
        if constexpr (!std::is_same_v<aod::MFTTracks, TTracksAssoc>) {    // if NOT TPC-MFT case -> TPC-TPC case
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
        // if (processMc) {
        if constexpr (std::is_same_v<McParticles, TTracksTrig> || std::is_same_v<McParticles, TTracksAssoc>) {
          if (!isAcceptedMftMcParticle(track2)) {
            continue;
          }
        }

        // if constexpr (std::is_same_v<McParticles, TTracksAssoc>) {
        //   registry.fill(HIST("MC/Gen/TpcMft/HfHadron/SameEvent/hEtaMFT"), track2.eta());
        // }

        float eta2 = track2.eta();
        float pt2 = track2.pt();
        float phi2 = track2.phi();
        o2::math_utils::bringTo02Pi(phi2);

        //  TODO: add getter for NUE associated efficiency here

        //  TODO: add pair cuts on phi*

        float deltaPhi = phi1 - phi2;
        //  set range of delta phi in (-pi/2 , 3/2*pi)
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

        // IF EVERYTHING WORKS WITH THE REASSOCIATED MFT TRACKS, I WILL HAVE TO CHANGE HOW THOSE FUNCTIONS ARE FILLED TOO
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
          // if constexpr (std::is_same_v<FilteredCollisionsWSelMult, TCollisions>) { // If DATA
          if constexpr (!std::is_same_v<aod::MFTTracks, TTracksAssoc>) {    // IF TPC-TPC case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-TPC D0-h
              fillTpcTpcHfChSameEventAssociatedQa(multiplicity, track2, true);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-TPC Lc-h
              fillTpcTpcHfChSameEventAssociatedQa(multiplicity, track2, false);
            }
            // No if condition if it is h-h, because it would be the same plots than for the trigger particle
          } else {                                                          // IF TPC-MFT case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-MFT D0-h
              fillTpcMftHfChSameEventAssociatedQa(multiplicity, track2, true);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-MFT Lc-h
              fillTpcMftHfChSameEventAssociatedQa(multiplicity, track2, false);
            } else { // IF NEITHER D0 NOR LC -> TPC-MFT h-h
              fillTpcMftChChSameEventQa(multiplicity, track2, false);
            } // end of if condition for TPC-TPC or TPC-MFT case
          }
          //} else {                                                        // If MC (add cases later)
          // fillTpcTpcChChSameEventQaMc(multiplicityTracks2, vz, tracks1);
          //}
        }

        if (sameEvent && (loopCounter == 1) && std::is_same_v<aod::MFTTracks, TTracksAssoc>) {
          // FILL USUAL MFT DISTRIBUTIONS
          registry.fill(HIST("Data/TpcMft/kCFStepAll/hEta"), eta2);
          registry.fill(HIST("Data/TpcMft/kCFStepAll/hPhi"), phi2);
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

    //
    // TRIGGER PARTICLE
    //
    for (const auto& track1 : tracks1) {

      loopCounter++;

      float eta1 = track1.eta();
      float pt1 = track1.pt();
      float phi1 = track1.phi();
      o2::math_utils::bringTo02Pi(phi1);

      //  TODO: add getter for NUE trigger efficiency here

      //  calculating inv. mass to be filled into the container below
      //  Note: this is needed only in case of HF-hadron correlations
      // TO DO ? Add one more if condition if its MC ?
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
          // Should add D0 bar ?
        } else { // If Lc
          invmass = hfHelper.invMassLcToPKPi(track1);
          // Should add Lc bar ? (maybe not its the same mass right ?)
        }
      }

      //// Selections for MC GENERATED
      // if constexpr (std::is_same_v<McParticles2ProngMatched, TTracksTrig> || std::is_same_v<McParticles3ProngMatched, TTracksTrig>) {
      //   //  TODO: Check how to put this into a Filter -> Pretty sure it cannot be a filter
      //   if (!isAcceptedMcCandidate(track1)) {
      //     continue;
      //   }
      //   fillingHFcontainer = true;
      //   if constexpr (std::is_same_v<McParticles2ProngMatched, TTracksTrig>) { // If D0
      //     invmass = o2::constants::physics::MassD0;
      //   } else { // If Lc
      //     invmass = o2::constants::physics::MassLambdaCPlus;
      //   }
      // }

      //  fill single-track distributions
      if (!fillingHFcontainer) { // if not HF-h case
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      // FILL QA PLOTS for trigger particle
      if (sameEvent && (cutAmbiguousTracks == false)) {
        if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) {
          fillTpcMftHfChSameEventCandidateQa(multiplicity, track1, true);
        } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
          fillTpcMftHfChSameEventCandidateQa(multiplicity, track1, false);
        } else {
          fillTpcMftChChSameEventQa(multiplicity, track1, true);
        }
      }

      //
      // ASSOCIATED PARTICLE
      //
      for (const auto& track2 : tracks2) {

        // Fill QA plot for all MFT tracks () (only if cutAmbiguousTracks is false to avoid double counting)
        if (!cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
          registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
        }

        // const auto& reassociatedMftTrack = track2.mfttrack();
        //  No one uses const and auto& here, so I will follow

        auto reassociatedMftTrack = track2.template mfttrack_as<aod::MFTTracks>();

        if (!isAcceptedMftTrack(reassociatedMftTrack, false)) {
          continue;
        }

        // Fill QA plot for MFT tracks after physical selection (eta + clusters)
        if (!cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
          registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);
          registry.fill(HIST("Data/TpcMft/hReassociationMftTracks"), ReassociationMftTracks::NotReassociatedMftTracks);
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
            registry.fill(HIST("Data/TpcMft/hReassociationMftTracks"), ReassociationMftTracks::ReassociatedMftTracks);
          }
        }

        //  case of h-h correlations where the two types of tracks are the same
        //  this avoids autocorrelations and double counting of particle pairs
        // if constexpr (std::is_same_v<TTracksAssoc, TTracksTrig>) {
        //  if (track1.index() <= reassociatedMftTrack.index()) {
        //    continue;
        //  }
        //}

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

        //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
        // if (processMc) {
        if constexpr (std::is_same_v<McParticles, TTracksTrig> || std::is_same_v<McParticles, TTracksAssoc>) {
          if (!isAcceptedMftMcParticle(reassociatedMftTrack)) {
            continue;
          }
        }

        // if constexpr (std::is_same_v<McParticles, TTracksAssoc>) {
        //   registry.fill(HIST("MC/Gen/TpcMft/HfHadron/SameEvent/hEtaMFT"), reassociatedMftTrack.eta());
        // }

        float eta2 = reassociatedMftTrack.eta();
        float pt2 = reassociatedMftTrack.pt();
        float phi2 = reassociatedMftTrack.phi();
        o2::math_utils::bringTo02Pi(phi2);

        //  TODO: add getter for NUE associated efficiency here

        //  TODO: add pair cuts on phi*

        float deltaPhi = phi1 - phi2;
        //  set range of delta phi in (-pi/2 , 3/2*pi)
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);

        // IF EVERYTHING WORKS WITH THE REASSOCIATED MFT TRACKS, I WILL HAVE TO CHANGE HOW THOSE FUNCTIONS ARE FILLED TOO
        if (!fillingHFcontainer) {
          //  fill pair correlations
          target->getPairHist()->Fill(step, eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);
        } else {
          target->getPairHist()->Fill(step, eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ, invmass,
                                      triggerWeight * associatedWeight);
        }

        // FILL QA PLOTS for associated particle
        if (sameEvent && (loopCounter == 1) && (cutAmbiguousTracks == false)) {
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) {
            fillTpcMftHfChSameEventAssociatedQa(multiplicity, reassociatedMftTrack, true);
          } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
            fillTpcMftHfChSameEventAssociatedQa(multiplicity, reassociatedMftTrack, false);
          } else {
            fillTpcMftChChSameEventQa(multiplicity, reassociatedMftTrack, false);
          }
        }

        // QA plots for basic MFT distributions for non-ambiguous tracks only (kCFStepTracked)
        if (cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
          // FILL USUAL MFT DISTRIBUTIONS
          registry.fill(HIST("Data/TpcMft/kCFStepTracked/hEta"), eta2);
          registry.fill(HIST("Data/TpcMft/kCFStepTracked/hPhi"), phi2);
        }

      } // end of loop over tracks2
    } // end of loop over tracks 1
  }

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelationsFV0(TTarget target, CorrelationContainer::CFStep step,
                           TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                           float multiplicity, float posZ, bool sameEvent)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    // To avoid filling associated tracks QA many times
    //  I fill it only for the first trigger track of the collision
    auto loopCounter = 0;

    //
    // TRIGGER PARTICLE
    //
    for (auto const& track1 : tracks1) {

      loopCounter++;

      float eta1 = track1.eta();
      float pt1 = track1.pt();
      float phi1 = track1.phi();
      o2::math_utils::bringTo02Pi(phi1);

      //  TODO: add getter for NUE trigger efficiency here

      //  calculating inv. mass to be filled into the container below
      //  Note: this is needed only in case of HF-hadron correlations
      // TO DO ? Add one more if condition if its MC ?
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
        if (!isAcceptedMcCandidate<step>(track1)) {
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
        if (processMc == false) {                                           // If DATA
          if constexpr (!std::is_same_v<aod::MFTTracks, TTracksTrig>) {     // If not aod::MFTTracks as trigger -> TPC-FV0a correlations
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-FV0a D0-h
              fillTpcFv0aHfChSameEventCandidateQa(multiplicity, track1, true);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-FV0a Lc-h
              fillTpcFv0aHfChSameEventCandidateQa(multiplicity, track1, false);
            } else { // IF NEITHER D0 NOR LC -> TPC-MFT h-h
              fillTpcFv0aChChSameEventTriggerQa(multiplicity, track1);
            }
          } else { // If aod::MFTTracks as trigger -> MFT-FV0a (non reassoc/ambiguous) correlations
            fillMftFv0aChChSameEventQa(multiplicity, track1);
            // registry.fill(HIST("Data/TpcMft/kCFStepAll/hEta"), eta2);
            // registry.fill(HIST("Data/TpcMft/kCFStepAll/hPhi"), phi2);
          }
        }
        //} else {                                                                    // If MC (add cases later)
        //  if constexpr (!std::is_same_v<FilteredMftTracksMcLabels, TTracksAssoc>) { // IF TPC-TPC case
        //    fillTpcTpcChChSameEventQaMc(multiplicity, track1);
        //  }
      } // end of if condition to fill QA plots for trigger particle

      //
      // ASSOCIATED PARTICLE
      //
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

        // IF EVERYTHING WORKS WITH THE REASSOCIATED MFT TRACKS, I WILL HAVE TO CHANGE HOW THOSE FUNCTIONS ARE FILLED TOO
        if (!fillingHFcontainer) {
          //  fill pair correlations
          target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);
        } else {
          target->getPairHist()->Fill(step, eta1 - eta2, pt1, pt1, multiplicity, deltaPhi, posZ, invmass,
                                      triggerWeight * associatedWeight);
        }

        // FILL QA PLOTS for associated particle
        if (sameEvent && (loopCounter == 1) && (step == CorrelationContainer::kCFStepReconstructed)) {
          // if constexpr (std::is_same_v<FilteredCollisionsWSelMult, TCollisions>) { // If DATA
          if constexpr (!std::is_same_v<aod::MFTTracks, TTracksTrig>) {     // If not aod::MFTTracks as trigger -> TPC-FV0a correlations
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-FV0a D0-h
              fillTpcFv0aHfChSameEventAssociatedQa(multiplicity, eta2, phi2, true);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-FV0a Lc-h
              fillTpcFv0aHfChSameEventAssociatedQa(multiplicity, eta2, phi2, false);
            } else { // IF NEITHER D0 NOR LC -> ch. part. - ch. part
              fillFv0aQa(multiplicity, eta2, phi2, true);
            }
          } else { // If aod::MFTTracks as trigger -> MFT-FV0a (non reassoc/ambiguous) correlations
            fillFv0aQa(multiplicity, eta2, phi2, false);
          }
          //} else {                                                        // If MC (add cases later)
          // fillTpcTpcChChSameEventQaMc(multiplicityTracks2, vz, tracks1);
          //}
        } // end of if condition to fill QA plots for associated particle

      } // end of loop over FV0 channel indices
    } // end of loop over tracks 1
  }

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelationsFV0ReassociatedMftTracks(TTarget target, CorrelationContainer::CFStep step,
                                                TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                                                float multiplicity, float posZ, bool sameEvent, bool cutAmbiguousTracks)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    // To avoid filling associated tracks QA many times
    //  I fill it only for the first trigger track of the collision
    auto loopCounter = 0;

    //
    // TRIGGER PARTICLE
    //
    for (auto const& track1 : tracks1) {
      loopCounter++;

      // TO-DO (if useful) : adapt this but or Mft-FV0a
      // Fill QA plot for all MFT tracks () (only if cutAmbiguousTracks is false to avoid double counting)
      // if (!cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
      //  registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
      //}

      auto reassociatedMftTrack = track1.template mfttrack_as<aod::MFTTracks>();

      if (!isAcceptedMftTrack(reassociatedMftTrack, false)) {
        continue;
      }

      // TO-DO (if useful) : adapt this but or Mft-FV0a
      // Fill QA plot for MFT tracks after physical selection (eta + clusters)
      // if (!cutAmbiguousTracks && sameEvent && (loopCounter == 1)) {
      //  registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);
      //  registry.fill(HIST("Data/TpcMft/hReassociationMftTracks"), ReassociationMftTracks::NotReassociatedMftTracks);
      //}

      // We check if the track is ambiguous or non-ambiguous (QA plots are filled in isAmbiguousMftTrack)
      // Fill plots only if cutAmbiguousTracks is false (to avoid double counting)
      if (isAmbiguousMftTrack(track1, false)) {
        // If the MFT track is ambiguous we may cut or not on the ambiguous track
        if (cutAmbiguousTracks) {
          continue;
        }
      }

      // if (reassociatedMftTrack.collisionId() != track2.bestCollisionId()) {
      //   if (sameEvent && (loopCounter == 1)) {
      //     registry.fill(HIST("Data/TpcMft/hReassociationMftTracks"), ReassociationMftTracks::ReassociatedMftTracks);
      //   }
      // }

      float eta1 = reassociatedMftTrack.eta();
      float pt1 = reassociatedMftTrack.pt();
      float phi1 = reassociatedMftTrack.phi();
      o2::math_utils::bringTo02Pi(phi1);

      //  TODO: add getter for NUE trigger efficiency here

      target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);

      // FILL QA PLOTS for trigger particle
      if (sameEvent && (step == CorrelationContainer::kCFStepReconstructed)) {
        if (processMc == false) { // If DATA
          fillMftFv0aChChSameEventQa(multiplicity, reassociatedMftTrack);
        }
      } // end of if condition to fill QA plots for trigger particle

      //
      // ASSOCIATED PARTICLE
      //
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
          fillFv0aQa(multiplicity, eta2, phi2, false);
        } // end of if condition to fill QA plots for associated particle

      } // end of loop over FV0 channel indices
    } // end of loop over tracks 1
  }

  // ===============================================================================================================================================================================
  //      mixCollisions for RECONSTRUCTED events
  // ===============================================================================================================================================================================

  template <typename TCollisions, typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisions(TCollisions const& collisions, CorrelationContainer::CFStep step,
                     TTracksTrig const& tracks1, TTracksAssoc const& tracks2,
                     TLambda getPartsSize,
                     OutputObj<CorrelationContainer>& corrContainer)
  {
    // The first one that I call "Data" should work for data and mc rec
    using BinningTypeData = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::collision::PosZ, decltype(getPartsSize)>;

    BinningTypeData binningWithTracksSize{{getPartsSize}, {binsMixingVertex, binsMixingMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<TCollisions, TTracksTrig, TTracksAssoc, BinningTypeData> pair{binningWithTracksSize, nMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if constexpr (!std::is_same_v<FilteredMcCollisions, TCollisions>) { // if NOT MC -> do collision cut
        if (!(isAcceptedCollision(collision1, false))) {
          continue;
        }
        if (!(isAcceptedCollision(collision2, false))) {
          continue;
        }
      }

      // auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      // int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicityTracks1 = getPartsSize(collision1);

      /*
      if constexpr (std::is_same_v<FilteredCollisionsWSelMultMcLabels, TCollisions>) { // If MC
        registry.fill(HIST("MC/Rec/TpcTpc/ChPartChPart/MixedEvent/hEventCountMixing"), bin);
      } else {                                                                                                              // If not MC
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {                                                       // IF TPC-MFT case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-MFT HF-h
            registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing"), bin);
          } else { // IF h-h case -> TPC-MFT h-h case
            registry.fill(HIST("Data/TpcMft/ChPartChPart/MixedEvent/hEventCountMixing"), bin);
          }
        } else {                                                                                                            // IF TPC-TPC case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-TPC HF-h
            registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
          } else { // IF h-h case -> TPC-TPC h-h case
            registry.fill(HIST("Data/TpcTpc/ChPartChPart/MixedEvent/hEventCountMixing"), bin);
          }
        } // end of if condition for TPC-TPC or TPC-MFT case
      }
      */

      corrContainer->fillEvent(multiplicityTracks1, step);
      fillCorrelations(corrContainer, step, tracks1, tracks2, multiplicityTracks1, collision1.posZ(), false);
    }
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
    Pair<TCollisions, TTracksTrig, TTracksAssoc, BinningTypeData> pair{binningWithTracksSize, nMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if constexpr (!std::is_same_v<FilteredMcCollisions, TCollisions>) { // if NOT MC -> do collision cut
        if (!(isAcceptedCollision(collision1, false))) {
          continue;
        }
        if (!(isAcceptedCollision(collision2, false))) {
          continue;
        }
      }

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicityTracks1 = getPartsSize(collision1);



      if constexpr (std::is_same_v<FilteredCollisionsWSelMultMcLabels, TCollisions>) { // If MC
        registry.fill(HIST("MC/Rec/TpcTpc/ChPartChPart/MixedEvent/hEventCountMixing"), bin);
      } else {                                                                                                              // If not MC
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {                                              // IF TPC-MFT case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-MFT HF-h
            registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing"), bin);
          } else { // IF h-h case -> TPC-MFT h-h case
            registry.fill(HIST("Data/TpcMft/ChPartChPart/MixedEvent/hEventCountMixing"), bin);
          }
        } else {                                                                                                            // IF TPC-TPC case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-TPC HF-h
            registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
          } else { // IF h-h case -> TPC-TPC h-h case
            registry.fill(HIST("Data/TpcTpc/ChPartChPart/MixedEvent/hEventCountMixing"), bin);
          }
        } // end of if condition for TPC-TPC or TPC-MFT case
      }


      corrContainer->fillEvent(multiplicityTracks1, step);
      fillCorrelationsReassociatedMftTracks(corrContainer, step, tracks1, tracks2, multiplicityTracks1, collision1.posZ(), false, cutAmbiguousTracks, field );
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

    BinningTypeMcTruth binningWithTracksSize{{getPartsSize}, {binsMixingVertex, binsMixingMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<TMcCollisions, TTracksTrig, TTracksAssoc, BinningTypeMcTruth> pair{binningWithTracksSize, nMixedEvents, -1, mcCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      // auto binningValues = binningWithTracksSize.getBinningValues(collision1, mcCollisions);
      //  int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicity = getPartsSize(collision1); // get multiplicity of charged hadrons, which is used for slicing in mixing

      // TO BE DONE : ADD ONE MORE IF CONDITION TO FILL THE MC CASE
      // TODO : FILL NEW PLOTS FOR MCTRUTH ONLY
      // registry.fill(HIST("MC/Gen/TpcTpc/ChPartChPart/MixedEvent/hEventCountMixing"), bin);

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
    const auto multiplicity = collision.multNTracksPV();
    // registry.fill(HIST("Data/TpcTpc/ChPartChPart/SameEvent/hMultiplicity"), multiplicity);
    // registry.fill(HIST("Data/TpcTpc/ChPartChPart/SameEvent/hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcTpc/ChPartChPart/SameEvent/hEventCountSame"), bin);

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
    if (doReferenceFlow)
      fillEventSelectionPlots = false;

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEventCountSame"), bin);

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
    if (doReferenceFlow)
      fillEventSelectionPlots = false;

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEventCountSame"), bin);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, tracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcLcCh, "DATA : Process same-event correlations for TPC-TPC Lc-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT h-h case
  // =====================================

  void processSameTpcMftChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             FilteredTracksWDcaSel const& tracks,
                             aod::MFTTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/ChPartChPart/SameEvent/hEventCountSame"), bin);

    // I use kCFStepAll for running my code with all MFTTracks were the reassociation process was not applied
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChCh, "DATA : Process same-event correlations for TPC-MFT h-h case", false);

  void processSameTpcMftChChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         FilteredTracksWDcaSel const& tracks,
                                         aod::MFTTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    registry.fill(HIST("Data/TpcMft/hNTracks"), tracks.size());
    registry.fill(HIST("Data/TpcMft/hNMftTracks"), mftTracks.size());
    registry.fill(HIST("Data/TpcMft/hNBestCollisionFwd"), reassociatedMftTracks.size());

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/ChPartChPart/SameEvent/hEventCountSame"), bin);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChChReassociated, "DATA : Process same-event correlations for TPC-MFT h-h case reassociated", false);

  void processSameTpcMftChChNonAmbiguous(FilteredCollisionsWSelMult::iterator const& collision,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         FilteredTracksWDcaSel const& tracks,
                                         aod::MFTTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    registry.fill(HIST("Data/TpcMft/hNTracks"), tracks.size());
    registry.fill(HIST("Data/TpcMft/hNMftTracks"), mftTracks.size());
    registry.fill(HIST("Data/TpcMft/hNBestCollisionFwd"), reassociatedMftTracks.size());

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/ChPartChPart/SameEvent/hEventCountSame"), bin);

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
                             aod::MFTTracks const& mftTracks)
  {
    auto fillEventSelectionPlots = true;

    // When doing reference flow, two cases are used (HF-h, h-h) and thus eventSelectionPlots was filled twice
    if (doReferenceFlow)
      fillEventSelectionPlots = false;

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hEventCountSame"), bin);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftD0Ch, "DATA : Process same-event correlations for TPC-MFT D0-h case", false);

  void processSameTpcMftD0ChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         HfCandidatesSelD0 const& candidates,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         aod::MFTTracks const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    const auto multiplicity = collision.multNTracksPV();

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
                             aod::MFTTracks const& mftTracks)
  {
    auto fillEventSelectionPlots = true;

    // When doing reference flow, two cases are used (HF-h, h-h) and thus eventSelectionPlots was filled twice
    if (doReferenceFlow)
      fillEventSelectionPlots = false;

    if (!(isAcceptedCollision(collision, fillEventSelectionPlots))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hEventCountSame"), bin);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftLcCh, "DATA : Process same-event correlations for TPC-MFT Lc-h case", false);

  void processSameTpcMftLcChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         HfCandidatesSelLc const& candidates,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         aod::MFTTracks const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/ChPartChPart/SameEvent/hEventCountSame"), bin);

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
                              aod::FV0As const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFV0(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, fv0, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFv0aChCh, "DATA : Process same-event correlations for TPC-FV0-A h-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FV0A D0 - Ch. Part
  // =====================================

  void processSameTpcFv0aD0Ch(FilteredCollisionsWSelMult::iterator const& collision,
                              HfCandidatesSelD0 const& candidates,
                              aod::FV0As const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFV0(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, fv0, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFv0aD0Ch, "DATA : Process same-event correlations for TPC-FV0-A D0-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-FV0A Lc - Ch. Part
  // =====================================

  void processSameTpcFv0aLcCh(FilteredCollisionsWSelMult::iterator const& collision,
                              HfCandidatesSelLc const& candidates,
                              aod::FV0As const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFV0(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, fv0, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcFv0aLcCh, "DATA : Process same-event correlations for TPC-FV0-A Lc-h case", false);

  // =====================================
  //    DATA : process same event correlations: MFT-FV0A Ch. Part. - Ch. Part
  // =====================================

  void processSameMftFv0aChCh(FilteredCollisionsWSelMult::iterator const& collision,
                              aod::MFTTracks const& mftTracks,
                              aod::FV0As const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFV0(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, mftTracks, fv0, multiplicity, collision.posZ(), true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFv0aChCh, "DATA : Process same-event correlations for MFT-FV0-A h-h case", false);

  void processSameMftFv0aChChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                          soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                          aod::MFTTracks const& /*mftTracks*/,
                                          aod::FV0As const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFV0ReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, reassociatedMftTracks, fv0, multiplicity, collision.posZ(), true, false);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFv0aChChReassociated, "DATA : Process same-event correlations for MFT-FV0a h-h case reassociated", false);

  void processSameMftFv0aChChNonAmbiguous(FilteredCollisionsWSelMult::iterator const& collision,
                                          soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                          aod::MFTTracks const& /*mftTracks*/,
                                          aod::FV0As const&)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();

    if (collision.has_foundFV0()) {
      const auto& fv0 = collision.foundFV0();

      sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsFV0ReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, reassociatedMftTracks, fv0, multiplicity, collision.posZ(), true, true);
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processSameMftFv0aChChNonAmbiguous, "DATA : Process same-event correlations for MFT-FV0a h-h non-ambiguous case", false);

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

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(mcCollision.posZ(), multiplicity));
    // registry.fill(HIST("MC/Gen/TpcMft/HfHadron/SameEvent/2Prong/hEventCountSame"), bin);

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

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(mcCollision.posZ(), multiplicity));
    // registry.fill(HIST("MC/Gen/TpcMft/HfHadron/SameEvent/2Prong/hEventCountSame"), bin);

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
    //  we want to group collisions based on charged-track multiplicity
    // auto getTracksSize = [&tracks, this](FilteredCollisionsWSelMult::iterator const& col) {
    //   auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
    //   auto size = associatedTracks.size();
    //   return size;
    //  };

    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    // mixCollisions(collisions, tracks, tracks, getTracksSize, mixedEvent);
    mixCollisions(collisions, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, tracks, getMultiplicity, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChCh, "DATA : Process mixed-event correlations for TPC-TPC h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-TPC HF-h case for D0
  // =====================================

  void processMixedTpcTpcD0Ch(FilteredCollisionsWSelMult const& collisions,
                              FilteredTracksWDcaSel const& tracks,
                              HfCandidatesSelD0 const& candidates)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, tracks, getMultiplicity, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcD0Ch, "DATA : Process mixed-event correlations for TPC-TPC D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-TPC HF-h case for Lc
  // =====================================

  void processMixedTpcTpcLcCh(FilteredCollisionsWSelMult const& collisions,
                              FilteredTracksWDcaSel const& tracks,
                              HfCandidatesSelLc const& candidates)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, tracks, getMultiplicity, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcLcCh, "DATA : Process mixed-event correlations for TPC-TPC Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT h-h case
  // =====================================

  void processMixedTpcMftChCh(FilteredCollisionsWSelMult const& collisions,
                              FilteredTracksWDcaSel const& tracks,
                              aod::MFTTracks const& mftTracks)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, tracks, mftTracks, getMultiplicity, mixedEvent);
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
                              aod::MFTTracks const& mftTracks,
                              FilteredTracksWDcaSel const& /*tracks*/)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, mftTracks, getMultiplicity, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftD0Ch, "DATA : Process mixed-event correlations for TPC-MFT D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT HF-h case
  // =====================================

  void processMixedTpcMftLcCh(FilteredCollisionsWSelMult const& collisions,
                              HfCandidatesSelLc const& candidates,
                              aod::MFTTracks const& mftTracks)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, CorrelationContainer::kCFStepReconstructed, candidates, mftTracks, getMultiplicity, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftLcCh, "DATA : Process mixed-event correlations for TPC-MFT Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A ch part. - ch. part. case
  // =====================================

  void processMixedTpcFv0aChCh(FilteredCollisionsWSelMult const& collisions,
                               FilteredTracksWDcaSel const& tracks,
                               aod::FV0As const&)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {binsMixingVertex, binsMixingMultiplicity}, true};

    for (auto const& [collision1, collision2] : soa::selfCombinations(binningOnVtxAndMult, nMixedEvents, -1, collisions, collisions)) {

      if (!isAcceptedCollision(collision1) || !isAcceptedCollision(collision2)) {
        continue;
      }

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      if (collision1.has_foundFV0() && collision2.has_foundFV0()) {

        const auto multiplicity = getMultiplicity(collision1);

        auto slicedTriggerTracks = tracks.sliceBy(perColTracks, collision1.globalIndex());
        const auto& fv0 = collision2.foundFV0();

        mixedEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
        fillCorrelationsFV0(mixedEvent, CorrelationContainer::CFStep::kCFStepReconstructed, slicedTriggerTracks, fv0, multiplicity, collision1.posZ(), false);
      }
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFv0aChCh, "DATA : Process mixed-event correlations for TPC-FV0-A h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A D0 - ch. part. case
  // =====================================

  void processMixedTpcFv0aD0Ch(FilteredCollisionsWSelMult const& collisions,
                               HfCandidatesSelD0 const& candidates,
                               aod::FV0As const&)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {binsMixingVertex, binsMixingMultiplicity}, true};

    for (auto const& [collision1, collision2] : soa::selfCombinations(binningOnVtxAndMult, nMixedEvents, -1, collisions, collisions)) {

      if (!isAcceptedCollision(collision1) || !isAcceptedCollision(collision2)) {
        continue;
      }

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      if (collision1.has_foundFV0() && collision2.has_foundFV0()) {

        const auto multiplicity = getMultiplicity(collision1);

        auto slicedTriggerCandidates = candidates.sliceBy(perColD0s, collision1.globalIndex());
        const auto& fv0 = collision2.foundFV0();

        mixedEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
        fillCorrelationsFV0(mixedEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, slicedTriggerCandidates, fv0, multiplicity, collision1.posZ(), false);
      }
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFv0aD0Ch, "DATA : Process mixed-event correlations for TPC-FV0-A D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A Lc - ch. part. case
  // =====================================

  void processMixedTpcFv0aLcCh(FilteredCollisionsWSelMult const& collisions,
                               HfCandidatesSelLc const& candidates,
                               aod::FV0As const&)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {binsMixingVertex, binsMixingMultiplicity}, true};

    for (auto const& [collision1, collision2] : soa::selfCombinations(binningOnVtxAndMult, nMixedEvents, -1, collisions, collisions)) {

      if (!isAcceptedCollision(collision1) || !isAcceptedCollision(collision2)) {
        continue;
      }

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      if (collision1.has_foundFV0() && collision2.has_foundFV0()) {

        const auto multiplicity = getMultiplicity(collision1);

        auto slicedTriggerCandidates = candidates.sliceBy(perColLcs, collision1.globalIndex());
        const auto& fv0 = collision2.foundFV0();

        mixedEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
        fillCorrelationsFV0(mixedEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, slicedTriggerCandidates, fv0, multiplicity, collision1.posZ(), false);
      }
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcFv0aLcCh, "DATA : Process mixed-event correlations for TPC-FV0-A Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-FV0-A ch part. - ch. part. case
  // =====================================

  void processMixedMftFv0aChCh(FilteredCollisionsWSelMult const& collisions,
                               aod::MFTTracks const& mftTracks,
                               aod::FV0As const&)
  {
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {binsMixingVertex, binsMixingMultiplicity}, true};

    for (auto const& [collision1, collision2] : soa::selfCombinations(binningOnVtxAndMult, nMixedEvents, -1, collisions, collisions)) {

      if (!isAcceptedCollision(collision1) || !isAcceptedCollision(collision2)) {
        continue;
      }

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      if (collision1.has_foundFV0() && collision2.has_foundFV0()) {

        const auto multiplicity = getMultiplicity(collision1);
        auto slicedTriggerMftTracks = mftTracks.sliceBy(perColMftTracks, collision1.globalIndex());
        const auto& fv0 = collision2.foundFV0();

        mixedEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
        fillCorrelationsFV0(mixedEvent, CorrelationContainer::CFStep::kCFStepReconstructed, slicedTriggerMftTracks, fv0, multiplicity, collision1.posZ(), false);
      }
    }
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedMftFv0aChCh, "DATA : Process mixed-event correlations for Mft-FV0-A h-h case", false);

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
    if (centralityBinsForMc) {
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
