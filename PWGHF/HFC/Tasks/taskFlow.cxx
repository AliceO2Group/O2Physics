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
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsPid.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
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

#include <cstdlib>
#include <string>
#include <tuple>
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

enum EventSelectionStep {
  AllEvents = 0,
  AfterEventSelection,
  NEventSelectionSteps
};

// static constexpr float kPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};

struct HfTaskFlow {

  //  configurables for processing options

  Configurable<bool> centralityBinsForMc{"centralityBinsForMc", false, "false = OFF, true = ON for data like multiplicity/centrality bins for MC steps"};
  Configurable<float> mftMaxDCAxy{"mftMaxDCAxy", 2.0f, "Cut on dcaXY for MFT tracks"};
  Configurable<bool> doReferenceFlow{"doReferenceFlow", false, "Flag to know if reference flow should be done"};
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
  // Service<o2::ccdb::BasicCCDBManager> ccdb;
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

  Preslice<aod::MFTTracks> perCol = o2::aod::fwdtrack::collisionId;

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
  ConfigurableAxis axisEtaMft{"axisEtaMft", {48, -4, -2}, "eta axis for MFT histograms"};
  ConfigurableAxis axisEtaTpc{"axisEtaTpc", {48, -1, 1}, "eta axis for TPC histograms"};
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

  //  =========================
  //      init()
  //  =========================
  void init(InitContext&)
  {
    const int nBinsMix = axisMultiplicity->size() * axisVertex->size();

    //  =========================
    //      Event histograms
    // TO-DO : do i have to separate event histograms between DATA and MC ?
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

    for (int iBin = 0; iBin < NMftAmbiguitySteps; iBin++) {
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

    //  =========================
    //      DATA : histograms for TPC-TPC h-h case
    //  =========================

    // DATA : event histograms for TPC-TPC h-h same event
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {axisVertex}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    // DATA : associated particles histograms for TPC-TPC h-h same event
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1D, {axisEtaTpc}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTpc}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {axisMultiplicity, axisEtaTpc, axisPhi}});

    // Katarina had this :
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hVzEta", "eta vs. Vz", {HistType::kTH2F, {axisEtaTpc, axisVertex}});

    // DATA : event mixing histograms for TPC-TPC h-h mixed event
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    //  =========================
    //      DATA : histograms for TPC-TPC HF-h case for 2PRONG
    //  =========================

    // DATA : event histograms for TPC-TPC HF-h same event
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hEta", "eta", {HistType::kTH1D, {axisEtaTpc}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hPhi", "phi", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTpc}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {axisMultiplicity, axisEtaTpc, axisPhi}});

    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEventCountSame", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEta", "eta", {HistType::kTH1D, {axisEtaTpc}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPhi", "phi", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    // DATA : trigger particles (candidates) histograms for TPC-TPC D0-h same event
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtCandidate", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisMass, axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPtVsMult", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {axisMass, axisPt, axisMultiplicity}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {axisMass}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEtaCandVsPt", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {axisEtaTpc, axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, axisPt}});

    //  =========================
    //      DATA : histograms for TPC-TPC HF-h case for 3PRONG
    //  ===================

    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEventCountSame", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPt", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisMass, axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMass", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {axisMass}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPtVsMult", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {axisMass, axisPt, axisMultiplicity}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1D, {axisEtaTpc}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {axisEtaTpc, axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {axisPhi, axisPt}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, axisPt}});

    //  =========================
    //      DATA : histograms for TPC-MFT h-h case
    //  =========================

    // DATA : trigger particles (TPC tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiTPC", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {axisMultiplicity, axisEtaTpc, axisPhi}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaTPC", "etaTPC", {HistType::kTH1D, {axisEtaTpc}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPhiTPC", "phiTPC", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPtTPC", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hYieldsTPC", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTpc}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hMultiplicityTPC", "multiplicity;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}});

    // DATA : associated particles (MFT tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {axisMultiplicity, axisEtaMft, axisPhi}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaMFT", "etaMFT", {HistType::kTH1D, {axisEtaMft}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPhiMFT", "phiMFT", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPtMFT", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hYieldsMFT", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaMft}});

    // DATA : histograms for TPC-MFT h-h event mixing for events QA
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    //  =========================
    //      DATA : histograms for TPC-MFT HF-h case FOR 2PRONG
    //  =========================

    // DATA : trigger particles (candidates) histograms for TPC-MFT HF-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEventCountSame", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaPhiCandidate", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {axisMultiplicity, axisEtaMft, axisPhi}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaCandidate", "etaTPC", {HistType::kTH1D, {axisEtaTpc}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPhiCandidate", "phiTPC", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hYieldsCandidate", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTpc}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hMultiplicityCandidate", "multiplicity;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}});

    // DATA : trigger particles (candidates) histograms for TPC-MFT HF-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtCandidate", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisMass, axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {axisMass}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPtVsMult", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {axisMass, axisPt, axisMultiplicity}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaCandVsPt", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {axisEtaTpc, axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, axisPt}});

    // DATA : associated particles (MFT tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {axisMultiplicity, axisEtaMft, axisPhi}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaMFT", "etaMFT", {HistType::kTH1D, {axisEtaMft}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPhiMFT", "phiMFT", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPtMFT", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hYieldsMFT", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaMft}});

    // DATA : histograms for TPC-MFT h-h event mixing for events QA
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    //  =========================
    //      DATA : histograms for TPC-MFT HF-h case FOR 3PRONG
    //  =========================

    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEventCountSame", "bin", {HistType::kTH1D, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hYieldsCandidate", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTpc}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMultiplicityCandidate", "multiplicity;multiplicity;entries", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEtaPhiCandidate", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {axisMultiplicity, axisEtaTpc, axisPhi}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPt", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisMass, axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMass", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1D, {axisMass}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPtVsMult", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {axisMass, axisPt, axisMultiplicity}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1D, {axisEtaTpc}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1D, {axisPhi}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {axisEtaTpc, axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {axisPhi, axisPt}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, axisPt}});

    //  =========================
    //      DATA : histograms for MFT reassociated distributions (to plot all MFT tracks, and only non-ambiguous MFT trakcs)
    //             ressociated tracks distribution can be found in the normal eta and phi distributions of each individual cases
    //  =========================

    // All MFT tracks
    registry.add("Data/TpcMft/kCFStepAll/hEta", "eta", {HistType::kTH1D, {axisEtaMft}});
    registry.add("Data/TpcMft/kCFStepAll/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

    // Only non-ambiguous MFT tracks
    registry.add("Data/TpcMft/kCFStepTracked/hEta", "eta", {HistType::kTH1D, {axisEtaMft}});
    registry.add("Data/TpcMft/kCFStepTracked/hPhi", "phi", {HistType::kTH1D, {axisPhi}});

    //  =========================
    //      MC : histograms for TPC-TPC h-h case
    //  =========================

    // MC reconstructed

    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {axisVertex}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    // Katarina had this :
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary", "hMultiplicityPrimary", {HistType::kTH1F, {axisMultiplicity}});
    //  histograms for MC associated particles
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {axisEtaTpc}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {axisPhi}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTpc}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {axisMultiplicity, axisEtaTpc, axisPhi}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hNtracks", "hNtracks", {HistType::kTH1F, {axisMultiplicity}});
    //  histograms for MC particles in event mixing
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {axisVertex}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {axisEtaTpc}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {axisPhi}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {axisMultiplicity}});

    // MC Truth

    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {axisVertex}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    // Katarina had this :
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary", "hMultiplicityPrimary", {HistType::kTH1F, {axisMultiplicity}});
    //  histograms for MC associated particles
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {axisEtaTpc}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {axisPhi}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {axisMultiplicity, axisPt, axisEtaTpc}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {axisMultiplicity, axisEtaTpc, axisPhi}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hNtracks", "hNtracks", {HistType::kTH1F, {axisMultiplicity}});
    //  histograms for MC particles in event mixing
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {axisMultiplicity}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {axisVertex}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {axisPt}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {axisEtaTpc}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {axisPhi}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {axisMultiplicity}});

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

    // initialization of correlation containers for data
    sameEvent.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
    mixedEvent.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));
    sameEventHf.setObject(new CorrelationContainer("sameEventHf", "sameEventHf", corrAxis, effAxis, userAxis));
    mixedEventHf.setObject(new CorrelationContainer("mixedEventHf", "mixedEventHf", corrAxis, effAxis, userAxis));

    // initialization of correlation containes for monte-carlo
    sameEventHfMc.setObject(new CorrelationContainer("sameEventHfMc", "sameEventHfMc", corrAxis, effAxis, userAxis));
    mixedEventHfMc.setObject(new CorrelationContainer("mixedEventHfMc", "mixedEventHfMc", corrAxis, effAxis, userAxis));
  } // End of init() function

  // =========================
  //      Quality Assesment plots for Same Event
  // =========================

  // ---- DATA : TPC-TPC h-h Same Event QA ----
  template <typename TTrack>
  void fillTpcTpcChChSameEventQa(float multiplicity, TTrack const& track)
  {
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hPt"), track.pt());
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hEta"), track.eta());
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hPhi"), track.phi());
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hYields"), multiplicity, track.pt(), track.eta());
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hEtaPhi"), multiplicity, track.eta(), track.phi());
  }

  template <typename TTrack>
  void fillTpcTpcHfChSameEventQa(float multiplicity, TTrack const& track)
  {
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/hPt"), track.pt());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/hEta"), track.eta());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/hPhi"), track.phi());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/hYields"), multiplicity, track.pt(), track.eta());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/hEtaPhi"), multiplicity, track.eta(), track.phi());
  }

  // ---- MC REC : TPC-TPC h-h Same Event QA ----
  template <typename TTrack>
  void fillTpcTpcChChSameEventQaMc(float multiplicity, TTrack const& track)
  {
    if constexpr (std::is_same_v<FilteredTracksWDcaSelMC, TTrack>) { // if MC Rec
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPt"), track.pt());
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEta"), track.eta());
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPhi"), track.phi());
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hYields"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEtaPhi"), multiplicity, track.eta(), track.phi());
    } else { // if MC Gen
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPt"), track.pt());
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEta"), track.eta());
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPhi"), track.phi());
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hYields"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEtaPhi"), multiplicity, track.eta(), track.phi());
    }
  }

  // ---- DATA : TPC-MFT h-h Same Event QA ----
  template <typename TTrack>
  void fillTpcMftChChSameEventQa(float multiplicity, TTrack const& track, bool isTPC)
  {
    float phi = track.phi();
    o2::math_utils::bringTo02Pi(phi);

    if (isTPC) { // trigger hadron from TPC
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaTPC"), track.eta());
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPhiTPC"), phi);
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiTPC"), multiplicity, track.eta(), phi);
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPtTPC"), track.pt());
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hYieldsTPC"), multiplicity, track.pt(), track.eta());
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hMultiplicityTPC"), multiplicity);
      // add multiplicity plot?
    } else { // associated hadron from MFT
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaMFT"), track.eta());
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPhiMFT"), phi);
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiMFT"), multiplicity, track.eta(), phi);
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPtMFT"), track.pt());
      registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hYieldsMFT"), multiplicity, track.pt(), track.eta());
    }
  }

  // ---- DATA : TPC-MFT HF-h Same Event QA ----

  template <typename TTrack>
  void fillTpcMftHfChSameEventQa(float multiplicity, TTrack const& track)
  {
    // Used to fill QA plots for associated track from MFT when doing TPC-MFT HF-h correlations
    float phi = track.phi();
    o2::math_utils::bringTo02Pi(phi);

    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hEtaMFT"), track.eta());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hPhiMFT"), phi);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hEtaPhiMFT"), multiplicity, track.eta(), phi);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hPtMFT"), track.pt());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hYieldsMFT"), multiplicity, track.pt(), track.eta());
    // add plot for multiplicity ?
  }

  // ---- DATA : TPC-TPC HF-h Same Event (Candidates) QA ----
  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTrack>
  void fillTpcTpcD0CandidateQa(float multiplicity, TTrack const& candidate)
  {
    float phi = candidate.phi();
    o2::math_utils::bringTo02Pi(phi);
    auto eta = candidate.eta();
    auto pt = candidate.pt();

    if (candidate.isSelD0() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPtVsMult"), hfHelper.invMassD0ToPiK(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPt"), hfHelper.invMassD0ToPiK(candidate), pt);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMass"), hfHelper.invMassD0ToPiK(candidate));
    }
    if (candidate.isSelD0bar() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPtVsMult"), hfHelper.invMassD0barToKPi(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPt"), hfHelper.invMassD0barToKPi(candidate), pt);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMass"), hfHelper.invMassD0barToKPi(candidate));
    }

    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMultiplicity"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEta"), eta);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPhi"), phi);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtCandidate"), pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng1"), candidate.ptProng1());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEtaCandVsPt"), eta, pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), pt);
  }

  // ---- DATA : TPC-TPC HF-h Same Event (Candidates) QA ----
  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTrack>
  void fillTpcTpcLcCandidateQa(float multiplicity, TTrack const& candidate)
  {
    float phi = candidate.phi();
    o2::math_utils::bringTo02Pi(phi);

    auto pt = candidate.pt();
    auto ptProng0 = candidate.ptProng0();
    auto ptProng1 = candidate.ptProng1();
    auto ptProng2 = candidate.ptProng2();
    auto eta = candidate.eta();

    if (candidate.isSelLcToPKPi() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMass"), hfHelper.invMassLcToPKPi(candidate));
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPtVsMult"), hfHelper.invMassLcToPKPi(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPt"), hfHelper.invMassLcToPKPi(candidate), pt);
    }
    if (candidate.isSelLcToPiKP() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMass"), hfHelper.invMassLcToPiKP(candidate));
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPtVsMult"), hfHelper.invMassLcToPiKP(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPt"), hfHelper.invMassLcToPiKP(candidate), pt);
    }

    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMultiplicity"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPt"), pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng0"), ptProng0);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng1"), ptProng1);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng2"), ptProng2);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEta"), eta);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEtaVsPt"), eta, pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPhi"), phi);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPhiVsPt"), phi, pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hSelectionStatus"), candidate.isSelLcToPKPi(), pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hSelectionStatus"), candidate.isSelLcToPiKP(), pt);
  }

  // ---- DATA : TPC-MFT HF-h Same Event (Candidates) QA ----
  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTrack>
  void fillTpcMftD0CandidateQa(float multiplicity, TTrack const& candidate)
  {
    float phi = candidate.phi();
    auto pt = candidate.pt();
    o2::math_utils::bringTo02Pi(phi);

    if (candidate.isSelD0() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPtVsMult"), hfHelper.invMassD0ToPiK(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPt"), hfHelper.invMassD0ToPiK(candidate), pt);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMass"), hfHelper.invMassD0ToPiK(candidate));
    }
    if (candidate.isSelD0bar() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPtVsMult"), hfHelper.invMassD0barToKPi(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPt"), hfHelper.invMassD0barToKPi(candidate), pt);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMass"), hfHelper.invMassD0barToKPi(candidate));
    }

    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaCandidate"), candidate.eta());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hPhiCandidate"), phi);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMultiplicityCandidate"), multiplicity);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaPhiCandidate"), multiplicity, candidate.eta(), phi);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hYieldsCandidate"), multiplicity, pt, candidate.eta());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtCandidate"), pt);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtProng1"), candidate.ptProng1());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaCandVsPt"), candidate.eta(), pt);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), pt);
  }

  // ---- DATA : TPC-MFT HF-h Same Event (Candidates) QA ----
  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works

  template <typename TTrack>
  void fillTpcMftLcCandidateQa(float multiplicity, TTrack const& candidate)
  {
    auto pt = candidate.pt();
    auto ptProng0 = candidate.ptProng0();
    auto ptProng1 = candidate.ptProng1();
    auto ptProng2 = candidate.ptProng2();
    float phi = candidate.phi();
    o2::math_utils::bringTo02Pi(phi);

    if (candidate.isSelLcToPKPi() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMass"), hfHelper.invMassLcToPKPi(candidate));
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPtVsMult"), hfHelper.invMassLcToPKPi(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPt"), hfHelper.invMassLcToPKPi(candidate), pt);
    }
    if (candidate.isSelLcToPiKP() >= selectionFlagHf) {
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMass"), hfHelper.invMassLcToPiKP(candidate));
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPtVsMult"), hfHelper.invMassLcToPiKP(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPt"), hfHelper.invMassLcToPiKP(candidate), pt);
    }

    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hYieldsCandidate"), multiplicity, pt, candidate.eta());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hEtaPhiCandidate"), multiplicity, candidate.eta(), phi);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMultiplicityCandidate"), multiplicity);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hPt"), pt);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng0"), ptProng0);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng1"), ptProng1);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng2"), ptProng2);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hEta"), candidate.eta());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hEtaVsPt"), candidate.eta(), pt);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hPhi"), candidate.phi());
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hPhiVsPt"), candidate.phi(), pt);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hSelectionStatus"), candidate.isSelLcToPKPi(), pt);
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hSelectionStatus"), candidate.isSelLcToPiKP(), pt);
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
      if (std::abs(mcCandidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {

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
      if (std::abs(mcCandidate.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {

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
          // Should add D0 bar ?
        } else { // If Lc
          invmass = hfHelper.invMassLcToPKPi(track1);
          // Should add Lc bar ? (maybe not its the same mass right ?)
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
          // Should add D0 bar ?
        } else { // If Lc
          invmass = o2::constants::physics::MassLambdaCPlus;
          // Should add Lc bar ? (maybe not its the same mass right ?)
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
              fillTpcTpcD0CandidateQa(multiplicity, track1);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-TPC Lc-h
              fillTpcTpcLcCandidateQa(multiplicity, track1);
            } else { // IF NEITHER D0 NOR LC -> TPC-TPC h-h
              fillTpcTpcChChSameEventQa(multiplicity, track1);
            }
          } else {                                                          // IF TPC-MFT case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-MFT D0-h
              fillTpcMftD0CandidateQa(multiplicity, track1);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-MFT Lc-h
              fillTpcMftLcCandidateQa(multiplicity, track1);
            } else { // IF NEITHER D0 NOR LC -> TPC-MFT h-h
              fillTpcMftChChSameEventQa(multiplicity, track1, true);
            } // end of if condition for TPC-TPC or TPC-MFT case
          }
          // Maybe I won't need it for MC (first files are way lighter in MC, but also I need to loop over all tracks in MC GEN)
        } else {                                                            // If MC (add cases later)
          if constexpr (!std::is_same_v<MftTracksMcLabels, TTracksAssoc>) { // IF TPC-TPC case
            fillTpcTpcChChSameEventQaMc(multiplicity, track1);
          }
        }
      }

      //
      // ASSOCIATED PARTICLE
      //
      for (const auto& track2 : tracks2) {

        // apply cuts for MFT tracks
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {
          registry.fill(HIST("Data/TpcMft/hMftTracksSelection"), MftTrackSelectionStep::NoSelection);

          if (!isAcceptedMftTrack(track2, true)) {
            continue;
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
              fillTpcTpcHfChSameEventQa(multiplicity, track2);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-TPC Lc-h
              fillTpcTpcHfChSameEventQa(multiplicity, track2);
            }
            // No if condition if it is h-h, because it would be the same plots than for the trigger particle
          } else {                                                          // IF TPC-MFT case
            if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig>) { // IF D0 CASE -> TPC-MFT D0-h
              fillTpcMftHfChSameEventQa(multiplicity, track2);
            } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF LC CASE -> TPC-MFT Lc-h
              fillTpcMftHfChSameEventQa(multiplicity, track2);
            } else { // IF NEITHER D0 NOR LC -> TPC-MFT h-h
              fillTpcMftChChSameEventQa(multiplicity, track2, false);
            } // end of if condition for TPC-TPC or TPC-MFT case
          }
          //} else {                                                        // If MC (add cases later)
          // fillTpcTpcChChSameEventQaMc(multiplicityTracks2, vz, tracks1);
          //}
        }

        if (sameEvent && (loopCounter == 1)) {
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
      //     // Should add D0 bar ?
      //   } else { // If Lc
      //     invmass = o2::constants::physics::MassLambdaCPlus;
      //     // Should add Lc bar ? (maybe not its the same mass right ?)
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
          fillTpcMftD0CandidateQa(multiplicity, track1);
        } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
          fillTpcTpcLcCandidateQa(multiplicity, track1);
        } else {
          fillTpcMftChChSameEventQa(multiplicity, track1, true);
        }
      }

      //
      // ASSOCIATED PARTICLE
      //
      for (const auto& track2 : tracks2) {

        // Fill QA plot for all MFT tracks () (only if cutAmbiguousTracks is false to avoid double counting)
        if (!cutAmbiguousTracks) {
          registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AllMftTracks);
        }

        // const auto& reassociatedMftTrack = track2.mfttrack();
        //  No one uses const and auto& here, so I will follow
        auto reassociatedMftTrack = track2.template mfttrack_as<aod::MFTTracks>();

        if (!isAcceptedMftTrack(reassociatedMftTrack, false)) {
          continue;
        }

        // Fill QA plot for MFT tracks after physical selection (eta + clusters)
        if (!cutAmbiguousTracks) {
          registry.fill(HIST("Data/TpcMft/hAmbiguityOfMftTracks"), MftTrackAmbiguityStep::AfterTrackSelection);
        }

        // We check if the track is ambiguous or non-ambiguous (QA plots are filled in isAmbiguousMftTrack)
        // Fill plots only if cutAmbiguousTracks is false (to avoid double counting)
        if (isAmbiguousMftTrack(track2, !cutAmbiguousTracks)) {
          // If the MFT track is ambiguous we may cut or not on the ambiguous track
          if (cutAmbiguousTracks) {
            continue;
          }
        }

        if (reassociatedMftTrack.collisionId() != track2.bestCollisionId()) {
          // track.collision_as<CollwEv>().posZ()
          continue;
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
            fillTpcMftHfChSameEventQa(multiplicity, reassociatedMftTrack);
          } else if constexpr (std::is_same_v<HfCandidatesSelLc, TTracksTrig>) {
            fillTpcMftHfChSameEventQa(multiplicity, reassociatedMftTrack);
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

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicityTracks1 = getPartsSize(collision1);

      if constexpr (std::is_same_v<FilteredCollisionsWSelMultMcLabels, TCollisions>) { // If MC
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
      } else {                                                                                                              // If not MC
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {                                                       // IF TPC-MFT case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-MFT HF-h
            registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing"), bin);
          } else { // IF h-h case -> TPC-MFT h-h case
            registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hEventCountMixing"), bin);
          }
        } else {                                                                                                            // IF TPC-TPC case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-TPC HF-h
            registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
          } else { // IF h-h case -> TPC-TPC h-h case
            registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
          }
        } // end of if condition for TPC-TPC or TPC-MFT case
      }

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
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
      } else {                                                                                                              // If not MC
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {                                              // IF TPC-MFT case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-MFT HF-h
            registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing"), bin);
          } else { // IF h-h case -> TPC-MFT h-h case
            registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hEventCountMixing"), bin);
          }
        } else {                                                                                                            // IF TPC-TPC case
          if constexpr (std::is_same_v<HfCandidatesSelD0, TTracksTrig> || std::is_same_v<HfCandidatesSelLc, TTracksTrig>) { // IF HF-h case -> TPC-TPC HF-h
            registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
          } else { // IF h-h case -> TPC-TPC h-h case
            registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
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

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, mcCollisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicity = getPartsSize(collision1); // get multiplicity of charged hadrons, which is used for slicing in mixing

      // TO BE DONE : ADD ONE MORE IF CONDITION TO FILL THE MC CASE
      // TODO : FILL NEW PLOTS FOR MCTRUTH ONLY
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);

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
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hMultiplicity"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hEventCountSame"), bin);

    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    // TO-DO : add if condition for when we will implant corrected correlations (kCFStepReconstructed -> kCFStepCorrected)
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
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEventCountSame"), bin);

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
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEventCountSame"), bin);

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
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEventCountSame"), bin);

    // I use kCFStepAll for running my code with all MFTTracks were the reassociation process was not applied
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, mftTracks, multiplicity, collision.posZ(), true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    // sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    // fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed,  tracks, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);

    // I use kCFStepTracked for running my code with only non-ambiguous MFTTracks
    // This is the same as running with reassociatedMftTracks, but applying one more cut in the fillCorrelations function
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    // sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
    // fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepTracked,  tracks, reassociatedMftTracks, multiplicity, collision.posZ(), false, true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChCh, "DATA : Process same-event correlations for TPC-MFT h-h case", false);

  void processSameTpcMftChChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         FilteredTracksWDcaSel const& tracks,
                                         aod::MFTTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEventCountSame"), bin);

    // I use kCFStepAll for running my code with all MFTTracks were the reassociation process was not applied
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations(sameEvent, CorrelationContainer::CFStep::kCFStepAll, tracks, mftTracks, multiplicity, collision.posZ(), true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepReconstructed, tracks, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);

    // I use kCFStepTracked for running my code with only non-ambiguous MFTTracks
    // This is the same as running with reassociatedMftTracks, but applying one more cut in the fillCorrelations function
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
    fillCorrelationsReassociatedMftTracks(sameEvent, CorrelationContainer::CFStep::kCFStepTracked, tracks, reassociatedMftTracks, multiplicity, collision.posZ(), true, true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChChReassociated, "DATA : Process same-event correlations for TPC-MFT h-h case reassociated", false);

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
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hEventCountSame"), bin);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftD0Ch, "DATA : Process same-event correlations for TPC-MFT D0-h case", false);

  void processSameTpcMftD0ChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         HfCandidatesSelD0 const& candidates,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         aod::MFTTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEventCountSame"), bin);

    // I use kCFStepAll for running my code with all MFTTracks were the reassociation process was not applied
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepAll, candidates, mftTracks, multiplicity, collision.posZ(), true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);

    // I use kCFStepTracked for running my code with only non-ambiguous MFTTracks
    // This is the same as running with reassociatedMftTracks, but applying one more cut in the fillCorrelations function
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
    fillCorrelationsReassociatedMftTracks(sameEventHf, CorrelationContainer::CFStep::kCFStepTracked, candidates, reassociatedMftTracks, multiplicity, collision.posZ(), true, true);
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
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hEventCountSame"), bin);

    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftLcCh, "DATA : Process same-event correlations for TPC-MFT Lc-h case", false);

  void processSameTpcMftLcChReassociated(FilteredCollisionsWSelMult::iterator const& collision,
                                         HfCandidatesSelLc const& candidates,
                                         soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                         aod::MFTTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return; // when process function has iterator
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    // registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEventCountSame"), bin);

    // I use kCFStepAll for running my code with all MFTTracks were the reassociation process was not applied
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations(sameEventHf, CorrelationContainer::CFStep::kCFStepAll, candidates, mftTracks, multiplicity, collision.posZ(), true);

    // I use the step kCFStepReconstructed for reassociatedMftTracks (most likely the ones we will use in the end)
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReassociatedMftTracks(sameEventHf, CorrelationContainer::CFStep::kCFStepReconstructed, candidates, reassociatedMftTracks, multiplicity, collision.posZ(), true, false);

    // I use kCFStepTracked for running my code with only non-ambiguous MFTTracks
    // This is the same as running with reassociatedMftTracks, but applying one more cut in the fillCorrelations function
    // We don't fill "normal" QA plots with these tracks, only specific plots to compare with other type of MFTTracks
    sameEventHf->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
    fillCorrelationsReassociatedMftTracks(sameEventHf, CorrelationContainer::CFStep::kCFStepTracked, candidates, reassociatedMftTracks, multiplicity, collision.posZ(), true, true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftLcChReassociated, "DATA : Process same-event correlations for TPC-MFT Lc-h case reassociated", false);

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
