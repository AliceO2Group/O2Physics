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

#include <string>
#include <vector>

#include <TDirectory.h>
#include <TH1F.h>
#include <THn.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskFlow {

  //  configurables for processing options
  Configurable<bool> doReferenceFlow{"doReferenceFlow", false, "Flag to know if reference flow should be done"};
  Configurable<bool> processRun2{"processRun2", false, "Flag to run on Run 2 data"};
  Configurable<bool> processRun3{"processRun3", true, "Flag to run on Run 3 data"};
  Configurable<bool> processMc{"processMc", false, "Flag to run on MC"};
  Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per event"};
  //  configurables for collisions
  Configurable<float> zVertexMax{"zVertexMax", 7.0f, "Accepted z-vertex range"};
  //  configurables for TPC tracks
  Configurable<float> etaTpcTrackMax{"etaTpcTrackMax", 0.8f, "max. eta of TPC tracks"};
  Configurable<float> ptTpcTrackMin{"ptTpcTrackMin", 0.5f, "min. pT of TPC tracks"};
  //  configurables for HF candidates
  Configurable<float> etaCandidateMax{"etaCandidateMax", 0.8f, "max. eta of HF candidate"};
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for LambdaC"};
  Configurable<int> selectionFlagLcToPiKP{"selectionFlagLcToPiKP", 1, "Selection Flag for LambdaC bar"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  //  configurables for MFT tracks
  Configurable<double> etaMftTrackMax{"etaMftTrackMax", 0, "Maximum value for the eta of MFT tracks"};
  Configurable<double> etaMftTrackMin{"etaMftTrackMin", -5, "Minimum value for the eta of MFT tracks"};
  Configurable<int> nClustersMftTrack{"nClustersMftTrack", 5, "Minimum number of clusters for the reconstruction of MFT tracks"};

  HfHelper hfHelper;
  SliceCache cache;
  Service<o2::framework::O2DatabasePDG> pdg;

  // =========================
  //      using declarations : DATA
  // =========================

  using FilteredCollisionsWSelMult = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using HfCandidatesSelD0 = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using HfCandidatesSelLc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using TracksWDcaSel = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection>>;

  // =========================
  //      using declarations : MC
  // =========================

  // Even add McCollisions in the join ?
  // Kata adds subscribes to it but do not add it in the join
  // using FilteredCollisionsWSelMultMC = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::McCollisions>>;
  using FilteredCollisionsWSelMultMC = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults>>;
  using FilteredMcCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::MultMCExtras, aod::McCollsExtra>>;
  using FilteredMcParticles = soa::Filtered<aod::McParticles>;
  using TracksWDcaSelMC = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::McTrackLabels>>;

  // Remnants, need Katarina's info
  // using FilteredCollisionsWDcaSelMcLabels = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>>;
  // using FilteredTracksWDcaSelMcLabels = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>>;

  // =========================
  //      Filters & partitions : DATA
  // =========================

  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilterD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0) ||
                             (aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);

  Filter candidateFilterLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi) ||
                             (aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPiKP);

  //  Collision filters
  //  FIXME: The filter is applied also on the candidates! Beware!
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < zVertexMax;

  Filter trackFilter = (nabs(aod::track::eta) < etaTpcTrackMax) &&
                       (aod::track::pt > ptTpcTrackMin) &&
                       requireGlobalTrackWoPtEtaInFilter();

  // =========================
  //      Filters & partitions : MC
  // =========================

  // From Katarina's code, but not sure if I use it
  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < zVertexMax;

  // From Katarina's code
  Filter mcParticlesFilter = (nabs(aod::mcparticle::eta) < etaTpcTrackMax) &&
                             (aod::mcparticle::pt > ptTpcTrackMin); //&&
                                                                    //(aod::mcparticle::sign != 0)

  // =========================
  //      Preslice : DATA
  // =========================

  Preslice<aod::Tracks> dataPerCol = aod::track::collisionId;

  // =========================
  //      Preslice : MC
  // =========================

  Preslice<aod::McParticles> mcTruthPerCol = aod::mcparticle::mcCollisionId;
  // Do I have to adapt this preslice to MC ? How does it work exactly ?
  // Preslice<aod::Tracks> mcRecPerCol = aod::track::collisionId;

  //  configurables for containers
  ConfigurableAxis axisVertex{"axisVertex", {14, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity axis for histograms"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  //  TODO: flow of HF will need to be done vs. invariant mass, in the signal and side-band regions
  //        either 1) add invariant mass axis or 2) define several containers for different inv. mass regions
  //        Note: don't forget to check inv. mass separately for D0 and D0bar candidate
  ConfigurableAxis axisMass{"axisMass", {2, 1.7, 2.0}, "axis of invariant mass of HF candidates"};
  ConfigurableAxis binsMixingVertex{"binsMixingVertex", {14, -7, 7}, "vertex bins for event mixing"};
  ConfigurableAxis binsMixingMultiplicity{"binsMixingMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity bins for event mixing"};

  HistogramRegistry registry{"registry"};

  // Correlation containers used for data
  OutputObj<CorrelationContainer> sameEvent{"sameEvent"};
  OutputObj<CorrelationContainer> mixedEvent{"mixedEvent"};
  OutputObj<CorrelationContainer> sameEventHf{"sameEventHf"};
  OutputObj<CorrelationContainer> mixedEventHf{"mixedEventHf"};
  // OutputObj<CorrelationContainer> sameTPCTPCChCh{"sameTPCTPCChCh"};
  // OutputObj<CorrelationContainer> mixedTPCTPCChCh{"mixedTPCTPCChCh"};
  // OutputObj<CorrelationContainer> sameTPCTPCHfCh{"sameTPCTPCHfCh"};   // I still keep only one Correlation Container for HF, whether is D0 or Lc
  // OutputObj<CorrelationContainer> mixedTPCTPCHfCh{"mixedTPCTPCHfCh"}; // Because only one should be run at the same time
  // OutputObj<CorrelationContainer> sameTPCMFTChCh{"sameTPCMFTChCh"};
  // OutputObj<CorrelationContainer> mixedTPCMFTChCh{"mixedTPCMFTChCh"};
  // OutputObj<CorrelationContainer> sameTPCMFTHfCh{"sameTPCMFTHfCh"};   // I still keep only one Correlation Container for HF, whether is D0 or Lc
  // OutputObj<CorrelationContainer> mixedTPCMFTHfCh{"mixedTPCMFTHfCh"}; // Because only one should be run at the same time

  // Correlation containers used for Monte-Carlo
  // OutputObj<CorrelationContainer> sameTPCTPCChChMC{"sameTPCTPCChChMC"};
  // OutputObj<CorrelationContainer> mixedTPCTPCChChMC{"mixedTPCTPCChChMC"};

  //  =========================
  //      init()
  //  =========================
  void init(InitContext&)
  {
    //  =========================
    //      Event histograms
    // TO-DO : do i have to separate event histograms between DATA and MC ?
    //  =========================
    constexpr int kNBinsEvents = 2;
    registry.add("Data/hEventCounter", "hEventCounter", {HistType::kTH1F, {{kNBinsEvents, 0.5, 0.5 + kNBinsEvents}}});
    //  set axes of the event counter histogram
    std::string labels[kNBinsEvents];
    labels[0] = "all";
    labels[1] = "after Physics selection";

    const int nBinsMix = axisMultiplicity->size() * 14; // 14 bins for z-vertex

    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("Data/hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    //  =========================
    //      DATA : histograms for TPC-TPC h-h case
    //  =========================

    // DATA : event histograms for TPC-TPC h-h same event
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    // DATA : associated particles histograms for TPC-TPC h-h same event
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});

    // Katarina had this :
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hVzEta", "eta vs. Vz", {HistType::kTH2F, {{100, -4, 4, "#eta"}, {20, -10, 10, "Vz"}}});

    // DATA : event mixing histograms for TPC-TPC h-h mixed event
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    //  =========================
    //      DATA : histograms for TPC-TPC HF-h case for 2PRONG
    //  =========================

    // DATA : event histograms for TPC-TPC HF-h same event
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});

    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    // DATA : trigger particles (candidates) histograms for TPC-TPC h-h same event
    auto vbins = (std::vector<double>)binsPt;
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtCandidate", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEtaCandVsPt", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/2Prong/hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    //  =========================
    //      DATA : histograms for TPC-TPC HF-h case for 3PRONG
    //  ===================

    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPt", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMass", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPtVsMult", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {5000, 0., 10000.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/3Prong/hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    //  =========================
    //      DATA : histograms for TPC-MFT h-h case
    //  =========================

    // DATA : trigger particles (TPC tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiTPC", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaTPC", "etaTPC", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPhiTPC", "phiTPC", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPtTPC", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hYieldsTPC", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hMultiplicityTPC", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : associated particles (MFT tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPhiMFT", "phiMFT", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPtMFT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hYieldsMFT", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -5, 0, "#eta"}}});

    // DATA : histograms for TPC-MFT h-h event mixing for events QA
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    //  =========================
    //      DATA : histograms for TPC-MFT HF-h case FOR 2PRONG
    //  =========================

    // DATA : trigger particles (candidates) histograms for TPC-MFT HF-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaPhiCandidate", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaCandidate", "etaTPC", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPhiCandidate", "phiTPC", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hYieldsCandidate", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hMultiplicityCandidate", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : trigger particles (candidates) histograms for TPC-MFT HF-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtCandidate", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hEtaCandVsPt", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/2Prong/hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    // DATA : associated particles (MFT tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPhiMFT", "phiMFT", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPtMFT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hYieldsMFT", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -5, 0, "#eta"}}});

    // DATA : histograms for TPC-MFT h-h event mixing for events QA
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});

    //  =========================
    //      DATA : histograms for TPC-MFT HF-h case FOR 3PRONG
    //  =========================

    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hYieldsCandidate", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMultiplicityCandidate", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEtaPhiCandidate", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPt", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMass", "3-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPtVsMult", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; multiplicity", {HistType::kTH3F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {5000, 0., 10000.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/3Prong/hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    //  =========================
    //      MC : histograms for TPC-TPC h-h case
    //  =========================

    // MC reconstructed

    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    // Katarina had this :
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary", "hMultiplicityPrimary", {HistType::kTH1F, {{500, 0, 500}}});
    //  histograms for MC associated particles
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    //  histograms for MC particles in event mixing
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

    // MC Truth

    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    // Katarina had this :
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary", "hMultiplicityPrimary", {HistType::kTH1F, {{500, 0, 500}}});
    //  histograms for MC associated particles
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    //  histograms for MC particles in event mixing
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{nBinsMix + 2, -2.5, -0.5 + nBinsMix, "bin"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

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
    // sameTPCTPCChCh.setObject(new CorrelationContainer("sameTPCTPCChCh", "sameTPCTPCChCh", corrAxis, effAxis, {}));
    // mixedTPCTPCChCh.setObject(new CorrelationContainer("mixedTPCTPCChCh", "mixedTPCTPCChCh", corrAxis, effAxis, {}));
    // sameTPCTPCHfCh.setObject(new CorrelationContainer("sameTPCTPCHfCh", "sameTPCTPCHfCh", corrAxis, effAxis, userAxis));
    // mixedTPCTPCHfCh.setObject(new CorrelationContainer("mixedTPCTPCHfCh", "mixedTPCTPCHfCh", corrAxis, effAxis, userAxis));
    // sameTPCMFTChCh.setObject(new CorrelationContainer("sameTPCMFTChCh", "sameTPCMFTChCh", corrAxis, effAxis, {}));
    // mixedTPCMFTChCh.setObject(new CorrelationContainer("mixedTPCMFTChCh", "mixedTPCMFTChCh", corrAxis, effAxis, {}));
    // sameTPCMFTHfCh.setObject(new CorrelationContainer("sameTPCMFTHfCh", "sameTPCMFTHfCh", corrAxis, effAxis, userAxis));
    // mixedTPCMFTHfCh.setObject(new CorrelationContainer("mixedTPCMFTHfCh", "mixedTPCMFTHfCh", corrAxis, effAxis, userAxis));

    // initialization of correlation containes for monte-carlo
    // sameTPCTPCChChMC.setObject(new CorrelationContainer("sameTPCTPCChChMC", "sameTPCTPCChChMC", corrAxis, effAxis, {}));
    // mixedTPCTPCChChMC.setObject(new CorrelationContainer("mixedTPCTPCChChMC", "mixedTPCTPCChChMC", corrAxis, effAxis, {}));
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
    if constexpr (std::is_same_v<TracksWDcaSelMC, TTrack>) { // if MC Rec
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

    auto pt = candidate.pt();

    if (candidate.isSelD0() >= selectionFlagD0) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPt"), hfHelper.invMassD0ToPiK(candidate), pt);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMass"), hfHelper.invMassD0ToPiK(candidate));
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMassVsPt"), hfHelper.invMassD0barToKPi(candidate), pt);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMass"), hfHelper.invMassD0barToKPi(candidate));
    }

    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hMultiplicity"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEta"), candidate.eta());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPhi"), phi);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtCandidate"), pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hPtProng1"), candidate.ptProng1());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/2Prong/hEtaCandVsPt"), candidate.eta(), pt);
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

    if (candidate.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMass"), hfHelper.invMassLcToPKPi(candidate));
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPtVsMult"), hfHelper.invMassLcToPKPi(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPt"), hfHelper.invMassLcToPKPi(candidate), pt);
    }
    if (candidate.isSelLcToPiKP() >= selectionFlagLcToPiKP) {
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMass"), hfHelper.invMassLcToPiKP(candidate));
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPtVsMult"), hfHelper.invMassLcToPiKP(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMassVsPt"), hfHelper.invMassLcToPiKP(candidate), pt);
    }

    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hMultiplicity"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPt"), pt);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng0"), ptProng0);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng1"), ptProng1);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hPtProng2"), ptProng2);
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEta"), candidate.eta());
    registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/3Prong/hEtaVsPt"), candidate.eta(), pt);
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

    if (candidate.isSelD0() >= selectionFlagD0) {
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMassVsPt"), hfHelper.invMassD0ToPiK(candidate), pt);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/2Prong/hMass"), hfHelper.invMassD0ToPiK(candidate));
    }
    if (candidate.isSelD0bar() >= selectionFlagD0bar) {
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

    if (candidate.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMass"), hfHelper.invMassLcToPKPi(candidate));
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPtVsMult"), hfHelper.invMassLcToPKPi(candidate), pt, multiplicity);
      registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/3Prong/hMassVsPt"), hfHelper.invMassLcToPKPi(candidate), pt);
    }
    if (candidate.isSelLcToPiKP() >= selectionFlagLcToPiKP) {
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
  //      Cuts with functions
  // =========================

  //  FIXME: Some collisions are rejected here, what causes (part of) differences with the D0 task
  template <typename TCollision>
  bool isAcceptedCollision(TCollision const& collision, bool fillHistograms = false)
  {
    if (fillHistograms) {
      registry.fill(HIST("Data/hEventCounter"), 1);
    }

    if (processMc == false) {
      if (!collision.sel8()) {
        return false;
      }
    }

    if (fillHistograms) {
      registry.fill(HIST("Data/hEventCounter"), 2);
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
      return true;
    } else { // For now, that means we do D0
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        return false;
      }
      if (etaCandidateMax >= 0. && std::abs(etaCandidate) > etaCandidateMax) {
        return false;
      }
      return true;
    }
  }

  //  TODO: Check how to put this into a Filter
  // I tried to put it as a filter, but filters for normal TPC tracks also apply to MFT tracks I think
  // and it seems that they are not compatible
  template <typename TTrack>
  bool isAcceptedMftTrack(TTrack const& mftTrack)
  {
    // cut on the eta of MFT tracks
    if (mftTrack.eta() > etaMftTrackMax || mftTrack.eta() < etaMftTrackMin) {
      return false;
    }

    // cut on the number of clusters of the reconstructed MFT track
    if (mftTrack.nClusters() < nClustersMftTrack) {
      return false;
    }

    return true;
  }

  // I am not sure if to template McParticles is useful, I'll address this when doing the MC Gen case of HF-h correlations
  template <CorrelationContainer::CFStep step, typename TMcParticles>
  bool isMcParticleSelected(TMcParticles& mcParticles)
  {
    //  remove MC particles with charge = 0
    TParticlePDG* pdgparticle = pdg->GetParticle(mcParticles.pdgCode());
    if (pdgparticle != nullptr) {
      if (pdgparticle->Charge() == 0) {
        return false;
      }
    }

    //  MC particle has to be primary
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
      return mcParticles.isPhysicalPrimary();
    }
    return true;
  }

  // =========================
  //      Correlation functions
  // =========================

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelations(TTarget target, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, float multiplicity, float posZ, bool sameEvent)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    // To avoid filling associated tracks QA many times
    //  I fill it only for the first trigger track of the collision
    auto loopCounter = 0;

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
        //  TODO: Check how to put this into a Filter
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

      // From Katarina's code
      //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
      // if (processMc) {
      // NOTE : this version with FilteredMcParticles is only for MC truth
      if constexpr (std::is_same_v<FilteredMcParticles, TTracksTrig> || std::is_same_v<FilteredMcParticles, TTracksAssoc>) {
        if (!isMcParticleSelected<step>(track1)) {
          continue;
        }
        // TO-DO : add other if constexpr conditions when I will have more MC cases
      }

      //  fill single-track distributions
      if (!fillingHFcontainer) { // if not HF-h case
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(step, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      // FILL QA PLOTS for trigger particle
      if (sameEvent) {
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
        } else {                                                         // If MC (add cases later)
          if constexpr (!std::is_same_v<aod::MFTTracks, TTracksAssoc>) { // IF TPC-TPC case
            fillTpcTpcChChSameEventQaMc(multiplicity, track1);
          }
        }
      }

      for (const auto& track2 : tracks2) {

        // apply cuts for MFT tracks
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {
          if (!isAcceptedMftTrack(track2)) {
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
        if constexpr (std::is_same_v<FilteredMcParticles, TTracksTrig> || std::is_same_v<FilteredMcParticles, TTracksAssoc>) {
          if (!isMcParticleSelected<step>(track2)) {
            continue;
          }
          // Note : no need for HF if condition as this will always be normal track, but maybe for MFT
        }

        float eta2 = track2.eta();
        float pt2 = track2.pt();
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
        if (sameEvent && (loopCounter == 1)) {
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

      } // end of loop over tracks2
    } // end of loop over tracks 1
  }

  template <typename TCollisions, typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisions(TCollisions const& collisions, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
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

      if constexpr (std::is_same_v<FilteredCollisionsWSelMultMC, TCollisions>) { // If MC
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

      corrContainer->fillEvent(multiplicityTracks1, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(corrContainer, tracks1, tracks2, multiplicityTracks1, collision1.posZ(), false);
    }
  }

  // template <typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  // void mixCollisions(FilteredCollisionsWSelMult const& collisions, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  template <typename TCollisions, typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisionsMcTruth(TCollisions const& collisions, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  {
    using BinningTypeMcTruth = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::mccollision::PosZ, decltype(getPartsSize)>;

    BinningTypeMcTruth binningWithTracksSize{{getPartsSize}, {axisVertex, axisMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<TCollisions, TTracksTrig, TTracksAssoc, BinningTypeMcTruth> pair{binningWithTracksSize, nMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicity = tracks2.size(); // get multiplicity of charged hadrons, which is used for slicing in mixing
      // const auto vz = collision1.posZ();

      // TO BE DONE : ADD ONE MORE IF CONDITION TO FILL THE MC CASE
      // TODO : FILL NEW PLOTS FOR MCTRUTH ONLY
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);

      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
      fillCorrelations<CorrelationContainer::kCFStepAll>(corrContainer, tracks1, tracks2, multiplicity, collision1.posZ(), false);
    }
  }

  // =====================================
  //    DATA : process same event correlations: TPC-TPC h-h case
  // =====================================

  void processSameTpcTpcChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    //  the event histograms below are only filled for h-h case
    //  because there is a possibility of double-filling if more correlation
    //  options are ran at the same time
    //  temporary solution, since other correlation options always have to be ran with h-h, too
    //  TODO: rewrite it in a more intelligent way
    // const auto multiplicity = tracks.size();
    const auto multiplicity = collision.multNTracksPV();
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hMultiplicity"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hEventCountSame"), bin);

    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    // TO-DO : add if condition for when we will implant corrected correlations (kCFStepReconstructed -> kCFStepCorrected)
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameEvent, tracks, tracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcChCh, "DATA : Process same-event correlations for TPC-TPC h-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-TPC HF-h case for D0
  // =====================================

  void processSameTpcTpcD0Ch(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks,
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

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameEventHf, candidates, tracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcD0Ch, "DATA : Process same-event correlations for TPC-TPC D0-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-TPC HF-h case for Lc
  // =====================================

  void processSameTpcTpcLcCh(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks,
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

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameEventHf, candidates, tracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcLcCh, "DATA : Process same-event correlations for TPC-TPC Lc-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT h-h case
  // =====================================

  void processSameTpcMftChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks,
                             aod::MFTTracks const& mftTracks)
  {
    if (!(isAcceptedCollision(collision, true))) {
      return;
    }

    const auto multiplicity = collision.multNTracksPV();
    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEventCountSame"), bin);

    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameEvent, tracks, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChCh, "DATA : Process same-event correlations for TPC-MFT h-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT HF-h case for D0
  // =====================================

  void processSameTpcMftD0Ch(FilteredCollisionsWSelMult::iterator const& collision,
                             HfCandidatesSelD0 const& candidates,
                             TracksWDcaSel const& /*tracks*/,
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

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameEventHf, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftD0Ch, "DATA : Process same-event correlations for TPC-MFT D0-h case", false);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT HF-h case for Lc
  // =====================================

  void processSameTpcMftLcCh(FilteredCollisionsWSelMult::iterator const& collision,
                             HfCandidatesSelLc const& candidates,
                             TracksWDcaSel const& /*tracks*/,
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

    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameEventHf, candidates, mftTracks, multiplicity, collision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftLcCh, "DATA : Process same-event correlations for TPC-MFT Lc-h case", false);

  // =====================================
  //    MONTE-CARLO : process same event correlations: TPC-TPC h-h case
  // =====================================

  // Katarina's version = MC Truth
  void processSameTpcTpcChChmcGEN(FilteredMcCollisions::iterator const& mcCollision,
                                  FilteredMcParticles const& mcParticles)
  {

    // Not sure why to use this
    // if (collisions.size() == 0) {
    //  return;
    //}

    // if (!collision.has_mcCollision()) {
    //   LOGF(info, "No MC collision for this collision, skip...");
    //   return;
    // }

    // TODO : check if I have to get my multiplicity based on multNTracksPV or mcParticles.size()
    // const auto multiplicity = mcParticles.size(); // Note: these are all MC particles after selection (not only primary)
    const auto multiplicity = mcCollision.multMCPVz();
    registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicity"), multiplicity);
    registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hVtxZ"), mcCollision.posZ());

    //  fill correlations for all MC collisions
    // In Katka's code, the first time doing this does not fill the histograms, right now will be filled two times..
    // auto multPrimaryCharge0 = fillTpcTpcChChSameEventQaMc<CorrelationContainer::kCFStepAll>(multiplicity, mcParticles);
    sameEvent->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(sameEvent, mcParticles, mcParticles, multiplicity, mcCollision.posZ(), true);

    // TO-DO : fill correlation container for MC collisions that have a reconstructed collision
    // got rid of the second const auto for multPrimaryCharge0
    // This line below for sure induce that some plots are filled two times
    // multPrimaryCharge0 = fillTpcTpcChChSameEventQaMc<CorrelationContainer::kCFStepVertex>(multiplicity, mcParticles);
    // sameTPCTPCChChMC->fillEvent(multPrimaryCharge0, CorrelationContainer::kCFStepVertex);
    // fillCorrelations<CorrelationContainer::kCFStepVertex>(sameTPCTPCChChMC, mcParticles, mcParticles, multPrimaryCharge0, mcCollision.posZ(), true);
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcChChmcGEN, "MONTE-CARLO : Process same-event correlations for TPC-TPC h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations:TPC-TPC h-h case
  // =====================================
  // TO BECOME DATA & MC REC ?

  void processMixedTpcTpcChCh(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks)
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
    mixCollisions(collisions, tracks, tracks, getMultiplicity, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChCh, "DATA : Process mixed-event correlations for TPC-TPC h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-TPC HF-h case for D0
  // =====================================

  void processMixedTpcTpcD0Ch(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks,
                              HfCandidatesSelD0 const& candidates)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, candidates, tracks, getMultiplicity, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcD0Ch, "DATA : Process mixed-event correlations for TPC-TPC D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-TPC HF-h case for Lc
  // =====================================

  void processMixedTpcTpcLcCh(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks,
                              HfCandidatesSelLc const& candidates)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks, this](FilteredCollisionsWSelMult::iterator const& col) {
      // Still o2::aod::track::collisionId with HF ??? -> I don't think so
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache);
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, candidates, tracks, getTracksSize, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcLcCh, "DATA : Process mixed-event correlations for TPC-TPC Lc-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT h-h case
  // =====================================

  void processMixedTpcMftChCh(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks,
                              aod::MFTTracks const& mftTracks)
  {
    //  we want to group collisions based on charged-track multiplicity
    // auto getTracksSize = [&tracks, this](FilteredCollisionsWSelMult::iterator const& col) {
    //   auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache);
    //   auto size = associatedTracks.size();
    //   return size;
    // };

    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, tracks, mftTracks, getMultiplicity, mixedEvent);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftChCh, "DATA : Process mixed-event correlations for TPC-MFT h-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT HF-h case for D0
  // =====================================

  void processMixedTpcMftD0Ch(FilteredCollisionsWSelMult const& collisions,
                              HfCandidatesSelD0 const& candidates,
                              aod::MFTTracks const& mftTracks,
                              TracksWDcaSel const& /*tracks*/)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, candidates, mftTracks, getMultiplicity, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftD0Ch, "DATA : Process mixed-event correlations for TPC-MFT D0-h case", false);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT HF-h case
  // =====================================

  void processMixedTpcMftLcCh(FilteredCollisionsWSelMult const& collisions,
                              HfCandidatesSelLc const& candidates,
                              aod::MFTTracks const& mftTracks)
  {

    //  we want to group collisions based on charged-track multiplicity
    auto getMultiplicity = [](FilteredCollisionsWSelMult::iterator const& collision) {
      auto multiplicity = collision.numContrib();
      return multiplicity;
    };

    mixCollisions(collisions, candidates, mftTracks, getMultiplicity, mixedEventHf);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftLcCh, "DATA : Process mixed-event correlations for TPC-MFT Lc-h case", false);

  // =====================================
  //    MONTE-CARLO : process mixed event correlations: TPC-TPC h-h case
  // =====================================

  // MC gen
  void processMixedTpcTpcChChmcGEN(FilteredMcCollisions const& mcCollisions,
                                   FilteredMcParticles const& mcParticles)
  {
    // use normal index instead of globalIndex for MixedEvent ??

    //  we want to group collisions based on charged-track multiplicity
    // auto getTracksSize = [&mcParticles, this](FilteredMcCollisions::iterator const& mcCol) {
    //  auto associatedTracks = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId, mcCol.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
    //  auto size = associatedTracks.size();
    //  return size;
    //};

    auto getMultiplicity = [](FilteredMcCollisions::iterator const& mcCollision) {
      auto multiplicity = mcCollision.multMCPVz();
      return multiplicity;
    };

    mixCollisionsMcTruth(mcCollisions, mcParticles, mcParticles, getMultiplicity, mixedEvent);

    // TO-DO : mixed event for particles that have a reconstructed collision kCFStepVertex
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChChmcGEN, "MONTE-CARLO : Process mixed-event correlations for TPC-TPC h-h case", false);
}; // End of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskFlow>(cfgc)};
}
