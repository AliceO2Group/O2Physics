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
/// \author Katarina Krizkova Gajdosova <katarina.gajdosova@cern.ch>, CERN
/// \author Maja Kabus <maja.kabus@cern.ch>, CERN

#include <iostream>

#include <TDirectory.h>
#include <TH1F.h>
#include <THn.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
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
  Service<o2::framework::O2DatabasePDG> pdg;

  //  configurables for processing options
  Configurable<bool> processRun2{"processRun2", false, "Flag to run on Run 2 data"};
  Configurable<bool> processRun3{"processRun3", true, "Flag to run on Run 3 data"};
  Configurable<bool> processMc{"processMc", false, "Flag to run on MC"};
  Configurable<int> nMixedEvents{"nMixedEvents", 5, "Number of mixed events per event"};
  //  configurables for collisions
  Configurable<float> zVertexMax{"zVertexMax", 7.0f, "Accepted z-vertex range"};
  //  configurables for associated particles
  Configurable<float> etaTrackAssocMax{"etaTrackAssocMax", 0.8f, "max. eta of associated tracks"};
  Configurable<float> ptTrackAssocMin{"ptTrackAssocMin", 0.5f, "min. pT of associated tracks"};
  //  configurables for HF candidates
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;
  SliceCache cache;

  // =========================
  //      DATA : Filters & partitions
  // =========================

  //  Collision filters
  //  FIXME: The filter is applied also on the candidates! Beware!
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < zVertexMax;

  using FilteredCollisionsWSelMult = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using TracksWDcaSel = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection>>;
  using HfCandidatesSel = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;

  Filter trackFilter = (nabs(aod::track::eta) < etaTrackAssocMax) &&
                       (aod::track::pt > ptTrackAssocMin) &&
                       requireGlobalTrackWoPtEtaInFilter();

  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilter = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 ||
                           aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  Preslice<aod::Tracks> dataPerCol = aod::track::collisionId;
  // Do I have to adapt this preslice to MC ? How does it work exactly ?
  // Preslice<aod::Tracks> mcRecPerCol = aod::track::collisionId;
  Preslice<aod::McParticles> mcTruthPerCol = aod::mcparticle::mcCollisionId;

  // =========================
  //      MONTE-CARLO : Filters & partitions
  // =========================

  // Even add McCollisions in the join ?
  // Kata adds subscribes to it but do not add it in the join
  // using FilteredCollisionsWSelMultMC = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults, aod::McCollisions>>;
  using FilteredCollisionsWSelMultMC = soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::Mults>>;
  using TracksWDcaSelMC = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection, aod::McTrackLabels>>;

  // using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>>;

  // From Katarina's code, but not sure if I use it
  Filter mcCollisionFilter = nabs(aod::mccollision::posZ) < zVertexMax;
  using aodMcCollisions = soa::Filtered<aod::McCollisions>;

  // From Katarina's code
  Filter mcParticlesFilter = (nabs(aod::mcparticle::eta) < etaTrackAssocMax) &&
                             (aod::mcparticle::pt > ptTrackAssocMin); //&&
                                                                      //(aod::mcparticle::sign != 0)
  using aodMcParticles = soa::Filtered<aod::McParticles>;

  // Katarina had this in her code :
  //  Charged track filters
  /*
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackAssocMax) &&
                       (aod::track::pt > ptTrackAssocMin) &&
                       requireGlobalTrackWoPtEtaInFilter();
  */

  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>>;

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

  HistogramRegistry registry{"registry"};

  // Correlation containers used for data
  OutputObj<CorrelationContainer> sameTPCTPCChCh{"sameTPCTPCChCh"};
  OutputObj<CorrelationContainer> mixedTPCTPCChCh{"mixedTPCTPCChCh"};
  OutputObj<CorrelationContainer> sameTPCTPCHfCh{"sameTPCTPCHfCh"};
  OutputObj<CorrelationContainer> mixedTPCTPCHfCh{"mixedTPCTPCHfCh"};
  OutputObj<CorrelationContainer> sameTPCMFTChCh{"sameTPCMFTChCh"};
  OutputObj<CorrelationContainer> mixedTPCMFTChCh{"mixedTPCMFTChCh"};
  OutputObj<CorrelationContainer> sameTPCMFTHfCh{"sameTPCMFTHfCh"};
  OutputObj<CorrelationContainer> mixedTPCMFTHfCh{"mixedTPCMFTHfCh"};

  // Correlation containers used for Monte-Carlo
  OutputObj<CorrelationContainer> sameTPCTPCChChMC{"sameTPCTPCChChMC"};
  OutputObj<CorrelationContainer> mixedTPCTPCChChMC{"mixedTPCTPCChChMC"};

  //  =========================
  //      init()
  //  =========================
  void init(InitContext&)
  {
    //  =========================
    //      Event histograms
    // TO-DO : do i have to separate event histograms between DATA and MC ?
    //  =========================
    constexpr int kNBinsEvents = 3;
    registry.add("Data/hEventCounter", "hEventCounter", {HistType::kTH1F, {{kNBinsEvents, 0.5, 0.5 + kNBinsEvents}}});
    //  set axes of the event counter histogram
    std::string labels[kNBinsEvents];
    labels[0] = "all";
    labels[1] = "after trigger selection (Run 2)";
    labels[2] = "after Physics selection";

    const int maxMixBin = axisMultiplicity->size() * 14; // 14 bins for z-vertex

    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("Data/hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }

    //  =========================
    //      DATA : histograms for TPC-TPC h-h case
    //  =========================

    // DATA : event histograms for TPC-TPC h-h same event
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});

    // DATA : associated particles histograms for TPC-TPC h-h same event
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    // Katarina had this :
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hVzEta", "eta vs. Vz", {HistType::kTH2F, {{100, -4, 4, "#eta"}, {20, -10, 10, "Vz"}}});

    // DATA : event mixing histograms for TPC-TPC h-h mixed event
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : particles histograms for TPC-TPC h-h mixed event
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

    //  =========================
    //      DATA : histograms for TPC-TPC HF-h case
    //  =========================

    // DATA : event histograms for TPC-TPC HF-h same event
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hMultiplicityHFMixing", "hMultiplicityHFMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hVtxZHFMixing", "hVtxZHFMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hNtracksHFMixing", "hNtracksHFMixing", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : trigger particles (candidates) histograms for TPC-TPC h-h same event
    auto vbins = (std::vector<double>)binsPt;
    registry.add("Data/TpcTpc/HfHadron/hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcTpc/HfHadron/hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcTpc/HfHadron/hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcTpc/HfHadron/hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hDecLengthXY", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hEtaCand", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcTpc/HfHadron/hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    // DATA : trigger particles (candidates) histograms for TPC-TPC h-h mixed event
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hPtHFMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hEtaHFMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hPhiHFMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

    //  =========================
    //      DATA : histograms for TPC-MFT h-h case
    //  =========================

    // DATA : trigger particles (TPC tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiTPC", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaTPC", "etaTPC", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPhiTPC", "phiTPC", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPtTPC", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hYieldsTPC", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hNtracksTPC", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : associated particles (MFT tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hEtaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPhiMFT", "phiMFT", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hPtMFT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HadronHadron/SameEvent/hYieldsMFT", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hNtracksMFT", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});


    // DATA : histograms for TPC-MFT h-h event mixing for TPC tracks
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hMultiplicityMixingTPC", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hVtxZMixingTPC", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hPtMixingTPC", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hEtaMixingTPC", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hPhiMixingTPC", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hNtracksMixingTPC", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : histograms for TPC-MFT h-h event mixing for MFT tracks
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hMultiplicityMixingMFT", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hVtxZMixingMFT", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hPtMixingMFT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hEtaMixingMFT", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hPhiMixingMFT", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hNtracksMixingMFT", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : histograms for TPC-MFT h-h event mixing for events QA
    registry.add("Data/TpcMft/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});

    //  =========================
    //      DATA : histograms for TPC-MFT HF-h case
    //  =========================

    // DATA : trigger particles (candidates) histograms for TPC-MFT HF-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaPhiCandidate", "multiplicity vs eta vs phi in TPC", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaCandidate", "etaTPC", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPhiCandidate", "phiTPC", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPtCandidate", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hYieldsCandidate", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hNtracksCandidate", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : trigger particles (candidates) histograms for TPC-TPC h-h same event
    registry.add("Data/TpcMft/HfHadron/hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcMft/HfHadron/hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcMft/HfHadron/hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("Data/TpcMft/HfHadron/hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hDecLengthXY", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hEtaCand", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/TpcMft/HfHadron/hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    // DATA : associated particles (MFT tracks) histograms for TPC-MFT h-h same event
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hEtaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPhiMFT", "phiMFT", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hPtMFT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HfHadron/SameEvent/hYieldsMFT", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HfHadron/SameEvent/hNtracksMFT", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});


    // DATA : histograms for TPC-MFT h-h event mixing for candidates
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hMultiplicityMixingCandidate", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hVtxZMixingCandidate", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hPtMixingCandidate", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hEtaMixingCandidate", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hPhiMixingCandidate", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hNtracksMixingCandidate", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : histograms for TPC-MFT h-h event mixing for MFT tracks
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hMultiplicityMixingMFT", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hVtxZMixingMFT", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hPtMixingMFT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hEtaMixingMFT", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hPhiMixingMFT", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hNtracksMixingMFT", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});

    // DATA : histograms for TPC-MFT h-h event mixing for events QA
    registry.add("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});

    //  =========================
    //      MC : histograms for TPC-MFT h-h case
    //  =========================

    // MC reconstructed

    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});
    // Katarina had this :
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary", "hMultiplicityPrimary", {HistType::kTH1F, {{500, 0, 500}}});
    //  histograms for MC associated particles
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/SameEvent/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});
    //  histograms for MC particles in event mixing
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});

    // MC Truth

    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});
    // Katarina had this :
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary", "hMultiplicityPrimary", {HistType::kTH1F, {{500, 0, 500}}});
    //  histograms for MC associated particles
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/SameEvent/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});
    //  histograms for MC particles in event mixing
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});
    registry.add("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});

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
    sameTPCTPCChCh.setObject(new CorrelationContainer("sameTPCTPCChCh", "sameTPCTPCChCh", corrAxis, effAxis, {}));
    mixedTPCTPCChCh.setObject(new CorrelationContainer("mixedTPCTPCChCh", "mixedTPCTPCChCh", corrAxis, effAxis, {}));
    sameTPCTPCHfCh.setObject(new CorrelationContainer("sameTPCTPCHfCh", "sameTPCTPCHfCh", corrAxis, effAxis, userAxis));
    mixedTPCTPCHfCh.setObject(new CorrelationContainer("mixedTPCTPCHfCh", "mixedTPCTPCHfCh", corrAxis, effAxis, userAxis));
    sameTPCMFTChCh.setObject(new CorrelationContainer("sameTPCMFTChCh", "sameTPCMFTChCh", corrAxis, effAxis, {}));
    mixedTPCMFTChCh.setObject(new CorrelationContainer("mixedTPCMFTChCh", "mixedTPCMFTChCh", corrAxis, effAxis, {}));
    sameTPCMFTHfCh.setObject(new CorrelationContainer("sameTPCMFTHfCh", "sameTPCMFTHfCh", corrAxis, effAxis, {}));
    mixedTPCMFTHfCh.setObject(new CorrelationContainer("mixedTPCMFTHfCh", "mixedTPCMFTHfCh", corrAxis, effAxis, {}));

    // initialization of correlation containes for monte-carlo
    sameTPCTPCChChMC.setObject(new CorrelationContainer("sameTPCTPCChChMC", "sameTPCTPCChChMC", corrAxis, effAxis, userAxis));
    mixedTPCTPCChChMC.setObject(new CorrelationContainer("mixedTPCTPCChChMC", "mixedTPCTPCChChMC", corrAxis, effAxis, userAxis));
  } // End of init() function

  // =========================
  //    templates
  //  FIXME: Some collisions are rejected here, what causes (part of) differences with the D0 task
  // =========================
  template <typename TCollision>
  bool isCollisionSelected(TCollision const& collision, bool fillHistograms = false)
  {
    if (fillHistograms)
      registry.fill(HIST("Data/hEventCounter"), 1);

    if (processMc == false) {
      if (!collision.sel8()) {
        return false;
      }
    }

    if (fillHistograms)
      registry.fill(HIST("Data/hEventCounter"), 3);

    return true;
  }

  // =========================
  //      Quality Assesment plots
  // =========================

  // ---- DATA : TPC-TPC h-h Same Event QA ----
  template <typename TTracks>
  void fillTpcTpcChChSameEventQA(float multiplicity, TTracks const& tracks)
  {
    int Ntracks = 0;
    for (const auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hPt"), track1.pt());
      registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hEta"), track1.eta());
      registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hPhi"), track1.phi());
      registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hYields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hEtaPhi"), multiplicity, track1.eta(), track1.phi());
    }
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hNtracks"), Ntracks);
  }

  // ---- MC : TPC-TPC h-h Same Event QA ----
  // Changed quickly void for int type + added the return for a test
  template <CorrelationContainer::CFStep step, typename TTracks>
  int fillTpcTpcChChSameEventQAmc(float multiplicity, TTracks const& tracks)
  {
    int Ntracks = 0;
    for (const auto& track1 : tracks) {
      // TO-DO : DEBUG THIS

      //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
      if constexpr (std::is_same_v<aodMcParticles, TTracks>) {
        if (!isMcParticleSelected<step>(track1)) {
          continue;
        }
        // TO-DO : add other if constexpr conditions when I will have more MC cases
      }

      Ntracks++;
      if constexpr (std::is_same_v<TracksWDcaSelMC, TTracks>) { // if MC Rec
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPt"), track1.pt());
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEta"), track1.eta());
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hPhi"), track1.phi());
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hYields"), multiplicity, track1.pt(), track1.eta());
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEtaPhi"), multiplicity, track1.eta(), track1.phi());
      } else { // if MC Gen
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPt"), track1.pt());
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEta"), track1.eta());
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hPhi"), track1.phi());
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hYields"), multiplicity, track1.pt(), track1.eta());
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEtaPhi"), multiplicity, track1.eta(), track1.phi());
      }
    }
    if constexpr (std::is_same_v<TracksWDcaSelMC, TTracks>) { // if MC Rec
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hNtracks"), Ntracks);
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary"), Ntracks);
    } else { // if MC Gen
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hNtracks"), Ntracks);
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicityPrimary"), Ntracks);
    }
    return Ntracks;
  }

  // ---- DATA : TPC-TPC h-h Mixed Event QA ----
  template <typename TTracks>
  void fillTpcTpcChChMixedEventQA(float multiplicity, float vz, TTracks const& tracks)
  {
    registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing"), vz);

    int Ntracks = 0;
    for (const auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hPtMixing"), track1.pt());
      registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hEtaMixing"), track1.eta());
      registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hPhiMixing"), track1.phi());
    }
    registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing"), Ntracks);
  }

  // ---- MC : TPC-TPC h-h Mixed Event QA ----
  template <typename TTracks>
  void fillTpcTpcChChMixedEventQAmc(float multiplicity, float vz, TTracks const& tracks)
  {
    if constexpr (std::is_same_v<TracksWDcaSelMC, TTracks>) { // if MC Rec
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing"), multiplicity);
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing"), vz);
    } else { // if MC Gen
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing"), multiplicity);
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing"), vz);
    }

    int Ntracks = 0;
    for (const auto& track1 : tracks) {
      Ntracks++;
      if constexpr (std::is_same_v<TracksWDcaSelMC, TTracks>) { // if MC Rec
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPtMixing"), track1.pt());
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEtaMixing"), track1.eta());
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hPhiMixing"), track1.phi());
      } else { // if MC Gen
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPtMixing"), track1.pt());
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEtaMixing"), track1.eta());
        registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hPhiMixing"), track1.phi());
      }
    }
    if constexpr (std::is_same_v<TracksWDcaSelMC, TTracks>) { // if MC Rec
      registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing"), Ntracks);
    } else { // if MC Gen
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing"), Ntracks);
    }
  }

  // ---- DATA : TPC-TPC HF-h Mixed Event QA ----
  template <typename TTracks>
  void fillTpcTpcHfChMixedEventQA(float multiplicity, float vz, TTracks const& tracks)
  {
    registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hMultiplicityHFMixing"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hVtxZHFMixing"), vz);

    int Ntracks = 0;
    for (const auto& track1 : tracks) {
      Ntracks++;
      registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hPtHFMixing"), track1.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEtaHFMixing"), track1.eta());
      registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hPhiHFMixing"), track1.phi());
    }
    registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hNtracksHFMixing"), Ntracks);
  }

  // ---- DATA : TPC-MFT h-h Same Event QA ----
  template <typename TTracks>
  void fillTpcMftChChSameEventQA(float multiplicity, TTracks const& tracks)
  {
    int Ntracks = 0;
    bool isMFT = false;
    for (const auto& track1 : tracks) {
      if constexpr (std::is_same_v<aod::MFTTracks, TTracks>) { // if MFT tracks
        isMFT = true;
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaMFT"), track1.eta());
        float phi = track1.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPhiMFT"), phi);
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiMFT"), multiplicity, track1.eta(), phi);

        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPtMFT"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hYieldsMFT"), multiplicity, track1.pt(), track1.eta());
      }
      else {                                                  // if TPC tracks
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaTPC"), track1.eta());
        float phi = track1.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPhiTPC"), phi);
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hEtaPhiTPC"), multiplicity, track1.eta(), phi);

        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hPtTPC"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HadronHadron/SameEvent/hYieldsTPC"), multiplicity, track1.pt(), track1.eta());
      }
      if (isMFT) {
        registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hNtracksMFT"), Ntracks);
      } else {
        registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hNtracksTPC"), Ntracks);
      }
      
    }
  }

  // ---- DATA : TPC-MFT HF-h Same Event QA ----
  template <typename TTracks>
  void fillTpcMftHfChSameEventQA(float multiplicity, TTracks const& tracks)
  {
    int Ntracks = 0;
    bool isMFT = false;
    for (const auto& track1 : tracks) {
      if constexpr (std::is_same_v<aod::MFTTracks, TTracks>) { // if MFT tracks
        isMFT = true;
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hEtaMFT"), track1.eta());
        float phi = track1.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hPhiMFT"), phi);
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hEtaPhiMFT"), multiplicity, track1.eta(), phi);

        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hPtMFT"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hYieldsMFT"), multiplicity, track1.pt(), track1.eta());
      }
      else {                                                  // if TPC tracks
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hEtaCandidate"), track1.eta());
        float phi = track1.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hPhiCandidate"), phi);
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hEtaPhiCandidate"), multiplicity, track1.eta(), phi);

        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hPtCandidate"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HfHadron/SameEvent/hYieldsCandidate"), multiplicity, track1.pt(), track1.eta());
      }
      if (isMFT) {
        registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/hNtracksMFT"), Ntracks);
      } else {
        registry.fill(HIST("Data/TpcTpc/HfHadron/SameEvent/hNtracksCandidate"), Ntracks);
      }
      
    }
  }

  // ---- DATA : TPC-MFT h-h Mixed Event QA ----
  template <typename TTracks>
  void fillTpcMftChChMixedEventQA(float multiplicity, float vz, TTracks const& tracks)
  {
    if constexpr (std::is_same_v<aod::MFTTracks, TTracks>) { // if MFT tracks
      registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hMultiplicityMixingMFT"), multiplicity);
      registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hVtxZMixingMFT"), vz);

      int Ntracks = 0;
      for (const auto& track1 : tracks) {
        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hPtMixingMFT"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hEtaMixingMFT"), track1.eta());
        registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hPhiMixingMFT"), track1.phi());
      }
      registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hNtracksMixingMFT"), Ntracks);
    }
    else {                                                   // if TPC tracks
      registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hMultiplicityMixingTPC"), multiplicity);
      registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hVtxZMixingTPC"), vz);

      int Ntracks = 0;
      for (const auto& track1 : tracks) {
        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hPtMixingTPC"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hEtaMixingTPC"), track1.eta());
        registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hPhiMixingTPC"), track1.phi());
      }
      registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hNtracksMixingTPC"), Ntracks);
    }
  }

  // ---- DATA : TPC-MFT h-h Mixed Event QA ----
  template <typename TTracks>
  void fillTpcMftHfChMixedEventQA(float multiplicity, float vz, TTracks const& tracks)
  {
    if constexpr (std::is_same_v<aod::MFTTracks, TTracks>) { // if MFT tracks
      registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hMultiplicityMixingMFT"), multiplicity);
      registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hVtxZMixingMFT"), vz);

      int Ntracks = 0;
      for (const auto& track1 : tracks) {
        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hPtMixingMFT"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hEtaMixingMFT"), track1.eta());
        registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hPhiMixingMFT"), track1.phi());
      }
      registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hNtracksMixingMFT"), Ntracks);
    }
    else {                                                   // if candidate tracks
      registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hMultiplicityMixingCandidate"), multiplicity);
      registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hVtxZMixingCandidate"), vz);

      int Ntracks = 0;
      for (const auto& track1 : tracks) {
        Ntracks++;
        registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hPtMixingCandidate"), track1.pt());
        registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hEtaMixingCandidate"), track1.eta());
        registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hPhiMixingCandidate"), track1.phi());
      }
      registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hNtracksMixingCandidate"), Ntracks);
    }
  }

  //  TODO: Check how to put this into a Filter
  template <typename TTrack>
  bool isAcceptedCandidate(TTrack const& candidate)
  {
    if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
      return false;
    }
    if (yCandMax >= 0. && std::abs(hfHelper.yD0(candidate)) > yCandMax) {
      return false;
    }
    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTrack>
  bool isMcParticleSelected(TTrack& track)
  {
    //  remove MC particles with charge = 0
    int8_t sign = 0;
    TParticlePDG* pdgparticle = pdg->GetParticle(track.pdgCode());
    if (pdgparticle != nullptr) {
      sign = (pdgparticle->Charge() > 0) ? 1.0 : ((pdgparticle->Charge() < 0) ? -1.0 : 0.0);
    }
    if (sign == 0) {
      return false;
    }

    //  MC particle has to be primary
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
      return track.isPhysicalPrimary();
    }
    return true;
  }

  // ---- DATA : TPC-TPC HF-h Same Event (Candidates) QA ----
  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTracks>
  void fillTpcTpcCandidateQA(TTracks const& candidates)
  {
    for (const auto& candidate : candidates) {
      if (!isAcceptedCandidate(candidate)) {
        continue;
      }

      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("Data/TpcTpc/HfHadron/hMass"), hfHelper.invMassD0ToPiK(candidate), candidate.pt());
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("Data/TpcTpc/HfHadron/hMass"), hfHelper.invMassD0barToKPi(candidate), candidate.pt());
      }

      registry.fill(HIST("Data/TpcTpc/HfHadron/hPtCand"), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hCTS"), hfHelper.cosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hCt"), hfHelper.ctD0(candidate), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hEtaCand"), candidate.eta(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("Data/TpcTpc/HfHadron/hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    }
  }

  // ---- DATA : TPC-MFT HF-h Same Event (Candidates) QA ----
  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTracks>
  void fillTpcMftCandidateQA(TTracks const& candidates)
  {
    for (const auto& candidate : candidates) {
      if (!isAcceptedCandidate(candidate)) {
        continue;
      }

      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("Data/TpcMft/HfHadron/hMass"), hfHelper.invMassD0ToPiK(candidate), candidate.pt());
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("Data/TpcMft/HfHadron/hMass"), hfHelper.invMassD0barToKPi(candidate), candidate.pt());
      }

      registry.fill(HIST("Data/TpcMft/HfHadron/hPtCand"), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("Data/TpcMft/HfHadron/hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("Data/TpcMft/HfHadron/hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hDecLengthXY"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hCTS"), hfHelper.cosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hCt"), hfHelper.ctD0(candidate), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hEtaCand"), candidate.eta(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("Data/TpcMft/HfHadron/hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    }
  }

  // =========================
  //      Correlation functions
  // =========================

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelations(TTarget target, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, float multiplicity, float posZ)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    for (const auto& track1 : tracks1) {

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
      if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig>) {
        //  TODO: Check how to put this into a Filter
        if (!isAcceptedCandidate(track1)) {
          continue;
        }
        fillingHFcontainer = true;
        invmass = hfHelper.invMassD0ToPiK(track1);
      }

      // From Katarina's code
      //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
      // if (processMc) {
      // NOTE : this version with aodMcParticles is only for Katarina's way of doing MC
      if constexpr (std::is_same_v<aodMcParticles, TTracksTrig> || std::is_same_v<aodMcParticles, TTracksAssoc>) {
        if (!isMcParticleSelected<step>(track1)) {
          continue;
        }
        // TO-DO : add other if constexpr conditions when I will have more MC cases
      }

      //  fill single-track distributions
      if (!fillingHFcontainer) {
        target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, pt1, multiplicity, posZ, triggerWeight);
      } else {
        target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, pt1, multiplicity, posZ, invmass, triggerWeight);
      }

      for (const auto& track2 : tracks2) {

        //  case of h-h correlations where the two types of tracks are the same
        //  this avoids autocorrelations and double counting of particle pairs
        if constexpr (std::is_same_v<TTracksAssoc, TTracksTrig>) {
          if (track1.index() <= track2.index()) {
            continue;
          }
        }

        //  in case of HF-h correlations, remove candidate daughters from the pool of associated hadrons
        //  with which the candidate is being correlated (will not have to do it for TPC-MFT case)
        if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig>) {
          if ((track1.prong0Id() == track2.globalIndex()) || (track1.prong1Id() == track2.globalIndex())) {
            continue;
          }
        }

        //  in case of MC-generated, do additional selection on MCparticles : charge and isPhysicalPrimary
        // if (processMc) {
        if constexpr (std::is_same_v<aodMcParticles, TTracksTrig> || std::is_same_v<aodMcParticles, TTracksAssoc>) {
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
          target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                      eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ,
                                      triggerWeight * associatedWeight);
        } else {
          target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                      eta1 - eta2, pt2, pt1, multiplicity, deltaPhi, posZ, invmass,
                                      triggerWeight * associatedWeight);
        }
      }
    }
  }

  // template <typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  // void mixCollisions(FilteredCollisionsWSelMult const& collisions, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  template <typename TCollisions, typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisions(TCollisions const& collisions, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  {
    // The first one that I call "Data" should work for data and mc rec
    using BinningTypeData = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::collision::PosZ, decltype(getPartsSize)>;

    BinningTypeData binningWithTracksSize{{getPartsSize}, {axisVertex, axisMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<TCollisions, TTracksTrig, TTracksAssoc, BinningTypeData> pair{binningWithTracksSize, nMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if constexpr (!std::is_same_v<aodMcCollisions, TCollisions>) {
        if (!(isCollisionSelected(collision1, false))) {
          continue;
        }
        if (!(isCollisionSelected(collision2, false))) {
          continue;
        }
      }

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);


      const auto multiplicityTracks2 = tracks2.size(); // get multiplicity of charged hadrons, which is used for slicing in mixing
      const auto multiplicityTracks1 = tracks1.size(); // get multiplicity of charged hadrons, which is used for slicing in mixing
      const auto vz = collision1.posZ();

      if constexpr (std::is_same_v<FilteredCollisionsWSelMultMC, TCollisions>) { // If MC
        registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
        fillTpcTpcChChMixedEventQAmc(multiplicityTracks2, vz, tracks1);
        /*
        if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig>) {
          registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
          fillHFMixingQA(multiplicity, vz, tracks1);
        } else {
          registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
          fillMixingQA(multiplicity, vz, tracks1);
        }
        */
      } else {                                                        // If not MC
        if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig>) { // DATA :  If TPC-TPC Hf-h case
          registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
          fillTpcTpcHfChMixedEventQA(multiplicityTracks2, vz, tracks1);
        } else if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) { // DATA : If TPC-MFT h-h case
          registry.fill(HIST("Data/TpcMft/HadronHadron/MixedEvent/hEventCountMixing"), bin);
          fillTpcMftChChMixedEventQA(multiplicityTracks1, vz, tracks1); // TPC tracks
          fillTpcMftChChMixedEventQA(multiplicityTracks2, vz, tracks2); // MFT tracks
        } else if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig> && std::is_same_v<aod::MFTTracks, TTracksAssoc>) { // DATA : If TPC-MFT HF-h case
          registry.fill(HIST("Data/TpcMft/HfHadron/MixedEvent/hEventCountMixing"), bin);
          fillTpcMftHfChMixedEventQA(multiplicityTracks1, vz, tracks1); // Candidates
          fillTpcMftHfChMixedEventQA(multiplicityTracks2, vz, tracks2); // MFT tracks
        } else {                                                             // DATA : If TPC-TPC h-h case
          registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
          fillTpcTpcChChMixedEventQA(multiplicityTracks2, vz, tracks1);
        }
      }

      corrContainer->fillEvent(multiplicityTracks2, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(corrContainer, tracks1, tracks2, multiplicityTracks2, collision1.posZ());
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

      // added this to try to compile when doing mixed event with aodMcParticles and aodMcCollisions (MC truth)
      // TODO : GET RID OF THE COLLISION SELECTION FOR MC TRUTH
      if constexpr (!std::is_same_v<aodMcCollisions, TCollisions>) {
        if (!(isCollisionSelected(collision1, false))) {
          continue;
        }
        if (!(isCollisionSelected(collision2, false))) {
          continue;
        }
      }

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicity = tracks2.size(); // get multiplicity of charged hadrons, which is used for slicing in mixing
      const auto vz = collision1.posZ();

      // TO BE DONE : ADD ONE MORE IF CONDITION TO FILL THE MC CASE
      // TODO : FILL NEW PLOTS FOR MCTRUTH ONLY
      registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
      fillTpcTpcChChMixedEventQAmc(multiplicity, vz, tracks1);
      /*
      if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig>) {
        registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
        fillHFMixingQA(multiplicity, vz, tracks1);
      } else {
        registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
        fillMixingQA(multiplicity, vz, tracks1);
      }
      */

      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
      // TO-DO : probably put the step of the correlation container in the template of this function
      fillCorrelations<CorrelationContainer::kCFStepAll>(corrContainer, tracks1, tracks2, multiplicity, collision1.posZ());
    }
  }

  // =====================================
  //    DATA : process same event correlations: TPC-TPC h-h case
  // =====================================

  void processSameTpcTpcChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    //  the event histograms below are only filled for h-h case
    //  because there is a possibility of double-filling if more correlation
    //  options are ran at the same time
    //  temporary solution, since other correlation options always have to be ran with h-h, too
    //  TODO: rewrite it in a more intelligent way
    const auto multiplicity = tracks.size();
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hMultiplicity"), multiplicity);
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hVtxZ"), collision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(collision.posZ(), multiplicity));
    registry.fill(HIST("Data/TpcTpc/HadronHadron/SameEvent/hEventCountSame"), bin);

    sameTPCTPCChCh->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    fillTpcTpcChChSameEventQA(multiplicity, tracks);
    // TO-DO : add if condition for when we will implant corrected correlations (kCFStepReconstructed -> kCFStepCorrected)
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTPCTPCChCh, tracks, tracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcChCh, "DATA : Process same-event correlations for TPC-TPC h-h case", true);

  // =====================================
  //    DATA : process same event correlations: TPC-TPC HF-h case
  // =====================================

  void processSameTpcTpcHfCh(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks,
                             HfCandidatesSel const& candidates)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }
    const auto multiplicity = tracks.size();

    sameTPCTPCHfCh->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    fillTpcTpcCandidateQA(candidates);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTPCTPCHfCh, candidates, tracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcHfCh, "DATA : Process same-event correlations for TPC-TPC HF-h case", true);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT h-h case
  // =====================================

  void processSameTpcMftChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks,
                             aod::MFTTracks const& mftTracks)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    const auto multiplicityTPC = tracks.size();
    const auto multiplicityMFT = mftTracks.size();

    sameTPCMFTChCh->fillEvent(multiplicityTPC, CorrelationContainer::kCFStepReconstructed);
    fillTpcMftChChSameEventQA(multiplicityTPC, tracks);
    fillTpcMftChChSameEventQA(multiplicityMFT, mftTracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTPCMFTChCh, tracks, mftTracks, multiplicityTPC, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChCh, "DATA : Process same-event correlations for TPC-MFT h-h case", true);

  // =====================================
  //    DATA : process same event correlations: TPC-MFT HF-h case
  // =====================================

  void processSameTpcMftHfCh(FilteredCollisionsWSelMult::iterator const& collision,
                             HfCandidatesSel const& candidates,
                             aod::MFTTracks const& mftTracks)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    const auto multiplicityCandidates = candidates.size();
    const auto multiplicityMFT = mftTracks.size();

    sameTPCMFTHfCh->fillEvent(multiplicityCandidates, CorrelationContainer::kCFStepReconstructed);
    fillTpcMftCandidateQA(candidates);
    fillTpcMftHfChSameEventQA(multiplicityMFT, mftTracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTPCMFTHfCh, candidates, mftTracks, multiplicityCandidates, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftHfCh, "DATA : Process same-event correlations for TPC-MFT HF-h case", true);

  // =====================================
  //    MONTE-CARLO : process same event correlations: TPC-TPC h-h case
  // =====================================

  void processSameTpcTpcChChmcREC(FilteredCollisionsWSelMultMC::iterator const& mcCollision,
                                  TracksWDcaSelMC const& mcTracks)
  {

    // NEED TO COMMENT THIS
    if (!(isCollisionSelected(mcCollision, true))) {
      return;
    }

    const auto multiplicity = mcTracks.size();
    registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hMultiplicity"), multiplicity);
    registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hVtxZ"), mcCollision.posZ());

    BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    int bin = baseBinning.getBin(std::make_tuple(mcCollision.posZ(), multiplicity));
    registry.fill(HIST("MC/Rec/TpcTpc/HadronHadron/SameEvent/hEventCountSame"), bin);

    sameTPCTPCChChMC->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    fillTpcTpcChChSameEventQAmc<CorrelationContainer::kCFStepReconstructed>(multiplicity, mcTracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(sameTPCTPCChChMC, mcTracks, mcTracks, multiplicity, mcCollision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcChChmcREC, "MONTE-CARLO : Process same-event correlations for TPC-TPC h-h case", true);

  // Katarina's version = MC Truth
  void processSameTpcTpcChChmcGEN(aodMcCollisions::iterator const& mcCollision,
                                  // FilteredCollisionsWSelMultMC const& collisions,
                                  aodMcParticles const& mcParticles)
  {

    // if (!(isCollisionSelected(mcCollision, true))) {
    //   return;
    // }

    // Not sure why to use this
    // if (collisions.size() == 0) {
    //  return;
    //}

    // if (!collision.has_mcCollision()) {
    //   LOGF(warning, "No MC collision for this collision, skip...");
    //   return;
    // }

    const auto multiplicity = mcParticles.size(); // Note: these are all MC particles after selection (not only primary)
    registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hMultiplicity"), multiplicity);
    registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hVtxZ"), mcCollision.posZ());

    //  fill correlations for all MC collisions
    // In Katka's code, the first time doing this does not fill the histograms, right now will be filled two times..
    auto multPrimaryCharge0 = fillTpcTpcChChSameEventQAmc<CorrelationContainer::kCFStepAll>(multiplicity, mcParticles);
    sameTPCTPCChChMC->fillEvent(multPrimaryCharge0, CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(sameTPCTPCChChMC, mcParticles, mcParticles, multPrimaryCharge0, mcCollision.posZ());

    // NOT USED BY KATARINA APPARENTLY
    // BinningPolicyBase<2> baseBinning{{axisVertex, axisMultiplicity}, true};
    // int bin = baseBinning.getBin(std::make_tuple(mcCollision.posZ(), multiplicity));
    // registry.fill(HIST("MC/Gen/TpcTpc/HadronHadron/SameEvent/hEventCountSame"), bin);

    //  fill correlations for MC collisions that have a reconstructed collision
    // got rid of the second const auto for multPrimaryCharge0
    multPrimaryCharge0 = fillTpcTpcChChSameEventQAmc<CorrelationContainer::kCFStepVertex>(multiplicity, mcParticles);
    sameTPCTPCChChMC->fillEvent(multPrimaryCharge0, CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(sameTPCTPCChChMC, mcParticles, mcParticles, multPrimaryCharge0, mcCollision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcChChmcGEN, "MONTE-CARLO : Process same-event correlations for TPC-TPC h-h case", true);

  // =====================================
  //    DATA : process mixed event correlations:TPC-TPC h-h case
  // =====================================
  // TO BECOME DATA & MC REC ?

  void processMixedTpcTpcChCh(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks, this](FilteredCollisionsWSelMult::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, tracks, tracks, getTracksSize, mixedTPCTPCChCh);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChCh, "DATA : Process mixed-event correlations for TPC-TPC h-h case", true);

  // =====================================
  //    DATA : process mixed event correlations: TPC-TPC HF-h case
  // =====================================

  void processMixedTpcTpcHfCh(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks,
                              HfCandidatesSel const& candidates)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks, this](FilteredCollisionsWSelMult::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache);
      auto size = associatedTracks.size();
      return size;
    };

    /*
    auto getTracksSize = [&candidates, this](FilteredCollisionsWSelMult::iterator const& col) {
      auto associatedTracks = candidates.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache);
      auto size = associatedTracks.size();
      return size;
    };
    */

    mixCollisions(collisions, candidates, tracks, getTracksSize, mixedTPCTPCHfCh);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcHfCh, "DATA : Process mixed-event correlations for TPC-TPC HF-h case", true);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT h-h case
  // =====================================

  void processMixedTpcMftChCh(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks,
                              aod::MFTTracks const& mftTracks)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks, this](FilteredCollisionsWSelMult::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache);
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, tracks, mftTracks, getTracksSize, mixedTPCMFTChCh);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftChCh, "DATA : Process mixed-event correlations for TPC-MFT h-h case", true);

  // =====================================
  //    DATA : process mixed event correlations: TPC-MFT HF-h case
  // =====================================

  void processMixedTpcMftHfCh(FilteredCollisionsWSelMult const& collisions,
                              HfCandidatesSel const& candidates,
                              aod::MFTTracks const& mftTracks)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&candidates, this](FilteredCollisionsWSelMult::iterator const& col) {
      // Still o2::aod::track::collisionId with HF ???
      auto associatedTracks = candidates.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache);
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, candidates, mftTracks, getTracksSize, mixedTPCMFTHfCh);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftHfCh, "DATA : Process mixed-event correlations for TPC-MFT HF-h case", true);

  // =====================================
  //    MONTE-CARLO : process mixed event correlations: TPC-TPC h-h case
  // =====================================

  // MC rec
  void processMixedTpcTpcChChmcREC(FilteredCollisionsWSelMultMC const& mcCollisions,
                                   TracksWDcaSelMC const& mcTracks)
  {
    // use normal index instead of globalIndex for MixedEvent ??

    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&mcTracks, this](FilteredCollisionsWSelMultMC::iterator const& mcCol) {
      auto associatedTracks = mcTracks.sliceByCached(o2::aod::track::collisionId, mcCol.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(mcCollisions, mcTracks, mcTracks, getTracksSize, mixedTPCTPCChChMC);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChChmcREC, "MONTE-CARLO : Process mixed-event correlations for TPC-TPC h-h case", true);

  // MC gen
  void processMixedTpcTpcChChmcGEN(aodMcCollisions const& mcCollisions,
                                   aodMcParticles const& mcParticles,
                                   FilteredCollisionsWSelMultMC const& collisions)
  {
    // use normal index instead of globalIndex for MixedEvent ??

    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&mcParticles, this](aodMcCollisions::iterator const& mcCol) {
      auto associatedTracks = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId, mcCol.globalIndex(), this->cache); // it's cached, so slicing/grouping happens only once
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisionsMcTruth(mcCollisions, mcParticles, mcParticles, getTracksSize, mixedTPCTPCChChMC);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChChmcGEN, "MONTE-CARLO : Process mixed-event correlations for TPC-TPC h-h case", true);
}; // End of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskFlow>(cfgc)};
}
