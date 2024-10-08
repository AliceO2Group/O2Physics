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

#include <TDirectory.h>
#include <TH1F.h>
#include <THn.h>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
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

  using FilteredCollisionsWSelMult = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using TracksWDcaSel = soa::Filtered<soa::Join<aod::TracksWDca, aod::TrackSelection>>;
  using HfCandidatesSel = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;

  //  Collision filters
  //  FIXME: The filter is applied also on the candidates! Beware!
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < zVertexMax;
  //  Charged track filters
  Filter trackFilter = (nabs(aod::track::eta) < etaTrackAssocMax) &&
                       (aod::track::pt > ptTrackAssocMin) &&
                       requireGlobalTrackWoPtEtaInFilter();
  //  HF candidate filter
  //  TODO: use Partition instead of filter
  Filter candidateFilter = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;

  Preslice<aod::Tracks> perCol = aod::track::collisionId;

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

  OutputObj<CorrelationContainer> sameTPCTPCChCh{"sameTPCTPCChCh"};
  OutputObj<CorrelationContainer> mixedTPCTPCChCh{"mixedTPCTPCChCh"};
  OutputObj<CorrelationContainer> sameTPCTPCHfCh{"sameTPCTPCHfCh"};
  OutputObj<CorrelationContainer> mixedTPCTPCHfCh{"mixedTPCTPCHfCh"};
  OutputObj<CorrelationContainer> sameTPCMFTChCh{"sameTPCMFTChCh"};
  OutputObj<CorrelationContainer> mixedTPCMFTChCh{"mixedTPCMFTChCh"};

  //  =========================
  //      init()
  //  =========================
  void init(InitContext&)
  {
    //  EVENT HISTOGRAMS
    constexpr int kNBinsEvents = 3;
    registry.add("Data/hEventCounter", "hEventCounter", {HistType::kTH1F, {{kNBinsEvents, 0.5, 0.5 + kNBinsEvents}}});
    //  set axes of the event counter histogram
    std::string labels[kNBinsEvents];
    labels[0] = "all";
    labels[1] = "after trigger selection (Run 2)";
    labels[2] = "after Physics selection";
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("Data/hEventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hVtxZ", "hVtxZ", {HistType::kTH1F, {{400, -50, 50}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hNtracks", "hNtracks", {HistType::kTH1F, {{500, 0, 500}}});

    //  histograms for event mixing
    const int maxMixBin = axisMultiplicity->size() * 14; // 14 bins for z-vertex
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEventCountSame", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hMultiplicityMixing", "hMultiplicityMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hVtxZMixing", "hVtxZMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hNtracksMixing", "hNtracksMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hMultiplicityHFMixing", "hMultiplicityHFMixing", {HistType::kTH1F, {{500, 0, 500}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hVtxZHFMixing", "hVtxZHFMixing", {HistType::kTH1F, {{100, -10, 10}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hNtracksHFMixing", "hNtracksHFMixing", {HistType::kTH1F, {{500, 0, 500}}});

    //  TRACK HISTOGRAMS
    //  histograms for associated particles
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hYields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEtaPhi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPt", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hEta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/SameEvent/hPhi", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

    //  histograms for particles in event mixing
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hPtMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hEtaMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HadronHadron/MixedEvent/hPhiMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

    //  histograms for MFT tracks
    registry.add("Data/TpcMft/HadronHadron/hEtaPhiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, TwoPI, "#varphi"}}});
    registry.add("Data/TpcMft/HadronHadron/hEtaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcMft/HadronHadron/hPhiMFT", "phiMFT", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

    //  histograms for candidates
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

    //  histograms for candidates in event mixing
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hPtHFMixing", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hEtaHFMixing", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("Data/TpcTpc/HfHadron/MixedEvent/hPhiHFMixing", "phi", {HistType::kTH1F, {{100, 0, TwoPI, "#varphi"}}});

    //  set axes of the correlation container
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

    sameTPCTPCChCh.setObject(new CorrelationContainer("sameTPCTPCChCh", "sameTPCTPCChCh", corrAxis, effAxis, {}));
    mixedTPCTPCChCh.setObject(new CorrelationContainer("mixedTPCTPCChCh", "mixedTPCTPCChCh", corrAxis, effAxis, {}));
    sameTPCTPCHfCh.setObject(new CorrelationContainer("sameTPCTPCHfCh", "sameTPCTPCHfCh", corrAxis, effAxis, userAxis));
    mixedTPCTPCHfCh.setObject(new CorrelationContainer("mixedTPCTPCHfCh", "mixedTPCTPCHfCh", corrAxis, effAxis, userAxis));
    sameTPCMFTChCh.setObject(new CorrelationContainer("sameTPCMFTChCh", "sameTPCMFTChCh", corrAxis, effAxis, {}));
    mixedTPCMFTChCh.setObject(new CorrelationContainer("mixedTPCMFTChCh", "mixedTPCMFTChCh", corrAxis, effAxis, {}));
  }

  //  ---------------
  //    templates
  //  FIXME: Some collisions are rejected here, what causes (part of) differences with the D0 task
  //  ---------------
  template <typename TCollision>
  bool isCollisionSelected(TCollision const& collision, bool fillHistograms = false)
  {
    if (processRun2 == true) {
      //  Run 2: trigger selection for data case
      if (fillHistograms)
        registry.fill(HIST("Data/hEventCounter"), 1);
      if (!processMc) {
        if (!collision.alias_bit(kINT7)) {
          return false;
        }
      }
      //  Run 2: further offline selection
      if (fillHistograms)
        registry.fill(HIST("Data/hEventCounter"), 2);
      if (!collision.sel7()) {
        return false;
      }
      if (fillHistograms)
        registry.fill(HIST("Data/hEventCounter"), 3);
    } else {
      //  Run 3: selection
      if (fillHistograms)
        registry.fill(HIST("Data/hEventCounter"), 1);
      if (!collision.sel8()) {
        return false;
      }
      if (fillHistograms)
        registry.fill(HIST("Data/hEventCounter"), 3);
    }
    return true;
  }

  template <typename TTracks>
  void fillQA(float multiplicity, TTracks const& tracks)
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

  template <typename TTracks>
  void fillMixingQA(float multiplicity, float vz, TTracks const& tracks)
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

  template <typename TTracks>
  void fillHFMixingQA(float multiplicity, float vz, TTracks const& tracks)
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

  template <typename TTracks>
  void fillMFTQA(float multiplicity, TTracks const& tracks)
  {
    for (const auto& track1 : tracks) {
      registry.fill(HIST("Data/TpcMft/HadronHadron/hEtaMFT"), track1.eta());
      float phi = track1.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("Data/TpcMft/HadronHadron/hPhiMFT"), phi);
      registry.fill(HIST("Data/TpcMft/HadronHadron/hEtaPhiMFT"), multiplicity, track1.eta(), phi);
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

  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTracks>
  void fillCandidateQA(TTracks const& candidates)
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

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc>
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
        //  with which the candidate is being correlated
        if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig>) {
          if ((track1.prong0Id() == track2.globalIndex()) || (track1.prong1Id() == track2.globalIndex())) {
            continue;
          }
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

  template <typename TTracksTrig, typename TTracksAssoc, typename TLambda>
  void mixCollisions(FilteredCollisionsWSelMult const& collisions, TTracksTrig const& tracks1, TTracksAssoc const& tracks2, TLambda getPartsSize, OutputObj<CorrelationContainer>& corrContainer)
  {
    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getPartsSize)>, aod::collision::PosZ, decltype(getPartsSize)>;
    BinningType binningWithTracksSize{{getPartsSize}, {axisVertex, axisMultiplicity}, true};
    auto tracksTuple = std::make_tuple(tracks1, tracks2);
    Pair<FilteredCollisionsWSelMult, TTracksTrig, TTracksAssoc, BinningType> pair{binningWithTracksSize, nMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pair) {

      if (!(isCollisionSelected(collision1, false))) {
        continue;
      }
      if (!(isCollisionSelected(collision2, false))) {
        continue;
      }

      auto binningValues = binningWithTracksSize.getBinningValues(collision1, collisions);
      int bin = binningWithTracksSize.getBin(binningValues);

      const auto multiplicity = tracks2.size(); // get multiplicity of charged hadrons, which is used for slicing in mixing
      const auto vz = collision1.posZ();

      if constexpr (std::is_same_v<HfCandidatesSel, TTracksTrig>) {
        registry.fill(HIST("Data/TpcTpc/HfHadron/MixedEvent/hEventCountHFMixing"), bin);
        fillHFMixingQA(multiplicity, vz, tracks1);
      } else {
        registry.fill(HIST("Data/TpcTpc/HadronHadron/MixedEvent/hEventCountMixing"), bin);
        fillMixingQA(multiplicity, vz, tracks1);
      }

      corrContainer->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations(corrContainer, tracks1, tracks2, multiplicity, collision1.posZ());
    }
  }

  // =====================================
  //    process same event correlations: h-h case
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

    fillQA(multiplicity, tracks);
    fillCorrelations(sameTPCTPCChCh, tracks, tracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcChCh, "Process same-event correlations for TPC-TPC h-h case", true);

  // =====================================
  //    process same event correlations: HF-h case
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

    fillCandidateQA(candidates);
    fillCorrelations(sameTPCTPCHfCh, candidates, tracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcTpcHfCh, "Process same-event correlations for TPC-TPC HF-h case", true);

  // =====================================
  //    process same event correlations: h-MFT case
  // =====================================
  void processSameTpcMftChCh(FilteredCollisionsWSelMult::iterator const& collision,
                             TracksWDcaSel const& tracks,
                             aod::MFTTracks const& mfttracks)
  {
    if (!(isCollisionSelected(collision, true))) {
      return;
    }

    const auto multiplicity = tracks.size();

    sameTPCMFTChCh->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillMFTQA(multiplicity, mfttracks);
    fillCorrelations(sameTPCMFTChCh, tracks, mfttracks, multiplicity, collision.posZ());
  }
  PROCESS_SWITCH(HfTaskFlow, processSameTpcMftChCh, "Process same-event correlations for TPC-MFT h-h case", true);

  // =====================================
  //    process mixed event correlations: h-h case
  // =====================================
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
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcChCh, "Process mixed-event correlations for TPC-TPC h-h case", true);

  // =====================================
  //    process mixed event correlations: h-HF case
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

    mixCollisions(collisions, candidates, tracks, getTracksSize, mixedTPCTPCHfCh);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcTpcHfCh, "Process mixed-event correlations for TPC-TPC HF-h case", true);

  // =====================================
  //    process mixed event correlations: h-MFT case
  // =====================================
  void processMixedTpcMftChCh(FilteredCollisionsWSelMult const& collisions,
                              TracksWDcaSel const& tracks,
                              aod::MFTTracks const& mfttracks)
  {
    //  we want to group collisions based on charged-track multiplicity
    auto getTracksSize = [&tracks, this](FilteredCollisionsWSelMult::iterator const& col) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, col.globalIndex(), this->cache);
      auto size = associatedTracks.size();
      return size;
    };

    mixCollisions(collisions, tracks, mfttracks, getTracksSize, mixedTPCMFTChCh);
  }
  PROCESS_SWITCH(HfTaskFlow, processMixedTpcMftChCh, "Process mixed-event correlations for TPC-MFT h-h case", true);
}; // End of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskFlow>(cfgc)};
}
