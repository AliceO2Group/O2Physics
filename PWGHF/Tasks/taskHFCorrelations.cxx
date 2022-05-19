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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/GroupSlicer.h"
#include "Framework/StepTHn.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "DataFormatsParameters/GRPObject.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>
#include <THn.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::analysis::hf_cuts_d0_topik;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct TaskHfCorrelations {

  HistogramRegistry registry{"registry"};
  OutputObj<CorrelationContainer> sameTPCTPCCh{"sameEventTPCTPCChHadrons"};
  OutputObj<CorrelationContainer> sameTPCMFTCh{"sameEventTPCMFTChHadrons"};
  OutputObj<CorrelationContainer> sameHF{"sameEventHFHadrons"};
  OutputObj<CorrelationContainer> mixedTPCTPCCh{"mixedEventTPCTPCChHadrons"};

  //  configurables for processing options
  Configurable<bool> processTPCTPChh{"processTPCTPChh", "true", "Flag to process TPC-TPC h-h correlations"};
  Configurable<bool> processTPCMFThh{"processTPCMFThh", "true", "Flag to process TPC-MFT h-h correlations"};
  Configurable<bool> processHFHadrons{"processHFHadrons", "false", "Flag to process HF-h correlations"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  //  configurables for collisions
  Configurable<float> cfgCutVertex{"cfgCutVertex", 7.0f, "Accepted z-vertex range"};

  //  configurables for associated particles
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.5f, "Minimum pT for tracks"};

  //  configurables for HF candidates
  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_d0_topik::pTBins_v}, "pT bin limits"};

  //  configurables for containers
  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity axis for histograms"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  //  Collision filters
  Filter collisionVtxZFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;

  using BinningType = BinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M>;

  //  Charged track filters
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) &&
                       (aod::track::pt > cfgCutPt) &&
                       ((aod::track::isGlobalTrack == (uint8_t) true) ||
                        (aod::track::isGlobalTrackSDD == (uint8_t) true));
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended, aod::TrackSelection>>;

  //  HF candidate filter
  Filter candidateFilter = (aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 ||
                            aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar);
  using hfCandidates = soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>>;

  //  =========================
  //      init()
  //  =========================
  void init(o2::framework::InitContext&)
  {
    registry.add("eventCounter", "eventCounter", {HistType::kTH1F, {{2, 0.5, 2.5}}});
    //  set axes of the event counter histogram
    const int nBins = 2;
    std::string labels[nBins];
    labels[0] = "all";
    labels[1] = "after sel8";
    for (int iBin = 0; iBin < nBins; iBin++) {
      registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    registry.add("hMultiplicity", "hMultiplicity", {HistType::kTH1F, {{200, 0, 200}}});

    //  histograms for event mixing
    const int maxMixBin = axisMultiplicity->size() * axisVertex->size();
    registry.add("eventcount", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});

    //  histograms for associated particles
    registry.add("yields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * M_PI, "#varphi"}}});
    registry.add("pT", "pT", {HistType::kTH1F, {{100, 0, 10, "p_{T}"}}});
    registry.add("eta", "eta", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("phi", "phi", {HistType::kTH1F, {{100, 0, 2 * M_PI, "#varphi"}}});

    //  histograms for MFT tracks
    registry.add("etaphiMFT", "multiplicity vs eta vs phi in MFT", {HistType::kTH3F, {{200, 0, 200, "multiplicity"}, {100, -2, 2, "#eta"}, {200, 0, 2 * M_PI, "#varphi"}}});
    registry.add("etaMFT", "etaMFT", {HistType::kTH1F, {{100, -4, 4, "#eta"}}});
    registry.add("phiMFT", "phiMFT", {HistType::kTH1F, {{100, -1 * M_PI, M_PI, "#varphi"}}});

    //  histograms for candidates
    auto vbins = (std::vector<double>)bins;

    registry.add("hptcand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hptprong0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hptprong1", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0, 10.}}});
    registry.add("hmass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hselectionstatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    //  set axes of the correlation container
    std::vector<AxisSpec> axisList = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisEtaEfficiency, "#eta"},
                                      {axisPtEfficiency, "p_{T} (GeV/c)"},
                                      {axisVertexEfficiency, "z-vtx (cm)"}};
    sameTPCTPCCh.setObject(new CorrelationContainer("sameEventTPCTPCChHadrons", "sameEventTPCTPCChHadrons", axisList));
    sameTPCMFTCh.setObject(new CorrelationContainer("sameEventTPCMFTChHadrons", "sameEventTPCMFTChHadrons", axisList));
    sameHF.setObject(new CorrelationContainer("sameEventHFHadrons", "sameEventHFHadrons", axisList));
    mixedTPCTPCCh.setObject(new CorrelationContainer("mixedEventTPCTPCChHadrons", "mixedEventTPCTPCChHadrons", axisList));
  }

  //  ---------------
  //    templates
  //  ---------------
  template <typename TTracks>
  void fillQA(float multiplicity, TTracks tracks)
  {
    for (auto& track1 : tracks) {
      //  Temporary track selection here, until mixing is solved for filtered tables
      if (isAcceptedTrack(track1) == false) {
        continue;
      }
      registry.fill(HIST("pT"), track1.pt());
      registry.fill(HIST("eta"), track1.eta());
      registry.fill(HIST("phi"), track1.phi());
      registry.fill(HIST("yields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphi"), multiplicity, track1.eta(), track1.phi());
    }
  }

  template <typename TTracks>
  void fillMFTQA(float multiplicity, TTracks tracks)
  {
    for (auto& track1 : tracks) {
      registry.fill(HIST("etaMFT"), track1.eta());
      registry.fill(HIST("phiMFT"), track1.phi());
      registry.fill(HIST("etaphiMFT"), multiplicity, track1.eta(), track1.phi());
    }
  }

  template <typename TCollision>
  bool isCollisionSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return false;
    }
    return true;
  }

  //  TODO: Check how to put this into a Filter
  template <typename TTrack>
  bool isAcceptedCandidate(TTrack candidate)
  {
    if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
      return false;
    }
    if (cutYCandMax >= 0. && std::abs(YD0(candidate)) > cutYCandMax) {
      return false;
    }
    return true;
  }

  //  TODO: Note: we do not need all these plots since they are in D0 and Lc task -> remove it after we are sure this works
  template <typename TTracks>
  void fillCandidateQA(float multiplicity, TTracks candidates)
  {
    for (auto& candidate : candidates) {
      if (isAcceptedCandidate(candidate) == false) {
        continue;
      }
      registry.fill(HIST("hptcand"), candidate.pt());
      registry.fill(HIST("hptprong0"), candidate.ptProng0());
      registry.fill(HIST("hptprong1"), candidate.ptProng1());
      registry.fill(HIST("hmass"), InvMassD0(candidate), candidate.pt());
      registry.fill(HIST("hdeclength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hdeclengthxy"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hCTS"), CosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("hCt"), CtD0(candidate), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hselectionstatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    }
  }

  template <typename TTarget, typename TTracksTrig, typename TTracksAssoc>
  void fillCorrelations(TTarget target, TTracksTrig tracks1, TTracksAssoc tracks2, float multiplicity, float posZ)
  {
    auto triggerWeight = 1;
    auto associatedWeight = 1;

    for (auto& track1 : tracks1) {
      //  TODO: add getter for NUE trigger efficiency here

      //  TODO: Check how to put this into a Filter
      if constexpr (std::is_same_v<hfCandidates, TTracksTrig>) {
        if (isAcceptedCandidate(track1) == false) {
          continue;
        }
      }

      target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), multiplicity, posZ, triggerWeight);

      for (auto& track2 : tracks2) {
        if constexpr (std::is_same_v<TTracksAssoc, TTracksTrig>) { // case of h-h correlations where the two types of tracks are the same
          if (track1 == track2) {
            continue;
          }
        }

        //  TODO: add getter for NUE associated efficiency here

        //  TODO: add pair cuts on phi*

        float deltaPhi = track1.phi() - track2.phi();
        if (deltaPhi > 1.5f * PI) {
          deltaPhi -= TwoPI;
        }
        if (deltaPhi < -PIHalf) {
          deltaPhi += TwoPI;
        }

        target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                    track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ,
                                    triggerWeight * associatedWeight);
      }
    }
  }

  // =====================================
  //    process same event correlations
  // =====================================
  void processSame(aodCollisions::iterator const& collision,
                   aodTracks const& tracks,
                   aod::MFTTracks const& mfttracks,
                   hfCandidates const& candidates)
  {
    const auto multiplicity = collision.multFV0M(); //  multV0M ? (work on adding centrality selection for Run3 ongoing)
    registry.fill(HIST("hMultiplicity"), multiplicity);

    registry.fill(HIST("eventCounter"), 1);
    if (isCollisionSelected(collision) == false) {
      return;
    }
    registry.fill(HIST("eventCounter"), 2);

    if (processTPCTPChh == true) {
      fillQA(multiplicity, tracks);
      fillCorrelations(sameTPCTPCCh, tracks, tracks, multiplicity, collision.posZ());
    }

    if (processTPCMFThh == true) {
      fillMFTQA(multiplicity, mfttracks);
      fillCorrelations(sameTPCMFTCh, tracks, mfttracks, multiplicity, collision.posZ());
    }

    if (processHFHadrons == true) {
      fillCandidateQA(multiplicity, candidates);
      fillCorrelations(sameHF, candidates, tracks, multiplicity, collision.posZ());
    }
  }
  PROCESS_SWITCH(TaskHfCorrelations, processSame, "Process same event", true);

  // =====================================
  //    process mixed event correlations
  // =====================================
  //  TODO: these collisions do not contain trigger selection, because SameKindPair needs table, not iterator -> check if it can be solved
  void processMixed(aodCollisions& collisions,
                    aodTracks const& tracks)
  {
    BinningType pairBinning{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default)
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aodCollisions, aodTracks, BinningType> pair{pairBinning, cfgNoMixedEvents, -1, collisions, tracksTuple};

    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (isCollisionSelected(collision1) == false) {
        continue;
      }
      if (isCollisionSelected(collision2) == false) {
        continue;
      }

      const auto multiplicity = collision1.multFV0M();

      int bin = pairBinning.getBin({collision1.posZ(), multiplicity});
      registry.fill(HIST("eventcount"), bin);

      if (processTPCTPChh == true) {
        fillCorrelations(mixedTPCTPCCh, tracks1, tracks2, multiplicity, collision1.posZ());
      }
    }
  }
  PROCESS_SWITCH(TaskHfCorrelations, processMixed, "Process mixed event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TaskHfCorrelations>(cfgc)};
}
