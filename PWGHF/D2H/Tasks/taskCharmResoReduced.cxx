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

/// \file taskCharmResoReduced.cxx
/// \brief Charmed Resonances analysis task
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, University and INFN Torino

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

// #include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum DecayChannel : uint8_t {
  Ds1ToDstarK0s = 0,
  Ds2StarToDplusK0s,
  XcToDplusLambda,
  LambdaDminus
};

struct HfTaskCharmResoReduced {
  Configurable<float> ptMinReso{"ptMinReso", 5, "Discard events with smaller pT"};
  Configurable<bool> cutBeforeMixing{"cutBeforeMixing", false, "Apply pT cut to candidates before event mixing"};
  Configurable<bool> cutAfterMixing{"cutAfterMixing", false, "Apply pT cut to candidates after event mixing"};
  // Configurables axis for histos
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisPtProng0{"axisPtProng0", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "prong0 bach. #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisPtProng1{"axisPtProng1", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "prong1 bach. #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisInvMassReso{"axisInvMassReso", {200, 2.34, 2.74}, "inv. mass (DV_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng0{"axisInvMassProng0", {175, 1.70, 2.05}, "inv. mass (D) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng1{"axisInvMassProng1", {80, 0.46, 0.54}, "inv. mass ({V}_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisCosThetaStar{"axisCosThetaStar", {40, -1, 1}, "cos(#vartheta*)"};
  ConfigurableAxis axisBkgBdtScore{"axisBkgBdtScore", {100, 0, 1}, "bkg BDT Score"};
  ConfigurableAxis axisNonPromptBdtScore{"axisNonPromptBdtScore", {100, 0, 1}, "non-prompt BDT Score"};
  // Configurables for ME
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<int> numberEventsToSkip{"numberEventsToSkip", -1, "Number of events to Skip in ME process"};
  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 45., 60., 75., 95, 250}, "event multiplicity pools (PV contributors for now)"};
  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -4, -1, 1, 4, 10.0}, "z vertex position pools"};
  // ConfigurableAxis bzPoolBins{"bzPoolBins", {2, -10, 10}, "Bz of collision"};

  using ReducedResoWithMl = soa::Join<aod::HfCandCharmReso, aod::HfCharmResoMLs>;
  SliceCache cache;
  Preslice<ReducedResoWithMl> resoPerCollision = aod::hf_track_index_reduced::hfRedCollisionId;

  // Histogram Registry
  HistogramRegistry registry;

  // init
  void init(InitContext&)
  {
    registry.add("hMass", "Charm resonance candidates inv. mass", {HistType::kTH1F, {axisInvMassReso}});
    registry.add("hMassProng0", "D daughters inv. mass", {HistType::kTH1F, {axisInvMassProng0}});
    registry.add("hMassProng1", "V0 daughter inv. mass", {HistType::kTH1F, {axisInvMassProng1}});
    registry.add("hPt", "Charm resonance candidates pT", {HistType::kTH1F, {axisPt}});
    registry.add("hPtProng0", "D daughters pT", {HistType::kTH1F, {axisPtProng0}});
    registry.add("hPtProng1", "V0 daughter pT", {HistType::kTH1F, {axisPtProng1}});
    registry.add("hNPvCont", "Collision number of PV contributors ; N contrib ; entries", {HistType::kTH1F, {{100, 0, 250}}});
    registry.add("hZvert", "Collision Z Vtx ; z PV [cm] ; entries", {HistType::kTH1F, {{120, -12., 12.}}});
    registry.add("hBz", "Collision Bz ; Bz [T] ; entries", {HistType::kTH1F, {{20, -10., 10.}}});
    registry.add("hSparse", "THn for production studies with cosThStar and BDT scores", HistType::kTHnSparseF, {axisPt, axisPtProng0, axisPtProng1, axisInvMassReso, axisInvMassProng0, axisInvMassProng1, axisCosThetaStar, axisBkgBdtScore, axisNonPromptBdtScore});

    if (doprocessDs1DataMixedEvent || doprocessDs2StarDataMixedEvent) {
      registry.add("hNPvContCorr", "Collision number of PV contributors ; N contrib ; N contrib", {HistType::kTH2F, {{100, 0, 250}, {100, 0, 250}}});
      registry.add("hZvertCorr", "Collision Z Vtx ; z PV [cm] ; z PV [cm]", {HistType::kTH2F, {{120, -12., 12.}, {120, -12., 12.}}});
      registry.add("hMassProng0Corr", "D daughters inv. mass", {HistType::kTH2F, {axisInvMassProng0, axisInvMassProng0}});
      registry.add("hMassProng1Corr", "V0 daughter inv. mass", {HistType::kTH2F, {axisInvMassProng1, axisInvMassProng1}});
    }
  }

  // Fill histograms
  /// \tparam channel is the decay channel of the Resonance
  /// \param candidate is a candidate
  /// \param coll is a reduced collision
  template <DecayChannel channel, typename Cand, typename Coll>
  void fillHisto(const Cand& candidate, const Coll& collision)
  {
    // Collision properties
    registry.fill(HIST("hNPvCont"), collision.numContrib());
    registry.fill(HIST("hZvert"), collision.posZ());
    registry.fill(HIST("hBz"), collision.bz());
    // Candidate properties
    registry.fill(HIST("hMass"), candidate.invMass());
    registry.fill(HIST("hMassProng0"), candidate.invMassProng0());
    registry.fill(HIST("hMassProng1"), candidate.invMassProng1());
    registry.fill(HIST("hPt"), candidate.pt());
    registry.fill(HIST("hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1"), candidate.ptProng1());
    float cosThetaStar{0.};
    switch (channel) {
      case DecayChannel::Ds1ToDstarK0s:
        cosThetaStar = candidate.cosThetaStarDs1();
        break;
      case DecayChannel::Ds2StarToDplusK0s:
        cosThetaStar = candidate.cosThetaStarDs2Star();
        break;
      default:
        cosThetaStar = candidate.cosThetaStarXiC3055();
        break;
    }
    registry.fill(HIST("hSparse"), candidate.pt(), candidate.ptProng0(), candidate.ptProng1(), candidate.invMass(), candidate.invMassProng0(), candidate.invMassProng1(), cosThetaStar, candidate.mlScoreBkgProng0(), candidate.mlScoreNonpromptProng0());
  } // fillHisto

  // Process data
  /// \tparam channel is the decay channel of the Resonance
  /// \param Coll is the reduced collisions table
  /// \param Cand is the candidates table
  template <DecayChannel channel, typename Coll, typename Candidates>
  void processData(Coll const&, Candidates const& candidates)
  {
    for (const auto& cand : candidates) {
      if (cand.pt() < ptMinReso) {
        continue;
      }
      auto coll = cand.template hfRedCollision_as<Coll>();
      fillHisto<channel>(cand, coll);
    }
  }

  // Process data with Mixed Event
  /// \tparam channel is the decay channel of the Resonance
  /// \param Coll is the reduced collisions table
  /// \param Cand is the candidates table
  template <DecayChannel channel, typename Coll, typename Candidates>
  void processDataMixedEvent(Coll const& collisions, Candidates const& candidates)
  {
    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    auto candsTuple = std::make_tuple(candidates);
    SameKindPair<aod::HfRedCollisions, Candidates, BinningType> pairs{corrBinning, numberEventsMixed, numberEventsToSkip, collisions, candsTuple, &cache};
    for (const auto& [collision1, cands1, collision2, cands2] : pairs) {
      // For each couple of candidate resonances I can make 2 mixed candidates by swithching daughters
      for (const auto& [cand1, cand2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(cands1, cands2))) {
        if (cutBeforeMixing && (cand1.pt() < ptMinReso || cand2.pt() < ptMinReso)) {
          continue;
        }
        float ptME1 = RecoDecay::pt(cand1.pVectorProng0(), cand2.pVectorProng1());
        float invMassME1;
        float cosThetaStarME1;
        if (!cutAfterMixing || ptME1 > ptMinReso) {
          switch (channel) {
            case DecayChannel::Ds1ToDstarK0s:
              invMassME1 = RecoDecay::m(std::array{cand1.pVectorProng0(), cand2.pVectorProng1()}, std::array{o2::constants::physics::MassDStar, o2::constants::physics::MassK0Short});
              cosThetaStarME1 = RecoDecay::cosThetaStar(std::array{cand1.pVectorProng0(), cand2.pVectorProng1()}, std::array{o2::constants::physics::MassDStar, o2::constants::physics::MassK0Short}, invMassME1, 1);
              break;
            case DecayChannel::Ds2StarToDplusK0s:
              invMassME1 = RecoDecay::m(std::array{cand1.pVectorProng0(), cand2.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassK0Short});
              cosThetaStarME1 = RecoDecay::cosThetaStar(std::array{cand1.pVectorProng0(), cand2.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassK0Short}, invMassME1, 1);
              break;
            default:
              invMassME1 = RecoDecay::m(std::array{cand1.pVectorProng0(), cand2.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassLambda0});
              cosThetaStarME1 = RecoDecay::cosThetaStar(std::array{cand1.pVectorProng0(), cand2.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassLambda0}, invMassME1, 1);
              break;
          }
          registry.fill(HIST("hMass"), invMassME1);
          registry.fill(HIST("hPt"), ptME1);
          registry.fill(HIST("hNPvContCorr"), collision1.numContrib(), collision2.numContrib());
          registry.fill(HIST("hZvertCorr"), collision1.posZ(), collision2.posZ());
          registry.fill(HIST("hMassProng0Corr"), cand1.invMassProng0(), cand2.invMassProng0());
          registry.fill(HIST("hMassProng1Corr"), cand1.invMassProng1(), cand2.invMassProng1());
          registry.fill(HIST("hSparse"), ptME1, cand1.ptProng0(), cand2.ptProng1(), invMassME1, cand1.invMassProng0(), cand2.invMassProng1(), cosThetaStarME1, cand1.mlScoreBkgProng0(), cand1.mlScoreNonpromptProng0());
        }
        float ptME2 = RecoDecay::pt(cand2.pVectorProng0(), cand1.pVectorProng1());
        float invMassME2;
        float cosThetaStarME2;
        if (!cutAfterMixing || ptME2 > ptMinReso) {
          switch (channel) {
            case DecayChannel::Ds1ToDstarK0s:
              invMassME2 = RecoDecay::m(std::array{cand2.pVectorProng0(), cand1.pVectorProng1()}, std::array{o2::constants::physics::MassDStar, o2::constants::physics::MassK0Short});
              cosThetaStarME2 = RecoDecay::cosThetaStar(std::array{cand2.pVectorProng0(), cand1.pVectorProng1()}, std::array{o2::constants::physics::MassDStar, o2::constants::physics::MassK0Short}, invMassME2, 1);
              break;
            case DecayChannel::Ds2StarToDplusK0s:
              invMassME2 = RecoDecay::m(std::array{cand2.pVectorProng0(), cand1.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassK0Short});
              cosThetaStarME2 = RecoDecay::cosThetaStar(std::array{cand2.pVectorProng0(), cand1.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassK0Short}, invMassME2, 1);
              break;
            default:
              invMassME2 = RecoDecay::m(std::array{cand2.pVectorProng0(), cand1.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassLambda0});
              cosThetaStarME2 = RecoDecay::cosThetaStar(std::array{cand2.pVectorProng0(), cand1.pVectorProng1()}, std::array{o2::constants::physics::MassDPlus, o2::constants::physics::MassLambda0}, invMassME2, 1);
              break;
          }
          registry.fill(HIST("hMass"), invMassME2);
          registry.fill(HIST("hPt"), ptME2);
          registry.fill(HIST("hSparse"), ptME2, cand2.ptProng0(), cand1.ptProng1(), invMassME2, cand2.invMassProng0(), cand1.invMassProng1(), cosThetaStarME2, cand2.mlScoreBkgProng0(), cand2.mlScoreNonpromptProng0());
        }
      }
    }
  }

  // process functions

  void processDs1Data(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates)
  {
    processData<DecayChannel::Ds1ToDstarK0s>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs1Data, "Process data", true);

  void processDs1DataMixedEvent(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates)
  {
    processDataMixedEvent<DecayChannel::Ds1ToDstarK0s>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs1DataMixedEvent, "Process data with Event Mixing", false);

  void processDs2StarData(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates)
  {
    processData<DecayChannel::Ds2StarToDplusK0s>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs2StarData, "Process data", false);

  void processDs2StarDataMixedEvent(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates)
  {
    processDataMixedEvent<DecayChannel::Ds2StarToDplusK0s>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs2StarDataMixedEvent, "Process data with Event Mixing", false);

}; // struct HfTaskCharmResoReduced
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmResoReduced>(cfgc)};
}
