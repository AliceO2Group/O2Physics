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

  // process functions

  void processDs1Data(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates)
  {
    processData<DecayChannel::Ds1ToDstarK0s>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs1Data, "Process data", true);

  void processDs2StarData(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates)
  {
    processData<DecayChannel::Ds2StarToDplusK0s>(collisions, candidates);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs2StarData, "Process data", false);

}; // struct HfTaskCharmResoReduced
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmResoReduced>(cfgc)};
}
