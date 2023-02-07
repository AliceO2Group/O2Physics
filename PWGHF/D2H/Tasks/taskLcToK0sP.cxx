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

/// \file taskLcToK0sP.cxx
/// \brief Lc -> K0S+p analysis task
///
/// \author Chiara Zampolli, <Chiara.Zampolli@cern.ch>, CERN
///         Paul Buehler, <paul.buehler@oeaw.ac.at>, Vienna
///
/// based on taskD0.cxx, taskLc.cxx

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc;
using namespace o2::framework::expressions;

/// LcToK0sp analysis task
struct HfTaskLcToK0sP {
  Configurable<int> selectionFlagLcToK0sP{"selectionFlagLcToK0sP", 1, "Selection Flag for Lc"};
  Configurable<int> selectionFlagLcbarToK0sP{"selectionFlagLcbarToK0sP", 1, "Selection Flag for Lcbar"};
  Configurable<double> etaCandMax{"etaCandMax", -1., "max. cand. pseudorapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcToK0sP || aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcbarToK0sP);

  HistogramRegistry registry{"registry"};

  void init(InitContext& context)
  {
    // data
    registry.add("hMassVsPtCand", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{600, 1.98f, 2.58f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtCand", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
    registry.add("hPtBachVsPtCand", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtV0VsPtCand", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0BachVsPtCand", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0V0posVsPtCand", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0V0negVsPtCand", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0CPAVsPtCand", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98f, 1.0001f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaVsPtCand", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {{500, -2.0f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "cascade candidates;selection status;p_{T}", {HistType::kTH1F, {{5, -0.5f, 4.5f}}});
    // add MC histograms
    if (context.mOptions.get<bool>("processMc")) {
      registry.add("hPtRecSig", "cascade candidates (MC);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hPtRecBg", "cascade candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hPtGen", "cascade (MC);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hPtGenSig", "cascade candidates (MC);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hCPAVsPtCandRecSig", "cascade candidates (matched);cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hCPAVsPtCandRecBg", "cascade candidates (unmatched);cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hV0CPAVsPtCandRecSig", "cascade candidates (matched);v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98, 1.0001}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hV0CPAVsPtCandRecBg", "cascade candidates (unmatched);v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98, 1.0001}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hEtaVsPtCandRecSig", "cascade candidates (matched);#it{#eta};p_{T}", {HistType::kTH2F, {{500, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hEtaVsPtCandRecBg", "cascade candidates (unmatched);#it{#eta};p_{T}", {HistType::kTH2F, {{500, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hEtaGen", "MC particles (MC);#it{#eta};entries", {HistType::kTH1F, {{500, -2., 2.}}});
    }
  }

  void
    process(soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP>> const& candidates)
  {
    // Printf("Candidates: %d", candidates.size());
    for (auto& candidate : candidates) {
      /*
      // no such selection for LcK0sp for now - it is the only cascade
      if (!(candidate.hfflag() & 1 << D0ToPiK)) {
        continue;
      }
      */
      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        // Printf("Candidate: eta rejection: %g", candidate.eta());
        continue;
      }

      auto ptCand = candidate.pt();
      registry.fill(HIST("hMassVsPtCand"), invMassLcToK0sP(candidate), ptCand);
      registry.fill(HIST("hPtCand"), ptCand);
      registry.fill(HIST("hPtBachVsPtCand"), candidate.ptProng0(), ptCand);
      registry.fill(HIST("hPtV0VsPtCand"), candidate.ptProng1(), ptCand);
      registry.fill(HIST("hd0BachVsPtCand"), candidate.impactParameter0(), ptCand);
      registry.fill(HIST("hd0V0posVsPtCand"), candidate.dcapostopv(), ptCand);
      registry.fill(HIST("hd0V0negVsPtCand"), candidate.dcanegtopv(), ptCand);
      registry.fill(HIST("hV0CPAVsPtCand"), candidate.v0cosPA(), ptCand);
      registry.fill(HIST("hEtaVsPtCand"), candidate.eta(), ptCand);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelLcToK0sP());
    }
  }

  void processMc(soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP, aod::HfCandCascadeMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandCascadeMcGen> const& particlesMC,
                 aod::BigTracksMC const& tracks)
  {
    // MC rec.
    // Printf("MC Candidates: %d", candidates.size());
    for (auto& candidate : candidates) {
      if (etaCandMax >= 0. && std::abs(candidate.eta()) > etaCandMax) {
        // Printf("MC Rec.: eta rejection: %g", candidate.eta());
        continue;
      }
      auto ptCand = candidate.pt();
      if (std::abs(candidate.flagMcMatchRec()) == 1) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandCascadeMcGen>>(), pdg::Code::kLambdaCPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        registry.fill(HIST("hPtRecSig"), ptCand);              // rec. level pT
        registry.fill(HIST("hCPAVsPtCandRecSig"), candidate.cpa(), ptCand);
        registry.fill(HIST("hV0CPAVsPtCandRecSig"), candidate.v0cosPA(), ptCand);
        registry.fill(HIST("hEtaVsPtCandRecSig"), candidate.eta(), ptCand);
      } else {
        registry.fill(HIST("hPtRecBg"), ptCand);
        registry.fill(HIST("hCPAVsPtCandRecBg"), candidate.cpa(), ptCand);
        registry.fill(HIST("hV0CPAVsPtCandRecBg"), candidate.v0cosPA(), ptCand);
        registry.fill(HIST("hEtaVsPtCandRecBg"), candidate.eta(), ptCand);
      }
    }
    // MC gen.
    // Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (etaCandMax >= 0. && std::abs(particle.eta()) > etaCandMax) {
        // Printf("MC Gen.: eta rejection: %g", particle.eta());
        continue;
      }
      if (std::abs(particle.flagMcMatchGen()) == 1) {
        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(HfTaskLcToK0sP, processMc, "Process MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskLcToK0sP>(cfgc),
  };
}
