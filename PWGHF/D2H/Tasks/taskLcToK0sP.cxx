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
    registry.add("hMass", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {{600, 1.98f, 2.58f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtCand", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
    registry.add("hPtBach", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtV0", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Bach", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{500, -0.5f, 0.5f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0V0pos", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0V0neg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {{1000, -5.0f, 5.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hV0CPA", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98f, 1.0001f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "cascade candidates;candidate #it{#eta};p_{T}", {HistType::kTH2F, {{500, -2.0f, 2.0f}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "cascade candidates;selection status;p_{T}", {HistType::kTH1F, {{5, -0.5f, 4.5f}}});
    // add MC histograms
    if (context.mOptions.get<bool>("processMc")) {
      registry.add("hPtRecSig", "cascade candidates (MC);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hPtRecBg", "cascade candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hPtGen", "cascade (MC);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hPtGenSig", "cascade candidates (MC);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{300, 0.0f, 30.0f}}});
      registry.add("hCPARecSig", "cascade candidates (matched);cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hCPARecBg", "cascade candidates (unmatched);cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hV0CPARecSig", "cascade candidates (matched);v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98, 1.0001}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hV0CPARecBg", "cascade candidates (unmatched);v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {{500, 0.98, 1.0001}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hEtaRecSig", "cascade candidates (matched);#it{#eta};p_{T}", {HistType::kTH2F, {{500, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("hEtaRecBg", "cascade candidates (unmatched);#it{#eta};p_{T}", {HistType::kTH2F, {{500, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
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

      auto candPt = candidate.pt();
      registry.fill(HIST("hMass"), invMassLcToK0sP(candidate), candPt);
      registry.fill(HIST("hPtCand"), candPt);
      registry.fill(HIST("hPtBach"), candidate.ptProng0(), candPt);
      registry.fill(HIST("hPtV0"), candidate.ptProng1(), candPt);
      registry.fill(HIST("hd0Bach"), candidate.impactParameter0(), candPt);
      registry.fill(HIST("hd0V0pos"), candidate.dcapostopv(), candPt);
      registry.fill(HIST("hd0V0neg"), candidate.dcanegtopv(), candPt);
      registry.fill(HIST("hV0CPA"), candidate.v0cosPA(), candPt);
      registry.fill(HIST("hEta"), candidate.eta(), candPt);
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
      auto candPt = candidate.pt();
      if (std::abs(candidate.flagMcMatchRec()) == 1) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandCascadeMcGen>>(), pdg::Code::kLambdaCPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        registry.fill(HIST("hPtRecSig"), candPt);              // rec. level pT
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candPt);
        registry.fill(HIST("hV0CPARecSig"), candidate.v0cosPA(), candPt);
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candPt);
      } else {
        registry.fill(HIST("hPtRecBg"), candPt);
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candPt);
        registry.fill(HIST("hV0CPARecBg"), candidate.v0cosPA(), candPt);
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candPt);
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
