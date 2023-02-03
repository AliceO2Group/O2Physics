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

  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcToK0sP || aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcbarToK0sP);

  HistogramRegistry registry{
    "registry",
    {// data
     {"hMass", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 1.5f, 3.0f}}}},
     {"hPtCand", "cascade candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
     {"hPtBach", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
     {"hPtV0", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
     {"hd0Bach", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{100, -1.0f, 1.0f}}}},
     {"hd0V0pos", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"hd0V0neg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{200, -5.0f, 5.0f}}}},
     {"hV0CPA", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {{100, 0.98f, 1.0001f}}}},
     {"hEta", "cascade candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{100, -2.0f, 2.0f}}}},
     {"hSelectionStatus", "cascade candidates;selection status;entries", {HistType::kTH1F, {{5, -0.5f, 4.5f}}}}}};

  void init(InitContext& context)
  {
    // add MC histograms
    if (context.mOptions.get<bool>("processMc")) {
      registry.add("hPtRecSig", "cascade candidates (MC);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hPtRecBg", "cascade candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hPtGen", "cascade (MC);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hPtGenSig", "cascade candidates (MC);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
      registry.add("hCPARecSig", "cascade candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("hCPARecBg", "cascade candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
      registry.add("hV0CPARecSig", "cascade candidates (matched);v0 cosine of pointing angle;entries", {HistType::kTH1F, {{100, 0.98, 1.0001}}});
      registry.add("hV0CPARecBg", "cascade candidates (unmatched);v0 cosine of pointing angle;entries", {HistType::kTH1F, {{100, 0.98, 1.0001}}});
      registry.add("hEtaRecSig", "cascade candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("hEtaRecBg", "cascade candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
      registry.add("hEtaGen", "MC particles (MC);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}});
    }
  }

  void process(soa::Filtered<soa::Join<aod::HfCandCascExt, aod::HfSelLcToK0sP>> const& candidates)
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

      registry.fill(HIST("hMass"), invMassLcToK0sP(candidate));
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtBach"), candidate.ptProng0());
      registry.fill(HIST("hPtV0"), candidate.ptProng1());
      registry.fill(HIST("hd0Bach"), candidate.impactParameter0());
      registry.fill(HIST("hd0V0pos"), candidate.dcapostopv());
      registry.fill(HIST("hd0V0neg"), candidate.dcanegtopv());
      registry.fill(HIST("hV0CPA"), candidate.v0cosPA());
      registry.fill(HIST("hEta"), candidate.eta());
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
      if (std::abs(candidate.flagMcMatchRec()) == 1) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandCascadeMcGen>>(), pdg::Code::kLambdaCPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        registry.fill(HIST("hPtRecSig"), candidate.pt());      // rec. level pT
        registry.fill(HIST("hCPARecSig"), candidate.cpa());
        registry.fill(HIST("hV0CPARecSig"), candidate.v0cosPA());
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa());
        registry.fill(HIST("hV0CPARecBg"), candidate.v0cosPA());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
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
