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

/// \file taskLcParametrizedPid.cxx
/// \brief Lc analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_lc_to_p_k_pi;

/// Fills MC histograms.
struct HfTaskLcParametrizedPid {
  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc_parametrized_pid::isSelLcToPKPiNoPid == 1 || aod::hf_sel_candidate_lc_parametrized_pid::isSelLcToPiKPNoPid == 1);

  HistogramRegistry registry{
    "registry",
    {{"hMassGen", "3-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigBkgLcNoPid", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcNoPid", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcNoPid", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLc", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLc", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLc", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLcPerfectPid", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcPerfectPid", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcPerfectPid", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  void process(soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLcParametrizedPid, aod::HfCand3ProngMcRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particlesMC,
               aod::TracksWMc const& tracks)
  {
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        continue;
      }
      if (std::abs(yLc(candidate)) > 4.0) {
        continue;
      }

      auto massLc = invMassLcToPKPi(candidate);
      auto massLcSwap = invMassLcToPiKP(candidate);
      auto ptCandidate = candidate.pt();
      auto rapidityCandidate = std::abs(yLc(candidate));

      if (candidate.isSelLcToPKPiNoPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcNoPid"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcToPKPiPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcNoPid"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcNoPid"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcToPiKPNoPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcNoPid"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcToPiKPPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcNoPid"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcNoPid"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcToPKPi() == 1) {
        registry.fill(HIST("hMassSigBkgLc"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcToPKPiPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLc"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLc"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcToPiKP() == 1) {
        registry.fill(HIST("hMassSigBkgLc"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcToPiKPPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLc"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLc"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcToPKPiPerfectPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcPerfectPid"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << DecayType::LcToPKPi)) {
          registry.fill(HIST("hMassSigLcPerfectPid"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcPerfectPid"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcToPiKPPerfectPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcPerfectPid"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << DecayType::LcToPKPi)) {
          registry.fill(HIST("hMassSigLcPerfectPid"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcPerfectPid"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }
    }

    for (const auto& particle : particlesMC) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::LcToPKPi) {
        if (std::abs(RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > 4.0) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        registry.fill(HIST("hMassGen"), ptGen, std::abs(yGen));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTaskLcParametrizedPid>(cfgc));
  return workflow;
}
