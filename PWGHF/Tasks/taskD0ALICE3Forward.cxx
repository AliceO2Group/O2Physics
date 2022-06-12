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

/// \file taskD0.cxx
/// \brief D0 analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::analysis::hf_cuts_d0_topik;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Fills MC histograms.
struct TaskD0ALICE3ForwardMC {
  HistogramRegistry registry{
    "registry",
    {{"hMassGen", "2-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigBkgD0ForwardRICHPID", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0ForwardRICHPID", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0ForwardRICHPID", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  Filter filterSelectCandidates = (aod::hf_selcandidate_d0_ALICE3_Forward::isSelHFFFlag >= 1);

  void process(soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0CandidateALICE3Forward, aod::HfCandProng2MCRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandProng2MCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (std::abs(YD0(candidate)) > 4.0) {
        continue;
      }

      auto massD0 = InvMassD0(candidate);
      //auto massD0bar = InvMassD0bar(candidate);
      auto ptCandidate = candidate.pt();
      auto rapidityCandidate = std::abs(YD0(candidate));

      if (candidate.isSelD0FRICHPID() >= 1) {
        registry.fill(HIST("hMassSigBkgD0ForwardRICHPID"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0ForwardRICHPID"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0ForwardRICHPID"), massD0, ptCandidate, rapidityCandidate);
        }
      }
    }

    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::D0ToPiK) {
        if (std::abs(RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > 4.0) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        registry.fill(HIST("hMassGen"), ptGen, std::abs(yGen));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<TaskD0ALICE3ForwardMC>(cfgc, TaskName{"hf-task-d0-ALICE3-Forward-mc"}));
  return workflow;
}
