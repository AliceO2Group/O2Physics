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
//#include "Common/Core/TrackSelection.h"
//#include "Common/DataModel/TrackSelectionTables.h"

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
struct TaskD0parametrizedPIDMC {
  HistogramRegistry registry{
    "registry",
    {{"hMassGen", "2-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigBkgD0NoPID", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0NoPID", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0NoPID", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassReflBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0PerfectPID", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0PerfectPID", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0PerfectPID", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  //Configurable<double> centralitySelectionMin{"centralitySelectionMin", 0.0, "Lower boundary of centrality selection"};
  //Configurable<double> centralitySelectionMax{"centralitySelectionMax", 30000.0, "Higher boundary of centrality selection"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_d0_parametrizedPID::isSelD0NoPID >= 1 || aod::hf_selcandidate_d0_parametrizedPID::isSelD0barNoPID >= 1);

  void process(soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0CandidateparametrizedPID, aod::HfCandProng2MCRec>> const& candidates, soa::Join<aod::McParticles, aod::HfCandProng2MCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  // void process(const o2::aod::Collision& collision, soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0CandidateparametrizedPID, aod::HfCandProng2MCRec>> const& candidates, soa::Join<aod::McParticles, aod::HfCandProng2MCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    //float ncontributor = collision.numContrib();
    for (auto& candidate : candidates) {
      /*if (ncontributor<=centralitySelectionMin && ncontributor>centralitySelectionMax) {
	continue;
	}*/
      if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (std::abs(YD0(candidate)) > 4.0) {
        continue;
      }

      auto massD0 = InvMassD0(candidate);
      auto massD0bar = InvMassD0bar(candidate);
      auto ptCandidate = candidate.pt();
      auto rapidityCandidate = std::abs(YD0(candidate));

      if (candidate.isSelD0NoPID() >= 1) {
        registry.fill(HIST("hMassSigBkgD0NoPID"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0NoPID"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0NoPID"), massD0, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0() >= 1) {
        registry.fill(HIST("hMassSigBkgD0"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0"), massD0, ptCandidate, rapidityCandidate);
        } else {
          if (candidate.flagMCMatchRec() == -(1 << DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0"), massD0, ptCandidate, rapidityCandidate);
          }
          registry.fill(HIST("hMassBkgD0"), massD0, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0PerfectPID() >= 1) {
        registry.fill(HIST("hMassSigBkgD0PerfectPID"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0PerfectPID"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0PerfectPID"), massD0, ptCandidate, rapidityCandidate);
        }
      }
    }

    for (auto& particle : particlesMC) {
      /*if (ncontributor<=centralitySelectionMin && ncontributor>centralitySelectionMax) {
	continue;
      }*/
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::D0ToPiK) {
        if (std::abs(RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > 4.0) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        registry.fill(HIST("hMassGen"), ptGen, std::abs(yGen));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<TaskD0parametrizedPIDMC>(cfgc, TaskName{"hf-task-d0-parametrizedPID-mc"}));
  return workflow;
}
