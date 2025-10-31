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

/// \file taskD0Alice3Barrel.cxx
/// \brief D0 analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Fills MC histograms.
struct HfTaskD0Alice3Barrel {
  // Configurable<double> centralitySelectionMin{"centralitySelectionMin", 0.0, "Lower boundary of centrality selection"};
  // Configurable<double> centralitySelectionMax{"centralitySelectionMax", 0.0, "Higher boundary of centrality selection"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_d0_alice3_barrel::isSelHfFlag >= 1);

  HistogramRegistry registry{
    "registry",
    {{"hMassGen", "2-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigBkgD0NoPid", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0NoPid", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0NoPid", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0TofPid", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0TofPid", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0TofPid", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0RICHPID", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0RICHPID", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0RICHPID", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassReflBkgD0RICHPID", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0TofPlusRichPid", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0TofPlusRichPid", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0TofPlusRichPid", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0barTofPlusRichPid", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0barTofPlusRichPid", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0barTofPlusRichPid", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0PerfectPid", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0PerfectPid", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0PerfectPid", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  // void process(soa::Join<aod::Collisions, aod::CentV0Ms>::iterator const& collision, soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0Alice3Barrel, aod::HfCand2ProngMcRec>> const& candidates, soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles, aod::TracksWMc const& tracks)
  void process(soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0Alice3Barrel, aod::HfCand2ProngMcRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
               aod::TracksWMc const&)
  {
    // float centrality = collision.centV0M();
    for (const auto& candidate : candidates) {
      // if (centrality<=centralitySelectionMin && centrality>centralitySelectionMax) {
      // continue;
      // }
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (std::abs(HfHelper::yD0(candidate)) > 4.0) {
        continue;
      }

      auto massD0 = HfHelper::invMassD0ToPiK(candidate);
      auto massD0bar = HfHelper::invMassD0barToKPi(candidate);
      auto ptCandidate = candidate.pt();
      auto rapidityCandidate = std::abs(HfHelper::yD0(candidate));

      if (candidate.isSelD0NoPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0NoPid"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0NoPid"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0NoPid"), massD0, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0TofPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0TofPid"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0TofPid"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0TofPid"), massD0, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0RichPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0RICHPID"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0RICHPID"), massD0, ptCandidate, rapidityCandidate);
        } else {
          if (candidate.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0RICHPID"), massD0, ptCandidate, rapidityCandidate);
          }
          registry.fill(HIST("hMassBkgD0RICHPID"), massD0, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0TofPlusRichPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0TofPlusRichPid"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0TofPlusRichPid"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0TofPlusRichPid"), massD0, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0barTofPlusRichPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0barTofPlusRichPid"), massD0bar, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0barTofPlusRichPid"), massD0bar, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0barTofPlusRichPid"), massD0bar, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0PerfectPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0PerfectPid"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0PerfectPid"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0PerfectPid"), massD0, ptCandidate, rapidityCandidate);
        }
      }
    }

    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
        if (std::abs(RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0)) > 4.0) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0);
        registry.fill(HIST("hMassGen"), ptGen, std::abs(yGen));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTaskD0Alice3Barrel>(cfgc));
  return workflow;
}
