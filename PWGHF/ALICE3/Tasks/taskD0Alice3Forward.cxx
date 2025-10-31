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

/// \file taskD0Alice3Forward.cxx
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
struct HfTaskD0Alice3Forward {
  Filter filterSelectCandidates = (aod::hf_sel_candidate_d0_alice3_forward::isSelHfFlag >= 1);

  HistogramRegistry registry{
    "registry",
    {{"hMassGen", "2-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigBkgD0ForwardRICHPID", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0ForwardRICHPID", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0ForwardRICHPID", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  void process(soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0Alice3Forward, aod::HfCand2ProngMcRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
               aod::TracksWMc const&)
  {
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (std::abs(HfHelper::yD0(candidate)) > 4.0) {
        continue;
      }

      auto massD0 = HfHelper::invMassD0ToPiK(candidate);
      // auto massD0bar = HfHelper::invMassD0barToKPi(candidate);
      auto ptCandidate = candidate.pt();
      auto rapidityCandidate = std::abs(HfHelper::yD0(candidate));

      if (candidate.isSelD0FRichPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0ForwardRICHPID"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0ForwardRICHPID"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0ForwardRICHPID"), massD0, ptCandidate, rapidityCandidate);
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
  workflow.push_back(adaptAnalysisTask<HfTaskD0Alice3Forward>(cfgc));
  return workflow;
}
