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

/// \file taskLcAlice3.cxx
/// \brief Lc analysis task for ALICE 3
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
struct HfTaskLcAlice3 {
  Filter filterSelectCandidates = (aod::hf_sel_candidate_lc_alice3::isSelLcToPKPiNoPid == 1 || aod::hf_sel_candidate_lc_alice3::isSelLcToPiKPNoPid == 1);

  HistogramRegistry registry{
    "registry",
    {{"hMassGen", "3-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigBkgLcNoPid", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcNoPid", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcNoPid", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLcTofPid", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcTofPid", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcTofPid", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLcTofPlusRichPid", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcTofPlusRichPid", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcTofPlusRichPid", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLcPerfectPid", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcPerfectPid", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcPerfectPid", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  void process(soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLcAlice3, aod::HfCand3ProngMcRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles,
               aod::TracksWMc const&)
  {
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      if (std::abs(HfHelper::yLc(candidate)) > 4.0) {
        continue;
      }

      auto massLc = HfHelper::invMassLcToPKPi(candidate);
      auto massLcSwap = HfHelper::invMassLcToPiKP(candidate);
      auto ptCandidate = candidate.pt();
      auto rapidityCandidate = std::abs(HfHelper::yLc(candidate));

      if (candidate.isSelLcToPKPiNoPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcNoPid"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi) && candidate.isSelLcToPKPiPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcNoPid"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcNoPid"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcToPiKPNoPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcNoPid"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi) && candidate.isSelLcToPiKPPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcNoPid"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcNoPid"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcToPKPiTofPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcTofPid"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi) && candidate.isSelLcToPKPiPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcTofPid"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTofPid"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcToPiKPTofPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcTofPid"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi) && candidate.isSelLcToPiKPPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcTofPid"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTofPid"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcToPKPiTofPlusRichPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcTofPlusRichPid"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi) && candidate.isSelLcToPKPiPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcTofPlusRichPid"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTofPlusRichPid"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcToPiKPTofPlusRichPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcTofPlusRichPid"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi) && candidate.isSelLcToPiKPPerfectPid() == 1) {
          registry.fill(HIST("hMassSigLcTofPlusRichPid"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTofPlusRichPid"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcToPKPiPerfectPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcPerfectPid"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          registry.fill(HIST("hMassSigLcPerfectPid"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcPerfectPid"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcToPiKPPerfectPid() == 1) {
        registry.fill(HIST("hMassSigBkgLcPerfectPid"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          registry.fill(HIST("hMassSigLcPerfectPid"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcPerfectPid"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }
    }

    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        if (std::abs(RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus)) > 4.0) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus);
        registry.fill(HIST("hMassGen"), ptGen, std::abs(yGen));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTaskLcAlice3>(cfgc));
  return workflow;
}
