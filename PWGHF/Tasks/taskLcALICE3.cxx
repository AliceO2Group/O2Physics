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

/// \file taskLC.cxx
/// \brief LC analysis task
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
using namespace o2::aod::hf_cand_prong3;
using namespace o2::analysis::hf_cuts_lc_topkpi;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Fills MC histograms.
struct TaskLcALICE3MC {
  HistogramRegistry registry{
    "registry",
    {{"hMassGen", "3-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigBkgLcNoPID", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcNoPID", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcNoPID", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLcTOFPID", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcTOFPID", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcTOFPID", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLcTOFplusRICHPID", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcTOFplusRICHPID", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcTOFplusRICHPID", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgLcPerfectPID", "3-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigLcPerfectPID", "3-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgLcPerfectPID", "3-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{500, 1.6, 3.1}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  Filter filterSelectCandidates = (aod::hf_selcandidate_lc_ALICE3::isSelLcPKPiNoPID == 1 || aod::hf_selcandidate_lc_ALICE3::isSelLcPiKPNoPID == 1);
  void process(soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidateALICE3, aod::HfCandProng3MCRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandProng3MCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        continue;
      }
      if (std::abs(YLc(candidate)) > 4.0) {
        continue;
      }

      auto massLc = InvMassLcpKpi(candidate);
      auto massLcSwap = InvMassLcpiKp(candidate);
      auto ptCandidate = candidate.pt();
      auto rapidityCandidate = std::abs(YLc(candidate));

      if (candidate.isSelLcPKPiNoPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcNoPID"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcPKPiPerfectPID() == 1) {
          registry.fill(HIST("hMassSigLcNoPID"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcNoPID"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcPiKPNoPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcNoPID"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcPiKPPerfectPID() == 1) {
          registry.fill(HIST("hMassSigLcNoPID"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcNoPID"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcPKPiTOFPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcTOFPID"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcPKPiPerfectPID() == 1) {
          registry.fill(HIST("hMassSigLcTOFPID"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTOFPID"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcPiKPTOFPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcTOFPID"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcPiKPPerfectPID() == 1) {
          registry.fill(HIST("hMassSigLcTOFPID"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTOFPID"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcPKPiTOFplusRICHPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcTOFplusRICHPID"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcPKPiPerfectPID() == 1) {
          registry.fill(HIST("hMassSigLcTOFplusRICHPID"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTOFplusRICHPID"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcPiKPTOFplusRICHPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcTOFplusRICHPID"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi) && candidate.isSelLcPiKPPerfectPID() == 1) {
          registry.fill(HIST("hMassSigLcTOFplusRICHPID"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcTOFplusRICHPID"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelLcPKPiPerfectPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcPerfectPID"), massLc, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi)) {
          registry.fill(HIST("hMassSigLcPerfectPID"), massLc, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcPerfectPID"), massLc, ptCandidate, rapidityCandidate);
        }
      }
      if (candidate.isSelLcPiKPPerfectPID() == 1) {
        registry.fill(HIST("hMassSigBkgLcPerfectPID"), massLcSwap, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::LcToPKPi)) {
          registry.fill(HIST("hMassSigLcPerfectPID"), massLcSwap, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgLcPerfectPID"), massLcSwap, ptCandidate, rapidityCandidate);
        }
      }
    }

    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::LcToPKPi) {
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
  workflow.push_back(adaptAnalysisTask<TaskLcALICE3MC>(cfgc, TaskName{"hf-task-lc-ALICE3-mc"}));
  return workflow;
}
