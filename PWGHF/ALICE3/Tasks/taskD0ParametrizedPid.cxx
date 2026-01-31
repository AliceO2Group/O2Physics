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

/// \file taskD0ParametrizedPid.cxx
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
struct HfTaskD0ParametrizedPid {
  // Configurable<double> centralitySelectionMin{"centralitySelectionMin", 0.0, "Lower boundary of centrality selection"};
  // Configurable<double> centralitySelectionMax{"centralitySelectionMax", 30000.0, "Higher boundary of centrality selection"};

  using McParticlesHf = soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_d0_parametrized_pid::isSelD0NoPid >= 1 || aod::hf_sel_candidate_d0_parametrized_pid::isSelD0barNoPid >= 1);

  HistogramRegistry registry{
    "registry",
    {{"hGenPtVsY", "2-prong candidates (generated); #it{p}_{T}; #it{y}", {HistType::kTH2F, {{150, 0., 30.}, {8, 0, 4.0}}}},
     {"hGenPtCMSYCutNoDaughterEtaCut", "2-prong candidates (generated); #it{p}_{T}", {HistType::kTH1F, {{150, 0., 30.}}}},
     {"hGenPtFiducialYCutNoDaughterEtaCut", "2-prong candidates (generated); #it{p}_{T}", {HistType::kTH1F, {{150, 0., 30.}}}},
     {"hGenPtCMSYCutCMSDaughterEtaCut", "2-prong candidates (generated); #it{p}_{T}", {HistType::kTH1F, {{150, 0., 30.}}}},
     {"hGenPtFiducialYCutALICE2DaughterEtaCut", "2-prong candidates (generated); #it{p}_{T}", {HistType::kTH1F, {{150, 0., 30.}}}},
     {"hGenPtFiducialYCutCMSDaughterEtaCut", "2-prong candidates (generated); #it{p}_{T}", {HistType::kTH1F, {{150, 0., 30.}}}},
     {"hMassSigBkgD0NoPid", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0NoPid", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0NoPid", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassReflBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassSigBkgD0PerfectPid", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0, 4.0}}}},
     {"hMassSigD0PerfectPid", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}},
     {"hMassBkgD0PerfectPid", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {8, 0., 4.}}}}}};

  void process(soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0ParametrizedPid, aod::HfCand2ProngMcRec>> const& candidates,
               McParticlesHf const& mcParticles,
               aod::TracksWMc const&)
  // void process(const o2::aod::Collision& collision, soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0ParametrizedPid, aod::HfCand2ProngMcRec>> const& candidates, soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles, aod::TracksWMc const& tracks)
  {
    // float ncontributor = collision.numContrib();
    for (const auto& candidate : candidates) {
      // if (ncontributor<=centralitySelectionMin && ncontributor>centralitySelectionMax) {
      //   continue;
      // }
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

      if (candidate.isSelD0NoPid() >= 1) {
        registry.fill(HIST("hMassSigBkgD0NoPid"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0NoPid"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0NoPid"), massD0, ptCandidate, rapidityCandidate);
        }
      }

      if (candidate.isSelD0() >= 1) {
        registry.fill(HIST("hMassSigBkgD0"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0"), massD0, ptCandidate, rapidityCandidate);
        } else {
          if (candidate.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0"), massD0, ptCandidate, rapidityCandidate);
          }
          registry.fill(HIST("hMassBkgD0"), massD0, ptCandidate, rapidityCandidate);
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
      // if (ncontributor<=centralitySelectionMin && ncontributor>centralitySelectionMax) {
      //   continue;
      // }
      float maxFiducialY = 0.8;
      float minFiducialY = -0.8;
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
        if (std::abs(RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0)) > 4.0) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassD0);
        registry.fill(HIST("hGenPtVsY"), ptGen, std::abs(yGen));
        if (ptGen < 5.0) {
          maxFiducialY = -0.2 / 15 * ptGen * ptGen + 1.9 / 15 * ptGen + 0.5;
          minFiducialY = 0.2 / 15 * ptGen * ptGen - 1.9 / 15 * ptGen - 0.5;
        }
        auto etaProng0 = std::abs(particle.daughters_as<McParticlesHf>().begin().eta());
        auto etaProng1 = std::abs((particle.daughters_as<McParticlesHf>().begin() + 1).eta());
        if (std::abs(yGen) < 2.5) {
          registry.fill(HIST("hGenPtCMSYCutNoDaughterEtaCut"), ptGen);
        }

        if (yGen > minFiducialY && yGen < maxFiducialY) {
          registry.fill(HIST("hGenPtFiducialYCutNoDaughterEtaCut"), ptGen);
        }
        if (etaProng0 < 2.5 && etaProng1 < 2.5) {
          if (std::abs(yGen) < 2.5) {
            registry.fill(HIST("hGenPtCMSYCutCMSDaughterEtaCut"), ptGen);
          }
          if (yGen > minFiducialY && yGen < maxFiducialY) {
            registry.fill(HIST("hGenPtFiducialYCutCMSDaughterEtaCut"), ptGen);
            if (etaProng0 < 0.8 && etaProng1 < 0.8) {
              registry.fill(HIST("hGenPtFiducialYCutALICE2DaughterEtaCut"), ptGen);
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTaskD0ParametrizedPid>(cfgc));
  return workflow;
}
