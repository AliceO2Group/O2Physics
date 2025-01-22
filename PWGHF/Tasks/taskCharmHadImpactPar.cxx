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

/// \file taskCharmHadImpactPar.cxx
/// \brief Analysis task to produce impact-parameter distributions of charm hadrons
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum Channel : uint8_t {
  DplusToKPiPi = 0,
  DzeroToKPi,
  NChannels
};

struct HfTaskCharmHadImpactPar {
  Configurable<int> selectionFlag{"selectionFlag", 15, "Selection Flag for the considered charm hadron"};
  ConfigurableAxis axisPt{"axisPt", {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 10.f, 12.f, 16.f, 24.f, 36.f, 50.f}, "axis for pT of charm hadron"};
  ConfigurableAxis axisMass{"axisMass", {250, 1.65f, 2.15f}, "axis for mass of charm hadron"};
  ConfigurableAxis axisImpPar{"axisImpPar", {2000, -500.f, 500.f}, "axis for impact-parameter of charm hadron"};
  ConfigurableAxis axisMlScore0{"axisMlScore0", {100, 0.f, 1.f}, "axis for ML output score 0"};
  ConfigurableAxis axisMlScore1{"axisMlScore1", {100, 0.f, 1.f}, "axis for ML output score 1"};
  ConfigurableAxis axisMlScore2{"axisMlScore2", {100, 0.f, 1.f}, "axis for ML output score 2"};

  HfHelper hfHelper;

  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandDplusDataWithMl = soa::Filtered<soa::Join<CandDplusData, aod::HfMlDplusToPiKPi>>;
  using CandDzeroData = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using CandDzeroDataWithMl = soa::Filtered<soa::Join<CandDzeroData, aod::HfMlD0>>;

  Filter filterDplusFlag = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterDzeroFlag = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 4> doprocess{doprocessDplus, doprocessDplusWithMl, doprocessDzero, doprocessDzeroWithMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
    }
    if (doprocessDplus || doprocessDzero) {
      registry.add("hMassPtImpPar", ";#it{M} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});dca XY (#mum);", HistType::kTHnSparseF, {axisMass, axisPt, axisImpPar});
    } else if (doprocessDplusWithMl || doprocessDzeroWithMl) {
      registry.add("hMassPtImpPar", ";#it{M} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});dca XY (#mum);ML score 0;ML score 1; ML score 2;", HistType::kTHnSparseF, {axisMass, axisPt, axisImpPar, axisMlScore0, axisMlScore1, axisMlScore2});
    }
  }

  // Fill THnSparses for the ML analysis
  /// \param candidate is a particle candidate
  template <Channel channel, bool withMl, typename CCands>
  void fillSparse(const CCands& candidate)
  {
    std::vector<float> outputMl = {-999., -999., -999.};
    float invMass{-1.f};
    if constexpr (channel == Channel::DplusToKPiPi) { // D+ -> Kpipi
      invMass = hfHelper.invMassDplusToPiKPi(candidate);
      if constexpr (withMl) {
        for (auto iScore{0u}; iScore < candidate.mlProbDplusToPiKPi().size(); ++iScore) {
          outputMl[iScore] = candidate.mlProbDplusToPiKPi()[iScore];
        }
        registry.fill(HIST("hMassPtImpPar"), invMass, candidate.pt(), candidate.impactParameterXY(), outputMl[0], outputMl[1], outputMl[2]);
      } else {
        registry.fill(HIST("hMassPtImpPar"), invMass, candidate.pt(), candidate.impactParameterXY());
      }
    } else if constexpr (channel == Channel::DzeroToKPi) {
      if (candidate.isSelD0()) { // D0 -> Kpi
        invMass = hfHelper.invMassD0ToPiK(candidate);
        if constexpr (withMl) {
          for (auto iScore{0u}; iScore < candidate.mlProbD0().size(); ++iScore) {
            outputMl[iScore] = candidate.mlProbD0()[iScore];
          }
          registry.fill(HIST("hMassPtImpPar"), invMass, candidate.pt(), candidate.impactParameterXY(), outputMl[0], outputMl[1], outputMl[2]);
        } else {
          registry.fill(HIST("hMassPtImpPar"), invMass, candidate.pt(), candidate.impactParameterXY());
        }
      }
      if (candidate.isSelD0bar()) {
        invMass = hfHelper.invMassD0barToKPi(candidate);
        if constexpr (withMl) {
          for (auto iScore{0u}; iScore < candidate.mlProbD0bar().size(); ++iScore) {
            outputMl[iScore] = candidate.mlProbD0bar()[iScore];
          }
          registry.fill(HIST("hMassPtImpPar"), invMass, candidate.pt(), candidate.impactParameterXY(), outputMl[0], outputMl[1], outputMl[2]);
        } else {
          registry.fill(HIST("hMassPtImpPar"), invMass, candidate.pt(), candidate.impactParameterXY());
        }
      }
    }
  }

  /// \param candidates are reconstructed candidates
  template <Channel channel, bool withMl, typename CCands>
  void runAnalysis(const CCands& candidates)
  {
    for (auto const& candidate : candidates) {
      fillSparse<channel, withMl>(candidate);
    }
  }

  // process functions
  void processDplus(CandDplusData const& candidates)
  {
    runAnalysis<Channel::DplusToKPiPi, false>(candidates);
  }
  PROCESS_SWITCH(HfTaskCharmHadImpactPar, processDplus, "Process D+ w/o ML", false);

  void processDplusWithMl(CandDplusDataWithMl const& candidates)
  {
    runAnalysis<Channel::DplusToKPiPi, true>(candidates);
  }
  PROCESS_SWITCH(HfTaskCharmHadImpactPar, processDplusWithMl, "Process D+ with ML", true);

  void processDzero(CandDzeroData const& candidates)
  {
    runAnalysis<Channel::DzeroToKPi, false>(candidates);
  }
  PROCESS_SWITCH(HfTaskCharmHadImpactPar, processDzero, "Process D0 w/o ML", false);

  void processDzeroWithMl(CandDzeroDataWithMl const& candidates)
  {
    runAnalysis<Channel::DzeroToKPi, true>(candidates);
  }
  PROCESS_SWITCH(HfTaskCharmHadImpactPar, processDzeroWithMl, "Process D0 with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmHadImpactPar>(cfgc)};
}
