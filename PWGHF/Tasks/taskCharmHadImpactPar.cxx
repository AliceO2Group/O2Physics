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
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum Channel : uint8_t {
  DplusToKPiPi = 0,
  DzeroToKPi,
  NChannels
};

namespace o2::aod
{
namespace hf_charm_cand_lite
{
DECLARE_SOA_COLUMN(M, m, float);                                 //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                               //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                 //! Rapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                             //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(ImpactParameterXY, impactParameterXY, float); //! Dca XY of candidate
DECLARE_SOA_COLUMN(MlScoreBkg, mlScoreBkg, float);               //! ML score for background class
DECLARE_SOA_COLUMN(MlScorePrompt, mlScorePrompt, float);         //! ML Prompt score for prompt class
DECLARE_SOA_COLUMN(MlScoreNonPrompt, mlScoreNonPrompt, float);   //! ML Non Prompt score for non prompt class
} // namespace hf_charm_cand_lite

DECLARE_SOA_TABLE(HfCharmCandLites, "AOD", "HFCHARMCANDLITE", //! Table with some B+ properties
                  hf_charm_cand_lite::M,
                  hf_charm_cand_lite::Pt,
                  hf_charm_cand_lite::Y,
                  hf_charm_cand_lite::Phi,
                  hf_charm_cand_lite::ImpactParameterXY,
                  hf_charm_cand_lite::MlScoreBkg,
                  hf_charm_cand_lite::MlScorePrompt,
                  hf_charm_cand_lite::MlScoreNonPrompt);
} // namespace o2::aod

struct HfTaskCharmHadImpactPar {
  Produces<aod::HfCharmCandLites> hfCharmCandLite;

  Configurable<int> selectionFlag{"selectionFlag", 15, "Selection Flag for the considered charm hadron"};
  Configurable<int> fillLightTreeCandidate{"fillLightTreeCandidate", 0, "Flag to store charm hadron features"};
  ConfigurableAxis axisPt{"axisPt", {0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 10.f, 12.f, 16.f, 24.f, 36.f, 50.f}, "axis for pT of charm hadron"};
  ConfigurableAxis axisMass{"axisMass", {250, 1.65f, 2.15f}, "axis for mass of charm hadron"};
  ConfigurableAxis axisPhi{"axisPhi", {180, 0.f, 2 * PI}, "axis for azimuthal angle of charm hadron"};
  ConfigurableAxis axisY{"axisY", {20, -1.f, 1.f}, "axis for rapidity of charm hadron"};
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
      registry.add("hMassPtImpParPhiY", ";#it{M} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});dca XY (#mum); phi; y;", HistType::kTHnSparseF, {axisMass, axisPt, axisImpPar, axisPhi, axisY});
    } else if (doprocessDplusWithMl || doprocessDzeroWithMl) {
      registry.add("hMassPtImpParPhiY", ";#it{M} (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});dca XY (#mum); phi; y; ML score 0;ML score 1; ML score 2;", HistType::kTHnSparseF, {axisMass, axisPt, axisImpPar, axisPhi, axisY, axisMlScore0, axisMlScore1, axisMlScore2});
    }
  }

  // Fill THnSparses for the ML analysis
  /// \param candidate is a particle candidate
  template <Channel channel, bool withMl, typename CCands>
  void fillSparse(const CCands& candidate)
  {
    std::vector<float> outputMl = {-999., -999., -999.};
    float invMass{-1.f};
    float yCand{-999.f};
    if constexpr (channel == Channel::DplusToKPiPi) { // D+ -> Kpipi
      invMass = hfHelper.invMassDplusToPiKPi(candidate);
      yCand = hfHelper.yDplus(candidate);
      if constexpr (withMl) {
        for (auto iScore{0u}; iScore < candidate.mlProbDplusToPiKPi().size(); ++iScore) {
          outputMl[iScore] = candidate.mlProbDplusToPiKPi()[iScore];
        }
        registry.fill(HIST("hMassPtImpParPhiY"), invMass, candidate.pt(), candidate.impactParameterXY(), candidate.phi(), yCand, outputMl[0], outputMl[1], outputMl[2]);
      } else {
        registry.fill(HIST("hMassPtImpParPhiY"), invMass, candidate.pt(), candidate.impactParameterXY(), candidate.phi(), yCand);
      }
    } else if constexpr (channel == Channel::DzeroToKPi) {
      if (candidate.isSelD0()) { // D0 -> Kpi
        invMass = hfHelper.invMassD0ToPiK(candidate);
        yCand = hfHelper.yD0(candidate);
        if constexpr (withMl) {
          for (auto iScore{0u}; iScore < candidate.mlProbD0().size(); ++iScore) {
            outputMl[iScore] = candidate.mlProbD0()[iScore];
          }
          registry.fill(HIST("hMassPtImpParPhiY"), invMass, candidate.pt(), candidate.impactParameterXY(), candidate.phi(), yCand, outputMl[0], outputMl[1], outputMl[2]);
        } else {
          registry.fill(HIST("hMassPtImpParPhiY"), invMass, candidate.pt(), candidate.impactParameterXY(), candidate.phi(), yCand);
        }
      }
      if (candidate.isSelD0bar()) {
        invMass = hfHelper.invMassD0barToKPi(candidate);
        yCand = hfHelper.yD0(candidate);
        if constexpr (withMl) {
          for (auto iScore{0u}; iScore < candidate.mlProbD0bar().size(); ++iScore) {
            outputMl[iScore] = candidate.mlProbD0bar()[iScore];
          }
          registry.fill(HIST("hMassPtImpParPhiY"), invMass, candidate.pt(), candidate.impactParameterXY(), candidate.phi(), yCand, outputMl[0], outputMl[1], outputMl[2]);
        } else {
          registry.fill(HIST("hMassPtImpParPhiY"), invMass, candidate.pt(), candidate.impactParameterXY(), candidate.phi(), yCand);
        }
      }
    }
  }

  // Fill THnSparses for the ML analysis
  /// \param candidate is a particle candidate
  template <Channel channel, bool withMl, typename CCands>
  void fillTree(const CCands& candidate)
  {
    std::vector<float> outputMl = {-999., -999., -999.};
    float invMass{-1.f};
    float yCand{-999.f};
    if constexpr (channel == Channel::DplusToKPiPi) { // D+ -> Kpipi
      invMass = hfHelper.invMassDplusToPiKPi(candidate);
      yCand = hfHelper.yDplus(candidate);
      if constexpr (withMl) {
        for (auto iScore{0u}; iScore < candidate.mlProbDplusToPiKPi().size(); ++iScore) {
          outputMl[iScore] = candidate.mlProbDplusToPiKPi()[iScore];
        }
      }
    } else if constexpr (channel == Channel::DzeroToKPi) {
      if (candidate.isSelD0()) { // D0 -> Kpi
        invMass = hfHelper.invMassD0ToPiK(candidate);
        yCand = hfHelper.yD0(candidate);
        if constexpr (withMl) {
          for (auto iScore{0u}; iScore < candidate.mlProbD0().size(); ++iScore) {
            outputMl[iScore] = candidate.mlProbD0()[iScore];
          }
        }
      }
      if (candidate.isSelD0bar()) {
        invMass = hfHelper.invMassD0barToKPi(candidate);
        yCand = hfHelper.yD0(candidate);
        if constexpr (withMl) {
          for (auto iScore{0u}; iScore < candidate.mlProbD0bar().size(); ++iScore) {
            outputMl[iScore] = candidate.mlProbD0bar()[iScore];
          }
        }
      }
    }
    hfCharmCandLite(
      // Charm candidate meson features
      invMass,
      candidate.pt(),
      yCand,
      candidate.phi(),
      candidate.impactParameterXY(),
      outputMl[0],
      outputMl[1],
      outputMl[2]);
  }

  /// \param candidates are reconstructed candidates
  template <Channel channel, bool withMl, typename CCands>
  void runAnalysis(const CCands& candidates)
  {
    for (auto const& candidate : candidates) {
      fillSparse<channel, withMl>(candidate);
      if (fillLightTreeCandidate) {
        fillTree<channel, withMl>(candidate);
      }
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
