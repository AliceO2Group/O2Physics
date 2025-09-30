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

/// \file taskXic0ToXiPi.cxx
/// \brief Task for Ξc^0 → Ξ∓ π± Kf analysis
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University
/// \author Ran Tu <ran.tu@cern.ch>, Fudan University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Xic0 analysis task

struct HfTaskXic0ToXiPi {
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<bool> fillCent{"fillCent", false, "Flag to fill centrality information"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.8, "max. gen particle rapidity"};
  Configurable<double> yCandRecMax{"yCandRecMax", 0.8, "max. cand. rapidity"};

  HfHelper hfHelper;
  SliceCache cache;

  using TracksMc = soa::Join<aod::Tracks, aod::TracksIU, aod::McTrackLabels>;

  using Xic0CandsKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf>>;
  using Xic0CandsMcKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec>>;

  using Xic0CandsMlKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfMlToXiPiKf>>;
  using Xic0CandsMlMcKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfMlToXiPiKf, aod::HfXicToXiPiMCRec>>;

  using Xic0Gen = soa::Filtered<soa::Join<aod::McParticles, aod::HfXicToXiPiMCGen>>;

  using CollisionsWithEvSels = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Filter filterSelectXic0Candidates = aod::hf_sel_toxipi::resultSelections == true;
  Filter filterXicMatchedRec = nabs(aod::hf_cand_xic0_omegac0::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi));
  Filter filterXicMatchedGen = nabs(aod::hf_cand_xic0_omegac0::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi));
  Preslice<Xic0CandsKF> candXicKFPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
  Preslice<Xic0CandsMlKF> candXicKFMlPerCollision = aod::hf_cand_xic0_omegac0::collisionId;

  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  // ThnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {100, 0, 1}, "Prompt score bins"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 2.4, 3.1}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {1000, 0, 100}, "Cand. beauty mother pTB bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0, 20}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0, 100}, "Centrality bins"};
  ConfigurableAxis thnConfigAxisPtPion{"thnConfigAxisPtPion", {100, 0, 10}, "PtPion from Xic0 bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisMatchFlag{"thnConfigAxisMatchFlag", {15, -7.5, 7.5}, "Cand. MC Match Flag type"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen Pt D"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    std::array<bool, 8> doprocess{doprocessDataWithKFParticle, doprocessMcWithKFParticle, doprocessDataWithKFParticleMl, doprocessMcWithKFParticleMl, doprocessDataWithKFParticleFT0C, doprocessDataWithKFParticleMlFT0C, doprocessDataWithKFParticleFT0M, doprocessDataWithKFParticleMlFT0M};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "One and only one process function should be enabled at a time.");
    }

    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#Xi#pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisMatchFlag{thnConfigAxisMatchFlag, "MatchFlag"};
    const AxisSpec thnAxisGenPtD{thnConfigAxisGenPtD, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};

    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl) {
      std::vector<AxisSpec> const axesAcc = {thnAxisGenPtD, thnAxisGenPtB, thnAxisY, thnAxisOrigin, thnAxisNumPvContr};
      registry.add("hSparseAcc", "Thn for generated Xic0 from charm and beauty", HistType::kTHnSparseD, axesAcc);
      registry.get<THnSparse>(HIST("hSparseAcc"))->Sumw2();
    }

    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisY};
    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl) {
      axes.push_back(thnAxisPtB);
      axes.push_back(thnAxisOrigin);
      axes.push_back(thnAxisMatchFlag);
      axes.push_back(thnAxisNumPvContr);
    }
    if (applyMl) {
      const AxisSpec thnAxisPromptScore{thnConfigAxisPromptScore, "BDT score prompt."};
      axes.insert(axes.begin(), thnAxisPromptScore);
      registry.add("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsXic0Type", "Thn for Xic0 candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsXic0Type"))->Sumw2();
    } else {
      registry.add("hMassVsPtVsPtBVsYVsOriginVsXic0Type", "Thn for Xic0 candidates", HistType::kTHnSparseF, axes);
      registry.get<THnSparse>(HIST("hMassVsPtVsPtBVsYVsOriginVsXic0Type"))->Sumw2();
    }
    if (fillCent) {
      const AxisSpec thnAxisPromptScore{thnConfigAxisPromptScore, "BDT score prompt."};
      const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality."};
      const AxisSpec thnAxisPtPion{thnConfigAxisPtPion, "Pt of Pion from Xic0."};
      std::vector<AxisSpec> const axesWithBdtCent = {thnAxisPromptScore, thnAxisMass, thnAxisPt, thnAxisY, thnAxisCent, thnAxisPtPion, thnConfigAxisNumPvContr};
      std::vector<AxisSpec> const axesWithCent = {thnAxisMass, thnAxisPt, thnAxisY, thnAxisCent, thnAxisPtPion, thnConfigAxisNumPvContr};
      registry.add("hBdtScoreVsMassVsPtVsYVsCentVsPtPion", "Thn for Xic0 candidates with BDT&Cent&pTpi", HistType::kTHnSparseD, axesWithBdtCent);
      registry.add("hMassVsPtVsYVsCentVsPtPion", "Thn for Xic0 candidates with Cent&pTpi", HistType::kTHnSparseD, axesWithCent);
      registry.get<THnSparse>(HIST("hBdtScoreVsMassVsPtVsYVsCentVsPtPion"))->Sumw2();
      registry.get<THnSparse>(HIST("hMassVsPtVsYVsCentVsPtPion"))->Sumw2();
    }
  }

  template <bool ApplyMl, typename CandType, typename CollType>
  void processData(const CandType& candidates, CollType const&)
  {
    for (const auto& candidate : candidates) {
      if (candidate.resultSelections() != true) {
        continue;
      }
      if (yCandRecMax >= 0. && std::abs(candidate.kfRapXic()) > yCandRecMax) {
        continue;
      }

      if constexpr (ApplyMl) {
        registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsXic0Type"), candidate.mlProbToXiPi()[0], candidate.invMassCharmBaryon(), candidate.kfptXic(), candidate.kfRapXic());
      } else {
        registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsXic0Type"), candidate.invMassCharmBaryon(), candidate.kfptXic(), candidate.kfRapXic());
      }
    }
  }

  template <bool UseCentrality, bool ApplyMl, typename CandType, typename CollType>
  void processDataCent(const CandType& candidates, CollType const& collisions)
  {
    for (const auto& collision : collisions) {

      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = ApplyMl
                                    ? candidates.sliceBy(candXicKFMlPerCollision, thisCollId)
                                    : candidates.sliceBy(candXicKFPerCollision, thisCollId);
      // auto numPvContributors = collision.numContrib();

      for (const auto& candidate : groupedXicCandidates) {
        if (candidate.resultSelections() != true) {
          continue;
        }
        if (yCandRecMax >= 0. && std::abs(candidate.kfRapXic()) > yCandRecMax) {
          continue;
        }

        auto numPvContributors = candidate.template collision_as<CollType>().numContrib();
        float centrality = -999.f;
        if constexpr (UseCentrality) {
          auto const& collision = candidate.template collision_as<CollType>();
          centrality = o2::hf_centrality::getCentralityColl(collision);
        }
        double const kfptXic = RecoDecay::sqrtSumOfSquares(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
        double const kfptPiFromXic = RecoDecay::sqrtSumOfSquares(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());
        if constexpr (ApplyMl) {
          registry.fill(HIST("hBdtScoreVsMassVsPtVsYVsCentVsPtPion"),
                        candidate.mlProbToXiPi()[0],
                        candidate.invMassCharmBaryon(),
                        kfptXic,
                        candidate.kfRapXic(),
                        centrality,
                        kfptPiFromXic,
                        numPvContributors);
        } else {
          registry.fill(HIST("hMassVsPtVsYVsCentVsPtPion"),
                        candidate.invMassCharmBaryon(),
                        kfptXic,
                        candidate.kfRapXic(),
                        centrality,
                        kfptPiFromXic,
                        numPvContributors);
        }
      }
    }
  }

  template <bool ApplyMl, typename CandType, typename CollType>
  void processMc(const CandType& candidates,
                 Xic0Gen const& mcParticles,
                 TracksMc const&,
                 CollType const& collisions,
                 aod::McCollisions const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (candidate.resultSelections() != true) {
        continue;
      }
      if (yCandRecMax >= 0. && std::abs(candidate.kfRapXic()) > yCandRecMax) {
        continue;
      }

      auto numPvContributors = candidate.template collision_as<CollType>().numContrib();
      double const kfptXic = RecoDecay::sqrtSumOfSquares(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
      if constexpr (ApplyMl) {
        registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsXic0Type"),
                      candidate.mlProbToXiPi()[0],
                      candidate.invMassCharmBaryon(),
                      kfptXic,
                      candidate.kfRapXic(),
                      candidate.ptBhadMotherPart(),
                      candidate.originMcRec(),
                      candidate.flagMcMatchRec(),
                      numPvContributors);
      } else {
        registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsXic0Type"),
                      candidate.invMassCharmBaryon(),
                      kfptXic,
                      candidate.kfRapXic(),
                      candidate.ptBhadMotherPart(),
                      candidate.originMcRec(),
                      candidate.flagMcMatchRec(),
                      numPvContributors);
      }
    }

    // MC gen.
    for (const auto& particle : mcParticles) {
      if (yCandGenMax >= 0. && std::abs(particle.rapidityCharmBaryonGen()) > yCandGenMax) {
        continue;
      }

      auto ptGen = particle.pt();
      auto yGen = particle.rapidityCharmBaryonGen();

      unsigned maxNumContrib = 0;
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
      for (const auto& recCol : recoCollsPerMcColl) {
        maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
      }

      if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hSparseAcc"),
                      ptGen,
                      -1.,
                      yGen,
                      RecoDecay::OriginType::Prompt,
                      maxNumContrib);
      } else {
        float const ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
        registry.fill(HIST("hSparseAcc"),
                      ptGen,
                      ptGenB,
                      yGen,
                      RecoDecay::OriginType::NonPrompt,
                      maxNumContrib);
      }
    }
  }

  void processDataWithKFParticle(Xic0CandsKF const& candidates,
                                 CollisionsWithEvSels const& collisions)
  {
    processDataCent<false, false>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticle, "process HfTaskXic0ToXiPi  with KFParticle", true);

  void processDataWithKFParticleMl(Xic0CandsMlKF const& candidates,
                                   CollisionsWithEvSels const& collisions)
  {
    processDataCent<false, true>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleMl, "process HfTaskXic0ToXiPi  with KFParticle and ML selections", false);

  void processDataWithKFParticleFT0C(Xic0CandsKF const& candidates,
                                     CollisionsWithFT0C const& collisions)
  {
    processDataCent<true, false>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleFT0C, "process HfTaskXic0ToXiPi  with KFParticle and with FT0C centrality", false);

  void processDataWithKFParticleFT0M(Xic0CandsKF const& candidates,
                                     CollisionsWithFT0M const& collisions)
  {
    processDataCent<true, false>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleFT0M, "process HfTaskXic0ToXiPi  with KFParticle and with FT0M centrality", false);

  void processDataWithKFParticleMlFT0C(Xic0CandsMlKF const& candidates,
                                       CollisionsWithFT0C const& collisions)
  {
    processDataCent<true, true>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleMlFT0C, "process HfTaskXic0ToXiPi  with KFParticle and ML selections and with FT0C centrality", false);

  void processDataWithKFParticleMlFT0M(Xic0CandsMlKF const& candidates,
                                       CollisionsWithFT0M const& collisions)
  {
    processDataCent<true, true>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleMlFT0M, "process HfTaskXic0ToXiPi  with KFParticle and ML selections and with FT0M centrality", false);

  void processMcWithKFParticle(Xic0CandsMcKF const& xic0CandidatesMcKf,
                               Xic0Gen const& mcParticles,
                               TracksMc const& tracks,
                               CollisionsWithMcLabels const& collisions,
                               aod::McCollisions const& mcCollisions)
  {
    processMc<false>(xic0CandidatesMcKf, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processMcWithKFParticle, "Process MC with KFParticle", false);

  void processMcWithKFParticleMl(Xic0CandsMlMcKF const& xic0CandidatesMlMcKf,
                                 Xic0Gen const& mcParticles,
                                 TracksMc const& tracks,
                                 CollisionsWithMcLabels const& collisions,
                                 aod::McCollisions const& mcCollisions)
  {
    processMc<true>(xic0CandidatesMlMcKf, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processMcWithKFParticleMl, "Process MC with KFParticle and ML selections", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskXic0ToXiPi>(cfgc)};
}
