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
/// \brief Task for Ξc^0 → Ξ∓ π± analysis
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University
/// \author Ran Tu <ran.tu@cern.ch>, Fudan University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
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
#include <cstdlib>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Xic0 analysis task

struct HfTaskXic0ToXiPi {
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<bool> fillCent{"fillCent", false, "Flag to fill centrality information"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.8, "max. gen particle rapidity"};
  Configurable<double> yCandRecMax{"yCandRecMax", 0.8, "max. cand. rapidity"};

  SliceCache cache;

  using TracksMc = soa::Join<aod::Tracks, aod::TracksIU, aod::McTrackLabels>;

  using Xic0Cands = soa::Filtered<soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi>>;
  using Xic0CandsKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf>>;
  using Xic0CandsMc = soa::Filtered<soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfXicToXiPiMCRec>>;
  using Xic0CandsMcKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec>>;

  using Xic0CandsMl = soa::Filtered<soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfMlToXiPi>>;
  using Xic0CandsMlKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfMlToXiPi>>;
  using Xic0CandsMlMc = soa::Filtered<soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfMlToXiPi, aod::HfXicToXiPiMCRec>>;
  using Xic0CandsMlMcKF = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfMlToXiPi, aod::HfXicToXiPiMCRec>>;

  using Xic0Gen = soa::Filtered<soa::Join<aod::McParticles, aod::HfXicToXiPiMCGen>>;

  using CollisionsWithEvSels = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterSelectXic0Candidates = aod::hf_sel_toxipi::resultSelections == true;
  Filter filterXicMatchedRec = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi));
  Filter filterXicMatchedGen = nabs(aod::hf_cand_mc_flag::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::XiczeroToXiPi));
  Preslice<Xic0Cands> candXicPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
  Preslice<Xic0CandsKF> candXicKFPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
  Preslice<Xic0CandsMl> candXicMlPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
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
    std::array<bool, 16> doprocess{doprocessDataWithDCAFitter, doprocessDataWithDCAFitterMl, doprocessDataWithDCAFitterFT0C, doprocessDataWithDCAFitterFT0M, doprocessDataWithDCAFitterMlFT0C, doprocessDataWithDCAFitterMlFT0M,
                                   doprocessDataWithKFParticle, doprocessDataWithKFParticleMl, doprocessDataWithKFParticleFT0C, doprocessDataWithKFParticleFT0M, doprocessDataWithKFParticleMlFT0C, doprocessDataWithKFParticleMlFT0M,
                                   doprocessMcWithKFParticle, doprocessMcWithKFParticleMl, doprocessMcWithDCAFitter, doprocessMcWithDCAFitterMl};
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
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality percentile"};

    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl || doprocessMcWithDCAFitter || doprocessMcWithDCAFitterMl) {
      std::vector<AxisSpec> const axesAcc = {thnAxisGenPtD, thnAxisGenPtB, thnAxisY, thnAxisOrigin, thnAxisCent, thnAxisNumPvContr};
      registry.add("hSparseAcc", "Thn for generated Xic0 from charm and beauty", HistType::kTHnSparseD, axesAcc);
      registry.get<THnSparse>(HIST("hSparseAcc"))->Sumw2();

      registry.add("hSparseAccWithRecoColl", "Gen. Xic0 from charm and beauty (associated to a reco collision)", HistType::kTHnSparseD, axesAcc);
      registry.get<THnSparse>(HIST("hSparseAccWithRecoColl"))->Sumw2();

      registry.add("hNumRecoCollPerMcColl", "Number of reco collisions associated to a mc collision;Num. reco. coll. per Mc coll.;", {HistType::kTH1D, {{10, -0.5, 9.5}}});
    }

    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisY};
    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl || doprocessMcWithDCAFitter || doprocessMcWithDCAFitterMl) {
      axes.push_back(thnAxisPtB);
      axes.push_back(thnAxisOrigin);
      axes.push_back(thnAxisMatchFlag);
      axes.push_back(thnAxisCent);
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

  template <bool UseKfParticle, bool UseCentrality, bool ApplyMl, typename CandType, typename CollType>
  void processDataCent(const CandType& candidate, CollType const& collision)
  {
    if (candidate.resultSelections() != true) {
      return;
    }
    double yCharmBaryon;
    if constexpr (UseKfParticle) {
      yCharmBaryon = candidate.kfRapXic();
    } else {
      yCharmBaryon = candidate.y(o2::constants::physics::MassXiC0);
    }
    if (yCandRecMax >= 0. && std::abs(yCharmBaryon) > yCandRecMax) {
      return;
    }

    auto numPvContributors = collision.numContrib();
    float centrality = -999.f;
    if constexpr (UseCentrality) {
      centrality = o2::hf_centrality::getCentralityColl(collision);
    }
    double const ptXic = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
    double const ptPiFromXic = RecoDecay::pt(candidate.pxBachFromCharmBaryon(), candidate.pyBachFromCharmBaryon());
    if constexpr (ApplyMl) {
      registry.fill(HIST("hBdtScoreVsMassVsPtVsYVsCentVsPtPion"),
                    candidate.mlProbToXiPi()[0],
                    candidate.invMassCharmBaryon(),
                    ptXic,
                    yCharmBaryon,
                    centrality,
                    ptPiFromXic,
                    numPvContributors);
    } else {
      registry.fill(HIST("hMassVsPtVsYVsCentVsPtPion"),
                    candidate.invMassCharmBaryon(),
                    ptXic,
                    yCharmBaryon,
                    centrality,
                    ptPiFromXic,
                    numPvContributors);
    }
  }

  template <bool UseKfParticle, bool ApplyMl, typename CandType, typename CollType, typename McCollisionWithCents>
  void processMc(const CandType& candidates,
                 Xic0Gen const& mcParticles,
                 TracksMc const&,
                 CollType const& collisions,
                 McCollisionWithCents const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (candidate.resultSelections() != true) {
        continue;
      }
      double yCharmBaryon;
      if constexpr (UseKfParticle) {
        yCharmBaryon = candidate.kfRapXic();
      } else {
        yCharmBaryon = candidate.y(o2::constants::physics::MassXiC0);
      }
      if (yCandRecMax >= 0. && std::abs(yCharmBaryon) > yCandRecMax) {
        continue;
      }

      auto collision = candidate.template collision_as<CollType>();
      auto numPvContributors = collision.numContrib();
      float const mcCent = o2::hf_centrality::getCentralityColl(collision.template mcCollision_as<McCollisionWithCents>());

      double const ptXic = RecoDecay::pt(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
      if constexpr (ApplyMl) {
        registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsXic0Type"),
                      candidate.mlProbToXiPi()[0],
                      candidate.invMassCharmBaryon(),
                      ptXic,
                      yCharmBaryon,
                      candidate.ptBhadMotherPart(),
                      candidate.originMcRec(),
                      candidate.flagMcMatchRec(),
                      mcCent,
                      numPvContributors);
      } else {
        registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsXic0Type"),
                      candidate.invMassCharmBaryon(),
                      ptXic,
                      yCharmBaryon,
                      candidate.ptBhadMotherPart(),
                      candidate.originMcRec(),
                      candidate.flagMcMatchRec(),
                      mcCent,
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

      auto mcCollision = particle.template mcCollision_as<McCollisionWithCents>();
      unsigned maxNumContrib = 0;
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      for (const auto& recCol : recoCollsPerMcColl) {
        maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
      }
      float const mcCent = o2::hf_centrality::getCentralityColl(mcCollision);

      if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hSparseAcc"),
                      ptGen,
                      -1.,
                      yGen,
                      RecoDecay::OriginType::Prompt,
                      mcCent,
                      maxNumContrib);
      } else {
        float const ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
        registry.fill(HIST("hSparseAcc"),
                      ptGen,
                      ptGenB,
                      yGen,
                      RecoDecay::OriginType::NonPrompt,
                      mcCent,
                      maxNumContrib);
      }

      registry.fill(HIST("hNumRecoCollPerMcColl"), recoCollsPerMcColl.size());

      // fill sparse only for gen particles associated to a reconstructed collision
      if (recoCollsPerMcColl.size() >= 1) {
        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hSparseAccWithRecoColl"),
                        ptGen,
                        -1.,
                        yGen,
                        RecoDecay::OriginType::Prompt,
                        mcCent,
                        maxNumContrib);
        } else {
          float const ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          registry.fill(HIST("hSparseAccWithRecoColl"),
                        ptGen,
                        ptGenB,
                        yGen,
                        RecoDecay::OriginType::NonPrompt,
                        mcCent,
                        maxNumContrib);
        }
      }
    }
  }

  void processDataWithDCAFitter(Xic0Cands const& candidates,
                                CollisionsWithEvSels const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<false, false, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithDCAFitter, "process HfTaskXic0ToXiPi with DCAFitter", true);

  void processDataWithKFParticle(Xic0CandsKF const& candidates,
                                 CollisionsWithEvSels const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicKFPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<true, false, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticle, "process HfTaskXic0ToXiPi with KFParticle", true);

  void processDataWithDCAFitterMl(Xic0CandsMl const& candidates,
                                  CollisionsWithEvSels const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicMlPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<false, false, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithDCAFitterMl, "process HfTaskXic0ToXiPi with DCAFitter and ML selections", false);

  void processDataWithKFParticleMl(Xic0CandsMlKF const& candidates,
                                   CollisionsWithEvSels const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicKFMlPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<true, false, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleMl, "process HfTaskXic0ToXiPi with KFParticle and ML selections", false);

  void processDataWithDCAFitterFT0C(Xic0Cands const& candidates,
                                    CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<false, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithDCAFitterFT0C, "process HfTaskXic0ToXiPi with DCAFitter and with FT0C centrality", false);

  void processDataWithKFParticleFT0C(Xic0CandsKF const& candidates,
                                     CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicKFPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<true, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleFT0C, "process HfTaskXic0ToXiPi with KFParticle and with FT0C centrality", false);

  void processDataWithDCAFitterFT0M(Xic0Cands const& candidates,
                                    CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<false, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithDCAFitterFT0M, "process HfTaskXic0ToXiPi with DCAFitter and with FT0M centrality", false);

  void processDataWithKFParticleFT0M(Xic0CandsKF const& candidates,
                                     CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicKFPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<true, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleFT0M, "process HfTaskXic0ToXiPi with KFParticle and with FT0M centrality", false);

  void processDataWithDCAFitterMlFT0C(Xic0CandsMl const& candidates,
                                      CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicMlPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<false, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithDCAFitterMlFT0C, "process HfTaskXic0ToXiPi with DCAFitter and ML selections and with FT0C centrality", false);

  void processDataWithKFParticleMlFT0C(Xic0CandsMlKF const& candidates,
                                       CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicKFMlPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<true, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleMlFT0C, "process HfTaskXic0ToXiPi with KFParticle and ML selections and with FT0C centrality", false);

  void processDataWithDCAFitterMlFT0M(Xic0CandsMl const& candidates,
                                      CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicMlPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<false, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithDCAFitterMlFT0M, "process HfTaskXic0ToXiPi with DCAFitter and ML selections and with FT0M centrality", false);

  void processDataWithKFParticleMlFT0M(Xic0CandsMlKF const& candidates,
                                       CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedXicCandidates = candidates.sliceBy(candXicKFMlPerCollision, thisCollId);
      for (const auto& candidate : groupedXicCandidates) {
        processDataCent<true, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processDataWithKFParticleMlFT0M, "process HfTaskXic0ToXiPi with KFParticle and ML selections and with FT0M centrality", false);

  void processMcWithDCAFitter(Xic0CandsMc const& xic0CandidatesMc,
                              Xic0Gen const& mcParticles,
                              TracksMc const& tracks,
                              CollisionsWithMcLabels const& collisions,
                              McCollisionsCentFT0Ms const& mcCollisions)
  {
    processMc<false, false>(xic0CandidatesMc, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processMcWithDCAFitter, "Process MC with KFParticle", false);

  void processMcWithKFParticle(Xic0CandsMcKF const& xic0CandidatesMcKf,
                               Xic0Gen const& mcParticles,
                               TracksMc const& tracks,
                               CollisionsWithMcLabels const& collisions,
                               McCollisionsCentFT0Ms const& mcCollisions)
  {
    processMc<true, false>(xic0CandidatesMcKf, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processMcWithKFParticle, "Process MC with KFParticle", false);

  void processMcWithDCAFitterMl(Xic0CandsMlMc const& xic0CandidatesMlMc,
                                Xic0Gen const& mcParticles,
                                TracksMc const& tracks,
                                CollisionsWithMcLabels const& collisions,
                                McCollisionsCentFT0Ms const& mcCollisions)
  {
    processMc<false, true>(xic0CandidatesMlMc, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processMcWithDCAFitterMl, "Process MC with KFParticle and ML selections", false);

  void processMcWithKFParticleMl(Xic0CandsMlMcKF const& xic0CandidatesMlMcKf,
                                 Xic0Gen const& mcParticles,
                                 TracksMc const& tracks,
                                 CollisionsWithMcLabels const& collisions,
                                 McCollisionsCentFT0Ms const& mcCollisions)
  {
    processMc<true, true>(xic0CandidatesMlMcKf, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskXic0ToXiPi, processMcWithKFParticleMl, "Process MC with KFParticle and ML selections", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskXic0ToXiPi>(cfgc)};
}
