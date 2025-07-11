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

/// \file taskOmegac0ToOmegapi.cxx
/// \brief OmegaC0 analysis task
/// \author Yunfan Liu <yunfan.liu@cern.ch>, China University of Geosciences
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston

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
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
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

namespace o2::aod
{
namespace ml
{
// collision info
DECLARE_SOA_COLUMN(KfptPiFromOmegac, kfptPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfptOmegac, kfptOmegac, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(MlProbOmegac, mlProbOmegac, float);
DECLARE_SOA_COLUMN(Cent, cent, float);
} // namespace ml
DECLARE_SOA_TABLE(HfKfOmegacML, "AOD", "HFKFOMEGACML",
                  ml::InvMassCharmBaryon, ml::KfptOmegac, ml::KfptPiFromOmegac, ml::MlProbOmegac, ml::Cent);
} // namespace o2::aod

/// Omegac0 analysis task

struct HfTaskOmegac0ToOmegapi {

  Produces<o2::aod::HfKfOmegacML> kfCandMl;
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<bool> fillCent{"fillCent", false, "Flag to fill centrality information"};
  Configurable<bool> fillTree{"fillTree", false, "Fill TTree for local analysis.(Enabled only with ML)"};
  Configurable<bool> selectionFlagOmegac0{"selectionFlagOmegac0", true, "Select Omegac0 candidates"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};

  HfHelper hfHelper;
  SliceCache cache;

  using TracksMc = soa::Join<aod::Tracks, aod::TracksIU, aod::McTrackLabels>;

  using Omegac0Cands = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi>>;
  using Omegac0CandsKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfOmegacKf>>;
  using OmegaC0CandsMcKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfOmegacKf, aod::HfToOmegaPiMCRec>>;

  using Omegac0CandsMl = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfMlSelOmegacToOmegaPi>>;
  using Omegac0CandsMlKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfMlSelOmegacToOmegaPi, aod::HfOmegacKf>>;
  using Omegac0CandsMlMcKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfMlSelOmegacToOmegaPi, aod::HfOmegacKf, aod::HfToOmegaPiMCRec>>;

  using Omegac0Gen = soa::Filtered<soa::Join<aod::McParticles, aod::HfToOmegaPiMCGen>>;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Filter filterOmegaCToOmegaPiFlag = (aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi))) != static_cast<uint8_t>(0);
  Filter filterOmegaCMatchedRec = nabs(aod::hf_cand_xic0_omegac0::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi));
  Filter filterOmegaCMatchedGen = nabs(aod::hf_cand_xic0_omegac0::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi));
  Preslice<Omegac0CandsKF> candOmegacKFPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
  Preslice<Omegac0CandsMlKF> candOmegacKFMlPerCollision = aod::hf_cand_xic0_omegac0::collisionId;

  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  // ThnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {100, 0, 1}, "Prompt score bins"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 2.4, 3.1}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {1000, 0, 100}, "Cand. beauty mother pTB bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0, 20}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0, 100}, "Centrality bins"};
  ConfigurableAxis thnConfigAxisPtPion{"thnConfigAxisPtPion", {100, 0, 10}, "PtPion from Omegac0 bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisMatchFlag{"thnConfigAxisMatchFlag", {15, -7.5, 7.5}, "Cand. MC Match Flag type"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen Pt D"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    std::array<bool, 6> doprocess{doprocessDataWithKFParticle, doprocessDataWithKFParticleMl, doprocessDataWithKFParticleFT0C, doprocessDataWithKFParticleMlFT0C, doprocessDataWithKFParticleFT0M, doprocessDataWithKFParticleMlFT0M};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) > 1) {
      LOGP(fatal, "At most one data process function should be enabled at a time.");
    }

    std::array<bool, 2> doprocessMc{doprocessMcWithKFParticle, doprocessMcWithKFParticleMl};
    if (std::accumulate(doprocessMc.begin(), doprocessMc.end(), 0) > 1) {
      LOGP(fatal, "At most one MC process function should be enabled at a time.");
    }

    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0) + std::accumulate(doprocessMc.begin(), doprocessMc.end(), 0)) == 0) {
      LOGP(fatal, "At least one process function should be enabled.");
    }

    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#Omega#pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisMatchFlag{thnConfigAxisMatchFlag, "MatchFlag"};
    const AxisSpec thnAxisGenPtD{thnConfigAxisGenPtD, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};

    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl) {
      std::vector<AxisSpec> axesAcc = {thnAxisGenPtD, thnAxisGenPtB, thnAxisY, thnAxisOrigin, thnAxisNumPvContr};
      registry.add("hSparseAcc", "Thn for generated Omega0 from charm and beauty", HistType::kTHnSparseD, axesAcc);
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
      registry.add("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsOmegac0Type", "Thn for Omegac0 candidates", HistType::kTHnSparseD, axes);
      registry.get<THnSparse>(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsOmegac0Type"))->Sumw2();
    } else {
      registry.add("hMassVsPtVsPtBVsYVsOriginVsOmegac0Type", "Thn for Omegac0 candidates", HistType::kTHnSparseF, axes);
      registry.get<THnSparse>(HIST("hMassVsPtVsPtBVsYVsOriginVsOmegac0Type"))->Sumw2();
    }
    if (fillCent) {
      const AxisSpec thnAxisPromptScore{thnConfigAxisPromptScore, "BDT score prompt."};
      const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality."};
      const AxisSpec thnAxisPtPion{thnConfigAxisPtPion, "Pt of Pion from Omegac0."};
      std::vector<AxisSpec> axesWithBdtCent = {thnAxisPromptScore, thnAxisMass, thnAxisPt, thnAxisY, thnAxisCent, thnAxisPtPion, thnConfigAxisNumPvContr};
      std::vector<AxisSpec> axesWithCent = {thnAxisMass, thnAxisPt, thnAxisY, thnAxisCent, thnAxisPtPion, thnConfigAxisNumPvContr};
      registry.add("hBdtScoreVsMassVsPtVsYVsCentVsPtPion", "Thn for Omegac0 candidates with BDT&Cent&pTpi", HistType::kTHnSparseD, axesWithBdtCent);
      registry.add("hMassVsPtVsYVsCentVsPtPion", "Thn for Omegac0 candidates with Cent&pTpi", HistType::kTHnSparseD, axesWithCent);
      registry.get<THnSparse>(HIST("hBdtScoreVsMassVsPtVsYVsCentVsPtPion"))->Sumw2();
      registry.get<THnSparse>(HIST("hMassVsPtVsYVsCentVsPtPion"))->Sumw2();
    }
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    return o2::hf_centrality::getCentralityColl<Coll>(collision);
  }

  template <bool applyMl, typename CandType, typename CollType>
  void processData(const CandType& candidates, CollType const&)
  {
    for (const auto& candidate : candidates) {
      if (!(candidate.resultSelections() == true || (candidate.resultSelections() == false && !selectionFlagOmegac0))) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(candidate.kfRapOmegac()) > yCandRecoMax) {
        continue;
      }

      if constexpr (applyMl) {
        registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), candidate.mlProbOmegac()[0], candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac());
      } else {
        registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac());
      }
    }
  }

  template <bool applyMl, typename CandType, typename CollType>
  void processDataCent(const CandType& candidates, CollType const& collisions)
  {
    for (const auto& collision : collisions) {

      auto thisCollId = collision.globalIndex();
      auto groupedOmegacCandidates = applyMl
                                       ? candidates.sliceBy(candOmegacKFMlPerCollision, thisCollId)
                                       : candidates.sliceBy(candOmegacKFPerCollision, thisCollId);
      auto numPvContributors = collision.numContrib();

      for (const auto& candidate : groupedOmegacCandidates) {
        if (!(candidate.resultSelections() == true || (candidate.resultSelections() == false && !selectionFlagOmegac0))) {
          continue;
        }
        if (yCandRecoMax >= 0. && std::abs(candidate.kfRapOmegac()) > yCandRecoMax) {
          continue;
        }
        float cent = evaluateCentralityColl(collision);
        if constexpr (applyMl) {
          if (fillTree) {
            kfCandMl(candidate.invMassCharmBaryon(),
                     candidate.ptCharmBaryon(),
                     candidate.kfptPiFromOmegac(),
                     candidate.mlProbOmegac()[0],
                     cent);
          } else {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsYVsCentVsPtPion"),
                          candidate.mlProbOmegac()[0],
                          candidate.invMassCharmBaryon(),
                          candidate.ptCharmBaryon(),
                          candidate.kfRapOmegac(),
                          cent,
                          candidate.kfptPiFromOmegac(),
                          numPvContributors);
          }
        } else {
          registry.fill(HIST("hMassVsPtVsYVsCentVsPtPion"),
                        candidate.invMassCharmBaryon(),
                        candidate.ptCharmBaryon(),
                        candidate.kfRapOmegac(),
                        cent,
                        candidate.kfptPiFromOmegac(),
                        numPvContributors);
        }
      }
    }
  }

  template <bool applyMl, typename CandType, typename CollType>
  void processMc(const CandType& candidates,
                 Omegac0Gen const& mcParticles,
                 TracksMc const&,
                 CollType const& collisions,
                 aod::McCollisions const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (!(candidate.resultSelections() == true || (candidate.resultSelections() == false && !selectionFlagOmegac0))) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(candidate.kfRapOmegac()) > yCandRecoMax) {
        continue;
      }

      auto numPvContributors = candidate.template collision_as<CollType>().numContrib();

      if constexpr (applyMl) {
        registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), candidate.mlProbOmegac()[0], candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac(), candidate.ptBhadMotherPart(), candidate.originRec(), candidate.flagMcMatchRec(), numPvContributors);

      } else {
        registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac(), candidate.ptBhadMotherPart(), candidate.originRec(), candidate.flagMcMatchRec(), numPvContributors);
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

      if (particle.originGen() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hSparseAcc"), ptGen, -1., yGen, RecoDecay::OriginType::Prompt, maxNumContrib);
      } else {
        float ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
        registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, RecoDecay::OriginType::NonPrompt, maxNumContrib);
      }
    }
  }

  void processDataWithKFParticle(Omegac0CandsKF const& candidates,
                                 Collisions const& collisions)
  {
    processData<false>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticle, "process HfTaskOmegac0ToOmegapi with KFParticle", false);

  void processDataWithKFParticleMl(Omegac0CandsMlKF const& candidates,
                                   Collisions const& collisions)
  {
    processData<true>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticleMl, "process HfTaskOmegac0ToOmegapi with KFParticle and ML selections", false);

  void processDataWithKFParticleFT0C(Omegac0CandsKF const& candidates,
                                     CollisionsWithFT0C const& collisions)
  {
    processDataCent<false>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticleFT0C, "process HfTaskOmegac0ToOmegapi with KFParticle and with FT0C centrality", false);

  void processDataWithKFParticleMlFT0C(Omegac0CandsMlKF const& candidates,
                                       CollisionsWithFT0C const& collisions)
  {
    processDataCent<true>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticleMlFT0C, "process HfTaskOmegac0ToOmegapi with KFParticle and ML selections and with FT0C centrality", false);

  void processDataWithKFParticleFT0M(Omegac0CandsKF const& candidates,
                                     CollisionsWithFT0M const& collisions)
  {
    processDataCent<false>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticleFT0M, "process HfTaskOmegac0ToOmegapi with KFParticle and with FT0M centrality", false);

  void processDataWithKFParticleMlFT0M(Omegac0CandsMlKF const& candidates,
                                       CollisionsWithFT0M const& collisions)
  {
    processDataCent<true>(candidates, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticleMlFT0M, "process HfTaskOmegac0ToOmegapi with KFParticle and ML selections and with FT0M centrality", false);

  void processMcWithKFParticle(OmegaC0CandsMcKF const& omegaC0CandidatesMcKF,
                               Omegac0Gen const& mcParticles,
                               TracksMc const& tracks,
                               CollisionsWithMcLabels const& collisions,
                               aod::McCollisions const& mcCollisions)
  {
    processMc<false>(omegaC0CandidatesMcKF, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processMcWithKFParticle, "Process MC with KFParticle", false);

  void processMcWithKFParticleMl(Omegac0CandsMlMcKF const& omegac0CandidatesMlMcKF,
                                 Omegac0Gen const& mcParticles,
                                 TracksMc const& tracks,
                                 CollisionsWithMcLabels const& collisions,
                                 aod::McCollisions const& mcCollisions)
  {
    processMc<true>(omegac0CandidatesMlMcKF, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processMcWithKFParticleMl, "Process MC with KFParticle and ML selections", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskOmegac0ToOmegapi>(cfgc)};
}
