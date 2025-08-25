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
#include <Framework/AnalysisHelpers.h>
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
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(KfptOmegac, kfptOmegac, float);
DECLARE_SOA_COLUMN(KfptPiFromOmegac, kfptPiFromOmegac, float);
DECLARE_SOA_COLUMN(MlProbOmegac, mlProbOmegac, float);
DECLARE_SOA_COLUMN(Cent, cent, float);
} // namespace ml
DECLARE_SOA_TABLE(HfKfOmegacML, "AOD", "HFKFOMEGACML",
                  ml::InvMassCharmBaryon, ml::KfptOmegac, ml::KfptPiFromOmegac, ml::MlProbOmegac, ml::Cent);
} // namespace o2::aod

/// Omegac0 analysis task
struct HfTaskOmegac0ToOmegapi {
  Produces<o2::aod::HfKfOmegacML> kfCandMl;

  Configurable<bool> selectionFlagOmegac0{"selectionFlagOmegac0", true, "Select Omegac0 candidates"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "Max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "Max. cand. rapidity"};
  Configurable<bool> fillTree{"fillTree", false, "Fill tree for local analysis (enabled only with ML)"};

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

  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {100, 0, 1}, "Prompt score"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {700, 2.4, 3.1}, "Cand. inv. mass"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {500, 0, 50}, "Cand. beauty mother pT"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {500, 0, 50}, "Cand. pT"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0, 100}, "Centrality"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin"};
  ConfigurableAxis thnConfigAxisMatchFlag{"thnConfigAxisMatchFlag", {15, -7.5, 7.5}, "Cand. MC match flag"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen pT"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {500, 0, 50}, "Gen beauty mother pT"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    std::array<bool, 8> doprocess{doprocessDataWithKFParticle, doprocessDataWithKFParticleMl, doprocessDataWithKFParticleFT0C, doprocessDataWithKFParticleMlFT0C,
                                  doprocessDataWithKFParticleFT0M, doprocessDataWithKFParticleMlFT0M, doprocessMcWithKFParticle, doprocessMcWithKFParticleMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "One and only one process function should be enabled at a time.");
    }

    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#Omega#pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisMatchFlag{thnConfigAxisMatchFlag, "MC match flag"};
    const AxisSpec thnAxisGenPtD{thnConfigAxisGenPtD, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};
    const AxisSpec thnAxisPromptScore{thnConfigAxisPromptScore, "BDT score prompt"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};

    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisY};

    if (doprocessDataWithKFParticleFT0C || doprocessDataWithKFParticleMlFT0C || doprocessDataWithKFParticleFT0M || doprocessDataWithKFParticleMlFT0M) {
      axes.push_back(thnAxisCent);
      axes.push_back(thnConfigAxisNumPvContr);
    }

    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl) {
      std::vector<AxisSpec> axesMcGen = {thnAxisGenPtD, thnAxisGenPtB, thnAxisY, thnAxisOrigin, thnAxisNumPvContr};
      registry.add("hMcGen", "Thn for generated #Omega_{c}^{0} from charm and beauty", HistType::kTHnSparseD, axesMcGen);
      registry.get<THnSparse>(HIST("hMcGen"))->Sumw2();

      axes.push_back(thnAxisPtB);
      axes.push_back(thnAxisOrigin);
      axes.push_back(thnAxisMatchFlag);
      axes.push_back(thnAxisNumPvContr);
    }

    if (doprocessDataWithKFParticleMl || doprocessDataWithKFParticleMlFT0C || doprocessDataWithKFParticleMlFT0M || doprocessMcWithKFParticleMl) {
      axes.push_back(thnAxisPromptScore);
    }

    registry.add("hReco", "Thn for reco. #Omega_{c}^{0} candidates", HistType::kTHnSparseD, axes);
    registry.get<THnSparse>(HIST("hReco"))->Sumw2();
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
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac(), candidate.mlProbOmegac()[0]);
      } else {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac());
      }
    }
  }

  template <bool applyMl, typename CandType, typename CollType>
  void processDataCent(const CandType& candidates, CollType const& collisions)
  {
    for (const auto& collision : collisions) {

      auto thisCollId = collision.globalIndex();
      auto groupedOmegacCandidates = applyMl ? candidates.sliceBy(candOmegacKFMlPerCollision, thisCollId) : candidates.sliceBy(candOmegacKFPerCollision, thisCollId);
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
          registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac(),
                        cent, numPvContributors, candidate.mlProbOmegac()[0]);
          if (fillTree) {
            kfCandMl(candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfptPiFromOmegac(), candidate.mlProbOmegac()[0], cent);
          }
        } else {
          registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac(),
                        cent, numPvContributors);
        }
      }
    }
  }

  template <bool applyMl, typename CandType, typename CollType>
  void processMc(const CandType& candidates, Omegac0Gen const& mcParticles, TracksMc const&,
                 CollType const& collisions, aod::McCollisions const&)
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
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac(), candidate.ptBhadMotherPart(), candidate.originMcRec(), candidate.flagMcMatchRec(), numPvContributors, candidate.mlProbOmegac()[0]);

      } else {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), candidate.kfRapOmegac(), candidate.ptBhadMotherPart(), candidate.originMcRec(), candidate.flagMcMatchRec(), numPvContributors);
      }
    }

    // MC gen.
    for (const auto& particle : mcParticles) {
      if (yCandGenMax >= 0. && std::abs(particle.rapidityCharmBaryonGen()) > yCandGenMax) {
        continue;
      }

      auto ptGen = particle.pt();
      auto yGen = particle.rapidityCharmBaryonGen();

      int maxNumContrib = 0;
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
      for (const auto& recCol : recoCollsPerMcColl) {
        maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
      }

      if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hMcGen"), ptGen, -1., yGen, RecoDecay::OriginType::Prompt, maxNumContrib);
      } else {
        float ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
        registry.fill(HIST("hMcGen"), ptGen, ptGenB, yGen, RecoDecay::OriginType::NonPrompt, maxNumContrib);
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
