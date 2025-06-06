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

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Omegac0 analysis task

struct HfTaskOmegac0ToOmegapi {
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
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
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Filter filterOmegaCToOmegaPiFlag = (aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi))) != static_cast<uint8_t>(0);
  Filter filterOmegaCMatchedRec = nabs(aod::hf_cand_xic0_omegac0::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi));
  Filter filterOmegaCMatchedGen = nabs(aod::hf_cand_xic0_omegac0::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi));

  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  // ThnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {50, 0, 1}, "Prompt score bins"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 2.4, 3.1}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {1000, 0, 100}, "Cand. beauty mother pTB bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0, 20}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisMatchFlag{"thnConfigAxisMatchFlag", {15, -7.5, 7.5}, "Cand. MC Match Flag type"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen Pt D"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    std::array<bool, 4> doprocess{doprocessDataWithKFParticle, doprocessMcWithKFParticle, doprocessDataWithKFParticleMl, doprocessMcWithKFParticleMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "One and only one process function should be enabled at a time.");
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
