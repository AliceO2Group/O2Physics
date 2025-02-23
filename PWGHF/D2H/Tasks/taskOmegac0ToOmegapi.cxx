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
  Configurable<bool> selectionFlagOmegac0{"selectionFlagOmegac0", false, "Selection Flag for Omegac0 candidates"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};

  HfHelper hfHelper;
  SliceCache cache;
  using MyTracksWMc = soa::Join<aod::Tracks, aod::TracksIU, aod::McTrackLabels>;

  using Omegac0Candidates = soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi>;
  using Omegac0CandidatesKF = soa::Join<Omegac0Candidates, aod::HfOmegacKf>;
  using OmegaC0CandidatesMcKF = soa::Join<Omegac0CandidatesKF, aod::HfToOmegaPiMCRec>;

  using Omegac0CandidatesMl = soa::Join<Omegac0Candidates, aod::HfMlSelOmegacToOmegaPi>;
  using Omegac0CandidatesMlKF = soa::Join<Omegac0CandidatesMl, aod::HfOmegacKf>;
  using Omegac0CandidatesMlMcKF = soa::Join<Omegac0CandidatesMlKF, aod::HfToOmegaPiMCRec>;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  Partition<Omegac0CandidatesKF> selectedOmegac0CandidatesKF = aod::hf_sel_toomegapi::resultSelections && !selectionFlagOmegac0;
  Partition<Omegac0CandidatesMlKF> selectedOmegac0CandidatesMlKF = aod::hf_sel_toomegapi::resultSelections && !selectionFlagOmegac0;

  // ThnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {50, 0, 1}, "Prompt score bins"};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {120, 2.4, 3.1}, "Cand. inv-mass bins"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {1000, 0, 100}, "Cand. beauty mother pTB bins"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0, 20}, "Cand. pT bins"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin type"};
  ConfigurableAxis thnConfigAxisGenPtD{"thnConfigAxisGenPtD", {500, 0, 50}, "Gen Pt D"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    std::array<bool, 12> doprocess{doprocessDataWithKFParticle, doprocessMcWithKFParticle, doprocessDataWithKFParticleMl, doprocessMcWithKFParticleMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "One and only one process function should be enabled at a time.");
    }

    const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (#Omega #pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisGenPtD{thnConfigAxisGenPtD, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisGenPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of PV contributors"};

    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl) {
      std::vector<AxisSpec> axesAcc = {thnAxisGenPtD, thnAxisGenPtB, thnAxisY, thnAxisOrigin, thnAxisNumPvContr};
      registry.add("hSparseAcc", "Thn for generated Omega0 from charm and beauty", HistType::kTHnSparseD, axesAcc);
      registry.get<THnSparse>(HIST("hSparseAcc"))->Sumw2();
    }

    std::vector<AxisSpec> axes = {
      thnAxisMass,
      thnAxisPt,
      thnAxisY,
    };
    if (doprocessMcWithKFParticle || doprocessMcWithKFParticleMl) {
      axes.push_back(thnAxisPtB);
      axes.push_back(thnAxisOrigin);
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
      if (!(candidate.hfflag() & 1 << aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(candidate.kfRapOmegac()) > yCandRecoMax) {
        continue;
      }
      float massOmegac0;
      massOmegac0 = candidate.invMassCharmBaryon();
      auto rapidityCandidate = candidate.kfRapOmegac();
      auto ptCandidate = candidate.ptCharmBaryon();
      if constexpr (applyMl) {
        registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), candidate.mlProbOmegac()[0], massOmegac0, ptCandidate, rapidityCandidate);
      } else {
        registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), massOmegac0, ptCandidate, rapidityCandidate);
      }
    }
  }

  void processDataWithKFParticle(Omegac0CandidatesKF const&, Collisions const& collisions)
  {
    processData<false>(selectedOmegac0CandidatesKF, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticle, "process HfTaskOmegac0ToOmegapi with KFParticle", false);
  // TODO: add processKFParticle

  void processDataWithKFParticleMl(Omegac0CandidatesMlKF const&, Collisions const& collisions)
  {
    processData<true>(selectedOmegac0CandidatesMlKF, collisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processDataWithKFParticleMl, "process HfTaskOmegac0ToOmegapi with KFParticle and ML selections", false);
  // TODO: add processKFParticleMl

  template <bool applyMl, typename CandType, typename CollType>
  void processMc(const CandType& candidates,
                 soa::Join<aod::McParticles, aod::HfToOmegaPiMCGen> const& mcParticles,
                 MyTracksWMc const&,
                 CollType const& collisions,
                 aod::McCollisions const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(candidate.kfRapOmegac()) > yCandRecoMax) {
        continue;
      }
      auto collision = candidate.template collision_as<CollType>();
      auto numPvContributors = collision.numContrib();
      float massOmegac0;
      massOmegac0 = candidate.invMassCharmBaryon();
      auto ptCandidate = candidate.ptCharmBaryon();
      auto rapidityCandidate = candidate.kfRapOmegac();
      if (candidate.resultSelections() && !selectionFlagOmegac0)
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi)) {
          if constexpr (applyMl) {
            registry.fill(HIST("hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), candidate.mlProbOmegac()[0], massOmegac0, ptCandidate, rapidityCandidate, candidate.ptBhadMotherPart(), candidate.originRec(), numPvContributors);

          } else {
            registry.fill(HIST("hMassVsPtVsPtBVsYVsOriginVsOmegac0Type"), massOmegac0, ptCandidate, rapidityCandidate, candidate.ptBhadMotherPart(), candidate.originRec(), numPvContributors);
          }
        }
    }
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi) {
        if (yCandGenMax >= 0. && std::abs(particle.rapidityCharmBaryonGen()) > yCandGenMax) {
          continue;
        }
        float ptGenB = -1;
        auto ptGen = particle.pt();
        auto yGen = particle.rapidityCharmBaryonGen();

        unsigned maxNumContrib = 0;
        const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
        for (const auto& recCol : recoCollsPerMcColl) {
          maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
        }

        if (particle.originGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 1, maxNumContrib);

        } else {
          ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          registry.fill(HIST("hSparseAcc"), ptGen, ptGenB, yGen, 2, maxNumContrib);
        }
      }
    }
  }

  void processMcWithKFParticle(OmegaC0CandidatesMcKF const& omegaC0CandidatesMcKF,
                               soa::Join<aod::McParticles, aod::HfToOmegaPiMCGen> const& mcParticles,
                               MyTracksWMc const& tracks,
                               CollisionsWithMcLabels const& collisions,
                               aod::McCollisions const& mcCollisions)
  {
    processMc<false>(omegaC0CandidatesMcKF, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processMcWithKFParticle, "Process MC with KFParticle", false);
  // TODO: add the processMcWithKFParticle

  void processMcWithKFParticleMl(Omegac0CandidatesMlMcKF const& omegac0CandidatesMlMcKF,
                                 soa::Join<aod::McParticles, aod::HfToOmegaPiMCGen> const& mcParticles,
                                 MyTracksWMc const& tracks,
                                 CollisionsWithMcLabels const& collisions,
                                 aod::McCollisions const& mcCollisions)
  {
    processMc<true>(omegac0CandidatesMlMcKF, mcParticles, tracks, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegapi, processMcWithKFParticleMl, "Process MC with KFParticle and ML selections", false);
  // TODO: add the processMcWithKFParticleMl
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskOmegac0ToOmegapi>(cfgc)};
}
