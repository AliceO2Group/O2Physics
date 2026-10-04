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

/// \file taskOmegac0ToOmegaPiQa.cxx
/// \brief OmegaC0 analysis task
/// \author Yunfan Liu <yunfan.liu@cern.ch>, China University of Geosciences
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston
/// \author Maria Fernanda Torres Cabrera <maria.fernanda.torres.cabrera@cern.ch>, University of Houston

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
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
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Omegac0 analysis task
struct HfTaskOmegac0ToOmegaPiQa {

  Configurable<bool> selectionFlagOmegac0{"selectionFlagOmegac0", true, "Select Omegac0 candidates"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "Max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "Max. cand. rapidity"};

  SliceCache cache;

  using Omegac0Cands = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi>>;
  using Omegac0CandsKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfOmegacKf>>;
  using Omegac0CandsMc = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfToOmegaPiMCRec>>;
  using Omegac0CandsMcKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfOmegacKf, aod::HfToOmegaPiMCRec>>;
  using Omegac0CandsMl = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfMlSelOmegacToOmegaPi>>;
  using Omegac0CandsMlKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfMlSelOmegacToOmegaPi, aod::HfOmegacKf>>;
  using Omegac0CandsMlMc = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfMlSelOmegacToOmegaPi, aod::HfToOmegaPiMCRec>>;
  using Omegac0CandsMlMcKF = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfMlSelOmegacToOmegaPi, aod::HfOmegacKf, aod::HfToOmegaPiMCRec>>;

  using Omegac0Gen = soa::Filtered<soa::Join<aod::McParticles, aod::HfToOmegaPiMCGen>>;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  using McCollisionsWithFT0M = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterOmegaCToOmegaPiFlag = (aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi))) != static_cast<uint8_t>(0);
  Filter filterOmegaCMatchedRec = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi));
  Filter filterOmegaCMatchedGen = nabs(aod::hf_cand_mc_flag::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_xic0_omegac0::DecayType::OmegaczeroToOmegaPi));

  Preslice<Omegac0Cands> candOmegacPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
  Preslice<Omegac0CandsKF> candOmegacKFPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
  Preslice<Omegac0CandsMl> candOmegacMlPerCollision = aod::hf_cand_xic0_omegac0::collisionId;
  Preslice<Omegac0CandsMlKF> candOmegacKFMlPerCollision = aod::hf_cand_xic0_omegac0::collisionId;

  PresliceUnsorted<CollisionsWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {700, 2.4, 3.1}, "Cand. inv. mass"};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {500, 0, 50}, "Cand. pT"};
  ConfigurableAxis thnConfigAxisPtB{"thnConfigAxisPtB", {500, 0, 50}, "Cand. beauty mother pT"};
  ConfigurableAxis thnConfigAxisY{"thnConfigAxisY", {20, -1, 1}, "Cand. rapidity"};
  ConfigurableAxis thnConfigAxisOrigin{"thnConfigAxisOrigin", {3, -0.5, 2.5}, "Cand. origin"};
  ConfigurableAxis thnConfigAxisMatchFlag{"thnConfigAxisMatchFlag", {15, -7.5, 7.5}, "Cand. MC match flag"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Coll. num. PV contributors"};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0, 100}, "Coll. centrality precentile"};
  ConfigurableAxis thnConfigAxisPromptScore{"thnConfigAxisPromptScore", {100, 0, 1}, "Prompt score"};
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    std::array<bool, 20> doprocess{doprocessDataDCAFitter, doprocessDataDCAFitterMl, doprocessDataDCAFitterFT0C, doprocessDataDCAFitterMlFT0C,
                                   doprocessDataDCAFitterFT0M, doprocessDataDCAFitterMlFT0M,
                                   doprocessDataKFParticle, doprocessDataKFParticleMl, doprocessDataKFParticleFT0C, doprocessDataKFParticleMlFT0C,
                                   doprocessDataKFParticleFT0M, doprocessDataKFParticleMlFT0M,
                                   doprocessMcDCAFitter, doprocessMcDCAFitterMl,
                                   doprocessMcDCAFitterFT0M, doprocessMcDCAFitterMlFT0M,
                                   doprocessMcKFParticle, doprocessMcKFParticleMl,
                                   doprocessMcKFParticleFT0M, doprocessMcKFParticleMlFT0M};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "One and only one process function should be enabled at a time.");
    }

    const AxisSpec thnAxisMass{thnConfigAxisMass, "Inv. mass (#Omega#pi) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPtB{thnConfigAxisPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
    const AxisSpec thnAxisY{thnConfigAxisY, "y"};
    const AxisSpec thnAxisOrigin{thnConfigAxisOrigin, "Origin"};
    const AxisSpec thnAxisMatchFlag{thnConfigAxisMatchFlag, "MC match flag"};
    const AxisSpec thnAxisNumPvContr{thnConfigAxisNumPvContr, "Number of primary vtx. contributors"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality percentile"};
    const AxisSpec thnAxisCentMc{thnConfigAxisCent, "Centrality percentile (from gen. MC info)"};
    const AxisSpec thnAxisPromptScore{thnConfigAxisPromptScore, "BDT score prompt"};

    std::vector<AxisSpec> axes = {thnAxisMass, thnAxisPt, thnAxisY};
    std::vector<AxisSpec> axesMcGen = {thnAxisPt, thnAxisPtB, thnAxisY, thnAxisOrigin};

    if (doprocessDataDCAFitterFT0C || doprocessDataDCAFitterMlFT0C || doprocessDataDCAFitterFT0M || doprocessDataDCAFitterMlFT0M || doprocessDataKFParticleFT0C || doprocessDataKFParticleMlFT0C || doprocessDataKFParticleFT0M || doprocessDataKFParticleMlFT0M) {
      axes.push_back(thnAxisCent);
      axes.emplace_back(thnConfigAxisNumPvContr);
    }

    if (doprocessMcDCAFitterFT0M || doprocessMcDCAFitterMlFT0M || doprocessMcKFParticleFT0M || doprocessMcKFParticleMlFT0M) {
      axes.push_back(thnAxisCentMc);
      axes.emplace_back(thnConfigAxisNumPvContr);
      axesMcGen.push_back(thnAxisCentMc);
      axesMcGen.emplace_back(thnConfigAxisNumPvContr);
    }

    if (doprocessMcDCAFitter || doprocessMcDCAFitterMl || doprocessMcDCAFitterFT0M || doprocessMcDCAFitterMlFT0M || doprocessMcKFParticle || doprocessMcKFParticleMl || doprocessMcKFParticleFT0M || doprocessMcKFParticleMlFT0M) {
      registry.add("hMcGen", "Gen. #Omega_{c}^{0} from charm and beauty", HistType::kTHnSparseD, axesMcGen);
      registry.get<THnSparse>(HIST("hMcGen"))->Sumw2();

      if (doprocessMcDCAFitterFT0M || doprocessMcDCAFitterMlFT0M || doprocessMcKFParticleFT0M || doprocessMcKFParticleMlFT0M) {
        registry.add("hMcGenWithRecoColl", "Gen. #Omega_{c}^{0} from charm and beauty (associated to a reco collision)", HistType::kTHnSparseD, axesMcGen);
        registry.add("hNumRecoCollPerMcColl", "Number of reco collisions associated to a mc collision;Num. reco. coll. per Mc coll.;", {HistType::kTH1D, {{10, -0.5, 9.5}}});
        registry.get<THnSparse>(HIST("hMcGenWithRecoColl"))->Sumw2();
      }

      axes.push_back(thnAxisPtB);
      axes.push_back(thnAxisOrigin);
      axes.push_back(thnAxisMatchFlag);
    }

    if (doprocessDataKFParticleMl || doprocessDataKFParticleMlFT0C || doprocessDataKFParticleMlFT0M || doprocessDataDCAFitterMl || doprocessDataDCAFitterMlFT0C || doprocessDataDCAFitterMlFT0M || doprocessMcKFParticleMl || doprocessMcKFParticleMlFT0M || doprocessMcDCAFitterMl || doprocessMcDCAFitterMlFT0M) {
      axes.push_back(thnAxisPromptScore);
    }

    registry.add("hReco", "Reco. #Omega_{c}^{0} candidates", HistType::kTHnSparseD, axes);
    registry.get<THnSparse>(HIST("hReco"))->Sumw2();
  }

  template <bool UseKfParticle, bool ApplyMl, typename CandType>
  void processData(const CandType& candidates)
  {
    for (const auto& candidate : candidates) {
      if (!(candidate.resultSelections() == true || (candidate.resultSelections() == false && !selectionFlagOmegac0))) {
        continue;
      }

      double yOmegac{0.};
      if constexpr (UseKfParticle) {
        yOmegac = candidate.kfRapOmegac();
      } else {
        yOmegac = RecoDecay::y(std::array{candidate.pxCharmBaryon(), candidate.pyCharmBaryon(), candidate.pzCharmBaryon()}, o2::constants::physics::MassOmegaC0);
      }
      if (yCandRecoMax >= 0. && std::abs(yOmegac) > yCandRecoMax) {
        continue;
      }

      if constexpr (ApplyMl) {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac, candidate.mlProbOmegac()[0]);
      } else {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac);
      }
    }
  }

  template <bool UseKfParticle, bool UseCentrality, bool ApplyMl, typename CandType, typename CollType>
  void processDataCent(const CandType& candidate, CollType const& collision)
  {
    if (!(candidate.resultSelections() == true || (candidate.resultSelections() == false && !selectionFlagOmegac0))) {
      return;
    }

    double yOmegac{0.};
    if constexpr (UseKfParticle) {
      yOmegac = candidate.kfRapOmegac();
    } else {
      yOmegac = RecoDecay::y(std::array{candidate.pxCharmBaryon(), candidate.pyCharmBaryon(), candidate.pzCharmBaryon()}, o2::constants::physics::MassOmegaC0);
    }

    if (yCandRecoMax >= 0. && std::abs(yOmegac) > yCandRecoMax) {
      return;
    }

    auto numPvContributors = collision.numContrib();
    if constexpr (UseCentrality) {
      float const cent = o2::hf_centrality::getCentralityColl(collision);
      if constexpr (ApplyMl) {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac,
                      cent, numPvContributors, candidate.mlProbOmegac()[0]);
      } else {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac,
                      cent, numPvContributors);
      }
    } else {
      if constexpr (ApplyMl) {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac, candidate.mlProbOmegac()[0]);
      } else {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac);
      }
    }
  }

  template <bool UseKfParticle, bool ApplyMl, typename CandType>
  void processMc(const CandType& candidates, Omegac0Gen const& mcParticles)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (!(candidate.resultSelections() == true || (candidate.resultSelections() == false && !selectionFlagOmegac0))) {
        continue;
      }

      double yOmegac{0.};
      if constexpr (UseKfParticle) {
        yOmegac = candidate.kfRapOmegac();
      } else {
        yOmegac = RecoDecay::y(std::array{candidate.pxCharmBaryon(), candidate.pyCharmBaryon(), candidate.pzCharmBaryon()}, o2::constants::physics::MassOmegaC0);
      }

      if (yCandRecoMax >= 0. && std::abs(yOmegac) > yCandRecoMax) {
        continue;
      }

      if constexpr (ApplyMl) {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac, candidate.ptBhadMotherPart(), candidate.originMcRec(), candidate.flagMcMatchRec(), candidate.mlProbOmegac()[0]);
      } else {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac, candidate.ptBhadMotherPart(), candidate.originMcRec(), candidate.flagMcMatchRec());
      }
    }

    // MC gen.
    for (const auto& particle : mcParticles) {
      if (yCandGenMax >= 0. && std::abs(particle.rapidityCharmBaryonGen()) > yCandGenMax) {
        continue;
      }

      auto ptGen = particle.pt();
      auto yGen = particle.rapidityCharmBaryonGen();

      if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hMcGen"), ptGen, -1., yGen, RecoDecay::OriginType::Prompt);
      } else {
        float const ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
        registry.fill(HIST("hMcGen"), ptGen, ptGenB, yGen, RecoDecay::OriginType::NonPrompt);
      }
    }
  }

  template <bool UseKfParticle, bool ApplyMl, typename CandType, typename McCollisionWithCents>
  void processMcCent(const CandType& candidates, Omegac0Gen const& mcParticles,
                     CollisionsWithMcLabels const& collisions, McCollisionWithCents const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (!(candidate.resultSelections() == true || (candidate.resultSelections() == false && !selectionFlagOmegac0))) {
        continue;
      }

      double yOmegac{0.};
      if constexpr (UseKfParticle) {
        yOmegac = candidate.kfRapOmegac();
      } else {
        yOmegac = RecoDecay::y(std::array{candidate.pxCharmBaryon(), candidate.pyCharmBaryon(), candidate.pzCharmBaryon()}, o2::constants::physics::MassOmegaC0);
      }

      if (yCandRecoMax >= 0. && std::abs(yOmegac) > yCandRecoMax) {
        continue;
      }

      auto collision = candidate.template collision_as<CollisionsWithMcLabels>();
      uint16_t const numPvContributors = collision.numContrib();
      float const mcCent = o2::hf_centrality::getCentralityColl(collision.template mcCollision_as<McCollisionWithCents>());

      if constexpr (ApplyMl) {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac, mcCent, numPvContributors, candidate.ptBhadMotherPart(), candidate.originMcRec(), candidate.flagMcMatchRec(), candidate.mlProbOmegac()[0]);
      } else {
        registry.fill(HIST("hReco"), candidate.invMassCharmBaryon(), candidate.ptCharmBaryon(), yOmegac, mcCent, numPvContributors, candidate.ptBhadMotherPart(), candidate.originMcRec(), candidate.flagMcMatchRec());
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

      int maxNumContrib = 0;
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      for (const auto& recCol : recoCollsPerMcColl) {
        maxNumContrib = recCol.numContrib() > maxNumContrib ? recCol.numContrib() : maxNumContrib;
      }

      float const mcCent = o2::hf_centrality::getCentralityColl(mcCollision);

      if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hMcGen"), ptGen, -1., yGen, RecoDecay::OriginType::Prompt, mcCent, maxNumContrib);
      } else {
        float const ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
        registry.fill(HIST("hMcGen"), ptGen, ptGenB, yGen, RecoDecay::OriginType::NonPrompt, mcCent, maxNumContrib);
      }

      registry.fill(HIST("hNumRecoCollPerMcColl"), recoCollsPerMcColl.size());

      // fill sparse only for gen particles associated to a reconstructed collision
      if (recoCollsPerMcColl.size() >= 1) {
        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hMcGenWithRecoColl"), ptGen, -1., yGen, RecoDecay::OriginType::Prompt, mcCent, maxNumContrib);
        } else {
          float const ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          registry.fill(HIST("hMcGenWithRecoColl"), ptGen, ptGenB, yGen, RecoDecay::OriginType::NonPrompt, mcCent, maxNumContrib);
        }
      }
    }
  }

  ///////////////////////////////////
  ///    Data with DCAFitter       //
  ///////////////////////////////////

  void processDataDCAFitter(Omegac0Cands const& candidates)
  {
    processData<false, false>(candidates);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataDCAFitter, "process data with DCAFitter", false);

  void processDataDCAFitterMl(Omegac0CandsMl const& candidates)
  {
    processData<false, true>(candidates);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataDCAFitterMl, "process data with DCAFitter, ML selections", false);

  void processDataDCAFitterFT0C(Omegac0Cands const& candidates,
                                CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<false, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataDCAFitterFT0C, "process data with DCAFitter, FT0C centrality", false);

  void processDataDCAFitterMlFT0C(Omegac0CandsMl const& candidates,
                                  CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacMlPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<false, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataDCAFitterMlFT0C, "process data with DCAFitter, ML selections, FT0C centrality", false);

  void processDataDCAFitterFT0M(Omegac0Cands const& candidates,
                                CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<false, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataDCAFitterFT0M, "process data with DCAFitter, FT0M centrality", false);

  void processDataDCAFitterMlFT0M(Omegac0CandsMl const& candidates,
                                  CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacMlPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<false, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataDCAFitterMlFT0M, "process data with DCAFitter, ML selections, FT0M centrality", false);

  ////////////////////////////////////
  ///    Data with KFParticle       //
  ////////////////////////////////////

  void processDataKFParticle(Omegac0CandsKF const& candidates)
  {
    processData<true, false>(candidates);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataKFParticle, "process data with KFParticle", false);

  void processDataKFParticleMl(Omegac0CandsMlKF const& candidates)
  {
    processData<true, true>(candidates);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataKFParticleMl, "process data with KFParticle, ML selections", false);

  void processDataKFParticleFT0C(Omegac0CandsKF const& candidates,
                                 CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacKFPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<true, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataKFParticleFT0C, "process data with KFParticle, FT0C centrality", false);

  void processDataKFParticleMlFT0C(Omegac0CandsMlKF const& candidates,
                                   CollisionsWithFT0C const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacKFMlPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<true, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataKFParticleMlFT0C, "process data with KFParticle, ML selections, FT0C centrality", false);

  void processDataKFParticleFT0M(Omegac0CandsKF const& candidates,
                                 CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacKFPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<true, true, false>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataKFParticleFT0M, "process data with KFParticle, FT0M centrality", false);

  void processDataKFParticleMlFT0M(Omegac0CandsMlKF const& candidates,
                                   CollisionsWithFT0M const& collisions)
  {
    for (const auto& collision : collisions) {
      auto groupedOmegacCandidates = candidates.sliceBy(candOmegacKFMlPerCollision, collision.globalIndex());

      for (const auto& candidate : groupedOmegacCandidates) {
        processDataCent<true, true, true>(candidate, collision);
      }
    }
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processDataKFParticleMlFT0M, "process data with KFParticle, ML selections, FT0M centrality", false);

  ///////////////////////////////////
  ///    MC with DCAFitter         //
  ///////////////////////////////////

  void processMcDCAFitter(Omegac0CandsMc const& omegaC0CandidatesMc,
                          Omegac0Gen const& mcParticles)
  {
    processMc<false, false>(omegaC0CandidatesMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcDCAFitter, "Process MC with DCAFitter", false);

  void processMcDCAFitterMl(Omegac0CandsMlMc const& omegac0CandidatesMlMc,
                            Omegac0Gen const& mcParticles)
  {
    processMc<false, true>(omegac0CandidatesMlMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcDCAFitterMl, "Process MC with DCAFitter, ML selections", false);

  void processMcDCAFitterFT0M(Omegac0CandsMc const& omegaC0CandidatesMc,
                              Omegac0Gen const& mcParticles,
                              CollisionsWithMcLabels const& collisions,
                              McCollisionsWithFT0M const& mcCollisions)
  {
    processMcCent<false, false>(omegaC0CandidatesMc, mcParticles, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcDCAFitterFT0M, "Process MC with DCAFitter, FT0M centrality (from MC)", false);

  void processMcDCAFitterMlFT0M(Omegac0CandsMlMc const& omegac0CandidatesMlMc,
                                Omegac0Gen const& mcParticles,
                                CollisionsWithMcLabels const& collisions,
                                McCollisionsWithFT0M const& mcCollisions)
  {
    processMcCent<false, true>(omegac0CandidatesMlMc, mcParticles, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcDCAFitterMlFT0M, "Process MC with DCAFitter, ML selections, FT0M centrality (from MC)", false);

  ////////////////////////////////////
  ///    MC with KFParticle         //
  ////////////////////////////////////

  void processMcKFParticle(Omegac0CandsMcKF const& omegaC0CandidatesMcKF,
                           Omegac0Gen const& mcParticles)
  {
    processMc<true, false>(omegaC0CandidatesMcKF, mcParticles);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcKFParticle, "Process MC with KFParticle", false);

  void processMcKFParticleMl(Omegac0CandsMlMcKF const& omegac0CandidatesMlMcKF,
                             Omegac0Gen const& mcParticles)
  {
    processMc<true, true>(omegac0CandidatesMlMcKF, mcParticles);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcKFParticleMl, "Process MC with KFParticle, ML selections", false);

  void processMcKFParticleFT0M(Omegac0CandsMcKF const& omegaC0CandidatesMcKF,
                               Omegac0Gen const& mcParticles,
                               CollisionsWithMcLabels const& collisions,
                               McCollisionsWithFT0M const& mcCollisions)
  {
    processMcCent<true, false>(omegaC0CandidatesMcKF, mcParticles, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcKFParticleFT0M, "Process MC with KFParticle, FT0M centrality (from MC)", false);

  void processMcKFParticleMlFT0M(Omegac0CandsMlMcKF const& omegac0CandidatesMlMcKF,
                                 Omegac0Gen const& mcParticles,
                                 CollisionsWithMcLabels const& collisions,
                                 McCollisionsWithFT0M const& mcCollisions)
  {
    processMcCent<true, true>(omegac0CandidatesMlMcKF, mcParticles, collisions, mcCollisions);
  }
  PROCESS_SWITCH(HfTaskOmegac0ToOmegaPiQa, processMcKFParticleMlFT0M, "Process MC with KFParticle, ML selections, FT0M centrality (from MC)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskOmegac0ToOmegaPiQa>(cfgc)};
}
