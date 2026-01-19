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

// jet tutorial task for hands on tutorial session (09/11/2023)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "PWGHF/Utils/utilsMcGen.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct HFDebugTask {
  HistogramRegistry registry{"registry",
                             {{"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
                              {"h_mccollisions", "mc event status;event status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
                              {"h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_mcpLcjet_nTracks", "mcp Lc jet ntracks;n_{track};entries", {HistType::kTH1F, {{80, 0.0, 80.}}}},
                              {"h_mcpLcjet_angularity", "mcp Lc jet angularity;angularity;entries", {HistType::kTH1F, {{50, 0.0, 1.}}}}}};

  Produces<o2::aod::HfD0Bases> D0BaseTable;
  Produces<o2::aod::HfD0CollBases> D0CollBaseTable;
  Produces<o2::aod::HfD0CollIds> D0CollIdTable;
  Produces<o2::aod::HfD0McCollBases> D0McCollBaseTable;
  Produces<o2::aod::HfD0McCollIds> D0McCollIdTable;
  Produces<o2::aod::HfD0McRCollIds> D0McRCollIdTable;
  Produces<o2::aod::HfD0PBases> D0PBaseTable;
  Produces<o2::aod::HfD0PIds> D0PIdTable;
  Produces<o2::aod::HfD0Pars> D0ParTable;
  Produces<o2::aod::HfD0ParEs> D0ParETable;
  Produces<o2::aod::HfD0Sels> D0SelTable;
  Produces<o2::aod::HfD0Mls> D0MlTable;
  Produces<o2::aod::HfD0Ids> D0IdTable;
  Produces<o2::aod::HfD0Mcs> D0McTable;

  Produces<o2::aod::HfLcBases> LcBaseTable;
  Produces<o2::aod::HfLcCollBases> LcCollBaseTable;
  Produces<o2::aod::HfLcCollIds> LcCollIdTable;
  Produces<o2::aod::HfLcMcCollBases> LcMcCollBaseTable;
  Produces<o2::aod::HfLcMcCollIds> LcMcCollIdTable;
  Produces<o2::aod::HfLcMcRCollIds> LcMcRCollIdTable;
  Produces<o2::aod::HfLcPBases> LcPBaseTable;
  Produces<o2::aod::HfLcPIds> LcPIdTable;
  Produces<o2::aod::HfLcPars> LcParTable;
  Produces<o2::aod::HfLcParEs> LcParETable;
  Produces<o2::aod::HfLcSels> LcSelTable;
  Produces<o2::aod::HfLcMls> LcMlTable;
  Produces<o2::aod::HfLcIds> LcIdTable;
  Produces<o2::aod::HfLcMcs> LcMcTable;

  Produces<o2::aod::BkgChargedRhos> rhoTable;
  Produces<o2::aod::BkgChargedMcRhos> rhoMCTable;
  Produces<o2::aod::BkgD0Rhos> rhoD0Table;
  Produces<o2::aod::BkgD0McRhos> rhoD0MCTable;
  Produces<o2::aod::BkgLcRhos> rhoLcTable;
  Produces<o2::aod::BkgLcMcRhos> rhoLcMCTable;

  std::vector<float> pTHF;
  std::vector<float> etaHF;
  std::vector<float> phiHF;
  std::vector<float> yHF;
  const int nCycles = 400;
  int mcCollisionIdD0 = -1;
  int mcCollisionIdLc = -1;
  void init(o2::framework::InitContext&)
  {
    pTHF.reserve(nCycles);
    etaHF.reserve(nCycles);
    phiHF.reserve(nCycles);
    yHF.reserve(nCycles);

    for (int i = 0; i < nCycles; ++i) {
      pTHF.push_back(0.125 * (i + 1));
      etaHF.push_back(-0.7 + (1.4 / 400.0 * i));
      phiHF.push_back(0.0 + (6.28 / 400.0 * i));
      yHF.push_back(-0.7 + (1.4 / 400.0 * i));
    }
  }

  Preslice<aod::McParticles> perMcCollisionParticles = aod::mcparticle::mcCollisionId;

  void clearObjetcs(aod::Collisions const& collisions)
  {
    mcCollisionIdD0 = -1;
    mcCollisionIdLc = -1;
    for (auto i = 0; i < collisions.size(); i++) {
      rhoTable(0.0, 0.0);
    }
  }
  PROCESS_SWITCH(HFDebugTask, clearObjetcs, "clear", true);

  void processMcRhoTable(aod::McCollision const&)
  {
    rhoMCTable(0.0, 0.0);
  }
  PROCESS_SWITCH(HFDebugTask, processMcRhoTable, "mc rho table", false);

  void processD0MC(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&, aod::McParticles const& particles, soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks)
  {
    bool foundD0 = false;
    std::vector<aod::McParticle> D0Particles;
    auto const& mcCollision = collision.mcCollision();
    if (mcCollisionIdD0 > mcCollision.globalIndex()) {
      return;
    } else {
      mcCollisionIdD0 = mcCollision.globalIndex();
    }
    auto particleSlice = particles.sliceBy(perMcCollisionParticles, mcCollision.globalIndex());
    for (auto const& particle : particleSlice) {
      if (std::abs(particle.pdgCode()) == 421 && std::abs(particle.y()) < 0.9) {
        foundD0 = true;
        D0Particles.push_back(particle);
      }
    }
    if (!foundD0) {
      return;
    }

    if (tracks.size() == 0) {
      return;
    }

    for (auto const& D0Particle : D0Particles) {
      int daughter1Id = -1;
      int daughter2Id = -1;
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto const& particle = track.mcParticle();
        if (particle.has_mothers() && particle.mothers_first_as<aod::McParticles>().globalIndex() == D0Particle.globalIndex()) {
          if (daughter1Id == -1) {
            daughter1Id = track.globalIndex();
          } else {
            daughter2Id = track.globalIndex();
            break;
          }
        }
      }
      if (daughter1Id == -1 || daughter2Id == -1) {
        continue;
      }
      D0BaseTable(collision.globalIndex(), D0Particle.pt(), D0Particle.eta(), D0Particle.phi(), static_cast<float>(o2::constants::physics::MassD0), D0Particle.y());
      D0CollBaseTable(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), -1.0, -1.0, -1.0, -1.0, -1.0);
      D0CollIdTable(collision.globalIndex());
      D0ParTable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      D0ParETable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      D0SelTable(0.0);
      D0MlTable(std::vector<float>{0.9f, 0.1f, 0.8f});
      D0IdTable(collision.globalIndex(), daughter1Id, daughter2Id);
      D0McCollBaseTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), 1.0);
      D0McCollIdTable(mcCollision.globalIndex());
      std::vector<int32_t> collIDs;
      collIDs.push_back(D0CollBaseTable.lastIndex());
      D0McRCollIdTable(collIDs);
      D0PBaseTable(mcCollision.globalIndex(), D0Particle.pt(), D0Particle.eta(), D0Particle.phi(), D0Particle.y(), static_cast<int>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK), 0);
      D0PIdTable(mcCollision.globalIndex(), D0Particle.globalIndex());
      D0McTable(static_cast<int>(o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK), 0);
      rhoD0Table(0.0, 0.0);
      rhoD0MCTable(0.0, 0.0);
    }
  }
  PROCESS_SWITCH(HFDebugTask, processD0MC, "D0 mc", false);

  void processD0Data(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Tracks const& tracks)
  {
    if (tracks.size() == 0) {
      return;
    }
    D0BaseTable(collision.globalIndex(), pTHF[collision.globalIndex() % nCycles], etaHF[collision.globalIndex() % nCycles], phiHF[collision.globalIndex() % nCycles], static_cast<float>(o2::constants::physics::MassD0), yHF[collision.globalIndex() % nCycles]);
    D0CollBaseTable(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), -1.0, -1.0, -1.0, -1.0, -1.0);
    D0CollIdTable(collision.globalIndex());
    D0ParTable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    D0ParETable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    D0SelTable(0.0); // this might not work
    D0MlTable(std::vector<float>{0.9, 0.9, 0.9});
    rhoD0Table(0.0, 0.0);
    for (const auto& track : tracks) {
      D0IdTable(collision.globalIndex(), track.globalIndex(), track.globalIndex());
      break;
    }
  }
  PROCESS_SWITCH(HFDebugTask, processD0Data, "D0data", true);

  void processLcMC(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&, aod::McParticles const& particles, soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks)
  {
    bool foundLc = false;
    std::vector<aod::McParticle> LcParticles;
    auto const& mcCollision = collision.mcCollision();
    if (mcCollisionIdLc > mcCollision.globalIndex()) {
      return;
    } else {
      mcCollisionIdLc = mcCollision.globalIndex();
    }
    auto particleSlice = particles.sliceBy(perMcCollisionParticles, mcCollision.globalIndex());
    for (auto const& particle : particleSlice) {
      if (std::abs(particle.pdgCode()) == 4122 && std::abs(particle.y()) < 0.9) {
        foundLc = true;
        LcParticles.push_back(particle);
      }
    }
    if (!foundLc) {
      return;
    }

    if (tracks.size() == 0) {
      return;
    }

    for (auto const& LcParticle : LcParticles) {
      int daughter1Id = -1;
      int daughter2Id = -1;
      int daughter3Id = -1;
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        auto const& particle = track.mcParticle();
        if (particle.has_mothers() && particle.mothers_first_as<aod::McParticles>().globalIndex() == LcParticle.globalIndex()) {
          if (daughter1Id == -1) {
            daughter1Id = track.globalIndex();
          } else if (daughter2Id == -1) {
            daughter2Id = track.globalIndex();
          } else {
            daughter3Id = track.globalIndex();
            break;
          }
        }
      }
      if (daughter1Id == -1 || daughter2Id == -1 || daughter3Id == -1) {
        continue;
      }
      LcBaseTable(collision.globalIndex(), LcParticle.pt(), LcParticle.eta(), LcParticle.phi(), static_cast<float>(o2::constants::physics::MassLambdaCPlus), LcParticle.y());
      LcCollBaseTable(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), -1.0, -1.0, -1.0, -1.0, -1.0);
      LcCollIdTable(collision.globalIndex());
      LcParTable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      LcParETable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      LcSelTable(0.0);
      LcMlTable(std::vector<float>{0.9f, 0.1f, 0.8f});
      LcIdTable(collision.globalIndex(), daughter1Id, daughter2Id, daughter3Id);
      LcMcCollBaseTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), 1.0);
      LcMcCollIdTable(mcCollision.globalIndex());
      std::vector<int32_t> collIDs;
      collIDs.push_back(LcCollBaseTable.lastIndex());
      LcMcRCollIdTable(collIDs);
      LcPBaseTable(mcCollision.globalIndex(), LcParticle.pt(), LcParticle.eta(), LcParticle.phi(), LcParticle.y(), static_cast<int>(o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi), 0);
      LcPIdTable(mcCollision.globalIndex(), LcParticle.globalIndex());
      LcMcTable(static_cast<int>(o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi), 0, 0);
      rhoLcTable(0.0, 0.0);
      rhoLcMCTable(0.0, 0.0);
    }
  }
  PROCESS_SWITCH(HFDebugTask, processLcMC, "Lc mc", false);

  void processLcData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Tracks const& tracks)
  {
    if (tracks.size() == 0) {
      return;
    }
    LcBaseTable(collision.globalIndex(), pTHF[collision.globalIndex() % nCycles], etaHF[collision.globalIndex() % nCycles], phiHF[collision.globalIndex() % nCycles], static_cast<float>(o2::constants::physics::MassLambdaCPlus), yHF[collision.globalIndex() % nCycles]);
    LcCollBaseTable(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), -1.0, -1.0, -1.0, -1.0, -1.0);
    LcCollIdTable(collision.globalIndex());
    LcParTable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    LcParETable(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    LcSelTable(0.0); // this might not work
    LcMlTable(std::vector<float>{0.9, 0.9, 0.9});
    rhoLcTable(0.0, 0.0);
    for (const auto& track : tracks) {
      LcIdTable(collision.globalIndex(), track.globalIndex(), track.globalIndex(), track.globalIndex());
      break;
    }
  }
  PROCESS_SWITCH(HFDebugTask, processLcData, "Lcdata", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<HFDebugTask>(cfgc, TaskName{"hf-debug"})}; }
