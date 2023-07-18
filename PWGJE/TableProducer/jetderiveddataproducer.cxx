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

// task to produce a self contained data format for jet analyses from the full AO2D
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "EventFiltering/filterTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetDerivedDataProducerTask {
  Produces<aod::JCollisions> jetCollisionsTable;
  Produces<aod::JCollisionPIs> jetCollisionsParentIndexTable;
  Produces<aod::JMcCollisionLbs> jetMcCollisionsLabelTable;
  Produces<aod::JMcCollisions> jetMcCollisionsTable;
  Produces<aod::JMcCollisionPIs> jetMcCollisionsParentIndexTable;
  Produces<aod::JChTrigSels> jetChargedTriggerSelsTable;
  Produces<aod::JTracks> jetTracksTable;
  Produces<aod::JTrackPIs> jetTracksParentIndexTable;
  Produces<aod::JMcTrackLbs> jetMcTracksLabelTable;
  Produces<aod::JMcParticles> jetMcParticlesTable;
  Produces<aod::JMcParticlePIs> jetParticlesParentIndexTable;

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillJetCollisionsTable(T const& collision)
  {
    jetCollisionsTable(collision.posZ(), JetDerivedDataUtilities::setEventSelectionBit(collision), collision.alias());
    jetCollisionsParentIndexTable(collision.globalIndex());
  }

  template <typename T>
  void fillJetMcCollisionsTable(T const& McCollision)
  {
    jetMcCollisionsTable(McCollision.posZ(), McCollision.weight());
    jetMcCollisionsParentIndexTable(McCollision.globalIndex());
  }

  template <typename T, typename U>
  void fillJetTracksTable(T const& collision, U const& track)
  {
    jetTracksTable(collision.globalIndex(), track.pt(), track.eta(), track.phi(), JetDerivedDataUtilities::trackEnergy(track), JetDerivedDataUtilities::setTrackSelectionBit(track));
    jetTracksParentIndexTable(track.globalIndex());
  }

  template <typename T, typename U>
  void fillJetParticlesTable(T const& collision, U const& particle)
  {
    std::vector<int> mothersId;
    int daughtersId[2];
    if (particle.has_mothers()) {
      for (auto const& mother : particle.template mothers_as<aod::McParticles>()) {
        mothersId.push_back(mother.globalIndex());
      }
    }
    auto i = 0;
    if (particle.has_daughters()) {
      for (auto const& daughter : particle.template daughters_as<aod::McParticles>()) {
        if (i > 1) {
          break;
        }
        daughtersId[i] = daughter.globalIndex();
        i++;
      }
    }
    jetMcParticlesTable(collision.globalIndex(), particle.pt(), particle.eta(), particle.phi(), particle.y(), particle.e(), particle.pdgCode(), particle.getGenStatusCode(), particle.getHepMCStatusCode(), particle.isPhysicalPrimary(), mothersId, daughtersId);
    jetParticlesParentIndexTable(particle.globalIndex());
  }

  void processChargedData(soa::Join<aod::Collisions, aod::EvSels, aod::BcSels>::iterator const& collision,
                          soa::Join<aod::Tracks, aod::TrackSelection> const& tracks)
  {
    fillJetCollisionsTable(collision);
    for (auto& track : tracks) {
      fillJetTracksTable(collision, track);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processChargedData, "jet track table producer data", true);

  void processChargedTriggered(soa::Join<aod::Collisions, aod::JetFilters>::iterator const& collision)
  {
    jetChargedTriggerSelsTable(JetDerivedDataUtilities::setChargedTriggerSelectionBit(collision));
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processChargedTriggered, "jet track table producer data", false);

  void processChargedMCD(soa::Join<aod::Collisions, aod::EvSels, aod::BcSels, aod::McCollisionLabels>::iterator const& collision,
                         soa::Join<aod::Tracks, aod::TrackSelection, aod::McTrackLabels> const& tracks,
                         aod::McParticles const& particles)
  {
    fillJetCollisionsTable(collision);
    if (collision.has_mcCollision()) {
      jetMcCollisionsLabelTable(collision.mcCollision().globalIndex());
    } else {
      jetMcCollisionsLabelTable(-1);
    }
    for (auto& track : tracks) {
      fillJetTracksTable(collision, track);
      if (track.has_mcParticle()) {
        jetMcTracksLabelTable(track.mcParticle().globalIndex());
      } else {
        jetMcTracksLabelTable(-1);
      }
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processChargedMCD, "jet track table producer MC Det", false);

  void processChargedMCP(aod::McCollision const& McCollision,
                         aod::McParticles const& particles)
  {
    fillJetMcCollisionsTable(McCollision);
    for (auto& particle : particles) {
      fillJetParticlesTable(McCollision, particle);
    }
  }
  PROCESS_SWITCH(JetDerivedDataProducerTask, processChargedMCP, "jet track table producer MC Part", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetDerivedDataProducerTask>(cfgc, TaskName{"jettracks-creator"})};
}