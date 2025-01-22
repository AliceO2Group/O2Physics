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

/// \file heavyFlavourDefinition.cxx
/// \brief Task to produce a table joinable to the jet tables for a flavour definition on MC
/// \author Hanseo Park <hanseo.park@cern.ch>

#include <memory>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/trackUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableMCD, typename JetTableMCP, typename JetFlavourDefMCD, typename JetFlavourDefMCP>
struct HeavyFlavourDefinitionTask {

  Produces<JetFlavourDefMCD> flavourTableMCD;
  Produces<JetFlavourDefMCP> flavourTableMCP;

  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};
  Configurable<bool> searchUpToQuark{"searchUpToQuark", true, "Finding first mother in particles to quark"};

  using JetTracksMCD = soa::Join<aod::JetTracksMCD, aod::JTrackPIs>;
  Preslice<aod::JetParticles> particlesPerCollision = aod::jmcparticle::mcCollisionId;
  Preslice<soa::Join<aod::JMcParticles, aod::JMcParticlePIs>> particlesPerMcCollision = aod::jmcparticle::mcCollisionId;

  void init(InitContext const&)
  {
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processDummy, "Dummy process", true);

  void processMCDByConstituents(aod::JetCollision const& /*collision*/, JetTableMCD const& mcdjets, JetTracksMCD const& tracks, aod::JetParticles const& particles)
  {
    for (auto const& mcdjet : mcdjets) {
      int8_t origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, tracks, particles, maxDeltaR, searchUpToQuark);
      flavourTableMCD(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCDByConstituents, "Fill definition of flavour for mcd jets using constituents", false);

  void processMCDByDistance(soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>::iterator const& collision, soa::Join<JetTableMCD, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcdjets, soa::Join<JetTableMCP, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& /*mcpjets*/, aod::JetParticles const& particles) // it used only for charged jets now
  {
    for (auto const& mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      int8_t origin = -1;
      if (mcdjet.has_matchedJetGeo()) {
        for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<soa::Join<JetTableMCP, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>()) {
          if (searchUpToQuark) {
            origin = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
          } else {
            origin = jettaggingutilities::getJetFlavorHadron(mcpjet, particlesPerColl);
          }
        }
      } else {
        origin = JetTaggingSpecies::none;
      }
      flavourTableMCD(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCDByDistance, "Fill definition of flavour for mcd jets using distance of jet with particles", false);

  void processMCPByConstituents(JetTableMCP const& mcpjets, aod::JetParticles const& particles)
  {
    for (auto const& mcpjet : mcpjets) {
      int8_t origin = jettaggingutilities::mcpJetFromHFShower(mcpjet, particles, maxDeltaR, searchUpToQuark);
      flavourTableMCP(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCPByConstituents, "Fill definition of flavour for mcp jets using constituents", false);

  void processMCPByDistance(soa::Join<aod::JMcCollisions, aod::JMcCollisionPIs> const& /*mcCollisions*/, JetTableMCP const& mcpjets, aod::JetParticles const& particles)
  {
    for (auto const& mcpjet : mcpjets) {
      auto const particlesPerMcColl = particles.sliceBy(particlesPerMcCollision, mcpjet.globalIndex());
      int8_t origin = -1;
      if (searchUpToQuark) {
        origin = jettaggingutilities::getJetFlavor(mcpjet, particlesPerMcColl);
      } else {
        origin = jettaggingutilities::getJetFlavorHadron(mcpjet, particlesPerMcColl);
      }
      flavourTableMCP(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCPByDistance, "Fill definition of flavour for mcp jets using distance of jet with particles", false);
};

using JetHfDefinitionCharged = HeavyFlavourDefinitionTask<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCParticleLevelJetFlavourDef>;
// using JetHfDefinitionFull = HeavyFlavourDefinitionTask<soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>, aod::FullMCDetectorLevelJetFlavourDef, aod::FullMCParticleLevelJetFlavourDef>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetHfDefinitionCharged>(cfgc, TaskName{"jet-hf-definition-charged"})}; // o2-linter: disable=name/o2-task
}
