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

#include <TF1.h>
#include <TH1.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/trackUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename JetTableMCD, typename JetTableMCP, typename JetFlavourDefMCD, typename JetFlavourDefMCP>
struct HeavyFlavourDefinitionTask {

  Produces<JetFlavourDefMCD> flavourTableMCD;
  Produces<JetFlavourDefMCP> flavourTableMCP;

  Configurable<float> maxDeltaR{"maxDeltaR", 0.25, "maximum distance of jet axis from flavour initiating parton"};
  Configurable<bool> searchUpToQuark{"searchUpToQuark", true, "Finding first mother in particles to quark"};

  using JetTracksMCD = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;
  Preslice<aod::JetParticles> particlesPerCollision = aod::jmcparticle::mcCollisionId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processDummy, "Dummy process", true);

  void processMCD(aod::JetCollision const& /*collision*/, JetTableMCD const& mcdjets, JetTracksMCD const& jtracks, aod::JetParticles const& particles)
  {
    for (auto const& mcdjet : mcdjets) {
      int origin = jettaggingutilities::mcdJetFromHFShower(mcdjet, jtracks, particles, maxDeltaR, searchUpToQuark);
      flavourTableMCD(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCD, "Fill definition of flavour for mcd jets", true);

  void processMCDRun2(soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>::iterator const& collision, soa::Join<JetTableMCD, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& mcdjets, soa::Join<JetTableMCP, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& /*mcpjets*/, aod::JetParticles const& particles) // it used only for charged jets now
  {
    for (auto const& mcdjet : mcdjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, collision.mcCollisionId());
      int origin = -1;
      for (auto const& mcpjet : mcdjet.template matchedJetGeo_as<soa::Join<JetTableMCP, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>()) {
        if (searchUpToQuark) {
          origin = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
        } else {
          origin = jettaggingutilities::getJetFlavorHadron(mcpjet, particlesPerColl);
        }
      }
      flavourTableMCD(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCDRun2, "Fill definition of flavour for mcd jets (Run2)", false);

  void processMCP(JetTableMCP const& mcpjets, aod::JetParticles const& particles)
  {
    for (auto const& mcpjet : mcpjets) {
      int origin = jettaggingutilities::mcpJetFromHFShower(mcpjet, particles, maxDeltaR, searchUpToQuark);
      flavourTableMCP(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCP, "Fill definition of flavour for mcp jets", true);

  void processMCPRun2(JetTableMCP const& mcpjets, aod::JetParticles const& particles)
  {
    for (auto const& mcpjet : mcpjets) {
      auto const particlesPerColl = particles.sliceBy(particlesPerCollision, mcpjet.mcCollisionId());
      int origin = -1;
      if (searchUpToQuark) {
        origin = jettaggingutilities::getJetFlavor(mcpjet, particlesPerColl);
      } else {
        origin = jettaggingutilities::getJetFlavorHadron(mcpjet, particlesPerColl);
      }
      flavourTableMCP(origin);
    }
  }
  PROCESS_SWITCH(HeavyFlavourDefinitionTask, processMCPRun2, "Fill definition of flavour for mcp jets (Run2)", false);
};

using JetFlavourDefCharged = HeavyFlavourDefinitionTask<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCParticleLevelJetFlavourDef>;
using JetFlavourDefFull = HeavyFlavourDefinitionTask<soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>, soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>, aod::FullMCDetectorLevelJetFlavourDef, aod::FullMCParticleLevelJetFlavourDef>;
// using JetTaggerhfNeutral = HeavyFlavourDefinitionTask<soa::Join<aod::NeutralMCDetectorLevelJets, aod::NeutralMCDetectorLevelJetConstituents>, aod::NeutralMCDetectorLevelJetTags, aod::NeutralMCParticleLevelJetFlavourDef>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(
    adaptAnalysisTask<JetFlavourDefCharged>(cfgc));

  tasks.emplace_back(
    adaptAnalysisTask<JetFlavourDefFull>(cfgc));
  /*
    tasks.emplace_back(
      adaptAnalysisTask<JetFlavourDefNeutral>(cfgc,
                                                  SetDefaultProcesses{}));
  */
  return WorkflowSpec{tasks};
}
