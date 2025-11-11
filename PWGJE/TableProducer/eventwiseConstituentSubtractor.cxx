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

// Task to produce a table groupable with the jcollisions table which contains the event-wise constituent subtracted track list
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <fastjet/PseudoJet.hh>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct eventWiseConstituentSubtractorTask {
  Produces<aod::JTrackSubs> trackSubtractedTable;
  Produces<aod::JMcParticleSubs> particleSubtractedTable;
  Produces<aod::JTrackD0Subs> trackSubtractedD0Table;
  Produces<aod::JMcParticleD0Subs> particleSubtractedD0Table;
  Produces<aod::JTrackDplusSubs> trackSubtractedDplusTable;
  Produces<aod::JMcParticleDplusSubs> particleSubtractedDplusTable;
  Produces<aod::JTrackDsSubs> trackSubtractedDsTable;
  Produces<aod::JMcParticleDsSubs> particleSubtractedDsTable;
  Produces<aod::JTrackDstarSubs> trackSubtractedDstarTable;
  Produces<aod::JMcParticleDstarSubs> particleSubtractedDstarTable;
  Produces<aod::JTrackLcSubs> trackSubtractedLcTable;
  Produces<aod::JMcParticleLcSubs> particleSubtractedLcTable;
  Produces<aod::JTrackB0Subs> trackSubtractedB0Table;
  Produces<aod::JMcParticleB0Subs> particleSubtractedB0Table;
  Produces<aod::JTrackBplusSubs> trackSubtractedBplusTable;
  Produces<aod::JMcParticleBplusSubs> particleSubtractedBplusTable;
  Produces<aod::JTrackXicToXiPiPiSubs> trackSubtractedXicToXiPiPiTable;
  Produces<aod::JMcParticleXicToXiPiPiSubs> particleSubtractedXicToXiPiPiTable;
  Produces<aod::JTrackDielectronSubs> trackSubtractedDielectronTable;
  Produces<aod::JMcParticleDielectronSubs> particleSubtractedDielectronTable;

  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};

  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<bool> applyTrackingEfficiency{"applyTrackingEfficiency", {false}, "configurable to decide whether to apply artificial tracking efficiency (discarding tracks) in the collision analysed by this task"};
  Configurable<std::vector<double>> trackingEfficiencyPtBinning{"trackingEfficiencyPtBinning", {0., 10, 999.}, "pt binning of tracking efficiency array if applyTrackingEfficiency is true"};
  Configurable<std::vector<double>> trackingEfficiency{"trackingEfficiency", {1.0, 1.0}, "tracking efficiency array applied if applyTrackingEfficiency is true"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  Configurable<float> alpha{"alpha", 1.0, "exponent of transverse momentum in calculating the distance measure between pairs"};
  Configurable<float> rMax{"rMax", 0.24, "maximum distance of subtraction"};
  Configurable<float> eventEtaMax{"eventEtaMax", 0.9, "maximum pseudorapidity of event"};
  Configurable<bool> doRhoMassSub{"doRhoMassSub", true, "perfom mass subtraction as well"};

  JetBkgSubUtils eventWiseConstituentSubtractor;
  std::vector<fastjet::PseudoJet> inputParticles;
  std::vector<fastjet::PseudoJet> tracksSubtracted;
  int trackSelection = -1;

  std::string particleSelection;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    eventWiseConstituentSubtractor.setDoRhoMassSub(doRhoMassSub);
    eventWiseConstituentSubtractor.setConstSubAlphaRMax(alpha, rMax);
    eventWiseConstituentSubtractor.setMaxEtaEvent(eventEtaMax);

    if (applyTrackingEfficiency) {
      if (trackingEfficiencyPtBinning->size() < 2) {
        LOGP(fatal, "eventWiseConstituentSubtractor workflow: trackingEfficiencyPtBinning configurable should have at least two bin edges");
      }
      if (trackingEfficiency->size() + 1 != trackingEfficiencyPtBinning->size()) {
        LOGP(fatal, "eventWiseConstituentSubtractor workflow: trackingEfficiency configurable should have exactly one less entry than the number of bin edges set in trackingEfficiencyPtBinning configurable");
      }
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta >= trackEtaMin && aod::jmcparticle::eta <= trackEtaMax && aod::jmcparticle::phi >= trackPhiMin && aod::jmcparticle::phi <= trackPhiMax);

  template <typename T, typename U, typename V>
  void analyseHF(T const& tracks, U const& candidates, V& trackSubTable)
  {
    for (auto& candidate : candidates) {
      inputParticles.clear();
      tracksSubtracted.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning, &candidate);

      tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, candidate.rho(), candidate.rhoM());
      for (auto const& trackSubtracted : tracksSubtracted) {
        trackSubTable(candidate.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), jetderiveddatautilities::setSingleTrackSelectionBit(trackSelection));
      }
    }
  }

  template <typename T, typename U, typename V>
  void analyseHFMc(T const& particles, U const& candidates, V& particleSubTable)
  {
    for (auto& candidate : candidates) {
      inputParticles.clear();
      tracksSubtracted.clear();
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, 1, particles, pdgDatabase, &candidate); // currently only works for charged analyses

      tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, candidate.rho(), candidate.rhoM());
      for (auto const& trackSubtracted : tracksSubtracted) {
        particleSubTable(candidate.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), trackSubtracted.rap(), trackSubtracted.e(), 211, 0, static_cast<uint8_t>(o2::aod::mcparticle::enums::PhysicalPrimary)); // everything after phi is artificial and should not be used for analyses
      }
    }
  }

  void processCollisions(soa::Join<aod::JetCollisions, aod::BkgChargedRhos>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (skipMBGapEvents && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    inputParticles.clear();
    tracksSubtracted.clear();
    jetfindingutilities::analyseTracks<soa::Filtered<aod::JetTracks>, soa::Filtered<aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning);

    tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, collision.rho(), collision.rhoM());

    for (auto const& trackSubtracted : tracksSubtracted) {
      trackSubtractedTable(collision.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), jetderiveddatautilities::setSingleTrackSelectionBit(trackSelection));
    }
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processCollisions, "Fill table of subtracted tracks for collisions", true);

  void processMcCollisions(soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>::iterator const& mcCollision, soa::Filtered<aod::JetParticles> const& particles)
  {
    if (skipMBGapEvents && mcCollision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    inputParticles.clear();
    tracksSubtracted.clear();
    jetfindingutilities::analyseParticles<false, soa::Filtered<aod::JetParticles>, soa::Filtered<aod::JetParticles>::iterator>(inputParticles, particleSelection, 1, particles, pdgDatabase);

    tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, mcCollision.rho(), mcCollision.rhoM());

    for (auto const& trackSubtracted : tracksSubtracted) {
      particleSubtractedTable(mcCollision.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), trackSubtracted.rap(), trackSubtracted.e(), 211, 0, static_cast<uint8_t>(o2::aod::mcparticle::enums::PhysicalPrimary)); // everything after phi is artificial and should not be used for analyses
    }
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processMcCollisions, "Fill table of subtracted tracks for Mc collisions", false);

  void processD0Collisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesD0Data, aod::BkgD0Rhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedD0Table);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processD0Collisions, "Fill table of subtracted tracks for collisions with D0 candidates", false);

  void processD0McCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesD0MCP, aod::BkgD0McRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedD0Table);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processD0McCollisions, "Fill table of subtracted tracks for collisions with D0 MCP candidates", false);

  void processDplusCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesDplusData, aod::BkgDplusRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedDplusTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDplusCollisions, "Fill table of subtracted tracks for collisions with Dplus candidates", false);

  void processDplusMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesDplusMCP, aod::BkgDplusMcRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedDplusTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDplusMcCollisions, "Fill table of subtracted tracks for collisions with Dplus MCP candidates", false);

  void processDsCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesDsData, aod::BkgDsRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedDsTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDsCollisions, "Fill table of subtracted tracks for collisions with Ds candidates", false);

  void processDsMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesDsMCP, aod::BkgDsMcRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedDsTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDsMcCollisions, "Fill table of subtracted tracks for collisions with Ds MCP candidates", false);

  void processDstarCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesDstarData, aod::BkgDstarRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedDstarTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDstarCollisions, "Fill table of subtracted tracks for collisions with D* candidates", false);

  void processDstarMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesDstarMCP, aod::BkgDstarMcRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedDstarTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDstarMcCollisions, "Fill table of subtracted tracks for collisions with D* MCP candidates", false);

  void processLcCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesLcData, aod::BkgLcRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedLcTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processLcCollisions, "Fill table of subtracted tracks for collisions with Lc candidates", false);

  void processLcMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesLcMCP, aod::BkgLcMcRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedLcTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processLcMcCollisions, "Fill table of subtracted tracks for collisions with Lc MCP candidates", false);

  void processB0Collisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesB0Data, aod::BkgB0Rhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedB0Table);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processB0Collisions, "Fill table of subtracted tracks for collisions with B0 candidates", false);

  void processB0McCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesB0MCP, aod::BkgB0McRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedB0Table);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processB0McCollisions, "Fill table of subtracted tracks for collisions with B0 MCP candidates", false);

  void processBplusCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesBplusData, aod::BkgBplusRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedBplusTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processBplusCollisions, "Fill table of subtracted tracks for collisions with Bplus candidates", false);

  void processBplusMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesBplusMCP, aod::BkgBplusMcRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedBplusTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processBplusMcCollisions, "Fill table of subtracted tracks for collisions with Bplus MCP candidates", false);

  void processXicToXiPiPiCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesXicToXiPiPiData, aod::BkgXicToXiPiPiRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedXicToXiPiPiTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processXicToXiPiPiCollisions, "Fill table of subtracted tracks for collisions with XicToXiPiPi candidates", false);

  void processXicToXiPiPiMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesXicToXiPiPiMCP, aod::BkgXicToXiPiPiMcRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedXicToXiPiPiTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processXicToXiPiPiMcCollisions, "Fill table of subtracted tracks for collisions with XicToXiPiPi MCP candidates", false);

  void processDielectronCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesDielectronData, aod::BkgDielectronRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedDielectronTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDielectronCollisions, "Fill table of subtracted tracks for collisions with Dielectron candidates", false);

  void processDielectronMcCollisions(aod::JetMcCollision const&, soa::Filtered<aod::JetParticles> const& tracks, soa::Join<aod::CandidatesDielectronMCP, aod::BkgDielectronMcRhos> const& candidates)
  {
    analyseHFMc(tracks, candidates, particleSubtractedDielectronTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDielectronMcCollisions, "Fill table of subtracted tracks for collisions with Dielectron MCP candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<eventWiseConstituentSubtractorTask>(cfgc, TaskName{"subtractor-eventwiseconstituent"})}; }
