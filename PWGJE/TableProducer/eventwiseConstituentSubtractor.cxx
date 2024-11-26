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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct eventWiseConstituentSubtractorTask {
  Produces<aod::JTrackSubs> trackSubtractedTable;
  Produces<aod::JTrackD0Subs> trackSubtractedD0Table;
  Produces<aod::JTrackLcSubs> trackSubtractedLcTable;
  Produces<aod::JTrackBplusSubs> trackSubtractedBplusTable;
  Produces<aod::JTrackDielectronSubs> trackSubtractedDielectronTable;

  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<double> trackingEfficiency{"trackingEfficiency", 1.0, "tracking efficiency applied to jet finding"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> alpha{"alpha", 1.0, "exponent of transverse momentum in calculating the distance measure between pairs"};
  Configurable<float> rMax{"rMax", 0.24, "maximum distance of subtraction"};
  Configurable<float> eventEtaMax{"eventEtaMax", 0.9, "maximum pseudorapidity of event"};
  Configurable<bool> doRhoMassSub{"doRhoMassSub", true, "perfom mass subtraction as well"};

  JetBkgSubUtils eventWiseConstituentSubtractor;
  float bkgPhiMax_;
  std::vector<fastjet::PseudoJet> inputParticles;
  std::vector<fastjet::PseudoJet> tracksSubtracted;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    eventWiseConstituentSubtractor.setDoRhoMassSub(doRhoMassSub);
    eventWiseConstituentSubtractor.setConstSubAlphaRMax(alpha, rMax);
    eventWiseConstituentSubtractor.setMaxEtaEvent(eventEtaMax);
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);

  template <typename T, typename U, typename V>
  void analyseHF(T const& tracks, U const& candidates, V& trackSubtractedTable)
  {

    for (auto& candidate : candidates) {
      inputParticles.clear();
      tracksSubtracted.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, trackingEfficiency, std::optional{candidate});

      tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, candidate.rho(), candidate.rhoM());
      for (auto const& trackSubtracted : tracksSubtracted) {

        trackSubtractedTable(candidate.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), trackSubtracted.E(), jetderiveddatautilities::setSingleTrackSelectionBit(trackSelection));
      }
    }
  }

  void processCollisions(soa::Join<aod::JetCollisions, aod::BkgChargedRhos>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {

    inputParticles.clear();
    tracksSubtracted.clear();
    jetfindingutilities::analyseTracks<soa::Filtered<aod::JetTracks>, soa::Filtered<aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, trackingEfficiency);

    tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, collision.rho(), collision.rhoM());

    for (auto const& trackSubtracted : tracksSubtracted) {
      trackSubtractedTable(collision.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), trackSubtracted.E(), jetderiveddatautilities::setSingleTrackSelectionBit(trackSelection));
    }
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processCollisions, "Fill table of subtracted tracks for collisions", true);

  void processD0Collisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesD0Data, aod::BkgD0Rhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedD0Table);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processD0Collisions, "Fill table of subtracted tracks for collisions with D0 candidates", false);

  void processLcCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesLcData, aod::BkgLcRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedLcTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processLcCollisions, "Fill table of subtracted tracks for collisions with Lc candidates", false);

  void processBplusCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesBplusData, aod::BkgBplusRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedBplusTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processBplusCollisions, "Fill table of subtracted tracks for collisions with Bplus candidates", false);

  void processDielectronCollisions(aod::JetCollision const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Join<aod::CandidatesDielectronData, aod::BkgDielectronRhos> const& candidates)
  {
    analyseHF(tracks, candidates, trackSubtractedDielectronTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processDielectronCollisions, "Fill table of subtracted tracks for collisions with Dielectron candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<eventWiseConstituentSubtractorTask>(cfgc, TaskName{"subtractor-eventwiseconstituent"})}; }
