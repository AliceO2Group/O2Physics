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

  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for Lc->PKPi"};
  Configurable<int> selectionFlagLcToPiPK{"selectionFlagLcToPiPK", 1, "Selection Flag for Lc->PiPK"};
  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for B+"};

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

  Preslice<aod::BkgD0Rhos> perD0Candidate = aod::bkgd0::candidateId;
  Preslice<aod::BkgLcRhos> perLcCandidate = aod::bkglc::candidateId;
  Preslice<aod::BkgBplusRhos> perBplusCandidate = aod::bkgbplus::candidateId;

  template <typename T, typename U, typename V, typename M>
  void analyseHF(T const& tracks, U const& candidates, V const& bkgRhos, M& trackSubtractedTable)
  {

    for (auto& candidate : candidates) {

      auto const bkgRhosSliced = jethfutilities::slicedPerCandidate(bkgRhos, candidate, perD0Candidate, perLcCandidate, perBplusCandidate);
      auto const bkgRho = bkgRhosSliced.iteratorAt(0);

      inputParticles.clear();
      tracksSubtracted.clear();
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, std::optional{candidate});

      tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, bkgRho.rho(), bkgRho.rhoM());
      for (auto const& trackSubtracted : tracksSubtracted) {

        trackSubtractedTable(candidate.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), trackSubtracted.E(), jetderiveddatautilities::setSingleTrackSelectionBit(trackSelection));
      }
    }
  }

  void processCollisions(soa::Join<JetCollisions, aod::BkgChargedRhos>::iterator const& collision, soa::Filtered<JetTracks> const& tracks)
  {

    inputParticles.clear();
    tracksSubtracted.clear();
    jetfindingutilities::analyseTracks<soa::Filtered<JetTracks>, soa::Filtered<JetTracks>::iterator>(inputParticles, tracks, trackSelection);

    tracksSubtracted = eventWiseConstituentSubtractor.JetBkgSubUtils::doEventConstSub(inputParticles, collision.rho(), collision.rhoM());

    for (auto const& trackSubtracted : tracksSubtracted) {
      trackSubtractedTable(collision.globalIndex(), trackSubtracted.pt(), trackSubtracted.eta(), trackSubtracted.phi(), trackSubtracted.E(), jetderiveddatautilities::setSingleTrackSelectionBit(trackSelection));
    }
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processCollisions, "Fill table of subtracted tracks for collisions", true);

  void processD0Collisions(JetCollision const& collision, aod::BkgD0Rhos const& bkgRhos, soa::Filtered<JetTracks> const& tracks, CandidatesD0Data const& candidates)
  {
    analyseHF(tracks, candidates, bkgRhos, trackSubtractedD0Table);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processD0Collisions, "Fill table of subtracted tracks for collisions with D0 candidates", false);

  void processLcCollisions(JetCollision const& collision, aod::BkgLcRhos const& bkgRhos, soa::Filtered<JetTracks> const& tracks, CandidatesLcData const& candidates)
  {
    analyseHF(tracks, candidates, bkgRhos, trackSubtractedLcTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processLcCollisions, "Fill table of subtracted tracks for collisions with Lc candidates", false);

  void processBplusCollisions(JetCollision const& collision, aod::BkgBplusRhos const& bkgRhos, soa::Filtered<JetTracks> const& tracks, CandidatesBplusData const& candidates)
  {
    analyseHF(tracks, candidates, bkgRhos, trackSubtractedBplusTable);
  }
  PROCESS_SWITCH(eventWiseConstituentSubtractorTask, processBplusCollisions, "Fill table of subtracted tracks for collisions with Bplus candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<eventWiseConstituentSubtractorTask>(cfgc, TaskName{"subtractor-eventwiseconstituent"})}; }
