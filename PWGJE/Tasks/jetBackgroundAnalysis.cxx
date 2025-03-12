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

// jet finder QA task
//
/// \author Aimeric Landou <aimeric.landou@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <string>
#include <vector>
#include "TLorentzVector.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetBackgroundAnalysisTask {
  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range for collisions"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality for collisions"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality for collisions"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of collisions in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of collisions in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};

  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  // // If weights are implemented for MCD bkg checks, the cut on pTHatMax of jets should be introduced
  // Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  // Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};
  Configurable<float> randomConeLeadJetDeltaR{"randomConeLeadJetDeltaR", -99.0, "min distance between leading jet axis and random cone (RC) axis; if negative, min distance is set to automatic value of R_leadJet+R_RC "};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    // selection settings initialisation
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // histogram definitions

    if (doprocessRho) {
      registry.add("h2_centrality_ntracks", "; centrality; N_{tracks};", {HistType::kTH2F, {{1100, 0., 110.0}, {10000, 0.0, 10000.0}}});
      registry.add("h2_ntracks_rho", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}});
      registry.add("h2_ntracks_rhom", "; N_{tracks}; #it{rho}_{m} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {100, 0.0, 100.0}}});
      registry.add("h2_centrality_rho", "; centrality; #it{rho} (GeV/area);", {HistType::kTH2F, {{1100, 0., 110.}, {400, 0., 400.0}}});
      registry.add("h2_centrality_rhom", ";centrality; #it{rho}_{m} (GeV/area)", {HistType::kTH2F, {{1100, 0., 110.}, {100, 0., 100.0}}});
    }

    if (doprocessBkgFluctuationsData || doprocessBkgFluctuationsMCD) {
      registry.add("h2_centrality_rhorandomcone", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
      registry.add("h2_centrality_rhorandomconerandomtrackdirection", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
      registry.add("h2_centrality_rhorandomconewithoutleadingjet", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
      registry.add("h2_centrality_rhorandomconerandomtrackdirectionwithoutoneleadingjets", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
      registry.add("h2_centrality_rhorandomconerandomtrackdirectionwithouttwoleadingjets", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, {800, -400.0, 400.0}}});
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  template <typename TTracks, typename TJets>
  bool trackIsInJet(TTracks const& track, TJets const& jet)
  {
    for (auto const& constituentId : jet.tracksIds()) {
      if (constituentId == track.globalIndex()) {
        return true;
      }
    }
    return false;
  }

  template <typename TCollisions, typename TJets, typename TTracks>
  void bkgFluctuationsRandomCone(TCollisions const& collision, TJets const& jets, TTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    TRandom3 randomNumber(0);
    float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
    float randomConePhi = randomNumber.Uniform(0.0, 2 * M_PI);
    float randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
        float dEta = track.eta() - randomConeEta;
        if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomcone"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * collision.rho());

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles
    randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, 2 * M_PI) - randomConePhi, static_cast<float>(-M_PI)); // ignores actual phi of track
        float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;                                            // ignores actual eta of track
        if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirection"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * collision.rho());

    // removing the leading jet from the random cone
    if (jets.size() > 0) { // if there are no jets in the acceptance (from the jetfinder cuts) then there can be no leading jet
      float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
      float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

      bool jetWasInCone = false;
      while ((randomConeLeadJetDeltaR <= 0 && (TMath::Sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR)) || (randomConeLeadJetDeltaR > 0 && (TMath::Sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < randomConeLeadJetDeltaR))) {
        jetWasInCone = true;
        randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
        randomConePhi = randomNumber.Uniform(0.0, 2 * M_PI);
        dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
        dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
      }
      if (jetWasInCone) {
        randomConePt = 0.0;
        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
            float dEta = track.eta() - randomConeEta;
            if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconewithoutleadingjet"), collision.centrality(), randomConePt - M_PI * randomConeR * randomConeR * collision.rho());

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles, removing tracks from 2 leading jets
    double randomConePtWithoutOneLeadJet = 0;
    double randomConePtWithoutTwoLeadJet = 0;
    if (jets.size() > 1) { // if there are no jets, or just one, in the acceptance (from the jetfinder cuts) then one cannot find 2 leading jets
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, 2 * M_PI) - randomConePhi, static_cast<float>(-M_PI)); // ignores actual phi of track
          float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;                                            // ignores actual eta of track
          if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
            if (!trackIsInJet(track, jets.iteratorAt(0))) {
              randomConePtWithoutOneLeadJet += track.pt();
              if (!trackIsInJet(track, jets.iteratorAt(1))) {
                randomConePtWithoutTwoLeadJet += track.pt();
              }
            }
          }
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirectionwithoutoneleadingjets"), collision.centrality(), randomConePtWithoutOneLeadJet - M_PI * randomConeR * randomConeR * collision.rho());
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirectionwithouttwoleadingjets"), collision.centrality(), randomConePtWithoutTwoLeadJet - M_PI * randomConeR * randomConeR * collision.rho());
  }

  void processRho(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    int nTracks = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nTracks++;
      }
    }
    registry.fill(HIST("h2_centrality_ntracks"), collision.centrality(), nTracks);
    registry.fill(HIST("h2_ntracks_rho"), nTracks, collision.rho());
    registry.fill(HIST("h2_ntracks_rhom"), nTracks, collision.rhoM());
    registry.fill(HIST("h2_centrality_rho"), collision.centrality(), collision.rho());
    registry.fill(HIST("h2_centrality_rhom"), collision.centrality(), collision.rhoM());
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processRho, "QA for rho-area subtracted jets", false);

  void processBkgFluctuationsData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    bkgFluctuationsRandomCone(collision, jets, tracks);
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsData, "QA for random cone estimation of background fluctuations in data", false);

  void processBkgFluctuationsMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    bkgFluctuationsRandomCone(collision, jets, tracks);
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsMCD, "QA for random cone estimation of background fluctuations in mcd", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetBackgroundAnalysisTask>(cfgc, TaskName{"jet-background-analysis"})}; }
