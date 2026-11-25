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
//
/// \author Aimeric Landou <aimeric.landou@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Yubiao Wang <yubiao.wang@cern.ch>
/// \file JetBackgroundAnalysisTask.cxx
/// \brief This file contains the implementation for the Charged Jet v2 analysis in the ALICE experiment

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TRandom3.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetBackgroundAnalysisTask {
  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range for collisions"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality for collisions"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality for collisions"};
  Configurable<bool> checkCentFT0M{"checkCentFT0M", false, "0: centFT0C as default, 1: use centFT0M estimator"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of collisions in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of collisions in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events"};
  Configurable<bool> applyMCDSelectionToMCP{"applyMCDSelectionToMCP", false, "require the first matched reconstructed collision to pass MCD event selection in MCP processes"};
  Configurable<int> nBinsFluct{"nBinsFluct", 1000, "number of bins for flucuations axes"};
  Configurable<double> deltaPtMinMaxFluct{"deltaPtMinMaxFluct", 200, "maximum and minimum delta pt value for fluctuation axes"};

  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "particle status selection, as in the jet finder"};

  Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};
  Configurable<float> randomConeLeadJetDeltaR{"randomConeLeadJetDeltaR", -99.0, "min distance between leading jet axis and random cone (RC) axis; if negative, min distance is set to automatic value of R_leadJet+R_RC "};

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  TRandom3 randomNumber{};

  o2::framework::Service<o2::framework::O2DatabasePDG> pdgDatabase;

  void init(o2::framework::InitContext&)
  {
    // Data/MCD and MCP use separate histograms and can run together.
    // Within each level, weighted and unweighted processes share histograms
    // and must not run together. Weighted processes subscribe to outlier tables.
    if (doprocessRho && doprocessRhoMCDWeighted) {
      LOGF(fatal, "Disable processRho when enabling processRhoMCDWeighted");
    }
    if (doprocessRhoMCP && doprocessRhoMCPWeighted) {
      LOGF(fatal, "Enable only one of processRhoMCP and processRhoMCPWeighted");
    }
    if (doprocessBkgFluctuationsMCDWeighted && (doprocessBkgFluctuationsData || doprocessBkgFluctuationsMCD)) {
      LOGF(fatal, "Disable processBkgFluctuationsData and processBkgFluctuationsMCD when enabling processBkgFluctuationsMCDWeighted");
    }
    if (doprocessBkgFluctuationsMCP && doprocessBkgFluctuationsMCPWeighted) {
      LOGF(fatal, "Enable only one of processBkgFluctuationsMCP and processBkgFluctuationsMCPWeighted");
    }

    // selection settings initialisation
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(trackSelections);
    randomNumber.SetSeed(0);

    // Axes definitions
    AxisSpec bkgFluctuationsAxis = {nBinsFluct, -deltaPtMinMaxFluct, deltaPtMinMaxFluct, "#delta #it{p}_{T} (GeV/#it{c})"};

    // histogram definitions

    if (doprocessRho || doprocessRhoMCDWeighted) {
      const bool doSumw2 = doprocessRhoMCDWeighted;
      registry.add("h2_centrality_ntracks", "; centrality; N_{tracks};", {HistType::kTH2F, {{1100, 0., 110.0}, {10000, 0.0, 10000.0}}}, doSumw2);
      registry.add("h2_ntracks_rho", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}}, doSumw2);
      registry.add("h2_ntracks_rhom", "; N_{tracks}; #it{rho}_{m} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {100, 0.0, 100.0}}}, doSumw2);
      registry.add("h2_centrality_rho", "; centrality; #it{rho} (GeV/area);", {HistType::kTH2F, {{1100, 0., 110.}, {400, 0., 400.0}}}, doSumw2);
      registry.add("h2_centrality_rhom", ";centrality; #it{rho}_{m} (GeV/area)", {HistType::kTH2F, {{1100, 0., 110.}, {100, 0., 100.0}}}, doSumw2);
    }
    if (doprocessRhoMCP || doprocessRhoMCPWeighted) {
      const bool doSumw2 = doprocessRhoMCPWeighted;
      registry.add("h2_nparticles_rho_mcp", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}}, doSumw2);
      registry.add("h2_nparticles_rhom_mcp", "; N_{tracks}; #it{rho}_{m} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {100, 0.0, 100.0}}}, doSumw2);
    }

    if (doprocessBkgFluctuationsData || doprocessBkgFluctuationsMCD || doprocessBkgFluctuationsMCDWeighted) {
      const bool doSumw2 = doprocessBkgFluctuationsMCDWeighted;
      registry.add("h2_centrality_rhorandomcone", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}}, doSumw2);
      registry.add("h2_centrality_rhorandomconerandomtrackdirection", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}}, doSumw2);
      registry.add("h2_centrality_rhorandomconewithoutleadingjet", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}}, doSumw2);
      registry.add("h2_centrality_rhorandomconerandomtrackdirectionwithoutoneleadingjets", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}}, doSumw2);
      registry.add("h2_centrality_rhorandomconerandomtrackdirectionwithouttwoleadingjets", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}}, doSumw2);
    }

    if (doprocessBkgFluctuationsMCP || doprocessBkgFluctuationsMCPWeighted) {
      const bool doSumw2 = doprocessBkgFluctuationsMCPWeighted;
      registry.add("h_rhorandomcone", "; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH1F, {bkgFluctuationsAxis}}, doSumw2);
      registry.add("h_rhorandomconerandomtrackdirection", "; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH1F, {bkgFluctuationsAxis}}, doSumw2);
      registry.add("h_rhorandomconewithoutleadingjet", "; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH1F, {bkgFluctuationsAxis}}, doSumw2);
      registry.add("h_rhorandomconerandomtrackdirectionwithoutoneleadingjets", "; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH1F, {bkgFluctuationsAxis}}, doSumw2);
      registry.add("h_rhorandomconerandomtrackdirectionwithouttwoleadingjets", "; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH1F, {bkgFluctuationsAxis}}, doSumw2);
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter particleCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta > trackEtaMin && aod::jmcparticle::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut);
  Filter mcEventCuts = (nabs(aod::jmccollision::posZ) < vertexZCut);

  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  Preslice<aod::JetParticles> particlesPerMCPCollision = aod::jmcparticle::mcCollisionId;
  Preslice<aod::ChargedMCParticleLevelJets> jetsPerMCPCollision = aod::jet::mcCollisionId;

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
  void bkgFluctuationsRandomCone(TCollisions const& collision, TJets const& jets, TTracks const& tracks, float centrality, float eventWeight = 1.0f)
  {
    float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
    float randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
    float randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, (-o2::constants::math::PI));
        float dEta = track.eta() - randomConeEta;
        if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomcone"), centrality, randomConePt - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), eventWeight);

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles
    float randomConePtRandomTrackDirection = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, o2::constants::math::TwoPI) - randomConePhi, (-o2::constants::math::PI)); // ignores actual phi of track
        float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;                                                               // ignores actual eta of track
        if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          randomConePtRandomTrackDirection += track.pt();
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirection"), centrality, randomConePtRandomTrackDirection - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), eventWeight);

    // removing the leading jet from the random cone
    const bool hasLead = jets.size() >= 1;
    const bool hasSub = jets.size() >= 2;
    float randomConePtWithoutLeadingJet = randomConePt;
    if (hasLead) {
      float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, (-o2::constants::math::PI));
      float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

      bool jetWasInCone = false;
      while ((randomConeLeadJetDeltaR <= 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR)) || (randomConeLeadJetDeltaR > 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < randomConeLeadJetDeltaR))) {
        jetWasInCone = true;
        randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
        randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
        dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, (-o2::constants::math::PI));
        dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
      }
      if (jetWasInCone) {
        randomConePtWithoutLeadingJet = 0.0;
        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, (-o2::constants::math::PI));
            float dEta = track.eta() - randomConeEta;
            if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
              randomConePtWithoutLeadingJet += track.pt();
            }
          }
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconewithoutleadingjet"), centrality, randomConePtWithoutLeadingJet - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), eventWeight);

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles, removing tracks from 2 leading jets
    double randomConePtWithoutOneLeadJet = randomConePtRandomTrackDirection;
    double randomConePtWithoutTwoLeadJet = randomConePtRandomTrackDirection;
    if (hasLead) {
      randomConePtWithoutOneLeadJet = 0.0;
      randomConePtWithoutTwoLeadJet = 0.0;
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, o2::constants::math::TwoPI) - randomConePhi, (-o2::constants::math::PI)); // ignores actual phi of track
          float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;                                                               // ignores actual eta of track
          if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
            const bool inLead = trackIsInJet(track, jets.iteratorAt(0));
            const bool inSub = hasSub && trackIsInJet(track, jets.iteratorAt(1));
            if (!inLead) {
              randomConePtWithoutOneLeadJet += track.pt();
              if (!hasSub || !inSub) {
                randomConePtWithoutTwoLeadJet += track.pt();
              }
            }
          }
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirectionwithoutoneleadingjets"), centrality, randomConePtWithoutOneLeadJet - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), eventWeight);
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirectionwithouttwoleadingjets"), centrality, randomConePtWithoutTwoLeadJet - o2::constants::math::PI * randomConeR * randomConeR * collision.rho(), eventWeight);
  }

  template <typename TMcCollision, typename TJets, typename TTracks>
  void bkgFluctuationsRandomConeMCP(TMcCollision const& mcCollision, TJets const& jets, TTracks const& tracks, float eventWeight = 1.0f)
  {
    float randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
    float randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
    float randomConePt = 0;
    for (auto const& track : tracks) {
      if (!jetfindingutilities::isParticleSelected(track, particleSelections, static_cast<int>(JetType::charged), pdgDatabase)) {
        continue;
      }
      float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, (-o2::constants::math::PI));
      float dEta = track.eta() - randomConeEta;
      if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
        randomConePt += track.pt();
      }
    }
    registry.fill(HIST("h_rhorandomcone"), randomConePt - o2::constants::math::PI * randomConeR * randomConeR * mcCollision.rho(), eventWeight);

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles
    float randomConePtRandomTrackDirection = 0;
    for (auto const& track : tracks) {
      if (!jetfindingutilities::isParticleSelected(track, particleSelections, static_cast<int>(JetType::charged), pdgDatabase)) {
        continue;
      }
      float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, o2::constants::math::TwoPI) - randomConePhi, (-o2::constants::math::PI)); // ignores actual phi of track
      float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;                                                               // ignores actual eta of track
      if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
        randomConePtRandomTrackDirection += track.pt();
      }
    }
    registry.fill(HIST("h_rhorandomconerandomtrackdirection"), randomConePtRandomTrackDirection - o2::constants::math::PI * randomConeR * randomConeR * mcCollision.rho(), eventWeight);

    // removing the leading jet from the random cone
    const bool hasLead = jets.size() >= 1;
    const bool hasSub = jets.size() >= 2;
    float randomConePtWithoutLeadingJet = randomConePt;
    if (hasLead) {
      float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, (-o2::constants::math::PI));
      float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

      bool jetWasInCone = false;
      while ((randomConeLeadJetDeltaR <= 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + randomConeR)) || (randomConeLeadJetDeltaR > 0 && (std::sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < randomConeLeadJetDeltaR))) {
        jetWasInCone = true;
        randomConeEta = randomNumber.Uniform(trackEtaMin + randomConeR, trackEtaMax - randomConeR);
        randomConePhi = randomNumber.Uniform(0.0, o2::constants::math::TwoPI);
        dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, (-o2::constants::math::PI));
        dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
      }
      if (jetWasInCone) {
        randomConePtWithoutLeadingJet = 0.0;
        for (auto const& track : tracks) {
          if (!jetfindingutilities::isParticleSelected(track, particleSelections, static_cast<int>(JetType::charged), pdgDatabase)) {
            continue;
          }
          float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, (-o2::constants::math::PI));
          float dEta = track.eta() - randomConeEta;
          if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
            randomConePtWithoutLeadingJet += track.pt();
          }
        }
      }
    }
    registry.fill(HIST("h_rhorandomconewithoutleadingjet"), randomConePtWithoutLeadingJet - o2::constants::math::PI * randomConeR * randomConeR * mcCollision.rho(), eventWeight);

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles, removing tracks from 2 leading jets
    double randomConePtWithoutOneLeadJet = randomConePtRandomTrackDirection;
    double randomConePtWithoutTwoLeadJet = randomConePtRandomTrackDirection;
    if (hasLead) {
      randomConePtWithoutOneLeadJet = 0.0;
      randomConePtWithoutTwoLeadJet = 0.0;
      for (auto const& track : tracks) {
        if (!jetfindingutilities::isParticleSelected(track, particleSelections, static_cast<int>(JetType::charged), pdgDatabase)) {
          continue;
        }
        float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, o2::constants::math::TwoPI) - randomConePhi, (-o2::constants::math::PI)); // ignores actual phi of track
        float dEta = randomNumber.Uniform(trackEtaMin, trackEtaMax) - randomConeEta;                                                               // ignores actual eta of track
        if (std::sqrt(dEta * dEta + dPhi * dPhi) < randomConeR) {
          const bool inLead = trackIsInJet(track, jets.iteratorAt(0));
          const bool inSub = hasSub && trackIsInJet(track, jets.iteratorAt(1));
          if (!inLead) {
            randomConePtWithoutOneLeadJet += track.pt();
            if (!hasSub || !inSub) {
              randomConePtWithoutTwoLeadJet += track.pt();
            }
          }
        }
      }
    }
    registry.fill(HIST("h_rhorandomconerandomtrackdirectionwithoutoneleadingjets"), randomConePtWithoutOneLeadJet - o2::constants::math::PI * randomConeR * randomConeR * mcCollision.rho(), eventWeight);
    registry.fill(HIST("h_rhorandomconerandomtrackdirectionwithouttwoleadingjets"), randomConePtWithoutTwoLeadJet - o2::constants::math::PI * randomConeR * randomConeR * mcCollision.rho(), eventWeight);
  }

  template <bool isMCPMatched = false, typename TCollision>
  bool isCollisionSelected(TCollision const& collision)
  {
    // Also check the vertex here for the unfiltered matched-collision table.
    if constexpr (isMCPMatched) {
      if (!(std::abs(collision.posZ()) < vertexZCut)) {
        return false;
      }
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents)) {
      return false;
    }
    if (collision.trackOccupancyInTimeRange() < trackOccupancyInTimeRangeMin || trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return false;
    }
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
    return centrality >= centralityMin && centrality <= centralityMax;
  }

  template <bool isWeighted = false, typename TMcCollision, typename TCollisions>
  bool isMcCollisionSelected(TMcCollision const& mcCollision, TCollisions const& collisions)
  {
    if (!applyMCDSelectionToMCP) {
      return true;
    }

    const auto matchedCollisions = collisions.sliceBy(collisionsPerMCPCollision, mcCollision.globalIndex());
    if (matchedCollisions.size() == 0) {
      return false;
    }

    const auto collision = matchedCollisions.begin();
    if constexpr (isWeighted) {
      if (collision.isOutlier()) {
        return false;
      }
    }
    return isCollisionSelected<true>(collision);
  }

  template <typename TCollision, typename TTracks>
  void fillRhoHistograms(TCollision const& collision, TTracks const& tracks, float eventWeight = 1.0f)
  {
    if (!isCollisionSelected(collision)) {
      return;
    }
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
    int nTracks = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nTracks++;
      }
    }
    registry.fill(HIST("h2_centrality_ntracks"), centrality, nTracks, eventWeight);
    registry.fill(HIST("h2_ntracks_rho"), nTracks, collision.rho(), eventWeight);
    registry.fill(HIST("h2_ntracks_rhom"), nTracks, collision.rhoM(), eventWeight);
    registry.fill(HIST("h2_centrality_rho"), centrality, collision.rho(), eventWeight);
    registry.fill(HIST("h2_centrality_rhom"), centrality, collision.rhoM(), eventWeight);
  }

  template <typename TMcCollision, typename TParticles>
  void fillRhoMCPHistograms(TMcCollision const& mcCollision, TParticles const& particles, float eventWeight = 1.0f)
  {
    // Optional reconstructed-collision selection is applied in the MCP processes.
    int nParticles = 0;
    for (auto const& particle : particles) {
      if (!jetfindingutilities::isParticleSelected(particle, particleSelections, static_cast<int>(JetType::charged), pdgDatabase)) {
        continue;
      }
      nParticles++;
    }
    registry.fill(HIST("h2_nparticles_rho_mcp"), nParticles, mcCollision.rho(), eventWeight);
    registry.fill(HIST("h2_nparticles_rhom_mcp"), nParticles, mcCollision.rhoM(), eventWeight);
  }

  template <typename TCollision, typename TJets, typename TTracks>
  void analyseBkgFluctuations(TCollision const& collision, TJets const& jets, TTracks const& tracks, float eventWeight = 1.0f)
  {
    if (!isCollisionSelected(collision)) {
      return;
    }
    float centrality = checkCentFT0M ? collision.centFT0M() : collision.centFT0C();
    bkgFluctuationsRandomCone(collision, jets, tracks, centrality, eventWeight);
  }

  void processRho(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    fillRhoHistograms(collision, tracks);
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processRho, "QA for rho-area subtracted jets", true);

  void processRhoMCDWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionOutliers>>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (collision.isOutlier()) {
      return;
    }
    fillRhoHistograms(collision, tracks, collision.weight());
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processRhoMCDWeighted, "Weighted rho QA at the MCD level with outlier rejection", false);

  void processRhoMCP(soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>> const& mcCollisions, soa::Filtered<aod::JetParticles> const& particles, aod::JetCollisionsMCD const& collisions)
  {
    for (auto const& mcCollision : mcCollisions) {
      if (!isMcCollisionSelected(mcCollision, collisions)) {
        continue;
      }
      const auto particlesPerCollision = particles.sliceBy(particlesPerMCPCollision, mcCollision.globalIndex());
      fillRhoMCPHistograms(mcCollision, particlesPerCollision);
    }
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processRhoMCP, "QA for rho-area subtracted jets at the MCP level", false);

  void processRhoMCPWeighted(soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>> const& mcCollisions, soa::Filtered<aod::JetParticles> const& particles, soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers> const& collisions)
  {
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isOutlier() || !isMcCollisionSelected<true>(mcCollision, collisions)) {
        continue;
      }
      const auto particlesPerCollision = particles.sliceBy(particlesPerMCPCollision, mcCollision.globalIndex());
      fillRhoMCPHistograms(mcCollision, particlesPerCollision, mcCollision.weight());
    }
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processRhoMCPWeighted, "Weighted rho QA at the MCP level with outlier rejection", false);

  void processBkgFluctuationsData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    analyseBkgFluctuations(collision, jets, tracks);
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsData, "QA for random cone estimation of background fluctuations in data", false);

  void processBkgFluctuationsMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    analyseBkgFluctuations(collision, jets, tracks);
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsMCD, "QA for random cone estimation of background fluctuations in mcd", false);

  void processBkgFluctuationsMCDWeighted(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionOutliers>>::iterator const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (collision.isOutlier()) {
      return;
    }
    analyseBkgFluctuations(collision, jets, tracks, collision.weight());
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsMCDWeighted, "Weighted random cone background fluctuations at the MCD level with outlier rejection", false);

  void processBkgFluctuationsMCP(soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>> const& mcCollisions, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets, soa::Filtered<aod::JetParticles> const& particles, aod::JetCollisionsMCD const& collisions)
  {
    for (auto const& mcCollision : mcCollisions) {
      if (!isMcCollisionSelected(mcCollision, collisions)) {
        continue;
      }
      const auto jetsPerCollision = jets.sliceBy(jetsPerMCPCollision, mcCollision.globalIndex());
      const auto particlesPerCollision = particles.sliceBy(particlesPerMCPCollision, mcCollision.globalIndex());
      bkgFluctuationsRandomConeMCP(mcCollision, jetsPerCollision, particlesPerCollision);
    }
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsMCP, "QA for random cone estimation of background fluctuations in mcp", false);

  void processBkgFluctuationsMCPWeighted(soa::Filtered<soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos, aod::JMcCollisionOutliers>> const& mcCollisions, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets, soa::Filtered<aod::JetParticles> const& particles, soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers> const& collisions)
  {
    for (auto const& mcCollision : mcCollisions) {
      if (mcCollision.isOutlier() || !isMcCollisionSelected<true>(mcCollision, collisions)) {
        continue;
      }
      const auto jetsPerCollision = jets.sliceBy(jetsPerMCPCollision, mcCollision.globalIndex());
      const auto particlesPerCollision = particles.sliceBy(particlesPerMCPCollision, mcCollision.globalIndex());
      bkgFluctuationsRandomConeMCP(mcCollision, jetsPerCollision, particlesPerCollision, mcCollision.weight());
    }
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsMCPWeighted, "Weighted random cone background fluctuations at the MCP level with outlier rejection", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetBackgroundAnalysisTask>(cfgc)};
}
