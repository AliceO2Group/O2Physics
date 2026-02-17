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

// jet finder hf-hf bar task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>

#ifndef PWGJE_JETFINDERS_JETFINDERHFHFBAR_H_
#define PWGJE_JETFINDERS_JETFINDERHFHFBAR_H_

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h> // IWYU pragma: export

#include <THn.h>
#include <TMathBase.h>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <string>
#include <vector>

template <typename CandidateTableData, typename CandidateTableMCD, typename CandidateTableMCP, typename JetTracksSubTable, typename JetParticlesSubTable, typename JetTable, typename ConstituentTable, typename JetEvtWiseSubTable, typename ConstituentEvtWiseSubTable>
struct JetFinderHFHFBarTask {
  o2::framework::Produces<JetTable> jetsTable;
  o2::framework::Produces<ConstituentTable> constituentsTable;
  o2::framework::Produces<JetEvtWiseSubTable> jetsEvtWiseSubTable;
  o2::framework::Produces<ConstituentEvtWiseSubTable> constituentsEvtWiseSubTable;

  o2::framework::HistogramRegistry registry;

  // event level configurables
  o2::framework::Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  o2::framework::Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  o2::framework::Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  o2::framework::Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  o2::framework::Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  o2::framework::Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  o2::framework::Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};
  o2::framework::Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

  // track level configurables
  o2::framework::Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  o2::framework::Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  o2::framework::Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  o2::framework::Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  o2::framework::Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  o2::framework::Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  o2::framework::Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  o2::framework::Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // cluster level configurables
  o2::framework::Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  o2::framework::Configurable<float> clusterEtaMin{"clusterEtaMin", -0.71, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  o2::framework::Configurable<float> clusterEtaMax{"clusterEtaMax", 0.71, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  o2::framework::Configurable<float> clusterPhiMin{"clusterPhiMin", 1.39, "minimum cluster phi"};
  o2::framework::Configurable<float> clusterPhiMax{"clusterPhiMax", 3.27, "maximum cluster phi"};
  o2::framework::Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  o2::framework::Configurable<float> clusterTimeMin{"clusterTimeMin", -25., "minimum Cluster time (ns)"};
  o2::framework::Configurable<float> clusterTimeMax{"clusterTimeMax", 25., "maximum Cluster time (ns)"};
  o2::framework::Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  // HF candidate level configurables
  o2::framework::Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  o2::framework::Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  o2::framework::Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  o2::framework::Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  // HF candidiate selection configurables
  o2::framework::Configurable<bool> rejectBackgroundMCDCandidates{"rejectBackgroundMCDCandidates", false, "reject background HF candidates at MC detector level"};
  o2::framework::Configurable<bool> rejectIncorrectDecaysMCP{"rejectIncorrectDecaysMCP", true, "reject HF paticles decaying to the non-analysed decay channels at MC generator level"};

  // jet level configurables
  o2::framework::Configurable<std::vector<double>> jetRadius{"jetRadius", {0.4}, "jet resolution parameters"};
  o2::framework::Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  o2::framework::Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  o2::framework::Configurable<float> jetPhiMin{"jetPhiMin", -99.0, "minimum jet phi"};
  o2::framework::Configurable<float> jetPhiMax{"jetPhiMax", 99.0, "maximum jet phi"};
  o2::framework::Configurable<float> jetEWSPtMin{"jetEWSPtMin", 0.0, "minimum event-wise subtracted jet pT"};
  o2::framework::Configurable<float> jetEWSPtMax{"jetEWSPtMax", 1000.0, "maximum event-wise subtracted jet pT"};
  o2::framework::Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  o2::framework::Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  o2::framework::Configurable<int> jetTypeParticleLevel{"jetTypeParticleLevel", 1, "Type of stored jets. 0 = full, 1 = charged, 2 = neutral"};
  o2::framework::Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  o2::framework::Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  o2::framework::Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
  o2::framework::Configurable<int> ghostRepeat{"ghostRepeat", 1, "set to 0 to gain speed if you dont need area calculation"};
  o2::framework::Configurable<bool> DoTriggering{"DoTriggering", false, "used for the charged jet trigger to remove the eta constraint on the jet axis"};
  o2::framework::Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  o2::framework::Configurable<int> jetPtBinWidth{"jetPtBinWidth", 5, "used to define the width of the jetPt bins for the THnSparse"};
  o2::framework::Configurable<bool> fillTHnSparse{"fillTHnSparse", false, "switch to fill the THnSparse"};
  o2::framework::Configurable<double> jetExtraParam{"jetExtraParam", -99.0, "sets the _extra_param in fastjet"};

  o2::framework::Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  std::vector<int> eventSelectionBits;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  std::vector<int> triggerMaskBits;

  void init(o2::framework::InitContext const&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.jetPtMax = jetPtMax;
    jetFinder.phiMin = trackPhiMin;
    jetFinder.phiMax = trackPhiMax;
    if (trackPhiMin < -98.0) {
      jetFinder.phiMin = -1.0 * M_PI;
      jetFinder.phiMax = 2.0 * M_PI;
    }
    jetFinder.jetPhiMin = jetPhiMin;
    jetFinder.jetPhiMax = jetPhiMax;
    if (jetPhiMin < -98.0) {
      jetFinder.jetPhiMin = -1.0 * M_PI;
      jetFinder.jetPhiMax = 2.0 * M_PI;
    }
    jetFinder.jetEtaMin = jetEtaMin;
    jetFinder.jetEtaMax = jetEtaMax;
    if (jetEtaMin < -98.0) {
      jetFinder.jetEtaDefault = true;
    }
    jetFinder.algorithm = static_cast<fastjet::JetAlgorithm>(static_cast<int>(jetAlgorithm));
    jetFinder.recombScheme = static_cast<fastjet::RecombinationScheme>(static_cast<int>(jetRecombScheme));
    jetFinder.ghostArea = jetGhostArea;
    jetFinder.ghostRepeatN = ghostRepeat;
    if (DoTriggering) {
      jetFinder.isTriggering = true;
    }
    jetFinder.fastjetExtraParam = jetExtraParam;

    auto jetRadiiBins = (std::vector<double>)jetRadius;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (TMath::Abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
    }
    int jetPtMaxInt = static_cast<int>(jetPtMax);
    int jetPtMinInt = static_cast<int>(jetPtMin);
    jetPtMinInt = (jetPtMinInt / jetPtBinWidth) * jetPtBinWidth;
    jetPtMaxInt = ((jetPtMaxInt + jetPtBinWidth - 1) / jetPtBinWidth) * jetPtBinWidth;
    int jetPtBinNumber = (jetPtMaxInt - jetPtMinInt) / jetPtBinWidth;
    double jetPtMinDouble = static_cast<double>(jetPtMinInt);
    double jetPtMaxDouble = static_cast<double>(jetPtMaxInt);

    registry.add("hJet", "sparse for data or mcd jets", {o2::framework::HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
    registry.add("hJetMCP", "sparse for mcp jets", {o2::framework::HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
  }

  o2::aod::EMCALClusterDefinition clusterDefinition = o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  o2::framework::expressions::Filter collisionFilter = (nabs(o2::aod::jcollision::posZ) < vertexZCut && o2::aod::jcollision::centFT0M >= centralityMin && o2::aod::jcollision::centFT0M < centralityMax && o2::aod::jcollision::trackOccupancyInTimeRange <= trackOccupancyInTimeRangeMax);
  o2::framework::expressions::Filter mcCollisionFilter = (nabs(o2::aod::jmccollision::posZ) < vertexZCut);
  o2::framework::expressions::Filter trackCuts = (o2::aod::jtrack::pt >= trackPtMin && o2::aod::jtrack::pt < trackPtMax && o2::aod::jtrack::eta >= trackEtaMin && o2::aod::jtrack::eta <= trackEtaMax && o2::aod::jtrack::phi >= trackPhiMin && o2::aod::jtrack::phi <= trackPhiMax);
  o2::framework::expressions::Filter partCuts = (o2::aod::jmcparticle::pt >= trackPtMin && o2::aod::jmcparticle::pt < trackPtMax && o2::aod::jmcparticle::eta >= trackEtaMin && o2::aod::jmcparticle::eta <= trackEtaMax && o2::aod::jmcparticle::phi >= trackPhiMin && o2::aod::jmcparticle::phi <= trackPhiMax);
  o2::framework::expressions::Filter clusterFilter = (o2::aod::jcluster::definition == static_cast<int>(clusterDefinition) && o2::aod::jcluster::eta >= clusterEtaMin && o2::aod::jcluster::eta <= clusterEtaMax && o2::aod::jcluster::phi >= clusterPhiMin && o2::aod::jcluster::phi <= clusterPhiMax && o2::aod::jcluster::energy >= clusterEnergyMin && o2::aod::jcluster::time > clusterTimeMin && o2::aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && o2::aod::jcluster::isExotic != true));
  // o2::framework::expressions::Filter candidateCuts = (o2::aod::hfcand::pt >= candPtMin && o2::aod::hfcand::pt < candPtMax && o2::aod::hfcand::y >= candYMin && o2::aod::hfcand::y < candYMax);

  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perD0Candidate = o2::aod::bkgd0::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perD0McCandidate = o2::aod::bkgd0mc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDplusCandidate = o2::aod::bkgdplus::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDplusMcCandidate = o2::aod::bkgdplusmc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDsCandidate = o2::aod::bkgds::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDsMcCandidate = o2::aod::bkgdsmc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDstarCandidate = o2::aod::bkgdstar::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDstarMcCandidate = o2::aod::bkgdstarmc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perLcCandidate = o2::aod::bkglc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perLcMcCandidate = o2::aod::bkglcmc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perB0Candidate = o2::aod::bkgb0::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perB0McCandidate = o2::aod::bkgb0mc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perBplusCandidate = o2::aod::bkgbplus::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perBplusMcCandidate = o2::aod::bkgbplusmc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perXicToXiPiPiCandidate = o2::aod::bkgxictoxipipi::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perXicToXiPiPiMcCandidate = o2::aod::bkgxictoxipipimc::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDielectronCandidate = o2::aod::bkgdielectron::candidateId;
  o2::framework::PresliceOptional<o2::soa::Filtered<JetTracksSubTable>> perDielectronMcCandidate = o2::aod::bkgdielectronmc::candidateId;

  // function that generalically processes Data and reco level events
  template <bool isEvtWiseSub, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& tracks, V const& candidate, V const& candidateBar, M& jetsTableInput, N& constituentsTableInput, O& /*originalTracks*/, float minJetPt, float maxJetPt)
  {
    if (candidate.globalIndex() == candidateBar.globalIndex() || candidate.candidateSelFlag() == candidateBar.candidateSelFlag()) {
      return;
    }
    for (auto const& track : tracks) {
      if (jetcandidateutilities::isDaughterTrack(track, candidate) && jetcandidateutilities::isDaughterTrack(track, candidateBar)) {
        return;
      }
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections) || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    inputParticles.clear();

    if constexpr (jetcandidateutilities::isCandidate<V>()) {
      if (!jetfindingutilities::analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax) || !jetfindingutilities::analyseCandidate(inputParticles, candidateBar, candPtMin, candPtMax, candYMin, candYMax)) {
        return;
      }
    }

    if constexpr (jetcandidateutilities::isMcCandidate<V>()) {
      if (!jetfindingutilities::analyseCandidateMC(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, rejectBackgroundMCDCandidates) || !jetfindingutilities::analyseCandidateMC(inputParticles, candidateBar, candPtMin, candPtMax, candYMin, candYMax, rejectBackgroundMCDCandidates)) {
        return;
      }
    }
    jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, &candidate);

    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, collision, jetsTableInput, constituentsTableInput, registry.get<THn>(HIST("hJet")), fillTHnSparse, true);
  }

  // function that generalically processes gen level events
  template <bool isEvtWiseSub, typename T, typename U, typename V, typename M, typename N>
  void analyseMCP(T const& mcCollision, U const& particles, V const& candidate, V const& candidateBar, M& jetsTableInput, N& constituentsTableInput, int jetTypeParticleLevel, float minJetPt, float maxJetPt)
  {
    if (candidate.globalIndex() == candidateBar.globalIndex() || candidate.flagMcMatchGen() == candidateBar.flagMcMatchGen()) {
      return;
    }
    for (auto const& particle : particles) {
      if (jetcandidateutilities::isDaughterParticle(candidate.template mcParticle_as<U>(), particle.globalIndex()) || jetcandidateutilities::isDaughterParticle(candidateBar.template mcParticle_as<U>(), particle.globalIndex())) {
        return;
      }
    }
    if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections)) {
      return;
    }

    if (rejectIncorrectDecaysMCP && (!jetcandidateutilities::isMatchedCandidate(candidate) || !jetcandidateutilities::isMatchedCandidate(candidateBar))) { // is this even needed in the new derived format? it means any simulations run have to force the decay channel
      return;
    }

    inputParticles.clear();
    if (!jetfindingutilities::analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax) || !jetfindingutilities::analyseCandidate(inputParticles, candidateBar, candPtMin, candPtMax, candYMin, candYMax)) {
      return;
    }
    jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, jetTypeParticleLevel, particles, pdgDatabase, &candidate);

    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, mcCollision, jetsTableInput, constituentsTableInput, registry.get<THn>(HIST("hJetMCP")), fillTHnSparse, true);
  }

  void processDummy(o2::aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetFinderHFHFBarTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(o2::soa::Filtered<o2::aod::JetCollisions>::iterator const& collision, o2::soa::Filtered<o2::aod::JetTracks> const& tracks, CandidateTableData const& candidates)
  {
    for (auto candidateIterator = candidates.begin(); candidateIterator != candidates.end(); ++candidateIterator) {
      auto candidateBarIterator = candidateIterator;
      ++candidateBarIterator;
      for (; candidateBarIterator != candidates.end(); ++candidateBarIterator) {
        typename CandidateTableData::iterator const& candidate = *candidateIterator;
        typename CandidateTableData::iterator const& candidateBar = *candidateBarIterator;
        analyseCharged<false>(collision, tracks, candidate, candidateBar, jetsTable, constituentsTable, tracks, jetPtMin, jetPtMax);
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFHFBarTask, processChargedJetsData, "charged hf jet finding on data", false);

  void processChargedJetsMCD(o2::soa::Filtered<o2::aod::JetCollisions>::iterator const& collision, o2::soa::Filtered<o2::aod::JetTracks> const& tracks, CandidateTableMCD const& candidates)
  {
    for (auto candidateIterator = candidates.begin(); candidateIterator != candidates.end(); ++candidateIterator) {
      auto candidateBarIterator = candidateIterator;
      ++candidateBarIterator;
      for (; candidateBarIterator != candidates.end(); ++candidateBarIterator) {
        typename CandidateTableMCD::iterator const& candidate = *candidateIterator;
        typename CandidateTableMCD::iterator const& candidateBar = *candidateBarIterator;
        analyseCharged<false>(collision, tracks, candidate, candidateBar, jetsTable, constituentsTable, tracks, jetPtMin, jetPtMax);
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFHFBarTask, processChargedJetsMCD, "charged hf jet finding on MC detector level", false);

  void processChargedJetsMCP(o2::soa::Filtered<o2::aod::JetMcCollisions>::iterator const& mcCollision,
                             o2::soa::Filtered<o2::aod::JetParticles> const& particles,
                             CandidateTableMCP const& candidates)
  {
    for (auto candidateIterator = candidates.begin(); candidateIterator != candidates.end(); ++candidateIterator) {
      auto candidateBarIterator = candidateIterator;
      ++candidateBarIterator;
      for (; candidateBarIterator != candidates.end(); ++candidateBarIterator) {
        typename CandidateTableMCP::iterator const& candidate = *candidateIterator;
        typename CandidateTableMCP::iterator const& candidateBar = *candidateBarIterator;
        analyseMCP<false>(mcCollision, particles, candidate, candidateBar, jetsTable, constituentsTable, 1, jetPtMin, jetPtMax);
      }
    }
  }
  PROCESS_SWITCH(JetFinderHFHFBarTask, processChargedJetsMCP, "hf jet finding on MC particle level", false);
};

#endif // PWGJE_JETFINDERS_JETFINDERHFHFBAR_H_
