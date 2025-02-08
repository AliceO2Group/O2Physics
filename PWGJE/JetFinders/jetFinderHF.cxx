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

// jet finder hf task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Jochen Klein <jochen.klein@cern.ch>

#include <vector>

#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename CandidateTableData, typename CandidateTableMCD, typename CandidateTableMCP, typename JetTracksSubTable, typename JetParticlesSubTable, typename JetTable, typename ConstituentTable, typename JetEvtWiseSubTable, typename ConstituentEvtWiseSubTable>
struct JetFinderHFTask {
  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;
  Produces<JetEvtWiseSubTable> jetsEvtWiseSubTable;
  Produces<ConstituentEvtWiseSubTable> constituentsEvtWiseSubTable;

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<double> trackingEfficiency{"trackingEfficiency", 1.0, "tracking efficiency applied to jet finding"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // cluster level configurables
  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.71, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.71, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMin{"clusterPhiMin", 1.39, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 3.27, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.5, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -25., "minimum Cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 25., "maximum Cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  // HF candidate level configurables
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  // HF candidiate selection configurables
  Configurable<bool> rejectBackgroundMCDCandidates{"rejectBackgroundMCDCandidates", false, "reject background HF candidates at MC detector level"};
  Configurable<bool> rejectIncorrectDecaysMCP{"rejectIncorrectDecaysMCP", true, "reject HF paticles decaying to the non-analysed decay channels at MC generator level"};

  // jet level configurables
  Configurable<std::vector<double>> jetRadius{"jetRadius", {0.4}, "jet resolution parameters"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEWSPtMin{"jetEWSPtMin", 0.0, "minimum event-wise subtracted jet pT"};
  Configurable<float> jetEWSPtMax{"jetEWSPtMax", 1000.0, "maximum event-wise subtracted jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
  Configurable<int> jetTypeParticleLevel{"jetTypeParticleLevel", 1, "Type of stored jets. 0 = full, 1 = charged, 2 = neutral"};
  Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
  Configurable<int> ghostRepeat{"ghostRepeat", 1, "set to 0 to gain speed if you dont need area calculation"};
  Configurable<bool> DoTriggering{"DoTriggering", false, "used for the charged jet trigger to remove the eta constraint on the jet axis"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<int> jetPtBinWidth{"jetPtBinWidth", 5, "used to define the width of the jetPt bins for the THnSparse"};
  Configurable<bool> fillTHnSparse{"fillTHnSparse", false, "switch to fill the THnSparse"};
  Configurable<double> jetExtraParam{"jetExtraParam", -99.0, "sets the _extra_param in fastjet"};

  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  std::vector<int> eventSelectionBits;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  std::vector<int> triggerMaskBits;

  void init(InitContext const&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.jetPtMax = jetPtMax;
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
    std::vector<double> jetPtBins;
    int jetPtMaxInt = static_cast<int>(jetPtMax);
    int jetPtMinInt = static_cast<int>(jetPtMin);
    jetPtMinInt = (jetPtMinInt / jetPtBinWidth) * jetPtBinWidth;
    jetPtMaxInt = ((jetPtMaxInt + jetPtBinWidth - 1) / jetPtBinWidth) * jetPtBinWidth;
    int jetPtBinNumber = (jetPtMaxInt - jetPtMinInt) / jetPtBinWidth;
    double jetPtMinDouble = static_cast<double>(jetPtMinInt);
    double jetPtMaxDouble = static_cast<double>(jetPtMaxInt);

    registry.add("hJet", "sparse for data or mcd jets", {HistType::kTHnC, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
    registry.add("hJetMCP", "sparse for mcp jets", {HistType::kTHnC, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
  }

  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax && aod::jcollision::trackOccupancyInTimeRange <= trackOccupancyInTimeRangeMax && (!skipMBGapEvents || aod::jcollision::subGeneratorId != static_cast<int>(jetderiveddatautilities::JCollisionSubGeneratorId::mbGap)));
  Filter mcCollisionFilter = (!skipMBGapEvents || aod::jmccollision::subGeneratorId != static_cast<int>(jetderiveddatautilities::JCollisionSubGeneratorId::mbGap)); // should we add a posZ vtx cut here or leave it to analysers?
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta >= trackEtaMin && aod::jtrack::eta <= trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta >= trackEtaMin && aod::jmcparticle::eta <= trackEtaMax && aod::jmcparticle::phi >= trackPhiMin && aod::jmcparticle::phi <= trackPhiMax);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta >= clusterEtaMin && aod::jcluster::eta <= clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));
  // Filter candidateCuts = (aod::hfcand::pt >= candPtMin && aod::hfcand::pt < candPtMax && aod::hfcand::y >= candYMin && aod::hfcand::y < candYMax);

  PresliceOptional<soa::Filtered<JetTracksSubTable>> perD0Candidate = aod::bkgd0::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perD0McCandidate = aod::bkgd0mc::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perDplusCandidate = aod::bkgdplus::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perDplusMcCandidate = aod::bkgdplusmc::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perLcCandidate = aod::bkglc::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perLcMcCandidate = aod::bkglcmc::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perBplusCandidate = aod::bkgbplus::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perBplusMcCandidate = aod::bkgbplusmc::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perDielectronCandidate = aod::bkgdielectron::candidateId;
  PresliceOptional<soa::Filtered<JetTracksSubTable>> perDielectronMcCandidate = aod::bkgdielectronmc::candidateId;

  // function that generalically processes Data and reco level events
  template <bool isEvtWiseSub, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& tracks, V const& candidate, M& jetsTableInput, N& constituentsTableInput, O& /*originalTracks*/, float minJetPt, float maxJetPt)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    inputParticles.clear();

    if constexpr (jetcandidateutilities::isCandidate<V>()) {
      if (!jetfindingutilities::analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax)) {
        return;
      }
    }

    if constexpr (jetcandidateutilities::isMcCandidate<V>()) {
      if (!jetfindingutilities::analyseCandidateMC(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, rejectBackgroundMCDCandidates)) {
        return;
      }
    }
    if constexpr (isEvtWiseSub) {
      jetfindingutilities::analyseTracks<U, typename U::iterator>(inputParticles, tracks, trackSelection, trackingEfficiency);
    } else {
      jetfindingutilities::analyseTracks(inputParticles, tracks, trackSelection, trackingEfficiency, std::optional{candidate});
    }
    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, collision, jetsTableInput, constituentsTableInput, registry.get<THn>(HIST("hJet")), fillTHnSparse, true);
  }

  // function that generalically processes gen level events
  template <bool checkIsDaughter, typename T, typename U, typename V>
  void analyseMCP(T const& collision, U const& particles, V const& candidate, int jetTypeParticleLevel, float minJetPt, float maxJetPt)
  {
    if (rejectIncorrectDecaysMCP && !jetcandidateutilities::isMatchedCandidate(candidate)) { // is this even needed in the new derived format? it means any simulations run have to force the decay channel
      return;
    }

    inputParticles.clear();
    if (!jetfindingutilities::analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax)) {
      return;
    }
    if constexpr (checkIsDaughter) {
      jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, jetTypeParticleLevel, particles, pdgDatabase, std::optional{candidate});
    } else {
      jetfindingutilities::analyseParticles<false>(inputParticles, particleSelection, jetTypeParticleLevel, particles, pdgDatabase, std::optional{candidate});
    }
    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, collision, jetsTable, constituentsTable, registry.get<THn>(HIST("hJetMCP")), fillTHnSparse, true);
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetFinderHFTask, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, CandidateTableData const& candidates)
  {
    for (typename CandidateTableData::iterator const& candidate : candidates) { // why can the type not be auto?  try const auto
      analyseCharged<false>(collision, tracks, candidate, jetsTable, constituentsTable, tracks, jetPtMin, jetPtMax);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsData, "charged hf jet finding on data", false);

  void processChargedEvtWiseSubJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<JetTracksSubTable> const& tracks, CandidateTableData const& candidates)
  {
    for (typename CandidateTableData::iterator const& candidate : candidates) {
      analyseCharged<true>(collision, jetcandidateutilities::slicedPerCandidate(tracks, candidate, perD0Candidate, perDplusCandidate, perLcCandidate, perBplusCandidate, perDielectronCandidate), candidate, jetsEvtWiseSubTable, constituentsEvtWiseSubTable, tracks, jetEWSPtMin, jetEWSPtMax);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedEvtWiseSubJetsData, "charged hf jet finding on data with event-wise constituent subtraction", false);

  void processChargedJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, CandidateTableMCD const& candidates)
  {
    for (typename CandidateTableMCD::iterator const& candidate : candidates) {
      analyseCharged<false>(collision, tracks, candidate, jetsTable, constituentsTable, tracks, jetPtMin, jetPtMax);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsMCD, "charged hf jet finding on MC detector level", false);

  void processChargedEvtWiseSubJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<JetTracksSubTable> const& tracks, CandidateTableMCD const& candidates)
  {
    for (typename CandidateTableMCD::iterator const& candidate : candidates) {
      analyseCharged<true>(collision, jetcandidateutilities::slicedPerCandidate(tracks, candidate, perD0Candidate, perDplusCandidate, perLcCandidate, perBplusCandidate, perDielectronCandidate), candidate, jetsEvtWiseSubTable, constituentsEvtWiseSubTable, tracks, jetEWSPtMin, jetEWSPtMax);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedEvtWiseSubJetsMCD, "charged hf jet finding on MC detector level with event-wise constituent subtraction", false);

  void processChargedJetsMCP(soa::Filtered<aod::JetMcCollisions>::iterator const& collision,
                             soa::Filtered<aod::JetParticles> const& particles,
                             CandidateTableMCP const& candidates)
  {
    for (typename CandidateTableMCP::iterator const& candidate : candidates) {
      analyseMCP<true>(collision, particles, candidate, 1, jetPtMin, jetPtMax);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedJetsMCP, "hf jet finding on MC particle level", false);

  void processChargedEvtWiseSubJetsMCP(soa::Filtered<aod::JetMcCollisions>::iterator const& collision,
                                       soa::Filtered<JetParticlesSubTable> const& particles,
                                       CandidateTableMCP const& candidates)
  {
    for (typename CandidateTableMCP::iterator const& candidate : candidates) {
      analyseMCP<false>(collision, jetcandidateutilities::slicedPerCandidate(particles, candidate, perD0McCandidate, perDplusMcCandidate, perLcMcCandidate, perBplusMcCandidate, perDielectronMcCandidate), candidate, 1, jetPtMin, jetPtMax);
    }
  }
  PROCESS_SWITCH(JetFinderHFTask, processChargedEvtWiseSubJetsMCP, "hf jet finding on MC particle level", false);
};
