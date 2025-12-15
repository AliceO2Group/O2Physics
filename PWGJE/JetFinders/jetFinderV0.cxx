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

// jet finder V0 task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <typename CandidateTableData, typename CandidateTableMCD, typename CandidateTableMCP, typename JetTable, typename ConstituentTable>
struct JetFinderV0Task {

  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "decide to run over MB gap events or not"};
  Configurable<bool> applyRCTSelections{"applyRCTSelections", true, "decide to apply RCT selections"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<bool> applyTrackingEfficiency{"applyTrackingEfficiency", {false}, "configurable to decide whether to apply artificial tracking efficiency (discarding tracks) in jet finding"};
  Configurable<std::vector<double>> trackingEfficiencyPtBinning{"trackingEfficiencyPtBinning", {0., 10, 999.}, "pt binning of tracking efficiency array if applyTrackingEfficiency is true"};
  Configurable<std::vector<double>> trackingEfficiency{"trackingEfficiency", {1.0, 1.0}, "tracking efficiency array applied to jet finding if applyTrackingEfficiency is true"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // V0 candidate level configurables
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  Configurable<float> candPDG{"candPDG", 310, "candidate PDG for mass in clustering"};

  // jet level configurables
  Configurable<std::vector<double>> jetRadius{"jetRadius", {0.4}, "jet resolution parameters"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
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
  Configurable<bool> fillTHnSparse{"fillTHnSparse", true, "switch to fill the THnSparse"};
  Configurable<double> jetExtraParam{"jetExtraParam", -99.0, "sets the _extra_param in fastjet"};
  Configurable<bool> useV0SignalFlags{"useV0SignalFlags", true, "use V0 signal flags table"};
  Configurable<bool> saveJetsWithCandidatesOnly{"saveJetsWithCandidatesOnly", true, "only save jets if they contain a V0"};

  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  std::vector<int> eventSelectionBits;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  std::vector<int> triggerMaskBits;

  int candIndex;

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

    if (candPDG == 310) {
      candIndex = 0;
    }
    if (candPDG == 3122) {
      candIndex = 1;
    }

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

    registry.add("hJet", "sparse for data or mcd jets", {HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});
    registry.add("hJetMCP", "sparse for mcp jets", {HistType::kTHnD, {{jetRadiiBins, ""}, {jetPtBinNumber, jetPtMinDouble, jetPtMaxDouble}, {40, -1.0, 1.0}, {18, 0.0, 7.0}}});

    if (applyTrackingEfficiency) {
      if (trackingEfficiencyPtBinning->size() < 2) {
        LOGP(fatal, "jetFinderV0 workflow: trackingEfficiencyPtBinning configurable should have at least two bin edges");
      }
      if (trackingEfficiency->size() + 1 != trackingEfficiencyPtBinning->size()) {
        LOGP(fatal, "jetFinderV0 workflow: trackingEfficiency configurable should have exactly one less entry than the number of bin edges set in trackingEfficiencyPtBinning configurable");
      }
    }
  }

  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax && aod::jcollision::trackOccupancyInTimeRange <= trackOccupancyInTimeRangeMax);
  Filter mcCollisionFilter = (nabs(aod::jmccollision::posZ) < vertexZCut);
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta >= trackEtaMin && aod::jtrack::eta <= trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax && aod::jmcparticle::eta >= trackEtaMin && aod::jmcparticle::eta <= trackEtaMax && aod::jmcparticle::phi >= trackPhiMin && aod::jmcparticle::phi <= trackPhiMax);

  // function that generalically processes Data and reco level events
  template <typename T, typename U, typename V, typename M, typename N>
  void analyseCharged(T const& collision, U const& tracks, V const& candidates, M& jetsTableInput, N& constituentsTableInput, float minJetPt, float maxJetPt)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRCTSelections) || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    inputParticles.clear();
    if (!jetfindingutilities::analyseV0s(inputParticles, candidates, candPtMin, candPtMax, candYMin, candYMax, candIndex, useV0SignalFlags)) {
      if (saveJetsWithCandidatesOnly) {
        return;
      }
    }

    /*
        if constexpr (jethfutilities::isHFMcCandidate<V>()) {
          if (!jetfindingutilities::analyseCandidateMC(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, rejectBackgroundMCDCandidates)) {
            return;
          }
        }
        */
    jetfindingutilities::analyseTracksMultipleCandidates(inputParticles, tracks, trackSelection, applyTrackingEfficiency, trackingEfficiency, trackingEfficiencyPtBinning, candidates);

    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, collision, jetsTableInput, constituentsTableInput, registry.get<THn>(HIST("hJet")), fillTHnSparse, saveJetsWithCandidatesOnly);
  }

  template <typename T, typename U, typename V>
  void analyseMCP(T const& mcCollision, U const& particles, V const& candidates, int jetTypeParticleLevel, float minJetPt, float maxJetPt)
  {

    if (!jetderiveddatautilities::selectMcCollision(mcCollision, skipMBGapEvents, applyRCTSelections)) {
      return;
    }

    inputParticles.clear();
    if (!jetfindingutilities::analyseV0s(inputParticles, candidates, candPtMin, candPtMax, candYMin, candYMax, candIndex, useV0SignalFlags)) {
      if (saveJetsWithCandidatesOnly) {
        return;
      }
    }
    jetfindingutilities::analyseParticles<true>(inputParticles, particleSelection, jetTypeParticleLevel, particles, pdgDatabase, &candidates);
    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, mcCollision, jetsTable, constituentsTable, registry.get<THn>(HIST("hJetMCP")), fillTHnSparse, saveJetsWithCandidatesOnly);
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetFinderV0Task, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, CandidateTableData const& candidates)
  {
    analyseCharged(collision, tracks, candidates, jetsTable, constituentsTable, jetPtMin, jetPtMax);
  }
  PROCESS_SWITCH(JetFinderV0Task, processChargedJetsData, "charged hf jet finding on data", false);

  void processChargedJetsMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, CandidateTableMCD const& candidates)
  {
    analyseCharged(collision, tracks, candidates, jetsTable, constituentsTable, jetPtMin, jetPtMax);
  }
  PROCESS_SWITCH(JetFinderV0Task, processChargedJetsMCD, "charged hf jet finding on MC detector level", false);

  void processChargedJetsMCP(soa::Filtered<aod::JetMcCollisions>::iterator const& mcCollision,
                             soa::Filtered<aod::JetParticles> const& particles,
                             CandidateTableMCP const& candidates)
  {
    analyseMCP(mcCollision, particles, candidates, 1, jetPtMin, jetPtMax);
  }
  PROCESS_SWITCH(JetFinderV0Task, processChargedJetsMCP, "hf jet finding on MC particle level", false);
};
