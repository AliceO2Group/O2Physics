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

#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/Core/JetFindingUtilities.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename CandidateTableData, typename CandidateTableMCD, typename CandidateTableMCP, typename JetTable, typename ConstituentTable>
struct JetFinderV0Task {
  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // V0 candidate level configurables
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  Configurable<float> candPDGMass{"candPDGMass", 310, "candidate PDG for mass in clustering"};

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
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  int trackSelection = -1;
  int eventSelection = -1;
  std::string particleSelection;

  JetFinder jetFinder;
  std::vector<fastjet::PseudoJet> inputParticles;

  int candIndex;

  void init(InitContext const&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
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

    if (candPDGMass == 310) {
      candIndex = 0;
    }
    if (candPDGMass == 3122) {
      candIndex = 1;
    }
  }

  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  // Filter candidateCuts = (aod::hfcand::pt >= candPtMin && aod::hfcand::pt < candPtMax && aod::hfcand::y >= candYMin && aod::hfcand::y < candYMax);

  // function that generalically processes Data and reco level events
  template <typename T, typename U, typename V, typename M, typename N>
  void analyseCharged(T const& collision, U const& tracks, V const& candidates, M& jetsTableInput, N& constituentsTableInput, float minJetPt, float maxJetPt)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    inputParticles.clear();
    if (!jetfindingutilities::analyseV0s(inputParticles, candidates, candPtMin, candPtMax, candYMin, candYMax, candIndex)) {
      return;
    }

    /*
        if constexpr (jethfutilities::isHFMcCandidate<V>()) {
          if (!jetfindingutilities::analyseCandidateMC(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, rejectBackgroundMCDCandidates)) {
            return;
          }
        }
        */
    jetfindingutilities::analyseTracksMultipleCandidates(inputParticles, tracks, trackSelection, candidates);

    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, collision, jetsTableInput, constituentsTableInput); // will give all jets, even those without V0s
  }

  template <typename T, typename U, typename V>
  void analyseMCP(T const& collision, U const& particles, V const& candidates, int jetTypeParticleLevel, float minJetPt, float maxJetPt)
  {

    inputParticles.clear();
    if (!jetfindingutilities::analyseV0s(inputParticles, candidates, candPtMin, candPtMax, candYMin, candYMax, candIndex)) {
      return;
    }
    jetfindingutilities::analyseParticles(inputParticles, particleSelection, jetTypeParticleLevel, particles, pdgDatabase, std::optional{candidates});
    jetfindingutilities::findJets(jetFinder, inputParticles, minJetPt, maxJetPt, jetRadius, jetAreaFractionMin, collision, jetsTable, constituentsTable, true);
  }

  void processDummy(JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetFinderV0Task, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<JetTracks> const& tracks, CandidateTableData const& candidates)
  {
    analyseCharged(collision, tracks, candidates, jetsTable, constituentsTable, jetPtMin, jetPtMax);
  }
  PROCESS_SWITCH(JetFinderV0Task, processChargedJetsData, "charged hf jet finding on data", false);

  void processChargedJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<JetTracks> const& tracks, CandidateTableMCD const& candidates)
  {
    analyseCharged(collision, tracks, candidates, jetsTable, constituentsTable, jetPtMin, jetPtMax);
  }
  PROCESS_SWITCH(JetFinderV0Task, processChargedJetsMCD, "charged hf jet finding on MC detector level", false);

  void processChargedJetsMCP(JetMcCollision const& collision,
                             soa::Filtered<JetParticles> const& particles,
                             CandidateTableMCP const& candidates)
  {
    analyseMCP(collision, particles, candidates, 1, jetPtMin, jetPtMax);
  }
  PROCESS_SWITCH(JetFinderV0Task, processChargedJetsMCP, "hf jet finding on MC particle level", false);
};
