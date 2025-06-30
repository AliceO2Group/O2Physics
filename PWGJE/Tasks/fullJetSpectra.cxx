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
/// \file fullJetSpectra.cxx
/// \brief Task for full jet spectra studies in pp collisions.
/// \author Archita Rani Dash <archita.rani.dash@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Multiplicity.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include <math.h>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FullJetSpectra {

  HistogramRegistry registry;

  // MC Sample split configurables
  /*  Configurable<int> mcSplitSeed{"mcSplitSeed", 12345, "Seed for reproducible MC event splitting"};
  Configurable<float> mcClosureSplitFrac{"mcClosureSplitFrac", 0.2f, "Fraction of MC events for closure test (MCD)"};
  Configurable<bool> doMcClosure{"doMcClosure", false, "Enable random splitting for MC closure test"};
  */
  // Event configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<bool> doEMCALEventWorkaround{"doEMCALEventWorkaround", false, "apply the workaround to read the EMC trigger bit by requiring a cell content in the EMCAL"};
  Configurable<bool> doMBGapTrigger{"doMBGapTrigger", true, "set to true only when using MB-Gap Trigger JJ MC to reject MB events at the collision and track level"};
  // Configurable<bool> doMBMC{"doMBMC", false, "set to true only when using MB MC"};
  Configurable<bool> checkMcCollisionIsMatched{"checkMcCollisionIsMatched", false, "0: count whole MCcollisions, 1: select MCcollisions which only have their correspond collisions"};

  // Jet configurables
  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> jetpTMin{"jetpTMin", 20.0, "minimum jet pT"};
  Configurable<float> jetpTMax{"jetpTMax", 350., "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.3, "minimum jet eta"}; // each of these jet configurables are for the fiducial emcal cuts
  Configurable<float> jetEtaMax{"jetEtaMax", 0.3, "maximum jet eta"};  // for R = 0.4 (EMCAL eta acceptance: eta_jet = 0.7 - R)
  Configurable<float> jetPhiMin{"jetPhiMin", 1.80, "minimum jet phi"}; // phi_jet_min for R = 0.4 is 1.80
  Configurable<float> jetPhiMax{"jetPhiMax", 2.86, "maximum jet phi"}; // phi_jet_min for R = 0.4 is 2.86
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};

  // Track configurables
  Configurable<float> trackpTMin{"trackpTMin", 0.15, "minimum track pT"};
  Configurable<float> trackpTMax{"trackpTMax", 350., "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.7, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.7, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", 1.396, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 3.283, "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<std::string> eventSelections{"eventSelections", "selMCFull", "choose event selection"};
  Configurable<std::string> particleSelections{"particleSelections", "PhysicalPrimary", "set particle selections"};

  // Cluster configurables
  Configurable<std::string> clusterDefinitionS{"clusterDefinitionS", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"};
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};
  Configurable<float> clusterPhiMin{"clusterPhiMin", 1.396, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 3.283, "maximum cluster phi"};
  Configurable<float> clusterEnergyMin{"clusterEnergyMin", 0.3, "minimum cluster energy in EMCAL (GeV)"};
  Configurable<float> clusterTimeMin{"clusterTimeMin", -15., "minimum cluster time (ns)"};
  Configurable<float> clusterTimeMax{"clusterTimeMax", 15., "maximum cluster time (ns)"};
  Configurable<bool> clusterRejectExotics{"clusterRejectExotics", true, "Reject exotic clusters"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 4.0, "exponent of the event weight for the calculation of pTHeventSelectionBitsat"}; // 6 for MB MC and 4 for JJ MC
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum value of pTHat"};

  int trackSelection = -1;
  const float kJetAreaFractionMinThreshold = -98.0f;
  const float kLeadingConstituentPtMinThreshold = -98.0f;
  std::vector<int> eventSelectionBits;
  std::vector<bool> filledJetR;
  std::vector<double> jetRadiiValues;

  std::string particleSelection;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  // Random splitter instance
  /*  TRandom3 randGen;
  // float eventRandomValue = -1.0;  // default invalid
  // Cache to store random values per MC collision ID
  std::unordered_map<int64_t, float> mcCollisionRandomValues;
  */
  /*
  MC CLOSURE SPLITTING LOGIC -> still not working across different process functions. Not so trivial in O2Physics Framework!
  --------------------------
  • doMcClosure=true activates MC sample splitting.
  • Each event gets ONE random value in [0, 1), stored in eventRandomValue.
  • Events are split as:
  - ≤ mcClosureSplitFrac -> Closure (MCD)
  - > mcClosureSplitFrac -> Response (MCP + Matched)
  • This ensures mutually exclusive processing — NO double-counting.
  • eventRandomValue is reset to -1 after each event -> this is done by the `endOfEvent` defined at the end
  */

  // Add Collision Histograms' Bin Labels for clarity
  void labelCollisionHistograms(HistogramRegistry& registry)
  {
    if (doprocessDataTracks || doprocessMCTracks) {
      auto hCollisionsUnweighted = registry.get<TH1>(HIST("hCollisionsUnweighted"));
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(1, "allDetColl");
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(3, "MBRejectedDetEvents");
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(4, "EventsNotSatisfyingEventSelection");
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(5, "EMCreadoutDetEventsWithkTVXinEMC");
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(6, "AllRejectedEventsAfterEMCEventSelection");
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(7, "EMCAcceptedDetColl");
      hCollisionsUnweighted->GetXaxis()->SetBinLabel(8, "EMCAcceptedCollAfterTrackSel");
    }

    if (doprocessTracksWeighted) {
      auto hCollisionsWeighted = registry.get<TH1>(HIST("hCollisionsWeighted"));
      hCollisionsWeighted->GetXaxis()->SetBinLabel(1, "AllWeightedDetColl");
      hCollisionsWeighted->GetXaxis()->SetBinLabel(2, "WeightedCollWithVertexZ");
      hCollisionsWeighted->GetXaxis()->SetBinLabel(3, "MBRejectedDetEvents");
      hCollisionsWeighted->GetXaxis()->SetBinLabel(4, "EventsNotSatisfyingEventSelection");
      hCollisionsWeighted->GetXaxis()->SetBinLabel(5, "EMCreadoutDetJJEventsWithkTVXinEMC");
      hCollisionsWeighted->GetXaxis()->SetBinLabel(6, "AllRejectedEventsAfterEMCEventSelection");
      hCollisionsWeighted->GetXaxis()->SetBinLabel(7, "EMCAcceptedWeightedDetColl");
      hCollisionsWeighted->GetXaxis()->SetBinLabel(8, "EMCAcceptedWeightedCollAfterTrackSel");
    }

    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted) {
      auto hDetcollisionCounter = registry.get<TH1>(HIST("hDetcollisionCounter"));
      hDetcollisionCounter->GetXaxis()->SetBinLabel(1, "allDetColl");
      hDetcollisionCounter->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      hDetcollisionCounter->GetXaxis()->SetBinLabel(3, "RejectedDetCollWithOutliers");
      hDetcollisionCounter->GetXaxis()->SetBinLabel(4, "MBRejectedDetEvents");
      hDetcollisionCounter->GetXaxis()->SetBinLabel(5, "EventsNotSatisfyingEventSelection");
      hDetcollisionCounter->GetXaxis()->SetBinLabel(6, "EMCreadoutDetEventsWithkTVXinEMC");
      hDetcollisionCounter->GetXaxis()->SetBinLabel(7, "AllRejectedEventsAfterEMCEventSelection");
      hDetcollisionCounter->GetXaxis()->SetBinLabel(8, "EMCAcceptedDetColl");
    }

    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      auto hPartcollisionCounter = registry.get<TH1>(HIST("hPartcollisionCounter"));
      hPartcollisionCounter->GetXaxis()->SetBinLabel(1, "allMcColl");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(2, "McCollWithVertexZ");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(3, "PartCollWithSize>1");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(4, "RejectedPartCollForDetCollWithSize0");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(5, "RejectedPartCollWithOutliers");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(6, "MBRejectedPartEvents");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(7, "EMCreadoutDetJJEventsWithkTVXinEMC");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(8, "AllRejectedPartEventsAfterEMCEventSelection");
      hPartcollisionCounter->GetXaxis()->SetBinLabel(9, "EMCAcceptedPartColl");
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      auto hMatchedcollisionCounter = registry.get<TH1>(HIST("hMatchedcollisionCounter"));
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(1, "allDetColl");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(3, "RejectedDetCollWithOutliers");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(4, "RejectedPartCollWithOutliers");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(5, "EMCMBRejectedDetColl");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(6, "EventsNotSatisfyingEventSelection");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(7, "EMCreadoutDetJJEventsWithkTVXinEMC");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(8, "AllRejectedDetEventsAfterEMCEventSelection");
      hMatchedcollisionCounter->GetXaxis()->SetBinLabel(9, "EMCAcceptedDetColl");
    }

    if (doprocessMBCollisionsDATAWithMultiplicity || doprocessMBCollisionsWithMultiplicity || doprocessCollisionsWeightedWithMultiplicity) {
      auto hEventmultiplicityCounter = registry.get<TH1>(HIST("hEventmultiplicityCounter"));
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(1, "allDetColl");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(3, "MBRejectedDetEvents");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(4, "EventsNotSatisfyingEventSelection");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(5, "EMCreadoutDetEventsWithkTVXinEMC");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(6, "AllRejectedEventsAfterEMCEventSelection");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(7, "EMCAcceptedDetColl");
    }
  }

  // Add Bin Labels for the MC Split Event Counter
  /*  void labelMCSplitHistogram(HistogramRegistry& registry) {
  auto hSpliteventSelector = registry.get<TH1>(HIST("hSpliteventSelector"));
  hSpliteventSelector->GetXaxis()->SetBinLabel(1, "MCD");
  hSpliteventSelector->GetXaxis()->SetBinLabel(2, "MCP");
  hSpliteventSelector->GetXaxis()->SetBinLabel(3, "MatchedforRM");
}
*/
  void init(o2::framework::InitContext&)
  {
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    particleSelection = static_cast<std::string>(particleSelections);
    jetRadiiValues = (std::vector<double>)jetRadii;

    /*  if (doMcClosure) {
    // randGen.SetSeed(mcSplitSeed);
    // randGen.SetSeed(static_cast<UInt_t>(std::time(nullptr)));
    // int seed = mcSplitSeed >= 0 ? mcSplitSeed : static_cast<int>(std::time(nullptr));
    // randGen.SetSeed(seed);
    // LOGF(info, "MC closure seed = %d", seed);

    int seed = mcSplitSeed >= 0 ? mcSplitSeed : static_cast<int>(std::time(nullptr));
    randGen.SetSeed(seed);
    LOGF(info, "MC closure splitting enabled with seed = %d, split fraction = %.2f", seed, static_cast<float>(mcClosureSplitFrac));

    registry.add("hSpliteventSelector", "Random MC Split Selector;Split Type;Entries",{HistType::kTH1F, {{3, 0.0, 3.0}}}); // 0=MCD, 1=MCP, 2=RM

    //individual processes' event counters for sanity checks
    registry.add("h_MCD_splitevent_counter", "Events into MCD split", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    registry.add("h_MCP_splitevent_counter", "Events into MCP split", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    registry.add("h_Matched_splitevent_counter", "Events into Matched split", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    registry.add("hRandomValueDebug", "Random values for debugging;Random Value;Entries", {HistType::kTH1F, {{100, 0.0, 1.0}}});

    // DEBUG: Add counters for total events processed (before splitting)
    registry.add("h_MCD_total_events", "Total MCD events processed", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    registry.add("h_MCP_total_events", "Total MCP events processed", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    registry.add("h_Matched_total_events", "Total Matched events processed", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    registry.add("hMCCollisionIdDebug_MCP", "MC Collision Ids being processed", {HistType::kTH1F, {{100000, 0.0, 100000.0}}});

  }
  */
    for (std::size_t iJetRadius = 0; iJetRadius < jetRadiiValues.size(); iJetRadius++) {
      filledJetR.push_back(0.0);
    }
    auto jetRadiiBins = (std::vector<double>)jetRadii;
    if (jetRadiiBins.size() > 1) {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + (std::abs(jetRadiiBins[jetRadiiBins.size() - 1] - jetRadiiBins[jetRadiiBins.size() - 2])));
    } else {
      jetRadiiBins.push_back(jetRadiiBins[jetRadiiBins.size() - 1] + 0.1);
    }

    // Track QA histograms
    if (doprocessDataTracks || doprocessMCTracks || doprocessTracksWeighted) {
      registry.add("hCollisionsUnweighted", "event status; event status;entries", {HistType::kTH1F, {{12, 0., 12.0}}});

      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_track_energy", "track energy;Energy of tracks;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_track_energysum", "track energy sum;Sum of track energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});

      // Cluster QA histograms
      registry.add("h_clusterTime", "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}});
      registry.add("h_cluster_pt", "cluster pT;#it{p}_{T_cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_cluster_eta", "cluster #eta;#eta_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_cluster_phi", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_cluster_energy", "cluster energy;Energy of cluster;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_cluster_energysum", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});

      if (doprocessTracksWeighted) {
        registry.add("hCollisionsWeighted", "event status;event status;entries", {HistType::kTH1F, {{12, 0.0, 12.0}}});
      }
    }

    // Jet QA histograms
    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted) {

      registry.add("hDetcollisionCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.}}});

      registry.add("h_full_jet_pt", "#it{p}_{T,jet};#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_pt_pTHatcut", "#it{p}_{T,jet};#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_clusterTime", "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}});
      registry.add("h2_full_jet_nef", "#it{p}_{T,jet} vs nef at Det Level; #it{p}_{T,jet} (GeV/#it{c});nef", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}});
      registry.add("h2_full_jet_nef_rejected", "#it{p}_{T,jet} vs nef at Det Level for rejected events; #it{p}_{T,jet} (GeV/#it{c});nef", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}});

      registry.add("h_Detjet_ntracks", "#it{p}_{T,track};#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h2_full_jet_chargedconstituents", "Number of charged constituents at Det Level;#it{p}_{T,jet} (GeV/#it{c});N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h2_full_jet_neutralconstituents", "Number of neutral constituents at Det Level;#it{p}_{T,jet} (GeV/#it{c});N_{ne}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h_full_jet_chargedconstituents_pt", "track pT;#it{p}^{T,jet}_{track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_chargedconstituents_eta", "track #eta;#eta^{jet}_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_chargedconstituents_phi", "track #varphi;#varphi^{jet}_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_chargedconstituents_energy", "track energy;Energy of tracks;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_chargedconstituents_energysum", "track energy sum;Sum of track energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_neutralconstituents_pt", "cluster pT;#it{p}^{T,jet}_{cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_full_jet_neutralconstituents_eta", "cluster #eta;#eta^{jet}_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_neutralconstituents_phi", "cluster #varphi;#varphi^{jet}_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_neutralconstituents_energy", "cluster energy;Energy of cluster;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_neutralconstituents_energysum", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h2_full_jettrack_pt", "#it{p}_{T,jet} vs #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{350, 0., 350.}, {200, 0., 200.}}});
      registry.add("h2_full_jettrack_eta", "jet #eta vs jet_track #eta; #eta_{jet};#eta_{track}", {HistType::kTH2F, {{100, -1., 1.}, {500, -5., 5.}}});
      registry.add("h2_full_jettrack_phi", "jet #varphi vs jet_track #varphi; #varphi_{jet}; #varphi_{track}", {HistType::kTH2F, {{160, 0., 7.}, {160, -1., 7.}}});

      registry.add("h2_track_etaphi", "jet_track #eta vs jet_track #varphi; #eta_{track};#varphi_{track}", {HistType::kTH2F, {{500, -5., 5.}, {160, -1., 7.}}});
      registry.add("h2_jet_etaphi", "jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});
    }
    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      registry.add("hPartcollisionCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});

      registry.add("h_full_mcpjet_tablesize", "", {HistType::kTH1F, {{4, 0., 5.}}});
      registry.add("h_full_mcpjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}});
      registry.add("h_full_jet_pt_part", "jet pT;#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_eta_part", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_jet_phi_part", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h2_full_jet_nef_part", "#it{p}_{T,jet} vs nef at Part Level;#it{p}_{T,jet} (GeV/#it{c});nef", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}});

      registry.add("h_Partjet_ntracks", "#it{p}_{T,constituent};#it{p}_{T_constituent} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h2_full_jet_chargedconstituents_part", "Number of charged constituents at Part Level;#it{p}_{T,jet} (GeV/#it{c});N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h2_full_jet_neutralconstituents_part", "Number of neutral constituents at Part Level;#it{p}_{T,jet} (GeV/#it{c});N_{ne}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h_full_jet_neutralconstituents_pt_part", "#it{p}_{T} of neutral constituents at Part Level;#it{p}_{T,ne} (GeV/#it{c}); entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_neutralconstituents_eta_part", "#eta of neutral constituents at Part Level;#eta_{ne};entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_neutralconstituents_phi_part", "#varphi of neutral constituents at Part Level;#varphi_{ne};entries", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_jet_neutralconstituents_energy_part", "neutral constituents' energy;Energy of neutral constituents;entries", {HistType::kTH1F, {{400, 0., 400.}}});
      registry.add("h_full_jet_neutralconstituents_energysum_part", "neutral constituents' energy sum;Sum of neutral constituents' energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}});

      registry.add("h2_jettrack_pt_part", "#it{p}_{T,jet} vs #it{p}_{T_track}; #it{p}_{T_jet} (GeV/#it{c});#it{p}_{T_track} (GeV/#it{c})", {HistType::kTH2F, {{350, 0., 350.}, {200, 0., 200.}}});
      registry.add("h2_jettrack_eta_part", "jet #eta vs jet_track #eta; #eta_{jet};#eta_{track}", {HistType::kTH2F, {{100, -1., 1.}, {500, -5., 5.}}});
      registry.add("h2_jettrack_phi_part", "jet #varphi vs jet_track #varphi; #varphi_{jet}; #varphi_{track}", {HistType::kTH2F, {{160, 0., 7.}, {160, -1., 7.}}});

      registry.add("h2_track_etaphi_part", "jet_track #eta vs jet_track #varphi; #eta_{track};#varphi_{track}", {HistType::kTH2F, {{500, -5., 5.}, {160, -1., 7.}}});
      registry.add("h2_jet_etaphi_part", "jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});

      // registry.add("h_NOmcpemcalcollisions", "event status;entries", {HistType::kTH1F, {{100, 0., 100.}}});
      // registry.add("h_mcpemcalcollisions", "event status;entries", {HistType::kTH1F, {{100, 0., 100.}}});
      registry.add("h2_full_mcpjetOutsideFiducial_pt", "MCP jet outside EMC Fiducial Acceptance #it{p}_{T,part};#it{p}_{T,part} (GeV/c); Ncounts", {HistType::kTH2F, {{350, 0., 350.}, {10000, 0., 10000.}}});
      registry.add("h_full_mcpjetOutside_eta_part", "MCP jet #eta outside EMC Fiducial Acceptance;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_mcpjetOutside_phi_part", "MCP jet #varphi outside EMC Fiducial Acceptance;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h2_full_mcpjetInsideFiducial_pt", "MCP jet #it{p}_{T,part} inside EMC Fiducial Acceptance;#it{p}_{T,part} (GeV/c); Ncounts", {HistType::kTH2F, {{350, 0., 350.}, {10000, 0., 10000.}}});
      registry.add("h_full_mcpjetInside_eta_part", "MCP jet #eta inside EMC Fiducial Acceptance;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_mcpjetInside_phi_part", "MCP jet #varphi inside EMC Fiducial Acceptance;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      registry.add("hMatchedcollisionCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});

      registry.add("h_full_matchedmcdjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_matchedmcpjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}});
      registry.add("h_full_matchedmcdjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}});
      registry.add("h_full_matchedmcpjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}});
      registry.add("h_full_matchedmcdjet_eta", "Matched MCD jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_matchedmcdjet_phi", "Matched MCD jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_matchedmcpjet_eta", "Matched MCP jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("h_full_matchedmcpjet_phi", "Matched MCP jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}});
      registry.add("h_full_jet_deltaR", "Distance between matched Det Jet and Part Jet; #Delta R; entries", {HistType::kTH1F, {{100, 0., 1.}}});

      registry.add("h2_full_jet_energyscaleDet", "Jet Energy Scale (det); p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});

      registry.add("h2_matchedjet_etaphiDet", "Det jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});
      registry.add("h2_matchedjet_etaphiPart", "Part jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}});
      registry.add("h2_matchedjet_deltaEtaCorr", "Correlation between Det Eta and Part Eta; #eta_{jet,det}; #eta_{jet,part}", {HistType::kTH2F, {{100, -1., 1.}, {100, -1., 1.}}});
      registry.add("h2_matchedjet_deltaPhiCorr", "Correlation between Det Phi and Part Phi; #varphi_{jet,det}; #varphi_{jet,part}", {HistType::kTH2F, {{160, 0., 7.}, {160, 0., 7.}}});

      registry.add("h2_full_jet_energyscalePart", "Jet Energy Scale (part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h3_full_jet_energyscalePart", "R dependence of Jet Energy Scale (Part); #it{R}_{jet};p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_etaresolutionPart", ";p_{T,part} (GeV/c); (#eta_{jet,det} - #eta_{jet,part})/#eta_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {100, -1., 1.}}});
      registry.add("h2_full_jet_phiresolutionPart", ";p_{T,part} (GeV/c); (#varphi_{jet,det} - #varphi_{jet,part})/#varphi_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {160, -1., 7.}}});
      registry.add("h2_full_jet_energyscaleChargedPart", "Jet Energy Scale (charged part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_energyscaleNeutralPart", "Jet Energy Scale (neutral part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_energyscaleChargedVsFullPart", "Jet Energy Scale (charged part, vs. full jet pt); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_jet_energyscaleNeutralVsFullPart", "Jet Energy Scale (neutral part, vs. full jet pt); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}});
      registry.add("h2_full_fakemcdjets", "Fake MCD Jets; p_{T,det} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h2FullfakeMcpJets", "Fake MCP Jets; p_{T,part} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}});
      registry.add("h2_full_matchedmcpjet_pt", "Matched MCP jet in EMC Fiducial Acceptance #it{p}_{T,part};#it{p}_{T,part} (GeV/c); Ncounts", {HistType::kTH2F, {{350, 0., 350.}, {10000, 0., 10000.}}});

      // Response Matrix
      registry.add("h_full_jet_ResponseMatrix", "Full Jets Response Matrix; p_{T,det} (GeV/c); p_{T,part} (GeV/c)", {HistType::kTH2F, {{350, 0., 350.}, {350, 0., 350.}}});
    }

    if (doprocessCollisionsWeightedWithMultiplicity || doprocessMBCollisionsWithMultiplicity || doprocessMBCollisionsDATAWithMultiplicity) {
      registry.add("hEventmultiplicityCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
      registry.add("h_FT0Mults_occupancy", "", {HistType::kTH1F, {{3500, 0., 3500.}}});
      registry.add("h2_full_jet_FT0Amplitude", "; FT0C Amplitude; Counts", {HistType::kTH1F, {{3500, 0., 3500.}}});
      registry.add("h2_full_jet_jetpTDetVsFT0Mults", "; p_{T,det} (GeV/c); FT0C Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}});
      registry.add("h3_full_jet_jetpTDet_FT0Mults_nef", "; p_{T,det} (GeV/c); FT0C Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {3500, 0., 3500.}, {105, 0.0, 1.05}}});
    }

    // Label the histograms
    labelCollisionHistograms(registry);
    // labelMCSplitHistogram(registry);

  } // init

  // Get or generate random value for a specific MC collision
  /*  float getMCCollisionRandomValue(int64_t mcCollisionId) {
  if (!doMcClosure) return 0.0f;

  // Check if I already have a random value for this MC collision
  auto it = mcCollisionRandomValues.find(mcCollisionId);
  if (it != mcCollisionRandomValues.end()) {
  LOGF(debug, "Using cached random value %.4f for MC collision %lld", it->second, mcCollisionId);
  return it->second;
  }

  // Generate new random value for this MC collision
  float randomVal = randGen.Uniform(0.0, 1.0);
  mcCollisionRandomValues[mcCollisionId] = randomVal;

  // Debug histogram
  registry.fill(HIST("hRandomValueDebug"), randomVal);

  LOGF(info, "Generated NEW random value %.4f for MC collision %lld", randomVal, mcCollisionId);
  return randomVal;
  }
  */
  using EMCCollisionsData = o2::soa::Join<aod::JetCollisions, aod::JEMCCollisionLbs>;   // JetCollisions with EMCAL Collision Labels
  using EMCCollisionsMCD = o2::soa::Join<aod::JetCollisionsMCD, aod::JEMCCollisionLbs>; // where, JetCollisionsMCD = JetCollisions+JMcCollisionLbs

  using FullJetTableDataJoined = soa::Join<aod::FullJets, aod::FullJetConstituents>;
  using JetTableMCDJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>;
  using JetTableMCDWeightedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetEventWeights>;
  using JetTableMCPJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>;
  using JetTableMCPWeightedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetEventWeights>;

  using JetTableMCDMatchedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets>;
  using jetMcpPerMcCollision = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets>;

  using JetTableMCDMatchedWeightedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets, aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets, aod::FullMCDetectorLevelJetEventWeights>;
  using JetTableMCPMatchedWeightedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets, aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets, aod::FullMCParticleLevelJetEventWeights>;

  // Applying some cuts(filters) on collisions, tracks, clusters

  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  // Filter EMCeventCuts = (nabs(aod::collision::posZ) < vertexZCut && aod::collision::centrality >= centralityMin && aod::collision::centrality < centralityMax);
  Filter trackCuts = (aod::jtrack::pt >= trackpTMin && aod::jtrack::pt < trackpTMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));
  Preslice<jetMcpPerMcCollision> JetMCPPerMcCollision = aod::jet::mcCollisionId;
  PresliceUnsorted<soa::Filtered<aod::JetCollisionsMCD>> CollisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {

    if (jetAreaFractionMin > kJetAreaFractionMinThreshold) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    if (leadingConstituentPtMin > kLeadingConstituentPtMinThreshold) {
      bool isMinleadingConstituent = false;
      for (const auto& constituent : jet.template tracks_as<T>()) {
        if (constituent.pt() >= leadingConstituentPtMin) {
          isMinleadingConstituent = true;
          break;
        }
      }

      if (!isMinleadingConstituent) {
        return false;
      }
    }
    return true;
  }
  template <typename T>
  void fillJetHistograms(T const& jet, float weight = 1.0)
  {
    float neutralEnergy = 0.0;
    double sumtrackE = 0.0;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_full_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_full_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_full_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h2_jet_etaphi"), jet.eta(), jet.phi(), weight);

      for (const auto& cluster : jet.template clusters_as<aod::JetClusters>()) {
        registry.fill(HIST("h2_full_jet_neutralconstituents"), jet.pt(), jet.clustersIds().size(), weight);

        neutralEnergy += cluster.energy();
        double clusterpt = cluster.energy() / std::cosh(cluster.eta());
        registry.fill(HIST("h_full_jet_clusterTime"), cluster.time(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_pt"), clusterpt, weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_eta"), cluster.eta(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_phi"), cluster.phi(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_energy"), cluster.energy(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_energysum"), neutralEnergy, weight);
      }
      auto nef = neutralEnergy / jet.energy();
      registry.fill(HIST("h2_full_jet_nef"), jet.pt(), nef, weight);

      for (const auto& jettrack : jet.template tracks_as<aod::JetTracks>()) {
        sumtrackE += jettrack.energy();

        registry.fill(HIST("h_Detjet_ntracks"), jettrack.pt(), weight);
        registry.fill(HIST("h2_full_jet_chargedconstituents"), jet.pt(), jet.tracksIds().size(), weight);
        registry.fill(HIST("h2_full_jettrack_pt"), jet.pt(), jettrack.pt(), weight);
        registry.fill(HIST("h2_full_jettrack_eta"), jet.eta(), jettrack.eta(), weight);
        registry.fill(HIST("h2_full_jettrack_phi"), jet.phi(), jettrack.phi(), weight);

        registry.fill(HIST("h2_track_etaphi"), jettrack.eta(), jettrack.phi(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_pt"), jettrack.pt(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_eta"), jettrack.eta(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_phi"), jettrack.phi(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_energy"), jettrack.energy(), weight);
        registry.fill(HIST("h_full_jet_chargedconstituents_energysum"), sumtrackE, weight);
      }
    } // jet.r()
  }

  // check for nef distribution for rejected events
  template <typename T>
  void fillRejectedJetHistograms(T const& jet, float weight = 1.0)
  {
    float neutralEnergy = 0.0;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      for (const auto& cluster : jet.template clusters_as<aod::JetClusters>()) {
        neutralEnergy += cluster.energy();
      }
      auto nef = neutralEnergy / jet.energy();
      registry.fill(HIST("h2_full_jet_nef_rejected"), jet.pt(), nef, weight);
    } // jet.r()
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {
    float neutralEnergy = 0.0;
    int neutralconsts = 0;
    int chargedconsts = 0;
    int mcpjetOutsideFid = 0;
    int mcpjetInsideFid = 0;

    auto isInFiducial = [&](auto const& jet) {
      return jet.eta() >= jetEtaMin && jet.eta() <= jetEtaMax &&
             jet.phi() >= jetPhiMin && jet.phi() <= jetPhiMax;
    };

    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_full_mcpjet_tablesize"), jet.size(), weight);
      registry.fill(HIST("h_full_mcpjet_ntracks"), jet.tracksIds().size(), weight);
      registry.fill(HIST("h_full_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_full_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_full_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h2_jet_etaphi_part"), jet.eta(), jet.phi(), weight);

      if (!isInFiducial(jet)) {
        // jet is outside
        mcpjetOutsideFid++;
        registry.fill(HIST("h2_full_mcpjetOutsideFiducial_pt"), jet.pt(), mcpjetOutsideFid, weight);
        registry.fill(HIST("h_full_mcpjetOutside_eta_part"), jet.eta(), weight);
        registry.fill(HIST("h_full_mcpjetOutside_phi_part"), jet.phi(), weight);
      } else {
        // jet is inside
        mcpjetInsideFid++;
        registry.fill(HIST("h2_full_mcpjetInsideFiducial_pt"), jet.pt(), mcpjetInsideFid, weight);
        registry.fill(HIST("h_full_mcpjetInside_eta_part"), jet.eta(), weight);
        registry.fill(HIST("h_full_mcpjetInside_phi_part"), jet.phi(), weight);
      }

      for (const auto& constituent : jet.template tracks_as<aod::JetParticles>()) {
        auto pdgParticle = pdgDatabase->GetParticle(constituent.pdgCode());
        if (pdgParticle->Charge() == 0) {
          neutralconsts++;
          neutralEnergy += constituent.e(); // neutral jet constituents at particle level
          double clusterpt = constituent.e() / std::cosh(constituent.eta());
          registry.fill(HIST("h2_full_jet_neutralconstituents_part"), jet.pt(), neutralconsts, weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_pt_part"), clusterpt, weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_eta_part"), constituent.eta(), weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_phi_part"), constituent.phi(), weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_energy_part"), constituent.e(), weight);
          registry.fill(HIST("h_full_jet_neutralconstituents_energysum_part"), neutralEnergy, weight);

        } else {
          chargedconsts++;
          registry.fill(HIST("h2_full_jet_chargedconstituents_part"), jet.pt(), chargedconsts, weight); // charged jet constituents at particle level
          registry.fill(HIST("h2_jettrack_pt_part"), jet.pt(), constituent.pt(), weight);
          registry.fill(HIST("h2_jettrack_eta_part"), jet.eta(), constituent.eta(), weight);
          registry.fill(HIST("h2_jettrack_phi_part"), jet.phi(), constituent.phi(), weight);
          registry.fill(HIST("h2_track_etaphi_part"), constituent.eta(), constituent.phi(), weight);
        }
      } // constituent loop
      auto nef = neutralEnergy / jet.energy();
      registry.fill(HIST("h2_full_jet_nef_part"), jet.pt(), nef, weight);
    } // jet.r()
  }

  template <typename T, typename U>
  void fillTrackHistograms(T const& tracks, U const& clusters, float weight = 1.0)
  {
    double sumtrackE = 0.0;

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      sumtrackE += track.energy();
      registry.fill(HIST("h_track_pt"), track.pt(), weight);
      registry.fill(HIST("h_track_eta"), track.eta(), weight);
      registry.fill(HIST("h_track_phi"), track.phi(), weight);
      registry.fill(HIST("h_track_energysum"), sumtrackE, weight);
    }
    double sumclusterE = 0.0;
    for (auto const& cluster : clusters) {
      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      sumclusterE += cluster.energy();

      registry.fill(HIST("h_clusterTime"), cluster.time(), weight);
      registry.fill(HIST("h_cluster_pt"), clusterpt, weight);
      registry.fill(HIST("h_cluster_eta"), cluster.eta(), weight);
      registry.fill(HIST("h_cluster_phi"), cluster.phi(), weight);
      registry.fill(HIST("h_cluster_energy"), cluster.energy(), weight);
      registry.fill(HIST("h_cluster_energysum"), sumclusterE, weight);
    }
  }

  template <typename T, typename U>
  void fillMatchedHistograms(T const& jetBase, float weight = 1.0)
  {
    if (jetBase.has_matchedJetGeo()) { // geometrical jet matching only needed for pp - here,matching Base(Det.level) with Tag (Part. level) jets
      registry.fill(HIST("h_full_matchedmcdjet_tablesize"), jetBase.size(), weight);
      registry.fill(HIST("h_full_matchedmcdjet_ntracks"), jetBase.tracksIds().size(), weight);
      registry.fill(HIST("h2_matchedjet_etaphiDet"), jetBase.eta(), jetBase.phi(), weight);

      for (const auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        auto deltaEta = jetBase.eta() - jetTag.eta();
        auto deltaPhi = jetBase.phi() - jetTag.phi();
        auto deltaR = jetutilities::deltaR(jetBase, jetTag);

        registry.fill(HIST("h_full_jet_deltaR"), deltaR, weight);
        registry.fill(HIST("h_full_matchedmcpjet_tablesize"), jetTag.size(), weight);
        registry.fill(HIST("h_full_matchedmcpjet_ntracks"), jetTag.tracksIds().size(), weight);
        registry.fill(HIST("h2_matchedjet_etaphiPart"), jetTag.eta(), jetTag.phi(), weight);
        registry.fill(HIST("h2_matchedjet_deltaEtaCorr"), jetBase.eta(), jetTag.eta(), weight);
        registry.fill(HIST("h2_matchedjet_deltaPhiCorr"), jetBase.phi(), jetTag.phi(), weight);

        // JES for fulljets
        registry.fill(HIST("h2_full_jet_energyscaleDet"), jetBase.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h2_full_jet_energyscalePart"), jetTag.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_full_jet_energyscalePart"), jetBase.r() / 100.0, jetTag.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h2_full_jet_etaresolutionPart"), jetTag.pt(), deltaEta / jetTag.eta(), weight);
        registry.fill(HIST("h2_full_jet_phiresolutionPart"), jetTag.pt(), deltaPhi / jetTag.phi(), weight);

        // Response Matrix
        registry.fill(HIST("h_full_jet_ResponseMatrix"), jetBase.pt(), jetTag.pt(), weight); // MCD vs MCP jet pT
      } // jetTag
    } // jetBase
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(FullJetSpectra, processDummy, "dummy task", true);

  void processJetsData(soa::Filtered<EMCCollisionsData>::iterator const& collision, FullJetTableDataJoined const& jets, aod::JetTracks const&, aod::JetClusters const&)
  {
    bool eventAccepted = false;

    registry.fill(HIST("hDetcollisionCounter"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hDetcollisionCounter"), 1.5); // DetCollWithVertexZ

    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hDetcollisionCounter"), 2.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hDetcollisionCounter"), 3.5); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hDetcollisionCounter"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        registry.fill(HIST("hDetcollisionCounter"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
        eventAccepted = true;
      }
    }

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) || !isAcceptedJet<aod::JetTracks>(jet)) {
          fillRejectedJetHistograms(jet, 1.0);
        }
      }
      registry.fill(HIST("hDetcollisionCounter"), 5.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hDetcollisionCounter"), 6.5); // EMCAcceptedDetColl

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsData, "Full Jets Data", false);

  void processJetsMCD(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDJoined const& jets, aod::JetTracks const&, aod::JetClusters const&)
  {
    bool eventAccepted = false;
    double weight = 1.0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    /*    if (doMcClosure) {
    // Count total events processed (before splitting decision)
    registry.fill(HIST("h_MCD_total_events"), 0.5);

    // DEBUG: Let's verify what collision IDs we're actually seeing
    LOGF(info, "[MCD DEBUG] Processing MC collision ID: %lld", collision.mcCollisionId());

    // Get random value for this MC collision
    float eventRandomValue = getMCCollisionRandomValue(collision.mcCollisionId());

    // MCD gets events with random value <= split fraction (20%)
    if (eventRandomValue > mcClosureSplitFrac) {
    LOGF(debug, "[MCD] Event REJECTED: rand = %.4f > split = %.2f (MC collision %d)",
    eventRandomValue, static_cast<float>(mcClosureSplitFrac), collision.mcCollisionId());
    return; // This event goes to MCP & Matched processes
  }

  LOGF(info, "[MCD] Event ACCEPTED: rand = %.4f <= split = %.2f (MC collision %d)",
  eventRandomValue, static_cast<float>(mcClosureSplitFrac), collision.mcCollisionId());

  registry.fill(HIST("hSpliteventSelector"), 0.5); // 20% Closure input for the measured spectra (reco)
  registry.fill(HIST("h_MCD_splitevent_counter"), 0.5);
  }
  */
    registry.fill(HIST("hDetcollisionCounter"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hDetcollisionCounter"), 1.5); // DetCollWithVertexZ

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) { // for MCD jets only to remove outliers; setting pTHatMaxMCD = 1 improves purity
        registry.fill(HIST("hDetcollisionCounter"), 2.5);               // RejectedDetCollWithOutliers
        return;
      }
      // this cut only to be used for calculating Jet Purity and not for Response Matrix
      // this is mainly applied to remove all high weight jets causing big fluctuations
      if (jet.pt() > 1 * pTHat) {
        registry.fill(HIST("h_full_jet_pt_pTHatcut"), jet.pt(), weight);
      }
    }
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hDetcollisionCounter"), 3.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hDetcollisionCounter"), 4.5); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hDetcollisionCounter"), 5.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hDetcollisionCounter"), 5.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) || !isAcceptedJet<aod::JetTracks>(jet)) {
          fillRejectedJetHistograms(jet, 1.0);
        }
      }
      registry.fill(HIST("hDetcollisionCounter"), 6.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hDetcollisionCounter"), 7.5); // EMCAcceptedDetColl

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCD, "Full Jets at Detector Level", false);

  void processJetsMCDWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDWeightedJoined const& jets, aod::JMcCollisions const&, aod::JetTracks const&, aod::JetClusters const&)
  {
    bool eventAccepted = false;
    double pTHat = 10. / (std::pow(collision.mcCollision().weight(), 1.0 / pTHatExponent));

    /*  if (doMcClosure) {
    // Count total events processed (before splitting decision)
    registry.fill(HIST("h_MCD_total_events"), 0.5);

    // DEBUG: Let's verify what collision IDs we're actually seeing
    LOGF(info, "[MCD DEBUG] Processing MC collision ID: %lld", collision.mcCollisionId());

    // Get random value for this MC collision
    float eventRandomValue = getMCCollisionRandomValue(collision.mcCollisionId());

    // MCD gets events with random value <= split fraction (20%)
    if (eventRandomValue > mcClosureSplitFrac) {
    LOGF(debug, "[MCD] Event REJECTED: rand = %.4f > split = %.2f (MC collision %d)",
    eventRandomValue, static_cast<float>(mcClosureSplitFrac), collision.mcCollisionId());
    return; // This event goes to MCP & Matched processes
  }

  LOGF(info, "[MCD] Event ACCEPTED: rand = %.4f <= split = %.2f (MC collision %d)",
  eventRandomValue, static_cast<float>(mcClosureSplitFrac), collision.mcCollisionId());

  registry.fill(HIST("hSpliteventSelector"), 0.5); // 20% Closure input for the measured spectra (reco)
  registry.fill(HIST("h_MCD_splitevent_counter"), 0.5);
  }
  */
    registry.fill(HIST("hDetcollisionCounter"), 0.5, collision.mcCollision().weight()); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hDetcollisionCounter"), 1.5, collision.mcCollision().weight()); // DetCollWithVertexZ
    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {                     // for MCD jets only to remove outliers; setting pTHatMaxMCD = 1 improves purity
        registry.fill(HIST("hDetcollisionCounter"), 2.5, collision.mcCollision().weight()); // RejectedDetCollWithOutliers
        return;
      }
      // this cut only to be used for calculating Jet Purity and not for Response Matrix
      // this is mainly applied to remove all high weight jets causing big fluctuations
      if (jet.pt() > 1 * pTHat) {
        registry.fill(HIST("h_full_jet_pt_pTHatcut"), jet.pt(), collision.mcCollision().weight());
      }
    }
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hDetcollisionCounter"), 3.5, collision.mcCollision().weight()); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hDetcollisionCounter"), 4.5, collision.mcCollision().weight()); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hDetcollisionCounter"), 5.5, collision.mcCollision().weight()); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hDetcollisionCounter"), 5.5, collision.mcCollision().weight()); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) || !isAcceptedJet<aod::JetTracks>(jet)) {
          fillRejectedJetHistograms(jet, collision.mcCollision().weight());
        }
      }
      registry.fill(HIST("hDetcollisionCounter"), 6.5, collision.mcCollision().weight()); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hDetcollisionCounter"), 7.5, collision.mcCollision().weight()); // EMCAcceptedDetColl

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      fillJetHistograms(jet, collision.mcCollision().weight());
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCDWeighted, "Full Jets at Detector Level on weighted events", false);

  void processJetsMCP(aod::JetMcCollision const& mccollision, JetTableMCPJoined const& jets, aod::JetParticles const&, soa::SmallGroups<EMCCollisionsMCD> const& collisions)
  {
    bool eventAccepted = false;
    double weight = 1.0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    /*  if (doMcClosure) {
    // Count total events processed (before splitting decision)
    registry.fill(HIST("h_MCP_total_events"), 0.5);

    // DEBUG: Let's verify what collision IDs we're actually seeing
    LOGF(info, "[MCP DEBUG] Processing MC collision ID: %lld", mccollision.globalIndex());

    // Get random value for this MC collision
    float eventRandomValue = getMCCollisionRandomValue(mccollision.globalIndex());

    // DEBUG: Track which MC collisions we're processing
    registry.fill(HIST("hMCCollisionIdDebug_MCP"), static_cast<float>(mccollision.globalIndex() % 100000));

    // MCP gets events with random value > split fraction (80%)
    if (eventRandomValue <= mcClosureSplitFrac) {
    LOGF(debug, "[MCP] Event REJECTED: rand = %.4f <= split = %.2f (MC collision %lld)",
    eventRandomValue, static_cast<float>(mcClosureSplitFrac), mccollision.globalIndex());
    return; // This event goes to MCD only
  }

  LOGF(info, "[MCP] Event ACCEPTED: rand = %.4f > split = %.2f (MC collision %lld)",
  eventRandomValue, static_cast<float>(mcClosureSplitFrac), mccollision.globalIndex());

  registry.fill(HIST("hSpliteventSelector"), 1.5); // remaining 80% input for MCP
  registry.fill(HIST("h_MCP_splitevent_counter"), 0.5);
  }
  */
    registry.fill(HIST("hPartcollisionCounter"), 0.5); // allMcColl
    if (std::fabs(mccollision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hPartcollisionCounter"), 1.5); // McCollWithVertexZ
    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("hPartcollisionCounter"), 2.5); // PartCollWithSize>1

    if (collisions.size() == 0) {
      registry.fill(HIST("hPartcollisionCounter"), 3.5); // RejectedPartCollForDetCollWithSize0
      return;
    }

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hPartcollisionCounter"), 4.5); // RejectedPartCollWithOutliers
        return;
      }
    }

    if (doMBGapTrigger && mccollision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      // Fill rejected MB events;
      registry.fill(HIST("hPartcollisionCounter"), 5.5); // MBRejectedPartEvents
      return;
    }

    for (auto const& collision : collisions) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        return;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          eventAccepted = true;
          if (collision.alias_bit(kTVXinEMC)) {
            registry.fill(HIST("hPartcollisionCounter"), 6.5); // EMCreadoutDetEventsWithkTVXinEMC
          }
        }
      } else {
        if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
          registry.fill(HIST("hPartcollisionCounter"), 6.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    }
    if (!eventAccepted) {
      registry.fill(HIST("hPartcollisionCounter"), 7.5); // AllRejectedPartEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hPartcollisionCounter"), 8.5); // EMCAcceptedWeightedPartColl

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        continue;
      }
      if (checkMcCollisionIsMatched) { // basically checks if the same collisions are generated at the Part level as those at the Det level
        auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, jet.mcCollisionId());
        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits)) {
          // Now here for every matched collision, I fill the corresponding jet histograms.
          fillMCPHistograms(jet);
        }
      } else {
        fillMCPHistograms(jet);
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCP, "Full Jets at Particle Level", false);

  void processJetsMCPWeighted(aod::JetMcCollision const& mccollision, JetTableMCPWeightedJoined const& jets, aod::JetParticles const&, soa::SmallGroups<EMCCollisionsMCD> const& collisions)
  {
    bool eventAccepted = false;
    float pTHat = 10. / (std::pow(mccollision.weight(), 1.0 / pTHatExponent));

    /*    if (doMcClosure) {
    // Count total events processed (before splitting decision)
    registry.fill(HIST("h_MCP_total_events"), 0.5);

    // DEBUG: Let's verify what collision IDs we're actually seeing
    LOGF(info, "[MCP DEBUG] Processing MC collision ID: %lld", mccollision.globalIndex());

    // Get random value for this MC collision
    float eventRandomValue = getMCCollisionRandomValue(mccollision.globalIndex());

    // DEBUG: Track which MC collisions we're processing
    registry.fill(HIST("hMCCollisionIdDebug_MCP"), static_cast<float>(mccollision.globalIndex() % 100000));

    // MCP gets events with random value > split fraction (80%)
    if (eventRandomValue <= mcClosureSplitFrac) {
    LOGF(debug, "[MCP] Event REJECTED: rand = %.4f <= split = %.2f (MC collision %lld)",
    eventRandomValue, static_cast<float>(mcClosureSplitFrac), mccollision.globalIndex());
    return; // This event goes to MCD only
  }

  LOGF(info, "[MCP] Event ACCEPTED: rand = %.4f > split = %.2f (MC collision %lld)",
  eventRandomValue, static_cast<float>(mcClosureSplitFrac), mccollision.globalIndex());

  registry.fill(HIST("hSpliteventSelector"), 1.5); // remaining 80% input for MCP
  registry.fill(HIST("h_MCP_splitevent_counter"), 0.5);
  }
  */
    registry.fill(HIST("hPartcollisionCounter"), 0.5, mccollision.weight()); // allMcColl
    if (std::fabs(mccollision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hPartcollisionCounter"), 1.5, mccollision.weight()); // McCollWithVertexZ
    if (collisions.size() < 1) {
      return;
    }
    registry.fill(HIST("hPartcollisionCounter"), 2.5, mccollision.weight()); // PartCollWithSize>1

    if (collisions.size() == 0) {
      registry.fill(HIST("hPartcollisionCounter"), 3.5, mccollision.weight()); // RejectedPartCollForDetCollWithSize0
      return;
    }

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hPartcollisionCounter"), 4.5, mccollision.weight()); // RejectedPartCollWithOutliers
        return;
      }
    }

    if (doMBGapTrigger && mccollision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      // Fill rejected MB events
      registry.fill(HIST("hPartcollisionCounter"), 5.5, mccollision.weight()); // MBRejectedPartEvents
      return;
    }

    for (auto const& collision : collisions) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        return;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          eventAccepted = true;
          if (collision.alias_bit(kTVXinEMC)) {
            registry.fill(HIST("hPartcollisionCounter"), 6.5, mccollision.weight()); // EMCreadoutDetJJEventsWithkTVXinEMC
          }
        }
      } else {
        if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
          registry.fill(HIST("hPartcollisionCounter"), 6.5, mccollision.weight()); // EMCreadoutDetJJEventsWithkTVXinEMC
        }
      }
    }
    if (!eventAccepted) {
      registry.fill(HIST("hPartcollisionCounter"), 7.5, mccollision.weight()); // AllRejectedPartEventsAfterEMCEventSelection
      return;
    }
    // Fill EMCAL JJ Part events
    registry.fill(HIST("hPartcollisionCounter"), 8.5, mccollision.weight()); // EMCAcceptedWeightedPartColl

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        return;
      }
      if (!isAcceptedJet<aod::JetParticles>(jet)) {
        return;
      }
      if (doMBGapTrigger && jet.eventWeight() == 1) {
        return;
      }

      if (checkMcCollisionIsMatched) {
        auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, jet.mcCollisionId());

        if (collisionspermcpjet.size() >= 1 && jetderiveddatautilities::selectCollision(collisionspermcpjet.begin(), eventSelectionBits)) {
          fillMCPHistograms(jet, jet.eventWeight());
        }
      } else {
        fillMCPHistograms(jet, jet.eventWeight());
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCPWeighted, "Full Jets at Particle Level on weighted events", false);

  void processJetsMCPMCDMatched(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDMatchedJoined const& mcdjets, jetMcpPerMcCollision const& mcpjets, aod::JMcCollisions const&, aod::JetTracks const&, aod::JetClusters const&, aod::JetParticles const&)
  {
    bool eventAccepted = false;
    int fakeMcdJet = 0;
    int fakeMcpJet = 0;
    double weight = 1.0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    const auto mcpJetsPerMcCollision = mcpjets.sliceBy(JetMCPPerMcCollision, collision.mcCollisionId());

    /*    if (doMcClosure) {
    // Count total events processed (before splitting decision)
    registry.fill(HIST("h_Matched_total_events"), 0.5);

    // Use consistent MC collision ID - same as MCD
    int64_t mcCollisionId = collision.mcCollisionId();

    // DEBUG: Let's verify what collision IDs we're actually seeing
    LOGF(info, "[Matched DEBUG] Processing MC collision ID: %lld", mcCollisionId);
    float eventRandomValue = getMCCollisionRandomValue(mcCollisionId);

    // Matched gets events with random value > split fraction (80%) - same as MCP
    if (eventRandomValue <= mcClosureSplitFrac) {
    LOGF(debug, "[Matched] Event REJECTED: rand = %.4f <= split = %.2f (MC collision %lld)",
    eventRandomValue, static_cast<float>(mcClosureSplitFrac), mcCollisionId);
    return; // This event goes to MCD only
  }

  LOGF(info, "[Matched] Event ACCEPTED: rand = %.4f > split = %.2f (MC collision %lld)",
  eventRandomValue, static_cast<float>(mcClosureSplitFrac), mcCollisionId);

  registry.fill(HIST("hSpliteventSelector"), 2.5); // Bin for Response Matrix
  registry.fill(HIST("h_Matched_splitevent_counter"), 0.5);
  }
  */
    registry.fill(HIST("hMatchedcollisionCounter"), 0.5); // allDetColl

    if (std::fabs(collision.posZ()) > vertexZCut) { // making double sure this condition is satisfied
      return;
    }
    registry.fill(HIST("hMatchedcollisionCounter"), 1.5); // DetCollWithVertexZ

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& mcdjet : mcdjets) {
      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hMatchedcollisionCounter"), 2.5); // RejectedDetCollWithOutliers
        return;
      }
    }
    // //outlier check for Part collisions: commenting out this for now otherwise this rejects all Det Colls
    // for (auto const& mcpjet : mcpjets) {
    //   if (mcpjet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
    //     registry.fill(HIST("hMatchedcollisionCounter"),3.5); //RejectedPartCollWithOutliers
    //     return;
    //   }
    // }

    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hMatchedcollisionCounter"), 4.5); // EMCMBRejectedDetColl
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hMatchedcollisionCounter"), 5.5); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hMatchedcollisionCounter"), 6.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hMatchedcollisionCounter"), 6.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }
    if (!eventAccepted) {
      registry.fill(HIST("hMatchedcollisionCounter"), 7.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hMatchedcollisionCounter"), 8.5); // EMCAcceptedDetColl

    for (const auto& mcdjet : mcdjets) {
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      // Check if MCD jet is within the EMCAL fiducial region; if not then flag it as a fake jet
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax || mcdjet.eta() < jetEtaMin || mcdjet.eta() > jetEtaMax) {
        fakeMcdJet++;
        registry.fill(HIST("h2_full_fakemcdjets"), mcdjet.pt(), fakeMcdJet, 1.0);
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<jetMcpPerMcCollision>()) {
        // apply emcal fiducial cuts to the matched particle level jets
        if (mcpjet.eta() > jetEtaMax || mcpjet.eta() < jetEtaMin || mcpjet.phi() > jetPhiMax || mcpjet.phi() < jetPhiMin) {
          fakeMcpJet++;
          registry.fill(HIST("h2_full_fakemcpjets"), mcpjet.pt(), fakeMcpJet, 1.0);
          continue;
        }
      } // mcpjet loop
      fillMatchedHistograms<JetTableMCDMatchedJoined::iterator, jetMcpPerMcCollision>(mcdjet);
    } // mcdjet loop
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCPMCDMatched, "Full Jet finder MCP matched to MCD", false);

  void processJetsMCPMCDMatchedWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDMatchedWeightedJoined const& mcdjets, JetTableMCPMatchedWeightedJoined const& mcpjets, aod::JMcCollisions const&, aod::JetTracks const&, aod::JetClusters const&, aod::JetParticles const&)
  {
    bool eventAccepted = false;
    int fakeMcdJet = 0;
    int fakeMcpJet = 0;
    int NPartJetFid = 0;
    float eventWeight = collision.mcCollision().weight();
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    const auto mcpJetsPerMcCollision = mcpjets.sliceBy(JetMCPPerMcCollision, collision.mcCollisionId());

    /*    if (doMcClosure) {
    // Count total events processed (before splitting decision)
    registry.fill(HIST("h_Matched_total_events"), 0.5);

    // Use consistent MC collision ID - same as MCD
    int64_t mcCollisionId = collision.mcCollisionId();

    // DEBUG: Let's verify what collision IDs we're actually seeing
    LOGF(info, "[Matched DEBUG] Processing MC collision ID: %lld", mcCollisionId);
    float eventRandomValue = getMCCollisionRandomValue(mcCollisionId);

    // Matched gets events with random value > split fraction (80%) - same as MCP
    if (eventRandomValue <= mcClosureSplitFrac) {
    LOGF(debug, "[Matched] Event REJECTED: rand = %.4f <= split = %.2f (MC collision %lld)",
    eventRandomValue, static_cast<float>(mcClosureSplitFrac), mcCollisionId);
    return; // This event goes to MCD only
  }

  LOGF(info, "[Matched] Event ACCEPTED: rand = %.4f > split = %.2f (MC collision %lld)",
  eventRandomValue, static_cast<float>(mcClosureSplitFrac), mcCollisionId);

  registry.fill(HIST("hSpliteventSelector"), 2.5); // Bin for Response Matrix
  registry.fill(HIST("h_Matched_splitevent_counter"), 0.5);
  }
  */
    registry.fill(HIST("hMatchedcollisionCounter"), 0.5, eventWeight); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {                    // making double sure this condition is satisfied
      return;
    }
    registry.fill(HIST("hMatchedcollisionCounter"), 1.5, eventWeight); // DetCollWithVertexZ

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& mcdjet : mcdjets) {
      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hMatchedcollisionCounter"), 2.5, eventWeight); // RejectedDetCollWithOutliers
        return;
      }
    }
    // outlier check for Part collisions: commenting out this for now otherwise this rejects all Det Colls
    //  for (auto const& mcpjet : mcpjets) {
    //    if (mcpjet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
    //      registry.fill(HIST("hMatchedcollisionCounter"),3.5, eventWeight); //RejectedPartCollWithOutliers
    //      return;
    //    }
    //  }

    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hMatchedcollisionCounter"), 4.5, eventWeight); // EMCMBRejectedDetColl
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hMatchedcollisionCounter"), 5.5, eventWeight); // EventsNotSatisfyingEventSelection
      return;
    }

    for (auto const& mcpjet : mcpJetsPerMcCollision) {
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) { // outlier rejection for MCP: Should I remove this cut as I'm already doing MC outlier rejection @L1071?
        return;
      }
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hMatchedcollisionCounter"), 6.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hMatchedcollisionCounter"), 6.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
      }
    }
    if (!eventAccepted) {
      registry.fill(HIST("hMatchedcollisionCounter"), 7.5, eventWeight); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hMatchedcollisionCounter"), 8.5, eventWeight); // EMCAcceptedDetColl

    for (const auto& mcdjet : mcdjets) {
      // Check if MCD jet is within the EMCAL fiducial region; if not then flag it as a fake jet
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax || mcdjet.eta() < jetEtaMin || mcdjet.eta() > jetEtaMax) {
        fakeMcdJet++;
        registry.fill(HIST("h2_full_fakemcdjets"), mcdjet.pt(), fakeMcdJet, eventWeight);
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }

      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedWeightedJoined>()) {
        // apply emcal fiducial cuts to the matched particle level jets - if the matched mcp jet lies outside of the EMCAL fiducial, flag it as a fake jet
        if (mcpjet.eta() > jetEtaMax || mcpjet.eta() < jetEtaMin || mcpjet.phi() > jetPhiMax || mcpjet.phi() < jetPhiMin) {
          fakeMcpJet++;
          registry.fill(HIST("h2FullfakeMcpJets"), mcpjet.pt(), fakeMcpJet, eventWeight);
          continue;
        } else {
          NPartJetFid++;
          // // If both MCD-MCP matched jet pairs are within the EMCAL fiducial region, fill these histos
          registry.fill(HIST("h2_full_matchedmcpjet_pt"), mcpjet.pt(), NPartJetFid, eventWeight);
          registry.fill(HIST("h_full_matchedmcpjet_eta"), mcpjet.eta(), eventWeight);
          registry.fill(HIST("h_full_matchedmcpjet_phi"), mcpjet.phi(), eventWeight);
        }
      } // mcpjet
      fillMatchedHistograms<JetTableMCDMatchedWeightedJoined::iterator, JetTableMCPMatchedWeightedJoined>(mcdjet, eventWeight);
    } // mcdjet
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCPMCDMatchedWeighted, "Full Jet finder MCP matched to MCD on weighted events", false);

  // Periodic cleanup to prevent unbounded memory growth
  /*void processCleanup(aod::Collision const&) {
  static int callCount = 0;
  callCount++;

  // Clean up cache every 50000 calls to prevent memory issues
  if (doMcClosure && callCount % 50000 == 0 && mcCollisionRandomValues.size() > 20000) {
  LOGF(info, "Cleaning up MC collision random values cache (size: %zu)", mcCollisionRandomValues.size());
  mcCollisionRandomValues.clear();

  // IMPROVEMENT: Add logging to verify our split ratios
  float mcdCount = registry.get<TH1>(HIST("h_MCD_splitevent_counter"))->GetBinContent(1);
  float mcpCount = registry.get<TH1>(HIST("h_MCP_splitevent_counter"))->GetBinContent(1);
  float matchedCount = registry.get<TH1>(HIST("h_Matched_splitevent_counter"))->GetBinContent(1);

  float totalEvents = mcdCount + mcpCount; // MCP and Matched should be the same, so don't double count
  float actualSplitFrac = totalEvents > 0 ? mcdCount / totalEvents : 0.0f;

  LOGF(info, "Current split statistics: MCD=%.1f, MCP=%.1f, Matched=%.1f", mcdCount, mcpCount, matchedCount);
  LOGF(info, "Actual split fraction: %.3f (target: %.3f)", actualSplitFrac, static_cast<float>(mcClosureSplitFrac));
  }
  }
  PROCESS_SWITCH(FullJetSpectra, processCleanup, "Periodic cleanup", true);
  */
  void processDataTracks(soa::Filtered<EMCCollisionsData>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;

    registry.fill(HIST("hCollisionsUnweighted"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hCollisionsUnweighted"), 1.5); // DetCollWithVertexZ

    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hCollisionsUnweighted"), 2.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hCollisionsUnweighted"), 3.5); // EventsNotSatisfyingEventSelection
      return;
    }
    // needed for the workaround to access EMCAL trigger bits. - This is needed for the MC productions in which the EMC trigger bits are missing. (MB MC LHC24f3, for ex.)
    // It first requires for atleast a cell in EMCAL to have energy content.
    // Once it finds a cell content,
    // it then checks if the collision is not an ambiguous collision (i.e. it has to be a unique collision = no bunch pile up)
    // If all of these conditions are satisfied then it checks for the required trigger bit in EMCAL.
    // For LHC22o, since the EMCAL didn't have hardware triggers, one would only require MB trigger (kTVXinEMC) in the EMCAL.

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hCollisionsUnweighted"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      // Check if EMCAL was readout with the MB trigger(kTVXinEMC) fired. If not then reject the event and exit the function.
      // This is the default check for the simulations with proper trigger flags not requiring the above workaround.
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hCollisionsUnweighted"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hCollisionsUnweighted"), 5.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hCollisionsUnweighted"), 6.5); // EMCAcceptedDetColl

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      // Fill Accepted events histos
      fillTrackHistograms(tracks, clusters, 1.0);
    }
    registry.fill(HIST("hCollisionsUnweighted"), 7.5); // EMCAcceptedCollAfterTrackSel
  }
  PROCESS_SWITCH(FullJetSpectra, processDataTracks, "Full Jet tracks for Data", false);

  void processMCTracks(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    double weight = 1.0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));

    registry.fill(HIST("hCollisionsUnweighted"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hCollisionsUnweighted"), 1.5); // DetCollWithVertexZ

    // for (auto const& track : tracks) {
    if (pTHat < pTHatAbsoluteMin) { // Track outlier rejection: should this be for every track iteration or for every collision?
      return;
    }
    // }
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hCollisionsUnweighted"), 2.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hCollisionsUnweighted"), 3.5); // EventsNotSatisfyingEventSelection
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hCollisionsUnweighted"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hCollisionsUnweighted"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hCollisionsUnweighted"), 5.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hCollisionsUnweighted"), 6.5); // EMCAcceptedDetColl

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      // Fill Accepted events histos
      fillTrackHistograms(tracks, clusters, 1.0);
    }
    registry.fill(HIST("hCollisionsUnweighted"), 7.5); // EMCAcceptedCollAfterTrackSel
  }
  PROCESS_SWITCH(FullJetSpectra, processMCTracks, "Full Jet tracks for MC", false);

  void processTracksWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision,
                             aod::JMcCollisions const&,
                             soa::Filtered<aod::JetTracks> const& tracks,
                             soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    float eventWeight = collision.mcCollision().weight();
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));

    registry.fill(HIST("hCollisionsWeighted"), 0.5, eventWeight); // AllWeightedDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hCollisionsWeighted"), 1.5, eventWeight); // WeightedCollWithVertexZ

    // for (auto const& track : tracks) {
    if (pTHat < pTHatAbsoluteMin) { // Track outlier rejection: should this be for every track iteration or for every collision?
      return;
    }
    // }
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hCollisionsWeighted"), 2.5, eventWeight); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hCollisionsWeighted"), 3.5, eventWeight); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doMBGapTrigger && eventWeight == 1) {
      registry.fill(HIST("hCollisionsWeighted"), 2.5, eventWeight); // MBRejectedDetEvents
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        fillTrackHistograms(tracks, clusters, eventWeight);
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hCollisionsWeighted"), 4.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hCollisionsWeighted"), 4.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hCollisionsWeighted"), 5.5, eventWeight); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hCollisionsWeighted"), 6.5); // EMCAcceptedWeightedDetColl

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      // Fill Accepted events histos
      fillTrackHistograms(tracks, clusters, eventWeight);
    }
    registry.fill(HIST("hCollisionsWeighted"), 7.5, eventWeight); // EMCAcceptedWeightedCollAfterTrackSel
  }
  PROCESS_SWITCH(FullJetSpectra, processTracksWeighted, "Full Jet tracks weighted", false);

  void processCollisionsWeightedWithMultiplicity(soa::Filtered<soa::Join<EMCCollisionsMCD, aod::FT0Mults>>::iterator const& collision, JetTableMCDWeightedJoined const& mcdjets, aod::JMcCollisions const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    float eventWeight = collision.mcCollision().weight();
    float neutralEnergy = 0.0;

    registry.fill(HIST("hEventmultiplicityCounter"), 0.5, eventWeight); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      registry.fill(HIST("hEventmultiplicityCounter"), 1.5, eventWeight); // DetCollWithVertexZ
      return;
    }
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hEventmultiplicityCounter"), 2.5, eventWeight); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hEventmultiplicityCounter"), 3.5, eventWeight); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doMBGapTrigger && eventWeight == 1) {
      registry.fill(HIST("hEventmultiplicityCounter"), 2.5, eventWeight); // MBRejectedDetEvents
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        fillTrackHistograms(tracks, clusters, eventWeight);
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hEventmultiplicityCounter"), 4.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hEventmultiplicityCounter"), 4.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hEventmultiplicityCounter"), 5.5, eventWeight); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hCollisionsWeighted"), 6.5); // EMCAcceptedWeightedDetColl

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 7.5, eventWeight); // EMCAcceptedWeightedCollAfterTrackSel
    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multiplicity(), eventWeight);

    for (auto const& mcdjet : mcdjets) {
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) { // MCD jets outlier rejection
        return;
      }
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      registry.fill(HIST("h2_full_jet_jetpTDetVsFT0Mults"), mcdjet.pt(), collision.multiplicity(), eventWeight);

      for (auto const& cluster : clusters) {
        neutralEnergy += cluster.energy();
      }
      auto nef = neutralEnergy / mcdjet.energy();
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef"), mcdjet.pt(), collision.multiplicity(), nef, eventWeight);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processCollisionsWeightedWithMultiplicity, "Weighted Collisions for Full Jets Multiplicity Studies", false);

  void processMBCollisionsWithMultiplicity(soa::Filtered<soa::Join<EMCCollisionsMCD, aod::FT0Mults>>::iterator const& collision, JetTableMCDJoined const& mcdjets, aod::JMcCollisions const&, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    float pTHat = 10. / (std::pow(1.0, 1.0 / pTHatExponent));
    float neutralEnergy = 0.0;

    registry.fill(HIST("hEventmultiplicityCounter"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 1.5); // DetCollWithVertexZ

    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hEventmultiplicityCounter"), 2.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hEventmultiplicityCounter"), 3.5); // EventsNotSatisfyingEventSelection
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        fillTrackHistograms(tracks, clusters, 1.0);
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hEventmultiplicityCounter"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hEventmultiplicityCounter"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hEventmultiplicityCounter"), 5.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 6.5); // EMCAcceptedDetColl

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 7.5); // EMCAcceptedCollAfterTrackSel
    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multiplicity());

    for (auto const& mcdjet : mcdjets) {
      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) { // MCD (Detector Level) Outlier Rejection
        return;
      }
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(mcdjet)) {
        continue;
      }
      registry.fill(HIST("h2_full_jet_jetpTDetVsFT0Mults"), mcdjet.pt(), collision.multiplicity(), 1.0);

      for (auto const& cluster : clusters) {
        neutralEnergy += cluster.energy();
      }
      auto nef = neutralEnergy / mcdjet.energy();
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef"), mcdjet.pt(), collision.multiplicity(), nef, 1.0);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processMBCollisionsWithMultiplicity, "MB MCD Collisions for Full Jets Multiplicity Studies", false);

  void processMBCollisionsDATAWithMultiplicity(soa::Filtered<soa::Join<EMCCollisionsData, aod::FT0Mults>>::iterator const& collision, FullJetTableDataJoined const& jets, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<aod::JetClusters> const& clusters)
  {
    bool eventAccepted = false;
    float neutralEnergy = 0.0;

    registry.fill(HIST("hEventmultiplicityCounter"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 1.5); // DetCollWithVertexZ
    if (doMBGapTrigger && collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hEventmultiplicityCounter"), 2.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hEventmultiplicityCounter"), 3.5); // EventsNotSatisfyingEventSelection
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        fillTrackHistograms(tracks, clusters, 1.0);
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hEventmultiplicityCounter"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hEventmultiplicityCounter"), 4.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hEventmultiplicityCounter"), 5.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 6.5); // EMCAcceptedDetColl

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 7.5); // EMCAcceptedCollAfterTrackSel
    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multiplicity());

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedJet<aod::JetTracks>(jet)) {
        continue;
      }
      registry.fill(HIST("h2_full_jet_jetpTDetVsFT0Mults"), jet.pt(), collision.multiplicity(), 1.0);

      for (auto const& cluster : clusters) {
        neutralEnergy += cluster.energy();
      }
      auto nef = neutralEnergy / jet.energy();
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef"), jet.pt(), collision.multiplicity(), nef, 1.0);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processMBCollisionsDATAWithMultiplicity, "MB DATA Collisions for Full Jets Multiplicity Studies", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FullJetSpectra>(cfgc)};
}
