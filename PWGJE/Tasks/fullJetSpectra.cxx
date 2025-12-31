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

/// TO DO : include histograms for cluster correction modes in the MC Mult processes.

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/EMCALClusterDefinition.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
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

#include <algorithm>
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
using ClusterWithCorrections = soa::Join<aod::JetClusters, aod::JClustersCorrectedEnergies>;

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
  Configurable<bool> doMBGapTrigger{"doMBGapTrigger", false, "set to true only when using MB-Gap Trigger JJ MC to reject MB events at the collision and track level"};

  // Software Trigger configurables
  Configurable<bool> doSoftwareTriggerSelection{"doSoftwareTriggerSelection", false, "set to true when using triggered datasets"};
  Configurable<std::string> triggerMasks{"triggerMasks", "fJetFullHighPt", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt"};

  // Jet configurables
  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.4, "resolution parameter for histograms without radius"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> jetpTMin{"jetpTMin", 20.0, "minimum jet pT"};
  Configurable<float> jetpTMax{"jetpTMax", 350., "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.3, "minimum jet eta for MCD"}; // each of these jet configurables are for the fiducial emcal cuts
  Configurable<float> jetEtaMax{"jetEtaMax", 0.3, "maximum jet eta for MCD"};  // for R = 0.4 (EMCAL eta acceptance: eta_jet = 0.7 - R)
  Configurable<float> jetPartEtaMin{"jetPartEtaMin", -0.3, "minimum jet eta for MCP with Fid"};
  Configurable<float> jetPartEtaMax{"jetPartEtaMax", 0.3, "maximum jet eta for MCP with Fid"};
  Configurable<float> jetNoFidPartEtaMin{"jetNoFidPartEtaMin", -0.9, "minimum jet eta for MCP w/o Fid"};
  Configurable<float> jetNoFidPartEtaMax{"jetNoFidPartEtaMax", 0.9, "maximum jet eta for MCP w/o Fid"};
  Configurable<float> emcalPhiMin{"emcalPhiMin", 1.3962634, "minimum emcal phi"};
  Configurable<float> emcalPhiMax{"emcalPhiMax", 3.2836100, "maximum emcal phi"};
  Configurable<float> jetPhiMin{"jetPhiMin", 1.8, "minimum emcal Fid phi"}; //For R = 0.4, jetPhiMin = emcalPhiMin + R
  Configurable<float> jetPhiMax{"jetPhiMax", 2.86, "maximum emcal Fid phi"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

  // Leading track and cluster pT configurables
  Configurable<float> minTrackPt{"minTrackPt", -99.0, "minimum pT selection on jet tracks"};
  Configurable<float> maxTrackPt{"maxTrackPt", 999.0, "maximum pT selection on jet tracks"};
  Configurable<float> minClusterPt{"minClusterPt", -99.0, "minimum pT selection on jet clusters"};
  Configurable<float> maxClusterPt{"maxClusterPt", 999.0, "maximum pT selection on jet clusters"};

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
  // const float kJetAreaFractionMinThreshold = -98.0f;
  const float kLeadingTrackPtMinThreshold = -98.0f;
  const float kLeadingTrackPtMaxThreshold = 9998.0f;
  const float kLeadingClusterPtMinThreshold = -98.0f;
  const float kLeadingClusterPtMaxThreshold = 9998.0f;

  std::vector<int> eventSelectionBits;
  std::vector<bool> filledJetR;
  std::vector<double> jetRadiiValues;
  std::vector<int> triggerMaskBits;
  std::string particleSelection;

  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Instantiate the Zorro processor for skimmed data and define an output object
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  bool doSumw2 = false;

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
    if (doprocessBCs) {
      auto hBCCounter = registry.get<TH1>(HIST("hBCCounter"));
      hBCCounter->GetXaxis()->SetBinLabel(1, "AllBC");
      hBCCounter->GetXaxis()->SetBinLabel(2, "BC+kTVX");
      hBCCounter->GetXaxis()->SetBinLabel(3, "BC+kTVX+NoTFB");
      hBCCounter->GetXaxis()->SetBinLabel(4, "BC+kTVX+NoTFB+NoITSROFB");
      hBCCounter->GetXaxis()->SetBinLabel(5, "CollinBC");
      hBCCounter->GetXaxis()->SetBinLabel(6, "CollinBC+kTVX");
      hBCCounter->GetXaxis()->SetBinLabel(7, "CollinBC+kTVX+Sel8");
      hBCCounter->GetXaxis()->SetBinLabel(8, "CollinBC+kTVX+Sel8Full");
      hBCCounter->GetXaxis()->SetBinLabel(9, "CollinBC+kTVX+Sel8Full+GoodZvtx");
      hBCCounter->GetXaxis()->SetBinLabel(10, "CollinBC+kTVX+Sel8Full+VtxZ+GoodZvtx");
    }

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

    if (doprocessJetsTriggeredData) {
      auto hDetTrigcollisionCounter = registry.get<TH1>(HIST("hDetTrigcollisionCounter"));
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(1, "allDetTrigColl");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(2, "DetTrigCollAfterZorroSelection");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(3, "DetTrigCollWithVertexZ");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(4, "EventsNotSatisfyingEvent+TriggerSelection");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(5, "OnlyFullJetHighPt+NoFullJetLowPt");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(6, "OnlyFullJetLowPt");
      // hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(7, "OnlyLowPt+NoMB");
      // hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(8, "OnlyMB");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(7, "FullJetHighPt+FullJetLowPt");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(8, "FullJetHighPt+MB");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(9, "FullJetLowPt+MB");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(10, "AllRejectedTrigOverlaps");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(11, "EMCAcceptedDetTrigCollAfterTrigOverlapChecks");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(12, "AllRejectedDetTrigEventsAfterEMCEventSelection");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(13, "EMCAcceptedDetTrigColl");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(14, "EMCAcceptedDetTrigCollWithLowChargedJetTriggers");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(15, "EMCAcceptedDetTrigCollWithHighChargedJetTriggers");
      hDetTrigcollisionCounter->GetXaxis()->SetBinLabel(16, "EMCAcceptedDetTrigCollWithLow+HighFullJetTriggers");
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

    if (doprocessJetsNoFidMCPMCDMatchedWeighted) {
      auto hMatchedNoFidcollisionCounter = registry.get<TH1>(HIST("hMatchedNoFidcollisionCounter"));
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(1, "allDetColl");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(3, "RejectedDetCollWithOutliers");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(4, "RejectedPartCollWithOutliers");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(5, "EMCMBRejectedDetColl");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(6, "EventsNotSatisfyingEventSelection");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(7, "EMCreadoutDetJJEventsWithkTVXinEMC");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(8, "AllRejectedDetEventsAfterEMCEventSelection");
      hMatchedNoFidcollisionCounter->GetXaxis()->SetBinLabel(9, "EMCAcceptedDetColl");
    }

    if (doprocessMBCollisionsDATAWithMultiplicity || doprocessMBMCDCollisionsWithMultiplicity) {
      auto hEventmultiplicityCounter = registry.get<TH1>(HIST("hEventmultiplicityCounter"));
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(1, "allDetColl");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(2, "DetCollWithVertexZ");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(3, "EventsNotSatisfyingEventSelection");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(4, "EMCreadoutDetEventsWithkTVXinEMC");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(5, "AllRejectedEventsAfterEMCEventSelection");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(6, "EMCAcceptedDetColl");
    }

    if (doprocessMCDCollisionsWeightedWithMultiplicity) {
      auto hEventmultiplicityCounter = registry.get<TH1>(HIST("hEventmultiplicityCounter"));
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(1, "allWeightedDetColl");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(2, "WeightedDetCollWithVertexZ");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(3, "MBRejectedDetEvents");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(4, "WeightedEventsNotSatisfyingEventSelection");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(5, "EMCreadoutWeightedDetEventsWithkTVXinEMC");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(6, "AllRejectedWeightedEventsAfterEMCEventSelection");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(7, "EMCAcceptedWeightedDetColl");
      hEventmultiplicityCounter->GetXaxis()->SetBinLabel(8, "EMCAcceptedWeightedCollAfterTrackSel");
    }

    if (doprocessMBMCPCollisionsWithMultiplicity) {
      auto hPartEventmultiplicityCounter = registry.get<TH1>(HIST("hPartEventmultiplicityCounter"));
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(1, "allMcColl");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(2, "McCollWithVertexZ");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(3, "RejectedPartCollWithOutliers");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(4, "MBRejectedPartEvents");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(5, "RejectedPartCollForDetCollWithSize0or<1");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(6, "AcceptedPartCollWithSize>=1");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(7, "EMCreadoutDetEventsWithkTVXinEMC");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(8, "AllRejectedPartEventsAfterEMCEventSelection");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(9, "EMCAcceptedPartColl");
    }

    if (doprocessMBMCPCollisionsWeightedWithMultiplicity) {
      auto hPartEventmultiplicityCounter = registry.get<TH1>(HIST("hPartEventmultiplicityCounter"));
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(1, "allWeightedMcColl");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(2, "WeightedMcCollWithVertexZ");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(3, "RejectedWeightedPartCollWithOutliers");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(4, "MBRejectedPartEvents");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(5, "RejectedWeightedPartCollForDetCollWithSize0or<1");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(6, "AcceptedWeightedPartCollWithSize>=1");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(7, "EMCreadoutDetEventsWithkTVXinEMC");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(8, "AllRejectedWeightedPartEventsAfterEMCEventSelection");
      hPartEventmultiplicityCounter->GetXaxis()->SetBinLabel(9, "EMCAcceptedWeightedPartColl");
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
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    particleSelection = static_cast<std::string>(particleSelections);
    jetRadiiValues = (std::vector<double>)jetRadii;
    doSumw2 = doMBGapTrigger;
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

    // Sanity Log check
    if (doSumw2) {
      LOGF(info, "HistogramRegistry initialized with Sumw2 = ON (weighted JJ MC mode).");
    } else {
      LOGF(info, "HistogramRegistry initialized with Sumw2 = OFF (unweighted mode).");
    }

    if (doprocessBCs) {
      registry.add("hBCCounter", "", {HistType::kTH1F, {{11, 0.0, 11.}}}, doSumw2);
    }

    // Track QA histograms
    if (doprocessDataTracks || doprocessMCTracks || doprocessTracksWeighted) {
      registry.add("hCollisionsUnweighted", "event status; event status;entries", {HistType::kTH1F, {{12, 0., 12.0}}}, doSumw2);

      registry.add("h_track_pt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_track_eta", "track #eta;#eta_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_track_phi", "track #varphi;#varphi_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_track_energy", "track energy;Energy of tracks;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_track_energysum", "track energy sum;Sum of track energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);

      // Cluster QA histograms
      registry.add("h_clusterTime", "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}}, doSumw2);
      registry.add("h_cluster_pt", "cluster pT;#it{p}_{T_cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, doSumw2);
      registry.add("h_cluster_eta", "cluster #eta;#eta_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_cluster_phi", "cluster #varphi;#varphi_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_cluster_energysum_uncorr", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energysum_corr_oneTrack100", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energysum_corr_oneTrack70", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energysum_corr_allTracks100", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energysum_corr_allTracks70", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);

      registry.add("h_cluster_energy_uncorr", "Cluster Energy (uncorrected); E_{cluster} [GeV]; N_{clusters}", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energy_corr_oneTrack100", "Cluster Energy (HadCorr, 1track, 100%); E_{cluster}^{corr,1,100} [GeV]; N_{clusters}", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energy_corr_oneTrack70", "Cluster Energy (HadCorr, 1track, 70%); E_{cluster}^{corr,1,70} [GeV]; N_{clusters}", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energy_corr_allTracks100", "Cluster Energy (HadCorr, alltracks, 100%); E_{cluster}^{corr,all,100} [GeV]; N_{clusters}", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_cluster_energy_corr_allTracks70", "Cluster Energy (HadCorr, alltracks, 70%); E_{cluster}^{corr,all,70} [GeV]; N_{clusters}", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);

      if (doprocessTracksWeighted) {
        registry.add("hCollisionsWeighted", "event status;event status;entries", {HistType::kTH1F, {{12, 0.0, 12.0}}}, doSumw2);
      }
    }

    // Jet QA histograms
    if (doprocessJetsData || doprocessJetsMCD || doprocessJetsMCDWeighted || doprocessJetsTriggeredData) {

      registry.add("hDetcollisionCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.}}}, doSumw2);

      registry.add("h_full_jet_pt", "#it{p}_{T,jet};#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_jet_pt_pTHatcut", "#it{p}_{T,jet};#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_jet_phi", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_full_jet_clusterTime", "Time of cluster", HistType::kTH1F, {{500, -250, 250, "#it{t}_{cls} (ns)"}}, doSumw2);
      registry.add("h2_full_jet_nef", "#it{p}_{T,jet} vs nef at Det Level; #it{p}_{T,jet} (GeV/#it{c});nef", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}}, doSumw2);
      registry.add("h2_full_jet_nef_rejected", "#it{p}_{T,jet} vs nef at Det Level for rejected events; #it{p}_{T,jet} (GeV/#it{c});nef", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}}, doSumw2);

      registry.add("h_Detjet_ntracks", "#it{p}_{T,track};#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h2_full_jet_chargedconstituents", "Number of charged constituents at Det Level;#it{p}_{T,jet} (GeV/#it{c});N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h2_full_jet_neutralconstituents", "Number of neutral constituents at Det Level;#it{p}_{T,jet} (GeV/#it{c});N_{ne}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h_full_jet_chargedconstituents_pt", "track pT;#it{p}^{T,jet}_{track} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_jet_chargedconstituents_eta", "track #eta;#eta^{jet}_{track};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_jet_chargedconstituents_phi", "track #varphi;#varphi^{jet}_{track};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_full_jet_chargedconstituents_energy", "track energy;Energy of tracks;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_full_jet_chargedconstituents_energysum", "track energy sum;Sum of track energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);

      registry.add("h_full_jet_neutralconstituents_pt_uncorr", "cluster pT;#it{p}^{T,jet}_{cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_pt_corr_oneTrack100", "cluster pT;#it{p}^{T,jet}_{cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_pt_corr_oneTrack70", "cluster pT;#it{p}^{T,jet}_{cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_pt_corr_allTracks100", "cluster pT;#it{p}^{T,jet}_{cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_pt_corr_allTracks70", "cluster pT;#it{p}^{T,jet}_{cluster} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}, doSumw2);

      registry.add("h_full_jet_neutralconstituents_eta", "cluster #eta;#eta^{jet}_{cluster};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_phi", "cluster #varphi;#varphi^{jet}_{cluster};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energy", "cluster energy;Energy of cluster;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energysum", "cluster energy sum;Sum of cluster energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h2_full_jettrack_pt", "#it{p}_{T,jet} vs #it{p}_{T,track}; #it{p}_{T,jet} (GeV/#it{c});#it{p}_{T,track} (GeV/#it{c})", {HistType::kTH2F, {{350, 0., 350.}, {200, 0., 200.}}}, doSumw2);
      registry.add("h2_full_jettrack_eta", "jet #eta vs jet_track #eta; #eta_{jet};#eta_{track}", {HistType::kTH2F, {{100, -1., 1.}, {500, -1., 1.}}}, doSumw2);
      registry.add("h2_full_jettrack_phi", "jet #varphi vs jet_track #varphi; #varphi_{jet}; #varphi_{track}", {HistType::kTH2F, {{160, 0., 7.}, {160, -1., 7.}}}, doSumw2);

      registry.add("h2_track_etaphi", "jet_track #eta vs jet_track #varphi; #eta_{track};#varphi_{track}", {HistType::kTH2F, {{500, -1., 1.}, {160, -1., 7.}}}, doSumw2);
      registry.add("h2_jet_etaphi", "jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}}, doSumw2);

      // NEW: Jet constituent histograms for each hadronic correction mode
      registry.add("h_full_jet_neutralconstituents_energy_uncorr", "Jet neutral cluster energy (uncorr)", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energy_corr_oneTrack100", "Jet neutral cluster energy (corr 1track100)", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energy_corr_oneTrack70", "Jet neutral cluster energy (corr 1track70)", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energy_corr_allTracks100", "Jet neutral cluster energy (corr alltracks100)", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energy_corr_allTracks70", "Jet neutral cluster energy (corr alltracks70)", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      // Corrected NEF histograms for the corresponding correction mode
      registry.add("h2_full_jet_nef_uncorr", "Jet pT vs NEF (uncorr); p_{T,jet}; NEF", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 1.}}}, doSumw2);
      registry.add("h2_full_jet_nef_corr_oneTrack100", "Jet pT vs NEF (corr, 1track100); p_{T,jet}; NEF", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 1.}}}, doSumw2);
      registry.add("h2_full_jet_nef_corr_oneTrack70", "Jet pT vs NEF (corr, 1track70); p_{T,jet}; NEF", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 1.}}}, doSumw2);
      registry.add("h2_full_jet_nef_corr_allTracks100", "Jet pT vs NEF (corr, alltracks100); p_{T,jet}; NEF", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 1.}}}, doSumw2);
      registry.add("h2_full_jet_nef_corr_allTracks70", "Jet pT vs NEF (corr, alltracks70); p_{T,jet}; NEF", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 1.}}}, doSumw2);
    }
    if (doprocessJetsTriggeredData) {
      registry.add("hDetTrigcollisionCounter", "event status;;entries", {HistType::kTH1F, {{17, 0.0, 17.}}}, doSumw2);
      // registry.add("h2_full_jet_nef_rejected", "#it{p}_{T,jet} vs nef at Det Level for rejected events; #it{p}_{T,jet} (GeV/#it{c});nef", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}}, doSumw2);
    }
    if (doprocessJetsMCP || doprocessJetsMCPWeighted) {
      registry.add("hPartcollisionCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}, doSumw2);
      registry.add("hRecoMatchesPerMcCollision", "split vertices QA;;entries", {HistType::kTH1F, {{5, 0.0, 5.0}}}, doSumw2);

      registry.add("h_full_mcpjet_tablesize", "", {HistType::kTH1F, {{4, 0., 5.}}}, doSumw2);
      registry.add("h_full_mcpjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}}, doSumw2);
      registry.add("h_full_jet_pt_part", "jet pT;#it{p}_{T_jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_jet_eta_part", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_jet_phi_part", "jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h2_full_jet_nef_part", "#it{p}_{T,jet} vs nef at Part Level;#it{p}_{T,jet} (GeV/#it{c});nef", {HistType::kTH2F, {{350, 0., 350.}, {105, 0., 1.05}}}, doSumw2);

      registry.add("h_Partjet_ntracks", "#it{p}_{T,constituent};#it{p}_{T_constituent} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h2_full_jet_chargedconstituents_part", "Number of charged constituents at Part Level;#it{p}_{T,jet} (GeV/#it{c});N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h2_full_jet_neutralconstituents_part", "Number of neutral constituents at Part Level;#it{p}_{T,jet} (GeV/#it{c});N_{ne}", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_pt_part", "#it{p}_{T} of neutral constituents at Part Level;#it{p}_{T,ne} (GeV/#it{c}); entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_eta_part", "#eta of neutral constituents at Part Level;#eta_{ne};entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_phi_part", "#varphi of neutral constituents at Part Level;#varphi_{ne};entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energy_part", "neutral constituents' energy;Energy of neutral constituents;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);
      registry.add("h_full_jet_neutralconstituents_energysum_part", "neutral constituents' energy sum;Sum of neutral constituents' energy per event;entries", {HistType::kTH1F, {{400, 0., 400.}}}, doSumw2);

      registry.add("h2_jettrack_pt_part", "#it{p}_{T,jet} vs #it{p}_{T_track}; #it{p}_{T_jet} (GeV/#it{c});#it{p}_{T_track} (GeV/#it{c})", {HistType::kTH2F, {{350, 0., 350.}, {200, 0., 200.}}}, doSumw2);
      registry.add("h2_jettrack_eta_part", "jet #eta vs jet_track #eta; #eta_{jet};#eta_{track}", {HistType::kTH2F, {{100, -1., 1.}, {500, -1., 1.}}}, doSumw2);
      registry.add("h2_jettrack_phi_part", "jet #varphi vs jet_track #varphi; #varphi_{jet}; #varphi_{track}", {HistType::kTH2F, {{160, 0., 7.}, {160, -1., 7.}}}, doSumw2);

      registry.add("h2_track_etaphi_part", "jet_track #eta vs jet_track #varphi; #eta_{track};#varphi_{track}", {HistType::kTH2F, {{500, -1., 1.}, {160, -1., 7.}}}, doSumw2);
      registry.add("h2_jet_etaphi_part", "jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}}, doSumw2);

      registry.add("h2_full_mcpjetOutsideFiducial_pt", "MCP jet outside EMC Fiducial Acceptance #it{p}_{T,part};#it{p}_{T,part} (GeV/c); Ncounts", {HistType::kTH2F, {{350, 0., 350.}, {10000, 0., 10000.}}}, doSumw2);
      registry.add("h_full_mcpjetOutside_eta_part", "MCP jet #eta outside EMC Fiducial Acceptance;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_mcpjetOutside_phi_part", "MCP jet #varphi outside EMC Fiducial Acceptance;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h2_full_mcpjetInsideFiducial_pt", "MCP jet #it{p}_{T,part} inside EMC Fiducial Acceptance;#it{p}_{T,part} (GeV/c); Ncounts", {HistType::kTH2F, {{350, 0., 350.}, {10000, 0., 10000.}}}, doSumw2);
      registry.add("h_full_mcpjetInside_eta_part", "MCP jet #eta inside EMC Fiducial Acceptance;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_mcpjetInside_phi_part", "MCP jet #varphi inside EMC Fiducial Acceptance;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
    }

    if (doprocessJetsMCPMCDMatched || doprocessJetsMCPMCDMatchedWeighted) {
      registry.add("hMatchedcollisionCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}, doSumw2);
      registry.add("h_allMatchedPartJetsCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}, doSumw2);

      registry.add("h_full_matchedmcdjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_matchedmcpjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_matchedmcdjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}}, doSumw2);
      registry.add("h_full_matchedmcpjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}}, doSumw2);
      // registry.add("h_full_matchedmcdjet_eta", "Matched MCD jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      // registry.add("h_full_matchedmcdjet_phi", "Matched MCD jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_full_matchedmcpjet_eta", "Matched MCP jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_matchedmcpjet_phi", "Matched MCP jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_allMatchedPartJetsPt", "Matched MCP jet Pt;p_{T,part} (GeV/c);entries", {HistType::kTH1F, {{350, 0.0, 350.0}}}, doSumw2);
      registry.add("h_full_jet_deltaR", "Distance between matched Det Jet and Part Jet; #Delta R; entries", {HistType::kTH1F, {{100, 0., 1.}}}, doSumw2);

      registry.add("h2_full_jet_energyscaleDet", "Jet Energy Scale (det); p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);

      registry.add("h2_matchedjet_etaphiDet", "Det jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}}, doSumw2);
      registry.add("h2_matchedjet_etaphiPart", "Part jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}}, doSumw2);
      registry.add("h2_matchedjet_deltaEtaCorr", "Correlation between Det Eta and Part Eta; #eta_{jet,det}; #eta_{jet,part}", {HistType::kTH2F, {{100, -1., 1.}, {100, -1., 1.}}}, doSumw2);
      registry.add("h2_matchedjet_deltaPhiCorr", "Correlation between Det Phi and Part Phi; #varphi_{jet,det}; #varphi_{jet,part}", {HistType::kTH2F, {{160, 0., 7.}, {160, 0., 7.}}}, doSumw2);

      registry.add("h2_full_jet_energyscalePart", "Jet Energy Scale (part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      registry.add("h3_full_jet_energyscalePart", "R dependence of Jet Energy Scale (Part); #it{R}_{jet};p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      registry.add("h2_full_jet_etaresolutionPart", ";p_{T,part} (GeV/c); (#eta_{jet,det} - #eta_{jet,part})/#eta_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {100, -1., 1.}}}, doSumw2);
      registry.add("h2_full_jet_phiresolutionPart", ";p_{T,part} (GeV/c); (#varphi_{jet,det} - #varphi_{jet,part})/#varphi_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {160, -1., 7.}}}, doSumw2);
      // registry.add("h2_full_jet_energyscaleChargedPart", "Jet Energy Scale (charged part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      // registry.add("h2_full_jet_energyscaleNeutralPart", "Jet Energy Scale (neutral part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      // registry.add("h2_full_jet_energyscaleChargedVsFullPart", "Jet Energy Scale (charged part, vs. full jet pt); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      // registry.add("h2_full_jet_energyscaleNeutralVsFullPart", "Jet Energy Scale (neutral part, vs. full jet pt); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      registry.add("h2_full_fakemcdjets", "Fake MCD Jets; p_{T,det} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h2_full_fakemcpjets", "Fake MCP Jets; p_{T,part} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h2_full_matchedmcpjet_pt", "Matched MCP jet in EMC Fiducial Acceptance #it{p}_{T,part};#it{p}_{T,part} (GeV/c); Ncounts", {HistType::kTH2F, {{350, 0., 350.}, {10000, 0., 10000.}}}, doSumw2);

      // Response Matrix
      registry.add("h_full_jet_ResponseMatrix", "Full Jets Response Matrix; p_{T,det} (GeV/c); p_{T,part} (GeV/c)", {HistType::kTH2F, {{500, 0., 500.}, {500, 0., 500.}}}, doSumw2);
    }

    if (doprocessJetsNoFidMCPMCDMatchedWeighted) {
      registry.add("hMatchedNoFidcollisionCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}, doSumw2);
      registry.add("h_allMatchedNoFidPartJetsCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}, doSumw2);

      registry.add("h_full_NoFidmatchedmcdjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_NoFidmatchedmcpjet_tablesize", "", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_full_NoFidmatchedmcdjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}}, doSumw2);
      registry.add("h_full_NoFidmatchedmcpjet_ntracks", "", {HistType::kTH1F, {{200, -0.5, 200.}}}, doSumw2);
      registry.add("h_full_NoFidmatchedmcpjet_eta", "Matched No Fid MCP jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1., 1.}}}, doSumw2);
      registry.add("h_full_NoFidmatchedmcpjet_phi", "Matched No Fid MCP jet #varphi;#varphi_{jet};entries", {HistType::kTH1F, {{160, 0., 7.}}}, doSumw2);
      registry.add("h_allMatchedNoFidPartJetsPt", "Matched No Fid MCP jet Pt;p_{T,part} (GeV/c);entries", {HistType::kTH1F, {{350, 0.0, 350.0}}}, doSumw2);
      registry.add("h_full_jet_NoFiddeltaR", "Distance between matched Det Jet and Part Jet; #Delta R; entries", {HistType::kTH1F, {{100, 0., 1.}}}, doSumw2);

      registry.add("h2_full_jet_NoFidenergyscaleDet", "Jet Energy Scale (det); p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);

      registry.add("h2_NoFidmatchedjet_etaphiDet", "Det jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}}, doSumw2);
      registry.add("h2_NoFidmatchedjet_etaphiPart", "Part jet #eta vs jet #varphi; #eta_{jet};#varphi_{jet}", {HistType::kTH2F, {{100, -1., 1.}, {160, -1., 7.}}}, doSumw2);
      registry.add("h2_NoFidmatchedjet_deltaEtaCorr", "Correlation between Det Eta and Part Eta; #eta_{jet,det}; #eta_{jet,part}", {HistType::kTH2F, {{100, -1., 1.}, {100, -1., 1.}}}, doSumw2);
      registry.add("h2_NoFidmatchedjet_deltaPhiCorr", "Correlation between Det Phi and Part Phi; #varphi_{jet,det}; #varphi_{jet,part}", {HistType::kTH2F, {{160, 0., 7.}, {160, 0., 7.}}}, doSumw2);

      registry.add("h2_full_jet_NoFidenergyscalePart", "Jet Energy Scale (part); p_{T,part} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH2F, {{400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      registry.add("h3_full_jet_NoFidenergyscalePart", "R dependence of Jet Energy Scale (Part); #it{R}_{jet};p_{T,det} (GeV/c); (p_{T,det} - p_{T,part})/p_{T,part}", {HistType::kTH3F, {{jetRadiiBins, ""}, {400, 0., 400.}, {200, -1., 1.}}}, doSumw2);
      registry.add("h2_full_jet_NoFidetaresolutionPart", ";p_{T,part} (GeV/c); (#eta_{jet,det} - #eta_{jet,part})/#eta_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {100, -1., 1.}}}, doSumw2);
      registry.add("h2_full_jet_NoFidphiresolutionPart", ";p_{T,part} (GeV/c); (#varphi_{jet,det} - #varphi_{jet,part})/#varphi_{jet,part}", {HistType::kTH2F, {{400, 0., 400.}, {160, -1., 7.}}}, doSumw2);
      registry.add("h2_full_NoFidfakemcdjets", "Fake MCD Jets; p_{T,det} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h2_full_NoFidfakemcpjets", "Fake MCP Jets; p_{T,part} (GeV/c); NCounts", {HistType::kTH2F, {{350, 0., 350.}, {100, 0., 100.}}}, doSumw2);
      registry.add("h2_full_NoFidmatchedmcpjet_pt", "Matched No Fid MCP jet #it{p}_{T,part};#it{p}_{T,part} (GeV/c); Ncounts", {HistType::kTH2F, {{350, 0., 350.}, {10000, 0., 10000.}}}, doSumw2);

      // Response Matrix
      registry.add("h_full_jet_NoFidResponseMatrix", "Full Jets Response Matrix; p_{T,det} (GeV/c); p_{T,part} (GeV/c)", {HistType::kTH2F, {{500, 0., 500.}, {500, 0., 500.}}}, doSumw2);
    }

    if (doprocessMBCollisionsDATAWithMultiplicity || doprocessMBMCDCollisionsWithMultiplicity || doprocessMCDCollisionsWeightedWithMultiplicity) {
      registry.add("hEventmultiplicityCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}, doSumw2);
      registry.add("h_FT0Mults_occupancy", "", {HistType::kTH1F, {{3500, 0., 3500.}}}, doSumw2);

      registry.add("h_all_fulljet_Njets", "Full Jet Multiplicity (per Event)", {HistType::kTH1F, {{20, 0., 20.}}}, doSumw2);
      registry.add("h_Leading_full_jet_pt", "#it{p}_{T,leading jet};#it{p}_{T_leading jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h2_full_jet_leadingJetPt_vs_counts", ";#it{p}_{T_leading jet} (GeV/#it{c}); Counts", {HistType::kTH2F, {{350, 0., 350.}, {20, 0., 20.}}}, doSumw2);
      registry.add("h_SubLeading_full_jet_pt", "#it{p}_{T,leading jet};#it{p}_{T_leading jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h2_full_jet_subLeadingJetPt_vs_counts", ";#it{p}_{T_leading jet} (GeV/#it{c}); Counts", {HistType::kTH2F, {{350, 0., 350.}, {20, 0., 20.}}}, doSumw2);
      // Inside Jet Loop:
      // CASE 1:
      registry.add("h_all_fulljet_pt", "#it{p}_{T,fulljet};#it{p}_{T_fulljet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_all_fulljet_Nch", ";N_{ch};", {HistType::kTH1F, {{50, 0., 50.}}}, doSumw2);
      registry.add("h_all_fulljet_NEF_uncorr", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h_all_fulljet_NEF_corr_oneTrack100", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h_all_fulljet_NEF_corr_oneTrack70", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h_all_fulljet_NEF_corr_allTracks100", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h_all_fulljet_NEF_corr_allTracks70", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);

      registry.add("h2_all_fulljet_jetpTDet_vs_FT0Mults", "; p_{T,det} (GeV/c); FT0M Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}}, doSumw2);
      registry.add("h2_all_fulljet_jetpTDet_vs_Nch", ";#it{p}_{T_fulljet} (GeV/#it{c}); N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {50, 0., 50.}}}, doSumw2);
      registry.add("h3_full_jet_jetpTDet_FT0Mults_nef_uncorr", "; p_{T,det} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      registry.add("h3_full_jet_jetpTDet_FT0Mults_nef_corr_oneTrack100", "; p_{T,det} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      registry.add("h3_full_jet_jetpTDet_FT0Mults_nef_corr_oneTrack70", "; p_{T,det} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      registry.add("h3_full_jet_jetpTDet_FT0Mults_nef_corr_allTracks100", "; p_{T,det} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      registry.add("h3_full_jet_jetpTDet_FT0Mults_nef_corr_allTracks70", "; p_{T,det} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      // CASE 2:
      registry.add("h_leading_fulljet_pt", "#it{p}_{T,Leading fulljet};#it{p}_{T_Leadingfulljet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_leading_fulljet_Nch", ";N_{ch};", {HistType::kTH1F, {{50, 0., 50.}}}, doSumw2);
      registry.add("h_leading_fulljet_NEF", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h2_leading_fulljet_jetpTDet_vs_FT0Mults", ";Leading p_{T,det} (GeV/c); FT0M Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}}, doSumw2);
      registry.add("h2_leading_fulljet_jetpTDet_vs_Nch", ";#it{p}_{T_Leadingfulljet} (GeV/#it{c}); N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {50, 0., 50.}}}, doSumw2);
      registry.add("h3_leading_fulljet_jetpTDet_FT0Mults_nef", "; Leading p_{T,det} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      // CASE 3:
      registry.add("h_subleading_fulljet_pt", "#it{p}_{T,SubLeading fulljet};#it{p}_{T_SubLeadingfulljet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_subleading_fulljet_Nch", ";N_{ch};", {HistType::kTH1F, {{50, 0., 50.}}}, doSumw2);
      registry.add("h_subleading_fulljet_NEF", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h2_subleading_fulljet_jetpTDet_vs_FT0Mults", ";SubLeading p_{T,det} (GeV/c); FT0M Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}}, doSumw2);
      registry.add("h2_subleading_fulljet_jetpTDet_vs_Nch", ";#it{p}_{T_SubLeadingfulljet} (GeV/#it{c}); N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {50, 0., 50.}}}, doSumw2);
      registry.add("h3_subleading_fulljet_jetpTDet_FT0Mults_nef", "; SubLeading p_{T,det} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
    }

    if (doprocessMBMCPCollisionsWithMultiplicity || doprocessMBMCPCollisionsWeightedWithMultiplicity) {
      registry.add("hPartEventmultiplicityCounter", "event status;event status;entries", {HistType::kTH1F, {{11, 0.0, 11.0}}}, doSumw2);
      registry.add("hRecoMatchesPerMcCollisionMult", "split vertices QA;;entries", {HistType::kTH1F, {{5, 0.0, 5.0}}}, doSumw2);
      registry.add("hMCCollMatchedFT0Mult", "", {HistType::kTH1F, {{3500, 0., 3500.}}}, doSumw2);
      registry.add("hMCCollMatchedFT0Cent", "", {HistType::kTH1F, {{105, 0., 105.}}}, doSumw2);

      registry.add("h_all_fulljet_Njets_part", "Full Jet Multiplicity (per Event)", {HistType::kTH1F, {{20, 0., 20.}}}, doSumw2);
      registry.add("h_Leading_full_jet_pt_part", "#it{p}_{T,leading jet};#it{p}_{T_leading jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h2_full_jet_leadingJetPt_vs_counts_part", ";#it{p}_{T_leading jet} (GeV/#it{c}); Counts", {HistType::kTH2F, {{350, 0., 350.}, {20, 0., 20.}}}, doSumw2);
      registry.add("h_SubLeading_full_jet_pt_part", "#it{p}_{T,leading jet};#it{p}_{T_leading jet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h2_full_jet_subLeadingJetPt_vs_counts_part", ";#it{p}_{T_leading jet} (GeV/#it{c}); Counts", {HistType::kTH2F, {{350, 0., 350.}, {20, 0., 20.}}}, doSumw2);

      // Inside Jet Loop:
      // CASE 1:
      registry.add("h_all_fulljet_pt_part", "#it{p}_{T,fulljet};#it{p}_{T_fulljet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_all_fulljet_Nch_part", ";N_{ch};", {HistType::kTH1F, {{50, 0., 50.}}}, doSumw2);
      registry.add("h_all_fulljet_Nne_part", ";N_{ne};", {HistType::kTH1F, {{100, 0., 100.}}}, doSumw2);
      registry.add("h_all_fulljet_NEF_part", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h2_all_fulljet_jetpT_vs_FT0Mults_part", "; p_{T,part} (GeV/c); FT0M Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}}, doSumw2);
      registry.add("h2_all_fulljet_jetpT_vs_Nch_part", ";#it{p}_{T_fulljet} (GeV/#it{c}); N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {50, 0., 50.}}}, doSumw2);
      registry.add("h3_full_jet_jetpT_FT0Mults_nef_part", "; p_{T,part} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      // CASE 2:
      registry.add("h_leading_fulljet_pt_part", "#it{p}_{T,Leading fulljet};#it{p}_{T_Leadingfulljet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_leading_fulljet_Nch_part", ";N_{ch};", {HistType::kTH1F, {{50, 0., 50.}}}, doSumw2);
      registry.add("h_leading_fulljet_Nne_part", ";N_{ne};", {HistType::kTH1F, {{100, 0., 100.}}}, doSumw2);
      registry.add("h_leading_fulljet_NEF_part", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h2_leading_fulljet_jetpT_vs_FT0Mults_part", ";Leading p_{T,part} (GeV/c); FT0M Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}}, doSumw2);
      registry.add("h2_leading_fulljet_jetpT_vs_Nch_part", ";#it{p}_{T_Leadingfulljet} (GeV/#it{c}); N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {50, 0., 50.}}}, doSumw2);
      registry.add("h3_leading_fulljet_jetpT_FT0Mults_nef_part", "; Leading p_{T,part} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
      // CASE 3:
      registry.add("h_subleading_fulljet_pt_part", "#it{p}_{T,SubLeading fulljet};#it{p}_{T_SubLeadingfulljet} (GeV/#it{c});entries", {HistType::kTH1F, {{350, 0., 350.}}}, doSumw2);
      registry.add("h_subleading_fulljet_Nch_part", ";N_{ch};", {HistType::kTH1F, {{50, 0., 50.}}}, doSumw2);
      registry.add("h_subleading_fulljet_Nne_part", ";N_{ne};", {HistType::kTH1F, {{100, 0., 100.}}}, doSumw2);
      registry.add("h_subleading_fulljet_NEF_part", ";NEF;", {HistType::kTH1F, {{105, 0., 1.05}}}, doSumw2);
      registry.add("h2_subleading_fulljet_jetpT_vs_FT0Mults_part", ";SubLeading p_{T,part} (GeV/c); FT0M Multiplicity", {HistType::kTH2F, {{350, 0., 350.}, {3500, 0., 3500.}}}, doSumw2);
      registry.add("h2_subleading_fulljet_jetpT_vs_Nch_part", ";#it{p}_{T_SubLeadingfulljet} (GeV/#it{c}); N_{ch}", {HistType::kTH2F, {{350, 0., 350.}, {50, 0., 50.}}}, doSumw2);
      registry.add("h3_subleading_fulljet_jetpT_FT0Mults_nef_part", "; SubLeading p_{T,part} (GeV/c); FT0M Multiplicity, nef", {HistType::kTH3F, {{350, 0., 350.}, {50, 0., 50.}, {105, 0.0, 1.05}}}, doSumw2);
    }

    // Label the histograms
    labelCollisionHistograms(registry);
    // labelMCSplitHistogram(registry);
  } // init

  // Initialize CCDB access and histogram registry for Zorro processing
  template <typename BCType>
  void initCCDB(const BCType& bc)
  {
    if (doSoftwareTriggerSelection) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerMasks.value);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }
  }

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
  using EMCCollisionsData = o2::soa::Join<aod::JetCollisions, aod::JEMCCollisionLbs>; // JetCollisions with EMCAL Collision Labels
  using EMCCollisionsTriggeredData = o2::soa::Join<aod::JetCollisions, aod::JCollisionBCs, aod::JEMCCollisionLbs>;

  using EMCCollisionsMCD = o2::soa::Join<aod::JetCollisionsMCD, aod::JEMCCollisionLbs>; // where, JetCollisionsMCD = JetCollisions+JMcCollisionLbs

  using FullJetTableDataJoined = soa::Join<aod::FullJets, aod::FullJetConstituents>;
  using JetTableMCDJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents>;
  // using JetTableMCDWeightedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents, aod::FullMCDetectorLevelJetEventWeights>;
  using JetTableMCPJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents>;
  // using JetTableMCPWeightedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents, aod::FullMCParticleLevelJetEventWeights>;

  using JetTableMCDMatchedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents,
                                   aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets>;

  using JetTableMCPMatchedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents,
                               aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets>;

  // Commenting these out for now to avoid dependency of the task on JE EventWeights tables
  /*using JetTableMCDMatchedWeightedJoined = soa::Join<aod::FullMCDetectorLevelJets, aod::FullMCDetectorLevelJetConstituents,
                                           aod::FullMCDetectorLevelJetsMatchedToFullMCParticleLevelJets,
                                           aod::FullMCDetectorLevelJetEventWeights>;*/

  /*using JetTableMCPMatchedWeightedJoined = soa::Join<aod::FullMCParticleLevelJets, aod::FullMCParticleLevelJetConstituents,
                                           aod::FullMCParticleLevelJetsMatchedToFullMCDetectorLevelJets,
                                           aod::FullMCParticleLevelJetEventWeights>;*/

  // Applying some cuts(filters) on collisions, tracks, clusters

  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centFT0M >= centralityMin && aod::jcollision::centFT0M < centralityMax);
  // Filter EMCeventCuts = (nabs(aod::collision::posZ) < vertexZCut && aod::collision::centrality >= centralityMin && aod::collision::centrality < centralityMax);
  Filter trackCuts = (aod::jtrack::pt >= trackpTMin && aod::jtrack::pt < trackpTMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax && aod::jtrack::phi >= trackPhiMin && aod::jtrack::phi <= trackPhiMax);
  aod::EMCALClusterDefinition clusterDefinition = aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter clusterFilter = (aod::jcluster::definition == static_cast<int>(clusterDefinition) && aod::jcluster::eta > clusterEtaMin && aod::jcluster::eta < clusterEtaMax && aod::jcluster::phi >= clusterPhiMin && aod::jcluster::phi <= clusterPhiMax && aod::jcluster::energy >= clusterEnergyMin && aod::jcluster::time > clusterTimeMin && aod::jcluster::time < clusterTimeMax && (clusterRejectExotics && aod::jcluster::isExotic != true));
  Preslice<JetTableMCPMatchedJoined> JetMCPPerMcCollision = aod::jet::mcCollisionId;
  PresliceUnsorted<soa::Filtered<aod::JetCollisionsMCD>> CollisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  PresliceUnsorted<o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>> perFoundBC = aod::evsel::foundBCId;

  template <typename T, typename S, typename U>
  bool isAcceptedRecoJet(U const& jet, double& filteredTrackPt, double& filteredClusterPt)
  {
    //Reset filtered pT accumulators (for QA if needed)
    filteredTrackPt = 0.0;
    filteredClusterPt = 0.0;

    // --- Track cuts: ALL tracks must satisfy 0.15 <= pT <= 200 or 150 GeV/c---
    // if (leadingTrackPtMin > kLeadingTrackPtMinThreshold || leadingTrackPtMax < kLeadingTrackPtMaxThreshold) {
      bool hasValidTrack = false;
      for (const auto& constituent : jet.template tracks_as<T>()) {
        const float pt = constituent.pt();
        if ((minTrackPt > kLeadingTrackPtMinThreshold && pt < minTrackPt) ||
            (maxTrackPt < kLeadingTrackPtMaxThreshold && pt > maxTrackPt)) {
          continue; //SKIP this invalid track
        }
        filteredTrackPt += pt; //Accumulate valid track pT
        hasValidTrack = true; // At least one track exists (if needed)
      }
      // Reject jets without valid tracks (edge case)
      if (!hasValidTrack) {
        return false;
      }
    // }

    // --- Cluster cuts: ALL clusters must satisfy min <= pT <= max == 0.3 <= pT <= 250
    // if (leadingClusterPtMin > kLeadingClusterPtMinThreshold || leadingClusterPtMax < kLeadingClusterPtMaxThreshold) {
      bool hasValidCluster = false;
      for (const auto& cluster : jet.template clusters_as<S>()) {
        const double pt = cluster.energy() / std::cosh(cluster.eta());
        if ((minClusterPt > kLeadingClusterPtMinThreshold && pt < minClusterPt) ||
            (maxClusterPt < kLeadingClusterPtMaxThreshold && pt > maxClusterPt)) {
          continue; //SKIP this invalid cluster
        }
        filteredClusterPt += pt;
        hasValidCluster = true; // At least one cluster exists
      }
      // Reject jets without valid clusters (edge case)
      if (!hasValidCluster) {
        return false;
      }
    // }
    return true; //Valid Jet
  } // isAcceptedRecoJet ends

/*  template <typename T, typename U>
  bool isAcceptedPartJet(U const& jet)
  {
    // if (jetAreaFractionMin > kJetAreaFractionMinThreshold) {
    //   if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
    //     return false;
    //   }
    // }
    // track pt Min cut at the Part level: 0 < pT <= 200 or 150 GeV/c
    if (leadingTrackPtMin > kLeadingTrackPtMinThreshold || leadingTrackPtMax < kLeadingTrackPtMaxThreshold) {
      bool hasValidParticle = false;
      for (const auto& constituent : jet.template tracks_as<T>()) {
        const float pt = constituent.pt();
        // Reject if ANY particle fails min or max cut
        if ((leadingTrackPtMin > kLeadingTrackPtMinThreshold && pt < leadingTrackPtMin) ||
            (leadingTrackPtMax < kLeadingTrackPtMaxThreshold && pt > leadingTrackPtMax)) {
          return false;
        }
        hasValidParticle = true; // At least one track exists (if needed)
      }
      // Reject if no particle exist (edge case)
      if (leadingTrackPtMin > kLeadingTrackPtMinThreshold && !hasValidParticle) {
        return false;
      }
    }
    return true;
  }*/

  template <typename T>
  bool isInPhiAcceptance(T const& jet) const
  {
    const double twoPi = 2.0 * M_PI;
    // convert encoded radius to real R (radians)
    const double R = static_cast<double>(jet.r()) / 100.0;

    // emcalPhiMin/emcalPhiMax are configurables for emcal phi edges in radians, e.g. 1.3962634, 3.2836100
    double jetFidPhiMin = emcalPhiMin.value + R;
    double jetFidPhiMax = emcalPhiMax.value - R;

    // normalize to [0, 2pi)
    auto norm = [&](double a) {
      while (a < 0) a += twoPi;
      while (a >= twoPi) a -= twoPi;
      return a;
    };

    double phi = norm(jet.phi());
    jetFidPhiMin = norm(jetFidPhiMin);
    jetFidPhiMax = norm(jetFidPhiMax);

    if (jetFidPhiMin <= jetFidPhiMax) {
      // non-wrap case (EMCal default)
      return (phi >= jetFidPhiMin && phi <= jetFidPhiMax);
    } else {
      // wrap-around case (defensive)
      return (phi >= jetFidPhiMin || phi <= jetFidPhiMax);
    }
  }

  template <typename T>
  void fillJetHistograms(T const& jet, float weight = 1.0)
  {
    // float neutralEnergy = 0.0;
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_full_jet_pt"), jet.pt(), weight);
      registry.fill(HIST("h_full_jet_eta"), jet.eta(), weight);
      registry.fill(HIST("h_full_jet_phi"), jet.phi(), weight);
      registry.fill(HIST("h2_jet_etaphi"), jet.eta(), jet.phi(), weight);

      // Sums for each correction mode
      double neutralEnergy_uncorr = 0.0;
      double neutralEnergy_corr_oneTrack100 = 0.0;
      double neutralEnergy_corr_oneTrack70 = 0.0;
      double neutralEnergy_corr_allTracks100 = 0.0;
      double neutralEnergy_corr_allTracks70 = 0.0;

      for (const auto& cluster : jet.template clusters_as<ClusterWithCorrections>()) {
        registry.fill(HIST("h2_full_jet_neutralconstituents"), jet.pt(), jet.clustersIds().size(), weight);

        // Sum energies for NEF calculation for each correction mode
        neutralEnergy_uncorr += cluster.energy();
        neutralEnergy_corr_oneTrack100 += cluster.energyCorrectedOneTrack1();
        neutralEnergy_corr_oneTrack70 += cluster.energyCorrectedOneTrack2();
        neutralEnergy_corr_allTracks100 += cluster.energyCorrectedAllTracks1();
        neutralEnergy_corr_allTracks70 += cluster.energyCorrectedAllTracks2();

        // neutralEnergy += cluster.energy();
        double clusterpt_uncorr = cluster.energy() / std::cosh(cluster.eta());
        double clusterpt_corr_oneTrack100 = cluster.energyCorrectedOneTrack1() / std::cosh(cluster.eta());
        double clusterpt_corr_oneTrack70 = cluster.energyCorrectedOneTrack2() / std::cosh(cluster.eta());
        double clusterpt_corr_allTracks100 = cluster.energyCorrectedAllTracks1() / std::cosh(cluster.eta());
        double clusterpt_corr_allTracks70 = cluster.energyCorrectedAllTracks2() / std::cosh(cluster.eta());

        registry.fill(HIST("h_full_jet_clusterTime"), cluster.time(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_pt_uncorr"), clusterpt_uncorr, weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_pt_corr_oneTrack100"), clusterpt_corr_oneTrack100, weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_pt_corr_oneTrack70"), clusterpt_corr_oneTrack70, weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_pt_corr_allTracks100"), clusterpt_corr_allTracks100, weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_pt_corr_allTracks70"), clusterpt_corr_allTracks70, weight);

        registry.fill(HIST("h_full_jet_neutralconstituents_eta"), cluster.eta(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_phi"), cluster.phi(), weight);
        registry.fill(HIST("h_full_jet_neutralconstituents_energysum"), neutralEnergy_uncorr, weight);

        registry.fill(HIST("h_full_jet_neutralconstituents_energy_uncorr"), cluster.energy(), weight);

        if (cluster.energyCorrectedOneTrack1()) {
          registry.fill(HIST("h_full_jet_neutralconstituents_energy_corr_oneTrack100"), cluster.energyCorrectedOneTrack1(), weight);
        }
        if (cluster.energyCorrectedOneTrack2()) {
          registry.fill(HIST("h_full_jet_neutralconstituents_energy_corr_oneTrack70"), cluster.energyCorrectedOneTrack2(), weight);
        }
        if (cluster.energyCorrectedAllTracks1()) {
          registry.fill(HIST("h_full_jet_neutralconstituents_energy_corr_allTracks100"), cluster.energyCorrectedAllTracks1(), weight);
        }
        if (cluster.energyCorrectedAllTracks2()) {
          registry.fill(HIST("h_full_jet_neutralconstituents_energy_corr_allTracks70"), cluster.energyCorrectedAllTracks2(), weight);
        }
      }
      // auto nef = neutralEnergy / jet.energy();
      // registry.fill(HIST("h2_full_jet_nef"), jet.pt(), nef, weight);
      auto jetEnergy = jet.energy();

      auto nef_uncorr = neutralEnergy_uncorr / jetEnergy;
      auto nef_corr_oneTrack100 = neutralEnergy_corr_oneTrack100 / jetEnergy;
      auto nef_corr_oneTrack70 = neutralEnergy_corr_oneTrack70 / jetEnergy;
      auto nef_corr_allTracks100 = neutralEnergy_corr_allTracks100 / jetEnergy;
      auto nef_corr_allTracks70 = neutralEnergy_corr_allTracks70 / jetEnergy;

      registry.fill(HIST("h2_full_jet_nef_uncorr"), jet.pt(), nef_uncorr, weight);
      registry.fill(HIST("h2_full_jet_nef_corr_oneTrack100"), jet.pt(), nef_corr_oneTrack100, weight);
      registry.fill(HIST("h2_full_jet_nef_corr_oneTrack70"), jet.pt(), nef_corr_oneTrack70, weight);
      registry.fill(HIST("h2_full_jet_nef_corr_allTracks100"), jet.pt(), nef_corr_allTracks100, weight);
      registry.fill(HIST("h2_full_jet_nef_corr_allTracks70"), jet.pt(), nef_corr_allTracks70, weight);
      double sumtrackE = 0;
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
    if (jet.r() == round(selectedJetsRadius * 100.0f)) {
      float neutralEnergy = 0.0;
      for (const auto& cluster : jet.template clusters_as<ClusterWithCorrections>()) {
        neutralEnergy += cluster.energy();
      }
      auto nef = neutralEnergy / jet.energy();
      registry.fill(HIST("h2_full_jet_nef_rejected"), jet.pt(), nef, weight);
    } // jet.r()
  }

  template <typename T>
  void fillMCPHistograms(T const& jet, float weight = 1.0)
  {
    auto isInFiducial = [&](auto const& jet) { //For QA purposes only
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
        registry.fill(HIST("h2_full_mcpjetOutsideFiducial_pt"), jet.pt(), 1, weight);
        registry.fill(HIST("h_full_mcpjetOutside_eta_part"), jet.eta(), weight);
        registry.fill(HIST("h_full_mcpjetOutside_phi_part"), jet.phi(), weight);
      } else {
        // jet is inside
        registry.fill(HIST("h2_full_mcpjetInsideFiducial_pt"), jet.pt(), 1, weight);
        registry.fill(HIST("h_full_mcpjetInside_eta_part"), jet.eta(), weight);
        registry.fill(HIST("h_full_mcpjetInside_phi_part"), jet.phi(), weight);
      }

      float neutralEnergy = 0.0;
      int neutralconsts = 0;
      int chargedconsts = 0;
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
    double sumclusterE_uncorr = 0.0;
    double sumclusterE_corr_oneTrack100 = 0.0;
    double sumclusterE_corr_oneTrack70 = 0.0;
    double sumclusterE_corr_allTracks100 = 0.0;
    double sumclusterE_corr_allTracks70 = 0.0;

    for (auto const& cluster : clusters) {

      double clusterpt = cluster.energy() / std::cosh(cluster.eta());
      sumclusterE_uncorr += cluster.energy();
      sumclusterE_corr_oneTrack100 += cluster.energyCorrectedOneTrack1();
      sumclusterE_corr_oneTrack70 += cluster.energyCorrectedOneTrack2();
      sumclusterE_corr_allTracks100 += cluster.energyCorrectedAllTracks1();
      sumclusterE_corr_allTracks70 += cluster.energyCorrectedAllTracks2();

      registry.fill(HIST("h_cluster_energy_uncorr"), cluster.energy(), weight);
      if (cluster.energyCorrectedOneTrack1()) {
        registry.fill(HIST("h_cluster_energy_corr_oneTrack100"), cluster.energyCorrectedOneTrack1(), weight);
      }
      if (cluster.energyCorrectedOneTrack2()) {
        registry.fill(HIST("h_cluster_energy_corr_oneTrack70"), cluster.energyCorrectedOneTrack2(), weight);
      }
      if (cluster.energyCorrectedAllTracks1()) {
        registry.fill(HIST("h_cluster_energy_corr_allTracks100"), cluster.energyCorrectedAllTracks1(), weight);
      }
      if (cluster.energyCorrectedAllTracks2()) {
        registry.fill(HIST("h_cluster_energy_corr_allTracks70"), cluster.energyCorrectedAllTracks2(), weight);
      }

      registry.fill(HIST("h_clusterTime"), cluster.time(), weight);
      registry.fill(HIST("h_cluster_pt"), clusterpt, weight);
      registry.fill(HIST("h_cluster_eta"), cluster.eta(), weight);
      registry.fill(HIST("h_cluster_phi"), cluster.phi(), weight);
      registry.fill(HIST("h_cluster_energysum_uncorr"), sumclusterE_uncorr, weight);
      registry.fill(HIST("h_cluster_energysum_corr_oneTrack100"), sumclusterE_corr_oneTrack100, weight);
      registry.fill(HIST("h_cluster_energysum_corr_oneTrack70"), sumclusterE_corr_oneTrack70, weight);
      registry.fill(HIST("h_cluster_energysum_corr_allTracks100"), sumclusterE_corr_allTracks100, weight);
      registry.fill(HIST("h_cluster_energysum_corr_allTracks70"), sumclusterE_corr_allTracks70, weight);
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

  template <typename T, typename U>
  void fillMatchedNoFidHistograms(T const& jetBase, float weight = 1.0)
  {
    if (jetBase.has_matchedJetGeo()) { // geometrical jet matching only needed for pp - here,matching Base(Det.level) with Tag (Part. level) jets
      registry.fill(HIST("h_full_NoFidmatchedmcdjet_tablesize"), jetBase.size(), weight);
      registry.fill(HIST("h_full_NoFidmatchedmcdjet_ntracks"), jetBase.tracksIds().size(), weight);
      registry.fill(HIST("h2_NoFidmatchedjet_etaphiDet"), jetBase.eta(), jetBase.phi(), weight);

      for (const auto& jetTag : jetBase.template matchedJetGeo_as<std::decay_t<U>>()) {
        auto deltaEta = jetBase.eta() - jetTag.eta();
        auto deltaPhi = jetBase.phi() - jetTag.phi();
        auto deltaR = jetutilities::deltaR(jetBase, jetTag);

        registry.fill(HIST("h_full_jet_NoFiddeltaR"), deltaR, weight);
        registry.fill(HIST("h_full_NoFidmatchedmcpjet_tablesize"), jetTag.size(), weight);
        registry.fill(HIST("h_full_NoFidmatchedmcpjet_ntracks"), jetTag.tracksIds().size(), weight);
        registry.fill(HIST("h2_NoFidmatchedjet_etaphiPart"), jetTag.eta(), jetTag.phi(), weight);
        registry.fill(HIST("h2_NoFidmatchedjet_deltaEtaCorr"), jetBase.eta(), jetTag.eta(), weight);
        registry.fill(HIST("h2_NoFidmatchedjet_deltaPhiCorr"), jetBase.phi(), jetTag.phi(), weight);

        // JES for fulljets
        registry.fill(HIST("h2_full_jet_NoFidenergyscaleDet"), jetBase.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h2_full_jet_NoFidenergyscalePart"), jetTag.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h3_full_jet_NoFidenergyscalePart"), jetBase.r() / 100.0, jetTag.pt(), (jetBase.pt() - jetTag.pt()) / jetTag.pt(), weight);
        registry.fill(HIST("h2_full_jet_NoFidetaresolutionPart"), jetTag.pt(), deltaEta / jetTag.eta(), weight);
        registry.fill(HIST("h2_full_jet_NoFidphiresolutionPart"), jetTag.pt(), deltaPhi / jetTag.phi(), weight);

        // Response Matrix
        registry.fill(HIST("h_full_jet_NoFidResponseMatrix"), jetBase.pt(), jetTag.pt(), weight); // MCD vs MCP jet pT
      } // jetTag
    } // jetBase
  }


  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(FullJetSpectra, processDummy, "dummy task", true);

  void processBCs(o2::soa::Join<o2::aod::BCs, o2::aod::BcSels> const& bcs, o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels> const& collisions)
  {
    if (bcs.size() == 0) {
      return;
    }
    for (auto bc : bcs) {
      registry.fill(HIST("hBCCounter"), 0.5); // All BC
      if (bc.selection_bit(aod::evsel::EventSelectionFlags::kIsTriggerTVX)) {
        registry.fill(HIST("hBCCounter"), 1.5); // BC+TVX
        if (bc.selection_bit(aod::evsel::EventSelectionFlags::kNoTimeFrameBorder)) {
          registry.fill(HIST("hBCCounter"), 2.5); // BC+TVX+NoTFB
          if (bc.selection_bit(aod::evsel::EventSelectionFlags::kNoITSROFrameBorder)) {
            registry.fill(HIST("hBCCounter"), 3.5); // BC+TVX+NoTFB+NoITSROFB ----> this goes to Lumi i.e. hLumiAfterBCcuts in eventSelection task
          }
        }
      }
      auto collisionsInBC = collisions.sliceBy(perFoundBC, bc.globalIndex());
      for (auto collision : collisionsInBC) {
        registry.fill(HIST("hBCCounter"), 4.5); // CollinBC
        if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
          registry.fill(HIST("hBCCounter"), 5.5); // CollinBC+TVX
          if (collision.sel8()) {
            registry.fill(HIST("hBCCounter"), 6.5); // CollinBC+TVX+sel8
            if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
              registry.fill(HIST("hBCCounter"), 7.5); // CollinBC+TVX+sel8Full
              if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
                registry.fill(HIST("hBCCounter"), 8.5); // CollinBC+TVX+sel8Full+GoodZvtx
                if (std::fabs(collision.posZ()) < vertexZCut) {
                  registry.fill(HIST("hBCCounter"), 9.5); // CollinBC+TVX+sel8Full+VtxZ+GoodZvtx ----> this goes to my analysis task for jet events selection
                }
              }
            }
          }
        }
      } // collision loop
    } // bc loop
  }
  PROCESS_SWITCH(FullJetSpectra, processBCs, "BCs for 0 vertex QA", false);

  void processJetsData(soa::Filtered<EMCCollisionsData>::iterator const& collision, FullJetTableDataJoined const& jets,
                       aod::JetTracks const&, ClusterWithCorrections const&)
  {
    bool eventAccepted = false;

    registry.fill(HIST("hDetcollisionCounter"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hDetcollisionCounter"), 1.5); // DetCollWithVertexZ

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
        registry.fill(HIST("hDetcollisionCounter"), 5.5); // EMCreadoutDetEventsWithkTVXinEMC
        eventAccepted = true;
      }
    }

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        {
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
      // if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax) {
      //   continue;
      // }
      if (!isInPhiAcceptance(jet)) {  // Using the new phi acceptance function
          continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsData, "Full Jets Data", false);

  void processJetsTriggeredData(soa::Filtered<EMCCollisionsTriggeredData>::iterator const& collision, FullJetTableDataJoined const& jets,
                                aod::JetTracks const&, ClusterWithCorrections const&, aod::JBCs const&)
  {
    bool eventAccepted = false;

    registry.fill(HIST("hDetTrigcollisionCounter"), 0.5); // allDetTrigColl

    // Get BC info associated with the collision before applying any event selections
    auto bc = collision.bc_as<aod::JBCs>();
    // Initialize CCDB objects using the BC info
    initCCDB(bc);
    // If SoftwareTriggerSelection (i.e. skimming) is enabled, skip this event unless it passes Zorro selection
    if (doSoftwareTriggerSelection && !zorro.isSelected(bc.globalBC())) {
      return;
    }
    registry.fill(HIST("hDetTrigcollisionCounter"), 1.5); // DetTrigCollAfterZorroSelection

    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hDetTrigcollisionCounter"), 2.5); // DetTrigCollWithVertexZ

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || !jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 3.5); // EventsNotSatisfyingEvent+TriggerSelection
      return;
    }
    //- should this kTVX HW trigger be still in place? - Removing it for now; probably not needed if we are only interested in SW triggers
    /*if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
      // eventAccepted = true;
      registry.fill(HIST("hDetTrigcollisionCounter"), 4.5); // EMCreadoutDetTrigEventsWithkTVXinEMC
    }*/
    // split event selections based on selected triggers -
    //  make sure there're no trigger overlaps: when analysing JetFullHighPt-> check no JetFullLowPt and kTVXinEMC are fired
    //  when analysing JetFullLowPt, check kTVXinEMC isn't fired!
    // apply exclusive trigger bit selections for event selection
    //  use a flag and fill QA for trigger overlap -
    //  - how often you reject a higher pT trig because lower trigs were fired : 5 cases -> 2D hist as a funtn of jet pT
    //  - check how often the ChJet Trigs were fired for every fullJetTrig fired.(don't reject these events but only for QA)

    // Get trigger status
    bool hasFullJetHighPt = jetderiveddatautilities::selectTrigger(collision, jetderiveddatautilities::JTrigSel::JetFullHighPt);
    bool hasFullJetLowPt = jetderiveddatautilities::selectTrigger(collision, jetderiveddatautilities::JTrigSel::JetFullLowPt);
    bool hasMB = jetderiveddatautilities::selectTrigger(collision, jetderiveddatautilities::JTrigSel::EMCALReadout);

    //*****Step 1: Check for pure triggers (NO overlaps)*****

    // Case 1: hasFullJetHighPt && !hasFullJetLowPt && !hasMB : Pure FullJetHighPt
    //  i.e. for every JetFullHighPt trig that was fired, check the low triggers weren't fired
    if (hasFullJetHighPt && !hasFullJetLowPt) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 4.5); // FullJetHighPt+FullJetLowPt
      eventAccepted = true;
    }
    // Case 2: hasFullJetLowPt && !hasMB : Pure FullJetLowPt
    //  i.e. for every hasFullJetLowPt trig that was fired, check the MB trig wasn't fired
    if (hasFullJetLowPt) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 5.5); // FullJetLowPt
      eventAccepted = true;
    }
    // // Case 3: hasMB && !hasFullJetLowPt && !hasFullJetHighPt : Pure MB
    // //  i.e. for every MB trig that was fired, check the higher trigs weren't fired
    // if (hasMB && !hasFullJetLowPt && !hasFullJetHighPt) {
    //   registry.fill(HIST("hDetTrigcollisionCounter"), 7.5); // OnlyMB
    // }

    //*****Step 2: Check for trigger overlap cases (for QA):*****

    if (hasFullJetHighPt && hasFullJetLowPt) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 6.5); // FullJetHighPt+FullJetLowPt
      eventAccepted = true;
    }
    if (hasFullJetHighPt && hasMB) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 7.5); // FullJetHighPt+MB
      eventAccepted = true;
    }
    if (hasFullJetLowPt && hasMB) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 8.5); // FullJetLowPt+MB
      eventAccepted = true;
    }

    //*****Step 3: Reject ALL overlapping events by applying EXCLUSIVE Trigger Selections *****
    // Skip further processing if ANY overlaps exist
    // if ((hasFullJetHighPt && (hasFullJetLowPt || hasMB)) || (hasFullJetLowPt && hasMB)) {
    //   registry.fill(HIST("hDetTrigcollisionCounter"), 11.5); // AllRejectedTrigOverlaps
    //   return;
    // }
    if ((hasFullJetHighPt && hasFullJetLowPt)) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 9.5); // AllRejectedTrigOverlaps
      return;
    }
    registry.fill(HIST("hDetTrigcollisionCounter"), 10.5); // EMCAcceptedDetTrigCollAfterTrigOverlapChecks

    if (!eventAccepted) {
      for (auto const& jet : jets) {
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          fillRejectedJetHistograms(jet, 1.0);
        }
      }
      registry.fill(HIST("hDetTrigcollisionCounter"), 11.5); // AllRejectedDetTrigEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hDetTrigcollisionCounter"), 12.5); // EMCAcceptedDetTrigColl

    if (jetderiveddatautilities::selectTrigger(collision, jetderiveddatautilities::JTrigSel::JetChLowPt)) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 13.5); // EMCAcceptedDetTrigCollWithLowChargedJetTriggers
      eventAccepted = true;
    }
    if (jetderiveddatautilities::selectTrigger(collision, jetderiveddatautilities::JTrigSel::JetChHighPt)) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 14.5); // EMCAcceptedDetTrigCollWithHighChargedJetTriggers
      eventAccepted = true;
    }

    if (jetderiveddatautilities::selectTrigger(collision, jetderiveddatautilities::JTrigSel::JetFullLowPt) && jetderiveddatautilities::selectTrigger(collision, jetderiveddatautilities::JTrigSel::JetFullHighPt)) {
      registry.fill(HIST("hDetTrigcollisionCounter"), 15.5); // EMCAcceptedDetTrigCollWithLow+HighFullJetTriggers
      eventAccepted = true;
    }
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isInPhiAcceptance(jet)) {  // Using the new phi acceptance function
          continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsTriggeredData, "Full Jets Triggered Data", false);

  void processJetsMCD(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDJoined const& jets,
                      aod::JetTracks const&, ClusterWithCorrections const&)
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
    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hDetcollisionCounter"), 3.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hDetcollisionCounter"), 4.5); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
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
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
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
      if (!isInPhiAcceptance(jet)) {  // Using the new phi acceptance function
          continue;
      }
      fillJetHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCD, "Full Jets at Detector Level", false);

  void processJetsMCDWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDJoined const& jets, aod::JMcCollisions const&,
                              aod::JetTracks const&, ClusterWithCorrections const&)
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
    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hDetcollisionCounter"), 3.5, collision.mcCollision().weight()); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hDetcollisionCounter"), 4.5, collision.mcCollision().weight()); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
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
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
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
      if (!isInPhiAcceptance(jet)) {  // Using the new phi acceptance function
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

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hPartcollisionCounter"), 2.5); // RejectedPartCollWithOutliers
        return;
      }
    }

    if (doMBGapTrigger && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      // Fill rejected MB events;
      registry.fill(HIST("hPartcollisionCounter"), 3.5); // MBRejectedPartEvents
      return;
    }

    auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, mccollision.globalIndex());
    registry.fill(HIST("hRecoMatchesPerMcCollision"), collisionspermcpjet.size()); // for split vertices QA

    if (collisionspermcpjet.size() == 0 || collisionspermcpjet.size() < 1) {
      registry.fill(HIST("hPartcollisionCounter"), 4.5); // RejectedPartCollForDetCollWithSize0or<1
      return;
    }
    registry.fill(HIST("hPartcollisionCounter"), 5.5); // AcceptedPartCollWithSize>=1

    for (auto const& collision : collisionspermcpjet) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        continue;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          if (collision.alias_bit(kTVXinEMC)) {
            eventAccepted = true;
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
    registry.fill(HIST("hPartcollisionCounter"), 8.5); // EMCAcceptedPartColl

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetNoFidPartEtaMin, jetNoFidPartEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      fillMCPHistograms(jet);
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCP, "Full Jets at Particle Level", false);

  void processJetsMCPWeighted(aod::JetMcCollision const& mccollision, JetTableMCPJoined const& jets, aod::JetParticles const&, soa::SmallGroups<EMCCollisionsMCD> const& collisions)
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

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hPartcollisionCounter"), 2.5, mccollision.weight()); // RejectedPartCollWithOutliers
        return;
      }
    }

    if (doMBGapTrigger && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      // Fill rejected MB events
      registry.fill(HIST("hPartcollisionCounter"), 3.5, mccollision.weight()); // MBRejectedPartEvents
      return;
    }

    auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, mccollision.globalIndex());
    registry.fill(HIST("hRecoMatchesPerMcCollision"), collisionspermcpjet.size(), mccollision.weight()); // for split vertices QA

    if (collisionspermcpjet.size() == 0 || collisionspermcpjet.size() < 1) {
      registry.fill(HIST("hPartcollisionCounter"), 4.5, mccollision.weight()); // RejectedPartCollForDetCollWithSize0or<1
      return;
    }
    registry.fill(HIST("hPartcollisionCounter"), 5.5, mccollision.weight()); // AcceptedPartCollWithSize>=1

    for (auto const& collision : collisionspermcpjet) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        continue;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          if (collision.alias_bit(kTVXinEMC)) {
            eventAccepted = true;
            registry.fill(HIST("hPartcollisionCounter"), 6.5, mccollision.weight()); // EMCreadoutDetEventsWithkTVXinEMC
          }
        }
      } else {
        if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
          registry.fill(HIST("hPartcollisionCounter"), 6.5, mccollision.weight()); // EMCreadoutDetEventsWithkTVXinEMC
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
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetNoFidPartEtaMin, jetNoFidPartEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (doMBGapTrigger && mccollision.weight() == 1) {
        continue;
      }
      fillMCPHistograms(jet, mccollision.weight());
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCPWeighted, "Full Jets at Particle Level on weighted events", false);

  void processJetsMCPMCDMatched(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDMatchedJoined const& mcdjets, JetTableMCPMatchedJoined const& mcpjets, aod::JMcCollisions const&,
                                aod::JetTracks const&, ClusterWithCorrections const&, aod::JetParticles const&)
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

    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hMatchedcollisionCounter"), 4.5); // EMCMBRejectedDetColl
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hMatchedcollisionCounter"), 5.5); // EventsNotSatisfyingEventSelection
      return;
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
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

      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) ||
          !isInPhiAcceptance(mcdjet)) {
        fakeMcdJet++;
        registry.fill(HIST("h2_full_fakemcdjets"), mcdjet.pt(), fakeMcdJet, 1.0);
        continue;
      }
      // Check if MCD jet is within the EMCAL fiducial region; if not then flag it as a fake jet
      // if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax || mcdjet.eta() < jetEtaMin || mcdjet.eta() > jetEtaMax) {
      //   fakeMcdJet++;
      //   registry.fill(HIST("h2_full_fakemcdjets"), mcdjet.pt(), fakeMcdJet, 1.0);
      //   continue;
      // }

      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedJoined>()) {
        if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetPartEtaMin, jetPartEtaMax, trackEtaMin, trackEtaMax) ||
            !isInPhiAcceptance(mcpjet)) {
          fakeMcpJet++;
          registry.fill(HIST("h2_full_fakemcpjets"), mcpjet.pt(), fakeMcpJet, 1.0);
          continue;
        } else {
          fillMatchedHistograms<JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined>(mcdjet);
        }
      } // mcpjet loop
    } // mcdjet loop
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsMCPMCDMatched, "Full Jet finder MCP matched to MCD", false);

  void processJetsNoFidMCPMCDMatchedWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDMatchedJoined const& mcdjets, JetTableMCPMatchedJoined const& mcpjets, aod::JMcCollisions const&,
                                        aod::JetTracks const&, ClusterWithCorrections const&, aod::JetParticles const&)
  {
    bool eventAccepted = false;
    int fakeMcdJet = 0;
    int fakeMcpJet = 0;
    int NPartJetFid = 0;
    int allMatchedPartJets = 0;
    float eventWeight = collision.mcCollision().weight();
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
    const auto mcpJetsPerMcCollision = mcpjets.sliceBy(JetMCPPerMcCollision, collision.mcCollisionId());

    registry.fill(HIST("hMatchedNoFidcollisionCounter"), 0.5, eventWeight); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {                    // making double sure this condition is satisfied
      return;
    }
    registry.fill(HIST("hMatchedNoFidcollisionCounter"), 1.5, eventWeight); // DetCollWithVertexZ

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& mcdjet : mcdjets) {
      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hMatchedNoFidcollisionCounter"), 2.5, eventWeight); // RejectedDetCollWithOutliers
        return;
      }
    }

    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hMatchedNoFidcollisionCounter"), 4.5, eventWeight); // EMCMBRejectedDetColl
      return;
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hMatchedNoFidcollisionCounter"), 5.5, eventWeight); // EventsNotSatisfyingEventSelection
      return;
    }

    for (auto const& mcpjet : mcpJetsPerMcCollision) {
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) { // outlier rejection for MCP: Should I remove this cut as I'm already doing MC outlier rejection @L1071?
        return;
      }
    }
    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
          registry.fill(HIST("hMatchedNoFidcollisionCounter"), 6.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hMatchedNoFidcollisionCounter"), 6.5, eventWeight); // EMCreadoutDetJJEventsWithkTVXinEMC
      }
    }
    if (!eventAccepted) {
      registry.fill(HIST("hMatchedNoFidcollisionCounter"), 7.5, eventWeight); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hMatchedNoFidcollisionCounter"), 8.5, eventWeight); // EMCAcceptedDetColl

    for (const auto& mcdjet : mcdjets) {

      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) ||
          !isInPhiAcceptance(mcdjet)) {
        fakeMcdJet++;
        registry.fill(HIST("h2_full_NoFidfakemcdjets"), mcdjet.pt(), fakeMcdJet, eventWeight);
        continue;
      }

      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedJoined>()) {
        allMatchedPartJets++;
        registry.fill(HIST("h_allMatchedNoFidPartJetsPt"), mcpjet.pt(), eventWeight);

        //Not applying any emcal fiducial cuts in eta and phi on MCP jets when matching.
        //Keeping jet eta here open = |0.9| and no cut in phi at all.
        if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetNoFidPartEtaMin, jetNoFidPartEtaMax, trackEtaMin, trackEtaMax)) {
          fakeMcpJet++;
          registry.fill(HIST("h2_full_NoFidfakemcpjets"), mcpjet.pt(), fakeMcpJet, eventWeight);
          continue;
        } else {
          NPartJetFid++;
          // Fill matched histograms (including Response Matrix) for valid MCD-MCP pairs
          fillMatchedNoFidHistograms<JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined>(mcdjet, eventWeight);
          registry.fill(HIST("h2_full_NoFidmatchedmcpjet_pt"), mcpjet.pt(), NPartJetFid, eventWeight);
          registry.fill(HIST("h_full_NoFidmatchedmcpjet_eta"), mcpjet.eta(), eventWeight);
          registry.fill(HIST("h_full_NoFidmatchedmcpjet_phi"), mcpjet.phi(), eventWeight);
        }
      } // mcpjet
    // Fill the total matched particle jets histogram after processing all MCP jets for the MCD jet
    registry.fill(HIST("h_allMatchedNoFidPartJetsCounter"), allMatchedPartJets, eventWeight);
    } // mcdjet
  }
  PROCESS_SWITCH(FullJetSpectra, processJetsNoFidMCPMCDMatchedWeighted, "Full Jet finder No Fid MCP matched to MCD on weighted events", false);


  void processJetsMCPMCDMatchedWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDMatchedJoined const& mcdjets, JetTableMCPMatchedJoined const& mcpjets, aod::JMcCollisions const&,
                                        aod::JetTracks const&, ClusterWithCorrections const&, aod::JetParticles const&)
  {
    bool eventAccepted = false;
    int fakeMcdJet = 0;
    int fakeMcpJet = 0;
    int NPartJetFid = 0;
    int allMatchedPartJets = 0;
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

    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
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
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
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

      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax) ||
          !isInPhiAcceptance(mcdjet)) {
        fakeMcdJet++;
        registry.fill(HIST("h2_full_fakemcdjets"), mcdjet.pt(), fakeMcdJet, eventWeight);
        continue;
      }

      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<JetTableMCPMatchedJoined>()) {
        allMatchedPartJets++;
        registry.fill(HIST("h_allMatchedPartJetsPt"), mcpjet.pt(), eventWeight);

        if (!jetfindingutilities::isInEtaAcceptance(mcpjet, jetPartEtaMin, jetPartEtaMax, trackEtaMin, trackEtaMax) ||
            !isInPhiAcceptance(mcpjet)) {
          fakeMcpJet++;
          registry.fill(HIST("h2_full_fakemcpjets"), mcpjet.pt(), fakeMcpJet, eventWeight);
          continue;
        } else {
          NPartJetFid++;
          // Fill matched histograms (including Response Matrix) for valid MCD-MCP pairs
          fillMatchedHistograms<JetTableMCDMatchedJoined::iterator, JetTableMCPMatchedJoined>(mcdjet, eventWeight);
          // If both MCD-MCP matched jet pairs are within the EMCAL fiducial region, fill these kinematic histos
          registry.fill(HIST("h2_full_matchedmcpjet_pt"), mcpjet.pt(), NPartJetFid, eventWeight);
          registry.fill(HIST("h_full_matchedmcpjet_eta"), mcpjet.eta(), eventWeight);
          registry.fill(HIST("h_full_matchedmcpjet_phi"), mcpjet.phi(), eventWeight);
        }
      } // mcpjet
    // Fill the total matched particle jets histogram after processing all MCP jets for the MCD jet
    registry.fill(HIST("h_allMatchedPartJetsCounter"), allMatchedPartJets, eventWeight);
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
  void processDataTracks(soa::Filtered<EMCCollisionsData>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks,
                         ClusterWithCorrections const& clusters) // replaced "soa::Filtered<aod::JetClusters>" with ClusterWithCorrections to include the hadcorr tables
  {
    bool eventAccepted = false;

    registry.fill(HIST("hCollisionsUnweighted"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hCollisionsUnweighted"), 1.5); // DetCollWithVertexZ

    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
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
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
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
    }
    fillTrackHistograms(tracks, clusters, 1.0);
    registry.fill(HIST("hCollisionsUnweighted"), 7.5); // EMCAcceptedCollAfterTrackSel
  }
  PROCESS_SWITCH(FullJetSpectra, processDataTracks, "Full Jet tracks for Data", false);

  void processMCTracks(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks, soa::Filtered<ClusterWithCorrections> const& clusters)
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
    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hCollisionsUnweighted"), 2.5); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hCollisionsUnweighted"), 3.5); // EventsNotSatisfyingEventSelection
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
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
    }
    fillTrackHistograms(tracks, clusters, 1.0);
    registry.fill(HIST("hCollisionsUnweighted"), 7.5); // EMCAcceptedCollAfterTrackSel
  }
  PROCESS_SWITCH(FullJetSpectra, processMCTracks, "Full Jet tracks for MC", false);

  void processTracksWeighted(soa::Filtered<EMCCollisionsMCD>::iterator const& collision,
                             aod::JMcCollisions const&,
                             soa::Filtered<aod::JetTracks> const& tracks,
                             soa::Filtered<ClusterWithCorrections> const& clusters)
  {
    bool eventAccepted = false;
    float eventWeight = collision.mcCollision().weight();
    float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));

    registry.fill(HIST("hCollisionsWeighted"), 0.5, eventWeight); // AllWeightedDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hCollisionsWeighted"), 1.5, eventWeight); // WeightedCollWithVertexZ

    if (pTHat < pTHatAbsoluteMin) { // Track outlier rejection: should this be for every track iteration or for every collision?
      return;
    }

    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
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
        fillTrackHistograms(tracks, clusters, eventWeight);
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
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
    }
    fillTrackHistograms(tracks, clusters, eventWeight);
    registry.fill(HIST("hCollisionsWeighted"), 7.5, eventWeight); // EMCAcceptedWeightedCollAfterTrackSel
  }
  PROCESS_SWITCH(FullJetSpectra, processTracksWeighted, "Full Jet tracks weighted", false);

  void processMBCollisionsDATAWithMultiplicity(soa::Filtered<EMCCollisionsData>::iterator const& collision,
                                               FullJetTableDataJoined const& jets, aod::JetTracks const& /*tracks*/, ClusterWithCorrections const& /*clusters*/)
  {
    bool eventAccepted = false;

    registry.fill(HIST("hEventmultiplicityCounter"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 1.5); // DetCollWithVertexZ

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hEventmultiplicityCounter"), 2.5); // EventsNotSatisfyingEventSelection
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        if (collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
          registry.fill(HIST("hEventmultiplicityCounter"), 3.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hEventmultiplicityCounter"), 3.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hEventmultiplicityCounter"), 4.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 5.5); // EMCAcceptedDetColl

    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multFT0M());

    // Jet processing - NEW IMPLEMENTATION
    std::vector<typename std::decay_t<decltype(jets)>::iterator> selectedJets;
    // static int eventCounter = 0;
    int nJetsThisEvent = 0;
    // Debug output
    // std::cout << "===== Event " << ++eventCounter << " (Collision ID: " << collision.globalIndex() << ") =====" << std::endl;

    // Verify jet-collision association
    for (auto const& jet : jets) {
      // Declare variables to store filtered track/cluster pT
      double filteredTrackPt = 0.0;
      double filteredClusterPt = 0.0;
      if (jet.collisionId() != collision.globalIndex()) {
        LOGF(warn, "Jet with pT %.2f belongs to collision %d but processing collision %d", jet.pt(), jet.collisionId(), collision.globalIndex());
        continue;
      }

      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        continue;
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax)
        continue;
      if (!isAcceptedRecoJet<aod::JetTracks, ClusterWithCorrections>(jet, filteredTrackPt, filteredClusterPt))
        continue;

      selectedJets.push_back(jet);
      nJetsThisEvent++;
      // std::cout << "Selected jet pT: " << jet.pt() << " (collision ID: " << jet.collisionId() << ")" << std::endl;
    }
    // 1. Sort selected jets by pT before processing
    std::sort(selectedJets.begin(), selectedJets.end(),
              [](auto const& a, auto const& b) { return (*a).pt() > (*b).pt(); });

    // std::cout << "Number of selected jets: " << selectedJets.size() << std::endl;

    // // 2. Reset counters for each event
    // int numberOfChargedParticles = 0;
    // int totalNumberOfChargedParticles = 0;
    // int leadingJetCount = 0;
    //
    // std::cout << "Number of selected jets per event: " << selectedJets.size() << std::endl;
    // for (const auto& jetIter : selectedJets) {
    //   std::cout << "Jet pT of selectedJets: " << (*jetIter).pt() << std::endl;
    // }
    // //Checking Event Counter
    // static int eventCounter = 0;
    // std::cout << "===== Event " << ++eventCounter << " =====" << std::endl;
    // std::cout << "******************************************** " << std::endl;
    if (selectedJets.size() == 0) { // no jets = no leading jet
      return;
    } else {
      // Jet multiplicity per event
      registry.fill(HIST("h_all_fulljet_Njets"), selectedJets.size(), 1.0);

      // Select Leading Jet for N_ch calculation (for every leading jet that is found). There's always one leading jet per event!
      auto const& leadingJet = *selectedJets[0];
      auto const& leadingJetPt = leadingJet.pt(); // jet pT distribution of the leading jet
      // std::cout << "Leading Jet pT: " << leadingJetPt << std::endl;
      registry.fill(HIST("h_Leading_full_jet_pt"), leadingJetPt, 1.0);
      registry.fill(HIST("h2_full_jet_leadingJetPt_vs_counts"), leadingJetPt, nJetsThisEvent, 1.0);
    }

    if (selectedJets.size() > 1) {
      auto const& subLeadingJet = *selectedJets[1];
      auto const& subLeadingJetPt = subLeadingJet.pt(); // jet pT distribution of the subleading jet i.e. 2nd leading jet
      registry.fill(HIST("h_SubLeading_full_jet_pt"), subLeadingJetPt, 1.0);
      registry.fill(HIST("h2_full_jet_subLeadingJetPt_vs_counts"), subLeadingJetPt, nJetsThisEvent, 1.0);
    }
    // Process ALL selected jets (not just leading)
    for (size_t i = 0; i < selectedJets.size(); i++) {
      auto const& jet = *selectedJets[i];
      float jetPt = jet.pt();
      bool isLeading = (i == 0);
      bool isSubLeading = (i == 1 && selectedJets.size() > 1); // first sub-leading jet

      // Count charged particles(NCh) for this jet
      int numberOfChargedParticles = 0;
      for (const auto& jettrack : jet.tracks_as<aod::JetTracks>()) {
        if (jetderiveddatautilities::selectTrack(jettrack, trackSelection)) {
          numberOfChargedParticles++;
        } else {
          continue;
        }
      }

      // Calculate neutral energy fraction for this jet
      // float neutralEnergy = 0.0;
      double neutralEnergy_uncorr = 0.0;
      double neutralEnergy_corr_oneTrack100 = 0.0;
      double neutralEnergy_corr_oneTrack70 = 0.0;
      double neutralEnergy_corr_allTracks100 = 0.0;
      double neutralEnergy_corr_allTracks70 = 0.0;

      for (const auto& jetcluster : jet.clusters_as<ClusterWithCorrections>()) {
        // Sum energies for NEF calculation for each correction mode
        neutralEnergy_uncorr += jetcluster.energy();
        neutralEnergy_corr_oneTrack100 += jetcluster.energyCorrectedOneTrack1();
        neutralEnergy_corr_oneTrack70 += jetcluster.energyCorrectedOneTrack2();
        neutralEnergy_corr_allTracks100 += jetcluster.energyCorrectedAllTracks1();
        neutralEnergy_corr_allTracks70 += jetcluster.energyCorrectedAllTracks2();
      }
      auto jetEnergy = jet.energy();

      auto nef_uncorr = neutralEnergy_uncorr / jetEnergy;
      auto nef_corr_oneTrack100 = neutralEnergy_corr_oneTrack100 / jetEnergy;
      auto nef_corr_oneTrack70 = neutralEnergy_corr_oneTrack70 / jetEnergy;
      auto nef_corr_allTracks100 = neutralEnergy_corr_allTracks100 / jetEnergy;
      auto nef_corr_allTracks70 = neutralEnergy_corr_allTracks70 / jetEnergy;

      // CASE 1: Fill histograms for ALL selected jets
      registry.fill(HIST("h_all_fulljet_pt"), jetPt, 1.0);
      registry.fill(HIST("h_all_fulljet_Nch"), numberOfChargedParticles, 1.0);
      registry.fill(HIST("h_all_fulljet_NEF_uncorr"), nef_uncorr, 1.0);
      registry.fill(HIST("h_all_fulljet_NEF_corr_oneTrack100"), nef_corr_oneTrack100, 1.0);
      registry.fill(HIST("h_all_fulljet_NEF_corr_oneTrack70"), nef_corr_oneTrack70, 1.0);
      registry.fill(HIST("h_all_fulljet_NEF_corr_allTracks100"), nef_corr_allTracks100, 1.0);
      registry.fill(HIST("h_all_fulljet_NEF_corr_allTracks70"), nef_corr_allTracks70, 1.0);

      registry.fill(HIST("h2_all_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), 1.0);
      registry.fill(HIST("h2_all_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, 1.0);
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef_uncorr"), jetPt, collision.multFT0M(), nef_uncorr, 1.0);
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef_corr_oneTrack100"), jetPt, collision.multFT0M(), nef_corr_oneTrack100, 1.0);
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef_corr_oneTrack70"), jetPt, collision.multFT0M(), nef_corr_oneTrack70, 1.0);
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef_corr_allTracks100"), jetPt, collision.multFT0M(), nef_corr_allTracks100, 1.0);
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef_corr_allTracks70"), jetPt, collision.multFT0M(), nef_corr_allTracks70, 1.0);

      // CASE 2: Additional leading jet processing
      if (isLeading) {
        registry.fill(HIST("h_leading_fulljet_pt"), jetPt, 1.0);
        registry.fill(HIST("h_leading_fulljet_Nch"), numberOfChargedParticles, 1.0);
        registry.fill(HIST("h_leading_fulljet_NEF"), nef_uncorr, 1.0);
        registry.fill(HIST("h2_leading_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), 1.0);
        registry.fill(HIST("h2_leading_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, 1.0);
        registry.fill(HIST("h3_leading_fulljet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef_uncorr, 1.0);
      }

      // CASE 3: Additional first sub-leading jet processing
      if (isSubLeading) {
        registry.fill(HIST("h_subleading_fulljet_pt"), jetPt, 1.0);
        registry.fill(HIST("h_subleading_fulljet_Nch"), numberOfChargedParticles, 1.0);
        registry.fill(HIST("h_subleading_fulljet_NEF"), nef_uncorr, 1.0);
        registry.fill(HIST("h2_subleading_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), 1.0);
        registry.fill(HIST("h2_subleading_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, 1.0);
        registry.fill(HIST("h3_subleading_fulljet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef_uncorr, 1.0);
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processMBCollisionsDATAWithMultiplicity, "MB DATA Collisions for Full Jets Multiplicity Studies", false);

  void processMBMCDCollisionsWithMultiplicity(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDJoined const& mcdjets, aod::JMcCollisions const&, aod::JetTracks const& /*tracks*/, ClusterWithCorrections const& /*clusters*/)
  {
    bool eventAccepted = false;

    registry.fill(HIST("hEventmultiplicityCounter"), 0.5); // allDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 1.5); // DetCollWithVertexZ

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hEventmultiplicityCounter"), 2.5); // EventsNotSatisfyingEventSelection
      return;
    }

    if (doEMCALEventWorkaround) {
      if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
        eventAccepted = true;
        // fillTrackHistograms(tracks, clusters, 1.0);
        if (collision.alias_bit(kTVXinEMC)) {
          registry.fill(HIST("hEventmultiplicityCounter"), 3.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hEventmultiplicityCounter"), 3.5); // EMCreadoutDetEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hEventmultiplicityCounter"), 4.5); // AllRejectedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 5.5); // EMCAcceptedDetColl

    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multFT0M());

    std::vector<typename std::decay_t<decltype(mcdjets)>::iterator> selectedJets;
    // static int eventCounter = 0;
    int nJetsThisEvent = 0;
    // Debug output
    // std::cout << "===== Event " << ++eventCounter << " (Collision ID: " << collision.globalIndex() << ") =====" << std::endl;

    // Verify jet-collision association
    for (auto const& mcdjet : mcdjets) {
      // Declare variables to store filtered track/cluster pT
      double filteredTrackPt = 0.0;
      double filteredClusterPt = 0.0;
      if (mcdjet.collisionId() != collision.globalIndex()) {
        LOGF(warn, "Jet with pT %.2f belongs to collision %d but processing collision %d", mcdjet.pt(), mcdjet.collisionId(), collision.globalIndex());
        continue;
      }

      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        continue;
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax)
        continue;
      if (!isAcceptedRecoJet<aod::JetTracks, ClusterWithCorrections>(mcdjet, filteredTrackPt, filteredClusterPt))
        continue;

      selectedJets.push_back(mcdjet);
      nJetsThisEvent++;
      // std::cout << "Selected jet pT: " << jet.pt() << " (collision ID: " << jet.collisionId() << ")" << std::endl;
    }
    // 1. Sort selected jets by pT before processing
    std::sort(selectedJets.begin(), selectedJets.end(),
              [](auto const& a, auto const& b) { return (*a).pt() > (*b).pt(); });

    if (selectedJets.size() == 0) { // no jets = no leading jet
      return;
    } else {
      // Jet multiplicity per event
      registry.fill(HIST("h_all_fulljet_Njets"), selectedJets.size(), 1.0);

      // Select Leading Jet for N_ch calculation (for every leading jet that is found). There's always one leading jet per event!
      auto const& leadingJet = *selectedJets[0];
      auto const& leadingJetPt = leadingJet.pt(); // jet pT distribution of the leading jet
      // std::cout << "Leading Jet pT: " << leadingJetPt << std::endl;
      registry.fill(HIST("h_Leading_full_jet_pt"), leadingJetPt, 1.0);
      registry.fill(HIST("h2_full_jet_leadingJetPt_vs_counts"), leadingJetPt, nJetsThisEvent, 1.0);
    }

    if (selectedJets.size() > 1) {
      auto const& subLeadingJet = *selectedJets[1];
      auto const& subLeadingJetPt = subLeadingJet.pt(); // jet pT distribution of the subleading jet i.e. 2nd leading jet
      registry.fill(HIST("h_SubLeading_full_jet_pt"), subLeadingJetPt, 1.0);
      registry.fill(HIST("h2_full_jet_subLeadingJetPt_vs_counts"), subLeadingJetPt, nJetsThisEvent, 1.0);
    }
    // Process ALL selected jets (not just leading)
    for (size_t i = 0; i < selectedJets.size(); i++) {
      auto const& jet = *selectedJets[i];
      float jetPt = jet.pt();
      bool isLeading = (i == 0);
      bool isSubLeading = (i == 1 && selectedJets.size() > 1); // first sub-leading jet

      // Count charged particles(NCh) for this jet
      int numberOfChargedParticles = 0;
      for (const auto& jettrack : jet.tracks_as<aod::JetTracks>()) {
        if (jetderiveddatautilities::selectTrack(jettrack, trackSelection)) {
          numberOfChargedParticles++;
        } else {
          continue;
        }
      }

      // Calculate neutral energy fraction for this jet
      float neutralEnergy = 0.0;
      for (const auto& jetcluster : jet.clusters_as<ClusterWithCorrections>()) {
        neutralEnergy += jetcluster.energy();
      }
      float nef = neutralEnergy / jet.energy();

      // CASE 1: Fill histograms for ALL selected jets
      registry.fill(HIST("h_all_fulljet_pt"), jetPt, 1.0);
      registry.fill(HIST("h_all_fulljet_Nch"), numberOfChargedParticles, 1.0);
      registry.fill(HIST("h_all_fulljet_NEF"), nef, 1.0);
      registry.fill(HIST("h2_all_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), 1.0);
      registry.fill(HIST("h2_all_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, 1.0);
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef, 1.0);

      // CASE 2: Additional leading jet processing
      if (isLeading) {
        registry.fill(HIST("h_leading_fulljet_pt"), jetPt, 1.0);
        registry.fill(HIST("h_leading_fulljet_Nch"), numberOfChargedParticles, 1.0);
        registry.fill(HIST("h_leading_fulljet_NEF"), nef, 1.0);
        registry.fill(HIST("h2_leading_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), 1.0);
        registry.fill(HIST("h2_leading_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, 1.0);
        registry.fill(HIST("h3_leading_fulljet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef, 1.0);
      }

      // CASE 3: Additional first sub-leading jet processing
      if (isSubLeading) {
        registry.fill(HIST("h_subleading_fulljet_pt"), jetPt, 1.0);
        registry.fill(HIST("h_subleading_fulljet_Nch"), numberOfChargedParticles, 1.0);
        registry.fill(HIST("h_subleading_fulljet_NEF"), nef, 1.0);
        registry.fill(HIST("h2_subleading_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), 1.0);
        registry.fill(HIST("h2_subleading_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, 1.0);
        registry.fill(HIST("h3_subleading_fulljet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef, 1.0);
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processMBMCDCollisionsWithMultiplicity, "MB MCD Collisions for Full Jets Multiplicity Studies", false);

  void processMCDCollisionsWeightedWithMultiplicity(soa::Filtered<EMCCollisionsMCD>::iterator const& collision, JetTableMCDJoined const& mcdjets, aod::JMcCollisions const&, aod::JetTracks const& tracks, ClusterWithCorrections const& clusters)
  {
    bool eventAccepted = false;
    float eventWeight = collision.mcCollision().weight();

    registry.fill(HIST("hEventmultiplicityCounter"), 0.5, eventWeight); // allWeightedDetColl
    if (std::fabs(collision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 1.5, eventWeight); // WeightedDetCollWithVertexZ

    if (doMBGapTrigger && collision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      registry.fill(HIST("hEventmultiplicityCounter"), 2.5, eventWeight); // MBRejectedDetEvents
      return;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
      registry.fill(HIST("hEventmultiplicityCounter"), 3.5, eventWeight); // WeightedEventsNotSatisfyingEventSelection
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
          registry.fill(HIST("hEventmultiplicityCounter"), 4.5, eventWeight); // EMCreadoutWeightedDetJJEventsWithkTVXinEMC
        }
      }
    } else {
      if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
        eventAccepted = true;
        registry.fill(HIST("hEventmultiplicityCounter"), 4.5, eventWeight); // EMCreadoutWeightedDetJJEventsWithkTVXinEMC
      }
    }

    if (!eventAccepted) {
      registry.fill(HIST("hEventmultiplicityCounter"), 5.5, eventWeight); // AllRejectedWeightedDetEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 6.5, eventWeight); // EMCAcceptedWeightedDetColl
    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
    }
    registry.fill(HIST("hEventmultiplicityCounter"), 7.5, eventWeight); // EMCAcceptedWeightedCollAfterTrackSel

    registry.fill(HIST("h_FT0Mults_occupancy"), collision.multFT0M(), eventWeight);

    std::vector<typename std::decay_t<decltype(mcdjets)>::iterator> selectedJets;
    int nJetsThisEvent = 0;

    for (auto const& mcdjet : mcdjets) {
      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      // Declare variables to store filtered track/cluster pT
      double filteredTrackPt = 0.0;
      double filteredClusterPt = 0.0;

      if (mcdjet.collisionId() != collision.globalIndex()) {
        LOGF(warn, "Jet with pT %.2f belongs to collision %d but processing collision %d", mcdjet.pt(), mcdjet.collisionId(), collision.globalIndex());
        continue;
      }

      if (mcdjet.pt() > pTHatMaxMCD * pTHat || pTHat < pTHatAbsoluteMin) { // MCD jets outlier rejection
        return;
      }
      if (!jetfindingutilities::isInEtaAcceptance(mcdjet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (mcdjet.phi() < jetPhiMin || mcdjet.phi() > jetPhiMax) {
        continue;
      }
      if (!isAcceptedRecoJet<aod::JetTracks, ClusterWithCorrections>(mcdjet, filteredTrackPt, filteredClusterPt)) {
        continue;
      }
      selectedJets.push_back(mcdjet);
      nJetsThisEvent++;
    }
    // 1. Sort selected jets by pT before processing
    std::sort(selectedJets.begin(), selectedJets.end(),
              [](auto const& a, auto const& b) { return (*a).pt() > (*b).pt(); });

    if (selectedJets.size() == 0) { // no jets = no leading jet
      return;
    } else {
      // Jet multiplicity per event
      registry.fill(HIST("h_all_fulljet_Njets"), selectedJets.size(), eventWeight);

      // Select Leading Jet for N_ch calculation (for every leading jet that is found). There's always one leading jet per event!
      auto const& leadingJet = *selectedJets[0];
      auto const& leadingJetPt = leadingJet.pt(); // jet pT distribution of the leading jet
      registry.fill(HIST("h_Leading_full_jet_pt"), leadingJetPt, eventWeight);
      registry.fill(HIST("h2_full_jet_leadingJetPt_vs_counts"), leadingJetPt, nJetsThisEvent, eventWeight);
    }

    if (selectedJets.size() > 1) {
      auto const& subLeadingJet = *selectedJets[1];
      auto const& subLeadingJetPt = subLeadingJet.pt(); // jet pT distribution of the subleading jet i.e. 2nd leading jet
      registry.fill(HIST("h_SubLeading_full_jet_pt"), subLeadingJetPt, eventWeight);
      registry.fill(HIST("h2_full_jet_subLeadingJetPt_vs_counts"), subLeadingJetPt, nJetsThisEvent, eventWeight);
    }
    // Process ALL selected jets (not just leading)
    for (size_t i = 0; i < selectedJets.size(); i++) {
      auto const& jet = *selectedJets[i];
      float jetPt = jet.pt();
      bool isLeading = (i == 0);
      bool isSubLeading = (i == 1 && selectedJets.size() > 1); // first sub-leading jet

      // Count charged particles(NCh) for this jet
      int numberOfChargedParticles = 0;
      for (const auto& jettrack : jet.tracks_as<aod::JetTracks>()) {
        if (jetderiveddatautilities::selectTrack(jettrack, trackSelection)) {
          numberOfChargedParticles++;
        } else {
          continue;
        }
      }

      // Calculate neutral energy fraction for this jet
      float neutralEnergy = 0.0;
      for (const auto& jetcluster : jet.clusters_as<ClusterWithCorrections>()) {
        neutralEnergy += jetcluster.energy();
      }
      float nef = neutralEnergy / jet.energy();

      // CASE 1: Fill histograms for ALL selected jets
      registry.fill(HIST("h_all_fulljet_pt"), jetPt, eventWeight);
      registry.fill(HIST("h_all_fulljet_Nch"), numberOfChargedParticles, eventWeight);
      registry.fill(HIST("h_all_fulljet_NEF"), nef, eventWeight);
      registry.fill(HIST("h2_all_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), eventWeight);
      registry.fill(HIST("h2_all_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, eventWeight);
      registry.fill(HIST("h3_full_jet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef, eventWeight);

      // CASE 2: Additional leading jet processing
      if (isLeading) {
        registry.fill(HIST("h_leading_fulljet_pt"), jetPt, eventWeight);
        registry.fill(HIST("h_leading_fulljet_Nch"), numberOfChargedParticles, eventWeight);
        registry.fill(HIST("h_leading_fulljet_NEF"), nef, eventWeight);
        registry.fill(HIST("h2_leading_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), eventWeight);
        registry.fill(HIST("h2_leading_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, eventWeight);
        registry.fill(HIST("h3_leading_fulljet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef, eventWeight);
      }

      // CASE 3: Additional first sub-leading jet processing
      if (isSubLeading) {
        registry.fill(HIST("h_subleading_fulljet_pt"), jetPt, eventWeight);
        registry.fill(HIST("h_subleading_fulljet_Nch"), numberOfChargedParticles, eventWeight);
        registry.fill(HIST("h_subleading_fulljet_NEF"), nef, eventWeight);
        registry.fill(HIST("h2_subleading_fulljet_jetpTDet_vs_FT0Mults"), jetPt, collision.multFT0M(), eventWeight);
        registry.fill(HIST("h2_subleading_fulljet_jetpTDet_vs_Nch"), jetPt, numberOfChargedParticles, eventWeight);
        registry.fill(HIST("h3_subleading_fulljet_jetpTDet_FT0Mults_nef"), jetPt, collision.multFT0M(), nef, eventWeight);
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processMCDCollisionsWeightedWithMultiplicity, "Weighted MCD Collisions for Full Jets Multiplicity Studies", false);

  void processMBMCPCollisionsWithMultiplicity(aod::JetMcCollision const& mccollision,
                                              JetTableMCPJoined const& jets, aod::JetParticles const& /*particles*/,
                                              soa::SmallGroups<EMCCollisionsMCD> const& collisions)
  {
    bool eventAccepted = false;
    double weight = 1.0;
    float pTHat = 10. / (std::pow(weight, 1.0 / pTHatExponent));
    float mcpMult = 0.0f;
    float mcpCent = 0.0f;

    registry.fill(HIST("hPartEventmultiplicityCounter"), 0.5); // allMcColl
    if (std::fabs(mccollision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hPartEventmultiplicityCounter"), 1.5); // McCollWithVertexZ

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hPartEventmultiplicityCounter"), 2.5); // RejectedPartCollWithOutliers
        return;
      }
    }

    if (doMBGapTrigger && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      // Fill rejected MB events;
      registry.fill(HIST("hPartEventmultiplicityCounter"), 3.5); // MBRejectedPartEvents
      return;
    }

    // Perform MC Collision matching, i.e. match the current MC collision to its associated reco (MCD) collision
    //  to get the corresponding FT0M component at the particle level
    auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, mccollision.globalIndex());
    registry.fill(HIST("hRecoMatchesPerMcCollisionMult"), collisionspermcpjet.size()); // for split vertices QA

    if (collisionspermcpjet.size() == 0 || collisionspermcpjet.size() < 1) {
      registry.fill(HIST("hPartEventmultiplicityCounter"), 4.5); // RejectedPartCollForDetCollWithSize0or<1
      return;
    }
    registry.fill(HIST("hPartEventmultiplicityCounter"), 5.5); // AcceptedPartCollWithSize>1

    for (auto const& collision : collisionspermcpjet) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        continue;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          if (collision.alias_bit(kTVXinEMC)) {
            eventAccepted = true;
            registry.fill(HIST("hPartEventmultiplicityCounter"), 6.5); // EMCreadoutDetEventsWithkTVXinEMC
          }
        }
      } else {
        if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
          registry.fill(HIST("hPartEventmultiplicityCounter"), 6.5); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
      mcpMult += collision.multFT0M();
      mcpCent += collision.centFT0M();
    } // collision loop ends

    if (!eventAccepted) {
      registry.fill(HIST("hPartEventmultiplicityCounter"), 7.5); // AllRejectedPartEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hPartEventmultiplicityCounter"), 8.5); // EMCAcceptedPartColl
    registry.fill(HIST("hMCCollMatchedFT0Mult"), mcpMult);
    registry.fill(HIST("hMCCollMatchedFT0Cent"), mcpCent);

    std::vector<typename std::decay_t<decltype(jets)>::iterator> selectedJets;
    int nJetsThisEvent = 0;

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        continue;
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax)
        continue;
      // if (!isAcceptedPartJet<aod::JetParticles>(jet))
      //   continue;

      selectedJets.push_back(jet);
      nJetsThisEvent++;
    }
    // 1. Sort selected jets by pT before processing
    std::sort(selectedJets.begin(), selectedJets.end(),
              [](auto const& a, auto const& b) { return (*a).pt() > (*b).pt(); });

    if (selectedJets.size() == 0) { // no jets = no leading jet
      return;
    } else {
      // Jet multiplicity per event
      registry.fill(HIST("h_all_fulljet_Njets_part"), selectedJets.size(), 1.0);

      // Select Leading Jet for N_ch calculation (for every leading jet that is found). There's always one leading jet per event!
      auto const& leadingJet = *selectedJets[0];
      auto const& leadingJetPt = leadingJet.pt(); // jet pT distribution of the leading jet
      registry.fill(HIST("h_Leading_full_jet_pt_part"), leadingJetPt, 1.0);
      registry.fill(HIST("h2_full_jet_leadingJetPt_vs_counts_part"), leadingJetPt, nJetsThisEvent, 1.0);
    }

    if (selectedJets.size() > 1) {
      auto const& subLeadingJet = *selectedJets[1];
      auto const& subLeadingJetPt = subLeadingJet.pt(); // jet pT distribution of the subleading jet i.e. 2nd leading jet
      registry.fill(HIST("h_SubLeading_full_jet_pt_part"), subLeadingJetPt, 1.0);
      registry.fill(HIST("h2_full_jet_subLeadingJetPt_vs_counts_part"), subLeadingJetPt, nJetsThisEvent, 1.0);
    }
    // Process ALL selected jets (not just leading)
    for (size_t i = 0; i < selectedJets.size(); i++) {
      auto const& jet = *selectedJets[i];
      float jetPt = jet.pt();
      bool isLeading = (i == 0);
      bool isSubLeading = (i == 1 && selectedJets.size() > 1); // first sub-leading jet

      // Count charged particles(NCh) and neutral particles (NNe) for this jet
      int numberOfChargedParticles = 0;
      int numberOfNeutralParticles = 0;
      float neutralEnergy = 0.0f;
      for (const auto& constituent : jet.tracks_as<aod::JetParticles>()) {
        auto pdgParticle = pdgDatabase->GetParticle(constituent.pdgCode());
        if (pdgParticle->Charge() == 0) {
          numberOfNeutralParticles++;
          neutralEnergy += constituent.e();
        } else {
          numberOfChargedParticles++;
        }
      }
      float nef = neutralEnergy / jet.energy();

      // CASE 1: Fill histograms for ALL selected jets
      registry.fill(HIST("h_all_fulljet_pt_part"), jetPt, 1.0);
      registry.fill(HIST("h_all_fulljet_Nch_part"), numberOfChargedParticles, 1.0);
      registry.fill(HIST("h_all_fulljet_Nne_part"), numberOfNeutralParticles, 1.0);
      registry.fill(HIST("h_all_fulljet_NEF_part"), nef, 1.0);
      registry.fill(HIST("h2_all_fulljet_jetpT_vs_FT0Mults_part"), jetPt, mcpMult, 1.0);
      registry.fill(HIST("h2_all_fulljet_jetpT_vs_Nch_part"), jetPt, numberOfChargedParticles, 1.0);
      registry.fill(HIST("h3_full_jet_jetpT_FT0Mults_nef_part"), jetPt, mcpMult, nef, 1.0);

      // CASE 2: Additional leading jet processing
      if (isLeading) {
        registry.fill(HIST("h_leading_fulljet_pt_part"), jetPt, 1.0);
        registry.fill(HIST("h_leading_fulljet_Nch_part"), numberOfChargedParticles, 1.0);
        registry.fill(HIST("h_leading_fulljet_Nne_part"), numberOfNeutralParticles, 1.0);
        registry.fill(HIST("h_leading_fulljet_NEF_part"), nef, 1.0);
        registry.fill(HIST("h2_leading_fulljet_jetpT_vs_FT0Mults_part"), jetPt, mcpMult, 1.0);
        registry.fill(HIST("h2_leading_fulljet_jetpT_vs_Nch_part"), jetPt, numberOfChargedParticles, 1.0);
        registry.fill(HIST("h3_leading_fulljet_jetpT_FT0Mults_nef_part"), jetPt, mcpMult, nef, 1.0);
      }

      // CASE 3: Additional first sub-leading jet processing
      if (isSubLeading) {
        registry.fill(HIST("h_subleading_fulljet_pt_part"), jetPt, 1.0);
        registry.fill(HIST("h_subleading_fulljet_Nch_part"), numberOfChargedParticles, 1.0);
        registry.fill(HIST("h_subleading_fulljet_Nne_part"), numberOfNeutralParticles, 1.0);
        registry.fill(HIST("h_subleading_fulljet_NEF_part"), nef, 1.0);
        registry.fill(HIST("h2_subleading_fulljet_jetpT_vs_FT0Mults_part"), jetPt, mcpMult, 1.0);
        registry.fill(HIST("h2_subleading_fulljet_jetpT_vs_Nch_part"), jetPt, numberOfChargedParticles, 1.0);
        registry.fill(HIST("h3_subleading_fulljet_jetpT_FT0Mults_nef_part"), jetPt, mcpMult, nef, 1.0);
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processMBMCPCollisionsWithMultiplicity, "MB MCP Collisions for Full Jets Multiplicity Studies", false);

  void processMBMCPCollisionsWeightedWithMultiplicity(aod::JetMcCollision const& mccollision,
                                                      JetTableMCPJoined const& jets, aod::JetParticles const& /*particles*/,
                                                      soa::SmallGroups<EMCCollisionsMCD> const& collisions)
  {
    bool eventAccepted = false;
    float pTHat = 10. / (std::pow(mccollision.weight(), 1.0 / pTHatExponent));
    float mcpMult = 0.0f;
    float mcpCent = 0.0f;

    registry.fill(HIST("hPartEventmultiplicityCounter"), 0.5, mccollision.weight()); // allWeightedMcColl
    if (std::fabs(mccollision.posZ()) > vertexZCut) {
      return;
    }
    registry.fill(HIST("hPartEventmultiplicityCounter"), 1.5, mccollision.weight()); // WeightedMcCollWithVertexZ

    // outlier check: for every outlier jet, reject the whole event
    for (auto const& jet : jets) {
      if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
        registry.fill(HIST("hPartEventmultiplicityCounter"), 2.5, mccollision.weight()); // RejectedWeightedPartCollWithOutliers
        return;
      }
    }

    if (doMBGapTrigger && mccollision.getSubGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      // Fill rejected MB events;
      registry.fill(HIST("hPartEventmultiplicityCounter"), 3.5, mccollision.weight()); // MBRejectedPartEvents
      return;
    }

    // Perform MC Collision matching, i.e. match the current MC collision to its associated reco (MCD) collision
    //  to get the corresponding FT0M component at the particle level
    auto collisionspermcpjet = collisions.sliceBy(CollisionsPerMCPCollision, mccollision.globalIndex());
    registry.fill(HIST("hRecoMatchesPerMcCollisionMult"), collisionspermcpjet.size(), mccollision.weight()); // for split vertices QA

    if (collisionspermcpjet.size() == 0 || collisionspermcpjet.size() < 1) {
      registry.fill(HIST("hPartEventmultiplicityCounter"), 4.5, mccollision.weight()); // RejectedWeightedPartCollForDetCollWithSize0or<1
      return;
    }
    registry.fill(HIST("hPartEventmultiplicityCounter"), 5.5, mccollision.weight()); // AcceptedWeightedPartCollWithSize>1

    for (auto const& collision : collisionspermcpjet) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, doMBGapTrigger)) {
        continue;
      }
      if (doEMCALEventWorkaround) {
        if (collision.isEmcalReadout() && !collision.isAmbiguous()) { // i.e. EMCAL has a cell content
          if (collision.alias_bit(kTVXinEMC)) {
            eventAccepted = true;
            registry.fill(HIST("hPartEventmultiplicityCounter"), 6.5, mccollision.weight()); // EMCreadoutDetEventsWithkTVXinEMC
          }
        }
      } else {
        if (!collision.isAmbiguous() && jetderiveddatautilities::eventEMCAL(collision) && collision.alias_bit(kTVXinEMC)) {
          eventAccepted = true;
          registry.fill(HIST("hPartEventmultiplicityCounter"), 6.5, mccollision.weight()); // EMCreadoutDetEventsWithkTVXinEMC
        }
      }
      mcpMult += collision.multFT0M();
      mcpCent += collision.centFT0M();
    } // collision loop ends

    if (!eventAccepted) {
      registry.fill(HIST("hPartEventmultiplicityCounter"), 7.5, mccollision.weight()); // AllRejectedWeightedPartEventsAfterEMCEventSelection
      return;
    }
    registry.fill(HIST("hPartEventmultiplicityCounter"), 8.5, mccollision.weight()); // EMCAcceptedWeightedPartColl
    registry.fill(HIST("hMCCollMatchedFT0Mult"), mcpMult, mccollision.weight());
    registry.fill(HIST("hMCCollMatchedFT0Cent"), mcpCent, mccollision.weight());

    std::vector<typename std::decay_t<decltype(jets)>::iterator> selectedJets;
    int nJetsThisEvent = 0;

    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax))
        continue;
      if (jet.phi() < jetPhiMin || jet.phi() > jetPhiMax)
        continue;
      // if (!isAcceptedPartJet<aod::JetParticles>(jet))
      //   continue;

      selectedJets.push_back(jet);
      nJetsThisEvent++;
    }
    // 1. Sort selected jets by pT before processing
    std::sort(selectedJets.begin(), selectedJets.end(),
              [](auto const& a, auto const& b) { return (*a).pt() > (*b).pt(); });

    if (selectedJets.size() == 0) { // no jets = no leading jet
      return;
    } else {
      // Jet multiplicity per event
      registry.fill(HIST("h_all_fulljet_Njets_part"), selectedJets.size(), mccollision.weight());

      // Select Leading Jet for N_ch calculation (for every leading jet that is found). There's always one leading jet per event!
      auto const& leadingJet = *selectedJets[0];
      auto const& leadingJetPt = leadingJet.pt(); // jet pT distribution of the leading jet
      registry.fill(HIST("h_Leading_full_jet_pt_part"), leadingJetPt, mccollision.weight());
      registry.fill(HIST("h2_full_jet_leadingJetPt_vs_counts_part"), leadingJetPt, nJetsThisEvent, mccollision.weight());
    }

    if (selectedJets.size() > 1) {
      auto const& subLeadingJet = *selectedJets[1];
      auto const& subLeadingJetPt = subLeadingJet.pt(); // jet pT distribution of the subleading jet i.e. 2nd leading jet
      registry.fill(HIST("h_SubLeading_full_jet_pt_part"), subLeadingJetPt, mccollision.weight());
      registry.fill(HIST("h2_full_jet_subLeadingJetPt_vs_counts_part"), subLeadingJetPt, nJetsThisEvent, mccollision.weight());
    }
    // Process ALL selected jets (not just leading)
    for (size_t i = 0; i < selectedJets.size(); i++) {
      auto const& jet = *selectedJets[i];
      float jetPt = jet.pt();
      bool isLeading = (i == 0);
      bool isSubLeading = (i == 1 && selectedJets.size() > 1); // first sub-leading jet

      // Count charged particles(NCh) and neutral particles (NNe) for this jet
      int numberOfChargedParticles = 0;
      int numberOfNeutralParticles = 0;
      float neutralEnergy = 0.0f;
      for (const auto& constituent : jet.tracks_as<aod::JetParticles>()) {
        auto pdgParticle = pdgDatabase->GetParticle(constituent.pdgCode());
        if (pdgParticle->Charge() == 0) {
          numberOfNeutralParticles++;
          neutralEnergy += constituent.e();
        } else {
          numberOfChargedParticles++;
        }
      }
      float nef = neutralEnergy / jet.energy();

      // CASE 1: Fill histograms for ALL selected jets
      registry.fill(HIST("h_all_fulljet_pt_part"), jetPt, mccollision.weight());
      registry.fill(HIST("h_all_fulljet_Nch_part"), numberOfChargedParticles, mccollision.weight());
      registry.fill(HIST("h_all_fulljet_Nne_part"), numberOfNeutralParticles, mccollision.weight());
      registry.fill(HIST("h_all_fulljet_NEF_part"), nef, mccollision.weight());
      registry.fill(HIST("h2_all_fulljet_jetpT_vs_FT0Mults_part"), jetPt, mcpMult, mccollision.weight());
      registry.fill(HIST("h2_all_fulljet_jetpT_vs_Nch_part"), jetPt, numberOfChargedParticles, mccollision.weight());
      registry.fill(HIST("h3_full_jet_jetpT_FT0Mults_nef_part"), jetPt, mcpMult, nef, mccollision.weight());

      // CASE 2: Additional leading jet processing
      if (isLeading) {
        registry.fill(HIST("h_leading_fulljet_pt_part"), jetPt, mccollision.weight());
        registry.fill(HIST("h_leading_fulljet_Nch_part"), numberOfChargedParticles, mccollision.weight());
        registry.fill(HIST("h_leading_fulljet_Nne_part"), numberOfNeutralParticles, mccollision.weight());
        registry.fill(HIST("h_leading_fulljet_NEF_part"), nef, mccollision.weight());
        registry.fill(HIST("h2_leading_fulljet_jetpT_vs_FT0Mults_part"), jetPt, mcpMult, mccollision.weight());
        registry.fill(HIST("h2_leading_fulljet_jetpT_vs_Nch_part"), jetPt, numberOfChargedParticles, mccollision.weight());
        registry.fill(HIST("h3_leading_fulljet_jetpT_FT0Mults_nef_part"), jetPt, mcpMult, nef, mccollision.weight());
      }

      // CASE 3: Additional first sub-leading jet processing
      if (isSubLeading) {
        registry.fill(HIST("h_subleading_fulljet_pt_part"), jetPt, mccollision.weight());
        registry.fill(HIST("h_subleading_fulljet_Nch_part"), numberOfChargedParticles, mccollision.weight());
        registry.fill(HIST("h_subleading_fulljet_Nne_part"), numberOfNeutralParticles, mccollision.weight());
        registry.fill(HIST("h_subleading_fulljet_NEF_part"), nef, mccollision.weight());
        registry.fill(HIST("h2_subleading_fulljet_jetpT_vs_FT0Mults_part"), jetPt, mcpMult, mccollision.weight());
        registry.fill(HIST("h2_subleading_fulljet_jetpT_vs_Nch_part"), jetPt, numberOfChargedParticles, mccollision.weight());
        registry.fill(HIST("h3_subleading_fulljet_jetpT_FT0Mults_nef_part"), jetPt, mcpMult, nef, mccollision.weight());
      }
    }
  }
  PROCESS_SWITCH(FullJetSpectra, processMBMCPCollisionsWeightedWithMultiplicity, "MB MCP Weighted Collisions for Full Jets Multiplicity Studies", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FullJetSpectra>(cfgc)};
}
