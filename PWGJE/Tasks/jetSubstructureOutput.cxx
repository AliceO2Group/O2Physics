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

// jet substructure tree filling task (subscribing to jet finder hf and jet substructure tasks)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include <vector>
#include <algorithm>
#include <utility>
#include <map>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

struct JetSubstructureOutputTask {

  Produces<aod::CJetCOs> collisionOutputTableData;
  Produces<aod::CJetOs> jetOutputTableData;
  Produces<aod::CJetSSOs> jetSubstructureOutputTableData;
  Produces<aod::CJetMOs> jetMatchingOutputTableData;
  Produces<aod::CEWSJetCOs> collisionOutputTableDataSub;
  Produces<aod::CEWSJetOs> jetOutputTableDataSub;
  Produces<aod::CEWSJetSSOs> jetSubstructureOutputTableDataSub;
  Produces<aod::CEWSJetMOs> jetMatchingOutputTableDataSub;
  Produces<aod::CMCDJetCOs> collisionOutputTableMCD;
  Produces<aod::CMCDJetOs> jetOutputTableMCD;
  Produces<aod::CMCDJetSSOs> jetSubstructureOutputTableMCD;
  Produces<aod::CMCDJetMOs> jetMatchingOutputTableMCD;
  Produces<aod::CMCPJetCOs> collisionOutputTableMCP;
  Produces<aod::CMCPJetOs> jetOutputTableMCP;
  Produces<aod::CMCPJetSSOs> jetSubstructureOutputTableMCP;
  Produces<aod::CMCPJetMOs> jetMatchingOutputTableMCP;

  Configurable<float> jetPtMinData{"jetPtMinData", 0.0, "minimum jet pT cut for data jets"};
  Configurable<float> jetPtMinDataSub{"jetPtMinDataSub", 0.0, "minimum jet pT cut for eventwise constituent subtracted data jets"};
  Configurable<float> jetPtMinMCD{"jetPtMinMCD", 0.0, "minimum jet pT cut for mcd jets"};
  Configurable<float> jetPtMinMCP{"jetPtMinMCP", 0.0, "minimum jet pT cut for mcp jets"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track pseudorapidity"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track pseudorapidity"};

  std::map<int32_t, int32_t> jetMappingData;
  std::map<int32_t, int32_t> jetMappingDataSub;
  std::map<int32_t, int32_t> jetMappingMCD;
  std::map<int32_t, int32_t> jetMappingMCP;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  template <typename T, typename U, typename V>
  void fillJetTables(T const& jet, int32_t collisionIndex, U& jetOutputTable, V& jetSubstructureOutputTable, std::map<int32_t, int32_t>& jetMapping)
  {
    std::vector<float> energyMotherVec;
    std::vector<float> ptLeadingVec;
    std::vector<float> ptSubLeadingVec;
    std::vector<float> thetaVec;
    std::vector<float> pairPtVec;
    std::vector<float> pairEnergyVec;
    std::vector<float> pairThetaVec;
    auto energyMotherSpan = jet.energyMother();
    auto ptLeadingSpan = jet.ptLeading();
    auto ptSubLeadingSpan = jet.ptSubLeading();
    auto thetaSpan = jet.theta();
    auto pairPtSpan = jet.pairPt();
    auto pairEnergySpan = jet.pairEnergy();
    auto pairThetaSpan = jet.pairTheta();
    std::copy(energyMotherSpan.begin(), energyMotherSpan.end(), std::back_inserter(energyMotherVec));
    std::copy(ptLeadingSpan.begin(), ptLeadingSpan.end(), std::back_inserter(ptLeadingVec));
    std::copy(ptSubLeadingSpan.begin(), ptSubLeadingSpan.end(), std::back_inserter(ptSubLeadingVec));
    std::copy(thetaSpan.begin(), thetaSpan.end(), std::back_inserter(thetaVec));
    std::copy(pairPtSpan.begin(), pairPtSpan.end(), std::back_inserter(pairPtVec));
    std::copy(pairEnergySpan.begin(), pairEnergySpan.end(), std::back_inserter(pairEnergyVec));
    std::copy(pairThetaSpan.begin(), pairThetaSpan.end(), std::back_inserter(pairThetaVec));
    jetOutputTable(collisionIndex, collisionIndex, jet.pt(), jet.phi(), jet.eta(), jet.y(), jet.r(), jet.tracksIds().size()); // second collision index is a dummy coloumn mirroring the hf candidate
    jetSubstructureOutputTable(jetOutputTable.lastIndex(), energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, jet.nSub2DR(), jet.nSub1(), jet.nSub2(), pairPtVec, pairEnergyVec, pairThetaVec, jet.angularity(), jet.ptLeadingConstituent());
    jetMapping.insert(std::make_pair(jet.globalIndex(), jetOutputTable.lastIndex()));
  }

  template <bool isMCP, typename T, typename U, typename V, typename M, typename N>
  void analyseCharged(T const& collision, U const& jets, V& collisionOutputTable, M& jetOutputTable, N& jetSubstructureOutputTable, std::map<int32_t, int32_t>& jetMapping, float jetPtMin, float eventWeight)
  {
    int nJetInCollision = 0;
    int32_t collisionIndex = -1;
    for (const auto& jet : jets) {
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          if (nJetInCollision == 0) {
            float centrality = -1.0;
            uint8_t eventSel = 0.0;
            if constexpr (!isMCP) {
              centrality = collision.centrality();
              eventSel = collision.eventSel();
            }
            collisionOutputTable(collision.posZ(), centrality, eventSel, eventWeight);
            collisionIndex = collisionOutputTable.lastIndex();
          }
          nJetInCollision++;
          fillJetTables(jet, collisionIndex, jetOutputTable, jetSubstructureOutputTable, jetMapping);
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void analyseMatched(T const& jets, U const& /*jetsTag*/, std::map<int32_t, int32_t>& jetMapping, std::map<int32_t, int32_t>& jetTagMapping, V& matchingOutputTable, float jetPtMin)
  {
    std::vector<int> candMatching;
    for (const auto& jet : jets) {
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          if (jet.has_matchedJetGeo()) {
            for (auto& jetTagId : jet.matchedJetGeoIds()) {
              auto jetTagIndex = jetTagMapping.find(jetTagId);
              if (jetTagIndex != jetTagMapping.end()) {
                geoMatching.push_back(jetTagIndex->second);
              }
            }
          }
          if (jet.has_matchedJetPt()) {
            for (auto& jetTagId : jet.matchedJetPtIds()) {
              auto jetTagIndex = jetTagMapping.find(jetTagId);
              if (jetTagIndex != jetTagMapping.end()) {
                ptMatching.push_back(jetTagIndex->second);
              }
            }
          }
          int storedJetIndex = -1;
          auto jetIndex = jetMapping.find(jet.globalIndex());
          if (jetIndex != jetMapping.end()) {
            storedJetIndex = jetIndex->second;
          }
          matchingOutputTable(storedJetIndex, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }

  void processClearMaps(aod::JetCollisions const&)
  {
    jetMappingData.clear();
    jetMappingDataSub.clear();
    jetMappingMCD.clear();
    jetMappingMCP.clear();
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processClearMaps, "process function that clears all the maps in each dataframe", true);

  void processOutputData(aod::JetCollision const& collision,
                         soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::CJetSSs> const& jets)
  {
    analyseCharged<false>(collision, jets, collisionOutputTableData, jetOutputTableData, jetSubstructureOutputTableData, jetMappingData, jetPtMinData, 1.0);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputData, "jet substructure output Data", false);

  void processOutputDataSub(aod::JetCollision const& collision,
                            soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::CEWSJetSSs> const& jets)
  {
    analyseCharged<false>(collision, jets, collisionOutputTableDataSub, jetOutputTableDataSub, jetSubstructureOutputTableDataSub, jetMappingDataSub, jetPtMinDataSub, 1.0);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputDataSub, "jet substructure output event-wise subtracted Data", false);

  void processOutputMatchingData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets> const& jets,
                                 soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets> const& jetsSub)
  {
    analyseMatched(jets, jetsSub, jetMappingData, jetMappingDataSub, jetMatchingOutputTableData, jetPtMinData);
    analyseMatched(jetsSub, jets, jetMappingDataSub, jetMappingData, jetMatchingOutputTableDataSub, jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMatchingData, "jet matching output Data", false);

  void processOutputMCD(aod::JetCollisionMCD const& collision, aod::JetMcCollisions const&,
                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::CMCDJetSSs> const& jets)
  {
    analyseCharged<false>(collision, jets, collisionOutputTableMCD, jetOutputTableMCD, jetSubstructureOutputTableMCD, jetMappingMCD, jetPtMinMCD, collision.mcCollision().weight());
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCD, "jet substructure output MCD", false);

  void processOutputMCP(aod::JetMcCollision const& collision,
                        soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::CMCPJetSSs> const& jets)
  {
    analyseCharged<true>(collision, jets, collisionOutputTableMCP, jetOutputTableMCP, jetSubstructureOutputTableMCP, jetMappingMCP, jetPtMinMCP, collision.weight());
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCP, "jet substructure output MCP", false);

  void processOutputMatchingMC(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& jetsMCD,
                               soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& jetsMCP)
  {
    analyseMatched(jetsMCD, jetsMCP, jetMappingMCD, jetMappingMCP, jetMatchingOutputTableMCD, jetPtMinMCD);
    analyseMatched(jetsMCP, jetsMCD, jetMappingMCP, jetMappingMCD, jetMatchingOutputTableMCP, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMatchingMC, "jet matching output MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureOutputTask>(
    cfgc, TaskName{"jet-substructure-output"})};
}
