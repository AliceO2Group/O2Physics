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

// heavy-flavour jet substructure tree filling task (subscribing to jet finder hf and jet substructure hf tasks)
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename CandidateCollisionTable, typename CandidateTable, typename CandidateTableMCD, typename CandidateTableMCP, typename TracksSub, typename JetTableData, typename JetMatchedTableData, typename OutputCollisionTableData, typename OutputTableData, typename SubstructureOutputTableData, typename MatchingOutputTableData, typename JetTableMCD, typename OutputCollisionTableMCD, typename OutputTableMCD, typename SubstructureOutputTableMCD, typename MatchingOutputTableMCD, typename JetTableMCP, typename OutputCollisionTableMCP, typename OutputTableMCP, typename SubstructureOutputTableMCP, typename MatchingOutputTableMCP, typename JetTableDataSub, typename OutputCollisionTableDataSub, typename OutputTableDataSub, typename SubstructureOutputTableDataSub, typename MatchingOutputTableDataSub, typename CandidateCollisionOutputTable, typename CandidateOutputTable, typename CandidateParOutputTable, typename CandidateParExtraOutputTable, typename CandidateSelOutputTable, typename CandidateMCDOutputTable, typename CandidateMCPOutputTable>
struct JetSubstructureHFOutputTask {
  Produces<OutputCollisionTableData> collisionOutputTableData;
  Produces<OutputTableData> jetOutputTableData;
  Produces<SubstructureOutputTableData> jetSubstructureOutputTableData;
  Produces<MatchingOutputTableData> jetMatchingOutputTableData;
  Produces<OutputCollisionTableDataSub> collisionOutputTableDataSub;
  Produces<OutputTableDataSub> jetOutputTableDataSub;
  Produces<SubstructureOutputTableDataSub> jetSubstructureOutputTableDataSub;
  Produces<MatchingOutputTableDataSub> jetMatchingOutputTableDataSub;
  Produces<OutputCollisionTableMCD> collisionOutputTableMCD;
  Produces<OutputTableMCD> jetOutputTableMCD;
  Produces<SubstructureOutputTableMCD> jetSubstructureOutputTableMCD;
  Produces<MatchingOutputTableMCD> jetMatchingOutputTableMCD;
  Produces<OutputCollisionTableMCP> collisionOutputTableMCP;
  Produces<OutputTableMCP> jetOutputTableMCP;
  Produces<SubstructureOutputTableMCP> jetSubstructureOutputTableMCP;
  Produces<MatchingOutputTableMCP> jetMatchingOutputTableMCP;
  Produces<CandidateCollisionOutputTable> hfCollisionsTable;
  Produces<CandidateOutputTable> candidateTable;
  Produces<CandidateParOutputTable> candidateParsTable;
  Produces<CandidateParExtraOutputTable> candidateParExtrasTable;
  Produces<CandidateSelOutputTable> candidateSelsTable;
  Produces<CandidateMCDOutputTable> candidateMcsTable;
  Produces<CandidateMCPOutputTable> hfParticlesTable;

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track pseudorapidity"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track pseudorapidity"};

  std::map<int32_t, int32_t> candidateMapping;
  std::map<int32_t, int32_t> candidateCollisionMapping;

  std::map<int32_t, int32_t> jetMappingData;
  std::map<int32_t, int32_t> jetMappingDataSub;
  std::map<int32_t, int32_t> jetMappingMCD;
  std::map<int32_t, int32_t> jetMappingMCP;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  Filter jetSelection = aod::jet::pt >= jetPtMin;

  template <typename T, typename U, typename V, typename M>
  void fillTables(T const& jet, U const& cand, int32_t collisionIndex, int32_t candidateIndex, V& jetOutputTable, M& jetSubstructureOutputTable, std::map<int32_t, int32_t>& jetMapping)
  {
    std::vector<float> energyMotherVec;
    std::vector<float> ptLeadingVec;
    std::vector<float> ptSubLeadingVec;
    std::vector<float> thetaVec;
    auto energyMotherSpan = jet.energyMother();
    auto ptLeadingSpan = jet.ptLeading();
    auto ptSubLeadingSpan = jet.ptSubLeading();
    auto thetaSpan = jet.theta();
    std::copy(energyMotherSpan.begin(), energyMotherSpan.end(), std::back_inserter(energyMotherVec));
    std::copy(ptLeadingSpan.begin(), ptLeadingSpan.end(), std::back_inserter(ptLeadingVec));
    std::copy(ptSubLeadingSpan.begin(), ptSubLeadingSpan.end(), std::back_inserter(ptSubLeadingVec));
    std::copy(thetaSpan.begin(), thetaSpan.end(), std::back_inserter(thetaVec));
    jetOutputTable(collisionIndex, candidateIndex, jet.pt(), jet.phi(), jet.eta(), jet.r(), jet.tracksIds().size() + jet.hfcandidatesIds().size()); // here we take the decision to keep the collision index consistent with the JE framework in case it is later needed to join to other tables. The candidate Index however can be linked to the HF tables
    jetSubstructureOutputTable(jetOutputTable.lastIndex(), energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec);
    jetMapping.insert(std::make_pair(jet.globalIndex(), jetOutputTable.lastIndex()));
  }

  template <bool isMc, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& jets, V const& candidates, M& collisionOutputTable, N& jetOutputTable, O& jetSubstructureOutputTable, std::map<int32_t, int32_t>& jetMapping)
  {

    int nJetInCollision = 0;
    int32_t collisionIndex = -1;
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          auto candidate = jet.template hfcandidates_first_as<V>();
          int32_t candidateIndex = -1;
          auto candidateTableIndex = candidateMapping.find(candidate.globalIndex());
          if (candidateTableIndex != candidateMapping.end()) {
            candidateIndex = candidateTableIndex->second;
          }
          if constexpr (!isMc) {
            if (nJetInCollision == 0) {
              collisionOutputTable(collision.posZ(), collision.centrality(), collision.eventSel());
              collisionIndex = collisionOutputTable.lastIndex();
            }
            nJetInCollision++;
          }
          fillTables(jet, candidate, collisionIndex, candidateIndex, jetOutputTable, jetSubstructureOutputTable, jetMapping);
        }
      }
    }
  }

  template <bool isMc, typename T, typename U, typename V>
  void analyseCandidates(T const& jets, U const& candidateCollisions, V const& candidates)
  {

    int nJetInCollision = 0;
    int32_t candidateCollisionIndex = -1;
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {

          auto candidate = jet.template hfcandidates_first_as<V>();

          auto candidateTableIndex = candidateMapping.find(candidate.globalIndex());
          if (candidateTableIndex != candidateMapping.end()) {
            continue;
          }

          auto candidateCollision = jethfutilities::getCandidateCollision<typename V::iterator>(candidate, candidateCollisions);
          if (nJetInCollision == 0) {
            auto candidateCollisionTableIndex = candidateCollisionMapping.find(candidateCollision.globalIndex());
            if (candidateCollisionTableIndex == candidateCollisionMapping.end()) {
              jethfutilities::fillHFCollisionTable(candidateCollision, candidates, hfCollisionsTable, candidateCollisionIndex);
              candidateCollisionMapping.insert(std::make_pair(candidateCollision.globalIndex(), candidateCollisionIndex));
            } else {
              candidateCollisionIndex = candidateCollisionTableIndex->second;
            }
          }
          nJetInCollision++;
          int32_t candidateIndex = -1;
          jethfutilities::fillCandidateTable<isMc>(candidate, candidateCollisionIndex, candidateTable, candidateParsTable, candidateParExtrasTable, candidateSelsTable, candidateMcsTable, candidateIndex);
          candidateMapping.insert(std::make_pair(candidate.globalIndex(), candidateIndex));
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void analyseMatched(T const& jets, U const& jetsTag, std::map<int32_t, int32_t>& jetMapping, std::map<int32_t, int32_t>& jetTagMapping, V& matchingOutputTable)
  {
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          std::vector<int> candMatching;
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
          if (jet.has_matchedJetCand()) {
            for (auto& jetTagId : jet.matchedJetCandIds()) {
              auto jetTagIndex = jetTagMapping.find(jetTagId);
              if (jetTagIndex != jetTagMapping.end()) {
                candMatching.push_back(jetTagIndex->second);
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

  void processDummy(JetCollisions const& collisions) {}
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processOutputHFData(JetCollision const& collision,
                           JetTableData const& jets,
                           CandidateCollisionTable const& canidateCollisions,
                           CandidateTable const& candidates)
  {

    analyseCandidates<false>(jets, canidateCollisions, candidates);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputHFData, "hf output Data", false);

  void processOutputHFDataSub(JetCollision const& collision,
                              JetTableDataSub const& jets,
                              CandidateCollisionTable const& canidateCollisions,
                              CandidateTable const& candidates)
  {

    analyseCandidates<false>(jets, canidateCollisions, candidates);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputHFDataSub, "hf output Data eventwise constituent subtracted", false);

  void processOutputHFMCD(JetCollision const& collision,
                          JetTableMCD const& jets,
                          CandidateCollisionTable const& canidateCollisions,
                          CandidateTableMCD const& candidates)
  {

    analyseCandidates<true>(jets, canidateCollisions, candidates);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputHFMCD, "hf output MCD", false);

  void processOutputData(JetCollision const& collision,
                         JetTableData const& jets,
                         CandidateTable const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, collisionOutputTableData, jetOutputTableData, jetSubstructureOutputTableData, jetMappingData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputData, "hf jet substructure output Data", false);

  void processOutputDataSub(JetCollision const& collision,
                            JetTableDataSub const& jets,
                            CandidateTable const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, collisionOutputTableDataSub, jetOutputTableDataSub, jetSubstructureOutputTableDataSub, jetMappingDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputDataSub, "hf jet substructure output event-wise subtracted Data", false);

  void processOutputMatchingData(JetMatchedTableData const& jets,
                                 JetTableDataSub const& jetsSub)
  {
    analyseMatched(jets, jetsSub, jetMappingData, jetMappingDataSub, jetMatchingOutputTableData);
    analyseMatched(jetsSub, jets, jetMappingDataSub, jetMappingData, jetMatchingOutputTableDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMatchingData, "jet matching output Data", false);

  void processOutputMCD(JetCollision const& collision,
                        JetTableMCD const& jets,
                        CandidateTableMCD const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, collisionOutputTableMCD, jetOutputTableMCD, jetSubstructureOutputTableMCD, jetMappingMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMCD, "hf jet substructure output MCD", false);

  void processOutputMCP(JetMcCollision const& collision,
                        JetTableMCP const& jets,
                        CandidateTableMCP const& candidates)
  {
    analyseCharged<true>(collision, jets, candidates, collisionOutputTableMCP, jetOutputTableMCP, jetSubstructureOutputTableMCP, jetMappingMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMCP, "hf jet substructure output MCP", false);

  void processOutputMatchingMC(JetTableMCD const& jetsMCD,
                               JetTableMCP const& jetsMCP)
  {
    analyseMatched(jetsMCD, jetsMCP, jetMappingMCD, jetMappingMCP, jetMatchingOutputTableMCD);
    analyseMatched(jetsMCP, jetsMCD, jetMappingMCP, jetMappingMCD, jetMatchingOutputTableMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMatchingMC, "jet matching output MC", false);
};
using JetSubstructureOutputD0 = JetSubstructureHFOutputTask<aod::HfD0CollBases, CandidatesD0Data, CandidatesD0MCD, CandidatesD0MCP, aod::JTrackD0Subs, soa::Filtered<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0CJetSSs>>, soa::Filtered<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0CJetSSs, aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets>>, aod::D0CJetCOs, aod::D0CJetOs, aod::D0CJetSSOs, aod::D0CJetMOs, soa::Filtered<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0CMCDJetSSs, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>>, aod::D0CMCDJetCOs, aod::D0CMCDJetOs, aod::D0CMCDJetSSOs, aod::D0CMCDJetMOs, soa::Filtered<soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0CMCPJetSSs, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>>, aod::D0CMCPJetCOs, aod::D0CMCPJetOs, aod::D0CMCPJetSSOs, aod::D0CMCPJetMOs, soa::Filtered<soa::Join<aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents, aod::D0CEWSJetSSs, aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets>>, aod::D0CEWSJetCOs, aod::D0CEWSJetOs, aod::D0CEWSJetSSOs, aod::D0CEWSJetMOs, aod::StoredHfD0CollBase, aod::StoredHfD0Bases, aod::StoredHfD0Pars, aod::StoredHfD0ParEs, aod::StoredHfD0Sels, aod::StoredHfD0Mcs, aod::StoredHfD0PBases>;
// using JetSubstructureOutputLc = JetSubstructureHFOutputTask<aod::HfLcCollBases,CandidatesLcData, CandidatesLcMCD, CandidatesLcMCP, aod::JTrackLcSubs, soa::Filtered<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcCJetSSs>>, soa::Filtered<soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcCJetSSs, aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets>>, aod::LcCJetCOs, aod::LcCJetOs, aod::LcCJetSSOs, aod::LcCJetMOs, soa::Filtered<soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcCMCDJetSSs, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>>, aod::LcCMCDJetCOs, aod::LcMCDJetOs, aod::LcCMCDJetSSOs, aod::LcCMCDJetMOs, soa::Filtered<soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcCMCPJetSSs, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>>, aod::LcCMCPJetCOs, aod::LcCMCPJetOs, aod::LcCMCPJetSSOs, aod::LcCMCPJetMOs, soa::Filtered<soa::Join<aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents, aod::LcCEWSJetSSs, aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets>>, aod::LcCEWSJetCOs, aod::LcCEWSJetOs, aod::LcCEWSJetSSOs, aod::LcCEWSJetMOs, aod::StoredHfLcCollBase, aod::StoredHfLcBases, aod::StoredHfLcPars, aod::StoredHfLcParEs, aod::StoredHfLcSels, aod::StoredHfLcMcs, aod::StoredHfLcPBases>;
// using JetSubstructureOutputBplus = JetSubstructureHFOutputTask<aod::HfBplusCollBases,CandidatesBplusData, CandidatesBplusMCD, CandidatesBplusMCP, aod::JTrackBplusSubs, soa::Filtered<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs>>, soa::Filtered<soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs, aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets>>, aod::BplusCJetCOs, aod::BplusCJetOs, aod::BplusCJetSSOs, aod::BplusCJetMOs, soa::Filtered<soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusCMCDJetSSs, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>>, aod::BplusCMCDJetCOs, aod::BplusCMCDJetOs, aod::BplusCMCDJetSSOs,  aod::BplusCMCDJetMOs, soa::Filtered<soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusCMCPJetSSs, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>>, aod::BplusCMCPJetCOs, aod::BplusCMCPJetOs, aod::BplusCMCPJetSSOs, aod::BplusCMCPJetMOs, soa::Filtered<soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents, aod::BplusCEWSJetSSs, aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets>>, aod::BplusCEWSJetCOs, aod::BplusCEWSJetOs, aod::BplusCEWSJetSSOs, aod::BplusCEWSJetMOs, aod::StoredHfBplusCollBase, aod::StoredHfBplusBases, aod::StoredHfBplusPars, aod::StoredHfBplusParEs, aod::StoredHfBplusSels, aod::StoredHfBplusMcs, aod::StoredHfBplusPBases>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputD0>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructured0-output"}));
  // tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputLc>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructurelc-output"}));
  // tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputBplus>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructurebplus-output"}));
  return WorkflowSpec{tasks};
}
