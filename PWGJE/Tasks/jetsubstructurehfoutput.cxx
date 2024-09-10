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
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename CandidateCollisionTable, typename CandidateMcCollisionTable, typename CandidateTable, typename CandidateTableMCD, typename CandidateTableMCP, typename TracksSub, typename JetTableData, typename JetMatchedTableData, typename OutputCollisionTableData, typename OutputTableData, typename SubstructureOutputTableData, typename MatchingOutputTableData, typename JetTableMCD, typename OutputCollisionTableMCD, typename OutputTableMCD, typename SubstructureOutputTableMCD, typename MatchingOutputTableMCD, typename JetTableMCP, typename OutputCollisionTableMCP, typename OutputTableMCP, typename SubstructureOutputTableMCP, typename MatchingOutputTableMCP, typename JetTableDataSub, typename OutputCollisionTableDataSub, typename OutputTableDataSub, typename SubstructureOutputTableDataSub, typename MatchingOutputTableDataSub, typename CandidateCollisionOutputTable, typename CandidateOutputTable, typename CandidateParOutputTable, typename CandidateParExtraOutputTable, typename CandidateSelOutputTable, typename CandidateMlOutputTable, typename CandidateMCDOutputTable, typename CandidateMcCollisionOutputTable, typename CandidateMcCollisionMatchingOutputTable, typename CandidateMCPOutputTable>
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
  Produces<CandidateMlOutputTable> candidateMlsTable;
  Produces<CandidateMCDOutputTable> candidateMcsTable;
  Produces<CandidateMcCollisionOutputTable> hfMcCollisionsTable;
  Produces<CandidateMcCollisionMatchingOutputTable> hfMcCollisionsMatchingTable;
  Produces<CandidateMCPOutputTable> hfParticlesTable;

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
  std::map<int32_t, int32_t> candidateMapping;
  std::map<int32_t, int32_t> candidateMappingMCP;
  std::map<int32_t, int32_t> candidateCollisionMapping;
  std::map<int32_t, int32_t> candidateMcCollisionMapping;

  std::vector<double> jetRadiiValues;

  std::vector<bool> collisionFlag;
  std::vector<bool> mcCollisionFlag;

  PresliceUnsorted<soa::Join<JetCollisions, aod::JMcCollisionLbs>> CollisionsPerMcCollision = aod::jmccollisionlb::mcCollisionId;
  PresliceOptional<CollisionsD0> D0CollisionsPerCollision = aod::jd0indices::collisionId;
  PresliceOptional<CollisionsLc> LcCollisionsPerCollision = aod::jlcindices::collisionId;
  PresliceOptional<CollisionsDielectron> DielectronCollisionsPerCollision = aod::jdielectronindices::collisionId;
  PresliceOptional<soa::Join<McCollisionsD0, aod::HfD0McRCollIds>> D0McCollisionsPerMcCollision = aod::jd0indices::mcCollisionId;
  PresliceOptional<soa::Join<McCollisionsLc, aod::Hf3PMcRCollIds>> LcMcCollisionsPerMcCollision = aod::jlcindices::mcCollisionId;
  PresliceOptional<McCollisionsDielectron> DielectronMcCollisionsPerMcCollision = aod::jdielectronindices::mcCollisionId;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  template <typename T, typename U, typename V, typename M>
  void fillJetTables(T const& jet, U const& /*cand*/, int32_t collisionIndex, int32_t candidateIndex, V& jetOutputTable, M& jetSubstructureOutputTable, std::map<int32_t, int32_t>& jetMap)
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
    jetOutputTable(collisionIndex, candidateIndex, jet.pt(), jet.phi(), jet.eta(), jet.y(), jet.r(), jet.tracksIds().size() + jet.candidatesIds().size()); // here we take the decision to keep the collision index consistent with the JE framework in case it is later needed to join to other tables. The candidate Index however can be linked to the HF tables
    jetSubstructureOutputTable(jetOutputTable.lastIndex(), energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, jet.nSub2DR(), jet.nSub1(), jet.nSub2(), pairPtVec, pairEnergyVec, pairThetaVec, jet.angularity());
    jetMap.insert(std::make_pair(jet.globalIndex(), jetOutputTable.lastIndex()));
  }

  template <bool isMCP, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& jets, V const& /*candidates*/, M& collisionOutputTable, N& jetOutputTable, O& jetSubstructureOutputTable, std::map<int32_t, int32_t>& jetMap, std::map<int32_t, int32_t>& candidateMap, float jetPtMin)
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
          auto candidate = jet.template candidates_first_as<V>();
          int32_t candidateIndex = -1;
          auto candidateTableIndex = candidateMap.find(candidate.globalIndex());
          if (candidateTableIndex != candidateMap.end()) {
            candidateIndex = candidateTableIndex->second;
          }
          if constexpr (!isMCP) {
            if (nJetInCollision == 0) {
              collisionOutputTable(collision.posZ(), collision.centrality(), collision.eventSel());
              collisionIndex = collisionOutputTable.lastIndex();
            }
            nJetInCollision++;
          }
          fillJetTables(jet, candidate, collisionIndex, candidateIndex, jetOutputTable, jetSubstructureOutputTable, jetMap);
        }
      }
    }
  }

  template <bool isMCD, bool isMCP, typename T, typename U>
  void analyseCandidates(T const& jets, U const& /*candidates*/, std::map<int32_t, int32_t>& candidateMap, float jetPtMin)
  {
    for (const auto& jet : jets) {
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {

          auto candidate = jet.template candidates_first_as<U>();

          auto candidateTableIndex = candidateMap.find(candidate.globalIndex());
          if (candidateTableIndex != candidateMap.end()) {
            continue;
          }
          int32_t candidateCollisionIndex = -1;
          int32_t candidateIndex = -1;
          if constexpr (isMCP) {
            auto hfMcCollisionIndex = candidateMcCollisionMapping.find(jetcandidateutilities::getMcCandidateCollisionId(candidate));
            if (hfMcCollisionIndex != candidateMcCollisionMapping.end()) {
              candidateCollisionIndex = hfMcCollisionIndex->second;
            }
            jetcandidateutilities::fillCandidateMcTable(candidate, candidateCollisionIndex, hfParticlesTable, candidateIndex);
          } else {
            auto hfCollisionIndex = candidateCollisionMapping.find(jetcandidateutilities::getCandidateCollisionId(candidate));
            if (hfCollisionIndex != candidateCollisionMapping.end()) {
              candidateCollisionIndex = hfCollisionIndex->second;
            }
            jetcandidateutilities::fillCandidateTable<isMCD>(candidate, candidateCollisionIndex, candidateTable, candidateParsTable, candidateParExtrasTable, candidateSelsTable, candidateMlsTable, candidateMcsTable, candidateIndex);
          }
          candidateMap.insert(std::make_pair(candidate.globalIndex(), candidateIndex));
        }
      }
    }
  }

  template <typename T, typename U, typename V>
  void analyseMatched(T const& jets, U const& /*jetsTag*/, std::map<int32_t, int32_t>& jetMapping, std::map<int32_t, int32_t>& jetTagMapping, V& matchingOutputTable, float jetPtMin)
  {
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

  template <bool isMC, bool isMCPOnly, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename S>
  void analyseHFCollisions(T const& collisions, U const& mcCollisions, V const& hfCollisions, M const& hfMcCollisions, N const& jets, O const& jetsMCP, P const& candidates, S const& candidatesMCP, float jetPtMin, float jetPtMinMCP = 0.0)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size());
    std::fill(collisionFlag.begin(), collisionFlag.end(), false);

    mcCollisionFlag.clear();
    mcCollisionFlag.resize(mcCollisions.size());
    std::fill(mcCollisionFlag.begin(), mcCollisionFlag.end(), false);

    if constexpr (!isMCPOnly) {
      for (const auto& jet : jets) {
        if (jet.pt() < jetPtMin) {
          continue;
        }
        if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        for (const auto& jetRadiiValue : jetRadiiValues) {
          if (jet.r() == round(jetRadiiValue * 100.0f)) {
            collisionFlag[jet.collisionId()] = true;
            if constexpr (isMC) {
              auto mcCollisionId = jet.template collision_as<T>().mcCollisionId();
              if (mcCollisionId >= 0) {
                mcCollisionFlag[mcCollisionId] = true;
              }
            }
          }
        }
      }
    }
    if constexpr (isMC) {
      for (const auto& jetMCP : jetsMCP) {
        if (jetMCP.pt() < jetPtMinMCP) {
          continue;
        }
        if (!jetfindingutilities::isInEtaAcceptance(jetMCP, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
          continue;
        }
        for (const auto& jetRadiiValue : jetRadiiValues) {
          if (jetMCP.r() == round(jetRadiiValue * 100.0f)) {

            mcCollisionFlag[jetMCP.mcCollisionId()] = true;
            if constexpr (!isMCPOnly) {
              const auto collisionsPerMcCollision = collisions.sliceBy(CollisionsPerMcCollision, jetMCP.mcCollisionId());
              for (auto collision : collisionsPerMcCollision) {
                collisionFlag[collision.globalIndex()] = true;
              }
            }
          }
        }
      }
    }
    if constexpr (!isMCPOnly) {
      for (const auto& collision : collisions) {
        if (collisionFlag[collision.globalIndex()]) {
          const auto hfCollisionsPerCollision = jetcandidateutilities::slicedPerCandidateCollision(hfCollisions, candidates, collision, D0CollisionsPerCollision, LcCollisionsPerCollision, D0CollisionsPerCollision, DielectronCollisionsPerCollision); // add Bplus later
          int32_t candidateCollisionIndex = -1;
          for (const auto& hfCollisionPerCollision : hfCollisionsPerCollision) { // should only ever be one
            auto hfCollisionTableIndex = candidateCollisionMapping.find(hfCollisionPerCollision.globalIndex());
            if (hfCollisionTableIndex != candidateCollisionMapping.end()) {
              continue;
            }
            jetcandidateutilities::fillCandidateCollisionTable(hfCollisionPerCollision, candidates, hfCollisionsTable, candidateCollisionIndex);
            candidateCollisionMapping.insert(std::make_pair(hfCollisionPerCollision.globalIndex(), hfCollisionsTable.lastIndex()));
          }
        }
      }
    }
    if constexpr (isMC) {
      for (const auto& mcCollision : mcCollisions) {
        if (mcCollisionFlag[mcCollision.globalIndex()]) {
          const auto hfMcCollisionsPerMcCollision = jetcandidateutilities::slicedPerCandidateCollision(hfMcCollisions, candidatesMCP, mcCollision, D0McCollisionsPerMcCollision, LcMcCollisionsPerMcCollision, D0McCollisionsPerMcCollision, DielectronMcCollisionsPerMcCollision); // add Bplus later
          int32_t candidateMcCollisionIndex = -1;
          for (const auto& hfMcCollisionPerMcCollision : hfMcCollisionsPerMcCollision) { // should only ever be one
            auto hfMcCollisionTableIndex = candidateMcCollisionMapping.find(hfMcCollisionPerMcCollision.globalIndex());
            if (hfMcCollisionTableIndex != candidateMcCollisionMapping.end()) {
              continue;
            }
            jetcandidateutilities::fillCandidateMcCollisionTable(hfMcCollisionPerMcCollision, candidatesMCP, hfMcCollisionsTable, candidateMcCollisionIndex);
            candidateMcCollisionMapping.insert(std::make_pair(hfMcCollisionPerMcCollision.globalIndex(), hfMcCollisionsTable.lastIndex()));
            if constexpr (!isMCPOnly && (jethfutilities::isHFTable<P>() || jethfutilities::isHFMcTable<S>())) { // the matching of mcCollision to Collision is only done for HF tables
              std::vector<int32_t> hfCollisionIDs;
              for (auto const& hfCollisionPerMcCollision : hfMcCollisionPerMcCollision.template hfCollBases_as<V>()) { // if added for others this line needs to be templated per type
                auto hfCollisionIndex = candidateCollisionMapping.find(hfCollisionPerMcCollision.globalIndex());
                if (hfCollisionIndex != candidateCollisionMapping.end()) {
                  hfCollisionIDs.push_back(hfCollisionIndex->second);
                }
              }
              hfMcCollisionsMatchingTable(hfCollisionIDs);
            }
          }
        }
      }
    }
  }

  void processClearMaps(JetCollisions const&)
  {
    candidateMapping.clear();
    candidateCollisionMapping.clear();
    candidateMappingMCP.clear();
    jetMappingData.clear();
    jetMappingDataSub.clear();
    jetMappingMCD.clear();
    jetMappingMCP.clear();
    candidateCollisionMapping.clear();
    candidateMcCollisionMapping.clear();
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processClearMaps, "process function that clears all the maps in each dataframe", true);

  void processOutputCollisionsData(JetCollisions const& collisions,
                                   JetTableData const& jets,
                                   CandidateCollisionTable const& canidateCollisions,
                                   CandidateTable const& candidates)
  {
    analyseHFCollisions<false, false>(collisions, collisions, canidateCollisions, canidateCollisions, jets, jets, candidates, candidates, jetPtMinData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsData, "hf collision output data", false);

  void processOutputCollisionsDataSub(JetCollisions const& collisions,
                                      JetTableDataSub const& jets,
                                      CandidateCollisionTable const& canidateCollisions,
                                      CandidateTable const& candidates)
  {
    analyseHFCollisions<false, false>(collisions, collisions, canidateCollisions, canidateCollisions, jets, jets, candidates, candidates, jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsDataSub, "hf collision output data eventwise constituent subtracted", false);

  void processOutputCollisionsMc(soa::Join<JetCollisions, aod::JMcCollisionLbs> const& collisions,
                                 JetMcCollisions const& mcCollisions,
                                 JetTableMCD const& jetsMCD,
                                 JetTableMCP const& jetsMCP,
                                 CandidateCollisionTable const& canidateCollisions,
                                 CandidateMcCollisionTable const& canidateMcCollisions,
                                 CandidateTableMCD const& candidatesMCD,
                                 CandidateTableMCP const& candidatesMCP)
  {
    analyseHFCollisions<true, false>(collisions, mcCollisions, canidateCollisions, canidateMcCollisions, jetsMCD, jetsMCP, candidatesMCD, candidatesMCP, jetPtMinMCD, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsMc, "hf collision output MC", false);

  void processOutputCollisionsMCPOnly(JetMcCollisions const& mcCollisions,
                                      JetTableMCP const& jetsMCP,
                                      CandidateMcCollisionTable const& canidateMcCollisions,
                                      CandidateTableMCP const& candidatesMCP)
  {
    analyseHFCollisions<true, true>(mcCollisions, mcCollisions, canidateMcCollisions, canidateMcCollisions, jetsMCP, jetsMCP, candidatesMCP, candidatesMCP, 0.0, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsMCPOnly, "hf collision output MCP only", false);

  void processOutputCandidatesData(JetCollision const&,
                                   JetTableData const& jets,
                                   CandidateTable const& candidates)
  {
    analyseCandidates<false, false>(jets, candidates, candidateMapping, jetPtMinData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesData, "hf candidate output data", false);

  void processOutputCandidatesDataSub(JetCollision const&,
                                      JetTableDataSub const& jets,
                                      CandidateTable const& candidates)
  {
    analyseCandidates<false, false>(jets, candidates, candidateMapping, jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesDataSub, "hf candidate output data eventwise constituent subtracted", false);

  void processOutputCandidatesMCD(JetCollision const&,
                                  JetTableMCD const& jets,
                                  CandidateTableMCD const& candidates)
  {

    analyseCandidates<true, false>(jets, candidates, candidateMapping, jetPtMinMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesMCD, "hf candidate output MCD", false);

  void processOutputCandidatesMCP(JetMcCollision const&,
                                  JetTableMCP const& jets,
                                  CandidateTableMCP const& candidates)
  {
    analyseCandidates<false, true>(jets, candidates, candidateMappingMCP, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesMCP, "hf candidate output MCP", false);

  void processOutputJetsData(JetCollision const& collision,
                             JetTableData const& jets,
                             CandidateTable const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, collisionOutputTableData, jetOutputTableData, jetSubstructureOutputTableData, jetMappingData, candidateMapping, jetPtMinData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsData, "hf jet substructure output Data", false);

  void processOutputJetsDataSub(JetCollision const& collision,
                                JetTableDataSub const& jets,
                                CandidateTable const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, collisionOutputTableDataSub, jetOutputTableDataSub, jetSubstructureOutputTableDataSub, jetMappingDataSub, candidateMapping, jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsDataSub, "hf jet substructure output event-wise subtracted Data", false);

  void processOutputMatchingData(JetMatchedTableData const& jets,
                                 JetTableDataSub const& jetsSub)
  {
    analyseMatched(jets, jetsSub, jetMappingData, jetMappingDataSub, jetMatchingOutputTableData, jetPtMinData);
    analyseMatched(jetsSub, jets, jetMappingDataSub, jetMappingData, jetMatchingOutputTableDataSub, jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMatchingData, "jet matching output Data", false);

  void processOutputJetsMCD(JetCollision const& collision,
                            JetTableMCD const& jets,
                            CandidateTableMCD const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, collisionOutputTableMCD, jetOutputTableMCD, jetSubstructureOutputTableMCD, jetMappingMCD, candidateMapping, jetPtMinMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsMCD, "hf jet substructure output MCD", false);

  void processOutputJetsMCP(JetMcCollision const& collision,
                            JetTableMCP const& jets,
                            CandidateTableMCP const& candidates)
  {
    analyseCharged<true>(collision, jets, candidates, collisionOutputTableMCP, jetOutputTableMCP, jetSubstructureOutputTableMCP, jetMappingMCP, candidateMappingMCP, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsMCP, "hf jet substructure output MCP", false);

  void processOutputMatchingMC(JetTableMCD const& jetsMCD,
                               JetTableMCP const& jetsMCP)
  {
    analyseMatched(jetsMCD, jetsMCP, jetMappingMCD, jetMappingMCP, jetMatchingOutputTableMCD, jetPtMinMCD);
    analyseMatched(jetsMCP, jetsMCD, jetMappingMCP, jetMappingMCD, jetMatchingOutputTableMCP, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputMatchingMC, "jet matching output MC", false);
};
using JetSubstructureOutputD0 = JetSubstructureHFOutputTask<CollisionsD0, soa::Join<McCollisionsD0, aod::HfD0McRCollIds>, CandidatesD0Data, CandidatesD0MCD, CandidatesD0MCP, aod::JTrackD0Subs, soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0CJetSSs>, soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents, aod::D0CJetSSs, aod::D0ChargedJetsMatchedToD0ChargedEventWiseSubtractedJets>, aod::D0CJetCOs, aod::D0CJetOs, aod::D0CJetSSOs, aod::D0CJetMOs, soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::D0CMCDJetSSs, aod::D0ChargedMCDetectorLevelJetsMatchedToD0ChargedMCParticleLevelJets>, aod::D0CMCDJetCOs, aod::D0CMCDJetOs, aod::D0CMCDJetSSOs, aod::D0CMCDJetMOs, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::D0CMCPJetSSs, aod::D0ChargedMCParticleLevelJetsMatchedToD0ChargedMCDetectorLevelJets>, aod::D0CMCPJetCOs, aod::D0CMCPJetOs, aod::D0CMCPJetSSOs, aod::D0CMCPJetMOs, soa::Join<aod::D0ChargedEventWiseSubtractedJets, aod::D0ChargedEventWiseSubtractedJetConstituents, aod::D0CEWSJetSSs, aod::D0ChargedEventWiseSubtractedJetsMatchedToD0ChargedJets>, aod::D0CEWSJetCOs, aod::D0CEWSJetOs, aod::D0CEWSJetSSOs, aod::D0CEWSJetMOs, aod::StoredHfD0CollBase, aod::StoredHfD0Bases, aod::StoredHfD0Pars, aod::StoredHfD0ParEs, aod::StoredHfD0Sels, aod::StoredHfD0Mls, aod::StoredHfD0Mcs, aod::StoredHfD0McCollBases, aod::StoredHfD0McRCollIds, aod::StoredHfD0PBases>;
using JetSubstructureOutputLc = JetSubstructureHFOutputTask<CollisionsLc, soa::Join<McCollisionsLc, aod::Hf3PMcRCollIds>, CandidatesLcData, CandidatesLcMCD, CandidatesLcMCP, aod::JTrackLcSubs, soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcCJetSSs>, soa::Join<aod::LcChargedJets, aod::LcChargedJetConstituents, aod::LcCJetSSs, aod::LcChargedJetsMatchedToLcChargedEventWiseSubtractedJets>, aod::LcCJetCOs, aod::LcCJetOs, aod::LcCJetSSOs, aod::LcCJetMOs, soa::Join<aod::LcChargedMCDetectorLevelJets, aod::LcChargedMCDetectorLevelJetConstituents, aod::LcCMCDJetSSs, aod::LcChargedMCDetectorLevelJetsMatchedToLcChargedMCParticleLevelJets>, aod::LcCMCDJetCOs, aod::LcCMCDJetOs, aod::LcCMCDJetSSOs, aod::LcCMCDJetMOs, soa::Join<aod::LcChargedMCParticleLevelJets, aod::LcChargedMCParticleLevelJetConstituents, aod::LcCMCPJetSSs, aod::LcChargedMCParticleLevelJetsMatchedToLcChargedMCDetectorLevelJets>, aod::LcCMCPJetCOs, aod::LcCMCPJetOs, aod::LcCMCPJetSSOs, aod::LcCMCPJetMOs, soa::Join<aod::LcChargedEventWiseSubtractedJets, aod::LcChargedEventWiseSubtractedJetConstituents, aod::LcCEWSJetSSs, aod::LcChargedEventWiseSubtractedJetsMatchedToLcChargedJets>, aod::LcCEWSJetCOs, aod::LcCEWSJetOs, aod::LcCEWSJetSSOs, aod::LcCEWSJetMOs, aod::StoredHf3PCollBase, aod::StoredHf3PBases, aod::StoredHf3PPars, aod::StoredHf3PParEs, aod::StoredHf3PSels, aod::StoredHf3PMls, aod::StoredHf3PMcs, aod::StoredHf3PMcCollBases, aod::StoredHf3PMcRCollIds, aod::StoredHf3PPBases>;
// using JetSubstructureOutputBplus = JetSubstructureHFOutputTask<aod::HfBplusCollBases, soa::Join<aod::StoredHfBplusMcCollBases, aod::JBplusMcCollisionIds, aod::HfBplusMcRCollIds>, CandidatesBplusData, CandidatesBplusMCD, CandidatesBplusMCP, aod::JTrackBplusSubs, soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs>, soa::Join<aod::BplusChargedJets, aod::BplusChargedJetConstituents, aod::BplusCJetSSs, aod::BplusChargedJetsMatchedToBplusChargedEventWiseSubtractedJets>, aod::BplusCJetCOs, aod::BplusCJetOs, aod::BplusCJetSSOs, aod::BplusCJetMOs, soa::Join<aod::BplusChargedMCDetectorLevelJets, aod::BplusChargedMCDetectorLevelJetConstituents, aod::BplusCMCDJetSSs, aod::BplusChargedMCDetectorLevelJetsMatchedToBplusChargedMCParticleLevelJets>, aod::BplusCMCDJetCOs, aod::BplusCMCDJetOs, aod::BplusCMCDJetSSOs,  aod::BplusCMCDJetMOs, soa::Join<aod::BplusChargedMCParticleLevelJets, aod::BplusChargedMCParticleLevelJetConstituents, aod::BplusCMCPJetSSs, aod::BplusChargedMCParticleLevelJetsMatchedToBplusChargedMCDetectorLevelJets>, aod::BplusCMCPJetCOs, aod::BplusCMCPJetOs, aod::BplusCMCPJetSSOs, aod::BplusCMCPJetMOs, soa::Join<aod::BplusChargedEventWiseSubtractedJets, aod::BplusChargedEventWiseSubtractedJetConstituents, aod::BplusCEWSJetSSs, aod::BplusChargedEventWiseSubtractedJetsMatchedToBplusChargedJets>, aod::BplusCEWSJetCOs, aod::BplusCEWSJetOs, aod::BplusCEWSJetSSOs, aod::BplusCEWSJetMOs, aod::StoredHfBplusCollBase, aod::StoredHfBplusBases, aod::StoredHfBplusPars, aod::StoredHfBplusParEs, aod::StoredHfBplusSels, aod::StoredHfBplusMls, aod::StoredHfBplusMcs, aod::StoredHfBplusPBases>;
using JetSubstructureOutputDielectron = JetSubstructureHFOutputTask<CollisionsDielectron, McCollisionsDielectron, CandidatesDielectronData, CandidatesDielectronMCD, CandidatesDielectronMCP, aod::JTrackDielectronSubs, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents, aod::DielectronCJetSSs>, soa::Join<aod::DielectronChargedJets, aod::DielectronChargedJetConstituents, aod::DielectronCJetSSs, aod::DielectronChargedJetsMatchedToDielectronChargedEventWiseSubtractedJets>, aod::DielectronCJetCOs, aod::DielectronCJetOs, aod::DielectronCJetSSOs, aod::DielectronCJetMOs, soa::Join<aod::DielectronChargedMCDetectorLevelJets, aod::DielectronChargedMCDetectorLevelJetConstituents, aod::DielectronCMCDJetSSs, aod::DielectronChargedMCDetectorLevelJetsMatchedToDielectronChargedMCParticleLevelJets>, aod::DielectronCMCDJetCOs, aod::DielectronCMCDJetOs, aod::DielectronCMCDJetSSOs, aod::DielectronCMCDJetMOs, soa::Join<aod::DielectronChargedMCParticleLevelJets, aod::DielectronChargedMCParticleLevelJetConstituents, aod::DielectronCMCPJetSSs, aod::DielectronChargedMCParticleLevelJetsMatchedToDielectronChargedMCDetectorLevelJets>, aod::DielectronCMCPJetCOs, aod::DielectronCMCPJetOs, aod::DielectronCMCPJetSSOs, aod::DielectronCMCPJetMOs, soa::Join<aod::DielectronChargedEventWiseSubtractedJets, aod::DielectronChargedEventWiseSubtractedJetConstituents, aod::DielectronCEWSJetSSs, aod::DielectronChargedEventWiseSubtractedJetsMatchedToDielectronChargedJets>, aod::DielectronCEWSJetCOs, aod::DielectronCEWSJetOs, aod::DielectronCEWSJetSSOs, aod::DielectronCEWSJetMOs, aod::StoredReducedEvents, aod::StoredDielectrons, aod::JDielectron1Dummys, aod::JDielectron2Dummys, aod::JDielectron3Dummys, aod::JDielectron4Dummys, aod::JDielectron5Dummys, aod::StoredJDielectronMcCollisions, aod::JDielectron6Dummys, aod::StoredJDielectronMcs>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputD0>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-d0-output"}));
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputLc>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-lc-output"}));
  // tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputBplus>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-bplus-output"}));
  tasks.emplace_back(adaptAnalysisTask<JetSubstructureOutputDielectron>(cfgc, SetDefaultProcesses{}, TaskName{"jet-substructure-dielectron-output"}));

  return WorkflowSpec{tasks};
}
