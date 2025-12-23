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

#ifndef PWGJE_TASKS_JETSUBSTRUCTUREHFOUTPUT_H_
#define PWGJE_TASKS_JETSUBSTRUCTUREHFOUTPUT_H_

#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataHF.h"
#include "PWGJE/DataModel/JetSubstructure.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

// NB: runDataProcessing.h must be included after customize!

template <typename CandidateCollisionTable, typename CandidateMcCollisionTable, typename CandidateMcOnlyCollisionTable, typename CandidateTable, typename CandidateTableMCD, typename CandidateTableMCP, typename CandidateRhosTable, typename CandidateMCRhosTable, typename TracksSub, typename JetTableData, typename JetMatchedTableData, typename RecoilTableData, typename SplittingTableData, typename PairTableData, typename OutputCollisionTableData, typename OutputTableData, typename SubstructureOutputTableData, typename MatchingOutputTableData, typename RecoilOutputTableData, typename JetTableMCD, typename RecoilTableMCD, typename SplittingTableMCD, typename PairTableMCD, typename OutputCollisionTableMCD, typename OutputTableMCD, typename SubstructureOutputTableMCD, typename MatchingOutputTableMCD, typename RecoilOutputTableMCD, typename MatchingRecoilOutputTableMCD, typename JetTableMCP, typename JetTableMatchedMCP, typename RecoilTableMCP, typename SplittingTableMCP, typename PairTableMCP, typename OutputCollisionTableMCP, typename CandidateMcOnlyCollisionOutputTable, typename OutputTableMCP, typename SubstructureOutputTableMCP, typename MatchingOutputTableMCP, typename RecoilOutputTableMCP, typename MatchingRecoilOutputTableMCP, typename JetTableDataSub, typename SplittingTableDataSub, typename PairTableDataSub, typename OutputCollisionTableDataSub, typename OutputTableDataSub, typename SubstructureOutputTableDataSub, typename MatchingOutputTableDataSub, typename CandidateCollisionOutputTable, typename CandidateOutputTable, typename CandidateParOutputTable, typename CandidateParExtraOutputTable, typename CandidateParDaughterOutputTable, typename CandidateSelOutputTable, typename CandidateMlOutputTable, typename CandidateMlDaughterOutputTable, typename CandidateMCDOutputTable, typename CandidateMcCollisionOutputTable, typename CandidateMcCollisionMatchingOutputTable, typename CandidateMCPOutputTable>
struct JetSubstructureHFOutputTask {

  struct : o2::framework::ProducesGroup {
    o2::framework::Produces<OutputCollisionTableData> collisionOutputTableData;
    o2::framework::Produces<OutputTableData> jetOutputTableData;
    o2::framework::Produces<SubstructureOutputTableData> jetSubstructureOutputTableData;
    o2::framework::Produces<RecoilOutputTableData> jetRecoilOutputTableData;
    o2::framework::Produces<MatchingOutputTableData> jetMatchingOutputTableData;
    o2::framework::Produces<OutputCollisionTableDataSub> collisionOutputTableDataSub;
    o2::framework::Produces<OutputTableDataSub> jetOutputTableDataSub;
    o2::framework::Produces<SubstructureOutputTableDataSub> jetSubstructureOutputTableDataSub;
    o2::framework::Produces<MatchingOutputTableDataSub> jetMatchingOutputTableDataSub;
    o2::framework::Produces<OutputCollisionTableMCD> collisionOutputTableMCD;
    o2::framework::Produces<OutputTableMCD> jetOutputTableMCD;
    o2::framework::Produces<SubstructureOutputTableMCD> jetSubstructureOutputTableMCD;
    o2::framework::Produces<RecoilOutputTableMCD> jetRecoilOutputTableMCD;
    o2::framework::Produces<MatchingOutputTableMCD> jetMatchingOutputTableMCD;
    o2::framework::Produces<MatchingRecoilOutputTableMCD> jetRecoilMatchingOutputTableMCD;
    o2::framework::Produces<OutputCollisionTableMCP> collisionOutputTableMCP;
    o2::framework::Produces<CandidateMcOnlyCollisionOutputTable> hfMcOnlyCollisionsTable;
    o2::framework::Produces<OutputTableMCP> jetOutputTableMCP;
    o2::framework::Produces<SubstructureOutputTableMCP> jetSubstructureOutputTableMCP;
    o2::framework::Produces<RecoilOutputTableMCP> jetRecoilOutputTableMCP;
    o2::framework::Produces<MatchingOutputTableMCP> jetMatchingOutputTableMCP;
    o2::framework::Produces<MatchingRecoilOutputTableMCP> jetRecoilMatchingOutputTableMCP;
    o2::framework::Produces<CandidateCollisionOutputTable> hfCollisionsTable;
    o2::framework::Produces<CandidateOutputTable> candidateTable;
    o2::framework::Produces<CandidateParOutputTable> candidateParsTable;
    o2::framework::Produces<CandidateParExtraOutputTable> candidateParExtrasTable;
    o2::framework::Produces<CandidateParDaughterOutputTable> candidateParsDaughterTable;
    o2::framework::Produces<CandidateSelOutputTable> candidateSelsTable;
    o2::framework::Produces<CandidateMlOutputTable> candidateMlsTable;
    o2::framework::Produces<CandidateMlDaughterOutputTable> candidateMlsDaughterTable;
    o2::framework::Produces<CandidateMCDOutputTable> candidateMcsTable;
    o2::framework::Produces<CandidateMcCollisionOutputTable> hfMcCollisionsTable;
    o2::framework::Produces<CandidateMcCollisionMatchingOutputTable> hfMcCollisionsMatchingTable;
    o2::framework::Produces<CandidateMCPOutputTable> hfParticlesTable;
  } products;

  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<float> jetPtMinData{"jetPtMinData", 0.0, "minimum jet pT cut for data jets"};
    o2::framework::Configurable<float> jetPtMinDataSub{"jetPtMinDataSub", 0.0, "minimum jet pT cut for eventwise constituent subtracted data jets"};
    o2::framework::Configurable<float> jetPtMinMCD{"jetPtMinMCD", 0.0, "minimum jet pT cut for mcd jets"};
    o2::framework::Configurable<float> jetPtMinMCP{"jetPtMinMCP", 0.0, "minimum jet pT cut for mcp jets"};
    o2::framework::Configurable<float> recoilJetPtMinData{"recoilJetPtMinData", 0.0, "minimum jet pT cut for data recoil jets"};
    o2::framework::Configurable<float> recoilJetPtMinMCD{"recoilJetPtMinMCD", 0.0, "minimum jet pT cut for mcd recoil jets"};
    o2::framework::Configurable<float> recoilJetPtMinMCP{"recoilJetPtMinMCP", 0.0, "minimum jet pT cut for mcp recoil jets"};
    o2::framework::Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
    o2::framework::Configurable<std::vector<double>> recoilJetRadii{"recoilJetRadii", std::vector<double>{0.4}, "recoil jet resolution parameters"};
    o2::framework::Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
    o2::framework::Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
    o2::framework::Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track pseudorapidity"};
    o2::framework::Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track pseudorapidity"};
  } configs;

  // need to add selection on pThat to post processing

  std::map<int32_t, int32_t> jetMappingData;
  std::map<int32_t, int32_t> jetMappingDataSub;
  std::map<int32_t, int32_t> jetMappingMCD;
  std::map<int32_t, int32_t> jetMappingMCP;
  std::map<int32_t, int32_t> candidateMapping;
  std::map<int32_t, int32_t> candidateMappingMCP;
  std::map<int32_t, int32_t> candidateCollisionMapping;
  std::map<int32_t, int32_t> candidateMcCollisionMapping;
  std::map<int32_t, int32_t> recoilJetMappingData; // not doing anything yet
  std::map<int32_t, int32_t> recoilJetMappingMCD;
  std::map<int32_t, int32_t> recoilJetMappingMCP;

  std::vector<bool> candidateSelectionFlagsData;
  std::vector<bool> candidateSelectionFlagsMCD;
  std::vector<bool> candidateSelectionFlagsMCP;
  std::vector<bool> candidateRecoilSelectionFlagsData;
  std::vector<bool> candidateRecoilSelectionFlagsMCD;
  std::vector<bool> candidateRecoilSelectionFlagsMCP;

  std::vector<std::vector<int32_t>> splittingMatchesGeoVecVecData;
  std::vector<std::vector<int32_t>> splittingMatchesPtVecVecData;
  std::vector<std::vector<int32_t>> splittingMatchesHFVecVecData;
  std::vector<std::vector<int32_t>> splittingMatchesGeoVecVecDataSub;
  std::vector<std::vector<int32_t>> splittingMatchesPtVecVecDataSub;
  std::vector<std::vector<int32_t>> splittingMatchesHFVecVecDataSub;
  std::vector<std::vector<int32_t>> splittingMatchesGeoVecVecMCD;
  std::vector<std::vector<int32_t>> splittingMatchesPtVecVecMCD;
  std::vector<std::vector<int32_t>> splittingMatchesHFVecVecMCD;
  std::vector<std::vector<int32_t>> splittingMatchesGeoVecVecMCP;
  std::vector<std::vector<int32_t>> splittingMatchesPtVecVecMCP;
  std::vector<std::vector<int32_t>> splittingMatchesHFVecVecMCP;

  std::vector<std::vector<int32_t>> pairMatchesVecVecData;
  std::vector<std::vector<int32_t>> pairMatchesVecVecDataSub;
  std::vector<std::vector<int32_t>> pairMatchesVecVecMCD;
  std::vector<std::vector<int32_t>> pairMatchesVecVecMCP;

  std::vector<double> jetRadiiValues;
  std::vector<double> recoilJetRadiiValues;

  std::vector<bool> collisionFlag;
  std::vector<bool> mcCollisionFlag;

  void init(o2::framework::InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)configs.jetRadii;
    recoilJetRadiiValues = (std::vector<double>)configs.recoilJetRadii;
  }

  struct : o2::framework::PresliceGroup {
    o2::framework::PresliceUnsorted<o2::soa::Join<o2::aod::JetCollisions, o2::aod::JMcCollisionLbs>> CollisionsPerMcCollision = o2::aod::jmccollisionlb::mcCollisionId;
    o2::framework::PresliceOptional<CandidateCollisionTable> CandidateCollisionsPerCollision = o2::aod::jcandidateindices::collisionId;
    o2::framework::PresliceOptional<CandidateMcCollisionTable> CandidateMcCollisionsPerMcCollision = o2::aod::jcandidateindices::mcCollisionId;
    o2::framework::PresliceOptional<CandidateMcOnlyCollisionTable> CandidateMcCollisionsPerMcCollisionMCPOnly = o2::aod::jcandidateindices::mcCollisionId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedSPs, o2::aod::D0ChargedSPsMatchedToD0ChargedEventWiseSubtractedSPs>> D0SplittingsPerJetData = o2::aod::d0chargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedEventWiseSubtractedSPs, o2::aod::D0ChargedEventWiseSubtractedSPsMatchedToD0ChargedSPs>> D0SplittingsPerJetDataSub = o2::aod::d0chargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelSPs, o2::aod::D0ChargedMCDetectorLevelSPsMatchedToD0ChargedMCParticleLevelSPs>> D0SplittingsPerJetMCD = o2::aod::d0chargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedMCParticleLevelSPs, o2::aod::D0ChargedMCParticleLevelSPsMatchedToD0ChargedMCDetectorLevelSPs>> D0SplittingsPerJetMCP = o2::aod::d0chargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedPRs, o2::aod::D0ChargedPRsMatchedToD0ChargedEventWiseSubtractedPRs>> D0PairsPerJetData = o2::aod::d0chargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedEventWiseSubtractedPRs, o2::aod::D0ChargedEventWiseSubtractedPRsMatchedToD0ChargedPRs>> D0PairsPerJetDataSub = o2::aod::d0chargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedMCDetectorLevelPRs, o2::aod::D0ChargedMCDetectorLevelPRsMatchedToD0ChargedMCParticleLevelPRs>> D0PairsPerJetMCD = o2::aod::d0chargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::D0ChargedMCParticleLevelPRs, o2::aod::D0ChargedMCParticleLevelPRsMatchedToD0ChargedMCDetectorLevelPRs>> D0PairsPerJetMCP = o2::aod::d0chargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedSPs, o2::aod::DplusChargedSPsMatchedToDplusChargedEventWiseSubtractedSPs>> DplusSplittingsPerJetData = o2::aod::dpluschargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedEventWiseSubtractedSPs, o2::aod::DplusChargedEventWiseSubtractedSPsMatchedToDplusChargedSPs>> DplusSplittingsPerJetDataSub = o2::aod::dpluschargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelSPs, o2::aod::DplusChargedMCDetectorLevelSPsMatchedToDplusChargedMCParticleLevelSPs>> DplusSplittingsPerJetMCD = o2::aod::dpluschargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedMCParticleLevelSPs, o2::aod::DplusChargedMCParticleLevelSPsMatchedToDplusChargedMCDetectorLevelSPs>> DplusSplittingsPerJetMCP = o2::aod::dpluschargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedPRs, o2::aod::DplusChargedPRsMatchedToDplusChargedEventWiseSubtractedPRs>> DplusPairsPerJetData = o2::aod::dpluschargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedEventWiseSubtractedPRs, o2::aod::DplusChargedEventWiseSubtractedPRsMatchedToDplusChargedPRs>> DplusPairsPerJetDataSub = o2::aod::dpluschargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedMCDetectorLevelPRs, o2::aod::DplusChargedMCDetectorLevelPRsMatchedToDplusChargedMCParticleLevelPRs>> DplusPairsPerJetMCD = o2::aod::dpluschargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DplusChargedMCParticleLevelPRs, o2::aod::DplusChargedMCParticleLevelPRsMatchedToDplusChargedMCDetectorLevelPRs>> DplusPairsPerJetMCP = o2::aod::dpluschargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedSPs, o2::aod::DsChargedSPsMatchedToDsChargedEventWiseSubtractedSPs>> DsSplittingsPerJetData = o2::aod::dschargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedEventWiseSubtractedSPs, o2::aod::DsChargedEventWiseSubtractedSPsMatchedToDsChargedSPs>> DsSplittingsPerJetDataSub = o2::aod::dschargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedMCDetectorLevelSPs, o2::aod::DsChargedMCDetectorLevelSPsMatchedToDsChargedMCParticleLevelSPs>> DsSplittingsPerJetMCD = o2::aod::dschargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedMCParticleLevelSPs, o2::aod::DsChargedMCParticleLevelSPsMatchedToDsChargedMCDetectorLevelSPs>> DsSplittingsPerJetMCP = o2::aod::dschargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedPRs, o2::aod::DsChargedPRsMatchedToDsChargedEventWiseSubtractedPRs>> DsPairsPerJetData = o2::aod::dschargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedEventWiseSubtractedPRs, o2::aod::DsChargedEventWiseSubtractedPRsMatchedToDsChargedPRs>> DsPairsPerJetDataSub = o2::aod::dschargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedMCDetectorLevelPRs, o2::aod::DsChargedMCDetectorLevelPRsMatchedToDsChargedMCParticleLevelPRs>> DsPairsPerJetMCD = o2::aod::dschargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DsChargedMCParticleLevelPRs, o2::aod::DsChargedMCParticleLevelPRsMatchedToDsChargedMCDetectorLevelPRs>> DsPairsPerJetMCP = o2::aod::dschargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedSPs, o2::aod::DstarChargedSPsMatchedToDstarChargedEventWiseSubtractedSPs>> DstarSplittingsPerJetData = o2::aod::dstarchargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedEventWiseSubtractedSPs, o2::aod::DstarChargedEventWiseSubtractedSPsMatchedToDstarChargedSPs>> DstarSplittingsPerJetDataSub = o2::aod::dstarchargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelSPs, o2::aod::DstarChargedMCDetectorLevelSPsMatchedToDstarChargedMCParticleLevelSPs>> DstarSplittingsPerJetMCD = o2::aod::dstarchargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedMCParticleLevelSPs, o2::aod::DstarChargedMCParticleLevelSPsMatchedToDstarChargedMCDetectorLevelSPs>> DstarSplittingsPerJetMCP = o2::aod::dstarchargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedPRs, o2::aod::DstarChargedPRsMatchedToDstarChargedEventWiseSubtractedPRs>> DstarPairsPerJetData = o2::aod::dstarchargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedEventWiseSubtractedPRs, o2::aod::DstarChargedEventWiseSubtractedPRsMatchedToDstarChargedPRs>> DstarPairsPerJetDataSub = o2::aod::dstarchargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedMCDetectorLevelPRs, o2::aod::DstarChargedMCDetectorLevelPRsMatchedToDstarChargedMCParticleLevelPRs>> DstarPairsPerJetMCD = o2::aod::dstarchargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DstarChargedMCParticleLevelPRs, o2::aod::DstarChargedMCParticleLevelPRsMatchedToDstarChargedMCDetectorLevelPRs>> DstarPairsPerJetMCP = o2::aod::dstarchargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedSPs, o2::aod::LcChargedSPsMatchedToLcChargedEventWiseSubtractedSPs>> LcSplittingsPerJetData = o2::aod::lcchargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedEventWiseSubtractedSPs, o2::aod::LcChargedEventWiseSubtractedSPsMatchedToLcChargedSPs>> LcSplittingsPerJetDataSub = o2::aod::lcchargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedMCDetectorLevelSPs, o2::aod::LcChargedMCDetectorLevelSPsMatchedToLcChargedMCParticleLevelSPs>> LcSplittingsPerJetMCD = o2::aod::lcchargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedMCParticleLevelSPs, o2::aod::LcChargedMCParticleLevelSPsMatchedToLcChargedMCDetectorLevelSPs>> LcSplittingsPerJetMCP = o2::aod::lcchargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedPRs, o2::aod::LcChargedPRsMatchedToLcChargedEventWiseSubtractedPRs>> LcPairsPerJetData = o2::aod::lcchargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedEventWiseSubtractedPRs, o2::aod::LcChargedEventWiseSubtractedPRsMatchedToLcChargedPRs>> LcPairsPerJetDataSub = o2::aod::lcchargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedMCDetectorLevelPRs, o2::aod::LcChargedMCDetectorLevelPRsMatchedToLcChargedMCParticleLevelPRs>> LcPairsPerJetMCD = o2::aod::lcchargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::LcChargedMCParticleLevelPRs, o2::aod::LcChargedMCParticleLevelPRsMatchedToLcChargedMCDetectorLevelPRs>> LcPairsPerJetMCP = o2::aod::lcchargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedSPs, o2::aod::B0ChargedSPsMatchedToB0ChargedEventWiseSubtractedSPs>> B0SplittingsPerJetData = o2::aod::b0chargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedEventWiseSubtractedSPs, o2::aod::B0ChargedEventWiseSubtractedSPsMatchedToB0ChargedSPs>> B0SplittingsPerJetDataSub = o2::aod::b0chargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelSPs, o2::aod::B0ChargedMCDetectorLevelSPsMatchedToB0ChargedMCParticleLevelSPs>> B0SplittingsPerJetMCD = o2::aod::b0chargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedMCParticleLevelSPs, o2::aod::B0ChargedMCParticleLevelSPsMatchedToB0ChargedMCDetectorLevelSPs>> B0SplittingsPerJetMCP = o2::aod::b0chargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedPRs, o2::aod::B0ChargedPRsMatchedToB0ChargedEventWiseSubtractedPRs>> B0PairsPerJetData = o2::aod::b0chargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedEventWiseSubtractedPRs, o2::aod::B0ChargedEventWiseSubtractedPRsMatchedToB0ChargedPRs>> B0PairsPerJetDataSub = o2::aod::b0chargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedMCDetectorLevelPRs, o2::aod::B0ChargedMCDetectorLevelPRsMatchedToB0ChargedMCParticleLevelPRs>> B0PairsPerJetMCD = o2::aod::b0chargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::B0ChargedMCParticleLevelPRs, o2::aod::B0ChargedMCParticleLevelPRsMatchedToB0ChargedMCDetectorLevelPRs>> B0PairsPerJetMCP = o2::aod::b0chargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedSPs, o2::aod::BplusChargedSPsMatchedToBplusChargedEventWiseSubtractedSPs>> BplusSplittingsPerJetData = o2::aod::bpluschargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedEventWiseSubtractedSPs, o2::aod::BplusChargedEventWiseSubtractedSPsMatchedToBplusChargedSPs>> BplusSplittingsPerJetDataSub = o2::aod::bpluschargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelSPs, o2::aod::BplusChargedMCDetectorLevelSPsMatchedToBplusChargedMCParticleLevelSPs>> BplusSplittingsPerJetMCD = o2::aod::bpluschargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedMCParticleLevelSPs, o2::aod::BplusChargedMCParticleLevelSPsMatchedToBplusChargedMCDetectorLevelSPs>> BplusSplittingsPerJetMCP = o2::aod::bpluschargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedPRs, o2::aod::BplusChargedPRsMatchedToBplusChargedEventWiseSubtractedPRs>> BplusPairsPerJetData = o2::aod::bpluschargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedEventWiseSubtractedPRs, o2::aod::BplusChargedEventWiseSubtractedPRsMatchedToBplusChargedPRs>> BplusPairsPerJetDataSub = o2::aod::bpluschargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedMCDetectorLevelPRs, o2::aod::BplusChargedMCDetectorLevelPRsMatchedToBplusChargedMCParticleLevelPRs>> BplusPairsPerJetMCD = o2::aod::bpluschargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::BplusChargedMCParticleLevelPRs, o2::aod::BplusChargedMCParticleLevelPRsMatchedToBplusChargedMCDetectorLevelPRs>> BplusPairsPerJetMCP = o2::aod::bpluschargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedSPs, o2::aod::XicToXiPiPiChargedSPsMatchedToXicToXiPiPiChargedEventWiseSubtractedSPs>> XicToXiPiPiSplittingsPerJetData = o2::aod::xictoxipipichargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedEventWiseSubtractedSPs, o2::aod::XicToXiPiPiChargedEventWiseSubtractedSPsMatchedToXicToXiPiPiChargedSPs>> XicToXiPiPiSplittingsPerJetDataSub = o2::aod::xictoxipipichargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelSPs, o2::aod::XicToXiPiPiChargedMCDetectorLevelSPsMatchedToXicToXiPiPiChargedMCParticleLevelSPs>> XicToXiPiPiSplittingsPerJetMCD = o2::aod::xictoxipipichargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCParticleLevelSPs, o2::aod::XicToXiPiPiChargedMCParticleLevelSPsMatchedToXicToXiPiPiChargedMCDetectorLevelSPs>> XicToXiPiPiSplittingsPerJetMCP = o2::aod::xictoxipipichargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedPRs, o2::aod::XicToXiPiPiChargedPRsMatchedToXicToXiPiPiChargedEventWiseSubtractedPRs>> XicToXiPiPiPairsPerJetData = o2::aod::xictoxipipichargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedEventWiseSubtractedPRs, o2::aod::XicToXiPiPiChargedEventWiseSubtractedPRsMatchedToXicToXiPiPiChargedPRs>> XicToXiPiPiPairsPerJetDataSub = o2::aod::xictoxipipichargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCDetectorLevelPRs, o2::aod::XicToXiPiPiChargedMCDetectorLevelPRsMatchedToXicToXiPiPiChargedMCParticleLevelPRs>> XicToXiPiPiPairsPerJetMCD = o2::aod::xictoxipipichargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::XicToXiPiPiChargedMCParticleLevelPRs, o2::aod::XicToXiPiPiChargedMCParticleLevelPRsMatchedToXicToXiPiPiChargedMCDetectorLevelPRs>> XicToXiPiPiPairsPerJetMCP = o2::aod::xictoxipipichargedmcparticlelevelpair::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedSPs, o2::aod::DielectronChargedSPsMatchedToDielectronChargedEventWiseSubtractedSPs>> DielectronSplittingsPerJetData = o2::aod::dielectronchargedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedEventWiseSubtractedSPs, o2::aod::DielectronChargedEventWiseSubtractedSPsMatchedToDielectronChargedSPs>> DielectronSplittingsPerJetDataSub = o2::aod::dielectronchargedeventwisesubtractedsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelSPs, o2::aod::DielectronChargedMCDetectorLevelSPsMatchedToDielectronChargedMCParticleLevelSPs>> DielectronSplittingsPerJetMCD = o2::aod::dielectronchargedmcdetectorlevelsplitting::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedMCParticleLevelSPs, o2::aod::DielectronChargedMCParticleLevelSPsMatchedToDielectronChargedMCDetectorLevelSPs>> DielectronSplittingsPerJetMCP = o2::aod::dielectronchargedmcparticlelevelsplitting::jetId;

    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedPRs, o2::aod::DielectronChargedPRsMatchedToDielectronChargedEventWiseSubtractedPRs>> DielectronPairsPerJetData = o2::aod::dielectronchargedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedEventWiseSubtractedPRs, o2::aod::DielectronChargedEventWiseSubtractedPRsMatchedToDielectronChargedPRs>> DielectronPairsPerJetDataSub = o2::aod::dielectronchargedeventwisesubtractedpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedMCDetectorLevelPRs, o2::aod::DielectronChargedMCDetectorLevelPRsMatchedToDielectronChargedMCParticleLevelPRs>> DielectronPairsPerJetMCD = o2::aod::dielectronchargedmcdetectorlevelpair::jetId;
    o2::framework::PresliceOptional<o2::soa::Join<o2::aod::DielectronChargedMCParticleLevelPRs, o2::aod::DielectronChargedMCParticleLevelPRsMatchedToDielectronChargedMCDetectorLevelPRs>> DielectronPairsPerJetMCP = o2::aod::dielectronchargedmcparticlelevelpair::jetId;

  } preslices;

  template <bool isMCP, typename T, typename U>
  auto candidateMCCollisionSlicer(T const& McCollisionsPerMcCollision, U const& McCollisionsPerMcCollisionMCPOnly)
  {
    if constexpr (isMCP) {
      return McCollisionsPerMcCollisionMCPOnly;
    } else {
      return McCollisionsPerMcCollision;
    }
  }

  template <typename T>
  void fillSplittingMatchingVectors(T const& splittingsMatches, int jetIndex, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec)
  {
    for (auto const& splittingMatches : splittingsMatches) {
      auto splittingMatchesGeoSpan = splittingMatches.splittingMatchingGeo();
      auto splittingMatchesPtSpan = splittingMatches.splittingMatchingPt();
      auto splittingMatchesHFSpan = splittingMatches.splittingMatchingHF();
      std::copy(splittingMatchesGeoSpan.begin(), splittingMatchesGeoSpan.end(), std::back_inserter(splittingMatchesGeoVecVec[jetIndex]));
      std::copy(splittingMatchesPtSpan.begin(), splittingMatchesPtSpan.end(), std::back_inserter(splittingMatchesPtVecVec[jetIndex]));
      std::copy(splittingMatchesHFSpan.begin(), splittingMatchesHFSpan.end(), std::back_inserter(splittingMatchesHFVecVec[jetIndex]));
    }
  }

  template <typename T>
  void fillPairMatchingVectors(T const& pairsMatches, int jetIndex, std::vector<std::vector<int32_t>>& pairMatchesVecVec)
  {
    for (auto const& pairMatches : pairsMatches) {
      auto pairMatchesSpan = pairMatches.pairMatching();
      std::copy(pairMatchesSpan.begin(), pairMatchesSpan.end(), std::back_inserter(pairMatchesVecVec[jetIndex]));
    }
  }

  template <typename T, typename U, typename V>
  void fillJetTables(T const& jet, int32_t collisionIndex, std::vector<int32_t> const& candidatesIndices, U& jetOutputTable, V& jetSubstructureOutputTable, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, float rho, std::map<int32_t, int32_t>& jetMap)
  {
    std::vector<float> energyMotherVec;
    std::vector<float> ptLeadingVec;
    std::vector<float> ptSubLeadingVec;
    std::vector<float> thetaVec;
    std::vector<float> pairJetPtVec;
    std::vector<float> pairJetEnergyVec;
    std::vector<float> pairJetThetaVec;
    std::vector<float> pairJetPerpCone1PtVec;
    std::vector<float> pairJetPerpCone1EnergyVec;
    std::vector<float> pairJetPerpCone1ThetaVec;
    std::vector<float> pairPerpCone1PerpCone1PtVec;
    std::vector<float> pairPerpCone1PerpCone1EnergyVec;
    std::vector<float> pairPerpCone1PerpCone1ThetaVec;
    std::vector<float> pairPerpCone1PerpCone2PtVec;
    std::vector<float> pairPerpCone1PerpCone2EnergyVec;
    std::vector<float> pairPerpCone1PerpCone2ThetaVec;

    auto energyMotherSpan = jet.energyMother();
    auto ptLeadingSpan = jet.ptLeading();
    auto ptSubLeadingSpan = jet.ptSubLeading();
    auto thetaSpan = jet.theta();
    auto pairJetPtSpan = jet.pairJetPt();
    auto pairJetEnergySpan = jet.pairJetEnergy();
    auto pairJetThetaSpan = jet.pairJetTheta();
    auto pairJetPerpCone1PtSpan = jet.pairJetPerpCone1Pt();
    auto pairJetPerpCone1EnergySpan = jet.pairJetPerpCone1Energy();
    auto pairJetPerpCone1ThetaSpan = jet.pairJetPerpCone1Theta();
    auto pairPerpCone1PerpCone1PtSpan = jet.pairPerpCone1PerpCone1Pt();
    auto pairPerpCone1PerpCone1EnergySpan = jet.pairPerpCone1PerpCone1Energy();
    auto pairPerpCone1PerpCone1ThetaSpan = jet.pairPerpCone1PerpCone1Theta();
    auto pairPerpCone1PerpCone2PtSpan = jet.pairPerpCone1PerpCone2Pt();
    auto pairPerpCone1PerpCone2EnergySpan = jet.pairPerpCone1PerpCone2Energy();
    auto pairPerpCone1PerpCone2ThetaSpan = jet.pairPerpCone1PerpCone2Theta();

    std::copy(energyMotherSpan.begin(), energyMotherSpan.end(), std::back_inserter(energyMotherVec));
    std::copy(ptLeadingSpan.begin(), ptLeadingSpan.end(), std::back_inserter(ptLeadingVec));
    std::copy(ptSubLeadingSpan.begin(), ptSubLeadingSpan.end(), std::back_inserter(ptSubLeadingVec));
    std::copy(thetaSpan.begin(), thetaSpan.end(), std::back_inserter(thetaVec));
    std::copy(pairJetPtSpan.begin(), pairJetPtSpan.end(), std::back_inserter(pairJetPtVec));
    std::copy(pairJetEnergySpan.begin(), pairJetEnergySpan.end(), std::back_inserter(pairJetEnergyVec));
    std::copy(pairJetThetaSpan.begin(), pairJetThetaSpan.end(), std::back_inserter(pairJetThetaVec));
    std::copy(pairJetPerpCone1PtSpan.begin(), pairJetPerpCone1PtSpan.end(), std::back_inserter(pairJetPerpCone1PtVec));
    std::copy(pairJetPerpCone1EnergySpan.begin(), pairJetPerpCone1EnergySpan.end(), std::back_inserter(pairJetPerpCone1EnergyVec));
    std::copy(pairJetPerpCone1ThetaSpan.begin(), pairJetPerpCone1ThetaSpan.end(), std::back_inserter(pairJetPerpCone1ThetaVec));
    std::copy(pairPerpCone1PerpCone1PtSpan.begin(), pairPerpCone1PerpCone1PtSpan.end(), std::back_inserter(pairPerpCone1PerpCone1PtVec));
    std::copy(pairPerpCone1PerpCone1EnergySpan.begin(), pairPerpCone1PerpCone1EnergySpan.end(), std::back_inserter(pairPerpCone1PerpCone1EnergyVec));
    std::copy(pairPerpCone1PerpCone1ThetaSpan.begin(), pairPerpCone1PerpCone1ThetaSpan.end(), std::back_inserter(pairPerpCone1PerpCone1ThetaVec));
    std::copy(pairPerpCone1PerpCone2PtSpan.begin(), pairPerpCone1PerpCone2PtSpan.end(), std::back_inserter(pairPerpCone1PerpCone2PtVec));
    std::copy(pairPerpCone1PerpCone2EnergySpan.begin(), pairPerpCone1PerpCone2EnergySpan.end(), std::back_inserter(pairPerpCone1PerpCone2EnergyVec));
    std::copy(pairPerpCone1PerpCone2ThetaSpan.begin(), pairPerpCone1PerpCone2ThetaSpan.end(), std::back_inserter(pairPerpCone1PerpCone2ThetaVec));

    std::vector<int> splittingMatchesGeoVec;
    std::vector<int> splittingMatchesPtVec;
    std::vector<int> splittingMatchesHFVec;
    std::vector<int> pairMatchesVec;
    if (doprocessOutputSubstructureMatchingData || doprocessOutputSubstructureMatchingMC) {
      splittingMatchesGeoVec = splittingMatchesGeoVecVec[jet.globalIndex()];
      splittingMatchesPtVec = splittingMatchesPtVecVec[jet.globalIndex()];
      splittingMatchesHFVec = splittingMatchesHFVecVec[jet.globalIndex()];
      pairMatchesVec = pairMatchesVecVec[jet.globalIndex()];
    }

    jetOutputTable(collisionIndex, candidatesIndices, jet.pt(), jet.phi(), jet.eta(), jet.y(), jet.r(), jet.area(), rho, jet.perpConeRho(), jet.tracksIds().size() + jet.candidatesIds().size()); // here we take the decision to keep the collision index consistent with the JE framework in case it is later needed to join to other tables. The candidate Index however can be linked to the HF tables
    jetSubstructureOutputTable(jetOutputTable.lastIndex(), energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, jet.nSub2DR(), jet.nSub1(), jet.nSub2(), pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, jet.angularity(), jet.ptLeadingConstituent(), splittingMatchesGeoVec, splittingMatchesPtVec, splittingMatchesHFVec, pairMatchesVec);
    jetMap.insert(std::make_pair(jet.globalIndex(), jetOutputTable.lastIndex()));
  }

  template <bool isMCP, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& jets, V const& /*candidates*/, M& collisionOutputTable, N& jetOutputTable, O& jetSubstructureOutputTable, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, std::map<int32_t, int32_t>& jetMap, std::map<int32_t, int32_t>& candidateMap, float jetPtMin, float eventWeight)
  {

    int nJetInCollision = 0;
    int32_t collisionIndex = -1;
    float rho = 0.0;
    for (const auto& jet : jets) {
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          std::vector<int32_t> candidatesIndices;
          for (auto const& candidate : jet.template candidates_as<V>()) {
            auto candidateTableIndex = candidateMap.find(candidate.globalIndex());
            if (candidateTableIndex != candidateMap.end()) {
              candidatesIndices.push_back(candidateTableIndex->second);
            }
            rho = candidate.rho();
          }
          if (nJetInCollision == 0) {
            float centrality = -1.0;
            uint8_t eventSel = 0.0;
            if constexpr (!isMCP) {
              centrality = collision.centFT0M();
              eventSel = collision.eventSel();
            }
            collisionOutputTable(collision.posZ(), centrality, eventSel, eventWeight);
            collisionIndex = collisionOutputTable.lastIndex();
          }
          nJetInCollision++;
          fillJetTables(jet, collisionIndex, candidatesIndices, jetOutputTable, jetSubstructureOutputTable, splittingMatchesGeoVecVec, splittingMatchesPtVecVec, splittingMatchesHFVecVec, pairMatchesVecVec, rho, jetMap);
        }
      }
    }
  }

  template <typename T, typename U>
  void analyseRecoilCharged(T const& jet, U& jetOutputTable, std::map<int32_t, int32_t>& jetMap, std::map<int32_t, int32_t>& candidateMap, float jetPtMin)
  {
    if (jet.pt() < jetPtMin) {
      return;
    }
    if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
      return;
    }
    for (const auto& jetRadiiValue : recoilJetRadiiValues) {
      if (jet.r() == round(jetRadiiValue * 100.0f)) {
        int32_t candidateIndex = -1;
        auto candidateTableIndex = candidateMap.find(jet.candidateId());
        if (candidateTableIndex != candidateMap.end()) {
          candidateIndex = candidateTableIndex->second;
        }
        jetOutputTable(candidateIndex, jet.jetPt(), jet.jetPhi(), jet.jetEta(), jet.jetR(), jet.jetNConstituents());
        jetMap.insert(std::make_pair(jet.jetId(), jetOutputTable.lastIndex())); // this is filled with the standard jet Id and the recoil jet table position
      }
    }
  }

  template <typename T, typename U>
  void selectCandidates(T const& jets, U const& /*candidates*/, float jetPtMin, std::vector<bool>& candidateSelectionFlags)
  {
    for (const auto& jet : jets) {
      auto const& candidates = jet.template candidates_as<U>();
      if (jet.pt() < jetPtMin) {
        for (const auto& candidate : candidates) {
          candidateSelectionFlags[candidate.globalIndex()] = false;
        }
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
        for (const auto& candidate : candidates) {
          candidateSelectionFlags[candidate.globalIndex()] = false;
        }
        continue;
      }
      bool radiusSelected = false;
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          radiusSelected = true;
        }
      }
      if (!radiusSelected) {
        for (const auto& candidate : candidates) {
          candidateSelectionFlags[candidate.globalIndex()] = false;
        }
        continue;
      }
    }
  }

  template <typename T>
  void selectRecoilCandidates(T const& jets, float jetPtMin, std::vector<bool>& candidateSelectionFlags)
  {
    for (const auto& jet : jets) {
      auto candidateId = jet.candidateId();
      if (jet.pt() < jetPtMin) {
        candidateSelectionFlags[candidateId] = false;
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
        candidateSelectionFlags[candidateId] = false;
        continue;
      }
      bool radiusSelected = false;
      for (const auto& jetRadiiValue : recoilJetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          radiusSelected = true;
        }
      }
      if (!radiusSelected) {
        candidateSelectionFlags[candidateId] = false;
        continue;
      }
    }
  }

  template <bool isMCD, bool isMCP, typename T>
  void analyseCandidate(T const& candidate, std::map<int32_t, int32_t>& candidateMap, std::vector<bool> const& candidateSelectionFlags, std::vector<bool> const& candidateRecoilSelectionFlags)
  {
    if (!candidateSelectionFlags[candidate.globalIndex()] && !candidateRecoilSelectionFlags[candidate.globalIndex()]) {
      return;
    }
    int32_t candidateCollisionIndex = -1;
    if constexpr (isMCP) {
      auto hfMcCollisionIndex = candidateMcCollisionMapping.find(jetcandidateutilities::getMcCandidateCollisionId(candidate));
      if (hfMcCollisionIndex != candidateMcCollisionMapping.end()) {
        candidateCollisionIndex = hfMcCollisionIndex->second;
      }
      jetcandidateutilities::fillCandidateMcTable(candidate, candidateCollisionIndex, products.hfParticlesTable);
      candidateMap.insert(std::make_pair(candidate.globalIndex(), products.hfParticlesTable.lastIndex()));
    } else {
      auto hfCollisionIndex = candidateCollisionMapping.find(jetcandidateutilities::getCandidateCollisionId(candidate));
      if (hfCollisionIndex != candidateCollisionMapping.end()) {
        candidateCollisionIndex = hfCollisionIndex->second;
      }
      jetcandidateutilities::fillCandidateTable<isMCD>(candidate, candidateCollisionIndex, products.candidateTable, products.candidateParsTable, products.candidateParExtrasTable, products.candidateParsDaughterTable, products.candidateSelsTable, products.candidateMlsTable, products.candidateMlsDaughterTable, products.candidateMcsTable);
      candidateMap.insert(std::make_pair(candidate.globalIndex(), products.candidateTable.lastIndex()));
    }
  }

  template <typename CandidateTableType, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename R, typename S, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H, typename I, typename J, typename K>
  void analyseSubstructureMatched(T const& jets, U const& allSplittings, V const& allPairs, M const& D0SplittingsPerJet, N const DplusSplittingsPerJet, O const& DsSplittingsPerJet, P const DstarSplittingsPerJet, Q const& LcSplittingsPerJet, R const& B0SplittingsPerJet, S const& BplusSplittingsPerJet, A const& XicToXiPiPiSplittingsPerJet, B const& DielectronSplittingsPerJet, C const& D0PairsPerJet, D const DplusPairsPerJet, E const& DsPairsPerJet, F const& DstarPairsPerJet, G const& LcPairsPerJet, H const& B0PairsPerJet, I const& BplusPairsPerJet, J const& XicToXiPiPiPairsPerJet, K const& DielectronPairsPerJet, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, float jetPtMin)
  {
    for (const auto& jet : jets) {
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          auto splittings = jetcandidateutilities::slicedPerJet<CandidateTableType>(allSplittings, jet, D0SplittingsPerJet, DplusSplittingsPerJet, DsSplittingsPerJet, DstarSplittingsPerJet, LcSplittingsPerJet, B0SplittingsPerJet, BplusSplittingsPerJet, XicToXiPiPiSplittingsPerJet, DielectronSplittingsPerJet);
          fillSplittingMatchingVectors(splittings, jet.globalIndex(), splittingMatchesGeoVecVec, splittingMatchesPtVecVec, splittingMatchesHFVecVec);
          auto pairs = jetcandidateutilities::slicedPerJet<CandidateTableType>(allPairs, jet, D0PairsPerJet, DplusPairsPerJet, DsPairsPerJet, DstarPairsPerJet, LcPairsPerJet, B0PairsPerJet, BplusPairsPerJet, XicToXiPiPiPairsPerJet, DielectronPairsPerJet);
          fillPairMatchingVectors(pairs, jet.globalIndex(), pairMatchesVecVec);
        }
      }
    }
  }

  template <typename T, typename U>
  void analyseJetMatched(T const& jets, std::map<int32_t, int32_t>& jetMapping, std::map<int32_t, int32_t>& jetTagMapping, U& matchingOutputTable, float jetPtMin)
  {
    for (const auto& jet : jets) {
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
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

  template <typename T, typename U, typename V>
  void analyseRecoilJetMatched(T const& recoilJets, U const& /*jets*/, std::map<int32_t, int32_t>& jetMapping, std::map<int32_t, int32_t>& jetTagMapping, V& matchingOutputTable, float jetPtMin)
  {
    for (const auto& recoilJet : recoilJets) {
      if (recoilJet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(recoilJet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
        continue;
      }
      int storedJetIndex = -1;
      auto jetIndex = jetMapping.find(recoilJet.jetId());
      if (jetIndex != jetMapping.end()) {
        storedJetIndex = jetIndex->second;
      } else {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (recoilJet.r() == round(jetRadiiValue * 100.0f)) {
          auto const& jet = recoilJet.template jet_as<U>();
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
          matchingOutputTable(storedJetIndex, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }

  template <bool isMC, bool isMCPOnly, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseHFCollisions(T const& collisions, U const& mcCollisions, V const& hfCollisions, M const& hfMcCollisions, N const& candidates, O const& candidatesMCP, std::vector<bool> const& candidateSelectionFlags, std::vector<bool> const& candidateRecoilSelectionFlags, std::vector<bool> const& candidateMCPSelectionFlags, std::vector<bool> const& candidateMCPRecoilSelectionFlags)
  {
    collisionFlag.clear();
    collisionFlag.resize(collisions.size());
    std::fill(collisionFlag.begin(), collisionFlag.end(), false);

    mcCollisionFlag.clear();
    mcCollisionFlag.resize(mcCollisions.size());
    std::fill(mcCollisionFlag.begin(), mcCollisionFlag.end(), false);

    if constexpr (!isMCPOnly) {
      for (auto const& candidate : candidates) {

        if (candidateSelectionFlags[candidate.globalIndex()] || candidateRecoilSelectionFlags[candidate.globalIndex()]) {
          collisionFlag[candidate.collisionId()] = true;
          if constexpr (isMC) {
            auto mcCollisionId = candidate.template collision_as<T>().mcCollisionId();
            if (mcCollisionId >= 0) {
              mcCollisionFlag[mcCollisionId] = true;
            }
          }
        }
      }
    }
    if constexpr (isMC) {
      for (auto const& candidateMCP : candidatesMCP) {
        if (candidateMCPSelectionFlags[candidateMCP.globalIndex()] || candidateMCPRecoilSelectionFlags[candidateMCP.globalIndex()]) {
          mcCollisionFlag[candidateMCP.mcCollisionId()] = true;
          if constexpr (!isMCPOnly) {
            const auto collisionsPerMcCollision = collisions.sliceBy(preslices.CollisionsPerMcCollision, candidateMCP.mcCollisionId());
            for (auto collision : collisionsPerMcCollision) {
              collisionFlag[collision.globalIndex()] = true;
            }
          }
        }
      }
    }
    if constexpr (!isMCPOnly) {
      for (const auto& collision : collisions) {
        if (collisionFlag[collision.globalIndex()]) {
          const auto hfCollisionsPerCollision = hfCollisions.sliceBy(preslices.CandidateCollisionsPerCollision, collision.globalIndex());
          for (const auto& hfCollisionPerCollision : hfCollisionsPerCollision) { // should only ever be one
            auto hfCollisionTableIndex = candidateCollisionMapping.find(hfCollisionPerCollision.globalIndex());
            if (hfCollisionTableIndex != candidateCollisionMapping.end()) {
              continue;
            }
            jetcandidateutilities::fillCandidateCollisionTable(hfCollisionPerCollision, candidates, products.hfCollisionsTable);
            candidateCollisionMapping.insert(std::make_pair(hfCollisionPerCollision.globalIndex(), products.hfCollisionsTable.lastIndex()));
          }
        }
      }
    }
    if constexpr (isMC) {
      for (const auto& mcCollision : mcCollisions) {
        if (mcCollisionFlag[mcCollision.globalIndex()]) {
          const auto hfMcCollisionsPerMcCollision = hfMcCollisions.sliceBy(candidateMCCollisionSlicer<isMCPOnly>(preslices.CandidateMcCollisionsPerMcCollision, preslices.CandidateMcCollisionsPerMcCollisionMCPOnly), mcCollision.globalIndex());
          for (const auto& hfMcCollisionPerMcCollision : hfMcCollisionsPerMcCollision) { // should only ever be one
            auto hfMcCollisionTableIndex = candidateMcCollisionMapping.find(hfMcCollisionPerMcCollision.globalIndex());
            if (hfMcCollisionTableIndex != candidateMcCollisionMapping.end()) {
              continue;
            }
            jetcandidateutilities::fillCandidateMcCollisionTable(hfMcCollisionPerMcCollision, candidatesMCP, products.hfMcCollisionsTable);
            candidateMcCollisionMapping.insert(std::make_pair(hfMcCollisionPerMcCollision.globalIndex(), products.hfMcCollisionsTable.lastIndex()));
            if constexpr (!isMCPOnly && (jethfutilities::isHFTable<N>() || jethfutilities::isHFMcTable<O>())) { // the matching of mcCollision to Collision is only done for HF tables
              std::vector<int32_t> hfCollisionIDs;
              for (auto const& hfCollisionPerMcCollision : hfMcCollisionPerMcCollision.template hfCollBases_as<V>()) { // if added for others this line needs to be templated per type
                auto hfCollisionIndex = candidateCollisionMapping.find(hfCollisionPerMcCollision.globalIndex());
                if (hfCollisionIndex != candidateCollisionMapping.end()) {
                  hfCollisionIDs.push_back(hfCollisionIndex->second);
                }
              }
              products.hfMcCollisionsMatchingTable(hfCollisionIDs);
            }
          }
        }
      }
    }
  }

  void processClearMaps(CandidateTable const& candidates)
  {
    candidateMapping.clear();
    jetMappingData.clear();
    jetMappingDataSub.clear();
    jetMappingMCD.clear();
    recoilJetMappingData.clear();
    recoilJetMappingMCD.clear();
    candidateCollisionMapping.clear();
    if (doprocessOutputCandidatesData) {
      candidateSelectionFlagsData.clear();
      candidateSelectionFlagsData.resize(candidates.size(), true);
      candidateRecoilSelectionFlagsData.clear();
      candidateRecoilSelectionFlagsData.resize(candidates.size(), true);
    }
    if (doprocessOutputCandidatesMCD) {
      candidateSelectionFlagsMCD.clear();
      candidateSelectionFlagsMCD.resize(candidates.size(), true);
      candidateRecoilSelectionFlagsMCD.clear();
      candidateRecoilSelectionFlagsMCD.resize(candidates.size(), true);
    }
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processClearMaps, "process function that clears all the non-mcp maps in each dataframe", true);

  void processClearMapsMCP(o2::aod::JetMcCollisions const& mcCollisions, CandidateTableMCP const& candidates)
  {
    candidateMappingMCP.clear();
    jetMappingMCP.clear();
    recoilJetMappingMCP.clear();
    candidateMcCollisionMapping.clear();
    for (auto mcCollision : mcCollisions) {
      products.hfMcOnlyCollisionsTable(mcCollision.posZ(), mcCollision.accepted(), mcCollision.attempted(), mcCollision.xsectGen(), mcCollision.xsectErr(), mcCollision.weight());
    }
    candidateSelectionFlagsMCP.clear();
    candidateSelectionFlagsMCP.resize(candidates.size(), true);
    candidateRecoilSelectionFlagsMCP.clear();
    candidateRecoilSelectionFlagsMCP.resize(candidates.size(), true);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processClearMapsMCP, "process function that clears all the mcp maps in each dataframe", true);

  void processSelectCandidatesData(JetTableData const& jets, CandidateTable const& candidates)
  {
    selectCandidates(jets, candidates, configs.jetPtMinData, candidateSelectionFlagsData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processSelectCandidatesData, "select HF candidates for data", false);

  void processSelectRecoilCandidatesData(RecoilTableData const& jets)
  {
    selectRecoilCandidates(jets, configs.recoilJetPtMinData, candidateRecoilSelectionFlagsData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processSelectRecoilCandidatesData, "select HF candidates for recoil data", false);

  void processSelectCandidatesMCD(JetTableMCD const& jets, CandidateTableMCD const& candidates)
  {
    selectCandidates(jets, candidates, configs.jetPtMinMCD, candidateSelectionFlagsMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processSelectCandidatesMCD, "select HF candidates for mcd", false);

  void processSelectRecoilCandidatesMCD(RecoilTableMCD const& jets)
  {
    selectRecoilCandidates(jets, configs.recoilJetPtMinMCD, candidateRecoilSelectionFlagsMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processSelectRecoilCandidatesMCD, "select HF candidates for recoil mcd", false);

  void processSelectCandidatesMCP(JetTableMCP const& jets, CandidateTableMCP const& candidates)
  {
    selectCandidates(jets, candidates, configs.jetPtMinMCP, candidateSelectionFlagsMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processSelectCandidatesMCP, "select HF candidates for MCP", false);

  void processSelectRecoilCandidatesMCP(RecoilTableMCP const& jets)
  {
    selectRecoilCandidates(jets, configs.recoilJetPtMinMCP, candidateRecoilSelectionFlagsMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processSelectRecoilCandidatesMCP, "select HF candidates for recoil MCP", false);

  void processOutputCollisionsData(o2::aod::JetCollisions const& collisions,
                                   CandidateCollisionTable const& canidateCollisions,
                                   CandidateTable const& candidates)
  {
    analyseHFCollisions<false, false>(collisions, collisions, canidateCollisions, canidateCollisions, candidates, candidates, candidateSelectionFlagsData, candidateRecoilSelectionFlagsData, candidateSelectionFlagsData, candidateRecoilSelectionFlagsData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsData, "hf collision output data", false);

  void processOutputCollisionsDataSub(o2::aod::JetCollisions const& collisions,
                                      CandidateCollisionTable const& canidateCollisions,
                                      CandidateTable const& candidates)
  {
    analyseHFCollisions<false, false>(collisions, collisions, canidateCollisions, canidateCollisions, candidates, candidates, candidateSelectionFlagsData, candidateRecoilSelectionFlagsData, candidateSelectionFlagsData, candidateRecoilSelectionFlagsData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsDataSub, "hf collision output data eventwise constituent subtracted", false);

  void processOutputCollisionsMc(o2::soa::Join<o2::aod::JetCollisions, o2::aod::JMcCollisionLbs> const& collisions,
                                 o2::aod::JetMcCollisions const& mcCollisions,
                                 CandidateCollisionTable const& canidateCollisions,
                                 CandidateMcCollisionTable const& canidateMcCollisions,
                                 CandidateTableMCD const& candidatesMCD,
                                 CandidateTableMCP const& candidatesMCP)
  {
    analyseHFCollisions<true, false>(collisions, mcCollisions, canidateCollisions, canidateMcCollisions, candidatesMCD, candidatesMCP, candidateSelectionFlagsMCD, candidateRecoilSelectionFlagsMCD, candidateSelectionFlagsMCP, candidateRecoilSelectionFlagsMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsMc, "hf collision output MC", false);

  void processOutputCollisionsMCPOnly(o2::aod::JetMcCollisions const& mcCollisions,
                                      CandidateMcOnlyCollisionTable const& canidateMcCollisions,
                                      CandidateTableMCP const& candidatesMCP)
  {
    analyseHFCollisions<true, true>(mcCollisions, mcCollisions, canidateMcCollisions, canidateMcCollisions, candidatesMCP, candidatesMCP, candidateSelectionFlagsMCP, candidateRecoilSelectionFlagsMCP, candidateSelectionFlagsMCP, candidateRecoilSelectionFlagsMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsMCPOnly, "hf collision output MCP only", false);

  void processOutputCandidatesData(typename CandidateTable::iterator const& candidate)
  {
    analyseCandidate<false, false>(candidate, candidateMapping, candidateSelectionFlagsData, candidateRecoilSelectionFlagsData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesData, "hf candidate output data", false);

  void processOutputCandidatesMCD(typename CandidateTableMCD::iterator const& candidate)
  {

    analyseCandidate<true, false>(candidate, candidateMapping, candidateSelectionFlagsMCD, candidateRecoilSelectionFlagsMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesMCD, "hf candidate output MCD", false);

  void processOutputCandidatesMCP(typename CandidateTableMCP::iterator const& candidate)
  {
    analyseCandidate<false, true>(candidate, candidateMappingMCP, candidateSelectionFlagsMCP, candidateRecoilSelectionFlagsMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesMCP, "hf candidate output MCP", false);

  void processOutputSubstructureMatchingData(JetMatchedTableData const& jets,
                                             JetTableDataSub const& jetsSub,
                                             CandidateTable const&,
                                             SplittingTableData const& splittingsData,
                                             SplittingTableDataSub const& splittingsDataSub,
                                             PairTableData const& pairsData,
                                             PairTableDataSub const& pairsDataSub)
  {
    splittingMatchesGeoVecVecData.assign(jets.size(), {});
    splittingMatchesPtVecVecData.assign(jets.size(), {});
    splittingMatchesHFVecVecData.assign(jets.size(), {});
    pairMatchesVecVecData.assign(jets.size(), {});
    analyseSubstructureMatched<CandidateTable>(jets, splittingsData, pairsData, preslices.D0SplittingsPerJetData, preslices.DplusSplittingsPerJetData, preslices.DsSplittingsPerJetData, preslices.DstarSplittingsPerJetData, preslices.LcSplittingsPerJetData, preslices.B0SplittingsPerJetData, preslices.BplusSplittingsPerJetData, preslices.XicToXiPiPiSplittingsPerJetData, preslices.DielectronSplittingsPerJetData, preslices.D0PairsPerJetData, preslices.DplusPairsPerJetData, preslices.DsPairsPerJetData, preslices.DstarPairsPerJetData, preslices.LcPairsPerJetData, preslices.B0PairsPerJetData, preslices.BplusPairsPerJetData, preslices.XicToXiPiPiPairsPerJetData, preslices.DielectronPairsPerJetData, splittingMatchesGeoVecVecData, splittingMatchesPtVecVecData, splittingMatchesHFVecVecData, pairMatchesVecVecData, configs.jetPtMinData);
    splittingMatchesGeoVecVecDataSub.assign(jetsSub.size(), {});
    splittingMatchesPtVecVecDataSub.assign(jetsSub.size(), {});
    splittingMatchesHFVecVecDataSub.assign(jetsSub.size(), {});
    pairMatchesVecVecDataSub.assign(jetsSub.size(), {});
    analyseSubstructureMatched<CandidateTable>(jetsSub, splittingsDataSub, pairsDataSub, preslices.D0SplittingsPerJetDataSub, preslices.DplusSplittingsPerJetDataSub, preslices.DsSplittingsPerJetDataSub, preslices.DstarSplittingsPerJetDataSub, preslices.LcSplittingsPerJetDataSub, preslices.B0SplittingsPerJetDataSub, preslices.BplusSplittingsPerJetDataSub, preslices.XicToXiPiPiSplittingsPerJetDataSub, preslices.DielectronSplittingsPerJetDataSub, preslices.D0PairsPerJetDataSub, preslices.DplusPairsPerJetDataSub, preslices.DsPairsPerJetDataSub, preslices.DstarPairsPerJetDataSub, preslices.LcPairsPerJetDataSub, preslices.B0PairsPerJetDataSub, preslices.BplusPairsPerJetDataSub, preslices.XicToXiPiPiPairsPerJetDataSub, preslices.DielectronPairsPerJetDataSub, splittingMatchesGeoVecVecDataSub, splittingMatchesPtVecVecDataSub, splittingMatchesHFVecVecDataSub, pairMatchesVecVecDataSub, configs.jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputSubstructureMatchingData, "jet substructure matching output Data", false);

  void processOutputJetsData(o2::aod::JetCollision const& collision,
                             JetTableData const& jets,
                             o2::soa::Join<CandidateTable, CandidateRhosTable> const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, products.collisionOutputTableData, products.jetOutputTableData, products.jetSubstructureOutputTableData, splittingMatchesGeoVecVecData, splittingMatchesPtVecVecData, splittingMatchesHFVecVecData, pairMatchesVecVecData, jetMappingData, candidateMapping, configs.jetPtMinData, 1.0);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsData, "hf jet substructure output Data", false);

  void processOutputJetsDataSub(o2::aod::JetCollision const& collision,
                                JetTableDataSub const& jets,
                                o2::soa::Join<CandidateTable, CandidateRhosTable> const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, products.collisionOutputTableDataSub, products.jetOutputTableDataSub, products.jetSubstructureOutputTableDataSub, splittingMatchesGeoVecVecDataSub, splittingMatchesPtVecVecDataSub, splittingMatchesHFVecVecDataSub, pairMatchesVecVecDataSub, jetMappingDataSub, candidateMapping, configs.jetPtMinDataSub, 1.0);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsDataSub, "hf jet substructure output event-wise subtracted Data", false);

  void processOutputJetMatchingData(JetMatchedTableData const& jets,
                                    JetTableDataSub const& jetsSub)
  {
    analyseJetMatched(jets, jetMappingData, jetMappingDataSub, products.jetMatchingOutputTableData, configs.jetPtMinData);
    analyseJetMatched(jetsSub, jetMappingDataSub, jetMappingData, products.jetMatchingOutputTableDataSub, configs.jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetMatchingData, "jet matching output Data", false);

  void processOutputSubstructureMatchingMC(JetTableMCD const& jetsMCD,
                                           JetTableMatchedMCP const& jetsMCP,
                                           CandidateTableMCD const&,
                                           CandidateTableMCP const&,
                                           SplittingTableMCD const& splittingsMCD,
                                           SplittingTableMCP const& splittingsMCP,
                                           PairTableMCD const& pairsMCD,
                                           PairTableMCP const& pairsMCP)
  {
    splittingMatchesGeoVecVecMCD.assign(jetsMCD.size(), {});
    splittingMatchesPtVecVecMCD.assign(jetsMCD.size(), {});
    splittingMatchesHFVecVecMCD.assign(jetsMCD.size(), {});
    pairMatchesVecVecMCD.assign(jetsMCD.size(), {});
    analyseSubstructureMatched<CandidateTableMCD>(jetsMCD, splittingsMCD, pairsMCD, preslices.D0SplittingsPerJetMCD, preslices.DplusSplittingsPerJetMCD, preslices.DsSplittingsPerJetMCD, preslices.DstarSplittingsPerJetMCD, preslices.LcSplittingsPerJetMCD, preslices.B0SplittingsPerJetMCD, preslices.BplusSplittingsPerJetMCD, preslices.XicToXiPiPiSplittingsPerJetMCD, preslices.DielectronSplittingsPerJetMCD, preslices.D0PairsPerJetMCD, preslices.DplusPairsPerJetMCD, preslices.DsPairsPerJetMCD, preslices.DstarPairsPerJetMCD, preslices.LcPairsPerJetMCD, preslices.B0PairsPerJetMCD, preslices.BplusPairsPerJetMCD, preslices.XicToXiPiPiPairsPerJetMCD, preslices.DielectronPairsPerJetMCD, splittingMatchesGeoVecVecMCD, splittingMatchesPtVecVecMCD, splittingMatchesHFVecVecMCD, pairMatchesVecVecMCD, configs.jetPtMinMCD);
    splittingMatchesGeoVecVecMCP.assign(jetsMCP.size(), {});
    splittingMatchesPtVecVecMCP.assign(jetsMCP.size(), {});
    splittingMatchesHFVecVecMCP.assign(jetsMCP.size(), {});
    pairMatchesVecVecMCP.assign(jetsMCP.size(), {});
    analyseSubstructureMatched<CandidateTableMCP>(jetsMCP, splittingsMCP, pairsMCP, preslices.D0SplittingsPerJetMCP, preslices.DplusSplittingsPerJetMCP, preslices.DsSplittingsPerJetMCP, preslices.DstarSplittingsPerJetMCP, preslices.LcSplittingsPerJetMCP, preslices.B0SplittingsPerJetMCP, preslices.BplusSplittingsPerJetMCP, preslices.XicToXiPiPiSplittingsPerJetMCP, preslices.DielectronSplittingsPerJetMCP, preslices.D0PairsPerJetMCP, preslices.DplusPairsPerJetMCP, preslices.DsPairsPerJetMCP, preslices.DstarPairsPerJetMCP, preslices.LcPairsPerJetMCP, preslices.B0PairsPerJetMCP, preslices.BplusPairsPerJetMCP, preslices.XicToXiPiPiPairsPerJetMCP, preslices.DielectronPairsPerJetMCP, splittingMatchesGeoVecVecMCP, splittingMatchesPtVecVecMCP, splittingMatchesHFVecVecMCP, pairMatchesVecVecMCP, configs.jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputSubstructureMatchingMC, "jet substructure matching output MC", false);

  void processOutputJetsMCD(o2::aod::JetCollisionMCD const& collision,
                            JetTableMCD const& jets,
                            o2::soa::Join<CandidateTableMCD, CandidateRhosTable> const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, products.collisionOutputTableMCD, products.jetOutputTableMCD, products.jetSubstructureOutputTableMCD, splittingMatchesGeoVecVecMCD, splittingMatchesPtVecVecMCD, splittingMatchesHFVecVecMCD, pairMatchesVecVecMCD, jetMappingMCD, candidateMapping, configs.jetPtMinMCD, collision.weight());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsMCD, "hf jet substructure output MCD", false);

  void processOutputJetsMCP(o2::aod::JetMcCollision const& collision,
                            JetTableMCP const& jets,
                            o2::soa::Join<CandidateTableMCP, CandidateMCRhosTable> const& candidates)
  {
    analyseCharged<true>(collision, jets, candidates, products.collisionOutputTableMCP, products.jetOutputTableMCP, products.jetSubstructureOutputTableMCP, splittingMatchesGeoVecVecMCP, splittingMatchesPtVecVecMCP, splittingMatchesHFVecVecMCP, pairMatchesVecVecMCP, jetMappingMCP, candidateMappingMCP, configs.jetPtMinMCP, collision.weight());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsMCP, "hf jet substructure output MCP", false);

  void processOutputJetMatchingMC(JetTableMCD const& jetsMCD,
                                  JetTableMatchedMCP const& jetsMCP)
  {
    analyseJetMatched(jetsMCD, jetMappingMCD, jetMappingMCP, products.jetMatchingOutputTableMCD, configs.jetPtMinMCD);
    analyseJetMatched(jetsMCP, jetMappingMCP, jetMappingMCD, products.jetMatchingOutputTableMCP, configs.jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetMatchingMC, "jet matching output MC", false);

  void processOutputRecoilJetsData(typename RecoilTableData::iterator const& recoilJet)
  {
    analyseRecoilCharged(recoilJet, products.jetRecoilOutputTableData, recoilJetMappingData, candidateMapping, configs.recoilJetPtMinData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputRecoilJetsData, "hf recoil jet output Data", false);

  void processOutputRecoilJetsMCD(typename RecoilTableMCD::iterator const& recoilJet)
  {
    analyseRecoilCharged(recoilJet, products.jetRecoilOutputTableMCD, recoilJetMappingMCD, candidateMapping, configs.recoilJetPtMinMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputRecoilJetsMCD, "hf recoil jet output mcd", false);

  void processOutputRecoilJetsMCP(typename RecoilTableMCP::iterator const& recoilJet)
  {
    analyseRecoilCharged(recoilJet, products.jetRecoilOutputTableMCP, recoilJetMappingMCP, candidateMappingMCP, configs.recoilJetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputRecoilJetsMCP, "hf recoil jet output mcp", false);

  void processOutputRecoilJetMatchingMC(RecoilTableMCD const& recoilJetsMCD,
                                        RecoilTableMCP const& recoilJetsMCP,
                                        o2::soa::Join<o2::aod::ChargedMCDetectorLevelJets, o2::aod::ChargedMCDetectorLevelJetConstituents, o2::aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& jetsMCD,
                                        o2::soa::Join<o2::aod::ChargedMCParticleLevelJets, o2::aod::ChargedMCParticleLevelJetConstituents, o2::aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& jetsMCP)
  {
    analyseRecoilJetMatched(recoilJetsMCD, jetsMCD, recoilJetMappingMCD, recoilJetMappingMCP, products.jetRecoilMatchingOutputTableMCD, configs.recoilJetPtMinMCD);
    analyseRecoilJetMatched(recoilJetsMCP, jetsMCP, recoilJetMappingMCP, recoilJetMappingMCD, products.jetRecoilMatchingOutputTableMCP, configs.recoilJetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputRecoilJetMatchingMC, "recoil jet matching output MC", false);
};

#endif // PWGJE_TASKS_JETSUBSTRUCTUREHFOUTPUT_H_
