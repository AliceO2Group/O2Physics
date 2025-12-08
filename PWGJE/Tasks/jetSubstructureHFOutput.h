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

#include "Framework/ASoA.h"
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!

template <typename CandidateCollisionTable, typename CandidateMcCollisionTable, typename CandidateMcOnlyCollisionTable, typename CandidateTable, typename CandidateTableMCD, typename CandidateTableMCP, typename CandidateRhosTable, typename CandidateMCRhosTable, typename TracksSub, typename JetTableData, typename JetMatchedTableData, typename SplittingTableData, typename PairTableData, typename OutputCollisionTableData, typename OutputTableData, typename SubstructureOutputTableData, typename MatchingOutputTableData, typename JetTableMCD, typename SplittingTableMCD, typename PairTableMCD, typename OutputCollisionTableMCD, typename OutputTableMCD, typename SubstructureOutputTableMCD, typename MatchingOutputTableMCD, typename JetTableMCP, typename JetTableMatchedMCP, typename SplittingTableMCP, typename PairTableMCP, typename OutputCollisionTableMCP, typename CandidateMcOnlyCollisionOutputTable, typename OutputTableMCP, typename SubstructureOutputTableMCP, typename MatchingOutputTableMCP, typename JetTableDataSub, typename SplittingTableDataSub, typename PairTableDataSub, typename OutputCollisionTableDataSub, typename OutputTableDataSub, typename SubstructureOutputTableDataSub, typename MatchingOutputTableDataSub, typename CandidateCollisionOutputTable, typename CandidateOutputTable, typename CandidateParOutputTable, typename CandidateParExtraOutputTable, typename CandidateParDaughterOutputTable, typename CandidateSelOutputTable, typename CandidateMlOutputTable, typename CandidateMlDaughterOutputTable, typename CandidateMCDOutputTable, typename CandidateMcCollisionOutputTable, typename CandidateMcCollisionMatchingOutputTable, typename CandidateMCPOutputTable>
struct JetSubstructureHFOutputTask {

  struct : ProducesGroup {
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
    Produces<CandidateMcOnlyCollisionOutputTable> hfMcOnlyCollisionsTable;
    Produces<OutputTableMCP> jetOutputTableMCP;
    Produces<SubstructureOutputTableMCP> jetSubstructureOutputTableMCP;
    Produces<MatchingOutputTableMCP> jetMatchingOutputTableMCP;
    Produces<CandidateCollisionOutputTable> hfCollisionsTable;
    Produces<CandidateOutputTable> candidateTable;
    Produces<CandidateParOutputTable> candidateParsTable;
    Produces<CandidateParExtraOutputTable> candidateParExtrasTable;
    Produces<CandidateParDaughterOutputTable> candidateParsDaughterTable;
    Produces<CandidateSelOutputTable> candidateSelsTable;
    Produces<CandidateMlOutputTable> candidateMlsTable;
    Produces<CandidateMlDaughterOutputTable> candidateMlsDaughterTable;
    Produces<CandidateMCDOutputTable> candidateMcsTable;
    Produces<CandidateMcCollisionOutputTable> hfMcCollisionsTable;
    Produces<CandidateMcCollisionMatchingOutputTable> hfMcCollisionsMatchingTable;
    Produces<CandidateMCPOutputTable> hfParticlesTable;
  } products;

  struct : ConfigurableGroup {
    Configurable<float> jetPtMinData{"jetPtMinData", 0.0, "minimum jet pT cut for data jets"};
    Configurable<float> jetPtMinDataSub{"jetPtMinDataSub", 0.0, "minimum jet pT cut for eventwise constituent subtracted data jets"};
    Configurable<float> jetPtMinMCD{"jetPtMinMCD", 0.0, "minimum jet pT cut for mcd jets"};
    Configurable<float> jetPtMinMCP{"jetPtMinMCP", 0.0, "minimum jet pT cut for mcp jets"};
    Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
    Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
    Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};
    Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track pseudorapidity"};
    Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track pseudorapidity"};
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

  std::vector<bool> collisionFlag;
  std::vector<bool> mcCollisionFlag;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)configs.jetRadii;
  }

  struct : PresliceGroup {
    PresliceUnsorted<soa::Join<aod::JetCollisions, aod::JMcCollisionLbs>> CollisionsPerMcCollision = aod::jmccollisionlb::mcCollisionId;
    PresliceOptional<CandidateCollisionTable> CandidateCollisionsPerCollision = aod::jcandidateindices::collisionId;
    PresliceOptional<CandidateMcCollisionTable> CandidateMcCollisionsPerMcCollision = aod::jcandidateindices::mcCollisionId;
    PresliceOptional<CandidateMcOnlyCollisionTable> CandidateMcCollisionsPerMcCollisionMCPOnly = aod::jcandidateindices::mcCollisionId;

    PresliceOptional<soa::Join<aod::D0ChargedSPs, aod::D0ChargedSPsMatchedToD0ChargedEventWiseSubtractedSPs>> D0SplittingsPerJetData = aod::d0chargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::D0ChargedEventWiseSubtractedSPs, aod::D0ChargedEventWiseSubtractedSPsMatchedToD0ChargedSPs>> D0SplittingsPerJetDataSub = aod::d0chargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::D0ChargedMCDetectorLevelSPs, aod::D0ChargedMCDetectorLevelSPsMatchedToD0ChargedMCParticleLevelSPs>> D0SplittingsPerJetMCD = aod::d0chargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::D0ChargedMCParticleLevelSPs, aod::D0ChargedMCParticleLevelSPsMatchedToD0ChargedMCDetectorLevelSPs>> D0SplittingsPerJetMCP = aod::d0chargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::D0ChargedPRs, aod::D0ChargedPRsMatchedToD0ChargedEventWiseSubtractedPRs>> D0PairsPerJetData = aod::d0chargedpair::jetId;
    PresliceOptional<soa::Join<aod::D0ChargedEventWiseSubtractedPRs, aod::D0ChargedEventWiseSubtractedPRsMatchedToD0ChargedPRs>> D0PairsPerJetDataSub = aod::d0chargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::D0ChargedMCDetectorLevelPRs, aod::D0ChargedMCDetectorLevelPRsMatchedToD0ChargedMCParticleLevelPRs>> D0PairsPerJetMCD = aod::d0chargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::D0ChargedMCParticleLevelPRs, aod::D0ChargedMCParticleLevelPRsMatchedToD0ChargedMCDetectorLevelPRs>> D0PairsPerJetMCP = aod::d0chargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::DplusChargedSPs, aod::DplusChargedSPsMatchedToDplusChargedEventWiseSubtractedSPs>> DplusSplittingsPerJetData = aod::dpluschargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DplusChargedEventWiseSubtractedSPs, aod::DplusChargedEventWiseSubtractedSPsMatchedToDplusChargedSPs>> DplusSplittingsPerJetDataSub = aod::dpluschargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DplusChargedMCDetectorLevelSPs, aod::DplusChargedMCDetectorLevelSPsMatchedToDplusChargedMCParticleLevelSPs>> DplusSplittingsPerJetMCD = aod::dpluschargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::DplusChargedMCParticleLevelSPs, aod::DplusChargedMCParticleLevelSPsMatchedToDplusChargedMCDetectorLevelSPs>> DplusSplittingsPerJetMCP = aod::dpluschargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::DplusChargedPRs, aod::DplusChargedPRsMatchedToDplusChargedEventWiseSubtractedPRs>> DplusPairsPerJetData = aod::dpluschargedpair::jetId;
    PresliceOptional<soa::Join<aod::DplusChargedEventWiseSubtractedPRs, aod::DplusChargedEventWiseSubtractedPRsMatchedToDplusChargedPRs>> DplusPairsPerJetDataSub = aod::dpluschargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::DplusChargedMCDetectorLevelPRs, aod::DplusChargedMCDetectorLevelPRsMatchedToDplusChargedMCParticleLevelPRs>> DplusPairsPerJetMCD = aod::dpluschargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::DplusChargedMCParticleLevelPRs, aod::DplusChargedMCParticleLevelPRsMatchedToDplusChargedMCDetectorLevelPRs>> DplusPairsPerJetMCP = aod::dpluschargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::DsChargedSPs, aod::DsChargedSPsMatchedToDsChargedEventWiseSubtractedSPs>> DsSplittingsPerJetData = aod::dschargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DsChargedEventWiseSubtractedSPs, aod::DsChargedEventWiseSubtractedSPsMatchedToDsChargedSPs>> DsSplittingsPerJetDataSub = aod::dschargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DsChargedMCDetectorLevelSPs, aod::DsChargedMCDetectorLevelSPsMatchedToDsChargedMCParticleLevelSPs>> DsSplittingsPerJetMCD = aod::dschargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::DsChargedMCParticleLevelSPs, aod::DsChargedMCParticleLevelSPsMatchedToDsChargedMCDetectorLevelSPs>> DsSplittingsPerJetMCP = aod::dschargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::DsChargedPRs, aod::DsChargedPRsMatchedToDsChargedEventWiseSubtractedPRs>> DsPairsPerJetData = aod::dschargedpair::jetId;
    PresliceOptional<soa::Join<aod::DsChargedEventWiseSubtractedPRs, aod::DsChargedEventWiseSubtractedPRsMatchedToDsChargedPRs>> DsPairsPerJetDataSub = aod::dschargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::DsChargedMCDetectorLevelPRs, aod::DsChargedMCDetectorLevelPRsMatchedToDsChargedMCParticleLevelPRs>> DsPairsPerJetMCD = aod::dschargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::DsChargedMCParticleLevelPRs, aod::DsChargedMCParticleLevelPRsMatchedToDsChargedMCDetectorLevelPRs>> DsPairsPerJetMCP = aod::dschargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::DstarChargedSPs, aod::DstarChargedSPsMatchedToDstarChargedEventWiseSubtractedSPs>> DstarSplittingsPerJetData = aod::dstarchargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DstarChargedEventWiseSubtractedSPs, aod::DstarChargedEventWiseSubtractedSPsMatchedToDstarChargedSPs>> DstarSplittingsPerJetDataSub = aod::dstarchargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DstarChargedMCDetectorLevelSPs, aod::DstarChargedMCDetectorLevelSPsMatchedToDstarChargedMCParticleLevelSPs>> DstarSplittingsPerJetMCD = aod::dstarchargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::DstarChargedMCParticleLevelSPs, aod::DstarChargedMCParticleLevelSPsMatchedToDstarChargedMCDetectorLevelSPs>> DstarSplittingsPerJetMCP = aod::dstarchargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::DstarChargedPRs, aod::DstarChargedPRsMatchedToDstarChargedEventWiseSubtractedPRs>> DstarPairsPerJetData = aod::dstarchargedpair::jetId;
    PresliceOptional<soa::Join<aod::DstarChargedEventWiseSubtractedPRs, aod::DstarChargedEventWiseSubtractedPRsMatchedToDstarChargedPRs>> DstarPairsPerJetDataSub = aod::dstarchargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::DstarChargedMCDetectorLevelPRs, aod::DstarChargedMCDetectorLevelPRsMatchedToDstarChargedMCParticleLevelPRs>> DstarPairsPerJetMCD = aod::dstarchargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::DstarChargedMCParticleLevelPRs, aod::DstarChargedMCParticleLevelPRsMatchedToDstarChargedMCDetectorLevelPRs>> DstarPairsPerJetMCP = aod::dstarchargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::LcChargedSPs, aod::LcChargedSPsMatchedToLcChargedEventWiseSubtractedSPs>> LcSplittingsPerJetData = aod::lcchargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::LcChargedEventWiseSubtractedSPs, aod::LcChargedEventWiseSubtractedSPsMatchedToLcChargedSPs>> LcSplittingsPerJetDataSub = aod::lcchargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::LcChargedMCDetectorLevelSPs, aod::LcChargedMCDetectorLevelSPsMatchedToLcChargedMCParticleLevelSPs>> LcSplittingsPerJetMCD = aod::lcchargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::LcChargedMCParticleLevelSPs, aod::LcChargedMCParticleLevelSPsMatchedToLcChargedMCDetectorLevelSPs>> LcSplittingsPerJetMCP = aod::lcchargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::LcChargedPRs, aod::LcChargedPRsMatchedToLcChargedEventWiseSubtractedPRs>> LcPairsPerJetData = aod::lcchargedpair::jetId;
    PresliceOptional<soa::Join<aod::LcChargedEventWiseSubtractedPRs, aod::LcChargedEventWiseSubtractedPRsMatchedToLcChargedPRs>> LcPairsPerJetDataSub = aod::lcchargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::LcChargedMCDetectorLevelPRs, aod::LcChargedMCDetectorLevelPRsMatchedToLcChargedMCParticleLevelPRs>> LcPairsPerJetMCD = aod::lcchargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::LcChargedMCParticleLevelPRs, aod::LcChargedMCParticleLevelPRsMatchedToLcChargedMCDetectorLevelPRs>> LcPairsPerJetMCP = aod::lcchargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::B0ChargedSPs, aod::B0ChargedSPsMatchedToB0ChargedEventWiseSubtractedSPs>> B0SplittingsPerJetData = aod::b0chargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::B0ChargedEventWiseSubtractedSPs, aod::B0ChargedEventWiseSubtractedSPsMatchedToB0ChargedSPs>> B0SplittingsPerJetDataSub = aod::b0chargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::B0ChargedMCDetectorLevelSPs, aod::B0ChargedMCDetectorLevelSPsMatchedToB0ChargedMCParticleLevelSPs>> B0SplittingsPerJetMCD = aod::b0chargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::B0ChargedMCParticleLevelSPs, aod::B0ChargedMCParticleLevelSPsMatchedToB0ChargedMCDetectorLevelSPs>> B0SplittingsPerJetMCP = aod::b0chargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::B0ChargedPRs, aod::B0ChargedPRsMatchedToB0ChargedEventWiseSubtractedPRs>> B0PairsPerJetData = aod::b0chargedpair::jetId;
    PresliceOptional<soa::Join<aod::B0ChargedEventWiseSubtractedPRs, aod::B0ChargedEventWiseSubtractedPRsMatchedToB0ChargedPRs>> B0PairsPerJetDataSub = aod::b0chargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::B0ChargedMCDetectorLevelPRs, aod::B0ChargedMCDetectorLevelPRsMatchedToB0ChargedMCParticleLevelPRs>> B0PairsPerJetMCD = aod::b0chargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::B0ChargedMCParticleLevelPRs, aod::B0ChargedMCParticleLevelPRsMatchedToB0ChargedMCDetectorLevelPRs>> B0PairsPerJetMCP = aod::b0chargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::BplusChargedSPs, aod::BplusChargedSPsMatchedToBplusChargedEventWiseSubtractedSPs>> BplusSplittingsPerJetData = aod::bpluschargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::BplusChargedEventWiseSubtractedSPs, aod::BplusChargedEventWiseSubtractedSPsMatchedToBplusChargedSPs>> BplusSplittingsPerJetDataSub = aod::bpluschargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::BplusChargedMCDetectorLevelSPs, aod::BplusChargedMCDetectorLevelSPsMatchedToBplusChargedMCParticleLevelSPs>> BplusSplittingsPerJetMCD = aod::bpluschargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::BplusChargedMCParticleLevelSPs, aod::BplusChargedMCParticleLevelSPsMatchedToBplusChargedMCDetectorLevelSPs>> BplusSplittingsPerJetMCP = aod::bpluschargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::BplusChargedPRs, aod::BplusChargedPRsMatchedToBplusChargedEventWiseSubtractedPRs>> BplusPairsPerJetData = aod::bpluschargedpair::jetId;
    PresliceOptional<soa::Join<aod::BplusChargedEventWiseSubtractedPRs, aod::BplusChargedEventWiseSubtractedPRsMatchedToBplusChargedPRs>> BplusPairsPerJetDataSub = aod::bpluschargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::BplusChargedMCDetectorLevelPRs, aod::BplusChargedMCDetectorLevelPRsMatchedToBplusChargedMCParticleLevelPRs>> BplusPairsPerJetMCD = aod::bpluschargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::BplusChargedMCParticleLevelPRs, aod::BplusChargedMCParticleLevelPRsMatchedToBplusChargedMCDetectorLevelPRs>> BplusPairsPerJetMCP = aod::bpluschargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedSPs, aod::XicToXiPiPiChargedSPsMatchedToXicToXiPiPiChargedEventWiseSubtractedSPs>> XicToXiPiPiSplittingsPerJetData = aod::xictoxipipichargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedEventWiseSubtractedSPs, aod::XicToXiPiPiChargedEventWiseSubtractedSPsMatchedToXicToXiPiPiChargedSPs>> XicToXiPiPiSplittingsPerJetDataSub = aod::xictoxipipichargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedMCDetectorLevelSPs, aod::XicToXiPiPiChargedMCDetectorLevelSPsMatchedToXicToXiPiPiChargedMCParticleLevelSPs>> XicToXiPiPiSplittingsPerJetMCD = aod::xictoxipipichargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedMCParticleLevelSPs, aod::XicToXiPiPiChargedMCParticleLevelSPsMatchedToXicToXiPiPiChargedMCDetectorLevelSPs>> XicToXiPiPiSplittingsPerJetMCP = aod::xictoxipipichargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedPRs, aod::XicToXiPiPiChargedPRsMatchedToXicToXiPiPiChargedEventWiseSubtractedPRs>> XicToXiPiPiPairsPerJetData = aod::xictoxipipichargedpair::jetId;
    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedEventWiseSubtractedPRs, aod::XicToXiPiPiChargedEventWiseSubtractedPRsMatchedToXicToXiPiPiChargedPRs>> XicToXiPiPiPairsPerJetDataSub = aod::xictoxipipichargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedMCDetectorLevelPRs, aod::XicToXiPiPiChargedMCDetectorLevelPRsMatchedToXicToXiPiPiChargedMCParticleLevelPRs>> XicToXiPiPiPairsPerJetMCD = aod::xictoxipipichargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::XicToXiPiPiChargedMCParticleLevelPRs, aod::XicToXiPiPiChargedMCParticleLevelPRsMatchedToXicToXiPiPiChargedMCDetectorLevelPRs>> XicToXiPiPiPairsPerJetMCP = aod::xictoxipipichargedmcparticlelevelpair::jetId;

    PresliceOptional<soa::Join<aod::DielectronChargedSPs, aod::DielectronChargedSPsMatchedToDielectronChargedEventWiseSubtractedSPs>> DielectronSplittingsPerJetData = aod::dielectronchargedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DielectronChargedEventWiseSubtractedSPs, aod::DielectronChargedEventWiseSubtractedSPsMatchedToDielectronChargedSPs>> DielectronSplittingsPerJetDataSub = aod::dielectronchargedeventwisesubtractedsplitting::jetId;
    PresliceOptional<soa::Join<aod::DielectronChargedMCDetectorLevelSPs, aod::DielectronChargedMCDetectorLevelSPsMatchedToDielectronChargedMCParticleLevelSPs>> DielectronSplittingsPerJetMCD = aod::dielectronchargedmcdetectorlevelsplitting::jetId;
    PresliceOptional<soa::Join<aod::DielectronChargedMCParticleLevelSPs, aod::DielectronChargedMCParticleLevelSPsMatchedToDielectronChargedMCDetectorLevelSPs>> DielectronSplittingsPerJetMCP = aod::dielectronchargedmcparticlelevelsplitting::jetId;

    PresliceOptional<soa::Join<aod::DielectronChargedPRs, aod::DielectronChargedPRsMatchedToDielectronChargedEventWiseSubtractedPRs>> DielectronPairsPerJetData = aod::dielectronchargedpair::jetId;
    PresliceOptional<soa::Join<aod::DielectronChargedEventWiseSubtractedPRs, aod::DielectronChargedEventWiseSubtractedPRsMatchedToDielectronChargedPRs>> DielectronPairsPerJetDataSub = aod::dielectronchargedeventwisesubtractedpair::jetId;
    PresliceOptional<soa::Join<aod::DielectronChargedMCDetectorLevelPRs, aod::DielectronChargedMCDetectorLevelPRsMatchedToDielectronChargedMCParticleLevelPRs>> DielectronPairsPerJetMCD = aod::dielectronchargedmcdetectorlevelpair::jetId;
    PresliceOptional<soa::Join<aod::DielectronChargedMCParticleLevelPRs, aod::DielectronChargedMCParticleLevelPRsMatchedToDielectronChargedMCDetectorLevelPRs>> DielectronPairsPerJetMCP = aod::dielectronchargedmcparticlelevelpair::jetId;

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

  template <typename T, typename U, typename V, typename M>
  void fillJetTables(T const& jet, U const& /*cand*/, int32_t collisionIndex, int32_t candidateIndex, V& jetOutputTable, M& jetSubstructureOutputTable, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, float rho, std::map<int32_t, int32_t>& jetMap)
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

    jetOutputTable(collisionIndex, candidateIndex, jet.pt(), jet.phi(), jet.eta(), jet.y(), jet.r(), jet.area(), rho, jet.perpConeRho(), jet.tracksIds().size() + jet.candidatesIds().size()); // here we take the decision to keep the collision index consistent with the JE framework in case it is later needed to join to other tables. The candidate Index however can be linked to the HF tables
    jetSubstructureOutputTable(jetOutputTable.lastIndex(), energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, jet.nSub2DR(), jet.nSub1(), jet.nSub2(), pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, jet.angularity(), jet.ptLeadingConstituent(), splittingMatchesGeoVec, splittingMatchesPtVec, splittingMatchesHFVec, pairMatchesVec);
    jetMap.insert(std::make_pair(jet.globalIndex(), jetOutputTable.lastIndex()));
  }

  template <bool isMCP, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& jets, V const& /*candidates*/, M& collisionOutputTable, N& jetOutputTable, O& jetSubstructureOutputTable, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, std::map<int32_t, int32_t>& jetMap, std::map<int32_t, int32_t>& candidateMap, float jetPtMin, float eventWeight)
  {

    int nJetInCollision = 0;
    int32_t collisionIndex = -1;
    for (const auto& jet : jets) {
      if (jet.pt() < jetPtMin) {
        continue;
      }
      if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
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
          fillJetTables(jet, candidate, collisionIndex, candidateIndex, jetOutputTable, jetSubstructureOutputTable, splittingMatchesGeoVecVec, splittingMatchesPtVecVec, splittingMatchesHFVecVec, pairMatchesVecVec, candidate.rho(), jetMap);
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
      if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
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
      }
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
        if (!jetfindingutilities::isInEtaAcceptance(jet, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
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
        if (!jetfindingutilities::isInEtaAcceptance(jetMCP, configs.jetEtaMin, configs.jetEtaMax, configs.trackEtaMin, configs.trackEtaMax)) {
          continue;
        }
        for (const auto& jetRadiiValue : jetRadiiValues) {
          if (jetMCP.r() == round(jetRadiiValue * 100.0f)) {

            mcCollisionFlag[jetMCP.mcCollisionId()] = true;
            if constexpr (!isMCPOnly) {
              const auto collisionsPerMcCollision = collisions.sliceBy(preslices.CollisionsPerMcCollision, jetMCP.mcCollisionId());
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
            if constexpr (!isMCPOnly && (jethfutilities::isHFTable<P>() || jethfutilities::isHFMcTable<S>())) { // the matching of mcCollision to Collision is only done for HF tables
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

  void processClearMaps(aod::JetCollisions const&)
  {
    candidateMapping.clear();
    jetMappingData.clear();
    jetMappingDataSub.clear();
    jetMappingMCD.clear();
    candidateCollisionMapping.clear();
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processClearMaps, "process function that clears all the non-mcp maps in each dataframe", true);

  void processClearMapsMCP(aod::JetMcCollisions const& mcCollisions)
  {
    candidateMappingMCP.clear();
    jetMappingMCP.clear();
    candidateMcCollisionMapping.clear();
    for (auto mcCollision : mcCollisions) {
      products.hfMcOnlyCollisionsTable(mcCollision.posZ(), mcCollision.accepted(), mcCollision.attempted(), mcCollision.xsectGen(), mcCollision.xsectErr(), mcCollision.weight());
    }
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processClearMapsMCP, "process function that clears all the mcp maps in each dataframe", true);

  void processOutputCollisionsData(aod::JetCollisions const& collisions,
                                   JetTableData const& jets,
                                   CandidateCollisionTable const& canidateCollisions,
                                   CandidateTable const& candidates)
  {
    analyseHFCollisions<false, false>(collisions, collisions, canidateCollisions, canidateCollisions, jets, jets, candidates, candidates, configs.jetPtMinData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsData, "hf collision output data", false);

  void processOutputCollisionsDataSub(aod::JetCollisions const& collisions,
                                      JetTableDataSub const& jets,
                                      CandidateCollisionTable const& canidateCollisions,
                                      CandidateTable const& candidates)
  {
    analyseHFCollisions<false, false>(collisions, collisions, canidateCollisions, canidateCollisions, jets, jets, candidates, candidates, configs.jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsDataSub, "hf collision output data eventwise constituent subtracted", false);

  void processOutputCollisionsMc(soa::Join<aod::JetCollisions, aod::JMcCollisionLbs> const& collisions,
                                 aod::JetMcCollisions const& mcCollisions,
                                 JetTableMCD const& jetsMCD,
                                 JetTableMatchedMCP const& jetsMCP,
                                 CandidateCollisionTable const& canidateCollisions,
                                 CandidateMcCollisionTable const& canidateMcCollisions,
                                 CandidateTableMCD const& candidatesMCD,
                                 CandidateTableMCP const& candidatesMCP)
  {
    analyseHFCollisions<true, false>(collisions, mcCollisions, canidateCollisions, canidateMcCollisions, jetsMCD, jetsMCP, candidatesMCD, candidatesMCP, configs.jetPtMinMCD, configs.jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsMc, "hf collision output MC", false);

  void processOutputCollisionsMCPOnly(aod::JetMcCollisions const& mcCollisions,
                                      JetTableMCP const& jetsMCP,
                                      CandidateMcOnlyCollisionTable const& canidateMcCollisions,
                                      CandidateTableMCP const& candidatesMCP)
  {
    analyseHFCollisions<true, true>(mcCollisions, mcCollisions, canidateMcCollisions, canidateMcCollisions, jetsMCP, jetsMCP, candidatesMCP, candidatesMCP, 0.0, configs.jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCollisionsMCPOnly, "hf collision output MCP only", false);

  void processOutputCandidatesData(aod::JetCollision const&,
                                   JetTableData const& jets,
                                   CandidateTable const& candidates)
  {
    analyseCandidates<false, false>(jets, candidates, candidateMapping, configs.jetPtMinData);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesData, "hf candidate output data", false);

  void processOutputCandidatesDataSub(aod::JetCollision const&,
                                      JetTableDataSub const& jets,
                                      CandidateTable const& candidates)
  {
    analyseCandidates<false, false>(jets, candidates, candidateMapping, configs.jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesDataSub, "hf candidate output data eventwise constituent subtracted", false);

  void processOutputCandidatesMCD(aod::JetCollision const&,
                                  JetTableMCD const& jets,
                                  CandidateTableMCD const& candidates)
  {

    analyseCandidates<true, false>(jets, candidates, candidateMapping, configs.jetPtMinMCD);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputCandidatesMCD, "hf candidate output MCD", false);

  void processOutputCandidatesMCP(aod::JetMcCollision const&,
                                  JetTableMCP const& jets,
                                  CandidateTableMCP const& candidates)
  {
    analyseCandidates<false, true>(jets, candidates, candidateMappingMCP, configs.jetPtMinMCP);
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

  void processOutputJetsData(aod::JetCollision const& collision,
                             JetTableData const& jets,
                             soa::Join<CandidateTable, CandidateRhosTable> const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, products.collisionOutputTableData, products.jetOutputTableData, products.jetSubstructureOutputTableData, splittingMatchesGeoVecVecData, splittingMatchesPtVecVecData, splittingMatchesHFVecVecData, pairMatchesVecVecData, jetMappingData, candidateMapping, configs.jetPtMinData, 1.0);
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsData, "hf jet substructure output Data", false);

  void processOutputJetsDataSub(aod::JetCollision const& collision,
                                JetTableDataSub const& jets,
                                soa::Join<CandidateTable, CandidateRhosTable> const& candidates)
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

  void processOutputJetsMCD(aod::JetCollisionMCD const& collision,
                            JetTableMCD const& jets,
                            soa::Join<CandidateTableMCD, CandidateRhosTable> const& candidates)
  {
    analyseCharged<false>(collision, jets, candidates, products.collisionOutputTableMCD, products.jetOutputTableMCD, products.jetSubstructureOutputTableMCD, splittingMatchesGeoVecVecMCD, splittingMatchesPtVecVecMCD, splittingMatchesHFVecVecMCD, pairMatchesVecVecMCD, jetMappingMCD, candidateMapping, configs.jetPtMinMCD, collision.weight());
  }
  PROCESS_SWITCH(JetSubstructureHFOutputTask, processOutputJetsMCD, "hf jet substructure output MCD", false);

  void processOutputJetsMCP(aod::JetMcCollision const& collision,
                            JetTableMCP const& jets,
                            soa::Join<CandidateTableMCP, CandidateMCRhosTable> const& candidates)
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
};

#endif // PWGJE_TASKS_JETSUBSTRUCTUREHFOUTPUT_H_
