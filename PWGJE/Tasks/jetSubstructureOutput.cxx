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

#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

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
  Produces<aod::CMCPJetMCCOs> mcCollisionOutputTable;

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

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  Preslice<soa::Join<aod::ChargedSPs, aod::ChargedSPsMatchedToChargedEventWiseSubtractedSPs>> splittingsPerJetData = aod::chargedsplitting::jetId;
  Preslice<soa::Join<aod::ChargedEventWiseSubtractedSPs, aod::ChargedEventWiseSubtractedSPsMatchedToChargedSPs>> splittingsPerJetDataSub = aod::chargedeventwisesubtractedsplitting::jetId;
  Preslice<soa::Join<aod::ChargedMCDetectorLevelSPs, aod::ChargedMCDetectorLevelSPsMatchedToChargedMCParticleLevelSPs>> splittingsPerJetMCD = aod::chargedmcdetectorlevelsplitting::jetId;
  Preslice<soa::Join<aod::ChargedMCParticleLevelSPs, aod::ChargedMCParticleLevelSPsMatchedToChargedMCDetectorLevelSPs>> splittingsPerJetMCP = aod::chargedmcparticlelevelsplitting::jetId;

  Preslice<soa::Join<aod::ChargedPRs, aod::ChargedPRsMatchedToChargedEventWiseSubtractedPRs>> pairsPerJetData = aod::chargedpair::jetId;
  Preslice<soa::Join<aod::ChargedEventWiseSubtractedPRs, aod::ChargedEventWiseSubtractedPRsMatchedToChargedPRs>> pairsPerJetDataSub = aod::chargedeventwisesubtractedpair::jetId;
  Preslice<soa::Join<aod::ChargedMCDetectorLevelPRs, aod::ChargedMCDetectorLevelPRsMatchedToChargedMCParticleLevelPRs>> pairsPerJetMCD = aod::chargedmcdetectorlevelpair::jetId;
  Preslice<soa::Join<aod::ChargedMCParticleLevelPRs, aod::ChargedMCParticleLevelPRsMatchedToChargedMCDetectorLevelPRs>> pairsPerJetMCP = aod::chargedmcparticlelevelpair::jetId;

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
  void fillJetTables(T const& jet, int32_t collisionIndex, U& jetOutputTable, V& jetSubstructureOutputTable, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, float rho, std::map<int32_t, int32_t>& jetMapping)
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
    jetOutputTable(collisionIndex, collisionIndex, jet.pt(), jet.phi(), jet.eta(), jet.y(), jet.r(), jet.area(), rho, jet.perpConeRho(), jet.tracksIds().size()); // second collision index is a dummy coloumn mirroring the hf candidate
    jetSubstructureOutputTable(jetOutputTable.lastIndex(), energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec, jet.nSub2DR(), jet.nSub1(), jet.nSub2(), pairJetPtVec, pairJetEnergyVec, pairJetThetaVec, pairJetPerpCone1PtVec, pairJetPerpCone1EnergyVec, pairJetPerpCone1ThetaVec, pairPerpCone1PerpCone1PtVec, pairPerpCone1PerpCone1EnergyVec, pairPerpCone1PerpCone1ThetaVec, pairPerpCone1PerpCone2PtVec, pairPerpCone1PerpCone2EnergyVec, pairPerpCone1PerpCone2ThetaVec, jet.angularity(), jet.ptLeadingConstituent(), splittingMatchesGeoVec, splittingMatchesPtVec, splittingMatchesHFVec, pairMatchesVec);
    jetMapping.insert(std::make_pair(jet.globalIndex(), jetOutputTable.lastIndex()));
  }

  template <bool isMCP, typename T, typename U, typename V, typename M, typename N>
  void analyseCharged(T const& collision, U const& jets, V& collisionOutputTable, M& jetOutputTable, N& jetSubstructureOutputTable, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, std::map<int32_t, int32_t>& jetMapping, float jetPtMin, float eventWeight)
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
              centrality = collision.centFT0M();
              eventSel = collision.eventSel();
            }
            collisionOutputTable(collision.posZ(), centrality, eventSel, eventWeight);
            collisionIndex = collisionOutputTable.lastIndex();
          }
          nJetInCollision++;
          fillJetTables(jet, collisionIndex, jetOutputTable, jetSubstructureOutputTable, splittingMatchesGeoVecVec, splittingMatchesPtVecVec, splittingMatchesHFVecVec, pairMatchesVecVec, collision.rho(), jetMapping);
        }
      }
    }
  }

  template <typename T, typename U, typename V, typename M, typename N>
  void analyseSubstructureMatched(T const& jets, U const& allSplittings, V const& allPairs, M const& splittingsSlicer, N const& pairsSlicer, std::vector<std::vector<int32_t>>& splittingMatchesGeoVecVec, std::vector<std::vector<int32_t>>& splittingMatchesPtVecVec, std::vector<std::vector<int32_t>>& splittingMatchesHFVecVec, std::vector<std::vector<int32_t>>& pairMatchesVecVec, float jetPtMin)
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

          auto splittings = allSplittings.sliceBy(splittingsSlicer, jet.globalIndex());
          fillSplittingMatchingVectors(splittings, jet.globalIndex(), splittingMatchesGeoVecVec, splittingMatchesPtVecVec, splittingMatchesHFVecVec);

          auto pairs = allPairs.sliceBy(pairsSlicer, jet.globalIndex());
          fillPairMatchingVectors(pairs, jet.globalIndex(), pairMatchesVecVec);
        }
      }
    }
  }

  template <typename T, typename U>
  void analyseJetMatched(T const& jets, std::map<int32_t, int32_t>& jetMapping, std::map<int32_t, int32_t>& jetTagMapping, U& matchingOutputTable, float jetPtMin)
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
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processClearMaps, "process function that clears all the non-mcp maps in each dataframe", true);

  void processClearMapsMCP(aod::JetMcCollisions const& mcCollisions)
  {
    jetMappingMCP.clear();
    for (auto mcCollision : mcCollisions) {
      mcCollisionOutputTable(mcCollision.posZ(), mcCollision.accepted(), mcCollision.attempted(), mcCollision.xsectGen(), mcCollision.xsectErr(), mcCollision.weight());
    }
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processClearMapsMCP, "process function that clears all the mcp maps in each dataframe", true);

  void processOutputSubstructureMatchingData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets> const& jets,
                                             soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets> const& jetsSub,
                                             soa::Join<aod::ChargedSPs, aod::ChargedSPsMatchedToChargedEventWiseSubtractedSPs> const& splittingsData,
                                             soa::Join<aod::ChargedEventWiseSubtractedSPs, aod::ChargedEventWiseSubtractedSPsMatchedToChargedSPs> const& splittingsDataSub,
                                             soa::Join<aod::ChargedPRs, aod::ChargedPRsMatchedToChargedEventWiseSubtractedPRs> const& pairsData,
                                             soa::Join<aod::ChargedEventWiseSubtractedPRs, aod::ChargedEventWiseSubtractedPRsMatchedToChargedPRs> const& pairsDataSub)
  {
    splittingMatchesGeoVecVecData.assign(jets.size(), {});
    splittingMatchesPtVecVecData.assign(jets.size(), {});
    splittingMatchesHFVecVecData.assign(jets.size(), {});
    pairMatchesVecVecData.assign(jets.size(), {});
    analyseSubstructureMatched(jets, splittingsData, pairsData, splittingsPerJetData, pairsPerJetData, splittingMatchesGeoVecVecData, splittingMatchesPtVecVecData, splittingMatchesHFVecVecData, pairMatchesVecVecData, jetPtMinData);
    splittingMatchesGeoVecVecDataSub.assign(jetsSub.size(), {});
    splittingMatchesPtVecVecDataSub.assign(jetsSub.size(), {});
    splittingMatchesHFVecVecDataSub.assign(jetsSub.size(), {});
    pairMatchesVecVecDataSub.assign(jetsSub.size(), {});
    analyseSubstructureMatched(jetsSub, splittingsDataSub, pairsDataSub, splittingsPerJetDataSub, pairsPerJetDataSub, splittingMatchesGeoVecVecDataSub, splittingMatchesPtVecVecDataSub, splittingMatchesHFVecVecDataSub, pairMatchesVecVecDataSub, jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputSubstructureMatchingData, "substructure matching output Data", false);

  void processOutputData(soa::Join<aod::JetCollisions, aod::BkgChargedRhos>::iterator const& collision,
                         soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::CJetSSs> const& jets)
  {
    analyseCharged<false>(collision, jets, collisionOutputTableData, jetOutputTableData, jetSubstructureOutputTableData, splittingMatchesGeoVecVecData, splittingMatchesPtVecVecData, splittingMatchesHFVecVecData, pairMatchesVecVecData, jetMappingData, jetPtMinData, 1.0);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputData, "jet substructure output Data", false);

  void processOutputDataSub(soa::Join<aod::JetCollisions, aod::BkgChargedRhos>::iterator const& collision,
                            soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::CEWSJetSSs> const& jets)
  {
    analyseCharged<false>(collision, jets, collisionOutputTableDataSub, jetOutputTableDataSub, jetSubstructureOutputTableDataSub, splittingMatchesGeoVecVecDataSub, splittingMatchesPtVecVecDataSub, splittingMatchesHFVecVecDataSub, pairMatchesVecVecDataSub, jetMappingDataSub, jetPtMinDataSub, 1.0);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputDataSub, "jet substructure output event-wise subtracted Data", false);

  void processOutputJetMatchingData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets> const& jets,
                                    soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets> const& jetsSub)
  {
    analyseJetMatched(jets, jetMappingData, jetMappingDataSub, jetMatchingOutputTableData, jetPtMinData);
    analyseJetMatched(jetsSub, jetMappingDataSub, jetMappingData, jetMatchingOutputTableDataSub, jetPtMinDataSub);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputJetMatchingData, "jet matching output Data", false);

  void processOutputSubstructureMatchingMC(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& jetsMCD,
                                           soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& jetsMCP,
                                           soa::Join<aod::ChargedMCDetectorLevelSPs, aod::ChargedMCDetectorLevelSPsMatchedToChargedMCParticleLevelSPs> const& splittingsMCD,
                                           soa::Join<aod::ChargedMCParticleLevelSPs, aod::ChargedMCParticleLevelSPsMatchedToChargedMCDetectorLevelSPs> const& splittingsMCP,
                                           soa::Join<aod::ChargedMCDetectorLevelPRs, aod::ChargedMCDetectorLevelPRsMatchedToChargedMCParticleLevelPRs> const& pairsMCD,
                                           soa::Join<aod::ChargedMCParticleLevelPRs, aod::ChargedMCParticleLevelPRsMatchedToChargedMCDetectorLevelPRs> const& pairsMCP)
  {
    splittingMatchesGeoVecVecMCD.assign(jetsMCD.size(), {});
    splittingMatchesPtVecVecMCD.assign(jetsMCD.size(), {});
    splittingMatchesHFVecVecMCD.assign(jetsMCD.size(), {});
    pairMatchesVecVecMCD.assign(jetsMCD.size(), {});
    analyseSubstructureMatched(jetsMCD, splittingsMCD, pairsMCD, splittingsPerJetMCD, pairsPerJetMCD, splittingMatchesGeoVecVecMCD, splittingMatchesPtVecVecMCD, splittingMatchesHFVecVecMCD, pairMatchesVecVecMCD, jetPtMinMCD);
    splittingMatchesGeoVecVecMCP.assign(jetsMCP.size(), {});
    splittingMatchesPtVecVecMCP.assign(jetsMCP.size(), {});
    splittingMatchesHFVecVecMCP.assign(jetsMCP.size(), {});
    pairMatchesVecVecMCP.assign(jetsMCP.size(), {});
    analyseSubstructureMatched(jetsMCP, splittingsMCP, pairsMCP, splittingsPerJetMCP, pairsPerJetMCP, splittingMatchesGeoVecVecMCP, splittingMatchesPtVecVecMCP, splittingMatchesHFVecVecMCP, pairMatchesVecVecMCP, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputSubstructureMatchingMC, "substructure matching output MC", false);

  void processOutputMCD(soa::Join<aod::JetCollisionsMCD, aod::BkgChargedRhos>::iterator const& collision,
                        soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::CMCDJetSSs> const& jets)
  {
    analyseCharged<false>(collision, jets, collisionOutputTableMCD, jetOutputTableMCD, jetSubstructureOutputTableMCD, splittingMatchesGeoVecVecMCD, splittingMatchesPtVecVecMCD, splittingMatchesHFVecVecMCD, pairMatchesVecVecMCD, jetMappingMCD, jetPtMinMCD, collision.weight());
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCD, "jet substructure output MCD", false);

  void processOutputMCP(soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>::iterator const& collision,
                        soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::CMCPJetSSs> const& jets)
  {
    analyseCharged<true>(collision, jets, collisionOutputTableMCP, jetOutputTableMCP, jetSubstructureOutputTableMCP, splittingMatchesGeoVecVecMCP, splittingMatchesPtVecVecMCP, splittingMatchesHFVecVecMCP, pairMatchesVecVecMCP, jetMappingMCP, jetPtMinMCP, collision.weight());
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCP, "jet substructure output MCP", false);

  void processOutputJetMatchingMC(soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets> const& jetsMCD,
                                  soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets> const& jetsMCP)
  {
    analyseJetMatched(jetsMCD, jetMappingMCD, jetMappingMCP, jetMatchingOutputTableMCD, jetPtMinMCD);
    analyseJetMatched(jetsMCP, jetMappingMCP, jetMappingMCD, jetMatchingOutputTableMCP, jetPtMinMCP);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputJetMatchingMC, "jet matching output MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureOutputTask>(
    cfgc, TaskName{"jet-substructure-output"})};
}
