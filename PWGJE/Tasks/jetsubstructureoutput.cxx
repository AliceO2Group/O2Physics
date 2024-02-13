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

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

struct JetSubstructureOutputTask {

  Produces<aod::ChargedJetCollisionOutputs> collisionOutputTableData;
  Produces<aod::ChargedJetOutputs> jetOutputTableData;
  Produces<aod::ChargedJetSubstructureOutputs> jetSubstructureOutputTableData;
  Produces<aod::ChargedEventWiseSubtractedJetCollisionOutputs> collisionOutputTableDataSub;
  Produces<aod::ChargedEventWiseSubtractedJetOutputs> jetOutputTableDataSub;
  Produces<aod::ChargedEventWiseSubtractedJetSubstructureOutputs> jetSubstructureOutputTableDataSub;
  Produces<aod::ChargedMCDetectorLevelJetCollisionOutputs> collisionOutputTableMCD;
  Produces<aod::ChargedMCDetectorLevelJetOutputs> jetOutputTableMCD;
  Produces<aod::ChargedMCDetectorLevelJetSubstructureOutputs> jetSubstructureOutputTableMCD;
  Produces<aod::ChargedMCParticleLevelJetCollisionOutputs> collisionOutputTableMCP;
  Produces<aod::ChargedMCParticleLevelJetOutputs> jetOutputTableMCP;
  Produces<aod::ChargedMCParticleLevelJetSubstructureOutputs> jetSubstructureOutputTableMCP;

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track pseudorapidity"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track pseudorapidity"};

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;
  }

  Filter jetSelection = aod::jet::pt >= jetPtMin;

  template <typename T, typename U, typename V, typename M>
  void fillTables(T const& jet, int32_t collisionIndex, U& collisionOutputTable, V& jetOutputTable, M& jetSubstructureOutputTable, std::vector<int> geoMatching, std::vector<int> ptMatching, std::vector<int> candMatching)
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
    jetOutputTable(collisionIndex, -1, geoMatching, ptMatching, candMatching, jet.pt(), jet.phi(), jet.eta(), jet.r(), jet.tracks().size());
    jetSubstructureOutputTable(jetOutputTable.lastIndex(), energyMotherVec, ptLeadingVec, ptSubLeadingVec, thetaVec);
  }

  template <bool hasMatching, typename T, typename U, typename V, typename M, typename N, typename O>
  void analyseCharged(T const& collision, U const& jets, V const& jetsTag, M& collisionOutputTable, N& jetOutputTable, O& jetSubstructureOutputTable)
  {

    std::vector<int> candMatching{-1};
    int nJetInCollision = 0;
    int32_t collisionIndex = -1;
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          if constexpr (hasMatching) {
            if (jet.has_matchedJetGeo()) {
              for (auto& jetTag : jet.template matchedJetGeo_as<V>()) {
                geoMatching.push_back(jetTag.globalIndex());
              }
            }
            if (jet.has_matchedJetPt()) {
              for (auto& jetTag : jet.template matchedJetPt_as<V>()) {
                ptMatching.push_back(jetTag.globalIndex());
              }
            }
          }
          if (nJetInCollision == 0) {
            collisionOutputTable(collision.posZ(), collision.centrality(), collision.eventSel());
            collisionIndex = collisionOutputTable.lastIndex();
          }
          nJetInCollision++;
          fillTables(jet, collisionIndex, collisionOutputTable, jetOutputTable, jetSubstructureOutputTable, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }

  void processDummy(JetCollisions const& collisions) {}
  PROCESS_SWITCH(JetSubstructureOutputTask, processDummy, "Dummy process function turned on by default", true);

  void processOutputData(JetCollision const& collision,
                         soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetSubstructures>> const& jets,
                         JetTracks const& tracks)
  {
    analyseCharged<false>(collision, jets, jets, collisionOutputTableData, jetOutputTableData, jetSubstructureOutputTableData);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputData, "jet substructure output Data", false);

  void processOutputDataSub(JetCollision const& collision,
                            soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetSubstructures, aod::ChargedJetsMatchedToChargedEventWiseSubtractedJets>> const& jets,
                            soa::Filtered<soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents, aod::ChargedEventWiseSubtractedJetSubstructures, aod::ChargedEventWiseSubtractedJetsMatchedToChargedJets>> const& jetsSub,
                            JetTracks const& tracks)
  {
    analyseCharged<true>(collision, jets, jetsSub, collisionOutputTableData, jetOutputTableData, jetSubstructureOutputTableData);
    analyseCharged<true>(collision, jetsSub, jets, collisionOutputTableDataSub, jetOutputTableDataSub, jetSubstructureOutputTableDataSub);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputDataSub, "jet substructure output event-wise subtracted Data", false);

  void processOutputMCD(JetCollision const& collision,
                        soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetSubstructures, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& jets,
                        aod::ChargedMCParticleLevelJets const& jetsTag,
                        JetTracks const& tracks)
  {
    analyseCharged<true>(collision, jets, jetsTag, collisionOutputTableMCD, jetOutputTableMCD, jetSubstructureOutputTableMCD);
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCD, "jet substructure output MCD", false);

  void processOutputMCP(JetMcCollision const& collision,
                        soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetSubstructures, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>> const& jets,
                        aod::ChargedMCDetectorLevelJets const& jetsTag,
                        JetParticles const& particles)
  {
    std::vector<int> candMatching{-1};
    for (const auto& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      for (const auto& jetRadiiValue : jetRadiiValues) {
        if (jet.r() == round(jetRadiiValue * 100.0f)) {
          std::vector<int> geoMatching;
          std::vector<int> ptMatching;
          if (jet.has_matchedJetGeo()) {
            for (auto& jetTag : jet.template matchedJetGeo_as<aod::ChargedMCDetectorLevelJets>()) {
              geoMatching.push_back(jetTag.globalIndex());
            }
          }
          if (jet.has_matchedJetPt()) {
            for (auto& jetTag : jet.template matchedJetPt_as<aod::ChargedMCDetectorLevelJets>()) {
              ptMatching.push_back(jetTag.globalIndex());
            }
          }
          fillTables(jet, -1, collisionOutputTableMCP, jetOutputTableMCP, jetSubstructureOutputTableMCP, geoMatching, ptMatching, candMatching);
        }
      }
    }
  }
  PROCESS_SWITCH(JetSubstructureOutputTask, processOutputMCP, "jet substructure output MCP", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetSubstructureOutputTask>(
    cfgc, TaskName{"jet-substructure-output"})};
}
