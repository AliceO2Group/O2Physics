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

/// \file jetGnnDefinition.cxx
/// \brief Task to produce GNN track indices and MC prediction label definitions
///
/// \author Changhwan Choi <changhwan.choi@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace jet_gnn_definition
{
// Optional outputs for MC track-origin and vertex definition tables.
template <typename T>
struct OptionalProduces : Produces<T> {
};
template <>
struct OptionalProduces<void> {
};
} // namespace jet_gnn_definition

template <bool isMCD, typename JetTable, typename JetGnnDefTrackTable = void, typename JetGnnDefVertexTable = void>
struct JetGnnDefinitionTask {
  static constexpr float largeNegativeNumber = -98.f;
  static constexpr float largePositiveNumber = 9999.f;

  Produces<std::conditional_t<isMCD, aod::ChargedMCDetectorLevelJetGnnTracks, aod::ChargedJetGnnTracks>> gnnTrackTable;
  Produces<std::conditional_t<isMCD, aod::ChargedMCDetectorLevelJetGnnTrackIds, aod::ChargedJetGnnTrackIds>> gnnTrackIndicesTable;
  jet_gnn_definition::OptionalProduces<std::conditional_t<isMCD, JetGnnDefTrackTable, void>> gnnDefTrackTable;
  jet_gnn_definition::OptionalProduces<std::conditional_t<isMCD, JetGnnDefVertexTable, void>> gnnDefVertexTable;

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  std::vector<int> eventSelectionBits;

  // Track selection for GNN inputs
  Configurable<float> trackDcaXYMax{"trackDcaXYMax", 1, "maximum DCA xy acceptance for tracks [cm]"};
  Configurable<float> trackDcaZMax{"trackDcaZMax", 2, "maximum DCA z acceptance for tracks [cm]"};

  // Match the jetFilter and isAcceptedJet criteria in bjetTaggingGnn.cxx.
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 300.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.9, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.9, "maximum jet pseudorapidity"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<int> nJetConst{"nJetConst", 100, "maximum number of jet constituents to be used for ML evaluation"};
  Configurable<float> vtxRes{"vtxRes", 0.01, "Vertex position resolution (cluster size) for GNN vertex predictions (cm)"};
  Configurable<bool> searchUpToQuark{"searchUpToQuark", false, "Flag whether to search for track origins up to the quark level"};

  using JetTracksExt = std::conditional_t<isMCD, soa::Join<aod::JetTracksMCD, aod::JTrackExtras>, soa::Join<aod::JetTracks, aod::JTrackExtras>>;
  using AnalysisJets = JetTable;
  using CollisionsData = aod::JetCollisions;
  using CollisionsMcd = soa::Join<aod::JetCollisions, aod::JMcCollisionLbs, aod::JCollisionOutliers>;
  using McCollisions = soa::Join<aod::JetMcCollisions, aod::JMcCollisionOutliers>;
  using JetParticlesWithIds = soa::Join<aod::JetParticles, aod::JMcParticlePIs>;

  template <typename Collision>
  bool isAcceptedCollision(Collision const& collision)
  {
    if (!(std::abs(collision.posZ()) < vertexZCut) ||
        !jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return false;
    }
    if constexpr (isMCD) {
      if (collision.isOutlier() ||
          (collision.has_mcCollision() && collision.template mcCollision_as<McCollisions>().isOutlier())) {
        return false;
      }
    }
    return true;
  }

  template <typename Tracks, typename Jet>
  bool isAcceptedJet(Jet const& jet)
  {
    // Apply the analysis jetFilter as a predicate rather than dropping rows
    // with soa::Filtered: GnnTrackIds must remain joinable to every input jet.
    if (!(jet.pt() >= jetPtMin && jet.pt() < jetPtMax &&
          jet.eta() < jetEtaMax - jet.r() / 100.f &&
          jet.eta() > jetEtaMin + jet.r() / 100.f)) {
      return false;
    }

    if (jetAreaFractionMin > largeNegativeNumber &&
        jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
      return false;
    }

    bool checkConstituentMinPt = leadingConstituentPtMin > largeNegativeNumber;
    bool checkConstituentMaxPt = leadingConstituentPtMax < largePositiveNumber;
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      return true;
    }

    bool isMinLeadingConstituent = !checkConstituentMinPt;
    bool isMaxLeadingConstituent = true;
    for (auto const& constituent : jet.template tracks_as<Tracks>()) {
      double pt = constituent.pt();
      if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
        isMinLeadingConstituent = true;
      }
      if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
        isMaxLeadingConstituent = false;
      }
    }
    return isMinLeadingConstituent && isMaxLeadingConstituent;
  }

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections);
    constexpr int MaxGnnTracks = std::numeric_limits<uint8_t>::max() + 1;
    if ((doprocessMcd || doprocessData) && (nJetConst < 1 || nJetConst > MaxGnnTracks)) {
      LOGF(fatal, "nJetConst must be in [1, %d] to fit the uint8_t trackIndex", MaxGnnTracks);
    }
    if ((doprocessMcd && !isMCD) || (doprocessData && isMCD)) {
      LOGF(fatal, "Use processMcd for MC detector-level jets and processData for data jets");
    }
  }

  void processDummy(aod::JetCollisions::iterator const&)
  {
  }
  PROCESS_SWITCH(JetGnnDefinitionTask, processDummy, "Dummy process", true);

  template <typename GetTrackLabels>
  void fillTrackTables(AnalysisJets const& jets, JetTracksExt const& tracks, bool acceptedCollision, GetTrackLabels const& getTrackLabels)
  {
    for (auto const& jet : jets) {
      // Collision grouping must preserve the positional extension of Jets.
      if (gnnTrackIndicesTable.lastIndex() + 1 != jet.globalIndex()) {
        LOG(fatal) << "GNN definition requires jets ordered by collision with valid collision indices.";
      }
      std::vector<int32_t> gnnTrackIds;
      if (!acceptedCollision || jet.tracksIds().size() <= 1 || !isAcceptedJet<JetTracksExt>(jet)) {
        gnnTrackIndicesTable(gnnTrackIds);
        continue;
      }

      // Reuse the inference selection, including replacement of low-significance
      // tracks at nJetConst. Keep the constituent position for each model node.
      std::vector<std::vector<float>> trackFeatures;
      std::vector<std::size_t> constituentIndices;
      jettaggingutilities::analyzeJetTrackInfo4GNN(jet, tracks, trackFeatures, trackPtMin, trackDcaXYMax, trackDcaZMax, nJetConst, &constituentIndices);
      if (constituentIndices.size() <= 1) {
        gnnTrackIndicesTable(gnnTrackIds);
        continue;
      }

      std::vector<int> origins;
      std::vector<int> vertexIndices;
      if constexpr (isMCD) {
        getTrackLabels(jet, origins, vertexIndices);
      }

      auto const constituents = jet.template tracks_as<JetTracksExt>();
      gnnTrackIds.reserve(constituentIndices.size());
      for (std::size_t nodeIndex = 0; nodeIndex < constituentIndices.size(); ++nodeIndex) {
        auto const constituentIndex = constituentIndices[nodeIndex];
        auto const& track = constituents[constituentIndex];
        gnnTrackTable(jet.globalIndex(), track.globalIndex(), static_cast<uint8_t>(nodeIndex));
        if constexpr (isMCD) {
          // The backbone and truth table share exactly the same row order.
          gnnDefTrackTable(static_cast<uint8_t>(origins[constituentIndex]));
          if constexpr (!std::is_void_v<JetGnnDefVertexTable>) {
            gnnDefVertexTable(static_cast<int8_t>(vertexIndices[constituentIndex]));
          }
        }
        gnnTrackIds.push_back(gnnTrackTable.lastIndex());
      }
      gnnTrackIndicesTable(gnnTrackIds);
    }
  }

  // PROCESS_SWITCH compiles both entry points for each task specialization,
  // including the disabled one. Their table-specific code needs compile-time guards.
  void processData(CollisionsData::iterator const& collision, AnalysisJets const& jets, JetTracksExt const& tracks)
  {
    if constexpr (!isMCD) {
      bool acceptedCollision = isAcceptedCollision(collision);
      fillTrackTables(jets, tracks, acceptedCollision, [](auto const&, std::vector<int>&, std::vector<int>&) {});
    }
  }
  PROCESS_SWITCH(JetGnnDefinitionTask, processData, "Produce GNN track indices for data jets", false);

  void processMcd(CollisionsMcd::iterator const& collision, AnalysisJets const& jets, JetTracksExt const& tracks, McCollisions const&, JetParticlesWithIds const& particles, aod::McParticles const& originalParticles)
  {
    if constexpr (isMCD) {
      bool acceptedCollision = isAcceptedCollision(collision);
      fillTrackTables(jets, tracks, acceptedCollision, [&](auto const& jet, std::vector<int>& origins, std::vector<int>& vertexIndices) {
        // vertexClustering labels the full constituent array, as in the training
        // tree creator. Project these labels onto the selected GNN nodes.
        std::unordered_map<std::string, std::vector<int>> trackLabels{{"trkVtxIndex", {}}, {"trkOrigin", {}}};
        if (collision.has_mcCollision()) {
          jettaggingutilities::vertexClustering(collision.template mcCollision_as<McCollisions>(), jet, tracks, particles, originalParticles, trackLabels, searchUpToQuark, vtxRes, trackPtMin);
          origins = std::move(trackLabels.at("trkOrigin"));
          vertexIndices = std::move(trackLabels.at("trkVtxIndex"));
        } else {
          // No matched MC collision: unknown origin (0) and no vertex label (-1).
          origins.assign(jet.tracksIds().size(), 0);
          vertexIndices.assign(jet.tracksIds().size(), -1);
        }
      });
    }
  }
  PROCESS_SWITCH(JetGnnDefinitionTask, processMcd, "Produce GNN track indices, origin and vertex definitions for MC detector-level jets", false);
};

using JetGnnDefinitionData = JetGnnDefinitionTask<false, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>;
using JetGnnDefinitionMcd = JetGnnDefinitionTask<true, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedMCDetectorLevelJetTrackOriginDef, aod::ChargedMCDetectorLevelJetVertexIndexDef>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetGnnDefinitionData>(
      cfgc, SetDefaultProcesses{}, TaskName{"jet-gnn-definition-data"}), // o2-linter: disable=name/o2-task (custom device name for templated task alias)
    adaptAnalysisTask<JetGnnDefinitionMcd>(
      cfgc, SetDefaultProcesses{}, TaskName{"jet-gnn-definition-mcd"})}; // o2-linter: disable=name/o2-task (custom device name for templated task alias)
}
