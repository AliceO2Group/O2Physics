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

/// \file jetGnnInference.cxx
/// \brief Task to produce tables for GNN inference
///
/// \author Changhwan Choi <changhwan.choi@cern.ch>

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/MlResponseHfTagging.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Tools/ML/MlResponse.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>

#include <onnxruntime_cxx_api.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace jet_gnn_inference
{
// Optional outputs for track-origin and future vertex prediction tables.
template <typename T>
struct OptionalProduces : Produces<T> {
};
template <>
struct OptionalProduces<void> {
};
} // namespace jet_gnn_inference

template <bool isMCD, typename JetTable, typename JetTaggingTable, typename JetGnnTrackIdsTable, typename JetGnnTracksTable, typename JetGnnInfTrackTable = void, typename JetGnnInfVertexTable = void>
struct JetGnnInferenceTask {
  static constexpr std::array<double, 3> defaultCutsMl = {0.5, 0.5, 0.5};
  static constexpr float invalidGnnScore = -999.f;
  static constexpr float largeNegativeNumber = -98.f;
  static constexpr float largePositiveNumber = 9999.f;

  Produces<JetTaggingTable> taggingTable;
  jet_gnn_inference::OptionalProduces<JetGnnInfTrackTable> gnnTrackTable;
  jet_gnn_inference::OptionalProduces<JetGnnInfVertexTable> gnnVertexTable;

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

  // ML configuration
  Configurable<int> nJetConst{"nJetConst", 100, "maximum number of jet constituents to be used for ML evaluation"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};

  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{5., 1000.}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cuts_ml::CutSmaller, cuts_ml::CutNot, cuts_ml::CutNot}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {defaultCutsMl.data(), 1, 3, {"pT bin 0"}, {"b", "c", "lf"}}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", 3, "Number of classes in ML model"};

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/c/chchoi"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"gnnbjettagger_LHC25a2b_v241.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // GNN configuration
  Configurable<float> jetpTMin{"jetpTMin", 5., "minimum jet pT for QA histograms"};
  Configurable<float> dbMin{"dbMin", -15., "minimum GNN Db"};
  Configurable<float> dbMax{"dbMax", 15., "maximum GNN Db"};
  Configurable<double> fC{"fC", 0.018, "Parameter f_c for D_b calculation"};
  Configurable<double> trackFC{"trackFC", 0.018, "Charm fraction for the track-origin D_b calculation"};
  Configurable<int64_t> nJetFeat{"nJetFeat", 4, "Number of jet GNN input features"};
  Configurable<int64_t> nTrkFeat{"nTrkFeat", 8, "Number of track GNN input features"};
  Configurable<int64_t> nTrkOrigin{"nTrkOrigin", 5, "Number of track origin categories"};
  Configurable<std::vector<float>> transformFeatureJetMean{"transformFeatureJetMean",
                                                           std::vector<float>{-999},
                                                           "Mean values for each GNN input feature (jet)"};
  Configurable<std::vector<float>> transformFeatureJetStdev{"transformFeatureJetStdev",
                                                            std::vector<float>{-999},
                                                            "Stdev values for each GNN input feature (jet)"};
  Configurable<std::vector<float>> transformFeatureTrkMean{"transformFeatureTrkMean",
                                                           std::vector<float>{-999},
                                                           "Mean values for each GNN input feature (track)"};
  Configurable<std::vector<float>> transformFeatureTrkStdev{"transformFeatureTrkStdev",
                                                            std::vector<float>{-999},
                                                            "Stdev values for each GNN input feature (track)"};
  Configurable<std::string> tfFuncTypeGNN{"tfFuncTypeGNN", "asinh", "Transformation function type for GNN"};

  o2::analysis::MlResponseHfTagging<float> bMlResponse;

  using JetTracksExt = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using AnalysisJets = std::conditional_t<isMCD, soa::Join<JetTable, aod::ChargedMCDetectorLevelJetFlavourDef>, JetTable>;
  using Collisions = std::conditional_t<isMCD, soa::Join<aod::JetCollisions, aod::JMcCollisionLbs, aod::JCollisionOutliers>, aod::JetCollisions>;
  using McCollisions = soa::Join<aod::JetMcCollisions, aod::JMcCollisionOutliers>;
  // Data reuses the already subscribed JCollisions table; only MCD requires
  // a separate context to dereference the matched MC collision's outlier flag.
  using CollisionContext = std::conditional_t<isMCD, McCollisions, aod::JCollisions>;

  o2::analysis::GNNBjetAllocator tensorAlloc;

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
    // with soa::Filtered: taggingTable must remain joinable to every input jet.
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

  template <typename Jet>
  void fillHistograms(Jet const& jet, float db, std::vector<float> const* trackScores = nullptr)
  {
    if (jet.pt() < jetpTMin) {
      return;
    }
    float dbRange;
    if (db < dbMin) {
      dbRange = 0.5; // underflow
    } else if (db < dbMax) {
      dbRange = 1.5; // in range
    } else {
      dbRange = 2.5; // overflow
    }
    registry.fill(HIST("h2_count_db"), 3.5, dbRange); // incl jet
    if constexpr (isMCD) {
      switch (jet.origin()) {
        case 2:
          registry.fill(HIST("h_db_b"), db);
          registry.fill(HIST("h2_count_db"), 0.5, dbRange);
          break;
        case 1:
          registry.fill(HIST("h_db_c"), db);
          registry.fill(HIST("h2_count_db"), 1.5, dbRange);
          break;
        case 0:
        case 3:
          registry.fill(HIST("h_db_lf"), db);
          registry.fill(HIST("h2_count_db"), 2.5, dbRange);
          break;
        default:
          LOGF(debug, "GNN inference, unexpected origin value: %d (%d)", jet.origin(), jet.globalIndex());
      }
    }
    registry.fill(HIST("h2_pt_db"), jet.pt(), db);

    // The tag-only jet type has no GNN-track array. Compile this block only
    // for full inference, and require MC truth only for the flavour split.
    if constexpr (requires { jet.tracksGnnIds(); }) {
      if (trackScores) {
        auto const gnnTracks = jet.template tracksGnn_as<JetGnnTracksTable>();
        if (trackScores->size() != gnnTracks.size()) {
          LOG(fatal) << "Track QA requires one score per GNN backbone row.";
        }
        for (auto const& track : gnnTracks) {
          std::size_t nodeIndex = track.trackIndex();
          if (nodeIndex >= trackScores->size()) {
            LOG(fatal) << "GNN trackIndex is out of range for the track QA scores.";
          }
          float score = (*trackScores)[nodeIndex];
          if (score == invalidGnnScore || std::isnan(score)) {
            continue;
          }
          registry.fill(HIST("h_trackDb"), score);
          if constexpr (isMCD) {
            switch (track.origin()) {
              case 2:
                registry.fill(HIST("h_trackDb_b"), score);
                break;
              case 1:
                registry.fill(HIST("h_trackDb_c"), score);
                break;
              case 0:
              case 3:
              case 4:
                registry.fill(HIST("h_trackDb_lf"), score);
                break;
              default:
                LOGF(debug, "GNN inference, unexpected track origin value: %d", track.origin());
            }
          }
        }
      }
    }
  }

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections);
    if (doprocessGnnTag && doprocessGnnFull) {
      LOG(fatal) << "Enable only one of processGnnTag and processGnnFull.";
    }
    if constexpr (std::is_void_v<JetGnnInfTrackTable>) {
      if (doprocessGnnFull) {
        LOG(fatal) << "processGnnFull requires a task with a track prediction output.";
      }
    }
    if (doprocessGnnTag || doprocessGnnFull) {
      // analyzeJetTrackInfo4GNN supplies four jet features and eight track features.
      if (nJetFeat != 4 || nTrkFeat != 8 || nClassesMl != 3 || nJetConst <= 0) {
        LOG(fatal) << "GNN jet tagging requires nJetFeat=4, nTrkFeat=8, nClassesMl=3 and nJetConst>0.";
      }
      if (doprocessGnnFull && (nTrkOrigin != 5 || !std::isfinite(trackFC.value) || trackFC < 0. || trackFC > 1.)) {
        LOG(fatal) << "Track-origin inference requires nTrkOrigin=5 and trackFC in [0, 1].";
      }
      if (transformFeatureJetMean.value.size() != static_cast<std::size_t>(nJetFeat.value) ||
          transformFeatureJetStdev.value.size() != static_cast<std::size_t>(nJetFeat.value) ||
          transformFeatureTrkMean.value.size() != static_cast<std::size_t>(nTrkFeat.value) ||
          transformFeatureTrkStdev.value.size() != static_cast<std::size_t>(nTrkFeat.value)) {
        LOG(fatal) << "GNN normalization vectors must match nJetFeat and nTrkFeat.";
      }
      bMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        o2::ccdb::CcdbApi ccdbApi;
        ccdbApi.init(ccdbUrl);
        bMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        bMlResponse.setModelPathsLocal(onnxFileNames);
      }
      bMlResponse.init();

      tensorAlloc = o2::analysis::GNNBjetAllocator(nJetFeat.value, nTrkFeat.value, nClassesMl.value, nTrkOrigin.value, transformFeatureJetMean.value, transformFeatureJetStdev.value, transformFeatureTrkMean.value, transformFeatureTrkStdev.value, nJetConst, tfFuncTypeGNN.value);

      registry.add("h2_count_db", "#it{D}_{b} underflow/overflow;Jet flavour;#it{D}_{b} range", {HistType::kTH2F, {{4, 0., 4.}, {3, 0., 3.}}});
      auto h2CountDb = registry.get<TH2>(HIST("h2_count_db"));
      h2CountDb->GetXaxis()->SetBinLabel(1, "b-jet");
      h2CountDb->GetXaxis()->SetBinLabel(2, "c-jet");
      h2CountDb->GetXaxis()->SetBinLabel(3, "lf-jet");
      h2CountDb->GetXaxis()->SetBinLabel(4, "incl jet");
      h2CountDb->GetYaxis()->SetBinLabel(1, "underflow");
      h2CountDb->GetYaxis()->SetBinLabel(2, "in range");
      h2CountDb->GetYaxis()->SetBinLabel(3, "overflow");

      registry.add("h_db_b", "#it{D}_{b} b-jet;#it{D}_{b}", {HistType::kTH1F, {{50, -20., 30.}}});
      registry.add("h_db_c", "#it{D}_{b} c-jet;#it{D}_{b}", {HistType::kTH1F, {{50, -20., 30.}}});
      registry.add("h_db_lf", "#it{D}_{b} lf-jet;#it{D}_{b}", {HistType::kTH1F, {{50, -20., 30.}}});
      registry.add("h2_pt_db", "#it{p}_{T} vs. #it{D}_{b};#it{p}_{T}^{ch jet} (GeV/#it{c}^{2});#it{D}_{b}", {HistType::kTH2F, {{100, 0., 200.}, {50, -20., 30.}}});
    }
    if (doprocessGnnFull) {
      registry.add("h_trackDb", "#it{d}_{b} inclusive track;#it{d}_{b}", {HistType::kTH1F, {{50, -20., 30.}}});
      if constexpr (isMCD) {
        registry.add("h_trackDb_b", "#it{d}_{b} b-track;#it{d}_{b}", {HistType::kTH1F, {{50, -20., 30.}}});
        registry.add("h_trackDb_c", "#it{d}_{b} c-track;#it{d}_{b}", {HistType::kTH1F, {{50, -20., 30.}}});
        registry.add("h_trackDb_lf", "#it{d}_{b} lf-track;#it{d}_{b}", {HistType::kTH1F, {{50, -20., 30.}}});
      }
    }
  }

  template <typename Jet, typename Tracks>
  float evaluateJet(Jet const& jet, Tracks const& tracks, std::vector<float>* trackScores = nullptr)
  {
    std::vector<std::vector<float>> trkFeat;
    if (trackScores) {
      trackScores->clear();
      if constexpr (!std::is_void_v<JetGnnInfTrackTable> && requires { jet.tracksGnnIds(); }) {
        // The definition backbone owns track selection and model-node order.
        // Do not reapply inference-side cuts or select a different track subset.
        auto const gnnTracks = jet.template tracksGnn_as<JetGnnTracksTable>();
        trackScores->assign(gnnTracks.size(), invalidGnnScore);
        if (gnnTracks.size() <= 1) {
          return invalidGnnScore;
        }
        if (gnnTracks.size() > static_cast<std::size_t>(nJetConst.value)) {
          LOG(fatal) << "GNN backbone contains " << gnnTracks.size() << " tracks, exceeding nJetConst=" << nJetConst.value << ". Increase the inference allocator limit to match jetGnnDefinition.";
        }
        trkFeat.reserve(gnnTracks.size());
        for (std::size_t nodeIndex = 0; nodeIndex < gnnTracks.size(); ++nodeIndex) {
          auto const& node = gnnTracks[nodeIndex];
          if (node.trackIndex() != nodeIndex || node.jetId() != jet.globalIndex()) {
            LOG(fatal) << "GNN backbone node order or jet reference does not match the jet's tracksGnnIds().";
          }
          auto const track = node.template track_as<Tracks>();
          int sign = jettaggingutilities::getGeoSign(jet, track);
          trkFeat.emplace_back(std::vector<float>{track.pt(), track.phi(), track.eta(), static_cast<float>(track.sign()), std::abs(track.dcaXY()) * sign, track.sigmadcaXY(), std::abs(track.dcaZ()) * sign, track.sigmadcaZ()});
        }
      } else {
        LOG(fatal) << "Track-origin inference requires a task with a track prediction output.";
      }
    } else {
      if (jet.tracksIds().size() <= 1 || !isAcceptedJet<Tracks>(jet)) {
        return invalidGnnScore;
      }
      jettaggingutilities::analyzeJetTrackInfo4GNN(jet, tracks, trkFeat, trackPtMin, trackDcaXYMax, trackDcaZMax, nJetConst);
      if (trkFeat.size() <= 1) {
        LOGF(debug, "GNN inference, fewer than two selected tracks (%d)", jet.globalIndex());
        return invalidGnnScore;
      }
    }

    std::vector<float> jetFeat{jet.pt(), jet.phi(), jet.eta(), jet.mass()};
    std::vector<float> feat;
    std::vector<Ort::Value> gnnInput;
    tensorAlloc.getGNNInput(jetFeat, trkFeat, feat, gnnInput);

    float db;
    if (trackScores) {
      // Exported head order: edge_outputs, node_outputs, flav_preds.
      // Retrieve jet and track logits in one model evaluation; leave edges unused.
      auto const outputs = bMlResponse.getModelOutputs(gnnInput, 0);
      if (outputs.empty()) {
        LOG(error) << "GNN returned no outputs; keeping invalid jet and track scores.";
        return invalidGnnScore;
      }
      if (outputs.size() != 3 || outputs.back().size() != static_cast<std::size_t>(nClassesMl.value) ||
          outputs[1].size() != trkFeat.size() * static_cast<std::size_t>(nTrkOrigin.value)) {
        LOG(fatal) << "Expected GNN outputs [edges, N*5 track-origin logits, 3 jet logits].";
      }
      db = jettaggingutilities::getDb(outputs.back(), fC);
      auto const& nodeLogits = outputs[1];
      std::vector<float> trackLogits(nTrkOrigin.value);
      for (std::size_t nodeIndex = 0; nodeIndex < trkFeat.size(); ++nodeIndex) {
        for (std::size_t category = 0; category < trackLogits.size(); ++category) {
          trackLogits[category] = nodeLogits[nodeIndex * trackLogits.size() + category];
        }
        float score = jettaggingutilities::getTrackDb(trackLogits, trackFC);
        (*trackScores)[nodeIndex] = std::isnan(score) ? invalidGnnScore : score;
      }
    } else {
      auto const modelOutput = bMlResponse.getModelOutput(gnnInput, 0);
      db = jettaggingutilities::getDb(modelOutput, fC);
    }
    if (std::isnan(db)) {
      LOGF(debug, "GNN inference, Db is NaN (%d)", jet.globalIndex());
      return 999.f;
    }
    return db;
  }

  void processDummy(aod::JetCollisions::iterator const&)
  {
  }
  PROCESS_SWITCH(JetGnnInferenceTask, processDummy, "Dummy process", true);

  void processGnnTag(typename Collisions::iterator const& collision, AnalysisJets const& jets, JetTracksExt const& tracks, CollisionContext const&)
  {
    bool acceptedCollision = isAcceptedCollision(collision);
    // Preserve a one-to-one positional join with Jets, including jets rejected
    // by collision/jet selection or with fewer than two GNN input tracks.
    for (auto const& jet : jets) {
      if (taggingTable.lastIndex() + 1 != jet.globalIndex()) {
        LOG(fatal) << "GNN inference requires jets ordered by collision with valid collision indices.";
      }
      float db = acceptedCollision ? evaluateJet(jet, tracks) : invalidGnnScore;
      taggingTable(0, -1.f, db);
      if (acceptedCollision) {
        fillHistograms(jet, db);
      }
    }
  }
  PROCESS_SWITCH(JetGnnInferenceTask, processGnnTag, "Run GNN inference and fill jet tagging tables", false);

  void processGnnFull(typename Collisions::iterator const& collision, soa::Join<AnalysisJets, JetGnnTrackIdsTable> const& jets, JetTracksExt const& tracks, JetGnnTracksTable const&, CollisionContext const&)
  {
    if constexpr (!std::is_void_v<JetGnnInfTrackTable>) {
      bool acceptedCollision = isAcceptedCollision(collision);
      std::vector<float> trackScores;
      for (auto const& jet : jets) {
        if (taggingTable.lastIndex() + 1 != jet.globalIndex()) {
          LOG(fatal) << "GNN inference requires jets ordered by collision with valid collision indices.";
        }
        float db = invalidGnnScore;
        if (acceptedCollision) {
          db = evaluateJet(jet, tracks, &trackScores);
        } else {
          // A different definition-side selection may have left backbone rows.
          // Keep one invalid prediction per row instead of breaking the join.
          trackScores.assign(jet.tracksGnnIds().size(), invalidGnnScore);
        }
        taggingTable(0, -1.f, db);
        auto const nodes = jet.template tracksGnn_as<JetGnnTracksTable>();
        for (std::size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex) {
          // Thin prediction tables join positionally with the input backbone.
          if (nodes[nodeIndex].globalIndex() != gnnTrackTable.lastIndex() + 1) {
            LOG(fatal) << "GNN backbone rows are not in jet/node order; refusing a misaligned track prediction table.";
          }
          gnnTrackTable(trackScores[nodeIndex]);
        }
        if (acceptedCollision) {
          fillHistograms(jet, db, &trackScores);
        }
      }
    }
  }
  PROCESS_SWITCH(JetGnnInferenceTask, processGnnFull, "Run GNN inference and fill jet tagging and track-origin prediction tables", false);
};

// One task per jet type owns both outputs. Process switches select tag-only or
// full inference without registering duplicate producers of the jet tags.
using JetGnnInferenceData = JetGnnInferenceTask<false, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, aod::ChargedJetTags, aod::ChargedJetGnnTrackIds, aod::ChargedJetGnnTracks, aod::ChargedJetGnnTrackOriginTags>;
using JetGnnInferenceMcd = JetGnnInferenceTask<true, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedMCDetectorLevelJetTags, aod::ChargedMCDetectorLevelJetGnnTrackIds, soa::Join<aod::ChargedMCDetectorLevelJetGnnTracks, aod::ChargedMCDetectorLevelJetTrackOriginDef>, aod::ChargedMCDetectorLevelJetGnnTrackOriginTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetGnnInferenceData>(
      cfgc, SetDefaultProcesses{}, TaskName{"jet-gnn-inference-data"}), // o2-linter: disable=name/o2-task (custom device name for templated task alias)
    adaptAnalysisTask<JetGnnInferenceMcd>(
      cfgc, SetDefaultProcesses{}, TaskName{"jet-gnn-inference-mcd"})}; // o2-linter: disable=name/o2-task (custom device name for templated task alias)
}
