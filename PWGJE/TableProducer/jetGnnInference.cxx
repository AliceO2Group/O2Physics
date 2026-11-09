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

#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/MlResponseHfTagging.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Tools/ML/MlResponse.h"

#include <CCDB/CcdbApi.h>
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

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace jet_gnn_inference
{
// Optional outputs for future track-origin and vertex prediction tables.
template <typename T>
struct OptionalProduces : Produces<T> {
};
template <>
struct OptionalProduces<void> {
};
} // namespace jet_gnn_inference

template <bool isMCD, typename JetTable, typename JetTaggingTable, typename JetGnnInfTrackTable = void, typename JetGnnInfVertexTable = void>
struct JetGnnInferenceTask {
  static constexpr double DefaultCutsMl[1][3] = {{0.5, 0.5, 0.5}};

  Produces<JetTaggingTable> taggingTable;
  jet_gnn_inference::OptionalProduces<JetGnnInfTrackTable> gnnTrackTable;
  jet_gnn_inference::OptionalProduces<JetGnnInfVertexTable> gnnVertexTable;

  // Track selection for GNN inputs
  Configurable<float> trackDcaXYMax{"trackDcaXYMax", 1, "maximum DCA xy acceptance for tracks [cm]"};
  Configurable<float> trackDcaZMax{"trackDcaZMax", 2, "maximum DCA z acceptance for tracks [cm]"};

  // ML configuration
  Configurable<int> nJetConst{"nJetConst", 100, "maximum number of jet constituents to be used for ML evaluation"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};

  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{5., 1000.}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cuts_ml::CutSmaller, cuts_ml::CutNot, cuts_ml::CutNot}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {DefaultCutsMl[0], 1, 3, {"pT bin 0"}, {"b", "c", "lf"}}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", 3, "Number of classes in ML model"};

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/c/chchoi"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"gnnbjettagger_LHC25a2b_v241.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // GNN configuration
  Configurable<float> jetpTMin{"jetpTMin", 5., "minimum jet pT"};
  Configurable<float> dbMin{"dbMin", -15., "minimum GNN Db"};
  Configurable<float> dbMax{"dbMax", 15., "maximum GNN Db"};
  Configurable<double> fC{"fC", 0.018, "Parameter f_c for D_b calculation"};
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

  o2::analysis::GNNBjetAllocator tensorAlloc;

  template <typename Jet>
  void fillHistograms(Jet const& jet, float db)
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
  }

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    if (doprocessAlgorithmGNN) {
      // analyzeJetTrackInfo4GNN supplies four jet features and eight track features.
      if (nJetFeat != 4 || nTrkFeat != 8 || nClassesMl != 3 || nJetConst <= 0) {
        LOG(fatal) << "GNN jet tagging requires nJetFeat=4, nTrkFeat=8, nClassesMl=3 and nJetConst>0.";
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
  }

  template <typename Jet, typename Tracks>
  float evaluateJet(Jet const& jet, Tracks const& tracks)
  {
    std::vector<std::vector<float>> trkFeat;
    jettaggingutilities::analyzeJetTrackInfo4GNN(jet, tracks, trkFeat, trackPtMin, trackDcaXYMax, trackDcaZMax, nJetConst);
    if (trkFeat.empty()) {
      LOGF(debug, "GNN inference, no selected tracks (%d)", jet.globalIndex());
      return -999.f;
    }

    std::vector<float> jetFeat{jet.pt(), jet.phi(), jet.eta(), jet.mass()};
    std::vector<float> feat;
    std::vector<Ort::Value> gnnInput;
    tensorAlloc.getGNNInput(jetFeat, trkFeat, feat, gnnInput);

    auto modelOutput = bMlResponse.getModelOutput(gnnInput, 0);
    float db = jettaggingutilities::getDb(modelOutput, fC);
    if (std::isnan(db)) {
      LOGF(debug, "GNN inference, Db is NaN (%d)", jet.globalIndex());
      return 999.f;
    }
    return db;
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetGnnInferenceTask, processDummy, "Dummy process", true);

  void processAlgorithmGNN(AnalysisJets const& jets, JetTracksExt const& tracks)
  {
    // Preserve a one-to-one positional join with Jets, including jets without
    // selected tracks and jets below the QA pT threshold.
    for (auto const& jet : jets) {
      float db = evaluateJet(jet, tracks);
      taggingTable(0, -1.f, db);
      fillHistograms(jet, db);
    }
  }
  PROCESS_SWITCH(JetGnnInferenceTask, processAlgorithmGNN, "Run GNN inference and fill jet tagging tables", false);
};

using JetGnnTagInferenceDataCharged = JetGnnInferenceTask<false, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, aod::ChargedJetTags>;
using JetGnnTagInferenceMcdCharged = JetGnnInferenceTask<true, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedMCDetectorLevelJetTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetGnnTagInferenceDataCharged>(
      cfgc, SetDefaultProcesses{}, TaskName{"jet-gnn-tag-inference-data-charged"}), // o2-linter: disable=name/o2-task (custom device name for templated task alias)
    adaptAnalysisTask<JetGnnTagInferenceMcdCharged>(
      cfgc, SetDefaultProcesses{}, TaskName{"jet-gnn-tag-inference-mcd-charged"})}; // o2-linter: disable=name/o2-task (custom device name for templated task alias)
}
