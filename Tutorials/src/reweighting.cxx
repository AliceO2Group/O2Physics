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
///
/// \brief A task to apply ML reweighting
/// \author
/// \since

/*
 * command to run:
 * o2-analysis-track-propagation |\
 * o2-analysis-timestamp |\
 * o2-analysis-multiplicity-table |\
 * o2-analysis-event-selection |\
 * o2-analysistutorial-reweighting
 */

/// Uses model produced in train_model.pynb
/// Currently only works locally because model file needs to be accessed

// This workflow is used to create a flat tree for model training
// Use o2-aod-merger to combine dataframes in output AnalysisResults_trees.root

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/Multiplicity.h"
#include "TrainingTree.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec MultAxis = {301, -0.5, 300.5};
AxisSpec ForwardMultAxis = {3001, -0.5, 3000.5};
AxisSpec PtAxis = {2401, -0.005, 24.005};

namespace o2::aod
{
namespace weights
{
DECLARE_SOA_COLUMN(Weight, weight, float);
} // namespace weights
DECLARE_SOA_TABLE(Weights, "AOD", "WGHTS",
                  weights::Weight);
} // namespace o2::aod

/// function to extract collision properties for feeding the model
/// note that this version returns a fixed-size array which is not necessary and
/// for real-life applications the input array should be created dynamically
/// according to model input requirements
template <typename C, typename T>
std::array<float, 6> collect(C const& collision, T const& tracks)
{
  return {collision.posZ(), collision.posX(), collision.posY(), analysis::meanPt(tracks), static_cast<float>(tracks.size()), collision.multFT0M()};
}

/// Task to process collisions and create weighting information by applying the
/// model
struct CreateWeights {
  Configurable<float> centralEtaCut{"centralEtaCut", 0.8, "central eta limit"};
  /// Currently local and cvmfs files can be used
  /// Will be possible to load models from CCDB in the future
  Configurable<std::string> onnxModel{"onnxModel", "test_xgboost.onnx", "Path to model"};
  /// onnx runtime environment
  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "ml-reweight"};
  Produces<aod::Weights> w;

  Filter centralTracks = nabs(aod::track::eta) < centralEtaCut;

  /// onnx runtime session handle
  std::shared_ptr<Ort::Experimental::Session> onnxSession = nullptr;
  /// onnx runtime session options
  Ort::SessionOptions sessionOptions;
  /// input vectore
  std::vector<Ort::Value> inputML;
  /// input shapes characterisc (needs to be adjusted)
  std::vector<std::vector<int64_t>> inputShapes;

  void init(InitContext&)
  {
    auto path = (std::string)onnxModel;
    /// create session
    onnxSession = std::make_shared<Ort::Experimental::Session>(env, path, sessionOptions);
    /// adjust input shape to use row-by-row model application
    inputShapes = onnxSession->GetInputShapes();
    if (inputShapes[0][0] < 0) {
      LOG(warning) << "Model with negative input shape, setting it to 1.";
      inputShapes[0][0] = 1;
    }
  }

  void process(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, soa::Filtered<aod::Tracks> const& tracks)
  {
    /// get the input variables
    auto features = collect(collision, tracks);
    /// add an entry in input vector
    inputML.push_back(Ort::Experimental::Value::CreateTensor<float>(features.data(), features.size(), inputShapes[0]));
    /// run inference
    auto result = onnxSession->Run(onnxSession->GetInputNames(), inputML, onnxSession->GetOutputNames());
    /// extract scores
    auto scores = result[1].GetTensorMutableData<float>();
    LOGP(info, "Col {}: scores ({}, {})", collision.globalIndex(), scores[0], scores[1]);
    /// write model output (BDT score)
    w(scores[1]);
    /// clear input vector for next iteration
    inputML.clear();
  }
};

struct ConsumeWeights {
  Configurable<float> centralEtaCut{"centralEtaCut", 0.8, "central eta limit"};
  HistogramRegistry registry{
    "registry", //
    {
      {"Events/Vertex/Z", "; Z_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}},                            //
      {"Events/Vertex/X", "; X_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}},                            //
      {"Events/Vertex/Y", "; Y_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}},                            //
      {"Events/Multiplicity/Central", "; N_{trk}", {HistType::kTH1F, {MultAxis}}},                  //
      {"Events/Multiplicity/Forward", "; N (a.u.)", {HistType::kTH1F, {ForwardMultAxis}}},          //
      {"Events/Averages/Pt", "; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}},                       //
      {"Tracks/Pt", "; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}},                                //
      {"Weighted/Events/Vertex/Z", "; Z_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}},                   //
      {"Weighted/Events/Vertex/X", "; X_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}},                   //
      {"Weighted/Events/Vertex/Y", "; Y_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}},                   //
      {"Weighted/Events/Multiplicity/Central", "; N_{trk}", {HistType::kTH1F, {MultAxis}}},         //
      {"Weighted/Events/Multiplicity/Forward", "; N (a.u.)", {HistType::kTH1F, {ForwardMultAxis}}}, //
      {"Weighted/Events/Averages/Pt", "; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}},              //
      {"Weighted/Tracks/Pt", "; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}}                        //
    },                                                                                              //
  };

  Filter centralTracks = nabs(aod::track::eta) < centralEtaCut;

  void processWeighted(soa::Join<aod::Collisions, aod::Mults, aod::Weights>::iterator const& collision, soa::Filtered<aod::Tracks> const& tracks)
  {
    /// fill histograms with using BDT scores produced by previous task as weights

    for (auto& track : tracks) {
      if (isfinite(track.pt())) {
        registry.fill(HIST("Weighted/Tracks/Pt"), track.pt(), collision.weight());
      }
    }
    auto output = collect(collision, tracks);
    registry.fill(HIST("Weighted/Events/Vertex/Z"), output[0], collision.weight());
    registry.fill(HIST("Weighted/Events/Vertex/X"), output[1], collision.weight());
    registry.fill(HIST("Weighted/Events/Vertex/Y"), output[2], collision.weight());
    registry.fill(HIST("Weighted/Events/Multiplicity/Central"), output[4], collision.weight());
    registry.fill(HIST("Weighted/Events/Multiplicity/Forward"), output[5], collision.weight());
    registry.fill(HIST("Weighted/Events/Averages/Pt"), output[3], collision.weight());
  }

  PROCESS_SWITCH(ConsumeWeights, processWeighted, "Use weighted events", true);

  void processNormal(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, soa::Filtered<aod::Tracks> const& tracks)
  {
    /// fill histograms without weights

    for (auto& track : tracks) {
      if (isfinite(track.pt())) {
        registry.fill(HIST("Tracks/Pt"), track.pt());
      }
    }
    auto output = collect(collision, tracks);
    registry.fill(HIST("Events/Vertex/Z"), output[0]);
    registry.fill(HIST("Events/Vertex/X"), output[1]);
    registry.fill(HIST("Events/Vertex/Y"), output[2]);
    registry.fill(HIST("Events/Multiplicity/Central"), output[4]);
    registry.fill(HIST("Events/Multiplicity/Forward"), output[5]);
    registry.fill(HIST("Events/Averages/Pt"), output[3]);
  }

  PROCESS_SWITCH(ConsumeWeights, processNormal, "Use original events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CreateWeights>(cfgc), adaptAnalysisTask<ConsumeWeights>(cfgc)};
}
