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
//
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//  Lambdakzero ML selection task
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    gianni.shigeru.setoue.liveraro@cern.ch
//    romain.schotter@cern.ch
//    david.dobrigkeit.chinellato@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::ml;
using std::array;
using std::cout;
using std::endl;

// For original data loops
using V0OriginalDatas = soa::Join<aod::V0Indices, aod::V0Cores>;

// For derived data analysis
using V0DerivedDatas = soa::Join<aod::V0Cores, aod::V0Extras, aod::V0CollRefs>;

struct lambdakzeromlselection {
  o2::ml::OnnxModel lambda_bdt;
  o2::ml::OnnxModel antilambda_bdt;
  o2::ml::OnnxModel gamma_bdt;
  o2::ml::OnnxModel kzeroshort_bdt;

  std::map<std::string, std::string> metadata;

  Produces<aod::V0GammaMLScores> gammaMLSelections;           // optionally aggregate information from ML output for posterior analysis (derived data)
  Produces<aod::V0LambdaMLScores> lambdaMLSelections;         // optionally aggregate information from ML output for posterior analysis (derived data)
  Produces<aod::V0AntiLambdaMLScores> antiLambdaMLSelections; // optionally aggregate information from ML output for posterior analysis (derived data)
  Produces<aod::V0K0ShortMLScores> kzeroShortMLSelections;    // optionally aggregate information from ML output for posterior analysis (derived data)

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // ML inference
  Configurable<bool> PredictLambda{"PredictLambda", true, "Flag to enable or disable the loading of model"};
  Configurable<bool> PredictAntiLambda{"PredictAntiLambda", false, "Flag to enable or disable the loading of model"};
  Configurable<bool> PredictGamma{"PredictGamma", true, "Flag to enable or disable the loading of model"};
  Configurable<bool> PredictKZeroShort{"PredictKZeroShort", false, "Flag to enable or disable the loading of model"};
  Configurable<bool> fIsMC{"fIsMC", false, "If true, save additional MC info for analysis"};

  // Feature selection masks:

  //// Order: LambdaMass, AntiLambdaMass, GammaMass, KZeroShortMass, PT, Qt, Alpha, PosEta, NegEta, V0Eta
  Configurable<std::vector<int>> Kine_SelMap{"Kine_SelMap", std::vector<int>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "Mask to select basic kinematic features for ML Inference"};

  //// Order: Z, V0radius, PA, DCApostopv, DCAnegtopv, DCAV0daughters, DCAv0topv, PsiPair
  Configurable<std::vector<int>> Topo_SelMap{"Topo_SelMap", std::vector<int>{0, 1, 1, 1, 1, 1, 1, 0}, "Mask to select basic topological features for ML Inference"};

  //// Casting
  std::vector<int> CastKine_SelMap, CastTopo_SelMap, Feature_SelMask;

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> BDTLocalPathLambda{"BDTLocalPathLambda", "Lambda_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
  Configurable<std::string> BDTLocalPathAntiLambda{"BDTLocalPathAntiLambda", "AntiLambda_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
  Configurable<std::string> BDTLocalPathGamma{"BDTLocalPathGamma", "Gamma_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
  Configurable<std::string> BDTLocalPathKZeroShort{"BDTLocalPathKZeroShort", "KZeroShort_BDTModel.onnx", "(std::string) Path to the local .onnx file."};

  Configurable<std::string> BDTPathCCDB{"BDTPathCCDB", "Users/g/gsetouel/MLModels2", "Path on CCDB"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

  // Axis
  // base properties
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};

  int nCandidates = 0;
  void init(InitContext const&)
  {
    // Histograms
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});

    ccdb->setURL(ccdbUrl.value);
    // Retrieve the model from CCDB
    if (loadModelsFromCCDB) {
      ccdbApi.init(ccdbUrl);

      /// Fetching model for specific timestamp
      LOG(info) << "Fetching model for timestamp: " << timestampCCDB.value;

      if (PredictLambda) {
        bool retrieveSuccessLambda = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathLambda.value);
        if (retrieveSuccessLambda) {
          lambda_bdt.initModel(BDTLocalPathLambda.value, enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Lambda model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (PredictAntiLambda) {
        bool retrieveSuccessAntiLambda = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathAntiLambda.value);
        if (retrieveSuccessAntiLambda) {
          antilambda_bdt.initModel(BDTLocalPathAntiLambda.value, enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the AntiLambda model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (PredictGamma) {
        bool retrieveSuccessGamma = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathGamma.value);
        if (retrieveSuccessGamma) {
          gamma_bdt.initModel(BDTLocalPathGamma.value, enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Gamma model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (PredictKZeroShort) {
        bool retrieveSuccessKZeroShort = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathKZeroShort.value);
        if (retrieveSuccessKZeroShort) {
          kzeroshort_bdt.initModel(BDTLocalPathKZeroShort.value, enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the KZeroShort model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }
    } else {
      if (PredictLambda)
        lambda_bdt.initModel(BDTLocalPathLambda.value, enableOptimizations.value);
      if (PredictAntiLambda)
        antilambda_bdt.initModel(BDTLocalPathAntiLambda.value, enableOptimizations.value);
      if (PredictGamma)
        gamma_bdt.initModel(BDTLocalPathGamma.value, enableOptimizations.value);
      if (PredictKZeroShort)
        kzeroshort_bdt.initModel(BDTLocalPathKZeroShort.value, enableOptimizations.value);
    }

    /// Here the Configurables are passed to std::vectors
    CastKine_SelMap = (std::vector<int>)Kine_SelMap;
    CastTopo_SelMap = (std::vector<int>)Topo_SelMap;

    // Concatenate the selection masks
    Feature_SelMask.reserve(CastKine_SelMap.size() + CastTopo_SelMap.size()); // Reserve space to avoid reallocations
    Feature_SelMask.insert(Feature_SelMask.end(), CastKine_SelMap.begin(), CastKine_SelMap.end());
    Feature_SelMask.insert(Feature_SelMask.end(), CastTopo_SelMap.begin(), CastTopo_SelMap.end());
    LOG(info) << "Feature_SelMask size: " << Feature_SelMask.size();
  }

  template <typename T, typename U>
  std::vector<float> extractSelectedElements(const std::vector<T>& base_features, const std::vector<U>& Sel_mask)
  {
    std::vector<float> selected_elements;
    for (size_t i = 0; i < Sel_mask.size(); ++i) {
      if (Sel_mask[i] >= 1) { // If the mask value is true, select the corresponding element
        selected_elements.push_back(base_features[i]);
      }
    }
    return selected_elements;
  }

  // Process candidate and store properties in object
  template <typename TV0Object, typename T>
  void processCandidate(TV0Object const& cand, const std::vector<T>& Feature_SelMask)
  {
    // Select features
    std::vector<float> base_features{cand.mLambda(), cand.mAntiLambda(),
                                     cand.mGamma(), cand.mK0Short(),
                                     cand.pt(), static_cast<float>(cand.qtarm()), cand.alpha(),
                                     cand.positiveeta(), cand.negativeeta(), cand.eta(),
                                     cand.z(), cand.v0radius(), static_cast<float>(TMath::ACos(cand.v0cosPA())),
                                     cand.dcapostopv(), cand.dcanegtopv(), cand.dcaV0daughters(),
                                     cand.dcav0topv(), cand.psipair()};

    // Apply mask to select features
    std::vector<float> inputFeatures = extractSelectedElements(base_features, Feature_SelMask);

    // calculate classifier output
    if (PredictLambda) {
      float* LambdaProbability = lambda_bdt.evalModel(inputFeatures);
      lambdaMLSelections(LambdaProbability[1]);
    }
    if (PredictGamma) {
      float* GammaProbability = gamma_bdt.evalModel(inputFeatures);
      gammaMLSelections(GammaProbability[1]);
    }
    if (PredictAntiLambda) {
      float* AntiLambdaProbability = antilambda_bdt.evalModel(inputFeatures);
      antiLambdaMLSelections(AntiLambdaProbability[1]);
    }
    if (PredictKZeroShort) {
      float* KZeroShortProbability = kzeroshort_bdt.evalModel(inputFeatures);
      kzeroShortMLSelections(KZeroShortProbability[1]);
    }
  }

  void processDerivedData(aod::StraCollision const& coll, V0DerivedDatas const& v0s)
  {
    histos.fill(HIST("hEventVertexZ"), coll.posZ());
    for (auto& v0 : v0s) {
      nCandidates++;
      if (nCandidates % 50000 == 0) {
        LOG(info) << "Candidates processed: " << nCandidates;
      }
      processCandidate(v0, Feature_SelMask);
    }
  }
  void processStandardData(aod::Collision const& coll, V0OriginalDatas const& v0s)
  {
    histos.fill(HIST("hEventVertexZ"), coll.posZ());
    for (auto& v0 : v0s) {
      nCandidates++;
      if (nCandidates % 50000 == 0) {
        LOG(info) << "Candidates processed: " << nCandidates;
      }
      processCandidate(v0, Feature_SelMask);
    }
  }

  PROCESS_SWITCH(lambdakzeromlselection, processStandardData, "Process standard data", false);
  PROCESS_SWITCH(lambdakzeromlselection, processDerivedData, "Process derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdakzeromlselection>(cfgc)};
}
