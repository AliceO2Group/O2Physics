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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/ASoA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CCDB/BasicCCDBManager.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::ml;
using std::array;
using std::cout;
using std::endl;

//OnnxModel bdt;
OnnxModel lambda_bdt;
OnnxModel antilambda_bdt;
OnnxModel gamma_bdt;
OnnxModel kzeroshort_bdt;

// For original data loops
using V0OriginalDatas = soa::Join<aod::V0Indices, aod::V0Cores>;

// For derived data analysis
using V0DerivedDatas = soa::Join<aod::V0Cores, aod::V0Extras, aod::V0CollRefs>;

std::map<std::string, std::string> metadata;
std::map<std::string, std::string> headers;

struct lambdakzeromlselection{
  Produces<aod::V0GammaMLScores> gammaMLSelections; // optionally aggregate information from ML output for posterior analysis (derived data)
  Produces<aod::V0LambdaMLScores> lambdaMLSelections; // optionally aggregate information from ML output for posterior analysis (derived data)
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // ML inference
  Configurable<bool> PredictLambda{"PredictLambda", true, "Flag to enable or disable the loading of model"};
  Configurable<bool> PredictAntiLambda{"PredictAntiLambda", false, "Flag to enable or disable the loading of model"};
  Configurable<bool> PredictGamma{"PredictGamma", true, "Flag to enable or disable the loading of model"};
  Configurable<bool> PredictKZeroShort{"PredictKZeroShort", false, "Flag to enable or disable the loading of model"};
  Configurable<bool> fIsMC{"fIsMC", false, "If true, save additional MC info for analysis"};

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
      //headers = ccdbApi.retrieveHeaders(BDTPathCCDB.value, metadata, timestampCCDB.value);

      if (PredictLambda) {
        bool retrieveSuccessLambda = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathLambda.value);
        if (retrieveSuccessLambda) lambda_bdt.initModel(BDTLocalPathLambda.value, enableOptimizations.value);
        else{
        LOG(fatal) << "Error encountered while fetching/loading the Lambda model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";}
        }
      
      if (PredictAntiLambda) {
        bool retrieveSuccessAntiLambda = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathAntiLambda.value);
        if (retrieveSuccessAntiLambda) antilambda_bdt.initModel(BDTLocalPathAntiLambda.value, enableOptimizations.value);
        else{
        LOG(fatal) << "Error encountered while fetching/loading the AntiLambda model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";}
        }

      if (PredictGamma) {
        bool retrieveSuccessGamma = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathGamma.value);
        if (retrieveSuccessGamma) gamma_bdt.initModel(BDTLocalPathGamma.value, enableOptimizations.value);
        else{
          LOG(fatal) << "Error encountered while fetching/loading the Gamma model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";}
        }
      
      if (PredictKZeroShort) {
        bool retrieveSuccessKZeroShort = ccdbApi.retrieveBlob(BDTPathCCDB.value, ".", metadata, timestampCCDB.value, false, BDTLocalPathKZeroShort.value);
        if (retrieveSuccessKZeroShort) kzeroshort_bdt.initModel(BDTLocalPathKZeroShort.value, enableOptimizations.value);
        else{
          LOG(fatal) << "Error encountered while fetching/loading the KZeroShort model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";}
      }
    } 
    else {
      if (PredictLambda) lambda_bdt.initModel(BDTLocalPathLambda.value, enableOptimizations.value);
      if (PredictAntiLambda) antilambda_bdt.initModel(BDTLocalPathAntiLambda.value, enableOptimizations.value);
      if (PredictGamma) gamma_bdt.initModel(BDTLocalPathGamma.value, enableOptimizations.value);
      if (PredictKZeroShort) kzeroshort_bdt.initModel(BDTLocalPathKZeroShort.value, enableOptimizations.value); 
    }
  }

  // Process candidate and store properties in object
  template <typename TV0Object>
  void processCandidate(TV0Object const& cand)
  {    
    std::vector<float> inputFeatures{cand.pt(), static_cast<float>(cand.qtarm()), 
                                     cand.alpha(), cand.v0radius(), 
                                     cand.v0cosPA(), cand.dcaV0daughters(), 
                                     cand.dcapostopv(), cand.dcanegtopv()};

    // calculate classifier
    float* LambdaProbability = lambda_bdt.evalModel(inputFeatures);
    float* GammaProbability = gamma_bdt.evalModel(inputFeatures);
    //float* AntiLambdaProbability = antilambda_bdt.evalModel(inputFeatures); // WIP
    //float* KZeroShortProbability = kzeroshort_bdt.evalModel(inputFeatures); // WIP

    gammaMLSelections(GammaProbability[1]);
    lambdaMLSelections(LambdaProbability[1]);
  }

  void processDerivedData(aod::StraCollision const& coll, V0DerivedDatas const& v0s)
  {
    histos.fill(HIST("hEventVertexZ"), coll.posZ());
    for (auto& v0: v0s){ // looping over lambdas 
      processCandidate(v0);
    }
  }
  void processStandardData(aod::Collision const& coll, V0OriginalDatas const& v0s)
  {
    histos.fill(HIST("hEventVertexZ"), coll.posZ());
    for (auto& v0: v0s){ // looping over lambdas 
      processCandidate(v0);
    }
  }

  PROCESS_SWITCH(lambdakzeromlselection, processStandardData, "Process standard data", true);
  PROCESS_SWITCH(lambdakzeromlselection, processDerivedData, "Process derived data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdakzeromlselection>(cfgc)};
}
