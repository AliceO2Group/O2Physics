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

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

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
#include <TPDGCode.h>
#include <TDatabasePDG.h>
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

// For original data loops
using CascOriginalDatas = soa::Join<aod::CascIndices, aod::CascCores>;

// For derived data analysis
using CascDerivedDatas = soa::Join<aod::CascCores, aod::CascExtras, aod::CascCollRefs>;

struct cascademlselection {
  o2::ml::OnnxModel mlModelXiMinus;
  o2::ml::OnnxModel mlModelXiPlus;
  o2::ml::OnnxModel mlModelOmegaMinus;
  o2::ml::OnnxModel mlModelOmegaPlus;

  std::map<std::string, std::string> metadata;

  Produces<aod::CascXiMLScores> xiMLSelections;    // optionally aggregate information from ML output for posterior analysis (derived data)
  Produces<aod::CascOmMLScores> omegaMLSelections; // optionally aggregate information from ML output for posterior analysis (derived data)

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  // CCDB options
  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } ccdbConfigurations;

  // Machine learning evaluation for pre-selection and corresponding information generation
  struct : ConfigurableGroup {
    // ML classifiers: master flags to populate ML Selection tables
    Configurable<bool> calculateXiMinusScores{"mlConfigurations.calculateXiMinusScores", true, "calculate XiMinus ML scores"};
    Configurable<bool> calculateXiPlusScores{"mlConfigurations.calculateXiPlusScores", true, "calculate XiPlus ML scores"};
    Configurable<bool> calculateOmegaMinusScores{"mlConfigurations.calculateOmegaMinusScores", true, "calculate OmegaMinus ML scores"};
    Configurable<bool> calculateOmegaPlusScores{"mlConfigurations.calculateOmegaPlusScores", true, "calculate OmegaPlus ML scores"};

    // ML input for ML calculation
    Configurable<std::string> modelPathCCDB{"mlConfigurations.modelPathCCDB", "", "ML Model path in CCDB"};
    Configurable<int64_t> timestampCCDB{"mlConfigurations.timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"mlConfigurations.loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"mlConfigurations.enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

    // Local paths for test purposes
    Configurable<std::string> localModelPathXiMinus{"mlConfigurations.localModelPathXiMinus", "XiMinus_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
    Configurable<std::string> localModelPathXiPlus{"mlConfigurations.localModelPathXiPlus", "XiPlus_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
    Configurable<std::string> localModelPathOmegaMinus{"mlConfigurations.localModelPathOmegaMinus", "OmegaMinus_BDTModel.onnx", "(std::string) Path to the local .onnx file."};
    Configurable<std::string> localModelPathOmegaPlus{"mlConfigurations.localModelPathOmegaPlus", "OmegaPlus_BDTModel.onnx", "(std::string) Path to the local .onnx file."};

    // Thresholds for choosing to populate V0Cores tables with pre-selections
    Configurable<float> thresholdXiMinus{"mlConfigurations.thresholdXiMinus", -1.0f, "Threshold to keep XiMinus candidates"};
    Configurable<float> thresholdXiPlus{"mlConfigurations.thresholdXiPlus", -1.0f, "Threshold to keep XiPlus candidates"};
    Configurable<float> thresholdOmegaMinus{"mlConfigurations.thresholdOmegaMinus", -1.0f, "Threshold to keep OmegaMinus candidates"};
    Configurable<float> thresholdOmegaPlus{"mlConfigurations.thresholdOmegaPlus", -1.0f, "Threshold to keep OmegaPlus candidates"};
  } mlConfigurations;

  // Axis
  // base properties
  ConfigurableAxis vertexZ{"vertexZ", {30, -15.0f, 15.0f}, ""};

  int nCandidates = 0;

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    int64_t timeStampML = 0;
    if constexpr (requires { collision.timestamp(); }) { // we are in derived data
      if (mRunNumber == collision.runNumber()) {
        return;
      }
      mRunNumber = collision.runNumber();
      timeStampML = collision.timestamp();
    }
    if constexpr (requires { collision.template bc_as<aod::BCsWithTimestamps>(); }) { // we are in original data
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      if (mRunNumber == bc.runNumber()) {
        return;
      }
      mRunNumber = bc.runNumber();
      timeStampML = bc.timestamp();
    }

    // machine learning initialization if requested
    if (mlConfigurations.calculateXiMinusScores ||
        mlConfigurations.calculateXiPlusScores ||
        mlConfigurations.calculateOmegaMinusScores ||
        mlConfigurations.calculateOmegaPlusScores) {
      if (mlConfigurations.timestampCCDB.value != -1)
        timeStampML = mlConfigurations.timestampCCDB.value;
      LoadMachines(timeStampML);
    }
  }

  // function to load models for ML-based classifiers
  void LoadMachines(int64_t timeStampML)
  {
    if (mlConfigurations.loadModelsFromCCDB) {
      ccdbApi.init(ccdbConfigurations.ccdburl);
      LOG(info) << "Fetching cascade models for timestamp: " << timeStampML;

      if (mlConfigurations.calculateXiMinusScores) {
        bool retrieveSuccess = ccdbApi.retrieveBlob(mlConfigurations.modelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathXiMinus.value);
        if (retrieveSuccess) {
          mlModelXiMinus.initModel(mlConfigurations.localModelPathXiMinus.value, mlConfigurations.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the XiMinus model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (mlConfigurations.calculateXiPlusScores) {
        bool retrieveSuccess = ccdbApi.retrieveBlob(mlConfigurations.modelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathXiPlus.value);
        if (retrieveSuccess) {
          mlModelXiPlus.initModel(mlConfigurations.localModelPathXiPlus.value, mlConfigurations.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the XiPlus model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (mlConfigurations.calculateOmegaMinusScores) {
        bool retrieveSuccess = ccdbApi.retrieveBlob(mlConfigurations.modelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathOmegaMinus.value);
        if (retrieveSuccess) {
          mlModelOmegaMinus.initModel(mlConfigurations.localModelPathOmegaMinus.value, mlConfigurations.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the OmegaMinus model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }

      if (mlConfigurations.calculateOmegaPlusScores) {
        bool retrieveSuccess = ccdbApi.retrieveBlob(mlConfigurations.modelPathCCDB, ".", metadata, timeStampML, false, mlConfigurations.localModelPathOmegaPlus.value);
        if (retrieveSuccess) {
          mlModelOmegaPlus.initModel(mlConfigurations.localModelPathOmegaPlus.value, mlConfigurations.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the OmegaPlus model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      }
    } else {
      if (mlConfigurations.calculateXiMinusScores)
        mlModelXiMinus.initModel(mlConfigurations.localModelPathXiMinus.value, mlConfigurations.enableOptimizations.value);
      if (mlConfigurations.calculateXiPlusScores)
        mlModelXiPlus.initModel(mlConfigurations.localModelPathXiPlus.value, mlConfigurations.enableOptimizations.value);
      if (mlConfigurations.calculateOmegaMinusScores)
        mlModelOmegaMinus.initModel(mlConfigurations.localModelPathOmegaMinus.value, mlConfigurations.enableOptimizations.value);
      if (mlConfigurations.calculateOmegaPlusScores)
        mlModelOmegaPlus.initModel(mlConfigurations.localModelPathOmegaPlus.value, mlConfigurations.enableOptimizations.value);
    }
    LOG(info) << "Cascade ML Models loaded.";
  }

  void init(InitContext const&)
  {
    // Histograms
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});

    ccdb->setURL(ccdbConfigurations.ccdburl);
  }

  // Process candidate and store properties in object
  template <typename TCascObject>
  void processCandidate(TCascObject const& cand)
  {
    // Select features
    // FIXME THIS NEEDS ADJUSTING
    std::vector<float> inputFeatures{0.0f, 0.0f,
                                     0.0f, 0.0f};

    // calculate scores
    if (cand.sign() < 0) {
      if (mlConfigurations.calculateXiMinusScores) {
        float* xiMinusProbability = mlModelXiMinus.evalModel(inputFeatures);
        xiMLSelections(xiMinusProbability[1]);
      } else {
        xiMLSelections(-1);
      }
      if (mlConfigurations.calculateOmegaMinusScores) {
        float* omegaMinusProbability = mlModelOmegaMinus.evalModel(inputFeatures);
        omegaMLSelections(omegaMinusProbability[1]);
      } else {
        omegaMLSelections(-1);
      }
    }
    if (cand.sign() > 0) {
      if (mlConfigurations.calculateXiPlusScores) {
        float* xiPlusProbability = mlModelXiPlus.evalModel(inputFeatures);
        xiMLSelections(xiPlusProbability[1]);
      } else {
        xiMLSelections(-1);
      }
      if (mlConfigurations.calculateOmegaPlusScores) {
        float* omegaPlusProbability = mlModelOmegaPlus.evalModel(inputFeatures);
        omegaMLSelections(omegaPlusProbability[1]);
      } else {
        omegaMLSelections(-1);
      }
    }
  }

  void processDerivedData(soa::Join<aod::StraCollisions, aod::StraStamps>::iterator const& collision, CascDerivedDatas const& cascades)
  {
    initCCDB(collision);

    histos.fill(HIST("hEventVertexZ"), collision.posZ());
    for (auto& casc : cascades) {
      nCandidates++;
      if (nCandidates % 50000 == 0) {
        LOG(info) << "Candidates processed: " << nCandidates;
      }
      processCandidate(casc);
    }
  }
  void processStandardData(aod::Collision const& collision, CascOriginalDatas const& cascades)
  {
    initCCDB(collision);

    histos.fill(HIST("hEventVertexZ"), collision.posZ());
    for (auto& casc : cascades) {
      nCandidates++;
      if (nCandidates % 50000 == 0) {
        LOG(info) << "Candidates processed: " << nCandidates;
      }
      processCandidate(casc);
    }
  }

  PROCESS_SWITCH(cascademlselection, processStandardData, "Process standard data", false);
  PROCESS_SWITCH(cascademlselection, processDerivedData, "Process derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cascademlselection>(cfgc)};
}
