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
/// \file Cascade ML selection task
/// \brief Produces ML response table for cascade selection at analysis level, either when running over original data or derived data.
/// \author Gianni Shigeru Setoue Liveraro Catalano <gianni.shigeru.setoue.liveraro@cern.ch>, UNICAMP
/// \author Romain Schotter <romain.schotter@cern.ch>, Austrian Academy of Sciences
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, Austrian Academy of Sciences
//

#include "PWGLF/DataModel/LFStrangenessMLTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/CascadeMlResponse.h"

#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::ml;

// For original data loops
using CascOriginalDatas = soa::Join<aod::CascIndices, aod::CascCores, aod::CascBBs>;

// For derived data analysis
using CascDerivedDatas = soa::Join<aod::CascCores, aod::CascExtras, aod::CascCollRefs, aod::CascBBs>;

struct cascademlselection {
  o2::analysis::CascadeMlResponse<float> mlModelXiMinus;
  o2::analysis::CascadeMlResponse<float> mlModelXiPlus;
  o2::analysis::CascadeMlResponse<float> mlModelOmegaMinus;
  o2::analysis::CascadeMlResponse<float> mlModelOmegaPlus;

  // Custom grouping
  std::vector<std::vector<int>> cascadesGrouped;

  Produces<aod::CascXiMLScores> xiMLSelections;    // optionally aggregate information from ML output for posterior analysis (derived data)
  Produces<aod::CascOmMLScores> omegaMLSelections; // optionally aggregate information from ML output for posterior analysis (derived data)

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // BDT score histograms, indexed [pT bin][class], one set per particle species
  std::vector<std::vector<std::shared_ptr<TH1>>> histScoreXiMinus;
  std::vector<std::vector<std::shared_ptr<TH1>>> histScoreXiPlus;
  std::vector<std::vector<std::shared_ptr<TH1>>> histScoreOmegaMinus;
  std::vector<std::vector<std::shared_ptr<TH1>>> histScoreOmegaPlus;

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = -1;

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdbConfigurations";
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } ccdbConfigurations;

  // Machine learning evaluation for pre-selection and corresponding information generation
  struct : ConfigurableGroup {
    std::string prefix = "mlConfigurations";
    // ML classifiers: master flags to populate ML Selection tables
    Configurable<bool> calculateXiMinusScores{"mlConfigurations.calculateXiMinusScores", true, "calculate XiMinus ML scores"};
    Configurable<bool> calculateXiPlusScores{"mlConfigurations.calculateXiPlusScores", true, "calculate XiPlus ML scores"};
    Configurable<bool> calculateOmegaMinusScores{"mlConfigurations.calculateOmegaMinusScores", true, "calculate OmegaMinus ML scores"};
    Configurable<bool> calculateOmegaPlusScores{"mlConfigurations.calculateOmegaPlusScores", true, "calculate OmegaPlus ML scores"};

    // List and order of input features fed to the ONNX models; any subset/order of the names
    // registered in CascMlResponse::setAvailableInputFeatures can be used here
    Configurable<std::vector<std::string>> namesInputFeatures{"mlConfigurations.namesInputFeatures", std::vector<std::string>{"cascradius", "v0radius", "casccosPA", "v0cosPA", "dcapostopv", "dcanegtopv", "dcabachtopv", "dcacascdaughters", "dcaV0daughters", "dcav0topv", "bachBaryonCosPA", "bachBaryonDCAxyToPV"}, "Names (and order) of the input features to be used in the ML models"};

    // ML input for ML calculation: one model (ONNX file / CCDB path) per pT bin
    Configurable<std::vector<std::string>> modelPathsCCDBXiMinus{"mlConfigurations.modelPathsCCDBXiMinus", std::vector<std::string>{""}, "ML Model paths in CCDB for Xi-. One per pT bin."};
    Configurable<std::vector<std::string>> modelPathsCCDBXiPlus{"mlConfigurations.modelPathsCCDBXiPlus", std::vector<std::string>{""}, "ML Model paths in CCDB for Xi+. One per pT bin."};
    Configurable<std::vector<std::string>> modelPathsCCDBOmegaMinus{"mlConfigurations.modelPathsCCDBOmegaMinus", std::vector<std::string>{""}, "ML Model paths in CCDB for Omega-. One per pT bin."};
    Configurable<std::vector<std::string>> modelPathsCCDBOmegaPlus{"mlConfigurations.modelPathsCCDBOmegaPlus", std::vector<std::string>{""}, "ML Model paths in CCDB for Omega+. One per pT bin."};
    Configurable<int64_t> timestampCCDB{"mlConfigurations.timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"mlConfigurations.loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"mlConfigurations.enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

    // Local/cvmfs paths (also used as CCDB download destination filenames), one per pT bin
    Configurable<std::vector<std::string>> onnxFileNamesXiMinus{"mlConfigurations.onnxFileNamesXiMinus", std::vector<std::string>{"XiMinus_BDTModel.onnx"}, "(std::string) Paths to the local .onnx file. One per pT bin."};
    Configurable<std::vector<std::string>> onnxFileNamesXiPlus{"mlConfigurations.onnxFileNamesXiPlus", std::vector<std::string>{"XiPlus_BDTModel.onnx"}, "(std::string) Paths to the local .onnx file. One per pT bin."};
    Configurable<std::vector<std::string>> onnxFileNamesOmegaMinus{"mlConfigurations.onnxFileNamesOmegaMinus", std::vector<std::string>{"OmegaMinus_BDTModel.onnx"}, "(std::string) Paths to the local .onnx file. One per pT bin."};
    Configurable<std::vector<std::string>> onnxFileNamesOmegaPlus{"mlConfigurations.onnxFileNamesOmegaPlus", std::vector<std::string>{"OmegaPlus_BDTModel.onnx"}, "(std::string) Paths to the local .onnx file. One per pT bin."};

    // Binning
    Configurable<std::vector<double>> binsPtXiMinus{"binsPtXiMinus", std::vector<double>{0., 10.}, "pT bin limits for ML application for Xi-"};
    Configurable<std::vector<double>> binsPtXiPlus{"binsPtXiPlus", std::vector<double>{0., 10.}, "pT bin limits for ML application for Xi+"};
    Configurable<std::vector<double>> binsPtOmegaMinus{"binsPtOmegaMinus", std::vector<double>{0., 10.}, "pT bin limits for ML application for Omega-"};
    Configurable<std::vector<double>> binsPtOmegaPlus{"binsPtOmegaPlus", std::vector<double>{0., 10.}, "pT bin limits for ML application for Omega+"};

    // Number of classes in the ML models. Default is 2 (signal and background)
    Configurable<int> nClassesMlXiMinus{"nClassesMlXiMinus", 2, "Number of classes in ML model for Xi-"};
    Configurable<int> nClassesMlXiPlus{"nClassesMlXiPlus", 2, "Number of classes in ML model for Xi+"};
    Configurable<int> nClassesMlOmegaMinus{"nClassesMlOmegaMinus", 2, "Number of classes in ML model for Omega-"};
    Configurable<int> nClassesMlOmegaPlus{"nClassesMlOmegaPlus", 2, "Number of classes in ML model for Omega+"};
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
      loadMachines(timeStampML);
    }
  }

  // function to load models for ML-based classifiers
  void loadMachines(int64_t timeStampML)
  {
    auto loadModel = [&](bool doCalculate, o2::analysis::CascadeMlResponse<float>& model, std::vector<std::string> const& onnxFileNames, std::vector<std::string> const& pathsCCDB) {
      if (!doCalculate) {
        return;
      }
      if (mlConfigurations.loadModelsFromCCDB) {
        model.setModelPathsCCDB(onnxFileNames, ccdbApi, pathsCCDB, timeStampML);
      } else {
        model.setModelPathsLocal(onnxFileNames);
      }
      model.init(mlConfigurations.enableOptimizations.value);
    };

    if (mlConfigurations.loadModelsFromCCDB) {
      ccdbApi.init(ccdbConfigurations.ccdburl);
      LOG(info) << "Fetching cascade models for timestamp: " << timeStampML;
    }

    loadModel(mlConfigurations.calculateXiMinusScores, mlModelXiMinus, mlConfigurations.onnxFileNamesXiMinus, mlConfigurations.modelPathsCCDBXiMinus);
    loadModel(mlConfigurations.calculateXiPlusScores, mlModelXiPlus, mlConfigurations.onnxFileNamesXiPlus, mlConfigurations.modelPathsCCDBXiPlus);
    loadModel(mlConfigurations.calculateOmegaMinusScores, mlModelOmegaMinus, mlConfigurations.onnxFileNamesOmegaMinus, mlConfigurations.modelPathsCCDBOmegaMinus);
    loadModel(mlConfigurations.calculateOmegaPlusScores, mlModelOmegaPlus, mlConfigurations.onnxFileNamesOmegaPlus, mlConfigurations.modelPathsCCDBOmegaPlus);

    LOG(info) << "Cascade ML Models loaded.";
  }

  void init(InitContext const&)
  {
    // Histograms
    histos.add("hEventVertexZ", "hEventVertexZ", kTH1F, {vertexZ});

    ccdb->setURL(ccdbConfigurations.ccdburl);

    // builds a shape-valid but functionally unused cuts array: this task only stores
    // raw ML scores (thresholds are applied downstream), so no cut direction is used
    auto dummyCuts = [](int nBins, int nClasses) {
      std::vector<double> zeros(static_cast<std::size_t>(nBins) * static_cast<std::size_t>(nClasses), 0.);
      return LabeledArray<double>(zeros.data(), nBins, nClasses);
    };

    auto configureModel = [&](o2::analysis::CascadeMlResponse<float>& model, std::vector<double> const& binsPt, int nClasses) {
      int nBins = static_cast<int>(binsPt.size()) - 1;
      model.configure(binsPt, dummyCuts(nBins, nClasses), std::vector<int>(nClasses, o2::cuts_ml::CutNot), nClasses);
      model.cacheInputFeaturesIndices(mlConfigurations.namesInputFeatures);
    };

    auto bookScoreHistos = [&](std::vector<std::vector<std::shared_ptr<TH1>>>& target, std::vector<double> const& binsPt, int nClasses, std::string const& particleName) {
      int nBins = static_cast<int>(binsPt.size()) - 1;
      target.resize(nBins);
      for (int iBin = 0; iBin < nBins; iBin++) {
        target[iBin].resize(nClasses);
        for (int iClass = 0; iClass < nClasses; iClass++) {
          target[iBin][iClass] = histos.add<TH1>(Form("BDTScore/%s/pTbin%d/class%d", particleName.c_str(), iBin, iClass),
                                                 Form("%s BDT score, %.2f #leq p_{T} < %.2f GeV/c, class %d;BDT score;entries", particleName.c_str(), binsPt[iBin], binsPt[iBin + 1], iClass),
                                                 kTH1F, {{100, 0., 1.}});
        }
      }
    };

    if (mlConfigurations.calculateXiMinusScores) {
      configureModel(mlModelXiMinus, mlConfigurations.binsPtXiMinus, mlConfigurations.nClassesMlXiMinus);
      bookScoreHistos(histScoreXiMinus, mlConfigurations.binsPtXiMinus, mlConfigurations.nClassesMlXiMinus, "XiMinus");
    }
    if (mlConfigurations.calculateXiPlusScores) {
      configureModel(mlModelXiPlus, mlConfigurations.binsPtXiPlus, mlConfigurations.nClassesMlXiPlus);
      bookScoreHistos(histScoreXiPlus, mlConfigurations.binsPtXiPlus, mlConfigurations.nClassesMlXiPlus, "XiPlus");
    }
    if (mlConfigurations.calculateOmegaMinusScores) {
      configureModel(mlModelOmegaMinus, mlConfigurations.binsPtOmegaMinus, mlConfigurations.nClassesMlOmegaMinus);
      bookScoreHistos(histScoreOmegaMinus, mlConfigurations.binsPtOmegaMinus, mlConfigurations.nClassesMlOmegaMinus, "OmegaMinus");
    }
    if (mlConfigurations.calculateOmegaPlusScores) {
      configureModel(mlModelOmegaPlus, mlConfigurations.binsPtOmegaPlus, mlConfigurations.nClassesMlOmegaPlus);
      bookScoreHistos(histScoreOmegaPlus, mlConfigurations.binsPtOmegaPlus, mlConfigurations.nClassesMlOmegaPlus, "OmegaPlus");
    }
  }

  // Finds the pT bin matching the MlResponse internal convention (upper_bound over bin edges);
  // returns -1 if pt falls outside the configured range, since evaluating the model in that
  // case would otherwise trigger a LOG(fatal) inside MlResponse::getModelOutput
  int findPtBin(std::vector<double> const& binsPt, float pt)
  {
    if (pt < binsPt.front() || pt >= binsPt.back()) {
      return -1;
    }
    return static_cast<int>(std::distance(binsPt.begin(), std::upper_bound(binsPt.begin(), binsPt.end(), pt))) - 1;
  }

  // Evaluates one particle hypothesis' model for a candidate, fills the per-bin/per-class
  // score histograms, and returns the signal (class 1) score to be stored in the output table
  template <typename TCascObject, typename TCollision>
  std::vector<float> evaluateModel(o2::analysis::CascadeMlResponse<float>& model, std::vector<std::vector<std::shared_ptr<TH1>>>& scoreHistos, std::vector<double> const& binsPt, bool doCalculate, TCascObject const& casc, float pt, TCollision const& coll)
  {
    if (!doCalculate) {
      return {};
    }
    int iBin = findPtBin(binsPt, pt);
    if (iBin < 0) {
      return {};
    }

    auto inputFeatures = model.getInputFeatures(casc, coll);
    std::vector<float> output;
    model.isSelectedMl(inputFeatures, pt, output);

    for (std::size_t iClass = 0; iClass < output.size() && iClass < scoreHistos[iBin].size(); iClass++) {
      scoreHistos[iBin][iClass]->Fill(output[iClass]);
    }
    return output;
  }

  // Process candidate and store properties in object
  template <typename TCascObject, typename TCollision>
  void processCandidate(TCascObject const& casc, float pt, TCollision const& coll)
  {
    // calculate scores (cascades only ever carry sign +1 or -1)
    if (casc.sign() < 0) {
      xiMLSelections(evaluateModel(mlModelXiMinus, histScoreXiMinus, mlConfigurations.binsPtXiMinus, mlConfigurations.calculateXiMinusScores, casc, pt, coll));
      omegaMLSelections(evaluateModel(mlModelOmegaMinus, histScoreOmegaMinus, mlConfigurations.binsPtOmegaMinus, mlConfigurations.calculateOmegaMinusScores, casc, pt, coll));
    } else {
      xiMLSelections(evaluateModel(mlModelXiPlus, histScoreXiPlus, mlConfigurations.binsPtXiPlus, mlConfigurations.calculateXiPlusScores, casc, pt, coll));
      omegaMLSelections(evaluateModel(mlModelOmegaPlus, histScoreOmegaPlus, mlConfigurations.binsPtOmegaPlus, mlConfigurations.calculateOmegaPlusScores, casc, pt, coll));
    }
  }

  void processDerivedData(soa::Join<aod::StraCollisions, aod::StraStamps> const& collisions, CascDerivedDatas const& cascades)
  {
    // Custom grouping
    cascadesGrouped.clear();
    cascadesGrouped.resize(collisions.size());

    for (const auto& cascade : cascades) {
      cascadesGrouped[cascade.straCollisionId()].push_back(cascade.globalIndex());
    }

    for (const auto& collision : collisions) {
      initCCDB(collision);

      histos.fill(HIST("hEventVertexZ"), collision.posZ());
      for (std::size_t i = 0; i < cascadesGrouped[collision.globalIndex()].size(); i++) {
        auto casc = cascades.rawIteratorAt(cascadesGrouped[collision.globalIndex()][i]);
        processCandidate(casc, casc.pt(), collision);
      }
    }
  }
  void processStandardData(aod::Collisions const& collisions, CascOriginalDatas const& cascades)
  {
    // Custom grouping
    cascadesGrouped.clear();
    cascadesGrouped.resize(collisions.size());

    for (const auto& cascade : cascades) {
      cascadesGrouped[cascade.collisionId()].push_back(cascade.globalIndex());
    }

    for (const auto& collision : collisions) {
      initCCDB(collision);

      histos.fill(HIST("hEventVertexZ"), collision.posZ());
      for (std::size_t i = 0; i < cascadesGrouped[collision.globalIndex()].size(); i++) {
        auto casc = cascades.rawIteratorAt(cascadesGrouped[collision.globalIndex()][i]);
        processCandidate(casc, casc.pt(), collision);
      }
    }
  }

  PROCESS_SWITCH(cascademlselection, processStandardData, "Process standard data", false);
  PROCESS_SWITCH(cascademlselection, processDerivedData, "Process derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cascademlselection>(cfgc)};
}
