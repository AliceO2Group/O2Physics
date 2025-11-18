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

/// \file jetTaggerHF.cxx
/// \brief Task to produce a table joinable to the jet tables for hf jet tagging
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>
/// \author Hanseo Park <hanseo.park@cern.ch>
/// \author Hadi Hassan <hadi.hassan@cern.ch>, University of Jyväskylä

#include "MlResponse.h"

#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/MlResponseHfTagging.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include <CCDB/CcdbApi.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TH1.h>

#include <onnxruntime_cxx_api.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

template <bool isMCD, typename JetTable, typename SVIndicesTable, typename SVTable, typename JetTaggingTable>
struct JetTaggerHFTask {
  static constexpr double DefaultCutsMl[1][2] = {{0.5, 0.5}};

  Produces<JetTaggingTable> taggingTable;

  // configuration topological cut for track and sv
  Configurable<float> trackDcaXYMax{"trackDcaXYMax", 1, "minimum DCA xy acceptance for tracks [cm]"};
  Configurable<float> trackDcaZMax{"trackDcaZMax", 2, "minimum DCA z acceptance for tracks [cm]"};
  Configurable<float> prongsigmaLxyMax{"prongsigmaLxyMax", 100, "maximum sigma of decay length of prongs on xy plane"};
  Configurable<float> prongsigmaLxyzMax{"prongsigmaLxyzMax", 100, "maximum sigma of decay length of prongs on xyz plane"};
  Configurable<float> prongIPxyMin{"prongIPxyMin", 0.008, "maximum impact paramter of prongs on xy plane [cm]"};
  Configurable<float> prongIPxyMax{"prongIPxyMax", 1, "minimum impact parmeter of prongs on xy plane [cm]"};
  Configurable<float> prongChi2PCAMin{"prongChi2PCAMin", 4, "minimum Chi2 PCA of decay length of prongs"};
  Configurable<float> prongChi2PCAMax{"prongChi2PCAMax", 100, "maximum Chi2 PCA of decay length of prongs"};
  Configurable<float> svDispersionMax{"svDispersionMax", 1, "maximum dispersion of sv"};

  // configuration about IP method
  Configurable<bool> useJetProb{"useJetProb", false, "fill table for track counting algorithm"};
  Configurable<bool> trackProbQA{"trackProbQA", false, "fill track probability histograms separately for geometric positive and negative tracks for QA"};
  Configurable<int> numCount{"numCount", 3, "number of track counting"};
  Configurable<int> resoFuncMatching{"resoFuncMatching", 0, "matching parameters of resolution function as MC samble (0: custom, 1: custom & inc, 2: MB, 3: MB & inc, 4: JJ, 5: JJ & inc)"};
  Configurable<std::vector<std::string>> pathsCCDBforIPDataparameter{"pathsCCDBforIPDataparameter", std::vector<std::string>{"Users/l/leehy/LHC24g4/f_inclusive_0"}, "Paths for fitting parameters of resolution functions of data for IP method on CCDB"};
  Configurable<std::vector<std::string>> pathsCCDBforIPIncparameter{"pathsCCDBforIPIncparameter", std::vector<std::string>{"Users/l/leehy/LHC24g4/f_inclusive_0"}, "Paths for fitting parameters of resolution functions of inclusive for IP method on CCDB"};
  Configurable<std::vector<std::string>> pathsCCDBforIPBeautyparameter{"pathsCCDBforIPBeautyparameter", std::vector<std::string>{"Users/l/leehy/LHC24g4/f_inclusive_0"}, "Paths for fitting parameters of resolution functions of beauty for IP method on CCDB"};
  Configurable<std::vector<std::string>> pathsCCDBforIPCharmparameter{"pathsCCDBforIPCharmparameter", std::vector<std::string>{"Users/l/leehy/LHC24g4/f_inclusive_0"}, "Paths for fitting parameters of resolution functions of charm for IP method on CCDB"};
  Configurable<std::vector<std::string>> pathsCCDBforIPLfparameter{"pathsCCDBforIPLfparameter", std::vector<std::string>{"Users/l/leehy/LHC24g4/f_inclusive_0"}, "Paths for fitting parameters of resolution functions of light flavour for IP method on CCDB"};
  Configurable<bool> usepTcategorize{"usepTcategorize", false, "p_T categorize TF1 function with Inclusive jet"};
  Configurable<std::vector<float>> paramsResoFuncData{"paramsResoFuncData", std::vector<float>{-1.0}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7))"};
  Configurable<std::vector<float>> paramsResoFuncIncJetMC{"paramsResoFuncIncJetMC", std::vector<float>{-1.0}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncCharmJetMC{"paramsResoFuncCharmJetMC", std::vector<float>{-1.0}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncBeautyJetMC{"paramsResoFuncBeautyJetMC", std::vector<float>{-1.0}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncLfJetMC{"paramsResoFuncLfJetMC", std::vector<float>{-1.0}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<float> minSignImpXYSig{"minSignImpXYSig", -40.0, "minimum of signed impact parameter significance"};
  Configurable<float> tagPointForIP{"tagPointForIP", 2.5, "tagging working point for IP"};
  Configurable<float> tagPointForIPxyz{"tagPointForIPxyz", 2.5, "tagging working point for IP xyz"};
  // configuration about SV method
  Configurable<float> tagPointForSV{"tagPointForSV", 40, "tagging working point for SV"};
  Configurable<float> tagPointForSVxyz{"tagPointForSVxyz", 40, "tagging working point for SV xyz"};

  // ML configuration
  Configurable<int> nJetConst{"nJetConst", 10, "maximum number of jet consistuents to be used for ML evaluation"};
  Configurable<float> trackPtMin{"trackPtMin", 0.5, "minimum track pT"};
  Configurable<float> svPtMin{"svPtMin", 0.5, "minimum SV pT"};

  Configurable<float> svReductionFactor{"svReductionFactor", 1.0, "factor for how many SVs to keep"};

  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{5., 1000.}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{cuts_ml::CutSmaller, cuts_ml::CutNot}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {DefaultCutsMl[0], 1, 2, {"pT bin 0"}, {"score for default b-jet tagging", "uncer 1"}}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", 2, "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  Configurable<bool> useDb{"useDb", false, "Flag to use DB for ML model instead of the score"};

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/h/hahassan"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ML_bjets/Models/LHC24g4_70_200/model.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // GNN configuration
  Configurable<float> jetpTMin{"jetpTMin", 5., "minimum jet pT"};
  Configurable<float> dbMin{"dbMin", -10., "minimum GNN Db"};
  Configurable<float> dbMax{"dbMax", 20., "maximum GNN Db"};
  Configurable<double> fC{"fC", 0.018, "Parameter f_c for D_b calculation"};
  Configurable<int64_t> nJetFeat{"nJetFeat", 4, "Number of jet GNN input features"};
  Configurable<int64_t> nTrkFeat{"nTrkFeat", 13, "Number of track GNN input features"};
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

  // axis spec
  ConfigurableAxis binTrackProbability{"binTrackProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};

  o2::analysis::MlResponseHfTagging<float> bMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  using JetTracksExt = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using OriginalTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov, aod::TracksExtra>;

  bool useResoFuncFromIncJet = false;
  int maxOrder = -1;
  int resoFuncMatch = 0;

  std::unique_ptr<TF1> fSignImpXYSigData = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigCharmJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigBeautyJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigLfJetMC = nullptr;

  std::vector<std::vector<float>> vecParamsDataJetCCDB;
  std::vector<std::vector<float>> vecParamsIncJetMcCCDB;
  std::vector<std::vector<float>> vecParamsBeautyJetMcCCDB;
  std::vector<std::vector<float>> vecParamsCharmJetMcCCDB;
  std::vector<std::vector<float>> vecParamsLfJetMcCCDB;

  std::vector<std::unique_ptr<TF1>> vecfSignImpXYSigDataJetCCDB;
  std::vector<std::unique_ptr<TF1>> vecfSignImpXYSigIncJetMcCCDB;
  std::vector<std::unique_ptr<TF1>> vecfSignImpXYSigCharmJetMcCCDB;
  std::vector<std::unique_ptr<TF1>> vecfSignImpXYSigBeautyJetMcCCDB;
  std::vector<std::unique_ptr<TF1>> vecfSignImpXYSigLfJetMcCCDB;

  std::vector<uint16_t> decisionNonML;
  std::vector<float> scoreML;

  o2::analysis::GNNBjetAllocator tensorAlloc;

  template <typename T, typename U>
  float calculateJetProbability(int origin, T const& jet, U const& tracks, bool const& isMC = false)
  {
    float jetProb = -1.0;
    if (!isMC) {
      if (usepTcategorize) {
        jetProb = jettaggingutilities::getJetProbability(vecfSignImpXYSigDataJetCCDB, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
      } else {
        jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigData, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
      }
    } else {
      if (useResoFuncFromIncJet) {
        if (usepTcategorize) {
          jetProb = jettaggingutilities::getJetProbability(vecfSignImpXYSigIncJetMcCCDB, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
        } else {
          jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigIncJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
        }
      } else {
        if (origin == JetTaggingSpecies::charm) {
          if (usepTcategorize) {
            jetProb = jettaggingutilities::getJetProbability(vecfSignImpXYSigCharmJetMcCCDB, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
          } else {
            jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigCharmJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
          }
        } else if (origin == JetTaggingSpecies::beauty) {
          if (usepTcategorize) {
            jetProb = jettaggingutilities::getJetProbability(vecfSignImpXYSigBeautyJetMcCCDB, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
          } else {
            jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigBeautyJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
          }
        } else {
          if (usepTcategorize) {
            jetProb = jettaggingutilities::getJetProbability(vecfSignImpXYSigLfJetMcCCDB, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
          } else {
            jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigLfJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
          }
        }
      }
    }
    return jetProb;
  }

  template <typename T, typename U>
  void evaluateTrackProbQA(int origin, T const& jet, U const& /*tracks*/, bool const& isMC = true)
  {
    for (const auto& track : jet.template tracks_as<U>()) {
      if (!jettaggingutilities::trackAcceptanceWithDca(track, trackDcaXYMax, trackDcaZMax))
        continue;
      auto geoSign = jettaggingutilities::getGeoSign(jet, track);
      float probTrack = -1;
      if (!isMC) {
        probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigData, track, minSignImpXYSig);
        if (geoSign > 0)
          registry.fill(HIST("h_pos_track_probability"), probTrack);
        else
          registry.fill(HIST("h_neg_track_probability"), probTrack);
      } else {
        if (useResoFuncFromIncJet) {
          probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigIncJetMC, track, minSignImpXYSig);
        } else {
          if (origin == JetTaggingSpecies::charm) {
            probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigCharmJetMC, track, minSignImpXYSig);
          }
          if (origin == JetTaggingSpecies::beauty) {
            probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigBeautyJetMC, track, minSignImpXYSig);
          }
          if (origin == JetTaggingSpecies::lightflavour) {
            probTrack = jettaggingutilities::getTrackProbability(fSignImpXYSigLfJetMC, track, minSignImpXYSig);
          }
        }
        if (geoSign > 0)
          registry.fill(HIST("h2_pos_track_probability_flavour"), probTrack, origin);
        else
          registry.fill(HIST("h2_neg_track_probability_flavour"), probTrack, origin);
      }
    }
  }

  template <bool isMC, typename T, typename U>
  void fillTables(T const& jet, U const& tracks)
  {
    int origin = -1;
    float jetProb = -1.0;
    if constexpr (isMC) {
      origin = jet.origin();
    }
    if (useJetProb) {
      if constexpr (isMC) {
        jetProb = calculateJetProbability(origin, jet, tracks, isMC);
        if (trackProbQA) {
          evaluateTrackProbQA(origin, jet, tracks, isMC);
        }
      } else {
        jetProb = calculateJetProbability(0, jet, tracks, isMC);
        if (trackProbQA) {
          evaluateTrackProbQA(0, jet, tracks, isMC);
        }
      }
    }
    if (doprocessAlgorithmGNN) {
      if (jet.pt() >= jetpTMin) {
        float dbRange;
        if (scoreML[jet.globalIndex()] < dbMin) {
          dbRange = 0.5; // underflow
        } else if (scoreML[jet.globalIndex()] < dbMax) {
          dbRange = 1.5; // in range
        } else {
          dbRange = 2.5; // overflow
        }
        registry.fill(HIST("h2_count_db"), 3.5, dbRange); // incl jet
        if constexpr (isMC) {
          switch (origin) {
            case 2:
              registry.fill(HIST("h_db_b"), scoreML[jet.globalIndex()]);
              registry.fill(HIST("h2_count_db"), 0.5, dbRange); // b-jet
              break;
            case 1:
              registry.fill(HIST("h_db_c"), scoreML[jet.globalIndex()]);
              registry.fill(HIST("h2_count_db"), 1.5, dbRange); // c-jet
              break;
            case 0:
            case 3:
              registry.fill(HIST("h_db_lf"), scoreML[jet.globalIndex()]);
              registry.fill(HIST("h2_count_db"), 2.5, dbRange); // lf-jet
              break;
            default:
              LOGF(debug, "doprocessAlgorithmGNN, Unexpected origin value: %d (%d)", origin, jet.globalIndex());
          }
        }
        registry.fill(HIST("h2_pt_db"), jet.pt(), scoreML[jet.globalIndex()]);
      }
    }
    taggingTable(decisionNonML[jet.globalIndex()], jetProb, scoreML[jet.globalIndex()]);
  }

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    std::vector<float> vecParamsData;
    std::vector<float> vecParamsIncJetMC;
    std::vector<float> vecParamsCharmJetMC;
    std::vector<float> vecParamsBeautyJetMC;
    std::vector<float> vecParamsLfJetMC;

    std::vector<TF1*> resoFuncDataCCDB;
    std::vector<TF1*> resoFuncIncCCDB;
    std::vector<TF1*> resoFuncBeautyCCDB;
    std::vector<TF1*> resoFuncCharmCCDB;
    std::vector<TF1*> resoFuncLfCCDB;

    ccdbApi.init(ccdbUrl);

    std::map<std::string, std::string> metadata;
    resoFuncMatch = resoFuncMatching;

    const int mIPmethodResolutionFunctionSize = 7;

    auto loadCCDBforIP = [&](const std::vector<std::string>& paths, std::vector<TF1*>& targetVec, const std::string& name) {
      if (paths.size() != mIPmethodResolutionFunctionSize) {
        usepTcategorize.value = false;
        LOG(info) << name << " does not have 7 entries. Disabling pT categorization (usepTcategorize = false).";
        resoFuncMatch = 0;
        return;
      }
      for (int i = 0; i < mIPmethodResolutionFunctionSize; i++) {
        targetVec.push_back(ccdbApi.retrieveFromTFileAny<TF1>(paths[i], metadata, -1));
      }
    };

    if (usepTcategorize) {
      switch (resoFuncMatch) {
        case 6:
          loadCCDBforIP(pathsCCDBforIPIncparameter, resoFuncIncCCDB, "pathsCCDBforIPIncparameter");
          break;

        case 7:
          loadCCDBforIP(pathsCCDBforIPBeautyparameter, resoFuncBeautyCCDB, "pathsCCDBforIPBeautyparameter");
          loadCCDBforIP(pathsCCDBforIPCharmparameter, resoFuncCharmCCDB, "pathsCCDBforIPCharmparameter");
          loadCCDBforIP(pathsCCDBforIPLfparameter, resoFuncLfCCDB, "pathsCCDBforIPLfparameter");
          break;

        case 8:
          loadCCDBforIP(pathsCCDBforIPDataparameter, resoFuncDataCCDB, "pathsCCDBforIPDataparameter");
          break;

        default:
          LOG(info) << "resoFuncMatching is neither 6 nor 7, although usepTcategorize is set to true. Resetting resoFuncMatching to 0.";
          resoFuncMatch = 0;
          break;
      }
    }

    maxOrder = numCount + 1; // 0: untagged, >1 : N ordering
    const int mIPmethodNumOfParameters = 9;

    // Set up the resolution function
    switch (resoFuncMatch) {
      case 0:
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsCharmJetMC = (std::vector<float>)paramsResoFuncCharmJetMC;
        vecParamsBeautyJetMC = (std::vector<float>)paramsResoFuncBeautyJetMC;
        vecParamsLfJetMC = (std::vector<float>)paramsResoFuncLfJetMC;
        LOG(info) << "defined parameters of resolution function: custom";
        break;
      case 1:
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = (std::vector<float>)paramsResoFuncIncJetMC;
        useResoFuncFromIncJet = true;
        LOG(info) << "defined parameters of resolution function: custom & use inclusive distribution";
        break;
      case 2: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsCharmJetMC = {282119.753, -0.065, 0.893, 11.608, 0.945, 8.029, 0.131, 6.244, 0.027};
        vecParamsBeautyJetMC = {74901.583, -0.082, 0.874, 10.332, 0.941, 7.352, 0.097, 6.220, 0.022};
        vecParamsLfJetMC = {1539435.343, -0.061, 0.896, 13.272, 1.034, 5.884, 0.004, 7.843, 0.090};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, MB, LHC23d1k";
        break;
      case 3: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = {1908803.027, -0.059, 0.895, 13.467, 1.005, 8.867, 0.098, 6.929, 0.011};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, MB, LHC23d1k & use inclusive distribution";
        useResoFuncFromIncJet = true;
        break;
      case 4: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsCharmJetMC = {743719.121, -0.960, -0.240, 13.765, 1.314, 10.761, 0.293, 8.538, 0.052};
        vecParamsBeautyJetMC = {88888.418, 0.256, 1.003, 10.185, 0.740, 8.216, 0.147, 7.228, 0.040};
        vecParamsLfJetMC = {414860.372, -1.000, 0.285, 14.561, 1.464, 11.693, 0.339, 9.183, 0.052};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, JJ, weighted, LHC24g4";
        break;
      case 5: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        vecParamsIncJetMC = {2211391.862, 0.360, 1.028, 13.019, 0.650, 11.151, 0.215, 9.462, 0.044};
        LOG(info) << "defined parameters of resolution function: PYTHIA8, JJ, weighted, LHC24g4 & use inclusive distribution";
        useResoFuncFromIncJet = true;
        break;
      case 6: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        for (size_t j = 0; j < resoFuncIncCCDB.size(); j++) {
          std::vector<float> params;
          if (resoFuncIncCCDB[j]) {
            for (int i = 0; i < mIPmethodNumOfParameters; i++) {
              params.emplace_back(resoFuncIncCCDB[j]->GetParameter(i));
            }
          }
          vecParamsIncJetMcCCDB.emplace_back(params);
        }
        LOG(info) << "defined parameters of resolution function from CCDB";
        useResoFuncFromIncJet = true;
        break;
      case 7: // TODO
        vecParamsData = (std::vector<float>)paramsResoFuncData;
        for (size_t j = 0; j < resoFuncBeautyCCDB.size(); j++) {
          std::vector<float> params;
          if (resoFuncBeautyCCDB[j]) {
            for (int i = 0; i < mIPmethodNumOfParameters; i++) {
              params.emplace_back(resoFuncBeautyCCDB[j]->GetParameter(i));
            }
          }
          vecParamsBeautyJetMcCCDB.emplace_back(params);
        }
        for (size_t j = 0; j < resoFuncCharmCCDB.size(); j++) {
          std::vector<float> params;
          if (resoFuncCharmCCDB[j]) {
            for (int i = 0; i < mIPmethodNumOfParameters; i++) {
              params.emplace_back(resoFuncCharmCCDB[j]->GetParameter(i));
            }
          }
          vecParamsCharmJetMcCCDB.emplace_back(params);
        }
        for (size_t j = 0; j < resoFuncLfCCDB.size(); j++) {
          std::vector<float> params;
          if (resoFuncLfCCDB[j]) {
            for (int i = 0; i < mIPmethodNumOfParameters; i++) {
              params.emplace_back(resoFuncLfCCDB[j]->GetParameter(i));
            }
          }
          vecParamsLfJetMcCCDB.emplace_back(params);
        }
        LOG(info) << "defined parameters of resolution function from CCDB for each flavour";
        break;
      case 8:
        for (size_t j = 0; j < resoFuncDataCCDB.size(); j++) {
          std::vector<float> params;
          if (resoFuncDataCCDB[j]) {
            for (int i = 0; i < mIPmethodNumOfParameters; i++) {
              params.emplace_back(resoFuncDataCCDB[j]->GetParameter(i));
            }
          }
          vecParamsDataJetCCDB.emplace_back(params);
        }
        LOG(info) << "defined parameters of resolution function from CCDB for data";
        break;
      default:
        LOG(fatal) << "undefined parameters of resolution function. Fix it!";
        break;
    }

    fSignImpXYSigData = jettaggingutilities::setResolutionFunction(vecParamsData);
    fSignImpXYSigIncJetMC = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC);
    fSignImpXYSigCharmJetMC = jettaggingutilities::setResolutionFunction(vecParamsCharmJetMC);
    fSignImpXYSigBeautyJetMC = jettaggingutilities::setResolutionFunction(vecParamsBeautyJetMC);
    fSignImpXYSigLfJetMC = jettaggingutilities::setResolutionFunction(vecParamsLfJetMC);

    for (const auto& params : vecParamsDataJetCCDB) {
      vecfSignImpXYSigDataJetCCDB.emplace_back(jettaggingutilities::setResolutionFunction(params));
    }
    for (const auto& params : vecParamsIncJetMcCCDB) {
      vecfSignImpXYSigIncJetMcCCDB.emplace_back(jettaggingutilities::setResolutionFunction(params));
    }
    for (const auto& params : vecParamsBeautyJetMcCCDB) {
      vecfSignImpXYSigBeautyJetMcCCDB.emplace_back(jettaggingutilities::setResolutionFunction(params));
    }
    for (const auto& params : vecParamsCharmJetMcCCDB) {
      vecfSignImpXYSigCharmJetMcCCDB.emplace_back(jettaggingutilities::setResolutionFunction(params));
    }
    for (const auto& params : vecParamsLfJetMcCCDB) {
      vecfSignImpXYSigLfJetMcCCDB.emplace_back(jettaggingutilities::setResolutionFunction(params));
    }

    // Use QA for effectivness of track probability
    if (trackProbQA) {
      AxisSpec trackProbabilityAxis = {binTrackProbability, "Track proability"};
      AxisSpec jetFlavourAxis = {binJetFlavour, "Jet flavour"};
      if (doprocessFillTables) {
        if (isMCD) {
          registry.add("h2_pos_track_probability_flavour", "positive track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
          registry.add("h2_neg_track_probability_flavour", "negative track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
        } else {
          registry.add("h_pos_track_probability", "positive track probability", {HistType::kTH1F, {{trackProbabilityAxis}}});
          registry.add("h_neg_track_probability", "negative track probability", {HistType::kTH1F, {{trackProbabilityAxis}}});
        }
      }
    }

    if (doprocessAlgorithmML || doprocessAlgorithmGNN) {
      bMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        bMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        bMlResponse.setModelPathsLocal(onnxFileNames);
      }
      bMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      bMlResponse.init();
    }

    if (doprocessAlgorithmGNN) {
      tensorAlloc = o2::analysis::GNNBjetAllocator(nJetFeat.value, nTrkFeat.value, nClassesMl.value, nTrkOrigin.value, transformFeatureJetMean.value, transformFeatureJetStdev.value, transformFeatureTrkMean.value, transformFeatureTrkStdev.value, nJetConst);

      registry.add("h2_count_db", "#it{D}_{b} underflow/overflow;Jet flavour;#it{D}_{b} range", {HistType::kTH2F, {{4, 0., 4.}, {3, 0., 3.}}});
      auto h2CountDb = registry.get<TH2>(HIST("h2_count_db"));
      h2CountDb->GetXaxis()->SetBinLabel(1, "b-jet");
      h2CountDb->GetXaxis()->SetBinLabel(2, "c-jet");
      h2CountDb->GetXaxis()->SetBinLabel(3, "lf-jet");
      h2CountDb->GetXaxis()->SetBinLabel(4, "incl jet");
      h2CountDb->GetYaxis()->SetBinLabel(1, "underflow");
      h2CountDb->GetYaxis()->SetBinLabel(2, "in range");
      h2CountDb->GetYaxis()->SetBinLabel(3, "overflow");

      registry.add("h_db_b", "#it{D}_{b} b-jet;#it{D}_{b}", {HistType::kTH1F, {{50, -10., 35.}}});
      registry.add("h_db_c", "#it{D}_{b} c-jet;#it{D}_{b}", {HistType::kTH1F, {{50, -10., 35.}}});
      registry.add("h_db_lf", "#it{D}_{b} lf-jet;#it{D}_{b}", {HistType::kTH1F, {{50, -10., 35.}}});
      registry.add("h2_pt_db", "#it{p}_{T} vs. #it{D}_{b};#it{p}_{T}^{ch jet} (GeV/#it{c}^{2});#it{D}_{b}", {HistType::kTH2F, {{100, 0., 200.}, {50, -10., 35.}}});
    }
  }

  template <typename AnyJets, typename AnyTracks, typename SecondaryVertices>
  void analyzeJetAlgorithmML(AnyJets const& alljets, AnyTracks const& allTracks, SecondaryVertices const& allSVs)
  {
    for (const auto& analysisJet : alljets) {

      std::vector<jettaggingutilities::BJetTrackParams> tracksParams;
      std::vector<jettaggingutilities::BJetSVParams> svsParams;

      jettaggingutilities::analyzeJetSVInfo4ML(analysisJet, allTracks, allSVs, svsParams, svPtMin, svReductionFactor);
      jettaggingutilities::analyzeJetTrackInfo4ML(analysisJet, allTracks, allSVs, tracksParams, trackPtMin, trackDcaXYMax, trackDcaZMax);

      int nSVs = analysisJet.template secondaryVertices_as<SecondaryVertices>().size();

      jettaggingutilities::BJetParams jetparam = {analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), static_cast<int>(tracksParams.size()), static_cast<int>(nSVs), analysisJet.mass()};
      tracksParams.resize(nJetConst); // resize to the number of inputs of the ML
      svsParams.resize(nJetConst);    // resize to the number of inputs of the ML

      std::vector<float> output;

      if (bMlResponse.getInputShape().size() > 1) {
        auto inputML = bMlResponse.getInputFeatures2D(jetparam, tracksParams, svsParams);
        bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);
      } else {
        auto inputML = bMlResponse.getInputFeatures1D(jetparam, tracksParams, svsParams);
        bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);
      }

      if (bMlResponse.getOutputNodes() > 1) {
        auto mDb = [](const std::vector<float>& scores, float fC) {
          return std::log(scores[2] / (fC * scores[1] + (1 - fC) * scores[0]));
        };

        scoreML[analysisJet.globalIndex()] = useDb ? mDb(output, fC) : output[2]; // 2 is the b-jet index
      } else {
        scoreML[analysisJet.globalIndex()] = output[0];
      }
    }
  }

  template <typename AnyJets, typename AnyTracks>
  void analyzeJetAlgorithmMLnoSV(AnyJets const& alljets, AnyTracks const& allTracks)
  {
    for (const auto& analysisJet : alljets) {

      std::vector<jettaggingutilities::BJetTrackParams> tracksParams;
      std::vector<jettaggingutilities::BJetSVParams> svsParams;

      jettaggingutilities::analyzeJetTrackInfo4MLnoSV(analysisJet, allTracks, tracksParams, trackPtMin, trackDcaXYMax, trackDcaZMax);

      jettaggingutilities::BJetParams jetparam = {analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), static_cast<int>(tracksParams.size()), 0, analysisJet.mass()};
      tracksParams.resize(nJetConst); // resize to the number of inputs of the ML

      std::vector<float> output;

      if (bMlResponse.getInputShape().size() > 1) {
        auto inputML = bMlResponse.getInputFeatures2D(jetparam, tracksParams, svsParams);
        bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);
      } else {
        auto inputML = bMlResponse.getInputFeatures1D(jetparam, tracksParams, svsParams);
        bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);
      }

      scoreML[analysisJet.globalIndex()] = output[0];
    }
  }

  template <typename AnyJets, typename AnyTracks, typename AnyOriginalTracks>
  void analyzeJetAlgorithmGNN(AnyJets const& jets, AnyTracks const& tracks, AnyOriginalTracks const& origTracks)
  {
    for (const auto& jet : jets) {
      std::vector<std::vector<float>> trkFeat;
      jettaggingutilities::analyzeJetTrackInfo4GNN(jet, tracks, origTracks, trkFeat, trackPtMin, nJetConst);

      std::vector<float> jetFeat{jet.pt(), jet.phi(), jet.eta(), jet.mass()};

      if (trkFeat.size() > 0) {
        std::vector<float> feat;
        std::vector<Ort::Value> gnnInput;
        tensorAlloc.getGNNInput(jetFeat, trkFeat, feat, gnnInput);

        auto modelOutput = bMlResponse.getModelOutput(gnnInput, 0);
        scoreML[jet.globalIndex()] = jettaggingutilities::getDb(modelOutput, fC);
      } else {
        scoreML[jet.globalIndex()] = -999.;
        LOGF(debug, "doprocessAlgorithmGNN, trkFeat.size() <= 0 (%d)", jet.globalIndex());
      }
    }
  }

  void processDummy(aod::JetCollisions const&)
  {
  }
  PROCESS_SWITCH(JetTaggerHFTask, processDummy, "Dummy process", true);

  void processSetup(JetTable const& jets)
  {
    decisionNonML.clear();
    scoreML.clear();
    decisionNonML.resize(jets.size());
    scoreML.resize(jets.size());
  }
  PROCESS_SWITCH(JetTaggerHFTask, processSetup, "Setup initialization and size of jets for filling table", false);

  void processIP(JetTable const& jets, JetTracksExt const& tracks)
  {
    for (const auto& jet : jets) {
      uint16_t bit = jettaggingutilities::setTaggingIPBit(jet, tracks, trackDcaXYMax, trackDcaZMax, tagPointForIP);
      decisionNonML[jet.globalIndex()] |= bit;
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processIP, "Fill tagging decision for jet with the IP algorithm", false);

  void processSV(soa::Join<JetTable, SVIndicesTable> const& jets, SVTable const& prongs)
  {
    for (const auto& jet : jets) {
      uint16_t bit = jettaggingutilities::setTaggingSVBit(jet, prongs, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, svDispersionMax, tagPointForSV);
      decisionNonML[jet.globalIndex()] |= bit;
    }
  }
  PROCESS_SWITCH(JetTaggerHFTask, processSV, "Fill tagging decision for jets with the SV", false);

  void processAlgorithmML(soa::Join<JetTable, SVIndicesTable> const& allJets, JetTracksExt const& allTracks, SVTable const& allSVs)
  {
    analyzeJetAlgorithmML(allJets, allTracks, allSVs);
  }
  PROCESS_SWITCH(JetTaggerHFTask, processAlgorithmML, "Fill ML evaluation score for charged jets", false);

  void processAlgorithmMLnoSV(JetTable const& allJets, JetTracksExt const& allTracks)
  {
    analyzeJetAlgorithmMLnoSV(allJets, allTracks);
  }
  PROCESS_SWITCH(JetTaggerHFTask, processAlgorithmMLnoSV, "Fill ML evaluation score for charged jets but without using SVs", false);

  void processAlgorithmGNN(JetTable const& jets, JetTracksExt const& jtracks, OriginalTracks const& origTracks)
  {
    analyzeJetAlgorithmGNN(jets, jtracks, origTracks);
  }
  PROCESS_SWITCH(JetTaggerHFTask, processAlgorithmGNN, "Fill GNN evaluation score (D_b) for charged jets", false);

  void processFillTables(std::conditional_t<isMCD, soa::Join<JetTable, aod::ChargedMCDetectorLevelJetFlavourDef>, JetTable>::iterator const& jet, JetTracksExt const& tracks)
  {
    fillTables<isMCD>(jet, tracks);
  }
  PROCESS_SWITCH(JetTaggerHFTask, processFillTables, "Fill Tables for tagging decision, jet probability, and ML score on charged jets", false);
};

using JetTaggerhfDataCharged = JetTaggerHFTask<false, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, aod::DataSecondaryVertex3ProngIndices, aod::DataSecondaryVertex3Prongs, aod::ChargedJetTags>;
using JetTaggerhfMCDCharged = JetTaggerHFTask<true, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::MCDSecondaryVertex3ProngIndices, aod::MCDSecondaryVertex3Prongs, aod::ChargedMCDetectorLevelJetTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{
    adaptAnalysisTask<JetTaggerhfDataCharged>(cfgc, SetDefaultProcesses{}, TaskName{"jet-taggerhf-data-charged"}), // o2-linter: disable=name/o2-task
    adaptAnalysisTask<JetTaggerhfMCDCharged>(cfgc, SetDefaultProcesses{}, TaskName{"jet-taggerhf-mcd-charged"})};  // o2-linter: disable=name/o2-task
}
