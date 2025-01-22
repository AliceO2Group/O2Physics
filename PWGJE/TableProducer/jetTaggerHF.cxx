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

#include <string>
#include <memory>
#include <vector>

#include <TF1.h>
#include <TH1.h>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/trackUtilities.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "Tools/ML/MlResponse.h"

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
  Configurable<std::vector<float>> paramsResoFuncData{"paramsResoFuncData", std::vector<float>{1306800, -0.1049, 0.861425, 13.7547, 0.977967, 8.96823, 0.151595, 6.94499, 0.0250301}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7))"};
  Configurable<std::vector<float>> paramsResoFuncIncJetMC{"paramsResoFuncIncJetMC", std::vector<float>{1908803.027, -0.059, 0.895, 13.467, 1.005, 8.867, 0.098, 6.929, 0.011}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncCharmJetMC{"paramsResoFuncCharmJetMC", std::vector<float>{282119.753, -0.065, 0.893, 11.608, 0.945, 8.029, 0.131, 6.244, 0.027}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncBeautyJetMC{"paramsResoFuncBeautyJetMC", std::vector<float>{74901.583, -0.082, 0.874, 10.332, 0.941, 7.352, 0.097, 6.220, 0.022}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
  Configurable<std::vector<float>> paramsResoFuncLfJetMC{"paramsResoFuncLfJetMC", std::vector<float>{1539435.343, -0.061, 0.896, 13.272, 1.034, 5.884, 0.004, 7.843, 0.090}, "parameters of gaus(0)+expo(3)+expo(5)+expo(7)))"};
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

  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"Users/h/hahassan"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ML_bjets/Models/LHC24g4_70_200/model.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  Configurable<std::string> IPparameterPathsCCDB{"IPparameterPathsCCDB", "Users/l/leehy/LHC24g4/", "Paths for fitting parameters of resolution functions for IP method on CCDB"};
  Configurable<std::vector<int64_t>> IPtimestampCCDB{"IPtimestampCCDB", std::vector<int64_t>{1737027389227, 1737027391774, 1737027393668, 1737027395548, 1737027397505, 1737027399396, 1737027401294}, "timestamp of the resolution function for IP method used to query in CCDB"};
  Configurable<bool> usepTcategorize{"usepTcategorize", false, "p_T categorize TF1 function with Inclusive jet"};

  // axis spec
  ConfigurableAxis binTrackProbability{"binTrackProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};

  o2::analysis::MlResponse<float> bMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  using JetTracksExt = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;

  bool useResoFuncFromIncJet = false;
  int maxOrder = -1;
  int resoFuncMatch = 0;

  std::unique_ptr<TF1> fSignImpXYSigData = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigCharmJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigBeautyJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigLfJetMC = nullptr;

  std::vector<float> vecParamsIncJetMC_CCDB_0;
  std::vector<float> vecParamsIncJetMC_CCDB_1;
  std::vector<float> vecParamsIncJetMC_CCDB_2;
  std::vector<float> vecParamsIncJetMC_CCDB_3;
  std::vector<float> vecParamsIncJetMC_CCDB_4;
  std::vector<float> vecParamsIncJetMC_CCDB_5;
  std::vector<float> vecParamsIncJetMC_CCDB_6;

  std::unique_ptr<TF1> fSignImpXYSigIncJetMC_CCDB_0 = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC_CCDB_1 = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC_CCDB_2 = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC_CCDB_3 = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC_CCDB_4 = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC_CCDB_5 = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC_CCDB_6 = nullptr;

  std::vector<std::unique_ptr<TF1>> fSignImpXYSigIncJetMC_CCDB_vec;

  std::vector<uint16_t> decisionNonML;
  std::vector<float> scoreML;

  template <typename T, typename U>
  float calculateJetProbability(int origin, T const& jet, U const& tracks, bool const& isMC = false)
  {
    float jetProb = -1.0;
    if (!isMC) {
      jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigData, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
    } else {
      if (useResoFuncFromIncJet) {
        if (usepTcategorize) {
          jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigIncJetMC_CCDB_vec, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
        } else {
          jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigIncJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
        }
      } else {
        if (origin == JetTaggingSpecies::charm) {
          jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigCharmJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
        } else if (origin == JetTaggingSpecies::beauty) {
          jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigBeautyJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
        } else {
          jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigLfJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
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
        jetProb = calculateJetProbability(origin, jet, tracks);
        if (trackProbQA) {
          evaluateTrackProbQA(origin, jet, tracks, isMC);
        }
      } else {
        jetProb = calculateJetProbability(0, jet, tracks);
        if (trackProbQA) {
          evaluateTrackProbQA(0, jet, tracks, isMC);
        }
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

    TF1* CCDB_ResoFunc_0 = nullptr;
    TF1* CCDB_ResoFunc_1 = nullptr;
    TF1* CCDB_ResoFunc_2 = nullptr;
    TF1* CCDB_ResoFunc_3 = nullptr;
    TF1* CCDB_ResoFunc_4 = nullptr;
    TF1* CCDB_ResoFunc_5 = nullptr;
    TF1* CCDB_ResoFunc_6 = nullptr;

    ccdbApi.init(ccdbUrl);
    if (usepTcategorize) {
      std::map<std::string, std::string> metadata; // dummy meta data (will be updated)
      // fill the timestamp directly of each TF1 according to p_T track range
      CCDB_ResoFunc_0 = ccdbApi.retrieveFromTFileAny<TF1>(IPparameterPathsCCDB, metadata, IPtimestampCCDB->at(0)); // 0 < p_T < 0.5
      CCDB_ResoFunc_1 = ccdbApi.retrieveFromTFileAny<TF1>(IPparameterPathsCCDB, metadata, IPtimestampCCDB->at(1)); // 0.5 < p_T < 1
      CCDB_ResoFunc_2 = ccdbApi.retrieveFromTFileAny<TF1>(IPparameterPathsCCDB, metadata, IPtimestampCCDB->at(2)); // 1 < p_T < 2
      CCDB_ResoFunc_3 = ccdbApi.retrieveFromTFileAny<TF1>(IPparameterPathsCCDB, metadata, IPtimestampCCDB->at(3)); // 2 < p_T < 4
      CCDB_ResoFunc_4 = ccdbApi.retrieveFromTFileAny<TF1>(IPparameterPathsCCDB, metadata, IPtimestampCCDB->at(4)); // 4 < p_T < 6
      CCDB_ResoFunc_5 = ccdbApi.retrieveFromTFileAny<TF1>(IPparameterPathsCCDB, metadata, IPtimestampCCDB->at(5)); // 6 < p_T < 9
      CCDB_ResoFunc_6 = ccdbApi.retrieveFromTFileAny<TF1>(IPparameterPathsCCDB, metadata, IPtimestampCCDB->at(6)); // 9 < p_T
    }

    maxOrder = numCount + 1; // 0: untagged, >1 : N ordering

    // Set up the resolution function
    resoFuncMatch = resoFuncMatching;
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
        for (int i = 0; i < 9; i++) {
          vecParamsIncJetMC_CCDB_0.emplace_back(CCDB_ResoFunc_0->GetParameter(i));
          vecParamsIncJetMC_CCDB_1.emplace_back(CCDB_ResoFunc_1->GetParameter(i));
          vecParamsIncJetMC_CCDB_2.emplace_back(CCDB_ResoFunc_2->GetParameter(i));
          vecParamsIncJetMC_CCDB_3.emplace_back(CCDB_ResoFunc_3->GetParameter(i));
          vecParamsIncJetMC_CCDB_4.emplace_back(CCDB_ResoFunc_4->GetParameter(i));
          vecParamsIncJetMC_CCDB_5.emplace_back(CCDB_ResoFunc_5->GetParameter(i));
          vecParamsIncJetMC_CCDB_6.emplace_back(CCDB_ResoFunc_6->GetParameter(i));
        }
        LOG(info) << "defined parameters of resolution function from CCDB";
        useResoFuncFromIncJet = true;
      default:
        LOG(fatal) << "undefined parameters of resolution function. Fix it!";
        break;
    }

    fSignImpXYSigData = jettaggingutilities::setResolutionFunction(vecParamsData);
    fSignImpXYSigIncJetMC = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC);
    fSignImpXYSigCharmJetMC = jettaggingutilities::setResolutionFunction(vecParamsCharmJetMC);
    fSignImpXYSigBeautyJetMC = jettaggingutilities::setResolutionFunction(vecParamsBeautyJetMC);
    fSignImpXYSigLfJetMC = jettaggingutilities::setResolutionFunction(vecParamsLfJetMC);

    fSignImpXYSigIncJetMC_CCDB_0 = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC_CCDB_0);
    fSignImpXYSigIncJetMC_CCDB_1 = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC_CCDB_1);
    fSignImpXYSigIncJetMC_CCDB_2 = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC_CCDB_2);
    fSignImpXYSigIncJetMC_CCDB_3 = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC_CCDB_3);
    fSignImpXYSigIncJetMC_CCDB_4 = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC_CCDB_4);
    fSignImpXYSigIncJetMC_CCDB_5 = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC_CCDB_5);
    fSignImpXYSigIncJetMC_CCDB_6 = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC_CCDB_6);

    fSignImpXYSigIncJetMC_CCDB_vec.emplace_back(std::move(fSignImpXYSigIncJetMC_CCDB_0));
    fSignImpXYSigIncJetMC_CCDB_vec.emplace_back(std::move(fSignImpXYSigIncJetMC_CCDB_1));
    fSignImpXYSigIncJetMC_CCDB_vec.emplace_back(std::move(fSignImpXYSigIncJetMC_CCDB_2));
    fSignImpXYSigIncJetMC_CCDB_vec.emplace_back(std::move(fSignImpXYSigIncJetMC_CCDB_3));
    fSignImpXYSigIncJetMC_CCDB_vec.emplace_back(std::move(fSignImpXYSigIncJetMC_CCDB_4));
    fSignImpXYSigIncJetMC_CCDB_vec.emplace_back(std::move(fSignImpXYSigIncJetMC_CCDB_5));
    fSignImpXYSigIncJetMC_CCDB_vec.emplace_back(std::move(fSignImpXYSigIncJetMC_CCDB_6));

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

    if (doprocessAlgorithmML) {
      bMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        bMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        bMlResponse.setModelPathsLocal(onnxFileNames);
      }
      // bMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      bMlResponse.init();
    }
  }

  template <typename AnyJets, typename AnyTracks, typename SecondaryVertices>
  void analyzeJetAlgorithmML(AnyJets const& alljets, AnyTracks const& allTracks, SecondaryVertices const& allSVs)
  {
    for (const auto& analysisJet : alljets) {

      std::vector<jettaggingutilities::BJetTrackParams> tracksParams;
      std::vector<jettaggingutilities::BJetSVParams> svsParams;

      analyzeJetSVInfo4ML(analysisJet, allTracks, allSVs, svsParams, svPtMin, svReductionFactor);
      analyzeJetTrackInfo4ML(analysisJet, allTracks, allSVs, tracksParams, trackPtMin);

      int nSVs = analysisJet.template secondaryVertices_as<SecondaryVertices>().size();

      jettaggingutilities::BJetParams jetparam = {analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), static_cast<int>(tracksParams.size()), static_cast<int>(nSVs), analysisJet.mass()};
      tracksParams.resize(nJetConst); // resize to the number of inputs of the ML
      svsParams.resize(nJetConst);    // resize to the number of inputs of the ML

      auto inputML = getInputsForML(jetparam, tracksParams, svsParams, nJetConst);

      std::vector<float> output;
      // bool isSelectedMl = bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);
      bMlResponse.isSelectedMl(inputML, analysisJet.pt(), output);

      scoreML[analysisJet.globalIndex()] = output[0];
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
