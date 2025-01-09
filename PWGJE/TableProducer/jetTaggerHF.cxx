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

template <typename JetTable, typename SVIndicesTable, typename SVTable, typename JetTaggingTable>
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

  // axis spec
  ConfigurableAxis binTrackProbability{"binTrackProbability", {100, 0.f, 1.f}, ""};
  ConfigurableAxis binJetFlavour{"binJetFlavour", {6, -0.5, 5.5}, ""};

  o2::analysis::MlResponse<float> bMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  using JetTracksExt = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;

  bool useResoFuncFromIncJet = false;
  int maxOrder = -1;
  int resoFuncMatch = 0;
  std::vector<float> jetProb;
  std::unique_ptr<TF1> fSignImpXYSigData = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigIncJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigCharmJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigBeautyJetMC = nullptr;
  std::unique_ptr<TF1> fSignImpXYSigLfJetMC = nullptr;

  std::vector<uint16_t> decisionNonML;
  std::vector<float> scoreML;

  template <typename T, typename U>
  float calculateJetProbability(int origin, T const& jet, U const& tracks, bool const& isMC = true)
  {
    float jetProb = -1.0;
    if (!isMC) {
      jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigData, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
    } else {
      if (useResoFuncFromIncJet) {
        jetProb = jettaggingutilities::getJetProbability(fSignImpXYSigIncJetMC, jet, tracks, trackDcaXYMax, trackDcaZMax, minSignImpXYSig);
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

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    std::vector<float> vecParamsData;
    std::vector<float> vecParamsIncJetMC;
    std::vector<float> vecParamsCharmJetMC;
    std::vector<float> vecParamsBeautyJetMC;
    std::vector<float> vecParamsLfJetMC;

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
      default:
        LOG(fatal) << "undefined parameters of resolution function. Fix it!";
        break;
    }

    fSignImpXYSigData = jettaggingutilities::setResolutionFunction(vecParamsData);
    fSignImpXYSigIncJetMC = jettaggingutilities::setResolutionFunction(vecParamsIncJetMC);
    fSignImpXYSigCharmJetMC = jettaggingutilities::setResolutionFunction(vecParamsCharmJetMC);
    fSignImpXYSigBeautyJetMC = jettaggingutilities::setResolutionFunction(vecParamsBeautyJetMC);
    fSignImpXYSigLfJetMC = jettaggingutilities::setResolutionFunction(vecParamsLfJetMC);

    // Use QA for effectivness of track probability
    if (trackProbQA) {
      AxisSpec trackProbabilityAxis = {binTrackProbability, "Track proability"};
      AxisSpec jetFlavourAxis = {binJetFlavour, "Jet flavour"};
      if (doprocessFillTablesData) {
        registry.add("h_pos_track_probability", "positive track probability", {HistType::kTH1F, {{trackProbabilityAxis}}});
        registry.add("h_neg_track_probability", "negative track probability", {HistType::kTH1F, {{trackProbabilityAxis}}});
      }
      if (doprocessFillTablesMCD) {
        registry.add("h2_pos_track_probability_flavour", "positive track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
        registry.add("h2_neg_track_probability_flavour", "negative track probability", {HistType::kTH2F, {{trackProbabilityAxis}, {jetFlavourAxis}}});
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

  void processFillTablesData(JetTable::iterator const& jet, JetTracksExt const& tracks)
  {
    float jetProb = -1.0;
    if (useJetProb) {
      jetProb = calculateJetProbability(0, jet, tracks, false);
      if (trackProbQA) {
        evaluateTrackProbQA(0, jet, tracks, false);
      }
    }
    taggingTable(decisionNonML[jet.globalIndex()], jetProb, scoreML[jet.globalIndex()]);
  }
  PROCESS_SWITCH(JetTaggerHFTask, processFillTablesData, "Fill Tables for tagging decision, jet probability, and ML score on charged jets in data", false);

  void processFillTablesMCD(soa::Join<JetTable, aod::ChargedMCDetectorLevelJetFlavourDef>::iterator const& jet, JetTracksExt const& tracks)
  {
    float jetProb = -1.0;
    int origin = jet.origin();
    if (useJetProb) {
      jetProb = calculateJetProbability(origin, jet, tracks);
      if (trackProbQA) {
        evaluateTrackProbQA(origin, jet, tracks);
      }
    }
    taggingTable(decisionNonML[jet.globalIndex()], jetProb, scoreML[jet.globalIndex()]);
  }
  PROCESS_SWITCH(JetTaggerHFTask, processFillTablesMCD, "Fill Tables for tagging decision, jet probability, and ML score on charged jets in MC", false);
};

using JetTaggerhfDataCharged = JetTaggerHFTask<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>, aod::DataSecondaryVertex3ProngIndices, aod::DataSecondaryVertex3Prongs, aod::ChargedJetTags>;
using JetTaggerhfMCDCharged = JetTaggerHFTask<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::MCDSecondaryVertex3ProngIndices, aod::MCDSecondaryVertex3Prongs, aod::ChargedMCDetectorLevelJetTags>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{
    adaptAnalysisTask<JetTaggerhfDataCharged>(cfgc, TaskName{"jet-taggerhf-data-charged"}),
    adaptAnalysisTask<JetTaggerhfMCDCharged>(cfgc, TaskName{"jet-taggerhf-mcd-charged"})};
}
