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

/// \file candidateSelectorLcToK0sP.cxx
/// \brief Lc --> K0s+p selection task.
/// \note based on candidateSelectorD0.cxx
///
/// \author Chiara Zampolli <Chiara.Zampolli@cern.ch>, CERN
///         Daniel Samitz, <daniel.samitz@cern.ch>, Vienna

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsDebugLcToK0sP.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc;
using namespace o2::analysis::hf_cuts_lc_to_k0s_p;

using MyBigTracksBayes = soa::Join<aod::BigTracksPID, aod::pidBayesPr, aod::pidBayesEl, aod::pidBayesMu, aod::pidBayesKa, aod::pidBayesPi>;

struct HfCandidateSelectorLcToK0sP {
  Produces<aod::HfSelLcToK0sP> hfSelLcToK0sPCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  // PID
  Configurable<double> pPidThreshold{"pPidThreshold", 1.0, "Threshold to switch between low and high p TrackSelectors"};
  // TPC
  Configurable<double> nSigmaTpcMaxLowP{"nSigmaTpcMaxLowP", 2.0, "Max nSigma in TPC for bachelor at low p"};
  Configurable<double> nSigmaTpcCombinedMaxLowP{"nSigmaTpcCombinedMaxLowP", 2.0, "Max nSigma in TPC combined with TOF for bachelor at low p"};
  Configurable<double> nSigmaTpcMaxHighP{"nSigmaTpcMaxHighP", 0., "Max nSigma in TPC for bachelor at high p"};
  Configurable<double> nSigmaTpcCombinedMaxHighP{"nSigmaTpcCombinedMaxHighP", 9999., "Max nSigma in TPC combined with TOF for bachelor at high p"};
  // TOF
  Configurable<double> nSigmaTofMaxLowP{"nSigmaTofMaxLowP", 0., "Max nSigma in TOF for bachelor at low p"};
  Configurable<double> nSigmaTofCombinedMaxLowP{"nSigmaTofCombinedMaxLowP", 9999., "Max nSigma in TOF combined with TPC for bachelor at low p"};
  Configurable<double> nSigmaTofMaxHighP{"nSigmaTofMaxHighP", 3.0, "Max nSigma in TOF for bachelor at high p"};
  Configurable<double> nSigmaTofCombinedMaxHighP{"nSigmaTofCombinedMaxHighP", 3.0, "Max nSigma in TOF combined with TPC for bachelor at high p"};
  // Bayesian
  Configurable<double> probBayesMinLowP{"probBayesMinLowP", 0.8, "min. Bayes probability for bachelor at low p [%]"};
  Configurable<double> probBayesMinHighP{"probBayesMinHighP", 0.8, "min. Bayes probability for bachelor at high p [%]"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_k0s_p::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Lc candidate selection per pT bin"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<std::string>> modelPathsMl{"modelPathsMl", std::vector<std::string>{hf_cuts_ml::modelPaths}, "Paths of the ML models, one for each pT bin"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  Configurable<std::vector<std::string>> inputFeaturesML{"inputFeaturesML", {""}, "List of input features for the ML model"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> modelPathsCCDB{"modelPathsCCDB", "EventFiltering/PWGHF/BDTLcToK0sP", "Path on CCDB"};
  Configurable<std::vector<std::string>> onnxFilesCCDB{"onnxFilesCCDB", std::vector<std::string>{"ModelHandler_onnx_LcToK0sP.onnx"}, "ONNX file names on CCDB, for each pT bin"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponse<float> hfMlResponse;
  std::vector<bool> selectedInputFeatures;

  o2::ccdb::CcdbApi ccdbApi;

  std::vector<std::shared_ptr<TH1>> hModelScore;
  std::vector<std::shared_ptr<TH2>> hModelScoreVsPtCand;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (!doprocessWithStandardPID && !doprocessWithBayesPID) {
      LOGF(fatal, "Neither processWithStandardPID nor processWithBayesPID enabled. Please choose one.");
    }
    if (doprocessWithStandardPID && doprocessWithBayesPID) {
      LOGF(fatal, "Cannot enable processWithStandardPID and processWithBayesPID at the same time. Please choose one.");
    }

    if (applyMl) {
      hfMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl, modelPathsMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFilesCCDB, ccdbApi, modelPathsCCDB.value, timestampCCDB);
      }
      hfMlResponse.init();

      // load histograms for ML score
      AxisSpec axisScore = {100, 0, 1, "score"};
      AxisSpec axisBinsPt = {binsPtMl, "#it{p}_{T} (GeV/#it{c})"};
      for (int classMl = 0; classMl < nClassesMl; classMl++) {
        hModelScore.push_back(registry.add<TH1>(Form("hMlScoreClass%d", classMl), "Model core distribution for Lc;Model score;counts", HistType::kTH1F, {axisScore}));
        hModelScoreVsPtCand.push_back(registry.add<TH2>(Form("hMlScoreClass%dVsPtCand", classMl), "Model score distribution for Lc;Model score;counts", HistType::kTH2F, {axisScore, axisBinsPt}));
      }

      // set the selected input features
      initInputFeatures();
    }
  }

  // Check the list of selected input features for ML model
  void initInputFeatures()
  {
    std::string inputFeaturesAll[] = {
      "numContrib",
      "posX",
      "posY",
      "posZ",
      "xSecondaryVertex",
      "ySecondaryVertex",
      "zSecondaryVertex",
      "errorDecayLength",
      "errorDecayLengthXY",
      "chi2PCA",
      "rSecondaryVertex",
      "decayLength",
      "decayLengthXY",
      "decayLengthNormalised",
      "decayLengthXYNormalised",
      "impactParameterNormalised0",
      "ptProng0",
      "pProng0",
      "impactParameterNormalised1",
      "ptProng1",
      "pProng1",
      "pxProng0",
      "pyProng0",
      "pzProng0",
      "pxProng1",
      "pyProng1",
      "pzProng1",
      "impactParameter0",
      "impactParameter1",
      "errorImpactParameter0",
      "errorImpactParameter1",
      "v0X",
      "v0Y",
      "v0Z",
      "v0Radius",
      "v0CosPA",
      "v0MLambda",
      "v0MAntiLambda",
      "v0MK0Short",
      "v0MGamma",
      "v0CtK0Short",
      "v0CtLambda",
      "dcaV0Daughters",
      "pxPos",
      "pyPos",
      "pzPos",
      "ptV0Pos",
      "dcaPosToPV",
      "pxNeg",
      "pyNeg",
      "pzNeg",
      "ptV0Neg",
      "dcaNegToPV",
      "nSigmaTPCPr0",
      "nSigmaTOFPr0",
      "m",
      "pt",
      "p",
      "cpa",
      "cpaXY",
      "ct",
      "eta",
      "phi",
      "y",
      "e"};

    // check for each possible input feature if it is included in the list of selected input features or not
    for (const auto& inputFeature : inputFeaturesAll) {
      if (std::find(std::begin(inputFeaturesML.value), std::end(inputFeaturesML.value), inputFeature) != std::end(inputFeaturesML.value)) {
        selectedInputFeatures.push_back(true);
        LOG(info) << "Included \'" << inputFeature << "\' in list of ML input features.";
      } else {
        selectedInputFeatures.push_back(false);
      }
    }

    // check if all given input features are recongnized
    for (const auto& inputFeature : inputFeaturesML.value) {
      if (std::find(std::begin(inputFeaturesAll), std::end(inputFeaturesAll), inputFeature) == std::end(inputFeaturesAll)) {
        LOG(fatal) << "Can not find \'" << inputFeature << "\' in list of possible ML input features.";
      }
    }
  }

  // fill only the selcted input features into the the vector of ML input features
  template <typename T, typename U>
  std::vector<float> setInputFeatures(const T& candidate, const U& bach)
  {
    std::vector<float> inputFeatures;

    if (selectedInputFeatures[0]) {
      inputFeatures.push_back(bach.collision().numContrib());
    }
    if (selectedInputFeatures[1]) {
      inputFeatures.push_back(candidate.posX());
    }
    if (selectedInputFeatures[2]) {
      inputFeatures.push_back(candidate.posY());
    }
    if (selectedInputFeatures[3]) {
      inputFeatures.push_back(candidate.posZ());
    }
    if (selectedInputFeatures[4]) {
      inputFeatures.push_back(candidate.xSecondaryVertex());
    }
    if (selectedInputFeatures[5]) {
      inputFeatures.push_back(candidate.ySecondaryVertex());
    }
    if (selectedInputFeatures[6]) {
      inputFeatures.push_back(candidate.zSecondaryVertex());
    }
    if (selectedInputFeatures[7]) {
      inputFeatures.push_back(candidate.errorDecayLength());
    }
    if (selectedInputFeatures[8]) {
      inputFeatures.push_back(candidate.errorDecayLengthXY());
    }
    if (selectedInputFeatures[9]) {
      inputFeatures.push_back(candidate.chi2PCA());
    }
    if (selectedInputFeatures[10]) {
      inputFeatures.push_back(candidate.rSecondaryVertex());
    }
    if (selectedInputFeatures[11]) {
      inputFeatures.push_back(candidate.decayLength());
    }
    if (selectedInputFeatures[12]) {
      inputFeatures.push_back(candidate.decayLengthXY());
    }
    if (selectedInputFeatures[13]) {
      inputFeatures.push_back(candidate.decayLengthNormalised());
    }
    if (selectedInputFeatures[14]) {
      inputFeatures.push_back(candidate.decayLengthXYNormalised());
    }
    if (selectedInputFeatures[15]) {
      inputFeatures.push_back(candidate.impactParameterNormalised0());
    }
    if (selectedInputFeatures[16]) {
      inputFeatures.push_back(candidate.ptProng0());
    }
    if (selectedInputFeatures[17]) {
      inputFeatures.push_back(RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()));
    }
    if (selectedInputFeatures[18]) {
      inputFeatures.push_back(candidate.impactParameterNormalised1());
    }
    if (selectedInputFeatures[19]) {
      inputFeatures.push_back(candidate.ptProng1());
    }
    if (selectedInputFeatures[20]) {
      inputFeatures.push_back(RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()));
    }
    if (selectedInputFeatures[21]) {
      inputFeatures.push_back(candidate.pxProng0());
    }
    if (selectedInputFeatures[22]) {
      inputFeatures.push_back(candidate.pyProng0());
    }
    if (selectedInputFeatures[23]) {
      inputFeatures.push_back(candidate.pzProng0());
    }
    if (selectedInputFeatures[24]) {
      inputFeatures.push_back(candidate.pxProng1());
    }
    if (selectedInputFeatures[25]) {
      inputFeatures.push_back(candidate.pyProng1());
    }
    if (selectedInputFeatures[26]) {
      inputFeatures.push_back(candidate.pzProng1());
    }
    if (selectedInputFeatures[27]) {
      inputFeatures.push_back(candidate.impactParameter0());
    }
    if (selectedInputFeatures[28]) {
      inputFeatures.push_back(candidate.impactParameter1());
    }
    if (selectedInputFeatures[29]) {
      inputFeatures.push_back(candidate.errorImpactParameter0());
    }
    if (selectedInputFeatures[30]) {
      inputFeatures.push_back(candidate.errorImpactParameter1());
    }
    if (selectedInputFeatures[31]) {
      inputFeatures.push_back(candidate.v0x());
    }
    if (selectedInputFeatures[32]) {
      inputFeatures.push_back(candidate.v0y());
    }
    if (selectedInputFeatures[33]) {
      inputFeatures.push_back(candidate.v0z());
    }
    if (selectedInputFeatures[34]) {
      inputFeatures.push_back(candidate.v0radius());
    }
    if (selectedInputFeatures[35]) {
      inputFeatures.push_back(candidate.v0cosPA());
    }
    if (selectedInputFeatures[36]) {
      inputFeatures.push_back(candidate.mLambda());
    }
    if (selectedInputFeatures[37]) {
      inputFeatures.push_back(candidate.mAntiLambda());
    }
    if (selectedInputFeatures[38]) {
      inputFeatures.push_back(candidate.mK0Short());
    }
    if (selectedInputFeatures[39]) {
      inputFeatures.push_back(candidate.mGamma());
    }
    if (selectedInputFeatures[40]) {
      inputFeatures.push_back(o2::aod::hf_cand_casc::ctV0K0s(candidate));
    }
    if (selectedInputFeatures[41]) {
      inputFeatures.push_back(o2::aod::hf_cand_casc::ctV0Lambda(candidate));
    }
    if (selectedInputFeatures[42]) {
      inputFeatures.push_back(candidate.dcaV0daughters());
    }
    if (selectedInputFeatures[43]) {
      inputFeatures.push_back(candidate.pxpos());
    }
    if (selectedInputFeatures[44]) {
      inputFeatures.push_back(candidate.pypos());
    }
    if (selectedInputFeatures[45]) {
      inputFeatures.push_back(candidate.pzpos());
    }
    if (selectedInputFeatures[46]) {
      inputFeatures.push_back(candidate.ptV0Pos());
    }
    if (selectedInputFeatures[47]) {
      inputFeatures.push_back(candidate.dcapostopv());
    }
    if (selectedInputFeatures[48]) {
      inputFeatures.push_back(candidate.pxneg());
    }
    if (selectedInputFeatures[49]) {
      inputFeatures.push_back(candidate.pyneg());
    }
    if (selectedInputFeatures[50]) {
      inputFeatures.push_back(candidate.pzneg());
    }
    if (selectedInputFeatures[51]) {
      inputFeatures.push_back(candidate.ptV0Neg());
    }
    if (selectedInputFeatures[52]) {
      inputFeatures.push_back(candidate.dcanegtopv());
    }
    if (selectedInputFeatures[53]) {
      inputFeatures.push_back(bach.tpcNSigmaPr());
    }
    if (selectedInputFeatures[54]) {
      inputFeatures.push_back(bach.tofNSigmaPr());
    }
    if (selectedInputFeatures[55]) {
      inputFeatures.push_back(o2::aod::hf_cand_casc::invMassLcToK0sP(candidate));
    }
    if (selectedInputFeatures[56]) {
      inputFeatures.push_back(candidate.pt());
    }
    if (selectedInputFeatures[57]) {
      inputFeatures.push_back(candidate.p());
    }
    if (selectedInputFeatures[58]) {
      inputFeatures.push_back(candidate.cpa());
    }
    if (selectedInputFeatures[59]) {
      inputFeatures.push_back(candidate.cpaXY());
    }
    if (selectedInputFeatures[60]) {
      inputFeatures.push_back(o2::aod::hf_cand_3prong::ctLc(candidate));
    }
    if (selectedInputFeatures[61]) {
      inputFeatures.push_back(candidate.eta());
    }
    if (selectedInputFeatures[62]) {
      inputFeatures.push_back(candidate.phi());
    }
    if (selectedInputFeatures[63]) {
      inputFeatures.push_back(o2::aod::hf_cand_3prong::yLc(candidate));
    }
    if (selectedInputFeatures[64]) {
      inputFeatures.push_back(o2::aod::hf_cand_3prong::eLc(candidate));
    }

    return inputFeatures;
  }

  /// Conjugate independent toplogical cuts
  /// \param hfCandCascade is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& hfCandCascade)
  {
    auto candPt = hfCandCascade.pt();
    int ptBin = findBin(binsPt, candPt);
    if (ptBin == -1) {
      return false;
    }

    if (candPt < ptCandMin || candPt >= ptCandMax) {
      return false; // check that the candidate pT is within the analysis range
    }

    if (std::abs(hfCandCascade.mK0Short() - RecoDecay::getMassPDG(kK0Short)) > cuts->get(ptBin, "mK0s")) {
      return false; // mass of the K0s
    }

    if ((std::abs(hfCandCascade.mLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda")) || (std::abs(hfCandCascade.mAntiLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda"))) {
      return false; // mass of the Lambda
    }

    if (std::abs(hfCandCascade.mGamma() - RecoDecay::getMassPDG(kGamma)) < cuts->get(ptBin, "mGamma")) {
      return false; // mass of the Gamma
    }

    if (hfCandCascade.ptProng0() < cuts->get(ptBin, "ptBach")) {
      return false; // pt of the p
    }

    if (hfCandCascade.ptV0Pos() < cuts->get(ptBin, "ptV0Dau")) {
      return false; // pt of the pos K0 daughter
    }

    if (hfCandCascade.ptV0Neg() < cuts->get(ptBin, "ptV0Dau")) {
      return false; // pt of the neg K0 daughter
    }

    if (hfCandCascade.ptProng1() < cuts->get(ptBin, "ptV0")) {
      return false; // pt of the K0
    }

    if (std::abs(hfCandCascade.impactParameter0()) > cuts->get(ptBin, "d0Bach")) {
      return false; // d0 of the bachelor
    }

    /*
    if ((std::abs(hfCandCascade.dcapostopv()) > d0K0Cut[ptBin]) || (std::abs(hfCandCascade.dcanegtopv()) > d0K0Cut[ptBin])) {
      LOG(debug) << "v0 daugh cut failed, positive v0 daugh --> " << hfCandCascade.dcapostopv() << ", negative v0 daugh --> " << hfCandCascade.dcanegtopv() << " , cut --> " << d0K0Cut[ptBin];
      return false; // d0 of the K0s daughters
    }
    */

    if (std::abs(hfCandCascade.impactParameter1()) > cuts->get(ptBin, "d0V0")) {
      return false; // d0 of the v0
    }

    return true;
  }

  template <typename T>
  bool selectionStandardPID(const T& track)
  {
    TrackSelectorPID selectorProton = TrackSelectorPID(kProton);
    if (track.p() < pPidThreshold) {
      selectorProton.setRangeNSigmaTPC(-nSigmaTpcMaxLowP, nSigmaTpcMaxLowP);
      selectorProton.setRangeNSigmaTOF(-nSigmaTofMaxLowP, nSigmaTofMaxLowP);
      selectorProton.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxLowP, nSigmaTpcCombinedMaxLowP);
      selectorProton.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxLowP, nSigmaTofCombinedMaxLowP);
    } else {
      selectorProton.setRangeNSigmaTPC(-nSigmaTpcMaxHighP, nSigmaTpcMaxHighP);
      selectorProton.setRangeNSigmaTOF(-nSigmaTofMaxHighP, nSigmaTofMaxHighP);
      selectorProton.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxHighP, nSigmaTpcCombinedMaxHighP);
      selectorProton.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxHighP, nSigmaTofCombinedMaxHighP);
    }

    return selectorProton.getStatusTrackPIDTpcAndTof(track) == TrackSelectorPID::Status::PIDAccepted;
  }

  template <typename T>
  bool selectionBayesPID(const T& track)
  {
    if (!selectionStandardPID(track)) { // possibility to add some pre-selection before using Bayesian PID
      return false;
    }

    TrackSelectorPID selectorProton = TrackSelectorPID(kProton);
    if (track.p() < pPidThreshold) {
      selectorProton.setProbBayesMin(probBayesMinLowP);
    } else {
      selectorProton.setProbBayesMin(probBayesMinHighP);
    }

    return selectorProton.getStatusTrackBayesProbPID(track) == TrackSelectorPID::Status::PIDAccepted;
  }

  template <typename T, typename U>
  bool selectionMl(const T& hfCandCascade, const U& bach)
  {

    auto ptCand = hfCandCascade.pt();
    std::vector<float> inputFeatures = setInputFeatures(hfCandCascade, bach);
    std::vector<float> outputMl = {};

    bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCand, outputMl);

    for (int classMl = 0; classMl < nClassesMl; classMl++) {
      hModelScore[classMl]->Fill(outputMl[classMl]);
      hModelScoreVsPtCand[classMl]->Fill(outputMl[classMl], ptCand);
    }

    return isSelectedMl;
  }

  void processWithStandardPID(aod::HfCandCascade const& candidates, aod::BigTracksPID const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {                     // looping over cascade candidates
      const auto& bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor track

      statusLc = 0;

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)
      if (!selectionTopol(candidate)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      if (!selectionStandardPID(bach)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      if (applyMl && !selectionMl(candidate, bach)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      statusLc = 1;

      hfSelLcToK0sPCandidate(statusLc);
    }
  }
  PROCESS_SWITCH(HfCandidateSelectorLcToK0sP, processWithStandardPID, "Use standard PID for bachelor track", true);

  void processWithBayesPID(aod::HfCandCascade const& candidates, MyBigTracksBayes const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {                    // looping over cascade candidates
      const auto& bach = candidate.prong0_as<MyBigTracksBayes>(); // bachelor track

      statusLc = 0;

      if (!selectionTopol(candidate)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      if (!selectionBayesPID(bach)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      if (applyMl && !selectionMl(candidate, bach)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      statusLc = 1;

      hfSelLcToK0sPCandidate(statusLc);
    }
  }
  PROCESS_SWITCH(HfCandidateSelectorLcToK0sP, processWithBayesPID, "Use Bayesian PID for bachelor track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorLcToK0sP>(cfgc)};
}
