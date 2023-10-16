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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsDebugLcToK0sP.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

// possible input features for ML
enum MLInputFeatures {
  numContrib = 0,
  posX,
  posY,
  posZ,
  xSecondaryVertex,
  ySecondaryVertex,
  zSecondaryVertex,
  errorDecayLength,
  errorDecayLengthXY,
  chi2PCA,
  rSecondaryVertex,
  decayLength,
  decayLengthXY,
  decayLengthNormalised,
  decayLengthXYNormalised,
  impactParameterNormalised0,
  ptProng0,
  pProng0,
  impactParameterNormalised1,
  ptProng1,
  pProng1,
  pxProng0,
  pyProng0,
  pzProng0,
  pxProng1,
  pyProng1,
  pzProng1,
  impactParameter0,
  impactParameter1,
  errorImpactParameter0,
  errorImpactParameter1,
  v0X,
  v0Y,
  v0Z,
  v0Radius,
  v0CosPA,
  v0MLambda,
  v0MAntiLambda,
  v0MK0Short,
  v0MGamma,
  v0CtK0Short,
  v0CtLambda,
  dcaV0Daughters,
  pxPos,
  pyPos,
  pzPos,
  ptV0Pos,
  dcaPosToPV,
  pxNeg,
  pyNeg,
  pzNeg,
  ptV0Neg,
  dcaNegToPV,
  nSigmaTPCPr0,
  nSigmaTOFPr0,
  m,
  pt,
  p,
  cpa,
  cpaXY,
  ct,
  eta,
  phi,
  y,
  e,
  NInputFeatures
};

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
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_k0s_p::cuts[0], hf_cuts_lc_to_k0s_p::nBinsPt, hf_cuts_lc_to_k0s_p::nCutVars, hf_cuts_lc_to_k0s_p::labelsPt, hf_cuts_lc_to_k0s_p::labelsCutVar}, "Lc candidate selection per pT bin"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  Configurable<std::vector<std::string>> inputFeaturesML{"inputFeaturesML", {""}, "List of input features for the ML model"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> modelPathsCCDB{"modelPathsCCDB", "EventFiltering/PWGHF/BDTLcToK0sP", "Path on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_D0ToKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  HfHelper hfHelper;
  TrackSelectorPr selectorProtonLowP;
  TrackSelectorPr selectorProtonHighP;

  o2::analysis::HfMlResponse<float> hfMlResponse;
  std::vector<bool> selectedInputFeatures{std::vector<bool>(MLInputFeatures::NInputFeatures, false)};

  o2::ccdb::CcdbApi ccdbApi;

  std::vector<std::shared_ptr<TH1>> hModelScore;
  std::vector<std::shared_ptr<TH2>> hModelScoreVsPtCand;

  using TracksSel = soa::Join<aod::TracksWExtra, aod::TracksPidPr>;
  using TracksSelBayes = soa::Join<TracksSel, aod::pidBayesPr>;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (!doprocessWithStandardPID && !doprocessWithBayesPID) {
      LOGF(fatal, "Neither processWithStandardPID nor processWithBayesPID enabled. Please choose one.");
    }
    if (doprocessWithStandardPID && doprocessWithBayesPID) {
      LOGF(fatal, "Cannot enable processWithStandardPID and processWithBayesPID at the same time. Please choose one.");
    }

    selectorProtonLowP.setRangeNSigmaTpc(-nSigmaTpcMaxLowP, nSigmaTpcMaxLowP);
    selectorProtonLowP.setRangeNSigmaTof(-nSigmaTofMaxLowP, nSigmaTofMaxLowP);
    selectorProtonLowP.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMaxLowP, nSigmaTpcCombinedMaxLowP);
    selectorProtonLowP.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMaxLowP, nSigmaTofCombinedMaxLowP);
    selectorProtonLowP.setProbBayesMin(probBayesMinLowP);

    selectorProtonHighP.setRangeNSigmaTpc(-nSigmaTpcMaxHighP, nSigmaTpcMaxHighP);
    selectorProtonHighP.setRangeNSigmaTof(-nSigmaTofMaxHighP, nSigmaTofMaxHighP);
    selectorProtonHighP.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMaxHighP, nSigmaTpcCombinedMaxHighP);
    selectorProtonHighP.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMaxHighP, nSigmaTofCombinedMaxHighP);
    selectorProtonHighP.setProbBayesMin(probBayesMinHighP);

    if (applyMl) {
      hfMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB.value, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
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
    std::map<MLInputFeatures, std::string> inputFeatureNames{
      {numContrib, "numContrib"},
      {posX, "posX"},
      {posY, "posY"},
      {posZ, "posZ"},
      {xSecondaryVertex, "xSecondaryVertex"},
      {ySecondaryVertex, "ySecondaryVertex"},
      {zSecondaryVertex, "zSecondaryVertex"},
      {errorDecayLength, "errorDecayLength"},
      {errorDecayLengthXY, "errorDecayLengthXY"},
      {chi2PCA, "chi2PCA"},
      {rSecondaryVertex, "rSecondaryVertex"},
      {decayLength, "decayLength"},
      {decayLengthXY, "decayLengthXY"},
      {decayLengthNormalised, "decayLengthNormalised"},
      {decayLengthXYNormalised, "decayLengthXYNormalised"},
      {impactParameterNormalised0, "impactParameterNormalised0"},
      {ptProng0, "ptProng0"},
      {pProng0, "pProng0"},
      {impactParameterNormalised1, "impactParameterNormalised1"},
      {ptProng1, "ptProng1"},
      {pProng1, "pProng1"},
      {pxProng0, "pxProng0"},
      {pyProng0, "pyProng0"},
      {pzProng0, "pzProng0"},
      {pxProng1, "pxProng1"},
      {pyProng1, "pyProng1"},
      {pzProng1, "pzProng1"},
      {impactParameter0, "impactParameter0"},
      {impactParameter1, "impactParameter1"},
      {errorImpactParameter0, "errorImpactParameter0"},
      {errorImpactParameter1, "errorImpactParameter1"},
      {v0X, "v0X"},
      {v0Y, "v0Y"},
      {v0Z, "v0Z"},
      {v0Radius, "v0Radius"},
      {v0CosPA, "v0CosPA"},
      {v0MLambda, "v0MLambda"},
      {v0MAntiLambda, "v0MAntiLambda"},
      {v0MK0Short, "v0MK0Short"},
      {v0MGamma, "v0MGamma"},
      {v0CtK0Short, "v0CtK0Short"},
      {v0CtLambda, "v0CtLambda"},
      {dcaV0Daughters, "dcaV0Daughters"},
      {pxPos, "pxPos"},
      {pyPos, "pyPos"},
      {pzPos, "pzPos"},
      {ptV0Pos, "ptV0Pos"},
      {dcaPosToPV, "dcaPosToPV"},
      {pxNeg, "pxNeg"},
      {pyNeg, "pyNeg"},
      {pzNeg, "pzNeg"},
      {ptV0Neg, "ptV0Neg"},
      {dcaNegToPV, "dcaNegToPV"},
      {nSigmaTPCPr0, "nSigmaTPCPr0"},
      {nSigmaTOFPr0, "nSigmaTOFPr0"},
      {m, "m"},
      {pt, "pt"},
      {p, "p"},
      {cpa, "cpa"},
      {cpaXY, "cpaXY"},
      {ct, "ct"},
      {eta, "eta"},
      {phi, "phi"},
      {y, "y"},
      {e, "e"}};

    // check for each possible input feature if it is included in the list of selected input features or not
    for (const auto& inputFeature : inputFeatureNames) {
      if (std::find(std::begin(inputFeaturesML.value), std::end(inputFeaturesML.value), inputFeature.second) != std::end(inputFeaturesML.value)) {
        selectedInputFeatures[inputFeature.first] = true;
        LOG(info) << "Included \'" << inputFeature.second << "\' in list of ML input features.";
      } else {
        selectedInputFeatures[inputFeature.first] = false;
      }
    }

    // check if all given input features are recongnized
    for (const auto& inputFeature : inputFeaturesML.value) {
      bool found = false;
      for (const auto& inputFeatureName : inputFeatureNames) {
        if (inputFeatureName.second == inputFeature) {
          found = true;
          break;
        }
      }
      if (!found) {
        LOG(fatal) << "Can not find \'" << inputFeature << "\' in list of possible ML input features.";
      }
    }
  }

  // fill only the selcted input features into the the vector of ML input features
  template <typename T, typename U>
  std::vector<float> setInputFeatures(const T& candidate, const U& bach)
  {
    std::vector<float> inputFeatures;

    if (selectedInputFeatures[MLInputFeatures::numContrib]) {
      inputFeatures.push_back(bach.collision().numContrib());
    }
    if (selectedInputFeatures[MLInputFeatures::posX]) {
      inputFeatures.push_back(candidate.posX());
    }
    if (selectedInputFeatures[MLInputFeatures::posY]) {
      inputFeatures.push_back(candidate.posY());
    }
    if (selectedInputFeatures[MLInputFeatures::posZ]) {
      inputFeatures.push_back(candidate.posZ());
    }
    if (selectedInputFeatures[MLInputFeatures::xSecondaryVertex]) {
      inputFeatures.push_back(candidate.xSecondaryVertex());
    }
    if (selectedInputFeatures[MLInputFeatures::ySecondaryVertex]) {
      inputFeatures.push_back(candidate.ySecondaryVertex());
    }
    if (selectedInputFeatures[MLInputFeatures::zSecondaryVertex]) {
      inputFeatures.push_back(candidate.zSecondaryVertex());
    }
    if (selectedInputFeatures[MLInputFeatures::errorDecayLength]) {
      inputFeatures.push_back(candidate.errorDecayLength());
    }
    if (selectedInputFeatures[MLInputFeatures::errorDecayLengthXY]) {
      inputFeatures.push_back(candidate.errorDecayLengthXY());
    }
    if (selectedInputFeatures[MLInputFeatures::chi2PCA]) {
      inputFeatures.push_back(candidate.chi2PCA());
    }
    if (selectedInputFeatures[MLInputFeatures::rSecondaryVertex]) {
      inputFeatures.push_back(candidate.rSecondaryVertex());
    }
    if (selectedInputFeatures[MLInputFeatures::decayLength]) {
      inputFeatures.push_back(candidate.decayLength());
    }
    if (selectedInputFeatures[MLInputFeatures::decayLengthXY]) {
      inputFeatures.push_back(candidate.decayLengthXY());
    }
    if (selectedInputFeatures[MLInputFeatures::decayLengthNormalised]) {
      inputFeatures.push_back(candidate.decayLengthNormalised());
    }
    if (selectedInputFeatures[MLInputFeatures::decayLengthXYNormalised]) {
      inputFeatures.push_back(candidate.decayLengthXYNormalised());
    }
    if (selectedInputFeatures[MLInputFeatures::impactParameterNormalised0]) {
      inputFeatures.push_back(candidate.impactParameterNormalised0());
    }
    if (selectedInputFeatures[MLInputFeatures::ptProng0]) {
      inputFeatures.push_back(candidate.ptProng0());
    }
    if (selectedInputFeatures[MLInputFeatures::pProng0]) {
      inputFeatures.push_back(RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()));
    }
    if (selectedInputFeatures[MLInputFeatures::impactParameterNormalised1]) {
      inputFeatures.push_back(candidate.impactParameterNormalised1());
    }
    if (selectedInputFeatures[MLInputFeatures::ptProng1]) {
      inputFeatures.push_back(candidate.ptProng1());
    }
    if (selectedInputFeatures[MLInputFeatures::pProng1]) {
      inputFeatures.push_back(RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()));
    }
    if (selectedInputFeatures[MLInputFeatures::pxProng0]) {
      inputFeatures.push_back(candidate.pxProng0());
    }
    if (selectedInputFeatures[MLInputFeatures::pyProng0]) {
      inputFeatures.push_back(candidate.pyProng0());
    }
    if (selectedInputFeatures[MLInputFeatures::pzProng0]) {
      inputFeatures.push_back(candidate.pzProng0());
    }
    if (selectedInputFeatures[MLInputFeatures::pxProng1]) {
      inputFeatures.push_back(candidate.pxProng1());
    }
    if (selectedInputFeatures[MLInputFeatures::pyProng1]) {
      inputFeatures.push_back(candidate.pyProng1());
    }
    if (selectedInputFeatures[MLInputFeatures::pzProng1]) {
      inputFeatures.push_back(candidate.pzProng1());
    }
    if (selectedInputFeatures[MLInputFeatures::errorImpactParameter0]) {
      inputFeatures.push_back(candidate.impactParameter0());
    }
    if (selectedInputFeatures[MLInputFeatures::impactParameter1]) {
      inputFeatures.push_back(candidate.impactParameter1());
    }
    if (selectedInputFeatures[MLInputFeatures::errorImpactParameter0]) {
      inputFeatures.push_back(candidate.errorImpactParameter0());
    }
    if (selectedInputFeatures[MLInputFeatures::errorImpactParameter1]) {
      inputFeatures.push_back(candidate.errorImpactParameter1());
    }
    if (selectedInputFeatures[MLInputFeatures::v0X]) {
      inputFeatures.push_back(candidate.v0x());
    }
    if (selectedInputFeatures[MLInputFeatures::v0Y]) {
      inputFeatures.push_back(candidate.v0y());
    }
    if (selectedInputFeatures[MLInputFeatures::v0Z]) {
      inputFeatures.push_back(candidate.v0z());
    }
    if (selectedInputFeatures[MLInputFeatures::v0Radius]) {
      inputFeatures.push_back(candidate.v0radius());
    }
    if (selectedInputFeatures[MLInputFeatures::v0CosPA]) {
      inputFeatures.push_back(candidate.v0cosPA());
    }
    if (selectedInputFeatures[MLInputFeatures::v0MLambda]) {
      inputFeatures.push_back(candidate.mLambda());
    }
    if (selectedInputFeatures[MLInputFeatures::v0MAntiLambda]) {
      inputFeatures.push_back(candidate.mAntiLambda());
    }
    if (selectedInputFeatures[MLInputFeatures::v0MK0Short]) {
      inputFeatures.push_back(candidate.mK0Short());
    }
    if (selectedInputFeatures[MLInputFeatures::v0MGamma]) {
      inputFeatures.push_back(candidate.mGamma());
    }
    if (selectedInputFeatures[MLInputFeatures::v0CtK0Short]) {
      inputFeatures.push_back(hfHelper.ctV0K0s(candidate));
    }
    if (selectedInputFeatures[MLInputFeatures::v0CtK0Short]) {
      inputFeatures.push_back(hfHelper.ctV0Lambda(candidate));
    }
    if (selectedInputFeatures[MLInputFeatures::dcaV0Daughters]) {
      inputFeatures.push_back(candidate.dcaV0daughters());
    }
    if (selectedInputFeatures[MLInputFeatures::pxPos]) {
      inputFeatures.push_back(candidate.pxpos());
    }
    if (selectedInputFeatures[MLInputFeatures::pyPos]) {
      inputFeatures.push_back(candidate.pypos());
    }
    if (selectedInputFeatures[MLInputFeatures::pzPos]) {
      inputFeatures.push_back(candidate.pzpos());
    }
    if (selectedInputFeatures[MLInputFeatures::ptV0Pos]) {
      inputFeatures.push_back(candidate.ptV0Pos());
    }
    if (selectedInputFeatures[MLInputFeatures::dcaPosToPV]) {
      inputFeatures.push_back(candidate.dcapostopv());
    }
    if (selectedInputFeatures[MLInputFeatures::pxNeg]) {
      inputFeatures.push_back(candidate.pxneg());
    }
    if (selectedInputFeatures[MLInputFeatures::pyNeg]) {
      inputFeatures.push_back(candidate.pyneg());
    }
    if (selectedInputFeatures[MLInputFeatures::pzNeg]) {
      inputFeatures.push_back(candidate.pzneg());
    }
    if (selectedInputFeatures[MLInputFeatures::ptV0Neg]) {
      inputFeatures.push_back(candidate.ptV0Neg());
    }
    if (selectedInputFeatures[MLInputFeatures::dcaNegToPV]) {
      inputFeatures.push_back(candidate.dcanegtopv());
    }
    if (selectedInputFeatures[MLInputFeatures::nSigmaTPCPr0]) {
      inputFeatures.push_back(bach.tpcNSigmaPr());
    }
    if (selectedInputFeatures[MLInputFeatures::nSigmaTOFPr0]) {
      inputFeatures.push_back(bach.tofNSigmaPr());
    }
    if (selectedInputFeatures[MLInputFeatures::m]) {
      inputFeatures.push_back(hfHelper.invMassLcToK0sP(candidate));
    }
    if (selectedInputFeatures[MLInputFeatures::pt]) {
      inputFeatures.push_back(candidate.pt());
    }
    if (selectedInputFeatures[MLInputFeatures::p]) {
      inputFeatures.push_back(candidate.p());
    }
    if (selectedInputFeatures[MLInputFeatures::cpa]) {
      inputFeatures.push_back(candidate.cpa());
    }
    if (selectedInputFeatures[MLInputFeatures::cpaXY]) {
      inputFeatures.push_back(candidate.cpaXY());
    }
    if (selectedInputFeatures[MLInputFeatures::ct]) {
      inputFeatures.push_back(hfHelper.ctLc(candidate));
    }
    if (selectedInputFeatures[MLInputFeatures::eta]) {
      inputFeatures.push_back(candidate.eta());
    }
    if (selectedInputFeatures[MLInputFeatures::phi]) {
      inputFeatures.push_back(candidate.phi());
    }
    if (selectedInputFeatures[MLInputFeatures::y]) {
      inputFeatures.push_back(hfHelper.yLc(candidate));
    }
    if (selectedInputFeatures[MLInputFeatures::e]) {
      inputFeatures.push_back(hfHelper.eLc(candidate));
    }

    return inputFeatures;
  }

  /// Conjugate independent topological cuts
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

    if (std::abs(hfCandCascade.mK0Short() - o2::analysis::pdg::MassK0Short) > cuts->get(ptBin, "mK0s")) {
      return false; // mass of the K0s
    }

    if ((std::abs(hfCandCascade.mLambda() - o2::analysis::pdg::MassLambda0) < cuts->get(ptBin, "mLambda")) || (std::abs(hfCandCascade.mAntiLambda() - o2::analysis::pdg::MassLambda0) < cuts->get(ptBin, "mLambda"))) {
      return false; // mass of the Lambda
    }

    if (std::abs(hfCandCascade.mGamma() - o2::analysis::pdg::MassGamma) < cuts->get(ptBin, "mGamma")) {
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
    if (track.p() < pPidThreshold) {
      return selectorProtonLowP.statusTpcAndTof(track) == TrackSelectorPID::Accepted;
    } else {
      return selectorProtonHighP.statusTpcAndTof(track) == TrackSelectorPID::Accepted;
    }
  }

  template <typename T>
  bool selectionBayesPID(const T& track)
  {
    if (!selectionStandardPID(track)) { // possibility to add some pre-selection before using Bayesian PID
      return false;
    }

    if (track.p() < pPidThreshold) {
      return selectorProtonLowP.statusBayesProb(track) == TrackSelectorPID::Accepted;
    } else {
      return selectorProtonHighP.statusBayesProb(track) == TrackSelectorPID::Accepted;
    }
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

  void processWithStandardPID(aod::HfCandCascade const& candidates,
                              TracksSel const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {                     // looping over cascade candidates
      const auto& bach = candidate.prong0_as<TracksSel>();         // bachelor track

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

  void processWithBayesPID(aod::HfCandCascade const& candidates,
                           TracksSelBayes const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {                    // looping over cascade candidates
      const auto& bach = candidate.prong0_as<TracksSelBayes>();   // bachelor track

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
