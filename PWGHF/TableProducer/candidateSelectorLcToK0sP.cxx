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
///         Elisa Meninno, <elisa.meninno@cern.ch>, Vienna

#include <string>
#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/Core/HfMlResponseLcToK0sP.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

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
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_k0s_p::Cuts[0], hf_cuts_lc_to_k0s_p::NBinsPt, hf_cuts_lc_to_k0s_p::NCutVars, hf_cuts_lc_to_k0s_p::labelsPt, hf_cuts_lc_to_k0s_p::labelsCutVar}, "Lc candidate selection per pT bin"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTLc"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_LcToPK0S.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  HfHelper hfHelper;
  TrackSelectorPr selectorProtonLowP;
  TrackSelectorPr selectorProtonHighP;

  o2::analysis::HfMlResponseLcToK0sP<float> hfMlResponse;

  o2::ccdb::CcdbApi ccdbApi;

  std::vector<std::shared_ptr<TH1>> hModelScore;
  std::vector<std::shared_ptr<TH2>> hModelScoreVsPtCand;

  using TracksSel = soa::Join<aod::TracksWExtra, aod::TracksPidPr, aod::PidTpcTofFullPr>;
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
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();

      // load histograms for ML score
      AxisSpec axisScore = {100, 0, 1, "score"};
      AxisSpec axisBinsPt = {binsPtMl, "#it{p}_{T} (GeV/#it{c})"};
      for (int classMl = 0; classMl < nClassesMl; classMl++) {
        hModelScore.push_back(registry.add<TH1>(Form("hMlScoreClass%d", classMl), "Model score distribution for Lc;Model score;counts", HistType::kTH1F, {axisScore}));
        hModelScoreVsPtCand.push_back(registry.add<TH2>(Form("hMlScoreClass%dVsPtCand", classMl), "Model score distribution for Lc;Model score;counts", HistType::kTH2F, {axisScore, axisBinsPt}));
      }
    }
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

    if (std::abs(hfCandCascade.mK0Short() - o2::constants::physics::MassK0Short) > cuts->get(ptBin, "mK0s")) {
      return false; // mass of the K0s
    }

    if ((std::abs(hfCandCascade.mLambda() - o2::constants::physics::MassLambda0) < cuts->get(ptBin, "mLambda")) || (std::abs(hfCandCascade.mAntiLambda() - o2::constants::physics::MassLambda0) < cuts->get(ptBin, "mLambda"))) {
      return false; // mass of the Lambda
    }

    if (std::abs(hfCandCascade.mGamma() - o2::constants::physics::MassGamma) < cuts->get(ptBin, "mGamma")) {
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
    std::vector<float> inputFeatures = hfMlResponse.getInputFeatures(hfCandCascade, bach);
    std::vector<float> outputMl = {};

    bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCand, outputMl);

    for (int classMl = 0; classMl < nClassesMl; classMl++) {
      hModelScore[classMl]->Fill(outputMl[classMl]);
      hModelScoreVsPtCand[classMl]->Fill(outputMl[classMl], ptCand);
    }

    return isSelectedMl;
  }

  void processWithStandardPID(aod::HfCandCascade const& candidates,
                              TracksSel const&)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {             // looping over cascade candidates
      const auto& bach = candidate.prong0_as<TracksSel>(); // bachelor track

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
                           TracksSelBayes const&)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {                  // looping over cascade candidates
      const auto& bach = candidate.prong0_as<TracksSelBayes>(); // bachelor track

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
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorLcToK0sP>(cfgc)};
}
