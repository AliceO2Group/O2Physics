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

/// \file candidateSelectorBsToDsPiReduced.cxx
/// \brief Bs → Ds- π+ candidate selector
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseBsToDsPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorBsToDsPiReduced {
  Produces<aod::HfSelBsToDsPi> hfSelBsToDsPiCandidate; // table defined in CandidateSelectionTables.h
  Produces<aod::HfMlBsToDsPi> hfMlBsToDsPiCandidate;   // table defined in CandidateSelectionTables.h

  Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Enable PID
  Configurable<int> pionPidMethod{"pionPidMethod", 1, "PID selection method for the bachelor pion (0: none, 1: TPC or TOF, 2: TPC and TOF)"};
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::NotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  // TPC PID
  Configurable<float> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<float> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<float> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<float> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<float> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<float> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<float> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<float> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bs_to_ds_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bs_to_ds_pi::cuts[0], hf_cuts_bs_to_ds_pi::nBinsPt, hf_cuts_bs_to_ds_pi::nCutVars, hf_cuts_bs_to_ds_pi::labelsPt, hf_cuts_bs_to_ds_pi::labelsCutVar}, "Bs candidate selection per pT bin"};
  // D-meson ML cuts
  Configurable<std::vector<double>> binsPtDmesMl{"binsPtDmesMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "D-meson pT bin limits for ML cuts"};
  Configurable<LabeledArray<double>> cutsDmesMl{"cutsDmesMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsDmesCutScore}, "D-meson ML cuts per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // B0 ML inference
  Configurable<bool> applyBsMl{"applyBsMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtBsMl{"binsPtBsMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirBsMl{"cutDirBsMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsBsMl{"cutsBsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesBsMl{"nClassesBsMl", static_cast<int>(hf_cuts_ml::nCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_Bs/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_BsToDsPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponseBsToDsPi<float> hfMlResponse;
  std::vector<float> outputMl = {};
  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;
  HfHelper hfHelper;

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 2> doprocess{doprocessSelection, doprocessSelectionWithDmesMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function for data should be enabled at a time.");
    }

    if (pionPidMethod < 0 || pionPidMethod > 2) {
      LOGP(fatal, "Invalid PID option in configurable, please set 0 (no PID), 1 (TPC or TOF), or 2 (TPC and TOF)");
    }

    if (pionPidMethod) {
      selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
      selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
      selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
      selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
      selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
      selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    }

    if (activateQA) {
      constexpr int kNBinsSelections = 1 + SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      labels[1 + aod::SelectionStep::RecoMl] = "ML selection";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }

    if (applyBsMl) {
      hfMlResponse.configure(binsPtBsMl, cutsBsMl, cutDirBsMl, nClassesBsMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
    }
  }

  /// Main function to perform Bs candidate selection
  /// \param withDmesMl is the flag to use the table with ML scores for the Ds- daughter (only possible if present in the derived data)
  /// \param hfCandsBs Bs candidates
  /// \param pionTracks pion tracks
  /// \param configs config inherited from the charm-hadron data creator
  template <bool withDmesMl, typename Cands>
  void runSelection(Cands const& hfCandsBs,
                    TracksPion const&,
                    HfCandBsConfigs const&)
  {
    for (const auto& hfCandBs : hfCandsBs) {
      int statusBsToDsPi = 0;
      outputMl.clear();
      auto ptCandBs = hfCandBs.pt();

      SETBIT(statusBsToDsPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusBsToDsPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandBs);
      }

      // topological cuts
      if (!hfHelper.selectionBsToDsPiTopol(hfCandBs, cuts, binsPt)) {
        hfSelBsToDsPiCandidate(statusBsToDsPi);
        if (applyBsMl) {
          hfMlBsToDsPiCandidate(outputMl);
        }
        continue;
      }

      if constexpr (withDmesMl) { // we include it in the topological selections
        if (!hfHelper.selectionDmesMlScoresForBReduced(hfCandBs, cutsDmesMl, binsPtDmesMl)) {
          hfSelBsToDsPiCandidate(statusBsToDsPi);
          if (applyBsMl) {
            hfMlBsToDsPiCandidate(outputMl);
          }
          continue;
        }
      }

      SETBIT(statusBsToDsPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusBsToDsPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandBs);
      }

      // track-level PID selection
      auto trackPi = hfCandBs.template prong1_as<TracksPion>();
      if (pionPidMethod) {
        int pidTrackPi{TrackSelectorPID::Status::NotApplicable};
        if (pionPidMethod == 1) {
          pidTrackPi = selectorPion.statusTpcOrTof(trackPi);
        } else {
          pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
        }
        if (!hfHelper.selectionBsToDsPiPid(pidTrackPi, acceptPIDNotApplicable.value)) {
          hfSelBsToDsPiCandidate(statusBsToDsPi);
          if (applyBsMl) {
            hfMlBsToDsPiCandidate(outputMl);
          }
          continue;
        }
        SETBIT(statusBsToDsPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusBsToDsPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandBs);
        }
      }

      if (applyBsMl) {
        // Bs ML selections
        std::vector<float> inputFeatures = hfMlResponse.getInputFeatures<withDmesMl>(hfCandBs, trackPi);
        bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCandBs, outputMl);
        hfMlBsToDsPiCandidate(outputMl);

        if (!isSelectedMl) {
          hfSelBsToDsPiCandidate(statusBsToDsPi);
          continue;
        }
        SETBIT(statusBsToDsPi, SelectionStep::RecoMl); // RecoML = 3 --> statusBsToDsPi = 15 if pionPidMethod, 11 otherwise
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoMl, ptCandBs);
        }
      }

      hfSelBsToDsPiCandidate(statusBsToDsPi);
    }
  }

  void processSelection(HfRedCandBs const& hfCandsBs,
                        TracksPion const& pionTracks,
                        HfCandBsConfigs const& configs)
  {
    runSelection<false>(hfCandsBs, pionTracks, configs);
  } // processSelection

  PROCESS_SWITCH(HfCandidateSelectorBsToDsPiReduced, processSelection, "Process selection without ML scores of D mesons", true);

  void processSelectionWithDmesMl(soa::Join<HfRedCandBs, HfRedBsDsMls> const& hfCandsBs,
                                  TracksPion const& pionTracks,
                                  HfCandBsConfigs const& configs)
  {
    runSelection<true>(hfCandsBs, pionTracks, configs);
  } // processSelectionWithDmesMl

  PROCESS_SWITCH(HfCandidateSelectorBsToDsPiReduced, processSelectionWithDmesMl, "Process selection with ML scores of D mesons", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorBsToDsPiReduced>(cfgc)};
}
