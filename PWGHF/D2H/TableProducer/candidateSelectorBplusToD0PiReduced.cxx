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

/// \file candidateSelectorBplusToD0PiReduced.cxx
/// \brief B+ → D0bar π+ candidate selector
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari

#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseBplusToD0PiReduced.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorBplusToD0PiReduced {
  Produces<aod::HfSelBplusToD0Pi> hfSelBplusToD0PiCandidate; // table defined in CandidateSelectionTables.h
  Produces<aod::HfMlBplusToD0Pi> hfMlBplusToD0PiCandidate;   // table defined in CandidateSelectionTables.h

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Enable PID
  Configurable<int> pionPidMethod{"pionPidMethod", 1, "PID selection method for the bachelor pion (0: none, 1: TPC or TOF, 2: TPC and TOF)"};
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::NotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bplus_to_d0_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bplus_to_d0_pi::cuts[0], hf_cuts_bplus_to_d0_pi::nBinsPt, hf_cuts_bplus_to_d0_pi::nCutVars, hf_cuts_bplus_to_d0_pi::labelsPt, hf_cuts_bplus_to_d0_pi::labelsCutVar}, "B+ candidate selection per pT bin"};
  // D0-meson ML cuts
  Configurable<std::vector<double>> binsPtDmesMl{"binsPtDmesMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "D0-meson pT bin limits for ML cuts"};
  Configurable<LabeledArray<double>> cutsDmesMl{"cutsDmesMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsDmesCutScore}, "D0-meson ML cuts per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // B+ ML inference
  Configurable<bool> applyBplusMl{"applyBplusMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtBpMl{"binsPtBpMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirBpMl{"cutDirBpMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsBpMl{"cutsBpMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesBpMl{"nClassesBpMl", static_cast<int>(hf_cuts_ml::nCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_BPLUS/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_BPLUSToD0Pi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  // variable that will store the value of selectionFlagD (defined in dataCreatorD0PiReduced.cxx)
  int mySelectionFlagD0 = -1;
  int mySelectionFlagD0bar = -1;

  o2::analysis::HfMlResponseBplusToD0PiReduced<float> hfMlResponse;
  float outputMlNotPreselected = -1.;
  std::vector<float> outputMl = {};
  o2::ccdb::CcdbApi ccdbApi;

  HfHelper hfHelper;
  TrackSelectorPi selectorPion;

  HistogramRegistry registry{"registry"};

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;

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

    if (applyBplusMl) {
      hfMlResponse.configure(binsPtBpMl, cutsBpMl, cutDirBpMl, nClassesBpMl);
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

  /// Main function to perform B+ candidate creation
  /// \param withDmesMl is the flag to use the table with ML scores for the D- daughter (only possible if present in the derived data)
  /// \param hfCandsBp B+ candidates
  /// \param pionTracks pion tracks
  /// \param configs config inherited from the D0pi data creator
  template <bool withDmesMl, typename Cands>
  void runSelection(Cands const& hfCandsBp,
                    TracksPion const& /*pionTracks*/,
                    HfCandBpConfigs const& configs)
  {
    // get D0Pi creator configurable
    for (const auto& config : configs) {
      mySelectionFlagD0 = config.mySelectionFlagD0();
      mySelectionFlagD0bar = config.mySelectionFlagD0bar();
    }

    for (const auto& hfCandBp : hfCandsBp) {
      int statusBplus = 0;
      auto ptCandBplus = hfCandBp.pt();

      SETBIT(statusBplus, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusBplus = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandBplus);
      }

      // topological cuts
      if (!hfHelper.selectionBplusToD0PiTopol(hfCandBp, cuts, binsPt)) {
        hfSelBplusToD0PiCandidate(statusBplus);
        if (applyBplusMl) {
          hfMlBplusToD0PiCandidate(outputMlNotPreselected);
        }
        // LOGF(info, "B+ candidate selection failed at topology selection");
        continue;
      }

      if constexpr (withDmesMl) { // we include it in the topological selections
        if (!hfHelper.selectionDmesMlScoresForBReduced(hfCandBp, cutsDmesMl, binsPtDmesMl)) {
          hfSelBplusToD0PiCandidate(statusBplus);
          if (applyBplusMl) {
            hfMlBplusToD0PiCandidate(outputMlNotPreselected);
          }
          // LOGF(info, "B+ candidate selection failed at D0-meson ML selection");
          continue;
        }
      }

      SETBIT(statusBplus, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusBplus = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandBplus);
      }

      // track-level PID selection
      auto trackPi = hfCandBp.template prong1_as<TracksPion>();
      if (pionPidMethod) {
        int pidTrackPi{TrackSelectorPID::Status::NotApplicable};
        if (pionPidMethod == 1) {
          pidTrackPi = selectorPion.statusTpcOrTof(trackPi);
        } else {
          pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
        }
        if (!hfHelper.selectionBplusToD0PiPid(pidTrackPi, acceptPIDNotApplicable.value)) {
          // LOGF(info, "B+ candidate selection failed at PID selection");
          hfSelBplusToD0PiCandidate(statusBplus);
          if (applyBplusMl) {
            hfMlBplusToD0PiCandidate(outputMlNotPreselected);
          }
          continue;
        }
        SETBIT(statusBplus, SelectionStep::RecoPID); // RecoPID = 2 --> statusBplus = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandBplus);
        }
      }
      if (applyBplusMl) {
        // B+ ML selections
        std::vector<float> inputFeatures = hfMlResponse.getInputFeatures<withDmesMl>(hfCandBp, trackPi);
        bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCandBplus, outputMl);
        hfMlBplusToD0PiCandidate(outputMl[1]); // storing ML score for signal class

        if (!isSelectedMl) {
          hfSelBplusToD0PiCandidate(statusBplus);
          continue;
        }
        SETBIT(statusBplus, SelectionStep::RecoMl); // RecoML = 3 --> statusBplus = 15 if pionPidMethod, 11 otherwise
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoMl, ptCandBplus);
        }
      }

      hfSelBplusToD0PiCandidate(statusBplus);
      // LOGF(info, "B+ candidate selection passed all selections");
    }
  }

  void processSelection(HfRedCandBplus const& hfCandsBp,
                        TracksPion const& pionTracks,
                        HfCandBpConfigs const& configs)
  {
    runSelection<false>(hfCandsBp, pionTracks, configs);
  } // processSelection

  PROCESS_SWITCH(HfCandidateSelectorBplusToD0PiReduced, processSelection, "Process selection without ML scores of D mesons", true);

  void processSelectionWithDmesMl(soa::Join<HfRedCandBplus, HfRedBplusD0Mls> const& hfCandsBp,
                                  TracksPion const& pionTracks,
                                  HfCandBpConfigs const& configs)
  {
    runSelection<true>(hfCandsBp, pionTracks, configs);
  } // processSelectionWithDmesMl

  PROCESS_SWITCH(HfCandidateSelectorBplusToD0PiReduced, processSelectionWithDmesMl, "Process selection with ML scores of D mesons", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorBplusToD0PiReduced>(cfgc)};
}
