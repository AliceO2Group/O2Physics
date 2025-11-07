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

/// \file candidateSelectorLbToLcPiReduced.cxx
/// \brief Lb → Lc+ π- candidate selector
///
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseLbToLcPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsPid.h"

#include "Common/Core/TrackSelectorPID.h"

#include <CCDB/CcdbApi.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;
using namespace o2::aod::pid_tpc_tof_utils;

struct HfCandidateSelectorLbToLcPiReduced {
  Produces<aod::HfSelLbToLcPi> hfSelLbToLcPiCandidate; // table defined in CandidateSelectionTables.h
  Produces<aod::HfMlLbToLcPi> hfMlLbToLcPiCandidate;   // table defined in CandidateSelectionTables.h

  Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<float> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Enable PID
  Configurable<int> pionPidMethod{"pionPidMethod", PidMethod::TpcOrTof, "PID selection method for the bachelor pion (PidMethod::NoPid: none, PidMethod::TpcOrTof: TPC or TOF, PidMethod::TpcAndTof: TPC and TOF)"};
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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lb_to_lc_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lb_to_lc_pi::Cuts[0], hf_cuts_lb_to_lc_pi::NBinsPt, hf_cuts_lb_to_lc_pi::NCutVars, hf_cuts_lb_to_lc_pi::labelsPt, hf_cuts_lb_to_lc_pi::labelsCutVar}, "Lb candidate selection per pT bin"};
  // Lc ML cuts
  Configurable<std::vector<double>> binsPtLcMl{"binsPtLcMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "Lc pT bin limits for ML cuts"};
  Configurable<LabeledArray<double>> cutsLcMl{"cutsLcMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsDmesCutScore}, "Lc ML cuts per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // Lb ML inference
  Configurable<bool> applyLbMl{"applyLbMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtLbMl{"binsPtLbMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirLbMl{"cutDirLbMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsLbMl{"cutsLbMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesLbMl{"nClassesLbMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_Lb/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_LbToLcPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponseLbToLcPi<float> hfMlResponse;
  float outputMlNotPreselected = -1.;
  std::vector<float> outputMl;
  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 2> doprocess{doprocessSelection, doprocessSelectionWithLcMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function for data should be enabled at a time.");
    }

    if (pionPidMethod < 0 || pionPidMethod >= PidMethod::NPidMethods) {
      LOGP(fatal, "Invalid PID option in configurable, please set 0 (no PID), 1 (TPC or TOF), or 2 (TPC and TOF)");
    }

    if (pionPidMethod != PidMethod::NoPid) {
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

    if (applyLbMl) {
      hfMlResponse.configure(binsPtLbMl, cutsLbMl, cutDirLbMl, nClassesLbMl);
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

  /// Main function to perform Lb candidate selection
  /// \param withLcMl is the flag to use the table with ML scores for the Ds- daughter (only possible if present in the derived data)
  /// \param hfCandsLb Lb candidates
  /// \param pionTracks pion tracks
  /// \param configs config inherited from the charm-hadron data creator
  template <bool WithLcMl, typename Cands>
  void runSelection(Cands const& hfCandsLb,
                    TracksPion const&,
                    HfCandLbConfigs const&)
  {
    for (const auto& hfCandLb : hfCandsLb) {
      int statusLbToLcPi = 0;
      outputMl.clear();
      auto ptCandLb = hfCandLb.pt();

      SETBIT(statusLbToLcPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusLbToLcPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandLb);
      }

      // topological cuts
      if (!HfHelper::selectionLbToLcPiTopol(hfCandLb, cuts, binsPt)) {
        hfSelLbToLcPiCandidate(statusLbToLcPi);
        if (applyLbMl) {
          hfMlLbToLcPiCandidate(outputMlNotPreselected);
        }
        continue;
      }

      if constexpr (WithLcMl) { // we include it in the topological selections
        if (!HfHelper::selectionDmesMlScoresForBReduced(hfCandLb, cutsLcMl, binsPtLcMl)) {
          hfSelLbToLcPiCandidate(statusLbToLcPi);
          if (applyLbMl) {
            hfMlLbToLcPiCandidate(outputMlNotPreselected);
          }
          continue;
        }
      }

      SETBIT(statusLbToLcPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusLbToLcPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandLb);
      }

      // track-level PID selection
      auto trackPi = hfCandLb.template prong1_as<TracksPion>();
      if (pionPidMethod == PidMethod::TpcOrTof || pionPidMethod == PidMethod::TpcAndTof) {
        int pidTrackPi{TrackSelectorPID::Status::NotApplicable};
        if (pionPidMethod == PidMethod::TpcOrTof) {
          pidTrackPi = selectorPion.statusTpcOrTof(trackPi);
        } else if (pionPidMethod == PidMethod::TpcAndTof) {
          pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
        }
        if (!HfHelper::selectionLbToLcPiPid(pidTrackPi, acceptPIDNotApplicable.value)) {
          hfSelLbToLcPiCandidate(statusLbToLcPi);
          if (applyLbMl) {
            hfMlLbToLcPiCandidate(outputMlNotPreselected);
          }
          continue;
        }
        SETBIT(statusLbToLcPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusLbToLcPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandLb);
        }
      }

      if (applyLbMl) {
        // Lb ML selections
        std::vector<float> inputFeatures = hfMlResponse.getInputFeatures<WithLcMl>(hfCandLb, trackPi);
        bool const isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCandLb, outputMl);
        hfMlLbToLcPiCandidate(outputMl[1]);

        if (!isSelectedMl) {
          hfSelLbToLcPiCandidate(statusLbToLcPi);
          continue;
        }
        SETBIT(statusLbToLcPi, SelectionStep::RecoMl); // RecoML = 3 --> statusLbToLcPi = 15 if PidMethod, 11 otherwise
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoMl, ptCandLb);
        }
      }

      hfSelLbToLcPiCandidate(statusLbToLcPi);
    }
  }

  void processSelection(HfRedCandLb const& hfCandsLb,
                        TracksPion const& pionTracks,
                        HfCandLbConfigs const& configs)
  {
    runSelection<false>(hfCandsLb, pionTracks, configs);
  } // processSelection

  PROCESS_SWITCH(HfCandidateSelectorLbToLcPiReduced, processSelection, "Process selection without ML scores of Lc", true);

  void processSelectionWithLcMl(soa::Join<HfRedCandLb, HfRedLbLcMls> const& hfCandsLb,
                                TracksPion const& pionTracks,
                                HfCandLbConfigs const& configs)
  {
    runSelection<true>(hfCandsLb, pionTracks, configs);
  } // processSelectionwithLcMl

  PROCESS_SWITCH(HfCandidateSelectorLbToLcPiReduced, processSelectionWithLcMl, "Process selection with ML scores of Lc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorLbToLcPiReduced>(cfgc)};
}
