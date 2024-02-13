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

/// \file candidateSelectorBsToDsPi.cxx
/// \brief Bs → Ds- π+ candidate selector
/// \note adapted from candidateSelectorB0ToDPi.cxx
///
/// \author Phil Stahlhut <phil.lennart.stahlhut@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorBsToDsPi {
  Produces<aod::HfSelBsToDsPi> hfSelBsToDsPiCandidate; // table defined in CandidateSelectionTables.h
  Produces<aod::HfMlBsToDsPi> hfMlBsToDsPiCandidate;   // table defined in CandidateSelectionTables.h

  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bs_to_ds_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bs_to_ds_pi::cuts[0], hf_cuts_bs_to_ds_pi::nBinsPt, hf_cuts_bs_to_ds_pi::nCutVars, hf_cuts_bs_to_ds_pi::labelsPt, hf_cuts_bs_to_ds_pi::labelsCutVar}, "Bs candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTBs"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFilesCCDB", std::vector<std::string>{"ModelHandler_onnx_BsToDsPi.onnx"}, "ONNX file names on CCDB for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // check if selectionFlagDs (defined in candidateCreatorBs.cxx) and usePid configurables are in sync
  bool selectionFlagDsAndUsePidInSync = true;

  o2::analysis::HfMlResponse<float> hfMlResponse;
  std::vector<float> outputMl = {};

  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;
  HfHelper hfHelper;

  using TracksPidWithSel = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
    if (usePid) {
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
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }

    if (applyMl) {
      hfMlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.init();
      outputMl.assign(((std::vector<int>)cutDirMl).size(), -1.f); // dummy value for ML output
    }

    int selectionFlagDs = -1;
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name.compare("hf-candidate-creator-bs") == 0) {
        for (const auto& option : device.options) {
          if (option.name.compare("selectionFlagDs") == 0) {
            selectionFlagDs = option.defaultValue.get<int>();
            LOGF(info, "selectionFlagDs = %d", selectionFlagDs);
          }
        }
      }
    }

    if (usePid && !TESTBIT(selectionFlagDs, SelectionStep::RecoPID)) {
      selectionFlagDsAndUsePidInSync = false;
      LOG(warning) << "PID selections required on Bs daughters (usePid=true) but no PID selections on Ds candidates were required a priori (selectionFlagDs<7). Set selectionFlagDs=7 in hf-candidate-creator-bs";
    }
    if (!usePid && TESTBIT(selectionFlagDs, SelectionStep::RecoPID)) {
      selectionFlagDsAndUsePidInSync = false;
      LOG(warning) << "No PID selections required on Bs daughters (usePid=false) but PID selections on Ds candidates were required a priori (selectionFlagDs=7). Set selectionFlagDs<7 in hf-candidate-creator-bs";
    }
  }

  void process(aod::HfCandBs const& hfCandsBs,
               TracksPidWithSel const&)
  {
    for (const auto& hfCandBs : hfCandsBs) {
      int statusBsToDsPi = 0;
      auto ptCandBs = hfCandBs.pt();

      // check if flagged as Bs → Ds π
      if (!TESTBIT(hfCandBs.hfflag(), hf_cand_bs::DecayType::BsToDsPi)) {
        hfSelBsToDsPiCandidate(statusBsToDsPi);
        if (applyMl) {
          hfMlBsToDsPiCandidate(outputMl);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCandBs);
        }
        continue;
      }
      SETBIT(statusBsToDsPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusBsToDsPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandBs);
      }

      // topological cuts
      if (!hfHelper.selectionBsToDsPiTopol(hfCandBs, cuts, binsPt)) {
        hfSelBsToDsPiCandidate(statusBsToDsPi);
        if (applyMl) {
          hfMlBsToDsPiCandidate(outputMl);
        }
        continue;
      }
      SETBIT(statusBsToDsPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusBsToDsPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandBs);
      }

      // checking if selectionFlagDs and usePid are in sync
      if (!selectionFlagDsAndUsePidInSync) {
        hfSelBsToDsPiCandidate(statusBsToDsPi);
        if (applyMl) {
          hfMlBsToDsPiCandidate(outputMl);
        }
        continue;
      }
      // track-level PID selection
      if (usePid) {
        auto trackPi = hfCandBs.prong1_as<TracksPidWithSel>();
        int pidTrackPi = selectorPion.statusTpcAndTof(trackPi);
        if (!hfHelper.selectionBsToDsPiPid(pidTrackPi, acceptPIDNotApplicable.value)) {
          hfSelBsToDsPiCandidate(statusBsToDsPi);
          if (applyMl) {
            hfMlBsToDsPiCandidate(outputMl);
          }
          continue;
        }
        SETBIT(statusBsToDsPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusBsToDsPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandBs);
        }
      }

      // ML selections
      if (applyMl) {
        std::vector<float> inputFeatures{hfCandBs.cpa(),
                                         hfCandBs.cpaXY(),
                                         hfCandBs.decayLength(),
                                         hfCandBs.decayLengthXY(),
                                         hfCandBs.chi2PCA(),
                                         hfCandBs.impactParameter0(),
                                         hfCandBs.impactParameter1(),
                                         hfCandBs.maxNormalisedDeltaIP(),
                                         hfCandBs.impactParameterProduct()};

        bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCandBs, outputMl);
        hfMlBsToDsPiCandidate(outputMl);

        if (!isSelectedMl) {
          hfSelBsToDsPiCandidate(statusBsToDsPi);
          continue;
        }
        SETBIT(statusBsToDsPi, aod::SelectionStep::RecoMl);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, ptCandBs);
        }
      }

      hfSelBsToDsPiCandidate(statusBsToDsPi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorBsToDsPi>(cfgc)};
}
