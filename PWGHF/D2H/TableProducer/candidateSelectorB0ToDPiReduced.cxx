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

/// \file candidateSelectorB0ToDPiReduced.cxx
/// \brief B0 → D- π+ candidate selector
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseB0ToDPi.h"
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

struct HfCandidateSelectorB0ToDPiReduced {
  Produces<aod::HfSelB0ToDPi> hfSelB0ToDPiCandidate; // table defined in CandidateSelectionTables.h
  Produces<aod::HfMlB0ToDPi> hfMlB0ToDPiCandidate;   // table defined in CandidateSelectionTables.h

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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_b0_to_d_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_b0_to_d_pi::Cuts[0], hf_cuts_b0_to_d_pi::NBinsPt, hf_cuts_b0_to_d_pi::NCutVars, hf_cuts_b0_to_d_pi::labelsPt, hf_cuts_b0_to_d_pi::labelsCutVar}, "B0 candidate selection per pT bin"};
  // D-meson ML cuts
  Configurable<std::vector<double>> binsPtDmesMl{"binsPtDmesMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "D-meson pT bin limits for ML cuts"};
  Configurable<LabeledArray<double>> cutsDmesMl{"cutsDmesMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsDmesCutScore}, "D-meson ML cuts per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // B0 ML inference
  Configurable<bool> applyB0Ml{"applyB0Ml", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtB0Ml{"binsPtB0Ml", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirB0Ml{"cutDirB0Ml", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsB0Ml{"cutsB0Ml", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesB0Ml{"nClassesB0Ml", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_B0/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_B0ToDPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  // variable that will store the value of selectionFlagD (defined in dataCreatorDplusPiReduced.cxx)
  int mySelectionFlagD = -1;

  o2::analysis::HfMlResponseB0ToDPi<float, true> hfMlResponse;
  float outputMlNotPreselected = -1.;
  std::vector<float> outputMl;
  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;

  using TracksBachPion = soa::Join<HfRedTracks, HfRedTracksPid>;
  using TracksSoftPions = soa::Join<aod::HfRedSoftPiBases, aod::HfRedSoftPiCov, aod::HfRedSoftPiPid>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    std::array<bool, 4> doprocess{doprocessSelectionDplusPi, doprocessSelectionDplusPiWithDmesMl, doprocessSelectionDstarPi, doprocessSelectionDstarPiWithDmesMl};
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

    if (applyB0Ml) {
      hfMlResponse.configure(binsPtB0Ml, cutsB0Ml, cutDirB0Ml, nClassesB0Ml);
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

  /// Utility function to retrieve the bach pion track
  /// from the B0 candidate in the D*-pi decay channel
  /// \param candidate is the B0 candidate
  /// \return bach pion track
  template <IsB0ToDstarPiChannel T1>
  auto getTrackBachPi(const T1& candidate)
  {
    return candidate.template prongBachPi_as<TracksBachPion>();
  }

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the B0 candidate
  /// \param prongBachPi is the candidate's bachelor pion prong
  /// \note this method is used for B0 → D*- π+ candidates with D meson ML scores
  template <bool WithDmesMl, IsB0ToDstarPiChannel T1, typename T2>
  auto getMlInputFeatures(const T1& candB0, const T2& prongBachPi)
  {
    auto prongSoftPi = candB0.template prongSoftPi_as<TracksSoftPions>();
    if constexpr (WithDmesMl) {
      return hfMlResponse.getInputFeaturesDStarPi<true>(candB0, prongBachPi, prongSoftPi);
    } else {
      return hfMlResponse.getInputFeaturesDStarPi<false>(candB0, prongBachPi, prongSoftPi);
    }
  }

  /// Utility function to retrieve the bach pion track
  /// from the B0 candidate in the D-pi decay channel
  /// \param candidate is the B0 candidate
  /// \return bach pion track
  template <typename T1>
  auto getTrackBachPi(const T1& candidate)
  {
    return candidate.template prong1_as<TracksBachPion>();
  }

  /// Method to get the input features vector needed for ML inference
  /// \param candB0 is the B0 candidate
  /// \param prongBachPi is the candidate's bachelor pion prong
  /// \note this method is used for B0 → D- π+ candidates with D meson ML scores
  template <bool WithDmesMl, typename T1, typename T2>
  auto getMlInputFeatures(const T1& candB0, const T2& prongBachPi)
  {
    if constexpr (WithDmesMl) {
      return hfMlResponse.getInputFeatures<true>(candB0, prongBachPi);
    } else {
      return hfMlResponse.getInputFeatures<false>(candB0, prongBachPi);
    }
  }

  /// Main function to perform B0 candidate creation
  /// \param withDmesMl is the flag to use the table with ML scores for the D- daughter (only possible if present in the derived data)
  /// \param hfCandsB0 B0 candidates
  /// \param pionTracks pion tracks
  /// \param configs config inherited from the Dpi data creator
  template <bool WithDmesMl, typename Cands>
  void runSelection(Cands const& hfCandsB0,
                    TracksBachPion const&,
                    HfCandB0Configs const& configs)
  {
    // get DplusPi creator configurable
    for (const auto& config : configs) {
      mySelectionFlagD = config.mySelectionFlagD();
    }

    for (const auto& hfCandB0 : hfCandsB0) {
      int statusB0ToDPi = 0;
      auto ptCandB0 = hfCandB0.pt();

      SETBIT(statusB0ToDPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusB0ToDPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandB0);
      }

      // topological cuts
      if (!HfHelper::selectionB0ToDPiTopol(hfCandB0, cuts, binsPt)) {
        hfSelB0ToDPiCandidate(statusB0ToDPi);
        if (applyB0Ml) {
          hfMlB0ToDPiCandidate(outputMlNotPreselected);
        }
        // LOGF(info, "B0 candidate selection failed at topology selection");
        continue;
      }

      if constexpr (WithDmesMl) { // we include it in the topological selections
        if (!HfHelper::selectionDmesMlScoresForBReduced(hfCandB0, cutsDmesMl, binsPtDmesMl)) {
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          if (applyB0Ml) {
            hfMlB0ToDPiCandidate(outputMlNotPreselected);
          }
          // LOGF(info, "B0 candidate selection failed at D-meson ML selection");
          continue;
        }
      }

      SETBIT(statusB0ToDPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusB0ToDPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandB0);
      }

      auto trackBachPi = getTrackBachPi(hfCandB0);
      if (pionPidMethod == PidMethod::TpcOrTof || pionPidMethod == PidMethod::TpcAndTof) {
        int pidTrackBachPi{TrackSelectorPID::Status::NotApplicable};
        if (pionPidMethod == PidMethod::TpcOrTof) {
          pidTrackBachPi = selectorPion.statusTpcOrTof(trackBachPi);
        } else if (pionPidMethod == PidMethod::TpcAndTof) {
          pidTrackBachPi = selectorPion.statusTpcAndTof(trackBachPi);
        }
        if (!HfHelper::selectionB0ToDPiPid(pidTrackBachPi, acceptPIDNotApplicable.value)) {
          // LOGF(info, "B0 candidate selection failed at PID selection");
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          if (applyB0Ml) {
            hfMlB0ToDPiCandidate(outputMlNotPreselected);
          }
          continue;
        }
        SETBIT(statusB0ToDPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusB0ToDPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandB0);
        }
      }
      if (applyB0Ml) {
        // B0 ML selections
        std::vector<float> inputFeatures = getMlInputFeatures<WithDmesMl>(hfCandB0, trackBachPi);
        bool const isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCandB0, outputMl);
        hfMlB0ToDPiCandidate(outputMl[1]); // storing ML score for signal class

        if (!isSelectedMl) {
          hfSelB0ToDPiCandidate(statusB0ToDPi);
          continue;
        }
        SETBIT(statusB0ToDPi, SelectionStep::RecoMl); // RecoML = 3 --> statusB0ToDPi = 15 if pionPidMethod, 11 otherwise
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoMl, ptCandB0);
        }
      }

      hfSelB0ToDPiCandidate(statusB0ToDPi);
      // LOGF(info, "B0 candidate selection passed all selections");
    }
  }

  void processSelectionDplusPi(HfRedCandB0 const& hfCandsB0,
                               TracksBachPion const& pionTracks,
                               HfCandB0Configs const& configs)
  {
    runSelection<false>(hfCandsB0, pionTracks, configs);
  } // processSelectionDplusPi

  PROCESS_SWITCH(HfCandidateSelectorB0ToDPiReduced, processSelectionDplusPi, "Process selection DplusPi without ML scores of D mesons", true);

  void processSelectionDplusPiWithDmesMl(soa::Join<HfRedCandB0, HfRedB0DpMls> const& hfCandsB0,
                                         TracksBachPion const& pionTracks,
                                         HfCandB0Configs const& configs)
  {
    runSelection<true>(hfCandsB0, pionTracks, configs);
  } // processSelectionDplusPiWithDmesMl

  PROCESS_SWITCH(HfCandidateSelectorB0ToDPiReduced, processSelectionDplusPiWithDmesMl, "Process selection DplusPi with ML scores of D mesons", false);

  void processSelectionDstarPi(HfRedCandB0DStar const& hfCandsB0,
                               TracksBachPion const& pionTracks,
                               HfCandB0Configs const& configs,
                               TracksSoftPions const& /*softPions*/)
  {
    runSelection<false>(hfCandsB0, pionTracks, configs);
  } // processSelectionDstarPi

  PROCESS_SWITCH(HfCandidateSelectorB0ToDPiReduced, processSelectionDstarPi, "Process selection DstarPi without ML scores of D mesons", false);

  void processSelectionDstarPiWithDmesMl(soa::Join<HfRedCandB0DStar, HfRedB0DpMls> const& hfCandsB0,
                                         TracksBachPion const& pionTracks,
                                         HfCandB0Configs const& configs,
                                         TracksSoftPions const& /*softPions*/)
  {
    runSelection<true>(hfCandsB0, pionTracks, configs);
  } // processSelectionDstarPiWithDmesMl

  PROCESS_SWITCH(HfCandidateSelectorB0ToDPiReduced, processSelectionDstarPiWithDmesMl, "Process selection DstarPi with ML scores of D mesons", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorB0ToDPiReduced>(cfgc)};
}
