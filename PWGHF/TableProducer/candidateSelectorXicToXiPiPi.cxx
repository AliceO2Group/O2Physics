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

/// \file candidateSelectorXicToXiPiPi.cxx
/// \brief Ξc± → Ξ∓ π± π± candidate selector
///
/// \author Phil Lennart Stahlhut <phil.lennart.stahlhut@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Utils/utilsAnalysis.h" // findBin function
#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorXicToXiPiPi {
  Produces<aod::HfSelXicToXiPiPi> hfSelXicToXiPiPiCandidate;
  Produces<aod::HfMlXicToXiPiPi> hfMlXicToXiPiPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
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
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_xi_pi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_xi_pi_pi::cuts[0], hf_cuts_xic_to_xi_pi_pi::nBinsPt, hf_cuts_xic_to_xi_pi_pi::nCutVars, hf_cuts_xic_to_xi_pi_pi::labelsPt, hf_cuts_xic_to_xi_pi_pi::labelsCutVar}, "Bs candidate selection per pT bin"};
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
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTXic"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFilesCCDB", std::vector<std::string>{"ModelHandler_onnx_XicTOXiPiPi.onnx"}, "ONNX file names on CCDB for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};


  o2::analysis::HfMlResponse<float> hfMlResponse;
  std::vector<float> outputMl = {};

  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;

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
  }

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T1>
  bool selectionTopol(const T1& hfCandXic)
  {
    auto candpT = hfCandXic.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    // check candidate mass is within a defined mass window
    if (std::abs(hfCandXic.invMassXic() - o2::constants::physics::MassXiCPlus) > cuts->get(pTBin, "m")) {
      return false;
    }

    // cosine of pointing angle
    if (hfCandXic.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    // cosine of pointing angle XY
    if (hfCandXic.cpaXY() <= cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }

    // candidate maximum decay length
    if (hfCandXic.decayLength() > cuts->get(pTBin, "max decay length")) {
      return false;
    }

    // candidate maximum decay length XY
    if (hfCandXic.decayLengthXY() > cuts->get(pTBin, "max decay length XY")) {
      return false;
    }

    // candidate chi2PC
    if (hfCandXic.chi2PCA() > cuts->get(pTBin, "chi2PCA")) {
      return false;
    }

    // maximum DCA of daughters
    if ((std::abs(hfCandXic.impactParameter0()) > cuts->get(pTBin, "max impParXY Xi")) ||
        (std::abs(hfCandXic.impactParameter1()) > cuts->get(pTBin, "max impParXY Pi0")) ||
        (std::abs(hfCandXic.impactParameter2()) > cuts->get(pTBin, "max impParXY Pi1"))) {
      return false;
    }

    // cut on daughter pT
    if (hfCandXic.ptProng0() < cuts->get(pTBin, "pT Xi")  ||
        hfCandXic.ptProng1() < cuts->get(pTBin, "pT Pi0") || 
        hfCandXic.ptProng2() < cuts->get(pTBin, "pT Pi1")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPi0 PID status of trackPi0 (prong1 of Xic candidate)
  /// \param pidTrackPi1 PID status of trackPi1 (prong2 of Xic candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prong1 and prong 2 of Xic candidate pass all selections
  template <typename T1 = int, typename T2 = int, typename T3 = bool>
  bool selectionPid(const T1& pidTrackPi0, const T2& pidTrackPi1, const T3& acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && pidTrackPi0 != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi0 == TrackSelectorPID::Rejected) {
      return false;
    }
    if (!acceptPIDNotApplicable && pidTrackPi1 != TrackSelectorPID::Accepted) {
      return false;
    }
    if (acceptPIDNotApplicable && pidTrackPi1 == TrackSelectorPID::Rejected) {
      return false;
    }

    return true;
  }

  void process(aod::HfCandXic const& hfCandsXic,
               TracksPidWithSel const&)
  {
    for (const auto& hfCandXic : hfCandsXic) {
      int statusXicToXiPiPi = 0;
      auto ptCandXic = hfCandXic.pt();

      // check if flagged as Ξc± → Ξ∓ π± π± 
      if (!TESTBIT(hfCandXic.hfflag(), hf_cand_xictoxipipi::DecayType::XicToXiPiPi)) {
        hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
        if (applyMl) {
          hfMlXicToXiPiPiCandidate(outputMl);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCandXic);
        }
        continue;
      }
      SETBIT(statusXicToXiPiPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusXicToXiPiPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandXic);
      }

      // topological cuts
      if (!selectionTopol(hfCandXic)) {
        hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
        if (applyMl) {
          hfMlXicToXiPiPiCandidate(outputMl);
        }
        continue;
      }
      SETBIT(statusXicToXiPiPi, SelectionStep::RecoTopol); // RecoTopol = 1 --> statusXicToXiPiPi = 3
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoTopol, ptCandXic);
      }

      // track-level PID selection
      if (usePid) {
        auto trackPi0 = hfCandXic.pi0_as<TracksPidWithSel>();
        auto trackPi1 = hfCandXic.pi1_as<TracksPidWithSel>();
        int pidTrackPi0 = selectorPion.statusTpcAndTof(trackPi0);
        int pidTrackPi1 = selectorPion.statusTpcAndTof(trackPi1);
        if (!selectionPid(pidTrackPi0, pidTrackPi1, acceptPIDNotApplicable.value)) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          if (applyMl) {
            hfMlXicToXiPiPiCandidate(outputMl);
          }
          continue;
        }
        SETBIT(statusXicToXiPiPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusXicToXiPiPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandXic);
        }
      }

      // ML selections
      if (applyMl) {
        std::vector<float> inputFeatures{};

        bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCandXic, outputMl);
        hfMlXicToXiPiPiCandidate(outputMl);

        if (!isSelectedMl) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          continue;
        }
        SETBIT(statusXicToXiPiPi, aod::SelectionStep::RecoMl);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, ptCandXic);
        }
      }

      hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorXicToXiPiPi>(cfgc)};
}
