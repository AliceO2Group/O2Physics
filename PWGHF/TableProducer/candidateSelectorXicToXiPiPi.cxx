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
/// \author Phil Lennart Stahlhut <phil.lennart.stahlhut@cern.ch>, Heidelberg University
/// \author Jaeyoon Cho <jaeyoon.cho@cern.ch>, Inha University

#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfMlResponseXicToXiPiPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h" // findBin function

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorXicToXiPiPi {
  Produces<aod::HfSelXicToXiPiPi> hfSelXicToXiPiPiCandidate;
  Produces<aod::HfMlXicToXiPiPi> hfMlXicToXiPiPiCandidate;

  Configurable<float> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<float> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_xi_pi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_xi_pi_pi::cuts[0], hf_cuts_xic_to_xi_pi_pi::nBinsPt, hf_cuts_xic_to_xi_pi_pi::nCutVars, hf_cuts_xic_to_xi_pi_pi::labelsPt, hf_cuts_xic_to_xi_pi_pi::labelsCutVar}, "Xicplus candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
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
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTXicToXiPiPi"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_XicToXiPiPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponseXicToXiPiPi<float> hfMlResponse;
  std::vector<float> outputMlXicToXiPiPi = {};
  o2::ccdb::CcdbApi ccdbApi;
  TrackSelectorPi selectorPion;
  TrackSelectorPr selectorProton;

  using TracksPidWithSel = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    if (usePid) {
      // pion
      selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
      selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
      selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
      selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
      selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
      selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
      // proton
      selectorProton.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
      selectorProton.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
      selectorProton.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
      selectorProton.setRangePtTof(ptPidTofMin, ptPidTofMax);
      selectorProton.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
      selectorProton.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    }

    if (activateQA) {
      constexpr int kNBinsSelections = 1 + SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      labels[1 + SelectionStep::RecoMl] = "Skims & Topological & PID & ML selections";
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
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
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
    if (std::abs(hfCandXic.invMassXicPlus() - o2::constants::physics::MassXiCPlus) > cuts->get(pTBin, "m")) {
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
    if (hfCandXic.ptProng0() < cuts->get(pTBin, "pT Xi") ||
        hfCandXic.ptProng1() < cuts->get(pTBin, "pT Pi0") ||
        hfCandXic.ptProng2() < cuts->get(pTBin, "pT Pi1")) {
      return false;
    }

    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPi0   PID status of trackPi0 (prong1 of Xic candidate)
  /// \param pidTrackPi1   PID status of trackPi1 (prong2 of Xic candidate)
  /// \param pidTrackPr    PID status of trackPr (positive daughter of V0 candidate)
  /// \param pidTrackPiLam PID status of trackPiLam (negative daughter of V0 candidate)
  /// \param pidTrackPiXi  PID status of trackPiXi (Bachelor of cascade candidate)
  /// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
  /// \return true if prongs of Xic candidate pass all selections
  bool selectionPid(TrackSelectorPID::Status const pidTrackPi0,
                    TrackSelectorPID::Status const pidTrackPi1,
                    TrackSelectorPID::Status const pidTrackPr,
                    TrackSelectorPID::Status const pidTrackPiLam,
                    TrackSelectorPID::Status const pidTrackPiXi,
                    bool const acceptPIDNotApplicable)
  {
    if (!acceptPIDNotApplicable && (pidTrackPi0 != TrackSelectorPID::Accepted || pidTrackPi1 != TrackSelectorPID::Accepted || pidTrackPr != TrackSelectorPID::Accepted || pidTrackPiLam != TrackSelectorPID::Accepted || pidTrackPiXi != TrackSelectorPID::Accepted)) {
      return false;
    }
    if (acceptPIDNotApplicable && (pidTrackPi0 == TrackSelectorPID::Rejected || pidTrackPi1 == TrackSelectorPID::Rejected || pidTrackPr == TrackSelectorPID::Rejected || pidTrackPiLam == TrackSelectorPID::Rejected || pidTrackPiXi == TrackSelectorPID::Rejected)) {
      return false;
    }
    return true;
  }

  void process(aod::HfCandXic const& hfCandsXic,
               TracksPidWithSel const&)
  {
    for (const auto& hfCandXic : hfCandsXic) {
      int statusXicToXiPiPi = 0;

      outputMlXicToXiPiPi.clear();

      auto ptCandXic = hfCandXic.pt();

      if (activateQA) {
        registry.fill(HIST("hSelections"), 1, ptCandXic);
      }

      // No hfFlag -> by default skim selected
      SETBIT(statusXicToXiPiPi, SelectionStep::RecoSkims); // RecoSkims = 0 --> statusXicToXiPiPi = 1
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoSkims, ptCandXic);
      }

      // topological cuts
      if (!selectionTopol(hfCandXic)) {
        hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
        if (applyMl) {
          hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
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
        auto trackV0PosDau = hfCandXic.posTrack_as<TracksPidWithSel>();
        auto trackV0NegDau = hfCandXic.negTrack_as<TracksPidWithSel>();
        auto trackPiFromXi = hfCandXic.bachelor_as<TracksPidWithSel>();
        // assign proton and pion hypothesis to V0 daughters
        auto trackPr = trackV0PosDau;
        auto trackPiFromLam = trackV0NegDau;
        if (hfCandXic.sign() < 0) {
          trackPr = trackV0NegDau;
          trackPiFromLam = trackV0PosDau;
        }
        // PID info
        TrackSelectorPID::Status pidTrackPi0 = selectorPion.statusTpcAndTof(trackPi0);
        TrackSelectorPID::Status pidTrackPi1 = selectorPion.statusTpcAndTof(trackPi1);
        TrackSelectorPID::Status pidTrackPr = selectorProton.statusTpcAndTof(trackPr);
        TrackSelectorPID::Status pidTrackPiLam = selectorPion.statusTpcAndTof(trackPiFromLam);
        TrackSelectorPID::Status pidTrackPiXi = selectorPion.statusTpcAndTof(trackPiFromXi);

        if (!selectionPid(pidTrackPi0, pidTrackPi1, pidTrackPr, pidTrackPiLam, pidTrackPiXi, acceptPIDNotApplicable.value)) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          if (applyMl) {
            hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
          }
          continue;
        }
        SETBIT(statusXicToXiPiPi, SelectionStep::RecoPID); // RecoPID = 2 --> statusXicToXiPiPi = 7
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoPID, ptCandXic);
        }
      }

      // ML selection
      if (applyMl) {
        bool isSelectedMlXicToXiPiPi = false;
        std::vector<float> inputFeaturesXicToXiPiPi = hfMlResponse.getInputFeatures(hfCandXic);

        isSelectedMlXicToXiPiPi = hfMlResponse.isSelectedMl(inputFeaturesXicToXiPiPi, ptCandXic, outputMlXicToXiPiPi);

        hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);

        if (!isSelectedMlXicToXiPiPi) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          continue;
        }
        SETBIT(statusXicToXiPiPi, aod::SelectionStep::RecoMl);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + SelectionStep::RecoMl, ptCandXic);
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
