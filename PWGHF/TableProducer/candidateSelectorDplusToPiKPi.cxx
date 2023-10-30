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

/// \file candidateSelectorDplusToPiKPi.cxx
/// \brief D± → π± K∓ π± selection task
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponse.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Struct to extend TracksPid tables
struct HfCandidateSelectorDplusToPiKPiExpressions {
  Spawns<aod::TracksPidPiExt> rowTracksPidFullPi;
  Spawns<aod::TracksPidKaExt> rowTracksPidFullKa;
};

/// Struct for applying Dplus to piKpi selection cuts
struct HfCandidateSelectorDplusToPiKPi {
  Produces<aod::HfSelDplusToPiKPi> hfSelDplusToPiKPiCandidate;
  Produces<aod::HfMlDplusToPiKPi> hfMlDplusToPiKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 1., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  // PID option
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::NotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 3., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 3., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_dplus_to_pi_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_dplus_to_pi_k_pi::cuts[0], hf_cuts_dplus_to_pi_k_pi::nBinsPt, hf_cuts_dplus_to_pi_k_pi::nCutVars, hf_cuts_dplus_to_pi_k_pi::labelsPt, hf_cuts_dplus_to_pi_k_pi::labelsCutVar}, "Dplus candidate selection per pT bin"};
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
  Configurable<std::string> modelPathsCCDB{"modelPathsCCDB", "EventFiltering/PWGHF/BDTD0", "Path on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_D0ToKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  o2::analysis::HfMlResponse<float> hfMlResponse;
  std::vector<float> outputMl = {};
  o2::ccdb::CcdbApi ccdbApi;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  HfHelper hfHelper;

  using TracksSel = soa::Join<aod::TracksWExtra, aod::TracksPidPiExt, aod::TracksPidKaExt>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorKaon = selectorPion;

    if (activateQA) {
      constexpr int kNBinsSelections = 1 + aod::SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + aod::SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + aod::SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + aod::SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      labels[1 + aod::SelectionStep::RecoMl] = "ML selection";
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
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB.value, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.init();
      outputMl.assign(((std::vector<int>)cutDirMl).size(), -1.f); // dummy value for ML output
    }
  }

  /// Candidate selections
  /// \param candidate is candidate
  /// \param trackPion1 is the first track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \param trackPion2 is the second track with the pion hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selection(const T1& candidate, const T2& trackPion1, const T2& trackKaon, const T2& trackPion2)
  {
    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }
    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT > ptCandMax) {
      return false;
    }
    // cut on daughter pT
    if (trackPion1.pt() < cuts->get(pTBin, "pT Pi") || trackKaon.pt() < cuts->get(pTBin, "pT K") || trackPion2.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    // invariant-mass cut
    if (std::abs(hfHelper.invMassDplusToPiKPi(candidate) - o2::analysis::pdg::MassDPlus) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (candidate.decayLength() < cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }
    if (std::abs(candidate.maxNormalisedDeltaIP()) > cuts->get(pTBin, "max normalized deltaIP")) {
      return false;
    }
    return true;
  }

  /// Apply PID selection
  /// \param pidTrackPos1Pion is the PID status of trackPos1 (prong0 of D candidate)
  /// \param pidTrackNegKaon is the PID status of trackNeg (prong1 of D candidate)
  /// \param pidTrackPos2Pion is the PID status of trackPos2 (prong2 of D candidate)
  /// \return true if prongs pass all selections
  template <typename T = int>
  bool selectionPID(const T& pidTrackPos1Pion, const T& pidTrackNegKaon, const T& pidTrackPos2Pion)
  {
    if (!acceptPIDNotApplicable &&
        (pidTrackPos1Pion != TrackSelectorPID::Accepted ||
         pidTrackNegKaon != TrackSelectorPID::Accepted ||
         pidTrackPos2Pion != TrackSelectorPID::Accepted)) {
      return false;
    }
    if (acceptPIDNotApplicable &&
        (pidTrackPos1Pion == TrackSelectorPID::Rejected ||
         pidTrackNegKaon == TrackSelectorPID::Rejected ||
         pidTrackPos2Pion == TrackSelectorPID::Rejected)) {
      return false;
    }

    return true;
  }

  void process(aod::HfCand3Prong const& candidates,
               TracksSel const&)
  {
    // looping over 3-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag:
      auto statusDplusToPiKPi = 0;

      auto ptCand = candidate.pt();

      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::DplusToPiKPi)) {
        hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
        if (applyMl) {
          hfMlDplusToPiKPiCandidate(outputMl);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCand);
        }
        continue;
      }
      SETBIT(statusDplusToPiKPi, aod::SelectionStep::RecoSkims);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoSkims, ptCand);
      }

      auto trackPos1 = candidate.prong0_as<TracksSel>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksSel>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksSel>(); // positive daughter (negative for the antiparticles)

      // topological selection
      if (!selection(candidate, trackPos1, trackNeg, trackPos2)) {
        hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
        if (applyMl) {
          hfMlDplusToPiKPiCandidate(outputMl);
        }
        continue;
      }
      SETBIT(statusDplusToPiKPi, aod::SelectionStep::RecoTopol);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoTopol, ptCand);
      }

      // track-level PID selection
      int pidTrackPos1Pion = selectorPion.statusTpcAndTof(trackPos1);
      int pidTrackNegKaon = selectorKaon.statusTpcAndTof(trackNeg);
      int pidTrackPos2Pion = selectorPion.statusTpcAndTof(trackPos2);

      if (!selectionPID(pidTrackPos1Pion, pidTrackNegKaon, pidTrackPos2Pion)) { // exclude D±
        hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
        if (applyMl) {
          hfMlDplusToPiKPiCandidate(outputMl);
        }
        continue;
      }
      SETBIT(statusDplusToPiKPi, aod::SelectionStep::RecoPID);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, ptCand);
      }

      if (applyMl) {
        // ML selections
        std::vector<float> inputFeatures{candidate.ptProng0(),
                                         candidate.impactParameter0(),
                                         candidate.ptProng1(),
                                         candidate.impactParameter1(),
                                         candidate.ptProng2(),
                                         candidate.impactParameter2(),
                                         candidate.decayLength(),
                                         candidate.decayLengthXYNormalised(),
                                         candidate.cpa(),
                                         candidate.cpaXY(),
                                         candidate.maxNormalisedDeltaIP(),
                                         trackPos1.tpcTofNSigmaPi(),
                                         trackPos1.tpcTofNSigmaKa(),
                                         trackNeg.tpcTofNSigmaPi(),
                                         trackNeg.tpcTofNSigmaKa(),
                                         trackPos2.tpcTofNSigmaPi(),
                                         trackPos2.tpcTofNSigmaKa()};

        bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCand, outputMl);
        hfMlDplusToPiKPiCandidate(outputMl);

        if (!isSelectedMl) {
          hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
          continue;
        }
        SETBIT(statusDplusToPiKPi, aod::SelectionStep::RecoMl);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, ptCand);
        }
      }

      hfSelDplusToPiKPiCandidate(statusDplusToPiKPi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorDplusToPiKPiExpressions>(cfgc),
    adaptAnalysisTask<HfCandidateSelectorDplusToPiKPi>(cfgc)};
}
