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

/// \file candidateSelectorDsToKKPi.cxx
/// \brief Ds± → K± K∓ π± selection task
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Universita and INFN Torino
/// \author Stefano Politano <stefano.politano@cern.ch>, Politecnico and INFN Torino

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseDsToKKPi.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Struct to extend TracksPid tables
struct HfCandidateSelectorDsToKKPiExpressions {
  Spawns<aod::TracksPidPiExt> rowTracksPidFullPi;
  Spawns<aod::TracksPidKaExt> rowTracksPidFullKa;
};

/// Struct for applying Ds to KKpi selection cuts
struct HfCandidateSelectorDsToKKPi {
  Produces<aod::HfSelDsToKKPi> hfSelDsToKKPiCandidate;
  Produces<aod::HfMlDsToKKPi> hfMlDsToKKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 1., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC"};
  //  TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_ds_to_k_k_pi::cuts[0], hf_cuts_ds_to_k_k_pi::nBinsPt, hf_cuts_ds_to_k_k_pi::nCutVars, hf_cuts_ds_to_k_k_pi::labelsPt, hf_cuts_ds_to_k_k_pi::labelsCutVar}, "Ds candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTDs"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_DsToKKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  HfHelper hfHelper;
  o2::analysis::HfMlResponseDsToKKPi<float> hfMlResponse;
  std::vector<float> outputMlDsToKKPi = {};
  std::vector<float> outputMlDsToPiKK = {};
  o2::ccdb::CcdbApi ccdbApi;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDca, aod::TracksPidPiExt, aod::TracksPidKaExt>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
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
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
    }
  }

  /// Candidate selections independent from the daugther-mass hypothesis
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T1>
  bool selection(const T1& candidate)
  {
    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    if (candpT < ptCandMin || candpT > ptCandMax) { // check that the candidate pT is within the analysis range
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
    if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "impact parameter XY")) {
      return false;
    }
    return true;
  }

  /// Candidate selections for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param trackKaon1 is the first track with the kaon hypothesis
  /// \param trackKaon2 is the second track with the kaon hypothesis
  /// \param trackPion is the track with the pion hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionKKPi(const T1& candidate, const T2& trackKaon1, const T2& trackKaon2, const T2& trackPion)
  {
    int pTBin = findBin(binsPt, candidate.pt());
    if (pTBin == -1) {
      return false;
    }

    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (std::abs(hfHelper.invMassDsToKKPi(candidate) - o2::constants::physics::MassDS) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (hfHelper.deltaMassPhiDsToKKPi(candidate) > cuts->get(pTBin, "deltaM Phi")) {
      return false;
    }
    if (std::abs(hfHelper.cos3PiKDsToKKPi(candidate)) < cuts->get(pTBin, "cos^3 theta_PiK")) {
      return false;
    }
    return true;
  }

  /// Candidate selections for the PiKK daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon1 is the first track with the kaon hypothesis
  /// \param trackKaon2 is the second track with the kaon hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionPiKK(const T1& candidate, const T2& trackPion, const T2& trackKaon1, const T2& trackKaon2)
  {
    int pTBin = findBin(binsPt, candidate.pt());
    if (pTBin == -1) {
      return false;
    }

    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (std::abs(hfHelper.invMassDsToPiKK(candidate) - o2::constants::physics::MassDS) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (hfHelper.deltaMassPhiDsToPiKK(candidate) > cuts->get(pTBin, "deltaM Phi")) {
      return false;
    }
    if (std::abs(hfHelper.cos3PiKDsToPiKK(candidate)) < cuts->get(pTBin, "cos^3 theta_PiK")) {
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
      auto statusDsToKKPi = 0;
      auto statusDsToPiKK = 0;

      outputMlDsToKKPi.clear();
      outputMlDsToPiKK.clear();

      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::DsToKKPi)) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMlDsToKKPi, outputMlDsToPiKK);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, candidate.pt());
        }
        continue;
      }
      SETBIT(statusDsToKKPi, aod::SelectionStep::RecoSkims);
      SETBIT(statusDsToPiKK, aod::SelectionStep::RecoSkims);
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoSkims, candidate.pt());
      }

      auto trackPos1 = candidate.prong0_as<TracksSel>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksSel>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksSel>(); // positive daughter (negative for the antiparticles)

      // topological selections
      if (!selection(candidate)) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMlDsToKKPi, outputMlDsToPiKK);
        }
        continue;
      }

      bool topolDsToKKPi = selectionKKPi(candidate, trackPos1, trackNeg, trackPos2);
      bool topolDsToPiKK = selectionPiKK(candidate, trackPos1, trackNeg, trackPos2);
      if (!topolDsToKKPi && !topolDsToPiKK) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMlDsToKKPi, outputMlDsToPiKK);
        }
        continue;
      }
      if (topolDsToKKPi) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoTopol);
      }
      if (topolDsToPiKK) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoTopol);
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoTopol, candidate.pt());
      }

      // track-level PID selection
      int pidTrackPos1Pion = selectorPion.statusTpcOrTof(trackPos1);
      int pidTrackPos1Kaon = selectorKaon.statusTpcOrTof(trackPos1);
      int pidTrackPos2Pion = selectorPion.statusTpcOrTof(trackPos2);
      int pidTrackPos2Kaon = selectorKaon.statusTpcOrTof(trackPos2);
      int pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg);

      bool pidDsToKKPi = !(pidTrackPos1Kaon == TrackSelectorPID::Rejected ||
                           pidTrackNegKaon == TrackSelectorPID::Rejected ||
                           pidTrackPos2Pion == TrackSelectorPID::Rejected);

      bool pidDsToPiKK = !(pidTrackPos1Pion == TrackSelectorPID::Rejected ||
                           pidTrackNegKaon == TrackSelectorPID::Rejected ||
                           pidTrackPos2Kaon == TrackSelectorPID::Rejected);

      if (!pidDsToKKPi && !pidDsToPiKK) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        if (applyMl) {
          hfMlDsToKKPiCandidate(outputMlDsToKKPi, outputMlDsToPiKK);
        }
        continue;
      }
      if (topolDsToKKPi && pidDsToKKPi) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoPID);
      }
      if (topolDsToPiKK && pidDsToPiKK) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoPID);
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, candidate.pt());
      }

      if (applyMl) {
        // ML selections
        bool isSelectedMlDsToKKPi = false;
        bool isSelectedMlDsToPiKK = false;

        if (topolDsToKKPi && pidDsToKKPi) {
          std::vector<float> inputFeaturesDsToKKPi = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2, true);
          isSelectedMlDsToKKPi = hfMlResponse.isSelectedMl(inputFeaturesDsToKKPi, candidate.pt(), outputMlDsToKKPi);
        }
        if (topolDsToPiKK && pidDsToPiKK) {
          std::vector<float> inputFeaturesDsToPiKK = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2, false);
          isSelectedMlDsToPiKK = hfMlResponse.isSelectedMl(inputFeaturesDsToPiKK, candidate.pt(), outputMlDsToPiKK);
        }

        hfMlDsToKKPiCandidate(outputMlDsToKKPi, outputMlDsToPiKK);

        if (!isSelectedMlDsToKKPi && !isSelectedMlDsToPiKK) {
          hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
          continue;
        }
        if (isSelectedMlDsToKKPi) {
          SETBIT(statusDsToKKPi, aod::SelectionStep::RecoMl);
        }
        if (isSelectedMlDsToPiKK) {
          SETBIT(statusDsToPiKK, aod::SelectionStep::RecoMl);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, candidate.pt());
        }
      }

      hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorDsToKKPiExpressions>(cfgc),
    adaptAnalysisTask<HfCandidateSelectorDsToKKPi>(cfgc)};
}
