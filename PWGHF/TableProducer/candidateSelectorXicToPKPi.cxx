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

/// \file candidateSelectorXicToPKPi.cxx
/// \brief Ξc± → p± K∓ π± selection task
/// \note Inspired from candidateSelectorLc.cxx
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN PADOVA
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Cristina Terrevoli <cristina.terrevoli@cern.ch>, INFN BARI

#include <string>
#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseXicToPKPi.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Struct for applying Xic selection cuts
struct HfCandidateSelectorXicToPKPi {
  Produces<aod::HfSelXicToPKPi> hfSelXicToPKPiCandidate;
  Produces<aod::HfMlXicToPKPi> hfMlXicToPKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> usePid{"usePid", true, "Bool to use or not the PID at filtering level"};
  // Combined PID options
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Bool to decide how to combine TPC and TOF PID: true =  both (if present, only one otherwise); false = one is enough"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<double> decayLengthXYNormalisedMin{"decayLengthXYNormalisedMin", 3., "Min. normalised decay length XY"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_p_k_pi::Cuts[0], hf_cuts_xic_to_p_k_pi::NBinsPt, hf_cuts_xic_to_p_k_pi::NCutVars, hf_cuts_xic_to_p_k_pi::labelsPt, hf_cuts_xic_to_p_k_pi::labelsCutVar}, "Xic candidate selection per pT bin"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTXic"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_XicToPKPi.onnx"}, "ONNX file names for each pT bin (if not from    CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", true, "Flag to enable QA histogram"};

  o2::analysis::HfMlResponseXicToPKPi<float> hfMlResponse;
  std::vector<float> outputMlXicToPKPi = {};
  std::vector<float> outputMlXicToPiKP = {};
  o2::ccdb::CcdbApi ccdbApi;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  TrackSelectorPr selectorProton;
  HfHelper hfHelper;

  using TracksSel = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;

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
    selectorProton = selectorPion;

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

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    // candidate chi2PCA
    if (candidate.chi2PCA() > cuts->get(pTBin, "chi2PCA")) {
      return false;
    }

    // candidate decay length
    if (candidate.decayLength() <= cuts->get(pTBin, "decay length")) {
      return false;
    }

    // candidate decay length XY
    if (candidate.decayLengthXY() <= cuts->get(pTBin, "decLengthXY")) {
      return false;
    }

    // candidate normalized decay length XY
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normDecLXY")) {
      return false;
    }

    // candidate normalised decay length (Inspired from Lc selector)
    if (candidate.decayLengthXYNormalised() < decayLengthXYNormalisedMin) {
      return false;
    }

    // candidate ct
    if (hfHelper.ctXic(candidate) > cuts->get(pTBin, "ct")) {
      return false;
    }

    // candidate impact parameter XY
    if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "impParXY")) {
      return false;
    }

    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackProton is the track with the proton hypothesis
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackProton, const T2& trackKaon, const T2& trackPion)
  {

    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // cut on daughter pT
    if (trackProton.pt() < cuts->get(pTBin, "pT p") || trackKaon.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    if (trackProton.globalIndex() == candidate.prong0Id()) {
      if (std::abs(hfHelper.invMassXicToPKPi(candidate) - o2::constants::physics::MassXiCPlus) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(hfHelper.invMassXicToPiKP(candidate) - o2::constants::physics::MassXiCPlus) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    return true;
  }

  void process(aod::HfCand3Prong const& candidates,
               TracksSel const&)
  {
    // looping over 3-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      auto statusXicToPKPi = 0;
      auto statusXicToPiKP = 0;

      outputMlXicToPKPi.clear();
      outputMlXicToPiKP.clear();

      auto ptCand = candidate.pt();

      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::XicToPKPi)) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        if (applyMl) {
          hfMlXicToPKPiCandidate(outputMlXicToPKPi, outputMlXicToPiKP);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCand);
        }
        continue;
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoSkims, ptCand);
      }
      SETBIT(statusXicToPKPi, aod::SelectionStep::RecoSkims);
      SETBIT(statusXicToPiKP, aod::SelectionStep::RecoSkims);

      auto trackPos1 = candidate.prong0_as<TracksSel>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksSel>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksSel>(); // positive daughter (negative for the antiparticles)

      // implement filter bit 4 cut - should be done before this task at the track selection level

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        if (applyMl) {
          hfMlXicToPKPiCandidate(outputMlXicToPKPi, outputMlXicToPiKP);
        }
        continue;
      }

      // conjugate-dependent topplogical selection for Xic

      bool topolXicToPKPi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolXicToPiKP = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolXicToPKPi && !topolXicToPiKP) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        if (applyMl) {
          hfMlXicToPKPiCandidate(outputMlXicToPKPi, outputMlXicToPiKP);
        }
        continue;
      }
      if (topolXicToPKPi) {
        SETBIT(statusXicToPKPi, aod::SelectionStep::RecoTopol);
      }
      if (topolXicToPiKP) {
        SETBIT(statusXicToPiKP, aod::SelectionStep::RecoTopol);
      }

      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoTopol, candidate.pt());
      }

      auto pidXicToPKPi = -1;
      auto pidXicToPiKP = -1;

      if (!usePid) {
        // PID non applied
        pidXicToPKPi = 1;
        pidXicToPiKP = 1;
      } else {
        // track-level PID selection
        TrackSelectorPID::Status pidTrackPos1Proton = TrackSelectorPID::Accepted;
        TrackSelectorPID::Status pidTrackPos2Proton = TrackSelectorPID::Accepted;
        TrackSelectorPID::Status pidTrackPos1Pion = TrackSelectorPID::Accepted;
        TrackSelectorPID::Status pidTrackPos2Pion = TrackSelectorPID::Accepted;
        TrackSelectorPID::Status pidTrackNegKaon = TrackSelectorPID::Accepted;
        if (usePidTpcAndTof) {

          pidTrackPos1Proton = selectorProton.statusTpcAndTof(trackPos1);
          pidTrackPos2Proton = selectorProton.statusTpcAndTof(trackPos2);
          pidTrackPos1Pion = selectorPion.statusTpcAndTof(trackPos1);
          pidTrackPos2Pion = selectorPion.statusTpcAndTof(trackPos2);
          pidTrackNegKaon = selectorKaon.statusTpcAndTof(trackNeg);
        } else {
          pidTrackPos1Proton = selectorProton.statusTpcOrTof(trackPos1);
          pidTrackPos2Proton = selectorProton.statusTpcOrTof(trackPos2);
          pidTrackPos1Pion = selectorPion.statusTpcOrTof(trackPos1);
          pidTrackPos2Pion = selectorPion.statusTpcOrTof(trackPos2);
          pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg);
        }

        if (pidTrackPos1Proton == TrackSelectorPID::Accepted &&
            pidTrackNegKaon == TrackSelectorPID::Accepted &&
            pidTrackPos2Pion == TrackSelectorPID::Accepted) {
          pidXicToPKPi = 1; // accept XicpKpi
        } else if (pidTrackPos1Proton == TrackSelectorPID::Rejected ||
                   pidTrackNegKaon == TrackSelectorPID::Rejected ||
                   pidTrackPos2Pion == TrackSelectorPID::Rejected) {
          pidXicToPKPi = 0; // exclude XicpKpi
        }
        if (pidTrackPos2Proton == TrackSelectorPID::Accepted &&
            pidTrackNegKaon == TrackSelectorPID::Accepted &&
            pidTrackPos1Pion == TrackSelectorPID::Accepted) {
          pidXicToPiKP = 1; // accept XicpiKp
        } else if (pidTrackPos1Pion == TrackSelectorPID::Rejected ||
                   pidTrackNegKaon == TrackSelectorPID::Rejected ||
                   pidTrackPos2Proton == TrackSelectorPID::Rejected) {
          pidXicToPiKP = 0; // exclude XicpiKp
        }
      }

      if (pidXicToPKPi == 0 && pidXicToPiKP == 0) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        if (applyMl) {
          hfMlXicToPKPiCandidate(outputMlXicToPKPi, outputMlXicToPiKP);
        }
        continue;
      }

      if ((pidXicToPKPi == -1 || pidXicToPKPi == 1) && topolXicToPKPi) {
        SETBIT(statusXicToPKPi, aod::SelectionStep::RecoPID);
      }
      if ((pidXicToPiKP == -1 || pidXicToPiKP == 1) && topolXicToPiKP) {
        SETBIT(statusXicToPiKP, aod::SelectionStep::RecoPID);
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, candidate.pt());
      }

      if (applyMl) {
        // ML selections
        bool isSelectedMlXicToPKPi = false;
        bool isSelectedMlXicToPiKP = false;

        if (topolXicToPKPi && pidXicToPKPi) {
          std::vector<float> inputFeaturesXicToPKPi = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2, true);
          isSelectedMlXicToPKPi = hfMlResponse.isSelectedMl(inputFeaturesXicToPKPi, ptCand, outputMlXicToPKPi);
        }
        if (topolXicToPiKP && pidXicToPiKP) {
          std::vector<float> inputFeaturesXicToPiKP = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2, false);
          isSelectedMlXicToPiKP = hfMlResponse.isSelectedMl(inputFeaturesXicToPiKP, ptCand, outputMlXicToPiKP);
        }

        hfMlXicToPKPiCandidate(outputMlXicToPKPi, outputMlXicToPiKP);

        if (!isSelectedMlXicToPKPi && !isSelectedMlXicToPiKP) {
          hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
          continue;
        }

        if (isSelectedMlXicToPKPi) {
          SETBIT(statusXicToPKPi, aod::SelectionStep::RecoMl);
        }
        if (isSelectedMlXicToPiKP) {
          SETBIT(statusXicToPiKP, aod::SelectionStep::RecoMl);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, candidate.pt());
        }
      }

      hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorXicToPKPi>(cfgc)};
}
