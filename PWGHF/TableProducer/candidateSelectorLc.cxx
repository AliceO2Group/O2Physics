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

/// \file candidateSelectorLc.cxx
/// \brief Λc± → p± K∓ π± selection task
///
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, University and INFN SALERNO
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Grazia Luparello  <grazia.luparello@cern.ch>, INFN Trieste

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseLcToPKPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Struct to extend TracksPid tables
struct HfCandidateSelectorLcExpressions {
  Spawns<aod::TracksPidPrExt> rowTracksPidFullPr;
  Spawns<aod::TracksPidKaExt> rowTracksPidFullKa;
  Spawns<aod::TracksPidPiExt> rowTracksPidFullPi;
};

/// Struct for applying Lc selection cuts
struct HfCandidateSelectorLc {
  Produces<aod::HfSelLc> hfSelLcCandidate;
  Produces<aod::HfMlLcToPKPi> hfMlLcToPKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> usePid{"usePid", true, "Bool to use or not the PID based on nSigma cut at filtering level"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.1, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 1., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.5, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 2.5, "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // Bayesian PID
  Configurable<bool> usePidBayes{"usePidBayes", true, "Bool to use or not the PID based on Bayesian probability cut at filtering level"};
  Configurable<double> ptPidBayesMin{"ptPidBayesMin", 0., "Lower bound of track pT for Bayesian PID"};
  Configurable<double> ptPidBayesMax{"ptPidBayesMax", 100, "Upper bound of track pT for Bayesian PID"};
  // Combined PID options
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Bool to decide how to combine TPC and TOF PID: true = both (if present, only one otherwise); false = one is enough"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_p_k_pi::cuts[0], hf_cuts_lc_to_p_k_pi::nBinsPt, hf_cuts_lc_to_p_k_pi::nCutVars, hf_cuts_lc_to_p_k_pi::labelsPt, hf_cuts_lc_to_p_k_pi::labelsCutVar}, "Lc candidate selection per pT bin"};
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
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTLc"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_LcToPKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  HfHelper hfHelper;
  o2::analysis::HfMlResponseLcToPKPi<float> hfMlResponse;
  std::vector<float> outputMlLcToPKPi = {};
  std::vector<float> outputMlLcToPiKP = {};
  o2::ccdb::CcdbApi ccdbApi;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  TrackSelectorPr selectorProton;

  using TracksSel = soa::Join<aod::TracksWExtra,
                              aod::TracksPidPiExt, aod::TracksPidKaExt, aod::TracksPidPrExt,
                              aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidBayes>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorPion.setRangePtBayes(ptPidBayesMin, ptPidBayesMax);
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
    if (candidate.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      return false;
    }

    if (candidate.decayLength() <= cuts->get(pTBin, "decay length")) {
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
      if (std::abs(hfHelper.invMassLcToPKPi(candidate) - o2::constants::physics::MassLambdaCPlus) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(hfHelper.invMassLcToPiKP(candidate) - o2::constants::physics::MassLambdaCPlus) > cuts->get(pTBin, "m")) {
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

      // final selection flag
      auto statusLcToPKPi = 0;
      auto statusLcToPiKP = 0;

      outputMlLcToPKPi.clear();
      outputMlLcToPiKP.clear();

      auto ptCand = candidate.pt();

      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        if (applyMl) {
          hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);
        }
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCand);
        }
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoSkims, ptCand);
      }

      auto trackPos1 = candidate.prong0_as<TracksSel>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TracksSel>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TracksSel>(); // positive daughter (negative for the antiparticles)

      // implement filter bit 4 cut - should be done before this task at the track selection level

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        if (applyMl) {
          hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);
        }
        continue;
      }

      // conjugate-dependent topological selection for Lc

      bool topolLcToPKPi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcToPiKP = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolLcToPKPi && !topolLcToPiKP) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        if (applyMl) {
          hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);
        }
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoTopol, candidate.pt());
      }

      auto pidLcToPKPi = -1;
      auto pidLcToPiKP = -1;
      auto pidBayesLcToPKPi = -1;
      auto pidBayesLcToPiKP = -1;

      if (!usePid) {
        // PID non applied
        pidLcToPKPi = 1;
        pidLcToPiKP = 1;
      } else {
        // track-level PID selection
        int pidTrackPos1Proton = 999;
        int pidTrackPos2Proton = 999;
        int pidTrackPos1Pion = 999;
        int pidTrackPos2Pion = 999;
        int pidTrackNegKaon = 999;
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
          pidLcToPKPi = 1; // accept LcToPKPi
        } else if (pidTrackPos1Proton == TrackSelectorPID::Rejected ||
                   pidTrackNegKaon == TrackSelectorPID::Rejected ||
                   pidTrackPos2Pion == TrackSelectorPID::Rejected) {
          pidLcToPKPi = 0; // exclude LcToPKPi
        }
        if (pidTrackPos2Proton == TrackSelectorPID::Accepted &&
            pidTrackNegKaon == TrackSelectorPID::Accepted &&
            pidTrackPos1Pion == TrackSelectorPID::Accepted) {
          pidLcToPiKP = 1; // accept LcToPiKP
        } else if (pidTrackPos1Pion == TrackSelectorPID::Rejected ||
                   pidTrackNegKaon == TrackSelectorPID::Rejected ||
                   pidTrackPos2Proton == TrackSelectorPID::Rejected) {
          pidLcToPiKP = 0; // exclude LcToPiKP
        }
      }

      if (!usePidBayes) {
        // PID non applied
        pidBayesLcToPKPi = 1;
        pidBayesLcToPiKP = 1;
      } else {
        int pidBayesTrackPos1Proton = selectorProton.statusBayes(trackPos1);
        int pidBayesTrackPos2Proton = selectorProton.statusBayes(trackPos2);
        int pidBayesTrackPos1Pion = selectorPion.statusBayes(trackPos1);
        int pidBayesTrackPos2Pion = selectorPion.statusBayes(trackPos2);
        int pidBayesTrackNegKaon = selectorKaon.statusBayes(trackNeg);

        if (pidBayesTrackPos1Proton == TrackSelectorPID::Accepted &&
            pidBayesTrackNegKaon == TrackSelectorPID::Accepted &&
            pidBayesTrackPos2Pion == TrackSelectorPID::Accepted) {
          pidBayesLcToPKPi = 1; // accept LcToPKPi
        } else if (pidBayesTrackPos1Proton == TrackSelectorPID::Rejected ||
                   pidBayesTrackNegKaon == TrackSelectorPID::Rejected ||
                   pidBayesTrackPos2Pion == TrackSelectorPID::Rejected) {
          pidBayesLcToPKPi = 0; // exclude LcToPKPi
        }
        if (pidBayesTrackPos2Proton == TrackSelectorPID::Accepted &&
            pidBayesTrackNegKaon == TrackSelectorPID::Accepted &&
            pidBayesTrackPos1Pion == TrackSelectorPID::Accepted) {
          pidBayesLcToPiKP = 1; // accept LcToPiKP
        } else if (pidBayesTrackPos1Pion == TrackSelectorPID::Rejected ||
                   pidBayesTrackNegKaon == TrackSelectorPID::Rejected ||
                   pidBayesTrackPos2Proton == TrackSelectorPID::Rejected) {
          pidBayesLcToPiKP = 0; // exclude LcToPiKP
        }
      }

      if (pidLcToPKPi == 0 && pidLcToPiKP == 0) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        if (applyMl) {
          hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);
        }
        continue;
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, candidate.pt());
      }

      if (pidBayesLcToPKPi == 0 && pidBayesLcToPiKP == 0) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }

      bool isSelectedMlLcToPKPi = true;
      bool isSelectedMlLcToPiKP = true;
      if (applyMl) {
        // ML selections
        isSelectedMlLcToPKPi = false;
        isSelectedMlLcToPiKP = false;

        if ((pidLcToPKPi == -1 || pidLcToPKPi == 1) && (pidBayesLcToPKPi == -1 || pidBayesLcToPKPi == 1) && topolLcToPKPi) {
          std::vector<float> inputFeaturesLcToPKPi = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2);
          isSelectedMlLcToPKPi = hfMlResponse.isSelectedMl(inputFeaturesLcToPKPi, candidate.pt(), outputMlLcToPKPi);
        }
        if ((pidLcToPiKP == -1 || pidLcToPiKP == 1) && (pidBayesLcToPiKP == -1 || pidBayesLcToPiKP == 1) && topolLcToPiKP) {
          std::vector<float> inputFeaturesLcToPiKP = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2);
          isSelectedMlLcToPiKP = hfMlResponse.isSelectedMl(inputFeaturesLcToPiKP, candidate.pt(), outputMlLcToPiKP);
        }

        hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);

        if (!isSelectedMlLcToPKPi && !isSelectedMlLcToPiKP) {
          hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
          continue;
        }

        if (activateQA) {
          registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, candidate.pt());
        }
      }

      if ((pidLcToPKPi == -1 || pidLcToPKPi == 1) && (pidBayesLcToPKPi == -1 || pidBayesLcToPKPi == 1) && isSelectedMlLcToPKPi && topolLcToPKPi) {
        statusLcToPKPi = 1; // identified as LcToPKPi
      }
      if ((pidLcToPiKP == -1 || pidLcToPiKP == 1) && (pidBayesLcToPiKP == -1 || pidBayesLcToPiKP == 1) && isSelectedMlLcToPiKP && topolLcToPiKP) {
        statusLcToPiKP = 1; // identified as LcToPiKP
      }

      hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorLcExpressions>(cfgc),
    adaptAnalysisTask<HfCandidateSelectorLc>(cfgc)};
}
