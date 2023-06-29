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

/// \file candidateSelectorLcMl.cxx
/// \brief Λc± → p± K∓ π± selection task using BDT
///
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, University and INFN SALERNO
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Maja Kabus <maja.kabus@cern.ch>, CERN, Warsaw University of Technology

#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/trackUtilities.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Tools/ML/model.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_lc_to_p_k_pi;
using namespace hf_cuts_bdt_multiclass;
using namespace o2::ml;

/// Struct for applying Lc selection cuts
struct HfCandidateSelectorLcMl {
  Produces<aod::HfSelLc> hfSelLcCandidate;

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

  Configurable<double> cpaMin{"cpaMin", 0.95, "Lower bound of candidate CPA"};
  Configurable<double> maxDeltaMass{"maxDeltaMass", 0.5, "Max difference of inv mass compared to PDG Lc mass"};

  // ONNX BDT
  Configurable<bool> applyML{"applyML", false, "Flag to enable or disable ML application"};
  Configurable<std::string> onnxFileLcToPiKPConf{"onnxFileLcToPiKPConf", "/cvmfs/alice.cern.ch/data/analysis/2022/vAN-20220818/PWGHF/o2/trigger/ModelHandler_onnx_LcToPKPi.onnx", "ONNX file for ML model for Lc+ candidates"};
  Configurable<LabeledArray<double>> thresholdBDTScoreLcToPiKP{"thresholdBDTScoreLcToPiKP", {hf_cuts_bdt_multiclass::cuts[0], hf_cuts_bdt_multiclass::nBinsPt, hf_cuts_bdt_multiclass::nCutBdtScores, hf_cuts_bdt_multiclass::labelsPt, hf_cuts_bdt_multiclass::labelsCutBdt}, "Threshold values for BDT output scores of Lc+ candidates"};

  o2::ccdb::CcdbApi ccdbApi;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> mlModelPathCCDB{"mlModelPathCCDB", "Analysis/PWGHF/ML/HFTrigger/Lc", "Path on CCDB"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB. Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  Configurable<int> activateQA{"activateQA", 0, "flag to enable QA histos (0 no QA, 1 basic QA, 2 extended QA)"};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  int dataTypeML;
  OnnxModel model;

  using TrksPID = soa::Join<aod::BigTracksPIDExtended, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidBayes>;

  void init(o2::framework::InitContext&)
  {
    AxisSpec bdtAxis{100, 0.f, 1.f};
    if (applyML && activateQA != 0) {
      registry.add<TH1>("hLcBDTScoreBkg", "BDT background score distribution for Lc;BDT background score;counts", HistType::kTH1F, {bdtAxis});
      registry.add<TH1>("hLcBDTScorePrompt", "BDT prompt score distribution for Lc;BDT prompt score;counts", HistType::kTH1F, {bdtAxis});
      registry.add<TH1>("hLcBDTScoreNonPrompt", "BDT nonprompt score distribution for Lc;BDT nonprompt score;counts", HistType::kTH1F, {bdtAxis});
    }

    ccdbApi.init(url);

    // init ONNX runtime session
    std::map<std::string, std::string> metadata;
    std::map<std::string, std::string> headers;
    bool retrieveSuccess = true;
    if (applyML) {
      if (onnxFileLcToPiKPConf.value == "") {
        LOG(error) << "Apply ML specified, but no name given to the local model file";
      }
      if (loadModelsFromCCDB && timestampCCDB > 0) {
        retrieveSuccess = ccdbApi.retrieveBlob(mlModelPathCCDB.value, ".", metadata, timestampCCDB.value, false, onnxFileLcToPiKPConf.value);
        headers = ccdbApi.retrieveHeaders(mlModelPathCCDB.value, metadata, timestampCCDB.value);
        model.initModel(onnxFileLcToPiKPConf.value, false, 1, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
      } else if (!loadModelsFromCCDB) {
        model.initModel(onnxFileLcToPiKPConf.value, false, 1);
      } else {
        LOG(error) << "Retrieving model based on current run number not implemented yet... But anyways it is just a single model!";
      }
      if (retrieveSuccess) {
        auto session = model.getSession();
        auto inputShapes = session->GetInputShapes();
        if (inputShapes[0][0] < 0) {
          LOGF(warning, "Model for Lc with negative input shape likely because converted with hummingbird, setting it to 1.");
          inputShapes[0][0] = 1;
        }
        std::vector<float> dummyInput(model.getNumInputNodes(), 1.);
        model.evalModel(dummyInput); // Init the model evaluations
        dataTypeML = session->GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetElementType();
      } else {
        LOG(fatal) << "Error encountered while fetching/loading the ML model from CCDB! Maybe the ML model doesn't exist yet for this runnumber/timestamp?";
      }
    }
  }

  /*
  /// Selection on goodness of daughter tracks
  /// \note should be applied at candidate selection
  /// \param track is daughter track
  /// \return true if track is good
  template <typename T>
  bool daughterSelection(const T& track)
  {
    if (track.tpcNClsFound() == 0) {
      return false; //is it clusters findable or found - need to check
    }
    return true;
  }
  */

  void
    process(aod::HfCand3Prong const& candidates, TrksPID const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorPion.setRangePtBayes(ptPidBayesMin, ptPidBayesMax);

    TrackSelectorPID selectorKaon(selectorPion);
    selectorKaon.setPDG(kKPlus);

    TrackSelectorPID selectorProton(selectorPion);
    selectorProton.setPDG(kProton);

    // looping over 3-prong candidates
    for (auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      auto statusLcToPKPi = 0;
      auto statusLcToPiKP = 0;

      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }

      auto trackPos1 = candidate.prong0_as<TrksPID>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<TrksPID>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<TrksPID>(); // positive daughter (negative for the antiparticles)

      /*
      // daughter track validity selection
      if (!daughterSelection(trackPos1) || !daughterSelection(trackNeg) || !daughterSelection(trackPos2)) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }
      */

      // implement filter bit 4 cut - should be done before this task at the track selection level

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
        int pidTrackPos1Proton = selectorProton.getStatusTrackPIDTpcOrTof(trackPos1);
        int pidTrackPos2Proton = selectorProton.getStatusTrackPIDTpcOrTof(trackPos2);
        int pidTrackPos1Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos1);
        int pidTrackPos2Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos2);
        int pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackNeg);

        if (pidTrackPos1Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidLcToPKPi = 1; // accept LcToPKPi
        } else if (pidTrackPos1Proton == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos2Pion == TrackSelectorPID::Status::PIDRejected) {
          pidLcToPKPi = 0; // exclude LcToPKPi
        }
        if (pidTrackPos2Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidLcToPiKP = 1; // accept LcToPiKP
        } else if (pidTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos2Proton == TrackSelectorPID::Status::PIDRejected) {
          pidLcToPiKP = 0; // exclude LcToPiKP
        }
      }

      if (!usePidBayes) {
        // PID non applied
        pidBayesLcToPKPi = 1;
        pidBayesLcToPiKP = 1;
      } else {
        int pidBayesTrackPos1Proton = selectorProton.getStatusTrackBayesPID(trackPos1);
        int pidBayesTrackPos2Proton = selectorProton.getStatusTrackBayesPID(trackPos2);
        int pidBayesTrackPos1Pion = selectorPion.getStatusTrackBayesPID(trackPos1);
        int pidBayesTrackPos2Pion = selectorPion.getStatusTrackBayesPID(trackPos2);
        int pidBayesTrackNegKaon = selectorKaon.getStatusTrackBayesPID(trackNeg);

        if (pidBayesTrackPos1Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidBayesLcToPKPi = 1; // accept LcToPKPi
        } else if (pidBayesTrackPos1Proton == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackPos2Pion == TrackSelectorPID::Status::PIDRejected) {
          pidBayesLcToPKPi = 0; // exclude LcToPKPi
        }
        if (pidBayesTrackPos2Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidBayesLcToPiKP = 1; // accept LcToPiKP
        } else if (pidBayesTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackPos2Proton == TrackSelectorPID::Status::PIDRejected) {
          pidBayesLcToPiKP = 0; // exclude LcToPiKP
        }
      }

      if (pidLcToPKPi == 0 && pidLcToPiKP == 0) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }

      if (pidBayesLcToPKPi == 0 && pidBayesLcToPiKP == 0) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }

      if ((pidLcToPKPi == -1 || pidLcToPKPi == 1) && (pidBayesLcToPKPi == -1 || pidBayesLcToPKPi == 1)) {
        statusLcToPKPi = 1; // identified as LcToPKPi
      }
      if ((pidLcToPiKP == -1 || pidLcToPiKP == 1) && (pidBayesLcToPiKP == -1 || pidBayesLcToPiKP == 1)) {
        statusLcToPiKP = 1; // identified as LcToPiKP
      }

      if (candidate.cpa() <= cpaMin) {
        statusLcToPKPi = 0;
        statusLcToPiKP = 0;
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }

      std::array<float, 3> pVecPos1 = {trackPos1.px(), trackPos1.py(), trackPos1.pz()};
      std::array<float, 3> pVecNeg = {trackNeg.px(), trackNeg.py(), trackNeg.pz()};
      std::array<float, 3> pVecPos2 = {trackPos2.px(), trackPos2.py(), trackPos2.pz()};
      const float massPi = RecoDecay::getMassPDG(kPiPlus);
      const float massK = RecoDecay::getMassPDG(kKPlus);
      const float massProton = RecoDecay::getMassPDG(kProton);
      const float massLc = RecoDecay::getMassPDG(o2::analysis::pdg::kLambdaCPlus);
      if (statusLcToPiKP == 1) {
        auto invMassLcToPiKP = RecoDecay::m(std::array{pVecPos1, pVecNeg, pVecPos2}, std::array{massPi, massK, massProton});
        if (std::abs(invMassLcToPiKP - massLc) >= maxDeltaMass && candidate.pt() < 10) {
          statusLcToPiKP = 0;
        }
      }
      if (statusLcToPKPi == 1) {
        auto invMassLcToPKPi = RecoDecay::m(std::array{pVecPos1, pVecNeg, pVecPos2}, std::array{massProton, massK, massPi});
        if (std::abs(invMassLcToPKPi - massLc) >= maxDeltaMass && candidate.pt() < 10) {
          statusLcToPKPi = 0;
        }
      }

      if ((statusLcToPiKP == 1 || statusLcToPKPi == 1) && applyML) {
        auto trackParPos1 = getTrackPar(trackPos1);
        auto trackParNeg = getTrackPar(trackNeg);
        auto trackParPos2 = getTrackPar(trackPos2);
        std::vector<float> inputFeaturesF{trackParPos1.getPt(), trackPos1.dcaXY(), trackPos1.dcaZ(), trackParNeg.getPt(), trackNeg.dcaXY(), trackNeg.dcaZ(), trackParPos2.getPt(), trackPos2.dcaXY(), trackPos2.dcaZ()};
        std::vector<double> inputFeaturesD{trackParPos1.getPt(), trackPos1.dcaXY(), trackPos1.dcaZ(), trackParNeg.getPt(), trackNeg.dcaXY(), trackNeg.dcaZ(), trackParPos2.getPt(), trackPos2.dcaXY(), trackPos2.dcaZ()};
        float scores[3] = {-1.f, -1.f, -1.f};
        if (dataTypeML == 1) {
          auto scoresRaw = model.evalModel(inputFeaturesF);
          for (int iScore = 0; iScore < 3; ++iScore) {
            scores[iScore] = scoresRaw[iScore];
          }
        } else if (dataTypeML == 11) {
          auto scoresRaw = model.evalModel(inputFeaturesD);
          for (int iScore = 0; iScore < 3; ++iScore) {
            scores[iScore] = scoresRaw[iScore];
          }
        } else {
          LOG(error) << "Error running model inference for Lc: Unexpected input data type.";
        }
        if (scores[0] > thresholdBDTScoreLcToPiKP.value.get(0u, "BDTbkg")) {
          // background
          statusLcToPKPi = 0;
          statusLcToPiKP = 0;
        }
        // This is an equivalent to the cut above but it depends on the thresholds set
        if (scores[1] <= thresholdBDTScoreLcToPiKP.value.get(0u, "BDTprompt") &&
            scores[2] <= thresholdBDTScoreLcToPiKP.value.get(0u, "BDTnonprompt")) {
          statusLcToPKPi = 0;
          statusLcToPiKP = 0;
        }
        if (scores[1] > thresholdBDTScoreLcToPiKP.value.get(0u, "BDTprompt")) {
          // prompt
        }
        if (scores[2] > thresholdBDTScoreLcToPiKP.value.get(0u, "BDTnonprompt")) {
          // non-prompt
          // NOTE: Can be both prompt and non-prompt!
        }
        if (activateQA != 0) {
          registry.fill(HIST("hLcBDTScoreBkg"), scores[0]);
          registry.fill(HIST("hLcBDTScorePrompt"), scores[1]);
          registry.fill(HIST("hLcBDTScoreNonPrompt"), scores[2]);
        }
      }

      hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorLcMl>(cfgc)};
}
