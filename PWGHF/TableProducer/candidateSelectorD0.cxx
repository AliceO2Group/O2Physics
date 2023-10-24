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

/// \file candidateSelectorD0.cxx
/// \brief D0(bar) → π± K∓ selection task
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN
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

/// Struct for applying D0 selection cuts
struct HfCandidateSelectorD0 {
  Produces<aod::HfSelD0> hfSelD0Candidate;
  Produces<aod::HfMlD0> hfMlD0Candidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // AND logic for TOF+TPC PID (as in Run2)
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Use AND logic for TPC and TOF PID"};
  // selecting only background candidates
  Configurable<bool> keepOnlySidebandCandidates{"keepOnlySidebandCandidates", false, "Select only sideband candidates, for studying background cut variable distributions"};
  Configurable<double> distanceFromD0MassForSidebands{"distanceFromD0MassForSidebands", 0.15, "Minimum distance from nominal D0 mass value for sideband region"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_d0_to_pi_k::cuts[0], hf_cuts_d0_to_pi_k::nBinsPt, hf_cuts_d0_to_pi_k::nCutVars, hf_cuts_d0_to_pi_k::labelsPt, hf_cuts_d0_to_pi_k::labelsCutVar}, "D0 candidate selection per pT bin"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::cuts[0], hf_cuts_ml::nBinsPt, hf_cuts_ml::nCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int8_t> nClassesMl{"nClassesMl", (int8_t)hf_cuts_ml::nCutScores, "Number of classes in ML model"};
  Configurable<bool> enableDebugMl{"enableDebugMl", false, "Flag to enable histograms to monitor BDT application"};
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

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidKa>;

  // Define histograms
  AxisSpec axisMassDmeson{200, 1.7f, 2.1f};
  AxisSpec axisBdtScore{100, 0.f, 1.f};
  AxisSpec axisSelStatus{2, -0.5f, 1.5f};
  HistogramRegistry registry{"registry"};

  void init(InitContext& initContext)
  {
    std::array<bool, 2> doprocess{doprocessWithDCAFitterN, doprocessWithKFParticle};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function can be enabled at a time.");
    }

    if (applyMl && enableDebugMl) {
      registry.add("DebugBdt/hBdtScore1VsStatus", ";BDT score;status", {HistType::kTH2F, {axisBdtScore, axisSelStatus}});
      registry.add("DebugBdt/hBdtScore2VsStatus", ";BDT score;status", {HistType::kTH2F, {axisBdtScore, axisSelStatus}});
      registry.add("DebugBdt/hBdtScore3VsStatus", ";BDT score;status", {HistType::kTH2F, {axisBdtScore, axisSelStatus}});
      registry.add("DebugBdt/hMassDmesonSel", ";#it{M}(D) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {axisMassDmeson}});
    }

    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorKaon = selectorPion;

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

  /// Conjugate-independent topological cuts
  /// \param reconstructionType is the reconstruction type (DCAFitterN or KFParticle)
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <int reconstructionType, typename T>
  bool selectionTopol(const T& candidate)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }
    // product of daughter impact parameters
    if (candidate.impactParameterProduct() > cuts->get(pTBin, "d0d0")) {
      return false;
    }
    // cosine of pointing angle
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    // cosine of pointing angle XY
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle xy")) {
      return false;
    }
    // normalised decay length in XY plane
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    // candidate DCA
    // if (candidate.chi2PCA() > cuts[pTBin][1]) return false;

    // candidate topological chi2 over ndf when using KFParticle, need to add this selection to the SelectorCuts.h
    // if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
    //   if (candidate.kfTopolChi2OverNdf() > cuts->get(pTBin, "topological chi2overndf as D0")) return false;
    // }
    if (std::abs(candidate.impactParameterNormalised0()) < cuts->get(pTBin, "norm dauImpParX") || std::abs(candidate.impactParameterNormalised1()) < cuts->get(pTBin, "norm dauImpParX")) {
      return false;
    }
    if (candidate.decayLength() < cuts->get(pTBin, "minimum decay length")) {
      return false;
    }
    if (candidate.decayLength() > cuts->get(pTBin, "max decay length")) {
      return false;
    }
    if (candidate.decayLengthXY() > cuts->get(pTBin, "decay length XY")) {
      return false;
    }

    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param reconstructionType is the reconstruction type (DCAFitterN or KFParticle)
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \note trackPion = positive and trackKaon = negative for D0 selection and inverse for D0bar
  /// \return true if candidate passes all cuts for the given Conjugate
  template <int reconstructionType, typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& trackKaon)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // invariant-mass cut
    float massD0, massD0bar;
    if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
      massD0 = candidate.kfGeoMassD0();
      massD0bar = candidate.kfGeoMassD0bar();
    } else {
      massD0 = hfHelper.invMassD0ToPiK(candidate);
      massD0bar = hfHelper.invMassD0barToKPi(candidate);
    }
    if (trackPion.sign() > 0) {
      if (std::abs(massD0 - o2::analysis::pdg::MassD0) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(massD0bar - o2::analysis::pdg::MassD0) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    // cut on daughter pT
    if (trackPion.pt() < cuts->get(pTBin, "pT Pi") || trackKaon.pt() < cuts->get(pTBin, "pT K")) {
      return false;
    }

    // cut on daughter DCA - need to add secondary vertex constraint here
    if (std::abs(trackPion.dcaXY()) > cuts->get(pTBin, "d0pi") || std::abs(trackKaon.dcaXY()) > cuts->get(pTBin, "d0K")) {
      return false;
    }

    // cut on cos(theta*)
    if (trackPion.sign() > 0) {
      if (std::abs(hfHelper.cosThetaStarD0(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    } else {
      if (std::abs(hfHelper.cosThetaStarD0bar(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    }

    // in case only sideband candidates have to be stored, additional invariant-mass cut
    if (keepOnlySidebandCandidates) {
      if (trackPion.sign() > 0) {
        if (std::abs(hfHelper.invMassD0ToPiK(candidate) - o2::analysis::pdg::MassD0) < distanceFromD0MassForSidebands) {
          return false;
        }
      } else {
        if (std::abs(hfHelper.invMassD0barToKPi(candidate) - o2::analysis::pdg::MassD0) < distanceFromD0MassForSidebands) {
          return false;
        }
      }
    }

    return true;
  }
  template <int reconstructionType, typename CandType>
  void processSel(CandType const& candidates,
                  TracksSel const&)
  {
    // looping over 2-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      int statusD0 = 0;
      int statusD0bar = 0;
      int statusHFFlag = 0;
      int statusTopol = 0;
      int statusCand = 0;
      int statusPID = 0;

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMl);
        }
        continue;
      }
      statusHFFlag = 1;

      auto ptCand = candidate.pt();
      auto trackPos = candidate.template prong0_as<TracksSel>(); // positive daughter
      auto trackNeg = candidate.template prong1_as<TracksSel>(); // negative daughter

      // conjugate-independent topological selection
      if (!selectionTopol<reconstructionType>(candidate)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMl);
        }
        continue;
      }
      statusTopol = 1;

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      // conjugate-dependent topological selection for D0
      bool topolD0 = selectionTopolConjugate<reconstructionType>(candidate, trackPos, trackNeg);
      // conjugate-dependent topological selection for D0bar
      bool topolD0bar = selectionTopolConjugate<reconstructionType>(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMl);
        }
        continue;
      }
      statusCand = 1;

      // track-level PID selection
      int pidTrackPosKaon = -1;
      int pidTrackPosPion = -1;
      int pidTrackNegKaon = -1;
      int pidTrackNegPion = -1;

      if (usePidTpcAndTof) {
        pidTrackPosKaon = selectorKaon.statusTpcAndTof(trackPos);
        pidTrackPosPion = selectorPion.statusTpcAndTof(trackPos);
        pidTrackNegKaon = selectorKaon.statusTpcAndTof(trackNeg);
        pidTrackNegPion = selectorPion.statusTpcAndTof(trackNeg);
      } else {
        pidTrackPosKaon = selectorKaon.statusTpcOrTof(trackPos);
        pidTrackPosPion = selectorPion.statusTpcOrTof(trackPos);
        pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg);
        pidTrackNegPion = selectorPion.statusTpcOrTof(trackNeg);
      }

      // int pidBayesTrackPos1Pion = selectorPion.statusBayes(trackPos);

      int pidD0 = -1;
      int pidD0bar = -1;

      if (pidTrackPosPion == TrackSelectorPID::Accepted &&
          pidTrackNegKaon == TrackSelectorPID::Accepted) {
        pidD0 = 1; // accept D0
      } else if (pidTrackPosPion == TrackSelectorPID::Rejected ||
                 pidTrackNegKaon == TrackSelectorPID::Rejected) {
        pidD0 = 0; // exclude D0
      }

      if (pidTrackNegPion == TrackSelectorPID::Accepted &&
          pidTrackPosKaon == TrackSelectorPID::Accepted) {
        pidD0bar = 1; // accept D0bar
      } else if (pidTrackNegPion == TrackSelectorPID::Rejected ||
                 pidTrackPosKaon == TrackSelectorPID::Rejected) {
        pidD0bar = 0; // exclude D0bar
      }

      if (pidD0 == 0 && pidD0bar == 0) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMl);
        }
        continue;
      }

      if ((pidD0 == -1 || pidD0 == 1) && topolD0) {
        statusD0 = 1; // identified as D0
      }
      if ((pidD0bar == -1 || pidD0bar == 1) && topolD0bar) {
        statusD0bar = 1; // identified as D0bar
      }
      statusPID = 1;

      if (applyMl) {
        // ML selections

        std::vector<float> inputFeatures{candidate.cpa(),
                                         candidate.cpaXY(),
                                         candidate.decayLength(),
                                         candidate.decayLengthXY(),
                                         candidate.impactParameter0(),
                                         candidate.impactParameter1(),
                                         candidate.impactParameterProduct()};

        bool isSelectedMl = hfMlResponse.isSelectedMl(inputFeatures, ptCand, outputMl);
        hfMlD0Candidate(outputMl);

        if (!isSelectedMl) {
          statusD0 = 0;
          statusD0bar = 0;
        }
        if (enableDebugMl) {
          registry.fill(HIST("DebugBdt/hBdtScore1VsStatus"), outputMl[0], statusD0);
          registry.fill(HIST("DebugBdt/hBdtScore1VsStatus"), outputMl[0], statusD0bar);
          registry.fill(HIST("DebugBdt/hBdtScore2VsStatus"), outputMl[1], statusD0);
          registry.fill(HIST("DebugBdt/hBdtScore2VsStatus"), outputMl[1], statusD0bar);
          registry.fill(HIST("DebugBdt/hBdtScore3VsStatus"), outputMl[2], statusD0);
          registry.fill(HIST("DebugBdt/hBdtScore3VsStatus"), outputMl[2], statusD0bar);
          if (statusD0 > 0) {
            registry.fill(HIST("DebugBdt/hMassDmesonSel"), hfHelper.invMassD0ToPiK(candidate));
          }
          if (statusD0bar > 0) {
            registry.fill(HIST("DebugBdt/hMassDmesonSel"), hfHelper.invMassD0barToKPi(candidate));
          }
        }
      }
      hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
    }
  }

  void processWithDCAFitterN(aod::HfCand2Prong const& candidates, TracksSel const& tracks)
  {
    processSel<aod::hf_cand::VertexerType::DCAFitter>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorD0, processWithDCAFitterN, "process candidates selection with DCAFitterN", true);

  void processWithKFParticle(soa::Join<aod::HfCand2Prong, aod::HfCand2ProngKF> const& candidates, TracksSel const& tracks)
  {
    processSel<aod::hf_cand::VertexerType::KfParticle>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorD0, processWithKFParticle, "process candidates selection with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorD0>(cfgc)};
}
