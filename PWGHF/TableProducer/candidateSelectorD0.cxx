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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseD0ToKPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/TrackSelectorPID.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
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

#include <array>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Struct for applying D0 selection cuts
struct HfCandidateSelectorD0 {
  Produces<aod::HfSelD0> hfSelD0Candidate;
  Produces<aod::HfMlD0> hfMlD0Candidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  Configurable<bool> usePid{"usePid", true, "Use PID selection"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  Configurable<bool> usePidTpcOnly{"usePidTpcOnly", false, "Only use TPC PID"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // AND logic for TOF+TPC PID (as in Run2)
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Use AND logic for TPC and TOF PID"};
  // ITS quality track cuts
  Configurable<int> itsNClustersFoundMin{"itsNClustersFoundMin", 0, "Minimum number of found ITS clusters"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 1e10f, "Maximum its fit chi2 per ITS cluster"};
  // TPC quality track cuts
  Configurable<int> tpcNClustersFoundMin{"tpcNClustersFoundMin", 0, "Minimum number of found TPC clusters"};
  Configurable<int> tpcNCrossedRowsMin{"tpcNCrossedRowsMin", 0, "Minimum number of crossed rows in TPC"};
  Configurable<float> tpcNCrossedRowsOverFindableClustersMin{"tpcNCrossedRowsOverFindableClustersMin", 0., "Minimum ratio crossed rows / findable clusters"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 1e10f, "Maximum TPC fit chi2 per TPC cluster"};
  // selecting only background candidates
  Configurable<bool> keepOnlySidebandCandidates{"keepOnlySidebandCandidates", false, "Select only sideband candidates, for studying background cut variable distributions"};
  Configurable<double> distanceFromD0MassForSidebands{"distanceFromD0MassForSidebands", 0.15, "Minimum distance from nominal D0 mass value for sideband region"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_d0_to_pi_k::Cuts[0], hf_cuts_d0_to_pi_k::NBinsPt, hf_cuts_d0_to_pi_k::NCutVars, hf_cuts_d0_to_pi_k::labelsPt, hf_cuts_d0_to_pi_k::labelsCutVar}, "D0 candidate selection per pT bin"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<bool> enableDebugMl{"enableDebugMl", false, "Flag to enable histograms to monitor BDT application"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTD0"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_D0ToKPi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  // Mass Cut for trigger analysis
  Configurable<bool> useTriggerMassCut{"useTriggerMassCut", false, "Flag to enable parametrize pT differential mass cut for triggered data"};

  o2::analysis::HfMlResponseD0ToKPi<float> hfMlResponse;
  std::vector<float> outputMlD0;
  std::vector<float> outputMlD0bar;
  o2::ccdb::CcdbApi ccdbApi;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  HfTrigger2ProngCuts hfTriggerCuts;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;

  // Define histograms
  AxisSpec axisMassDmeson{200, 1.7f, 2.1f};
  AxisSpec axisBdtScore{100, 0.f, 1.f};
  AxisSpec axisSelStatus{2, -0.5f, 1.5f};
  HistogramRegistry registry{"registry"};

  void init(InitContext&)
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
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
    }
  }

  /// Single track quality cuts
  /// \param track is track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedCandidateProng(const T& trackPos, const T& trackNeg)
  {
    if (!isSelectedTrackTpcQuality(trackPos, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value) ||
        !isSelectedTrackTpcQuality(trackNeg, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value)) {
      return false;
    }
    if (!isSelectedTrackItsQuality(trackPos, itsNClustersFoundMin.value, itsChi2PerClusterMax.value) ||
        !isSelectedTrackItsQuality(trackNeg, itsNClustersFoundMin.value, itsChi2PerClusterMax.value)) {
      return false;
    }
    return true;
  }

  /// Conjugate-independent topological cuts
  /// \param reconstructionType is the reconstruction type (DCAFitterN or KFParticle)
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <int ReconstructionType, typename T>
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
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "min norm decay length XY")) {
      return false;
    }
    // candidate DCA
    if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "DCA")) {
      return false;
    }

    // candidate topological chi2 over ndf when using KFParticle, need to add this selection to the SelectorCuts.h
    // if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
    //   if (candidate.kfTopolChi2OverNdf() > cuts->get(pTBin, "topological chi2overndf as D0")) return false;
    // }
    if (std::abs(candidate.impactParameterNormalised0()) < cuts->get(pTBin, "norm dauImpPar XY") || std::abs(candidate.impactParameterNormalised1()) < cuts->get(pTBin, "norm dauImpPar XY")) {
      return false;
    }
    if (candidate.decayLength() < cuts->get(pTBin, "min decay length")) {
      return false;
    }
    if (candidate.decayLength() > cuts->get(pTBin, "max decay length")) {
      return false;
    }
    if (candidate.decayLengthXY() > cuts->get(pTBin, "max decay length XY")) {
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
  template <int ReconstructionType, typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& trackKaon)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // invariant-mass cut
    float massD0, massD0bar;
    if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
      massD0 = candidate.kfGeoMassD0();
      massD0bar = candidate.kfGeoMassD0bar();
    } else {
      massD0 = HfHelper::invMassD0ToPiK(candidate);
      massD0bar = HfHelper::invMassD0barToKPi(candidate);
    }
    if (trackPion.sign() > 0) {
      if (std::abs(massD0 - o2::constants::physics::MassD0) > cuts->get(pTBin, "m")) {
        return false;
      }
      if (useTriggerMassCut && !isCandidateInMassRange(massD0, o2::constants::physics::MassD0, candidate.pt(), hfTriggerCuts)) {
        return false;
      }
    } else {
      if (std::abs(massD0bar - o2::constants::physics::MassD0) > cuts->get(pTBin, "m")) {
        return false;
      }
      if (useTriggerMassCut && !isCandidateInMassRange(massD0bar, o2::constants::physics::MassD0, candidate.pt(), hfTriggerCuts)) {
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
      if (std::abs(HfHelper::cosThetaStarD0(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    } else {
      if (std::abs(HfHelper::cosThetaStarD0bar(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    }

    // in case only sideband candidates have to be stored, additional invariant-mass cut
    if (keepOnlySidebandCandidates) {
      if (trackPion.sign() > 0) {
        if (std::abs(HfHelper::invMassD0ToPiK(candidate) - o2::constants::physics::MassD0) < distanceFromD0MassForSidebands) {
          return false;
        }
      } else {
        if (std::abs(HfHelper::invMassD0barToKPi(candidate) - o2::constants::physics::MassD0) < distanceFromD0MassForSidebands) {
          return false;
        }
      }
    }

    return true;
  }
  template <int ReconstructionType, typename CandType>
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

      outputMlD0.clear();
      outputMlD0bar.clear();

      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMlD0, outputMlD0bar);
        }
        continue;
      }
      statusHFFlag = 1;

      auto ptCand = candidate.pt();
      auto trackPos = candidate.template prong0_as<TracksSel>(); // positive daughter
      auto trackNeg = candidate.template prong1_as<TracksSel>(); // negative daughter

      // implement track quality selection for D0 daughters
      if (!isSelectedCandidateProng(trackPos, trackNeg)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMlD0, outputMlD0bar);
        }
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol<ReconstructionType>(candidate)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMlD0, outputMlD0bar);
        }
        continue;
      }
      statusTopol = 1;

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      // conjugate-dependent topological selection for D0
      bool const topolD0 = selectionTopolConjugate<ReconstructionType>(candidate, trackPos, trackNeg);
      // conjugate-dependent topological selection for D0bar
      bool const topolD0bar = selectionTopolConjugate<ReconstructionType>(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlD0Candidate(outputMlD0, outputMlD0bar);
        }
        continue;
      }
      statusCand = 1;

      if (usePid) {
        // track-level PID selection
        int pidTrackPosKaon = -1;
        int pidTrackPosPion = -1;
        int pidTrackNegKaon = -1;
        int pidTrackNegPion = -1;

        if (usePidTpcOnly) {
          /// kaon TPC PID positive daughter
          pidTrackPosKaon = selectorKaon.statusTpc(trackPos, candidate.nSigTpcKa0());
          /// pion TPC PID positive daughter
          pidTrackPosPion = selectorPion.statusTpc(trackPos, candidate.nSigTpcPi0());
          /// kaon TPC PID negative daughter
          pidTrackNegKaon = selectorKaon.statusTpc(trackNeg, candidate.nSigTpcKa1());
          /// pion TPC PID negative daughter
          pidTrackNegPion = selectorPion.statusTpc(trackNeg, candidate.nSigTpcPi1());
        } else if (usePidTpcAndTof) {
          /// kaon TPC, TOF PID positive daughter
          pidTrackPosKaon = selectorKaon.statusTpcAndTof(trackPos, candidate.nSigTpcKa0(), candidate.nSigTofKa0());
          /// pion TPC, TOF PID positive daughter
          pidTrackPosPion = selectorPion.statusTpcAndTof(trackPos, candidate.nSigTpcPi0(), candidate.nSigTofPi0());
          /// kaon TPC, TOF PID negative daughter
          pidTrackNegKaon = selectorKaon.statusTpcAndTof(trackNeg, candidate.nSigTpcKa1(), candidate.nSigTofKa1());
          /// pion TPC, TOF PID negative daughter
          pidTrackNegPion = selectorPion.statusTpcAndTof(trackNeg, candidate.nSigTpcPi1(), candidate.nSigTofPi1());
        } else {
          /// kaon TPC, TOF PID positive daughter
          pidTrackPosKaon = selectorKaon.statusTpcOrTof(trackPos, candidate.nSigTpcKa0(), candidate.nSigTofKa0());
          /// pion TPC, TOF PID positive daughter
          pidTrackPosPion = selectorPion.statusTpcOrTof(trackPos, candidate.nSigTpcPi0(), candidate.nSigTofPi0());
          /// kaon TPC, TOF PID negative daughter
          pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg, candidate.nSigTpcKa1(), candidate.nSigTofKa1());
          /// pion TPC, TOF PID negative daughter
          pidTrackNegPion = selectorPion.statusTpcOrTof(trackNeg, candidate.nSigTpcPi1(), candidate.nSigTofPi1());
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
            hfMlD0Candidate(outputMlD0, outputMlD0bar);
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
      } else {
        if (topolD0) {
          statusD0 = 1; // identified as D0
        }
        if (topolD0bar) {
          statusD0bar = 1; // identified as D0bar
        }
      }

      if (applyMl) {
        // ML selections
        bool isSelectedMlD0 = false;
        bool isSelectedMlD0bar = false;

        if (statusD0 > 0) {
          std::vector<float> inputFeaturesD0 = hfMlResponse.getInputFeatures(candidate, o2::constants::physics::kD0);
          isSelectedMlD0 = hfMlResponse.isSelectedMl(inputFeaturesD0, ptCand, outputMlD0);
        }
        if (statusD0bar > 0) {
          std::vector<float> inputFeaturesD0bar = hfMlResponse.getInputFeatures(candidate, o2::constants::physics::kD0Bar);
          isSelectedMlD0bar = hfMlResponse.isSelectedMl(inputFeaturesD0bar, ptCand, outputMlD0bar);
        }

        if (!isSelectedMlD0) {
          statusD0 = 0;
        }
        if (!isSelectedMlD0bar) {
          statusD0bar = 0;
        }

        hfMlD0Candidate(outputMlD0, outputMlD0bar);

        if (enableDebugMl) {
          if (isSelectedMlD0) {
            registry.fill(HIST("DebugBdt/hBdtScore1VsStatus"), outputMlD0[0], statusD0);
            registry.fill(HIST("DebugBdt/hBdtScore2VsStatus"), outputMlD0[1], statusD0);
            registry.fill(HIST("DebugBdt/hBdtScore3VsStatus"), outputMlD0[2], statusD0);
            registry.fill(HIST("DebugBdt/hMassDmesonSel"), HfHelper::invMassD0ToPiK(candidate));
          }
          if (isSelectedMlD0bar) {
            registry.fill(HIST("DebugBdt/hBdtScore1VsStatus"), outputMlD0bar[0], statusD0bar);
            registry.fill(HIST("DebugBdt/hBdtScore2VsStatus"), outputMlD0bar[1], statusD0bar);
            registry.fill(HIST("DebugBdt/hBdtScore3VsStatus"), outputMlD0bar[2], statusD0bar);
            registry.fill(HIST("DebugBdt/hMassDmesonSel"), HfHelper::invMassD0barToKPi(candidate));
          }
        }
      }
      hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
    }
  }

  void processWithDCAFitterN(aod::HfCand2ProngWPid const& candidates, TracksSel const& tracks)
  {
    processSel<aod::hf_cand::VertexerType::DCAFitter>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorD0, processWithDCAFitterN, "process candidates selection with DCAFitterN", true);

  void processWithKFParticle(soa::Join<aod::HfCand2ProngWPid, aod::HfCand2ProngKF> const& candidates, TracksSel const& tracks)
  {
    processSel<aod::hf_cand::VertexerType::KfParticle>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorD0, processWithKFParticle, "process candidates selection with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorD0>(cfgc)};
}
