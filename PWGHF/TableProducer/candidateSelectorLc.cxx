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

#include <string>
#include <vector>

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
  Configurable<double> ptPidBayesMin{"ptPidBayesMin", 0., "Lower bound of track pT for Bayesian PID"};
  Configurable<double> ptPidBayesMax{"ptPidBayesMax", 100, "Upper bound of track pT for Bayesian PID"};
  // Combined PID options
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Bool to decide how to combine TPC and TOF PID: true = both (if present, only one otherwise); false = one is enough"};
  // TPC quality track cuts
  Configurable<int> tpcNClustersFoundMin{"tpcNClustersFoundMin", 0, "min number of found TPC clusters"};
  Configurable<int> tpcNCrossedRowsMin{"tpcNCrossedRowsMin", 0, "min number of crossed rows in TPC"};
  Configurable<float> tpcNCrossedRowsOverFindableClustersMin{"tpcNCrossedRowsOverFindableClustersMin", 0., "min ratio crossed rows / findable clusters"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 1e10f, "max tpc fit chi2 per TPC cluster"};
  // ITS quality track cuts
  Configurable<int> itsNClustersFoundMin{"itsNClustersFoundMin", 0, "min. number of found ITS clusters"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 1e10f, "max its fit chi2 per ITS cluster"};
  // DCA track cuts
  Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for DCA XY/Z pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsSingleTrack{"cutsSingleTrack", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_p_k_pi::Cuts[0], hf_cuts_lc_to_p_k_pi::NBinsPt, hf_cuts_lc_to_p_k_pi::NCutVars, hf_cuts_lc_to_p_k_pi::labelsPt, hf_cuts_lc_to_p_k_pi::labelsCutVar}, "Lc candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
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
                              aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using TracksSelBayesPid = soa::Join<TracksSel, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidBayes>;

  HistogramRegistry registry{"registry"};

  double massK0Star892;

  void init(InitContext const&)
  {
    std::array<bool, 4> processes = {doprocessNoBayesPidWithDCAFitterN, doprocessBayesPidWithDCAFitterN, doprocessNoBayesPidWithKFParticle, doprocessBayesPidWithKFParticle};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }

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

    massK0Star892 = o2::constants::physics::MassK0Star892;
  }

  /// Single track quality cuts
  /// \param track is track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedCandidateProngQuality(const T& trackPos1, const T& trackNeg, const T& trackPos2)
  {
    if (!isSelectedTrackTpcQuality(trackPos1, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value) ||
        !isSelectedTrackTpcQuality(trackNeg, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value) ||
        !isSelectedTrackTpcQuality(trackPos2, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value)) {
      return false;
    }
    if (!isSelectedTrackItsQuality(trackPos1, itsNClustersFoundMin.value, itsChi2PerClusterMax.value) ||
        !isSelectedTrackItsQuality(trackNeg, itsNClustersFoundMin.value, itsChi2PerClusterMax.value) ||
        !isSelectedTrackItsQuality(trackPos2, itsNClustersFoundMin.value, itsChi2PerClusterMax.value)) {
      return false;
    }
    return true;
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

    // candidate decay length XY
    if (candidate.decayLengthXY() <= cuts->get(pTBin, "decLengthXY")) {
      return false;
    }

    // candidate normalized decay length XY
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normDecLXY")) {
      return false;
    }

    // candidate impact parameter XY
    if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "impParXY")) {
      return false;
    }

    if (!isSelectedCandidateProngDca(candidate)) {
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
  template <aod::hf_cand::VertexerType reconstructionType, typename T1, typename T2>
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

    // cut on Lc->pKpi, piKp mass values
    /// cut on the Kpi pair invariant mass, to study Lc->pK*(->Kpi)
    float massLc, massKPi;
    if constexpr (reconstructionType == aod::hf_cand::VertexerType::DCAFitter) {
      if (trackProton.globalIndex() == candidate.prong0Id()) {
        massLc = hfHelper.invMassLcToPKPi(candidate);
        massKPi = hfHelper.invMassKPiPairLcToPKPi(candidate);
      } else {
        massLc = hfHelper.invMassLcToPiKP(candidate);
        massKPi = hfHelper.invMassKPiPairLcToPiKP(candidate);
      }
    } else if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
      if (trackProton.globalIndex() == candidate.prong0Id()) {
        massLc = candidate.kfMassPKPi();
        massKPi = candidate.kfMassKPi();
      } else {
        massLc = candidate.kfMassPiKP();
        massKPi = candidate.kfMassPiK();
      }
    }

    if (std::abs(massLc - o2::constants::physics::MassLambdaCPlus) > cuts->get(pTBin, "m")) {
      return false;
    }

    /// cut on the Kpi pair invariant mass, to study Lc->pK*(->Kpi)
    const double cutMassKPi = cuts->get(pTBin, "mass (Kpi)");
    if (cutMassKPi > 0 && std::abs(massKPi - massK0Star892) > cutMassKPi) {
      return false;
    }

    return true;
  }

  /// Single-track dca_xy and dca_z cuts
  /// \param candidate is the Lc candidate
  /// \return true if all the prongs pass the selections
  template <typename T1>
  bool isSelectedCandidateProngDca(const T1& candidate)
  {
    return (isSelectedTrackDca(binsPtTrack, cutsSingleTrack, candidate.ptProng0(), candidate.impactParameter0(), candidate.impactParameterZ0()) &&
            isSelectedTrackDca(binsPtTrack, cutsSingleTrack, candidate.ptProng1(), candidate.impactParameter1(), candidate.impactParameterZ1()) &&
            isSelectedTrackDca(binsPtTrack, cutsSingleTrack, candidate.ptProng2(), candidate.impactParameter2(), candidate.impactParameterZ2()));
  }

  /// Apply PID selection
  /// \param pidTrackProton is the PID status of proton candidate track
  /// \param pidTrackKaon is the PID status of kaon candidate track
  /// \param pidTrackPion is the PID status of pion candidate track
  /// \return true if prongs pass all selections
  bool isSelectedPID(const TrackSelectorPID::Status pidTrackProton, const TrackSelectorPID::Status pidTrackKaon, const TrackSelectorPID::Status pidTrackPion)
  {
    if (pidTrackProton == TrackSelectorPID::Rejected ||
        pidTrackKaon == TrackSelectorPID::Rejected ||
        pidTrackPion == TrackSelectorPID::Rejected) {
      return false;
    }
    return true;
  }

  /// \brief function to apply Lc selections
  /// \param reconstructionType is the reconstruction type (DCAFitterN or KFParticle)
  /// \param candidates Lc candidate table
  /// \param tracks track table
  template <bool useBayesPid = false, aod::hf_cand::VertexerType reconstructionType, typename CandType, typename TTracks>
  void runSelectLc(CandType const& candidates, TTracks const&)
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

      auto trackPos1 = candidate.template prong0_as<TTracks>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.template prong1_as<TTracks>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.template prong2_as<TTracks>(); // positive daughter (negative for the antiparticles)

      // implement filter bit 4 cut - should be done before this task at the track selection level

      // track quality selection
      bool trackQualitySel = isSelectedCandidateProngQuality(trackPos1, trackNeg, trackPos2);
      if (!trackQualitySel) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        if (applyMl) {
          hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);
        }
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        if (applyMl) {
          hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);
        }
        continue;
      }

      // conjugate-dependent topological selection for Lc
      bool topolLcToPKPi = selectionTopolConjugate<reconstructionType>(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcToPiKP = selectionTopolConjugate<reconstructionType>(candidate, trackPos2, trackNeg, trackPos1);

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

      // PID not applied, accepted by default
      auto pidLcToPKPi = 1;
      auto pidLcToPiKP = 1;
      auto pidBayesLcToPKPi = 1;
      auto pidBayesLcToPiKP = 1;

      if (usePid) {
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

        if (!isSelectedPID(pidTrackPos1Proton, pidTrackNegKaon, pidTrackPos2Pion)) {
          pidLcToPKPi = 0; // reject LcToPKPi
        }
        if (!isSelectedPID(pidTrackPos2Proton, pidTrackNegKaon, pidTrackPos1Pion)) {
          pidLcToPiKP = 0; // accept LcToPiKP
        }
      }

      if constexpr (useBayesPid) {
        TrackSelectorPID::Status pidBayesTrackPos1Proton = selectorProton.statusBayes(trackPos1);
        TrackSelectorPID::Status pidBayesTrackPos2Proton = selectorProton.statusBayes(trackPos2);
        TrackSelectorPID::Status pidBayesTrackPos1Pion = selectorPion.statusBayes(trackPos1);
        TrackSelectorPID::Status pidBayesTrackPos2Pion = selectorPion.statusBayes(trackPos2);
        TrackSelectorPID::Status pidBayesTrackNegKaon = selectorKaon.statusBayes(trackNeg);

        if (!isSelectedPID(pidBayesTrackPos1Proton, pidBayesTrackNegKaon, pidBayesTrackPos2Pion)) {
          pidBayesLcToPKPi = 0; // reject LcToPKPi
        }

        if (!isSelectedPID(pidBayesTrackPos2Proton, pidBayesTrackNegKaon, pidBayesTrackPos1Pion)) {
          pidBayesLcToPiKP = 0; // reject LcToPiKP
        }
      }

      if ((pidLcToPKPi == 0 && pidLcToPiKP == 0) || (pidBayesLcToPKPi == 0 && pidBayesLcToPiKP == 0)) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        if (applyMl) {
          hfMlLcToPKPiCandidate(outputMlLcToPKPi, outputMlLcToPiKP);
        }
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, candidate.pt());
      }

      bool isSelectedMlLcToPKPi = true;
      bool isSelectedMlLcToPiKP = true;
      if (applyMl) {
        // ML selections
        isSelectedMlLcToPKPi = false;
        isSelectedMlLcToPiKP = false;

        if (pidLcToPKPi == 1 && pidBayesLcToPKPi == 1 && topolLcToPKPi) {
          std::vector<float> inputFeaturesLcToPKPi = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2, true);
          isSelectedMlLcToPKPi = hfMlResponse.isSelectedMl(inputFeaturesLcToPKPi, candidate.pt(), outputMlLcToPKPi);
        }
        if (pidLcToPiKP == 1 && pidBayesLcToPiKP == 1 && topolLcToPiKP) {
          std::vector<float> inputFeaturesLcToPiKP = hfMlResponse.getInputFeatures(candidate, trackPos1, trackNeg, trackPos2, false);
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

      if (pidLcToPKPi == 1 && pidBayesLcToPKPi == 1 && isSelectedMlLcToPKPi && topolLcToPKPi && trackQualitySel) {
        statusLcToPKPi = 1; // identified as LcToPKPi
      }
      if (pidLcToPiKP == 1 && pidBayesLcToPiKP == 1 && isSelectedMlLcToPiKP && topolLcToPiKP && trackQualitySel) {
        statusLcToPiKP = 1; // identified as LcToPiKP
      }

      hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
    }
  }

  /// \brief process function w/o Bayes PID with DCAFitterN
  /// \param candidates Lc candidate table
  /// \param tracks track table
  void processNoBayesPidWithDCAFitterN(aod::HfCand3Prong const& candidates,
                                       TracksSel const& tracks)
  {
    runSelectLc<false, aod::hf_cand::VertexerType::DCAFitter>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorLc, processNoBayesPidWithDCAFitterN, "Process Lc selection w/o Bayes PID with DCAFitterN", true);

  /// \brief process function with Bayes PID with DCAFitterN
  /// \param candidates Lc candidate table
  /// \param tracks track table with Bayes PID information
  void processBayesPidWithDCAFitterN(aod::HfCand3Prong const& candidates,
                                     TracksSelBayesPid const& tracks)
  {
    runSelectLc<true, aod::hf_cand::VertexerType::DCAFitter>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorLc, processBayesPidWithDCAFitterN, "Process Lc selection with Bayes PID with DCAFitterN", false);

  /// \brief process function w/o Bayes PID with KFParticle
  /// \param candidates Lc candidate table
  /// \param tracks track table
  void processNoBayesPidWithKFParticle(soa::Join<aod::HfCand3Prong, aod::HfCand3ProngKF> const& candidates,
                                       TracksSel const& tracks)
  {
    runSelectLc<false, aod::hf_cand::VertexerType::KfParticle>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorLc, processNoBayesPidWithKFParticle, "Process Lc selection w/o Bayes PID with KFParticle", false);

  /// \brief process function with Bayes PID with KFParticle
  /// \param candidates Lc candidate table
  /// \param tracks track table with Bayes PID information
  void processBayesPidWithKFParticle(soa::Join<aod::HfCand3Prong, aod::HfCand3ProngKF> const& candidates,
                                     TracksSelBayesPid const& tracks)
  {
    runSelectLc<true, aod::hf_cand::VertexerType::KfParticle>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorLc, processBayesPidWithKFParticle, "Process Lc selection with Bayes PID with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorLc>(cfgc)};
}
