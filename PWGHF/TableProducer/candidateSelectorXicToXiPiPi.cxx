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

#include "PWGHF/Core/HfMlResponseXicToXiPiPi.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

#include "Common/Core/TrackSelectorPID.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::analysis;

struct HfCandidateSelectorXicToXiPiPi {
  Produces<aod::HfSelXicToXiPiPi> hfSelXicToXiPiPiCandidate;
  Produces<aod::HfMlXicToXiPiPi> hfMlXicToXiPiPiCandidate;

  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_xi_pi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_xi_pi_pi::Cuts[0], hf_cuts_xic_to_xi_pi_pi::NBinsPt, hf_cuts_xic_to_xi_pi_pi::NCutVars, hf_cuts_xic_to_xi_pi_pi::labelsPt, hf_cuts_xic_to_xi_pi_pi::labelsCutVar}, "Xicplus candidate selection per pT bin"};
  Configurable<bool> fillHistogram{"fillHistogram", false, "Flag to filling of counter histogram"};
  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Switch for PID selection at track level"};
  Configurable<bool> useTpcPidOnly{"useTpcPidOnly", false, "Switch to use TPC PID only instead of TPC OR TOF)"};
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
  // TrackQualitySelection
  Configurable<bool> doTrackQualitySelection{"doTrackQualitySelection", true, "Switch to apply track quality selections on final state daughter particles"};
  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement"};
  Configurable<int> nTpcCrossedRowsMin{"nTpcCrossedRowsMin", 70, "Minimum number of TPC crossed rows requirement"};
  Configurable<double> tpcCrossedRowsOverFindableClustersRatioMin{"tpcCrossedRowsOverFindableClustersRatioMin", 0.8, "Minimum ratio TPC crossed rows over findable clusters requirement"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 4, "Maximum value of chi2 fit over TPC clusters"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pi <- charm baryon"};
  Configurable<int> nClustersItsInnBarrMin{"nClustersItsInnBarrMin", 1, "Minimum number of ITS clusters in inner barrel requirement for pi <- charm baryon"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 36, "Maximum value of chi2 fit over ITS clusters for pi <- charm baryon"};
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
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
  enum XicSelCriteria { All = 0,
                        Pt,
                        Mass,
                        Rapidity,
                        Eta,
                        EtaDaughters,
                        PtPionFromXicPlus,
                        Chi2SV,
                        MinDecayLength,
                        MaxInvMassXiPiPairs,
                        TpcTrackQualityXiDaughters,
                        TpcTrackQualityPiFromCharm,
                        ItsTrackQualityPiFromCharm,
                        PidSelected,
                        BdtSelected };

  using TracksExtraWPid = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::TracksPidPr>;

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

    if (fillHistogram) {
      registry.add("hSelCandidates", ";;entries", {HistType::kTH1F, {{15, -0.5, 14.5}}});
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + All, "All");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + Pt, "#it{p}_{T}");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + Mass, "#Delta M");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + Rapidity, "y");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + Eta, "#eta");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + EtaDaughters, "#eta final state daughters");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + PtPionFromXicPlus, "#it{p}_{T} (#pi #leftarrow #Xi_{c}^{+})");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + Chi2SV, "#chi^{2}_{SV}");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + MinDecayLength, "Decay length");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + MaxInvMassXiPiPairs, "M_{#Xi #pi}");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + TpcTrackQualityXiDaughters, "TPC track quality selection on #Xi daughters");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + TpcTrackQualityPiFromCharm, "TPC track quality selection on #pi #leftarrow #Xi_{c}^{+}");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + ItsTrackQualityPiFromCharm, "ITS track quality selection on #pi #leftarrow #Xi_{c}^{+}");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + PidSelected, "PID selection");
      registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(1 + BdtSelected, "BDT selection");
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
  bool isSelectedXic(const T1& hfCandXic,
                     const float& etaPi0,
                     const float& etaPi1,
                     const float& etaPiFromXi,
                     const float& etaV0PosDau,
                     const float& etaV0NegDau)
  {
    auto candpT = hfCandXic.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), Pt);
    }

    // check whether candidate mass is within a defined mass window
    if (std::abs(hfCandXic.invMassXicPlus() - o2::constants::physics::MassXiCPlus) > cuts->get(pTBin, "m")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), Mass);
    }

    // cut on candidate rapidity
    if (std::abs(hfCandXic.y(o2::constants::physics::MassXiCPlus)) > cuts->get(pTBin, "y")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), Rapidity);
    }

    // cut on candidate pseudorapidity
    if (std::abs(hfCandXic.eta()) > cuts->get(pTBin, "eta")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), Eta);
    }

    // cut on pseudorapidity of final state daughters
    if (std::abs(etaPi0) > cuts->get(pTBin, "eta Daughters") || std::abs(etaPi1) > cuts->get(pTBin, "eta Daughters") || std::abs(etaPiFromXi) > cuts->get(pTBin, "eta Daughters") || std::abs(etaV0PosDau) > cuts->get(pTBin, "eta Daughters") || std::abs(etaV0NegDau) > cuts->get(pTBin, "eta Daughters")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), EtaDaughters);
    }

    // cut on pion pT
    if (hfCandXic.ptProng1() < cuts->get(pTBin, "pT Pi0") || hfCandXic.ptProng2() < cuts->get(pTBin, "pT Pi1")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), PtPionFromXicPlus);
    }

    // cut on chi2 of secondary vertex
    if (hfCandXic.chi2PCA() > cuts->get(pTBin, "chi2SV")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), Chi2SV);
    }

    // cut on candidate decay length
    if (hfCandXic.decayLength() < cuts->get(pTBin, "min decay length") || hfCandXic.decayLengthXY() < cuts->get(pTBin, "min decay length XY")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), MinDecayLength);
    }

    // cut on invariant mass of Xi-pion pairs
    if (hfCandXic.invMassXiPi0() > cuts->get(pTBin, "max inv mass Xi-Pi0") || hfCandXic.invMassXiPi1() > cuts->get(pTBin, "max inv mass Xi-Pi1")) {
      return false;
    }
    if (fillHistogram) {
      registry.fill(HIST("hSelCandidates"), MaxInvMassXiPiPairs);
    }

    return true;
  }

  /// Apply PID selection
  /// \param statusPidPi0   PID status of trackPi0 (prong1 of Xic candidate)
  /// \param statusPidPi1   PID status of trackPi1 (prong2 of Xic candidate)
  /// \param statusPidPiXi  PID status of trackPiXi (Bachelor of cascade candidate)
  /// \param statusPidPrLam PID status of trackPr (positive daughter of V0 candidate)
  /// \param statusPidPiLam PID status of trackPiLam (negative daughter of V0 candidate)
  /// \param usePidTpcOnly  switch to check only TPC status
  /// \return true if prongs of Xic candidate pass all selections
  bool selectionPid(TrackSelectorPID::Status const statusPidPi0,
                    TrackSelectorPID::Status const statusPidPi1,
                    TrackSelectorPID::Status const statusPidPiXi,
                    TrackSelectorPID::Status const statusPidPrLam,
                    TrackSelectorPID::Status const statusPidPiLam,
                    bool const useTpcPidOnly)
  {
    if (useTpcPidOnly) {
      if ((statusPidPi0 != TrackSelectorPID::Accepted && statusPidPi0 != TrackSelectorPID::NotApplicable) || (statusPidPi1 != TrackSelectorPID::Accepted && statusPidPi1 != TrackSelectorPID::NotApplicable) || (statusPidPiXi != TrackSelectorPID::Accepted && statusPidPiXi != TrackSelectorPID::NotApplicable) || (statusPidPrLam != TrackSelectorPID::Accepted && statusPidPrLam != TrackSelectorPID::NotApplicable) || (statusPidPiLam != TrackSelectorPID::Accepted && statusPidPiLam != TrackSelectorPID::NotApplicable)) {
        return false;
      }
      return true;
    }
    if (statusPidPi0 == TrackSelectorPID::Rejected || statusPidPi1 == TrackSelectorPID::Rejected || statusPidPiXi == TrackSelectorPID::Rejected || statusPidPrLam == TrackSelectorPID::Rejected || statusPidPiLam == TrackSelectorPID::Rejected) {
      return false;
    }
    return true;
  }

  void process(aod::HfCandXic const& hfCandsXic,
               TracksExtraWPid const&)
  {
    for (const auto& hfCandXic : hfCandsXic) {
      int statusXicToXiPiPi = 0;
      if (applyMl) {
        outputMlXicToXiPiPi.clear();
      }

      float ptCandXic = hfCandXic.pt();
      auto trackPi0 = hfCandXic.pi0_as<TracksExtraWPid>();
      auto trackPi1 = hfCandXic.pi1_as<TracksExtraWPid>();
      auto trackPiFromXi = hfCandXic.bachelor_as<TracksExtraWPid>();
      auto trackV0PosDau = hfCandXic.posTrack_as<TracksExtraWPid>();
      auto trackV0NegDau = hfCandXic.negTrack_as<TracksExtraWPid>();

      // Succesful reconstruction
      SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoTotal); // RecoTotal = 0 --> statusXicToXiPiPi += 1

      // kinematic and topological selection
      if (!isSelectedXic(hfCandXic, trackPi0.eta(), trackPi1.eta(), trackPiFromXi.eta(), trackV0PosDau.eta(), trackV0NegDau.eta())) {
        hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
        if (applyMl) {
          hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
        }
        continue;
      }
      SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoKinTopol); // RecoKinTopol = 1 --> statusXicToXiPiPi += 2

      // track quality selection
      if (doTrackQualitySelection) {
        if (!isSelectedTrackTpcQuality(trackPiFromXi, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax) ||
            !isSelectedTrackTpcQuality(trackV0PosDau, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax) ||
            !isSelectedTrackTpcQuality(trackV0NegDau, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          if (applyMl) {
            hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
          }
          continue;
        }
        if (fillHistogram) {
          registry.fill(HIST("hSelCandidates"), TpcTrackQualityXiDaughters);
        }

        if (!isSelectedTrackTpcQuality(trackPi0, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax) ||
            !isSelectedTrackTpcQuality(trackPi1, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          if (applyMl) {
            hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
          }
          continue;
        }
        if (fillHistogram) {
          registry.fill(HIST("hSelCandidates"), TpcTrackQualityPiFromCharm);
        }

        if ((!isSelectedTrackItsQuality(trackPi0, nClustersItsMin, itsChi2PerClusterMax) || trackPi0.itsNClsInnerBarrel() < nClustersItsInnBarrMin) ||
            (!isSelectedTrackItsQuality(trackPi0, nClustersItsMin, itsChi2PerClusterMax) || trackPi1.itsNClsInnerBarrel() < nClustersItsInnBarrMin)) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          if (applyMl) {
            hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
          }
          continue;
        }
        if (fillHistogram) {
          registry.fill(HIST("hSelCandidates"), ItsTrackQualityPiFromCharm);
        }

        SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoTrackQuality); // RecoTrackQuality = 2 --> statusXicToXiPiPi += 4
      }
      if (!doTrackQualitySelection && fillHistogram) {
        registry.fill(HIST("hSelCandidates"), TpcTrackQualityXiDaughters);
        registry.fill(HIST("hSelCandidates"), TpcTrackQualityPiFromCharm);
        registry.fill(HIST("hSelCandidates"), ItsTrackQualityPiFromCharm);
      }

      // track-level PID selection
      if (usePid) {
        TrackSelectorPID::Status statusPidPi0 = TrackSelectorPID::NotApplicable;
        TrackSelectorPID::Status statusPidPi1 = TrackSelectorPID::NotApplicable;
        TrackSelectorPID::Status statusPidPiXi = TrackSelectorPID::NotApplicable;
        TrackSelectorPID::Status statusPidPrLam = TrackSelectorPID::NotApplicable;
        TrackSelectorPID::Status statusPidPiLam = TrackSelectorPID::NotApplicable;

        // assign proton and pion hypothesis to V0 daughters
        auto trackPr = trackV0PosDau;
        auto trackPiFromLam = trackV0NegDau;
        if (hfCandXic.sign() < 0) {
          trackPr = trackV0NegDau;
          trackPiFromLam = trackV0PosDau;
        }

        if (useTpcPidOnly) {
          statusPidPi0 = selectorPion.statusTpc(trackPi0);
          statusPidPi1 = selectorPion.statusTpc(trackPi1);
          statusPidPiXi = selectorPion.statusTpc(trackPiFromXi);
          statusPidPrLam = selectorProton.statusTpc(trackPr);
          statusPidPiLam = selectorPion.statusTpc(trackPiFromLam);
        } else {
          statusPidPi0 = selectorPion.statusTpcOrTof(trackPi0);
          statusPidPi1 = selectorPion.statusTpcOrTof(trackPi1);
          statusPidPiXi = selectorPion.statusTpcOrTof(trackPiFromXi);
          statusPidPrLam = selectorProton.statusTpcOrTof(trackPr);
          statusPidPiLam = selectorPion.statusTpcOrTof(trackPiFromLam);
        }

        if (!selectionPid(statusPidPi0, statusPidPi1, statusPidPiXi, statusPidPrLam, statusPidPiLam, useTpcPidOnly.value)) {
          hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
          if (applyMl) {
            hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
          }
          continue;
        }
        SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoPID); // RecoPID = 3 --> statusXicToXiPiPi += 8
        if (fillHistogram) {
          registry.fill(HIST("hSelCandidates"), PidSelected);
        }
      }
      if (!usePid && fillHistogram) {
        registry.fill(HIST("hSelCandidates"), PidSelected);
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
        SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoMl); // RecoPID = 4 --> statusXicToXiPiPi += 16
        if (fillHistogram) {
          registry.fill(HIST("hSelCandidates"), BdtSelected);
        }
      }
      if (!applyMl && fillHistogram) {
        registry.fill(HIST("hSelCandidates"), BdtSelected);
      }

      hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorXicToXiPiPi>(cfgc)};
}
