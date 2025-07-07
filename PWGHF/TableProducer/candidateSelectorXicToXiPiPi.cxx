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
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <Rtypes.h>

#include <cstdint>
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
                        BdtSelected,
                        NSelectionCriteria };

  using TracksExtraWPid = soa::Join<aod::TracksWExtra, aod::TracksPidPi, aod::TracksPidPr>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    if ((doprocessData + doprocessMc) != 1) {
      LOGP(fatal, "Enable exactly one process function at a time.");
    }

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
      std::string labels[NSelectionCriteria];
      labels[All] = "All";
      labels[Pt] = "#it{p}_{T}";
      labels[Mass] = "#Delta M";
      labels[Rapidity] = "y";
      labels[Eta] = "#eta";
      labels[EtaDaughters] = "#eta final state daughters";
      labels[PtPionFromXicPlus] = "#it{p}_{T} (#pi #leftarrow #Xi_{c}^{+})";
      labels[Chi2SV] = "#chi^{2}_{SV}";
      labels[MinDecayLength] = "Decay length";
      labels[MaxInvMassXiPiPairs] = "M_{#Xi #pi}";
      labels[TpcTrackQualityXiDaughters] = "TPC track quality selection on #Xi daughters";
      labels[TpcTrackQualityPiFromCharm] = "TPC track quality selection on #pi #leftarrow #Xi_{c}^{+}";
      labels[ItsTrackQualityPiFromCharm] = "ITS track quality selection on #pi #leftarrow #Xi_{c}^{+}";
      labels[PidSelected] = "PID selection";
      labels[BdtSelected] = "BDT selection";

      if (doprocessData) {
        registry.add("hSelCandidates", ";;entries", {HistType::kTH1F, {{NSelectionCriteria, -0.5f, +NSelectionCriteria - 0.5f}}});
        for (int iBin = 0; iBin < NSelectionCriteria; ++iBin) {
          registry.get<TH1>(HIST("hSelCandidates"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
        }
      } else if (doprocessMc) {
        registry.add("hSelCandidatesRecSig", ";;entries", {HistType::kTH1F, {{NSelectionCriteria, -0.5f, +NSelectionCriteria - 0.5f}}});
        registry.add("hSelCandidatesRecBkg", ";;entries", {HistType::kTH1F, {{NSelectionCriteria, -0.5f, +NSelectionCriteria - 0.5f}}});
        for (int iBin = 0; iBin < NSelectionCriteria; ++iBin) {
          registry.get<TH1>(HIST("hSelCandidatesRecSig"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
          registry.get<TH1>(HIST("hSelCandidatesRecBkg"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
        }
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

  /// Apply PID selection
  /// \param statusPidPi0   PID status of prong1 with pion hypothesis
  /// \param statusPidPi1   PID status of prong2 with pion hypothesis
  /// \param statusPidPiXi  PID status of bachelor track with pion hypothesis
  /// \param statusPidPrLam PID status of V0 daughter with proton hypothesis
  /// \param statusPidPiLam PID status of V0 daughter with pion hypothesis
  /// \param usePidTpcOnly  switch to check only TPC status
  /// \return true if prongs of Xic candidate pass all PID selections
  bool isSelectedPid(TrackSelectorPID::Status const statusPidPi0,
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

  /// Combine kinematic, topological, track quality and PID selections
  /// \param hfCandXic         Xic candidate
  /// \param statusXicToXiPiPi Flag to store selection status as defined in hf_sel_candidate_xic::XicToXiPiPiSelectionStep
  /// \param isMatchedSignal   Flag to indicate if the candidate is matched to a genereated XiCplus MC particle
  /// \return true if Xic candidate passes all selections, otherwise false
  template <bool isMc, typename XicCandidate>
  bool isSelectedXicToXiPiPiCandidateWoMl(XicCandidate const& hfCandXic,
                                          TracksExtraWPid const&,
                                          int& statusXicToXiPiPi,
                                          bool const isMatchedSignal = false)
  {
    // Successful reconstruction
    SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoTotal); // RecoTotal = 0 --> statusXicToXiPiPi += 1
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), All);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), All);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), All);
    }

    // Retrieve daughter tracks
    auto trackPi0 = hfCandXic.template pi0_as<TracksExtraWPid>();
    auto trackPi1 = hfCandXic.template pi1_as<TracksExtraWPid>();
    auto trackPiFromXi = hfCandXic.template bachelor_as<TracksExtraWPid>();
    auto trackV0PosDau = hfCandXic.template posTrack_as<TracksExtraWPid>();
    auto trackV0NegDau = hfCandXic.template negTrack_as<TracksExtraWPid>();

    ////////////////////////////////////////////////
    //     Kinematic and topological selection    //
    ////////////////////////////////////////////////

    // check whether candidate is in analyzed pT range
    auto ptCandXic = hfCandXic.pt();
    int pTBin = findBin(binsPt, ptCandXic);
    if (pTBin == -1) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), Pt);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), Pt);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), Pt);
    }
    // check whether candidate mass is within a defined mass window
    if (std::abs(hfCandXic.invMassXicPlus() - o2::constants::physics::MassXiCPlus) > cuts->get(pTBin, "m")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), Mass);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), Mass);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), Mass);
    }

    // cut on candidate rapidity
    if (std::abs(hfCandXic.y(o2::constants::physics::MassXiCPlus)) > cuts->get(pTBin, "y")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), Rapidity);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), Rapidity);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), Rapidity);
    }

    // cut on candidate pseudorapidity
    if (std::abs(hfCandXic.eta()) > cuts->get(pTBin, "eta")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), Eta);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), Eta);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), Eta);
    }

    // cut on pseudorapidity of final state daughters
    if (std::abs(trackPi0.eta()) > cuts->get(pTBin, "eta Daughters") || std::abs(trackPi1.eta()) > cuts->get(pTBin, "eta Daughters") || std::abs(trackPiFromXi.eta()) > cuts->get(pTBin, "eta Daughters") || std::abs(trackV0PosDau.eta()) > cuts->get(pTBin, "eta Daughters") || std::abs(trackV0NegDau.eta()) > cuts->get(pTBin, "eta Daughters")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), EtaDaughters);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), EtaDaughters);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), EtaDaughters);
    }

    // cut on pion pT
    if (hfCandXic.ptProng1() < cuts->get(pTBin, "pT Pi0") || hfCandXic.ptProng2() < cuts->get(pTBin, "pT Pi1")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), PtPionFromXicPlus);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), PtPionFromXicPlus);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), PtPionFromXicPlus);
    }

    // cut on chi2 of secondary vertex
    if (hfCandXic.chi2PCA() > cuts->get(pTBin, "chi2SV")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), Chi2SV);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), Chi2SV);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), Chi2SV);
    }

    // cut on candidate decay length
    if (hfCandXic.decayLength() < cuts->get(pTBin, "min decay length") || hfCandXic.decayLengthXY() < cuts->get(pTBin, "min decay length XY")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), MinDecayLength);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), MinDecayLength);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), MinDecayLength);
    }

    // cut on invariant mass of Xi-pion pairs
    if (hfCandXic.invMassXiPi0() > cuts->get(pTBin, "max inv mass Xi-Pi0") || hfCandXic.invMassXiPi1() > cuts->get(pTBin, "max inv mass Xi-Pi1")) {
      return false;
    }
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), MaxInvMassXiPiPairs);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), MaxInvMassXiPiPairs);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), MaxInvMassXiPiPairs);
    }

    // Successful kinematic and topological selection
    SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoKinTopol); // RecoKinTopol = 1 --> statusXicToXiPiPi += 2

    ////////////////////////////////////////////////
    //          Track quality selection           //
    ////////////////////////////////////////////////
    if (doTrackQualitySelection) {
      // TPC track quality selection on Xi daughters (bachelor pion and V0 daughters)
      if (!isSelectedTrackTpcQuality(trackPiFromXi, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax) ||
          !isSelectedTrackTpcQuality(trackV0PosDau, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax) ||
          !isSelectedTrackTpcQuality(trackV0NegDau, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
        return false;
      }
      if (fillHistogram && isMc && isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecSig"), TpcTrackQualityXiDaughters);
      } else if (fillHistogram && isMc && !isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecBkg"), TpcTrackQualityXiDaughters);
      } else if (fillHistogram && !isMc) {
        registry.fill(HIST("hSelCandidates"), TpcTrackQualityXiDaughters);
      }

      // TPC track quality selection on pions from charm baryon
      if (!isSelectedTrackTpcQuality(trackPi0, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax) ||
          !isSelectedTrackTpcQuality(trackPi1, nClustersTpcMin, nTpcCrossedRowsMin, tpcCrossedRowsOverFindableClustersRatioMin, tpcChi2PerClusterMax)) {
        return false;
      }
      if (fillHistogram && isMc && isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecSig"), TpcTrackQualityPiFromCharm);
      } else if (fillHistogram && isMc && !isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecBkg"), TpcTrackQualityPiFromCharm);
      } else if (fillHistogram && !isMc) {
        registry.fill(HIST("hSelCandidates"), TpcTrackQualityPiFromCharm);
      }

      // ITS track quality selection on pions from charm baryon
      if ((!isSelectedTrackItsQuality(trackPi0, nClustersItsMin, itsChi2PerClusterMax) || trackPi0.itsNClsInnerBarrel() < nClustersItsInnBarrMin) ||
          (!isSelectedTrackItsQuality(trackPi0, nClustersItsMin, itsChi2PerClusterMax) || trackPi1.itsNClsInnerBarrel() < nClustersItsInnBarrMin)) {
        return false;
      }
      if (fillHistogram && isMc && isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecSig"), ItsTrackQualityPiFromCharm);
      } else if (fillHistogram && isMc && !isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecBkg"), ItsTrackQualityPiFromCharm);
      } else if (fillHistogram && !isMc) {
        registry.fill(HIST("hSelCandidates"), ItsTrackQualityPiFromCharm);
      }

      // Successful track quality selection
      SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoTrackQuality); // RecoTrackQuality = 2 --> statusXicToXiPiPi += 4
    } else if (fillHistogram) {
      if (isMc && isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecSig"), TpcTrackQualityXiDaughters);
        registry.fill(HIST("hSelCandidatesRecSig"), TpcTrackQualityPiFromCharm);
        registry.fill(HIST("hSelCandidatesRecSig"), ItsTrackQualityPiFromCharm);
      } else if (isMc && !isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecBkg"), TpcTrackQualityXiDaughters);
        registry.fill(HIST("hSelCandidatesRecBkg"), TpcTrackQualityPiFromCharm);
        registry.fill(HIST("hSelCandidatesRecBkg"), ItsTrackQualityPiFromCharm);
      } else if (!isMc) {
        registry.fill(HIST("hSelCandidates"), TpcTrackQualityXiDaughters);
        registry.fill(HIST("hSelCandidates"), TpcTrackQualityPiFromCharm);
        registry.fill(HIST("hSelCandidates"), ItsTrackQualityPiFromCharm);
      }
    }

    ////////////////////////////////////////////////
    //                PID selection               //
    ////////////////////////////////////////////////
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

      if (!isSelectedPid(statusPidPi0, statusPidPi1, statusPidPiXi, statusPidPrLam, statusPidPiLam, useTpcPidOnly.value)) {
        return false;
      }

      // Successful PID selection
      SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoPID); // RecoPID = 3 --> statusXicToXiPiPi += 8
      if (fillHistogram && isMc && isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecSig"), PidSelected);
      } else if (fillHistogram && isMc && !isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecBkg"), PidSelected);
      } else if (fillHistogram && !isMc) {
        registry.fill(HIST("hSelCandidates"), PidSelected);
      }
    } else if (fillHistogram) {
      if (isMc && isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecSig"), PidSelected);
      } else if (isMc && !isMatchedSignal) {
        registry.fill(HIST("hSelCandidatesRecBkg"), PidSelected);
      } else {
        registry.fill(HIST("hSelCandidates"), PidSelected);
      }
    }

    return true;
  }

  /// Apply BDT selection
  /// \param hfCandXic         Xic candidate
  /// \param statusXicToXiPiPi Flag to store selection status as defined in hf_sel_candidate_xic::XicToXiPiPiSelectionStep
  /// \param isMatchedSignal   Flag to indicate if the candidate is matched to a genereated XiCplus MC particle
  template <bool isMc, typename XicCandidate>
  void isBdtSelected(XicCandidate const& hfCandXic,
                     int& statusXicToXiPiPi,
                     bool const isMatchedSignal = false)
  {
    bool isSelectedMlXicToXiPiPi = false;
    float ptCandXic = hfCandXic.pt();

    std::vector<float> inputFeaturesXicToXiPiPi = hfMlResponse.getInputFeatures(hfCandXic);
    isSelectedMlXicToXiPiPi = hfMlResponse.isSelectedMl(inputFeaturesXicToXiPiPi, ptCandXic, outputMlXicToXiPiPi);

    hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);

    if (!isSelectedMlXicToXiPiPi) {
      return;
    }

    // Successful ML selection
    SETBIT(statusXicToXiPiPi, hf_sel_candidate_xic::XicToXiPiPiSelectionStep::RecoMl); // RecoPID = 4 --> statusXicToXiPiPi += 16
    if (fillHistogram && isMc && isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecSig"), BdtSelected);
    } else if (fillHistogram && isMc && !isMatchedSignal) {
      registry.fill(HIST("hSelCandidatesRecBkg"), BdtSelected);
    } else if (fillHistogram && !isMc) {
      registry.fill(HIST("hSelCandidates"), BdtSelected);
    }
  }

  void processData(aod::HfCandXic const& hfCandsXic,
                   TracksExtraWPid const& tracks)
  {
    for (const auto& hfCandXic : hfCandsXic) {
      int statusXicToXiPiPi = 0;
      if (applyMl) {
        outputMlXicToXiPiPi.clear();
      }

      // Kinematic, topological, track quality and PID selections
      if (!isSelectedXicToXiPiPiCandidateWoMl<false>(hfCandXic, tracks, statusXicToXiPiPi)) {
        hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
        if (applyMl) {
          hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
        }
        continue;
      }

      // ML selection
      if (applyMl) {
        isBdtSelected<false>(hfCandXic, statusXicToXiPiPi);
      } else if (fillHistogram) {
        registry.fill(HIST("hSelCandidates"), BdtSelected);
      }

      hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
    }
  }
  PROCESS_SWITCH(HfCandidateSelectorXicToXiPiPi, processData, "Select candidates without MC matching information", true);

  void processMc(soa::Join<aod::HfCandXic, aod::HfCandXicMcRec> const& hfCandsXic,
                 TracksExtraWPid const& tracks)
  {
    for (const auto& hfCandXic : hfCandsXic) {
      int statusXicToXiPiPi = 0;
      if (applyMl) {
        outputMlXicToXiPiPi.clear();
      }

      bool isMatchedCandidate = false;
      if (hfCandXic.flagMcMatchRec() != int8_t(0)) {
        isMatchedCandidate = true;
      }

      // Kinematic, topological, track quality and PID selections
      if (!isSelectedXicToXiPiPiCandidateWoMl<true>(hfCandXic, tracks, statusXicToXiPiPi, isMatchedCandidate)) {
        hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
        if (applyMl) {
          hfMlXicToXiPiPiCandidate(outputMlXicToXiPiPi);
        }
        continue;
      }

      // ML selection
      if (applyMl) {
        isBdtSelected<true>(hfCandXic, statusXicToXiPiPi, isMatchedCandidate);
      } else if (fillHistogram && isMatchedCandidate) {
        registry.fill(HIST("hSelCandidatesRecSig"), BdtSelected);
      } else if (fillHistogram && !isMatchedCandidate) {
        registry.fill(HIST("hSelCandidatesRecBkg"), BdtSelected);
      }

      hfSelXicToXiPiPiCandidate(statusXicToXiPiPi);
    }
  }
  PROCESS_SWITCH(HfCandidateSelectorXicToXiPiPi, processMc, "Select candidates with MC matching information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorXicToXiPiPi>(cfgc)};
}
