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

/// \file candidateSelectorDstarToD0Pi.cxx
/// \brief Selection on D*± → D0(bar) π±  decay candidates
///
/// \author Deependra Sharma <deependra.sharma@cern.ch>, IITB
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <algorithm>
#include <string>
#include <vector>

// O2
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"
// O2Physics
#include "Common/Core/TrackSelectorPID.h"
// PWGHF
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseDstarToD0Pi.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsAnalysis.h"

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::analysis;
using namespace o2::framework;

// Struct to applying Dstar selection cuts
struct HfCandidateSelectorDstarToD0Pi {
  Produces<aod::HfSelDstarToD0Pi> hfSelDstarCandidate;
  Produces<aod::HfMlDstarToD0Pi> hfMlDstarCandidate;

  // Configurable specific to D0
  Configurable<double> ptD0CandMin{"ptD0CandMin", 0., "Minimum D0 candidate pT"};
  Configurable<double> ptD0CandMax{"ptD0CandMax", 50., "Maximum D0 candidate pT"};
  Configurable<std::vector<double>> binsPtD0{"binsPtD0", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits for D0"};
  Configurable<LabeledArray<double>> cutsD0{"cutsD0", {hf_cuts_d0_to_pi_k::Cuts[0], hf_cuts_d0_to_pi_k::NBinsPt, hf_cuts_d0_to_pi_k::NCutVars, hf_cuts_d0_to_pi_k::labelsPt, hf_cuts_d0_to_pi_k::labelsCutVar}, "D0 candidate selection per pT bin"};
  // Mass Cut for trigger analysis
  Configurable<bool> useTriggerMassCut{"useTriggerMassCut", false, "Flag to enable parametrize pT differential mass cut for triggered data"};

  // Configurable specific to Dstar
  Configurable<double> ptDstarCandMin{"ptDstarCandMin", 0., "Minimum Dstar candidate pT"};
  Configurable<double> ptDstarCandMax{"ptDstarCandMax", 50., "Maximum Dstar candidate pT"};
  Configurable<std::vector<double>> binsPtDstar{"binsPtDstar", std::vector<double>{hf_cuts_dstar_to_d0_pi::vecBinsPt}, "pT bin limits for Dstar"};
  Configurable<LabeledArray<double>> cutsDstar{"cutsDstar", {hf_cuts_dstar_to_d0_pi::Cuts[0], hf_cuts_dstar_to_d0_pi::NBinsPt, hf_cuts_dstar_to_d0_pi::NCutVars, hf_cuts_dstar_to_d0_pi::labelsPt, hf_cuts_dstar_to_d0_pi::labelsCutVar}, "Dstar candidate selection per pT bin"};

  // common Configurable
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Minimum track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Maximum track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};

  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Minimum track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 5., "Maximum track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};

  // AND logic for TOF+TPC PID
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Use AND logic for TPC and TOF PID"};

  // selecting only background candidates
  Configurable<bool> keepOnlySidebandCandidates{"keepOnlySidebandCandidates", false, "Select only sideband candidates, for studying background cut variable distributions"};
  Configurable<double> distanceFromDeltaMassForSidebands{"distanceFromDeltaMassForSidebands", 0.05, "Minimum distance from nominal (D*-D0) mass value for sideband region"};

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
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{""}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"Model.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  // PDG mass for kaon, pion and D0
  double massD0, massPi, massK;

  HfHelper hfHelper;
  o2::analysis::HfMlResponseDstarToD0Pi<float> hfMlResponse;
  std::vector<float> outputMlDstarToD0Pi = {};
  o2::ccdb::CcdbApi ccdbApi;

  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;

  using TracksSel = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;
  // using TracksSel = soa::Join<aod::Tracks, aod::TracksPidPi, aod::TracksPidKa>;
  using HfFullDstarCandidate = soa::Join<aod::HfD0FromDstar, aod::HfCandDstars>;

  AxisSpec axisBdtScore{100, 0.f, 1.f};
  AxisSpec axisSelStatus{2, -0.5f, 1.5f};
  HistogramRegistry registry{"registry"};

  HfTrigger2ProngCuts hfTriggerCuts;

  void init(InitContext&)
  {
    massPi = MassPiPlus;
    massK = MassKPlus;
    massD0 = MassD0;

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
      registry.add("QA/hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPtDstar, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("QA/hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }

      if (applyMl) {
        registry.add("QA/hBdtScore1VsStatus", ";BDT score", {HistType::kTH1F, {axisBdtScore}});
        registry.add("QA/hBdtScore2VsStatus", ";BDT score", {HistType::kTH1F, {axisBdtScore}});
        registry.add("QA/hBdtScore3VsStatus", ";BDT score", {HistType::kTH1F, {axisBdtScore}});
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

  /// Conjugate-independent topological cuts on D0
  /// @brief Topological selection on D0 candidate from Dstar
  /// @tparam T table iterator type of the candidate
  /// @param candidate candidate instance(object)
  /// @return true or false depending on selected or not
  template <typename T>
  bool selectionTopolD0(const T& candidate)
  {
    auto candpT = candidate.ptD0();
    auto binPt = findBin(binsPtD0, candpT);
    if (binPt == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptD0CandMin || candpT >= ptD0CandMax) {
      return false;
    }
    // product of daughter impact parameters
    if (candidate.impactParameterProductD0() > cutsD0->get(binPt, "d0d0")) {
      return false;
    }
    // cosine of pointing angle
    if (candidate.cpaD0() < cutsD0->get(binPt, "cos pointing angle")) {
      return false;
    }
    // cosine of pointing angle XY
    if (candidate.cpaXYD0() < cutsD0->get(binPt, "cos pointing angle xy")) {
      return false;
    }
    // normalised decay length in XY plane
    if (candidate.decayLengthXYNormalisedD0() < cutsD0->get(binPt, "min norm decay length XY")) {
      return false;
    }

    // Note: follwoing two cuts are not defined in namespace: hf_cuts_d0_to_pi_k of  SelectionCuts.h, while are defined in namespace: hf_cuts_dstar_to_d0_pi
    // Chi2PCA of secondary vertex reconstruction
    if (candidate.chi2PCAD0() > cutsDstar->get(binPt, "chi2PCA")) {
      return false;
    }
    if (std::abs(candidate.impactParameterXYD0()) > cutsD0->get(binPt, "DCA")) {
      return false;
    }
    // d0Prong0Normalised,1
    if (std::abs(candidate.impactParameterNormalised0()) < cutsDstar->get(binPt, "d0Prong0Normalised") || std::abs(candidate.impactParameterNormalised1()) < cutsDstar->get(binPt, "d0Prong1Normalised")) {
      return false;
    }

    // decay exponentail law, with tau = beta*gamma*ctau
    // decay length > ctau retains (1-1/e)

    if (candidate.decayLengthD0() < cutsD0->get(binPt, "min decay length")) {
      return false;
    }
    if (candidate.decayLengthD0() > cutsD0->get(binPt, "max decay length")) {
      return false;
    }
    if (candidate.decayLengthXYD0() > cutsD0->get(binPt, "max decay length XY")) {
      return false;
    }

    //.............Why is this if condition commented?
    if (candidate.decayLengthNormalisedD0() * candidate.decayLengthNormalisedD0() < 1.0) {
      // return false; // add back when getter fixed
    }
    return true;
  }

  /// Conjugate-independent topological cuts on Dstar
  /// @brief selection on Dstar candidate
  /// @tparam T table iterator type of the candidate
  /// @param candidate object
  /// @return true or false depending on selected or not
  template <typename T>
  bool selectionDstar(const T& candidate)
  {
    auto ptDstar = candidate.pt();
    auto binPt = findBin(binsPtDstar, ptDstar);
    if (binPt == -1) {
      return false;
    }
    // check that the candidate pT is within the analysis range
    if (ptDstar < ptDstarCandMin || ptDstar >= ptDstarCandMax) {
      return false;
    }
    // selction on DCA of softpi
    if (std::abs(candidate.impParamSoftPi()) > cutsDstar->get(binPt, "d0SoftPi")) {
      return false;
    }
    if (std::abs(candidate.normalisedImpParamSoftPi()) > cutsDstar->get(binPt, "d0SoftPiNormalised")) {
      return false;
    }
    // selection on pT of soft Pi
    if ((candidate.ptSoftPi() < cutsDstar->get(binPt, "ptSoftPiMin")) || (candidate.ptSoftPi() > cutsDstar->get(binPt, "ptSoftPiMax"))) {
      return false;
    }

    // selection on D0Candidate
    if (!selectionTopolD0(candidate)) {
      return false;
    }
    return true;
  }

  /// @brief Conjugate-dependent topological cuts
  /// @tparam T Table iterator type of candidate
  /// @param candidate candidate object
  /// @return true/false depending on success of selection
  template <typename T>
  bool selectionTopolConjugate(const T& candidate)
  {
    auto ptDstar = candidate.pt();
    auto binPt = findBin(binsPtDstar, ptDstar);
    if (binPt == -1) {
      return false;
    }
    auto prongSoftPi = candidate.template prongPi_as<TracksSel>();
    double mInvDstar = -999., mInvD0 = -999., mInvAntiDstar = -999., mInvD0Bar = -999.;

    if (prongSoftPi.sign() > 0.) { // Selection of D*+
      mInvDstar = candidate.invMassDstar();
      mInvD0 = candidate.invMassD0();
      if (std::abs(mInvD0 - massD0) > cutsD0->get(binPt, "m")) {
        return false;
      }
      if (useTriggerMassCut && !isCandidateInMassRange(mInvD0, massD0, candidate.ptD0(), hfTriggerCuts)) {
        return false;
      }
      // cut on daughter pT
      auto d0prong0 = candidate.template prong0_as<TracksSel>();
      auto d0prong1 = candidate.template prong1_as<TracksSel>();
      if (d0prong0.pt() < cutsD0->get(binPt, "pT Pi") || d0prong1.pt() < cutsD0->get(binPt, "pT K")) {
        return false;
      }
      // cut on difference of Dstar and D0 invariant mass
      if (std::abs(mInvDstar - mInvD0) > cutsDstar->get(binPt, "deltaMInvDstar")) {
        return false;
      }
      // cut on D0 daughter DCA - need to add secondary vertex constraint here
      if (std::abs(candidate.impactParameter0()) > cutsD0->get(binPt, "d0pi") || std::abs(candidate.impactParameter1()) > cutsD0->get(binPt, "d0K")) {
        return false;
      }
      // cut on cos(theta*)
      if (std::abs(candidate.cosThetaStarD0()) > cutsD0->get(binPt, "cos theta*")) {
        return false;
      }

    } else if (prongSoftPi.sign() < 0.) { // Selection of D*-
      mInvAntiDstar = candidate.invMassAntiDstar();
      mInvD0Bar = candidate.invMassD0Bar();
      if (std::abs(mInvD0Bar - massD0) > cutsD0->get(binPt, "m")) {
        return false;
      }
      if (useTriggerMassCut && !isCandidateInMassRange(mInvD0Bar, massD0, candidate.ptD0(), hfTriggerCuts)) {
        return false;
      }
      // cut on daughter pT
      auto d0prong0 = candidate.template prong0_as<TracksSel>();
      auto d0prong1 = candidate.template prong1_as<TracksSel>();
      if (d0prong0.pt() < cutsD0->get(binPt, "pT K") || d0prong1.pt() < cutsD0->get(binPt, "pT Pi")) {
        return false;
      }
      // cut on difference of Dstar and D0 invariant mass
      if (std::abs(mInvAntiDstar - mInvD0Bar) > cutsDstar->get(binPt, "deltaMInvDstar")) {
        return false;
      }
      // cut on D0 daughter DCA - need to add secondary vertex constraint here
      if (std::abs(candidate.impactParameter0()) > cutsD0->get(binPt, "d0K") || std::abs(candidate.impactParameter1()) > cutsD0->get(binPt, "d0pi")) {
        return false;
      }
      // cut on cos(theta*)
      if (std::abs(candidate.cosThetaStarD0Bar()) > cutsD0->get(binPt, "cos theta*")) {
        return false;
      }
    }

    // in case only sideband candidates have to be stored, additional invariant-mass cut
    if (keepOnlySidebandCandidates && prongSoftPi.sign() > 0.) {
      if (std::abs((mInvDstar - mInvD0) - massPi) < distanceFromDeltaMassForSidebands) {
        return false;
      }
    } else if (keepOnlySidebandCandidates && prongSoftPi.sign() < 0.) {
      if (std::abs((mInvAntiDstar - mInvD0Bar) - massPi) < distanceFromDeltaMassForSidebands) {
        return false;
      }
    }

    return true;
  }

  void process(TracksSel const&,
               HfFullDstarCandidate const& rowsDstarCand)
  {
    // LOG(info) << "selector called";
    for (const auto& candDstar : rowsDstarCand) {
      // final selection flag: false - rejected, true - accepted
      bool statusDstar = false, statusD0Flag = false, statusTopol = false, statusCand = false, statusPID = false;

      outputMlDstarToD0Pi.clear();
      auto ptCand = candDstar.pt();

      auto prongPi = candDstar.prongPi_as<TracksSel>();
      auto prong0 = candDstar.prong0_as<TracksSel>();
      auto prong1 = candDstar.prong1_as<TracksSel>();

      if (!TESTBIT(candDstar.hfflag(), aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlDstarCandidate(outputMlDstarToD0Pi);
        }
        if (activateQA) {
          registry.fill(HIST("QA/hSelections"), 1, ptCand);
        }
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("QA/hSelections"), 2 + aod::SelectionStep::RecoSkims, ptCand);
      }
      statusD0Flag = true;

      if (!selectionDstar(candDstar)) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlDstarCandidate(outputMlDstarToD0Pi);
        }
        continue;
      }
      statusTopol = true;
      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      // conjugate-dependent topological selection for Dstar
      bool topoDstar = selectionTopolConjugate(candDstar);
      if (!topoDstar) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlDstarCandidate(outputMlDstarToD0Pi);
        }
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("QA/hSelections"), 2 + aod::SelectionStep::RecoTopol, ptCand);
      }
      statusCand = true;

      // track-level PID selection
      int pidTrackPosKaon = -1;
      int pidTrackPosPion = -1;
      int pidTrackNegKaon = -1;
      int pidTrackNegPion = -1;

      if (usePidTpcAndTof) {
        pidTrackPosKaon = selectorKaon.statusTpcAndTof(candDstar.prong0_as<TracksSel>());
        pidTrackPosPion = selectorPion.statusTpcAndTof(candDstar.prong0_as<TracksSel>());
        pidTrackNegKaon = selectorKaon.statusTpcAndTof(candDstar.prong1_as<TracksSel>());
        pidTrackNegPion = selectorPion.statusTpcAndTof(candDstar.prong1_as<TracksSel>());
      } else {
        pidTrackPosKaon = selectorKaon.statusTpcOrTof(candDstar.prong0_as<TracksSel>());
        pidTrackPosPion = selectorPion.statusTpcOrTof(candDstar.prong0_as<TracksSel>());
        pidTrackNegKaon = selectorKaon.statusTpcOrTof(candDstar.prong1_as<TracksSel>());
        pidTrackNegPion = selectorPion.statusTpcOrTof(candDstar.prong1_as<TracksSel>());
      }

      int pidDstar = -1;
      if (candDstar.prongPi_as<TracksSel>().sign() > 0.) {
        if (pidTrackPosPion == TrackSelectorPID::Accepted && pidTrackNegKaon == TrackSelectorPID::Accepted) {
          pidDstar = 1; // accept D*+
        } else if (pidTrackPosPion == TrackSelectorPID::Rejected && pidTrackNegKaon == TrackSelectorPID::Rejected) {
          pidDstar = 0; // reject D*+
        }
      } else if (candDstar.prongPi_as<TracksSel>().sign() < 0.) {
        if (pidTrackNegPion == TrackSelectorPID::Accepted && pidTrackPosKaon == TrackSelectorPID::Accepted) {
          pidDstar = 1; // Accept D*-
        } else if (pidTrackNegPion == TrackSelectorPID::Rejected && pidTrackPosKaon == TrackSelectorPID::Rejected) {
          pidDstar = 0; // reject D*-
        }
      }

      if (pidDstar == 0) {
        hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
        if (applyMl) {
          hfMlDstarCandidate(outputMlDstarToD0Pi);
        }
        continue;
      }

      if ((pidDstar == -1 || pidDstar == 1) && topoDstar) {
        statusDstar = true; // identified as dstar
      }

      if (activateQA) {
        registry.fill(HIST("QA/hSelections"), 2 + aod::SelectionStep::RecoPID, ptCand);
      }
      statusPID = true;

      if (applyMl) {
        // ML selections
        bool isSelectedMlDstar = false;

        std::vector<float> inputFeatures = hfMlResponse.getInputFeatures(candDstar, prong0, prong1, prongPi);
        isSelectedMlDstar = hfMlResponse.isSelectedMl(inputFeatures, ptCand, outputMlDstarToD0Pi);

        hfMlDstarCandidate(outputMlDstarToD0Pi);
        if (activateQA) {
          registry.fill(HIST("QA/hBdtScore1VsStatus"), outputMlDstarToD0Pi[0]);
          registry.fill(HIST("QA/hBdtScore2VsStatus"), outputMlDstarToD0Pi[1]);
          registry.fill(HIST("QA/hBdtScore3VsStatus"), outputMlDstarToD0Pi[2]);
        }

        if (!isSelectedMlDstar) {
          statusDstar = false;
          hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
          continue;
        }

        if (activateQA) {
          registry.fill(HIST("QA/hSelections"), 2 + aod::SelectionStep::RecoMl, ptCand);
        }
      }
      hfSelDstarCandidate(statusDstar, statusD0Flag, statusTopol, statusCand, statusPID);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorDstarToD0Pi>(cfgc)};
}
