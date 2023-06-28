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

/// \file candidateSelectorLcToK0sP.cxx
/// \brief Lc --> K0s+p selection task.
/// \note based on candidateSelectorD0.cxx
///
/// \author Chiara Zampolli <Chiara.Zampolli@cern.ch>, CERN
///         Daniel Samitz, <daniel.samitz@cern.ch>, Vienna

#include "Common/Core/TrackSelectorPID.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsDebugLcToK0sP.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc;
using namespace o2::analysis::hf_cuts_lc_to_k0s_p;

using MyBigTracksBayes = soa::Join<aod::BigTracksPID, aod::pidBayesPr, aod::pidBayesEl, aod::pidBayesMu, aod::pidBayesKa, aod::pidBayesPi>;

struct HfCandidateSelectorLcToK0sP {
  Produces<aod::HfSelLcToK0sP> hfSelLcToK0sPCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};

  // PID
  Configurable<double> pPidThreshold{"pPidThreshold", 1.0, "Threshold to switch between low and high p TrackSelectors"};
  // TPC
  Configurable<double> nSigmaTpcMaxLowP{"nSigmaTpcMaxLowP", 2.0, "Max nSigma in TPC for bachelor at low p"};
  Configurable<double> nSigmaTpcCombinedMaxLowP{"nSigmaTpcCombinedMaxLowP", 2.0, "Max nSigma in TPC combined with TOF for bachelor at low p"};
  Configurable<double> nSigmaTpcMaxHighP{"nSigmaTpcMaxHighP", 0., "Max nSigma in TPC for bachelor at high p"};
  Configurable<double> nSigmaTpcCombinedMaxHighP{"nSigmaTpcCombinedMaxHighP", 9999., "Max nSigma in TPC combined with TOF for bachelor at high p"};
  // TOF
  Configurable<double> nSigmaTofMaxLowP{"nSigmaTofMaxLowP", 0., "Max nSigma in TOF for bachelor at low p"};
  Configurable<double> nSigmaTofCombinedMaxLowP{"nSigmaTofCombinedMaxLowP", 9999., "Max nSigma in TOF combined with TPC for bachelor at low p"};
  Configurable<double> nSigmaTofMaxHighP{"nSigmaTofMaxHighP", 3.0, "Max nSigma in TOF for bachelor at high p"};
  Configurable<double> nSigmaTofCombinedMaxHighP{"nSigmaTofCombinedMaxHighP", 3.0, "Max nSigma in TOF combined with TPC for bachelor at high p"};
  // Bayesian
  Configurable<double> probBayesMinLowP{"probBayesMinLowP", 0.8, "min. Bayes probability for bachelor at low p [%]"};
  Configurable<double> probBayesMinHighP{"probBayesMinHighP", 0.8, "min. Bayes probability for bachelor at high p [%]"};

  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_k0s_p::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Lc candidate selection per pT bin"};

  void init(InitContext&)
  {
    if (!doprocessWithStandardPID && !doprocessWithBayesPID) {
      LOGF(fatal, "Neither processWithStandardPID nor processWithBayesPID enabled. Please choose one.");
    }
    if (doprocessWithStandardPID && doprocessWithBayesPID) {
      LOGF(fatal, "Cannot enable processWithStandardPID and processWithBayesPID at the same time. Please choose one.");
    }
  }

  /// Conjugate independent toplogical cuts
  /// \param hfCandCascade is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& hfCandCascade)
  {
    auto candPt = hfCandCascade.pt();
    int ptBin = findBin(binsPt, candPt);
    if (ptBin == -1) {
      return false;
    }

    if (candPt < ptCandMin || candPt >= ptCandMax) {
      return false; // check that the candidate pT is within the analysis range
    }

    if (std::abs(hfCandCascade.mK0Short() - RecoDecay::getMassPDG(kK0Short)) > cuts->get(ptBin, "mK0s")) {
      return false; // mass of the K0s
    }

    if ((std::abs(hfCandCascade.mLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda")) || (std::abs(hfCandCascade.mAntiLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda"))) {
      return false; // mass of the Lambda
    }

    if (std::abs(hfCandCascade.mGamma() - RecoDecay::getMassPDG(kGamma)) < cuts->get(ptBin, "mGamma")) {
      return false; // mass of the Gamma
    }

    if (hfCandCascade.ptProng0() < cuts->get(ptBin, "ptBach")) {
      return false; // pt of the p
    }

    if (hfCandCascade.ptV0Pos() < cuts->get(ptBin, "ptV0Dau")) {
      return false; // pt of the pos K0 daughter
    }

    if (hfCandCascade.ptV0Neg() < cuts->get(ptBin, "ptV0Dau")) {
      return false; // pt of the neg K0 daughter
    }

    if (hfCandCascade.ptProng1() < cuts->get(ptBin, "ptV0")) {
      return false; // pt of the K0
    }

    if (std::abs(hfCandCascade.impactParameter0()) > cuts->get(ptBin, "d0Bach")) {
      return false; // d0 of the bachelor
    }

    /*
    if ((std::abs(hfCandCascade.dcapostopv()) > d0K0Cut[ptBin]) || (std::abs(hfCandCascade.dcanegtopv()) > d0K0Cut[ptBin])) {
      LOG(debug) << "v0 daugh cut failed, positive v0 daugh --> " << hfCandCascade.dcapostopv() << ", negative v0 daugh --> " << hfCandCascade.dcanegtopv() << " , cut --> " << d0K0Cut[ptBin];
      return false; // d0 of the K0s daughters
    }
    */

    if (std::abs(hfCandCascade.impactParameter1()) > cuts->get(ptBin, "d0V0")) {
      return false; // d0 of the v0
    }

    return true;
  }

  template <typename T>
  bool selectionStandardPID(const T& track)
  {
    TrackSelectorPID selectorProton = TrackSelectorPID(kProton);
    if (track.p() < pPidThreshold) {
      selectorProton.setRangeNSigmaTPC(-nSigmaTpcMaxLowP, nSigmaTpcMaxLowP);
      selectorProton.setRangeNSigmaTOF(-nSigmaTofMaxLowP, nSigmaTofMaxLowP);
      selectorProton.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxLowP, nSigmaTpcCombinedMaxLowP);
      selectorProton.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxLowP, nSigmaTofCombinedMaxLowP);
    } else {
      selectorProton.setRangeNSigmaTPC(-nSigmaTpcMaxHighP, nSigmaTpcMaxHighP);
      selectorProton.setRangeNSigmaTOF(-nSigmaTofMaxHighP, nSigmaTofMaxHighP);
      selectorProton.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMaxHighP, nSigmaTpcCombinedMaxHighP);
      selectorProton.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMaxHighP, nSigmaTofCombinedMaxHighP);
    }

    return selectorProton.getStatusTrackPIDTpcAndTof(track) == TrackSelectorPID::Status::PIDAccepted;
  }

  template <typename T>
  bool selectionBayesPID(const T& track)
  {
    if (!selectionStandardPID(track)) { // possibility to add some pre-selection before using Bayesian PID
      return false;
    }

    TrackSelectorPID selectorProton = TrackSelectorPID(kProton);
    if (track.p() < pPidThreshold) {
      selectorProton.setProbBayesMin(probBayesMinLowP);
    } else {
      selectorProton.setProbBayesMin(probBayesMinHighP);
    }

    return selectorProton.getStatusTrackBayesProbPID(track) == TrackSelectorPID::Status::PIDAccepted;
  }

  void processWithStandardPID(aod::HfCandCascade const& candidates, aod::BigTracksPID const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {                     // looping over cascade candidates
      const auto& bach = candidate.prong0_as<aod::BigTracksPID>(); // bachelor track

      statusLc = 0;

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)
      if (!selectionTopol(candidate)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      if (!selectionStandardPID(bach)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      statusLc = 1;

      hfSelLcToK0sPCandidate(statusLc);
    }
  }
  PROCESS_SWITCH(HfCandidateSelectorLcToK0sP, processWithStandardPID, "Use standard PID for bachelor track", true);

  void processWithBayesPID(aod::HfCandCascade const& candidates, MyBigTracksBayes const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {                    // looping over cascade candidates
      const auto& bach = candidate.prong0_as<MyBigTracksBayes>(); // bachelor track

      statusLc = 0;

      if (!selectionTopol(candidate)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      if (!selectionBayesPID(bach)) {
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      statusLc = 1;

      hfSelLcToK0sPCandidate(statusLc);
    }
  }
  PROCESS_SWITCH(HfCandidateSelectorLcToK0sP, processWithBayesPID, "Use Bayesian PID for bachelor track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorLcToK0sP>(cfgc)};
}
