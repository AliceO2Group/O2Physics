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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"
#include "PWGHF/Utils/utilsDebugLcToK0sP.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc;
using namespace o2::analysis::hf_cuts_lc_to_k0s_p;

//#define MY_DEBUG
#ifdef MY_DEBUG
#define MY_DEBUG_MSG(condition, cmd)
if (condition) {
  cmd;
}
using MyBigTracks = soa::Join<aod::BigTracksPID, aod::pidBayesPr, aod::pidBayesEl, aod::pidBayesMu, aod::pidBayesKa, aod::pidBayesPi, aod::McTrackLabels>;
#else
#define MY_DEBUG_MSG(condition, cmd)
using MyBigTracks = soa::Join<aod::BigTracksPID, aod::pidBayesPr, aod::pidBayesEl, aod::pidBayesMu, aod::pidBayesKa, aod::pidBayesPi>;
#endif

struct HfCandidateSelectorLcToK0sP {
  Produces<aod::HfSelLcToK0sP> hfSelLcToK0sPCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // PID
  Configurable<bool> usePidBach{"usePidBach", true, "Use PID on bachelor track"};
  Configurable<double> pPidThreshold{"pPidThreshold", 1.0, "Threshold to switch between low and high p TrackSelectors"};
  Configurable<double> nSigmaTpcMaxLowP{"nSigmaTpcMaxLowP", 2.0, "Max nSigam in TPC for bachelor at low p"};
  Configurable<double> nSigmaTpcMaxHighP{"nSigmaTpcMaxHighP", 9999., "Max nSigam in TPC for bachelor at high p"};
  Configurable<double> nSigmaTofMaxLowP{"nSigmaTofMaxLowP", 9999., "Max nSigam in TOF for bachelor at low p"};
  Configurable<double> nSigmaTofMaxHighP{"nSigmaTofMaxHighP", 3.0, "Max nSigam in TOF for bachelor at high p"};
  Configurable<bool> requireTofLowP{"requireTofLowP", false, "require TOF information for bachelor at low p"};
  Configurable<bool> requireTofHighP{"requireTofHighP", true, "require TOF information for bachelor at high p"};
  Configurable<double> probBayesMinLowP{"probBayesMinLowP", -1., "min. Bayes probability for bachelor at low p [%]"};
  Configurable<double> probBayesMinHighP{"probBayesMinHighP", -1., "min. Bayes probability for bachelor at high p [%]"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_k0s_p::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_k0s_p::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Lc candidate selection per pT bin"};

  // for debugging
#ifdef MY_DEBUG
  Configurable<std::vector<int>> indexK0Spos{"indexK0Spos", {729, 2866, 4754, 5457, 6891, 7824, 9243, 9810}, "indices of K0S positive daughters, for debug"};
  Configurable<std::vector<int>> indexK0Sneg{"indexK0Sneg", {730, 2867, 4755, 5458, 6892, 7825, 9244, 9811}, "indices of K0S negative daughters, for debug"};
  Configurable<std::vector<int>> indexProton{"indexProton", {717, 2810, 4393, 5442, 6769, 7793, 9002, 9789}, "indices of protons, for debug"};
#endif

  /// Selection on goodness of daughter tracks
  /// \note should be applied at candidate selection
  /// \param track is daughter track
  /// \return true if track is good
  template <typename T>
  bool daughterSelection(const T& track) // there could be more checks done on the bachelor here
  {
    // this is for now just a placeholder, in case we want to add extra checks
    return true;
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
      LOG(debug) << "cand pt (first check) cut failed: from cascade --> " << candPt << ", cut --> " << ptCandMax;
      return false; // check that the candidate pT is within the analysis range
    }

    if (std::abs(hfCandCascade.mK0Short() - RecoDecay::getMassPDG(kK0Short)) > cuts->get(ptBin, "mK0s")) {
      LOG(debug) << "massK0s cut failed: from v0 in cascade, K0s --> " << hfCandCascade.mK0Short() << ", in PDG K0s --> " << RecoDecay::getMassPDG(kK0Short) << ", cut --> " << cuts->get(ptBin, "mK0s");
      return false; // mass of the K0s
    }

    if ((std::abs(hfCandCascade.mLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda")) || (std::abs(hfCandCascade.mAntiLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin, "mLambda"))) {
      LOG(debug) << "mass L cut failed: from v0 in cascade, Lambda --> " << hfCandCascade.mLambda() << ", AntiLambda --> " << hfCandCascade.mAntiLambda() << ", in PDG, Lambda --> " << RecoDecay::getMassPDG(kLambda0) << ", cut --> " << cuts->get(ptBin, "mLambda");
      return false; // mass of the Lambda
    }

    if (std::abs(hfCandCascade.mGamma() - RecoDecay::getMassPDG(kGamma)) < cuts->get(ptBin, "mGamma")) {
      LOG(debug) << "mass gamma cut failed: from v0 in cascade, gamma --> " << hfCandCascade.mGamma() << ", cut --> " << cuts->get(ptBin, "mGamma");
      return false; // mass of the Gamma
    }

    if (hfCandCascade.ptProng0() < cuts->get(ptBin, "ptBach")) {
      LOG(debug) << "bach pt cut failed, from cascade --> " << hfCandCascade.ptProng0() << " , cut --> " << cuts->get(ptBin, "ptBach");
      return false; // pt of the p
    }

    if (hfCandCascade.ptV0Pos() < cuts->get(ptBin, "ptV0Dau")) {
      LOG(debug) << "v0 pos daugh pt cut failed, from cascade --> " << hfCandCascade.ptV0Pos() << ", cut --> " << cuts->get(ptBin, "ptV0Dau");
      return false; // pt of the K0
    }

    if (hfCandCascade.ptV0Neg() < cuts->get(ptBin, "ptV0Dau")) {
      LOG(debug) << "v0 neg daugh pt cut failed, from cascade --> " << hfCandCascade.ptV0Neg() << ", cut --> " << cuts->get(ptBin, "ptV0Dau");
      return false; // pt of the K0
    }

    if (hfCandCascade.ptProng1() < cuts->get(ptBin, "ptV0")) {
      LOG(debug) << "cand pt cut failed, from cascade --> " << hfCandCascade.ptProng1() << ", cut --> " << cuts->get(ptBin, "ptV0");
      return false; // pt of the Lc
    }

    if (std::abs(hfCandCascade.impactParameter0()) > cuts->get(ptBin, "d0Bach")) {
      LOG(debug) << "d0 bach cut failed, in cascade --> " << hfCandCascade.impactParameter0() << ", cut --> " << cuts->get(ptBin, "d0Bach");
      return false; // d0 of the bachelor
    }

    /*
    if ((std::abs(hfCandCascade.dcapostopv()) > d0K0Cut[ptBin]) || (std::abs(hfCandCascade.dcanegtopv()) > d0K0Cut[ptBin])) {
      LOG(debug) << "v0 daugh cut failed, positive v0 daugh --> " << hfCandCascade.dcapostopv() << ", negative v0 daugh --> " << hfCandCascade.dcanegtopv() << " , cut --> " << d0K0Cut[ptBin];
      return false; // d0 of the K0s daughters
    }
    */

    if (std::abs(hfCandCascade.impactParameter1()) > cuts->get(ptBin, "d0V0")) {
      LOG(debug) << "d0 v0 cut failed, in cascade --> " << hfCandCascade.impactParameter1() << ", cut --> " << cuts->get(ptBin, "d0V0");
      return false; // d0 of the v0
    }

    return true;
  }

  template <typename T>
  bool selectionPID(const T& track)
  {
    TrackSelectorPID selectorProton[2] = {TrackSelectorPID(kProton), TrackSelectorPID(kProton)};
    selectorProton[0].setRangeNSigmaTPC(-nSigmaTpcMaxLowP, nSigmaTpcMaxLowP);
    selectorProton[1].setRangeNSigmaTPC(-nSigmaTpcMaxHighP, nSigmaTpcMaxHighP);
    selectorProton[0].setRangeNSigmaTOF(-nSigmaTofMaxLowP, nSigmaTofMaxLowP);
    selectorProton[1].setRangeNSigmaTOF(-nSigmaTofMaxHighP, nSigmaTofMaxHighP);
    selectorProton[0].setProbBayesMin(probBayesMinLowP);
    selectorProton[1].setProbBayesMin(probBayesMinHighP);
    bool requireTof[2] = {requireTofLowP, requireTofHighP};

    int whichSelector;
    if (track.p() < pPidThreshold) {
      whichSelector = 0;
    } else {
      whichSelector = 1;
    }
    int pidProtonTPC = selectorProton[whichSelector].getStatusTrackPIDTPC(track);
    int pidProtonBayes = selectorProton[whichSelector].getStatusTrackBayesProbPID(track);
    int pidProtonTOF = -1;
    if (!track.hasTOF()) {
      if (requireTof[whichSelector]) {
        pidProtonTOF = TrackSelectorPID::Status::PIDRejected;
      } else {
        pidProtonTOF = TrackSelectorPID::Status::PIDAccepted;
      }
    } else {
      pidProtonTOF = selectorProton[whichSelector].getStatusTrackPIDTOF(track);
    }
    return ((pidProtonTPC == TrackSelectorPID::Status::PIDAccepted) && (pidProtonTOF == TrackSelectorPID::Status::PIDAccepted) && (pidProtonBayes == TrackSelectorPID::Status::PIDAccepted));
  }

  void process(aod::HfCandCascade const& candidates, MyBigTracks const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted

    for (const auto& candidate : candidates) {               // looping over cascade candidates
      const auto& bach = candidate.prong0_as<MyBigTracks>(); // bachelor track
#ifdef MY_DEBUG
      auto indexV0DaughPos = candidate.posTrack_as<MyBigTracks>().mcParticleId();
      auto indexV0DaughNeg = candidate.negTrack_as<MyBigTracks>().mcParticleId();
      auto indexBach = bach.mcParticleId();
      bool isLc = isLcK0SpFunc(indexBach, indexV0DaughPos, indexV0DaughNeg, indexProton, indexK0Spos, indexK0Sneg);
      bool isK0SfromLc = isK0SfromLcFunc(indexV0DaughPos, indexV0DaughNeg, indexK0Spos, indexK0Sneg);
#endif
      MY_DEBUG_MSG(isLc, printf("\n"); LOG(info) << "In selector: correct Lc found: proton --> " << indexBach << ", posTrack --> " << indexV0DaughPos << ", negTrack --> " << indexV0DaughNeg);
      // MY_DEBUG_MSG(isLc != 1, printf("\n"); LOG(info) << "In selector: wrong Lc found: proton --> " << indexBach << ", posTrack --> " << indexV0DaughPos << ", negTrack --> " << indexV0DaughNeg);

      statusLc = 0;

      // daughter track validity selection
      LOG(debug) << "daughterSelection(bach) = " << daughterSelection(bach);
      if (!daughterSelection(bach)) {
        MY_DEBUG_MSG(isLc, LOG(info) << "In selector: Lc rejected due to selections on bachelor");
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)
      LOG(debug) << "selectionTopol(candidate) = " << selectionTopol(candidate);
      if (!selectionTopol(candidate)) {
        MY_DEBUG_MSG(isLc, LOG(info) << "In selector: Lc rejected due to topological selections");
        hfSelLcToK0sPCandidate(statusLc);
        continue;
      }

      if (usePidBach) {
        LOG(debug) << "selectionPID(bach) = " << selectionPID(bach);
        if (!selectionPID(bach)) {
          MY_DEBUG_MSG(isLc, LOG(info) << "In selector: Lc rejected due to PID selection on bachelor");
          hfSelLcToK0sPCandidate(statusLc);
          continue;
        }
      }

      statusLc = 1;
      MY_DEBUG_MSG(isLc && pidProton == 1, LOG(info) << "In selector: Lc ACCEPTED");

      hfSelLcToK0sPCandidate(statusLc);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfcg)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorLcToK0sP>(cfcg));
  return workflow;
}
