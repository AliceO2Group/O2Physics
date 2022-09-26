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

/// \file HFLcK0spCandidateSelector.cxx
/// \brief Lc --> K0s+p selection task.
///
/// \author Chiara Zampolli <Chiara.Zampolli@cern.ch>, CERN
///         Daniel Samitz, <daniel.samitz@cern.ch>, Vienna

/// based on HFD0CandidateSelector.cxx

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "PWGHF/Utils/UtilsDebugLcK0Sp.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_casc;
using namespace o2::analysis::hf_cuts_lc_tok0sp;


#ifdef MY_DEBUG
#define MY_DEBUG_MSG(condition, cmd) \
  if (condition) {                   \
    cmd;                             \
  }
using MyBigTracks = soa::Join<aod::BigTracksPID, aod::McTrackLabels>;
#else
#define MY_DEBUG_MSG(condition, cmd)
using MyBigTracks = aod::BigTracksPID;
#endif

struct HFLcK0sPCandidateSelector {

  Produces<aod::HFSelLcK0sPCandidate> hfSelLcK0sPCandidate;

  Configurable<double> pTCandMin{"pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> pTCandMax{"pTCandMax", 50., "Upper bound of candidate pT"};

  // PID
  Configurable<double> applyPidTPCMinPt{"applyPidTPCMinPt", 4., "Lower bound of track pT to apply TPC PID"};
  Configurable<double> pidTPCMinPt{"pidTPCMinPt", 0., "Lower bound of track pT for TPC PID"};
  Configurable<double> pidTPCMaxPt{"pidTPCMaxPt", 100., "Upper bound of track pT for TPC PID"};
  Configurable<double> pidCombMaxP{"pidCombMaxP", 4., "Upper bound of track p to use TOF + TPC Bayes PID"};
  Configurable<double> nSigmaTPC{"nSigmaTPC", 3., "Nsigma cut on TPC only"};

  // track quality
  Configurable<double> TPCNClsFindablePIDCut{"TPCNClsFindablePIDCut", 50., "Lower bound of TPC findable clusters for good PID"};
  Configurable<bool> requireTPC{"requireTPC", true, "Flag to require a positive Number of found clusters in TPC"};

  //cuts
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_lc_tok0sp::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Lc_to_K0s_p_cuts", {hf_cuts_lc_tok0sp::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Lc candidate selection per pT bin"};

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
    int ptBin = findBin(pTBins,candPt);
    if (ptBin == -1) {
      return false;
    }

    if (candPt < pTCandMin || candPt >= pTCandMax) {
      LOG(debug) << "cand pt (first check) cut failed: from cascade --> " << candPt << ", cut --> " << pTCandMax;
      return false; //check that the candidate pT is within the analysis range
    }

    if (std::abs(hfCandCascade.mK0Short() - RecoDecay::getMassPDG(kK0Short)) > cuts->get(ptBin,"mK0s")) {
      LOG(debug) << "massK0s cut failed: from v0 in cascade, K0s --> " << hfCandCascade.mK0Short() << ", in PDG K0s --> " << RecoDecay::getMassPDG(kK0Short) << ", cut --> " << cuts->get(ptBin,"mK0s");
      return false; // mass of the K0s
    }

    if ((std::abs(hfCandCascade.mLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin,"mLambda")) || (std::abs(hfCandCascade.mAntiLambda() - RecoDecay::getMassPDG(kLambda0)) < cuts->get(ptBin,"mLambda"))) {
      LOG(debug) << "mass L cut failed: from v0 in cascade, Lambda --> " << hfCandCascade.mLambda() << ", AntiLambda --> " << hfCandCascade.mAntiLambda() << ", in PDG, Lambda --> " << RecoDecay::getMassPDG(kLambda0) << ", cut --> " << cuts->get(ptBin,"mLambda");
      return false; // mass of the Lambda
    }

    if (std::abs(InvMassGamma(hfCandCascade) - RecoDecay::getMassPDG(kGamma)) < cuts->get(ptBin,"mGamma")) {
      LOG(debug) << "mass gamma cut failed: from v0 in cascade, gamma --> " << InvMassGamma(hfCandCascade) << ", cut --> " << cuts->get(ptBin,"mGamma");
      return false; // mass of the Gamma
    }

    if (hfCandCascade.ptProng0() < cuts->get(ptBin,"ptBach")) {
      LOG(debug) << "bach pt cut failed, from cascade --> " << hfCandCascade.ptProng0() << " , cut --> " << cuts->get(ptBin,"ptBach");
      return false; // pt of the p
    }

    if (hfCandCascade.ptV0Pos() < cuts->get(ptBin,"ptV0Dau")) {
      LOG(debug) << "v0 pos daugh pt cut failed, from cascade --> " << hfCandCascade.ptV0Pos() << ", cut --> " << cuts->get(ptBin,"ptV0Dau");
      return false; // pt of the K0
    }

    if (hfCandCascade.ptV0Neg() < cuts->get(ptBin,"ptV0Dau")) {
      LOG(debug) << "v0 neg daugh pt cut failed, from cascade --> " << hfCandCascade.ptV0Neg() << ", cut --> " << cuts->get(ptBin,"ptV0Dau");
      return false; // pt of the K0
    }

    if (hfCandCascade.ptProng1() < cuts->get(ptBin,"ptV0")) {
      LOG(debug) << "cand pt cut failed, from cascade --> " << hfCandCascade.ptProng1() << ", cut --> " << cuts->get(ptBin,"ptV0");
      return false; // pt of the Lc
    }

    if (std::abs(hfCandCascade.impactParameter0()) > cuts->get(ptBin,"d0Bach")) {
      LOG(debug) << "d0 bach cut failed, in cascade --> " << hfCandCascade.impactParameter0() << ", cut --> " << cuts->get(ptBin,"d0Bach");
      return false; // d0 of the bachelor
    }

    /*
    if ((std::abs(hfCandCascade.dcapostopv()) > d0K0Cut[ptBin]) || (std::abs(hfCandCascade.dcanegtopv()) > d0K0Cut[ptBin])) {
      LOG(debug) << "v0 daugh cut failed, positive v0 daugh --> " << hfCandCascade.dcapostopv() << ", negative v0 daugh --> " << hfCandCascade.dcanegtopv() << " , cut --> " << d0K0Cut[ptBin];
      return false; // d0 of the K0s daughters
    }
    */

    if (std::abs(hfCandCascade.impactParameter1()) > cuts->get(ptBin,"d0V0")) {
      LOG(debug) << "d0 v0 cut failed, in cascade --> " << hfCandCascade.impactParameter1() << ", cut --> " << cuts->get(ptBin,"d0V0");
      return false; // d0 of the v0
    }

    return true;
  }

  /// Check if track is ok for TPC PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TPC PID
  template <typename T>
  bool validTPCPID(const T& track)
  {
    if (track.pt() < pidTPCMinPt || track.pt() >= pidTPCMaxPt) {
      LOG(debug) << "Bachelor pt is " << track.pt() << ", we trust TPC PID in [" << pidTPCMinPt << ", " << pidTPCMaxPt << "]";
      return false;
    }
    return true;
  }

  /// Check if we will use TPC PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TPC PID
  template <typename T>
  bool applyTPCPID(const T& track)
  {
    if (track.pt() < applyPidTPCMinPt) {
      LOG(debug) << "Bachelor pt is " << track.pt() << ", we apply TPC PID from " << applyPidTPCMinPt;
      return false;
    }
    LOG(debug) << "Bachelor pt is " << track.pt() << ", we apply TPC PID from " << applyPidTPCMinPt;
    return true;
  }

  /// Check if track is ok for TOF PID
  /// \param track is the track
  /// \note function to be expanded
  /// \return true if track is ok for TOF PID
  template <typename T>
  bool validCombPID(const T& track)
  {
    if (track.pt() > pidCombMaxP) { // is the pt sign used for the charge? If it is always positive, we should remove the abs
      return false;
    }
    return true;
  }

  /// Check if track is compatible with given TPC Nsigma cut for a given flavour hypothesis
  /// \param track is the track
  /// \param nPDG is the flavour hypothesis PDG number
  /// \param nSigmaCut is the nsigma threshold to test against
  /// \note nPDG=211 pion  nPDG=321 kaon
  /// \return true if track satisfies TPC PID hypothesis for given Nsigma cut
  template <typename T>
  bool selectionPIDTPC(const T& track, double nSigmaCut)
  {
    double nSigma = 100.0; //arbitarily large value
    nSigma = track.tpcNSigmaPr();
    LOG(debug) << "nSigma for bachelor = " << nSigma << ", cut at " << nSigmaCut;
    return std::abs(nSigma) < nSigmaCut;
  }

  /*
  /// Check if track is compatible with given TOF NSigma cut for a given flavour hypothesis
  /// \param track is the track
  /// \param nPDG is the flavour hypothesis PDG number
  /// \param nSigmaCut is the nSigma threshold to test against
  /// \note nPDG=211 pion  nPDG=321 kaon
  /// \return true if track satisfies TOF PID hypothesis for given NSigma cut
  template <typename T>
  bool selectionPIDTOF(const T& track, int nPDG, double nSigmaCut)
  {
    double nSigma = 100.0; //arbitarily large value
    nPDG = TMath::Abs(nPDG);
    if (nPDG == 111) {
      nSigma = track.tofNSigmaPi();
    } else if (nPDG == 321) {
      nSigma = track.tofNSigmaKa();
    } else {
      return false;
    }
    return nSigma < nSigmaCut;
  }
  */

  /// PID selection on daughter track
  /// \param track is the daughter track
  /// \param nPDG is the PDG code of the flavour hypothesis
  /// \note nPDG=211 pion  nPDG=321 kaon
  /// \return 1 if successful PID match, 0 if successful PID rejection, -1 if no PID info
  template <typename T>
  int selectionPID(const T& track)
  {
    int statusTPC = -1;
    //    int statusTOF = -1;

    if (!applyTPCPID(track)) {
      // we do not apply TPC PID in this range
      return 1;
    }

    if (validTPCPID(track)) {
      LOG(debug) << "We check the TPC PID now";
      if (!selectionPIDTPC(track, nSigmaTPC)) {
        statusTPC = 0;
        /*
        if (!selectionPIDTPC(track, nPDG, nSigmaTPCCombined)) {
          statusTPC = 0; //rejected by PID
        } else {
          statusTPC = 1; //potential to be acceepted if combined with TOF
        }
      } else {
        statusTPC = 2; //positive PID
      }
	*/
      } else {
        statusTPC = 1;
      }
    }

    return statusTPC;
    /*
    if (validTOFPID(track)) {
      if (!selectionPIDTOF(track, nPDG, nSigmaTOF)) {
        if (!selectionPIDTOF(track, nPDG, nSigmaTOFCombined)) {
          statusTOF = 0; //rejected by PID
        } else {
          statusTOF = 1; //potential to be acceepted if combined with TOF
        }
      } else {
        statusTOF = 2; //positive PID
      }
    } else {
      statusTOF = -1; //no PID info
    }

    if (statusTPC == 2 || statusTOF == 2) {
      return 1; //what if we have 2 && 0 ?
    } else if (statusTPC == 1 && statusTOF == 1) {
      return 1;
    } else if (statusTPC == 0 || statusTOF == 0) {
      return 0;
    } else {
      return -1;
    }
      */
  }

  void process(aod::HfCandCascade const& candidates, MyBigTracks const& tracks)
  {
    int statusLc = 0; // final selection flag : 0-rejected  1-accepted
    //bool topolLc = 0;
    int pidProton = -1;
    //int pidLc = -1;

    for (auto& candidate : candidates) { //looping over cascade candidates
      const auto& bach = candidate.index0_as<MyBigTracks>(); //bachelor track
#ifdef MY_DEBUG
      auto indexV0DaughPos = candidate.posTrack_as<MyBigTracks>().mcParticleId();
      auto indexV0DaughNeg = candidate.negTrack_as<MyBigTracks>().mcParticleId();
      auto indexBach = bach.mcParticleId();
      bool isLc = isLcK0SpFunc(indexBach, indexV0DaughPos, indexV0DaughNeg, indexProton, indexK0Spos, indexK0Sneg);
      bool isK0SfromLc = isK0SfromLcFunc(indexV0DaughPos, indexV0DaughNeg, indexK0Spos, indexK0Sneg);
#endif
      MY_DEBUG_MSG(isLc, printf("\n"); LOG(info) << "In selector: correct Lc found: proton --> " << indexBach << ", posTrack --> " << indexV0DaughPos << ", negTrack --> " << indexV0DaughNeg);
      //MY_DEBUG_MSG(isLc != 1, printf("\n"); LOG(info) << "In selector: wrong Lc found: proton --> " << indexBach << ", posTrack --> " << indexV0DaughPos << ", negTrack --> " << indexV0DaughNeg);

      statusLc = 0;
      /* // not needed for the Lc
      if (!(candidate.hfflag() & 1 << D0ToPiK)) {
        hfSelD0Candidate(statusLc);
        continue;
      }
      */

      //topolLc = true;
      pidProton = -1;

      // daughter track validity selection
      LOG(debug) << "daughterSelection(bach) = " << daughterSelection(bach);
      if (!daughterSelection(bach)) {
        MY_DEBUG_MSG(isLc, LOG(info) << "In selector: Lc rejected due to selections on bachelor");
        hfSelLcK0sPCandidate(statusLc);
        continue;
      }

      //implement filter bit 4 cut - should be done before this task at the track selection level
      //need to add special cuts (additional cuts on decay length and d0 norm)
      LOG(debug) << "selectionTopol(candidate) = " << selectionTopol(candidate);
      if (!selectionTopol(candidate)) {
        MY_DEBUG_MSG(isLc, LOG(info) << "In selector: Lc rejected due to topological selections");
        hfSelLcK0sPCandidate(statusLc);
        continue;
      }

      pidProton = selectionPID(bach);

      LOG(debug) << "pidProton = " << pidProton;

      if (pidProton == 1) {
        statusLc = 1;
      }

      MY_DEBUG_MSG(isLc && pidProton != 1, LOG(info) << "In selector: Lc rejected due to PID selections on bachelor");
      MY_DEBUG_MSG(isLc && pidProton == 1, LOG(info) << "In selector: Lc ACCEPTED");

      hfSelLcK0sPCandidate(statusLc);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfcg)
{
  return WorkflowSpec{
    adaptAnalysisTask<HFLcK0sPCandidateSelector>(cfcg, TaskName{"hf-lc-tok0sp-candidate-selector"})};
}
