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

#include "Common/Core/TrackSelectorPID.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_lc_to_p_k_pi;

/// Struct for applying Lc selection cuts
struct HfCandidateSelectorLc {
  Produces<aod::HfSelLc> hfSelLcCandidate;

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
  Configurable<bool> usePidBayes{"usePidBayes", true, "Bool to use or not the PID based on Bayesian probability cut at filtering level"};
  Configurable<double> ptPidBayesMin{"ptPidBayesMin", 0., "Lower bound of track pT for Bayesian PID"};
  Configurable<double> ptPidBayesMax{"ptPidBayesMax", 100, "Upper bound of track pT for Bayesian PID"};
  // Combined PID options
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Bool to decide how to combine TPC and TOF PID: true = both (if present, only one otherwise); false = one is enough"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_lc_to_p_k_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Lc candidate selection per pT bin"};

  using TrksPID = soa::Join<aod::BigTracksPID, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidBayes>;

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
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackProton is the track with the proton hypothesis
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2>
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

    if (trackProton.globalIndex() == candidate.prong0Id()) {
      if (std::abs(invMassLcToPKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(invMassLcToPiKP(candidate) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    return true;
  }

  void process(aod::HfCand3Prong const& candidates, TrksPID const&)
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

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }

      // conjugate-dependent topological selection for Lc

      bool topolLcToPKPi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcToPiKP = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolLcToPKPi && !topolLcToPiKP) {
        hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
        continue;
      }

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
        int pidTrackPos1Proton = 999;
        int pidTrackPos2Proton = 999;
        int pidTrackPos1Pion = 999;
        int pidTrackPos2Pion = 999;
        int pidTrackNegKaon = 999;
        if (usePidTpcAndTof) {
          pidTrackPos1Proton = selectorProton.getStatusTrackPIDTpcAndTof(trackPos1);
          pidTrackPos2Proton = selectorProton.getStatusTrackPIDTpcAndTof(trackPos2);
          pidTrackPos1Pion = selectorPion.getStatusTrackPIDTpcAndTof(trackPos1);
          pidTrackPos2Pion = selectorPion.getStatusTrackPIDTpcAndTof(trackPos2);
          pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcAndTof(trackNeg);
        } else {
          pidTrackPos1Proton = selectorProton.getStatusTrackPIDTpcOrTof(trackPos1);
          pidTrackPos2Proton = selectorProton.getStatusTrackPIDTpcOrTof(trackPos2);
          pidTrackPos1Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos1);
          pidTrackPos2Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos2);
          pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackNeg);
        }

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

      if ((pidLcToPKPi == -1 || pidLcToPKPi == 1) && (pidBayesLcToPKPi == -1 || pidBayesLcToPKPi == 1) && topolLcToPKPi) {
        statusLcToPKPi = 1; // identified as LcToPKPi
      }
      if ((pidLcToPiKP == -1 || pidLcToPiKP == 1) && (pidBayesLcToPiKP == -1 || pidBayesLcToPiKP == 1) && topolLcToPiKP) {
        statusLcToPiKP = 1; // identified as LcToPiKP
      }

      hfSelLcCandidate(statusLcToPKPi, statusLcToPiKP);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorLc>(cfgc)};
}
