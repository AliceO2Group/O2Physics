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

#include "Common/Core/TrackSelectorPID.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

/// Struct for applying D0 selection cuts
struct HfCandidateSelectorD0 {
  Produces<aod::HfSelD0> hfSelD0Candidate;

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
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_d0_to_pi_k::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "D0 candidate selection per pT bin"};

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

    // decay exponentail law, with tau = beta*gamma*ctau
    // decay length > ctau retains (1-1/e)
    if (std::abs(candidate.impactParameterNormalised0()) < 0.5 || std::abs(candidate.impactParameterNormalised1()) < 0.5) {
      return false;
    }
    double decayLengthCut = std::min((candidate.p() * 0.0066) + 0.01, cuts->get(pTBin, "minimum decay length"));
    if (candidate.decayLength() * candidate.decayLength() < decayLengthCut * decayLengthCut) {
      return false;
    }
    if (candidate.decayLength() > cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXY() > cuts->get(pTBin, "decay length XY")) {
      return false;
    }
    if (candidate.decayLengthNormalised() * candidate.decayLengthNormalised() < 1.0) {
      // return false; // add back when getter fixed
    }
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \note trackPion = positive and trackKaon = negative for D0 selection and inverse for D0bar
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& trackKaon)
  {
    auto candpT = candidate.pt();
    auto pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // invariant-mass cut
    if (trackPion.sign() > 0) {
      if (std::abs(invMassD0ToPiK(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(invMassD0barToKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "m")) {
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
      if (std::abs(cosThetaStarD0(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    } else {
      if (std::abs(cosThetaStarD0bar(candidate)) > cuts->get(pTBin, "cos theta*")) {
        return false;
      }
    }

    // in case only sideband candidates have to be stored, additional invariant-mass cut
    if (keepOnlySidebandCandidates) {
      if (trackPion.sign() > 0) {
        if (std::abs(invMassD0ToPiK(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) < distanceFromD0MassForSidebands) {
          return false;
        }
      } else {
        if (std::abs(invMassD0barToKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) < distanceFromD0MassForSidebands) {
          return false;
        }
      }
    }

    return true;
  }

  void process(aod::HfCand2Prong const& candidates, aod::BigTracksPIDExtended const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);

    TrackSelectorPID selectorKaon(selectorPion);
    selectorKaon.setPDG(kKPlus);

    // looping over 2-prong candidates
    for (auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      int statusD0 = 0;
      int statusD0bar = 0;
      int statusHFFlag = 0;
      int statusTopol = 0;
      int statusCand = 0;
      int statusPID = 0;

      if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusHFFlag = 1;

      auto trackPos = candidate.prong0_as<aod::BigTracksPIDExtended>(); // positive daughter
      auto trackNeg = candidate.prong1_as<aod::BigTracksPIDExtended>(); // negative daughter

      /*
      if (!daughterSelection(trackPos) || !daughterSelection(trackNeg)) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }
      */

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusTopol = 1;

      // implement filter bit 4 cut - should be done before this task at the track selection level
      // need to add special cuts (additional cuts on decay length and d0 norm)

      // conjugate-dependent topological selection for D0
      bool topolD0 = selectionTopolConjugate(candidate, trackPos, trackNeg);
      // conjugate-dependent topological selection for D0bar
      bool topolD0bar = selectionTopolConjugate(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        continue;
      }
      statusCand = 1;

      // track-level PID selection
      int pidTrackPosKaon = -1;
      int pidTrackPosPion = -1;
      int pidTrackNegKaon = -1;
      int pidTrackNegPion = -1;

      if (usePidTpcAndTof) {
        pidTrackPosKaon = selectorKaon.getStatusTrackPIDTpcAndTof(trackPos);
        pidTrackPosPion = selectorPion.getStatusTrackPIDTpcAndTof(trackPos);
        pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcAndTof(trackNeg);
        pidTrackNegPion = selectorPion.getStatusTrackPIDTpcAndTof(trackNeg);
      } else {
        pidTrackPosKaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackPos);
        pidTrackPosPion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos);
        pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackNeg);
        pidTrackNegPion = selectorPion.getStatusTrackPIDTpcOrTof(trackNeg);
      }

      int pidD0 = -1;
      int pidD0bar = -1;

      if (pidTrackPosPion == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted) {
        pidD0 = 1; // accept D0
      } else if (pidTrackPosPion == TrackSelectorPID::Status::PIDRejected ||
                 pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected) {
        pidD0 = 0; // exclude D0
      }

      if (pidTrackNegPion == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackPosKaon == TrackSelectorPID::Status::PIDAccepted) {
        pidD0bar = 1; // accept D0bar
      } else if (pidTrackNegPion == TrackSelectorPID::Status::PIDRejected ||
                 pidTrackPosKaon == TrackSelectorPID::Status::PIDRejected) {
        pidD0bar = 0; // exclude D0bar
      }

      if (pidD0 == 0 && pidD0bar == 0) {
        hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
        continue;
      }

      if ((pidD0 == -1 || pidD0 == 1) && topolD0) {
        statusD0 = 1; // identified as D0
      }
      if ((pidD0bar == -1 || pidD0bar == 1) && topolD0bar) {
        statusD0bar = 1; // identified as D0bar
      }
      statusPID = 1;
      hfSelD0Candidate(statusD0, statusD0bar, statusHFFlag, statusTopol, statusCand, statusPID);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorD0>(cfgc)};
}
