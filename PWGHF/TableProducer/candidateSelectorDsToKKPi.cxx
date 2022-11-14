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

/// \file candidateSelectorDsToKKPi.cxx
/// \brief Ds± → K± K∓ π± selection task
///
/// \author Stefano Politano <stefano.politano@cern.ch>, Politecnico and INFN Torino

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_ds_to_k_k_pi;

/// Struct for applying Ds to KKpi selection cuts
struct HfCandidateSelectorDsToKKPi {
  Produces<aod::HfSelDsToKKPi> hfSelDsToKKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 1., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC"};
  //  TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_ds_to_k_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_ds_to_k_k_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Ds candidate selection per pT bin"};

  /// Candidate selections
  /// \param candidate is candidate
  /// \param trackKaon1 is the first track with the kaon hypothesis
  /// \param trackKaon2 is the second track with the kaon hypothesis
  /// \param trackPion is the second track with the pion hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selection(const T1& candidate, const T2& trackKaon1, const T2& trackKaon2, const T2& trackPion)
  {
    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }
    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT > ptCandMax) {
      return false;
    }
    // cut on daughter pT
    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    // invariant-mass cut
    if (std::abs(invMassDsToKKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kDS)) > cuts->get(pTBin, "m") && (std::abs(invMassDsToPiKK(candidate) - RecoDecay::getMassPDG(pdg::Code::kDS)) > cuts->get(pTBin, "m"))) {
      return false;
    }
    // decay length cut
    if (candidate.decayLength() < cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    // cos. pointing angle cut
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }
    if (std::abs(candidate.maxNormalisedDeltaIP()) > cuts->get(pTBin, "max normalized deltaIP")) {
      return false;
    }
    return true;
  }

  void process(aod::HfCand3Prong const& candidates, aod::BigTracksPID const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);

    TrackSelectorPID selectorKaon(selectorPion);
    selectorKaon.setPDG(kKPlus);

    // looping over 3-prong candidates
    for (auto& candidate : candidates) {

      // final selection flag:
      auto statusDsToKKPi = 0;
      auto statusDsToPiKK = 0;

      if (!(candidate.hfflag() & 1 << DecayType::DsToKKPi)) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        continue;
      }
      SETBIT(statusDsToKKPi, aod::SelectionStep::RecoSkims);
      SETBIT(statusDsToPiKK, aod::SelectionStep::RecoSkims);

      auto trackPos1 = candidate.prong0_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<aod::BigTracksPID>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)

      // topological selection
      bool topoDsToKKPi = selection(candidate, trackPos1, trackNeg, trackPos2);
      bool topoDsToPiKK = selection(candidate, trackPos2, trackNeg, trackPos1);
      if ((!topoDsToKKPi) && (!topoDsToPiKK)) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        continue;
      }
      if (topoDsToKKPi) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoTopol);
      }
      if (topoDsToPiKK) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoTopol);
      }

      // track-level PID selection
      int pidTrackPos1Pion = selectorPion.getStatusTrackPIDAll(trackPos1);
      int pidTrackPos2Pion = selectorPion.getStatusTrackPIDAll(trackPos2);
      int pidTrackPos1Kaon = selectorKaon.getStatusTrackPIDAll(trackPos1);
      int pidTrackPos2Kaon = selectorKaon.getStatusTrackPIDAll(trackPos2);
      int pidTrackNegKaon = selectorKaon.getStatusTrackPIDAll(trackNeg);
      bool pidDsToKKPi = false;
      bool pidDsToPiKK = false;

      // excluding candidates with negative track rejected as K
      if (pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        continue;
      }
      // checking PID for Ds to KKPi hypothesis
      if (pidTrackPos1Kaon == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoPID); // accept DsKKPi
        pidDsToKKPi = true;
      }
      // checking PID for Ds to PiKK hypothesis
      if (pidTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackPos2Kaon == TrackSelectorPID::Status::PIDAccepted) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoPID); // accept DsPiKK
        pidDsToPiKK = true;
      }
      // both PID hypotheses rejected
      if (!pidDsToKKPi && !pidDsToPiKK) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        continue;
      }

      hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorDsToKKPi>(cfgc)};
}