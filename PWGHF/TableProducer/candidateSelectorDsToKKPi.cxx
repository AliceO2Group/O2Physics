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
/// \author Fabio Catalano <fabio.catalano@cern.ch>, Universita and INFN Torino
/// \author Stefano Politano <stefano.politano@cern.ch>, Politecnico and INFN Torino

#include "Common/Core/TrackSelectorPID.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

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

  /// Candidate selections independent from the daugther-mass hypothesis
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T1>
  bool selection(const T1& candidate)
  {
    auto candpT = candidate.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    if (candpT < ptCandMin || candpT > ptCandMax) { // check that the candidate pT is within the analysis range
      return false;
    }
    if (candidate.decayLength() < cuts->get(pTBin, "decay length")) {
      return false;
    }
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normalized decay length XY")) {
      return false;
    }
    if (candidate.cpa() < cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }
    if (candidate.cpaXY() < cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }
    if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "impact parameter XY")) {
      return false;
    }
    return true;
  }

  /// Candidate selections for the KKPi daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param trackKaon1 is the first track with the kaon hypothesis
  /// \param trackKaon2 is the second track with the kaon hypothesis
  /// \param trackPion is the track with the pion hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionKKPi(const T1& candidate, const T2& trackKaon1, const T2& trackKaon2, const T2& trackPion)
  {
    int pTBin = findBin(binsPt, candidate.pt());
    if (pTBin == -1) {
      return false;
    }

    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (std::abs(invMassDsToKKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kDS)) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (deltaMassPhiDsToKKPi(candidate) > cuts->get(pTBin, "deltaM Phi")) {
      return false;
    }
    if (std::abs(cos3PiKDsToKKPi(candidate)) < cuts->get(pTBin, "cos^3 theta_PiK")) {
      return false;
    }
    return true;
  }

  /// Candidate selections for the PiKK daugther-mass hypothesis
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon1 is the first track with the kaon hypothesis
  /// \param trackKaon2 is the second track with the kaon hypothesis
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2>
  bool selectionPiKK(const T1& candidate, const T2& trackPion, const T2& trackKaon1, const T2& trackKaon2)
  {
    int pTBin = findBin(binsPt, candidate.pt());
    if (pTBin == -1) {
      return false;
    }

    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (std::abs(invMassDsToPiKK(candidate) - RecoDecay::getMassPDG(pdg::Code::kDS)) > cuts->get(pTBin, "deltaM")) {
      return false;
    }
    if (deltaMassPhiDsToPiKK(candidate) > cuts->get(pTBin, "deltaM Phi")) {
      return false;
    }
    if (std::abs(cos3PiKDsToPiKK(candidate)) < cuts->get(pTBin, "cos^3 theta_PiK")) {
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

      // topological selections
      if (!selection(candidate)) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        continue;
      }

      bool topolDsToKKPi = selectionKKPi(candidate, trackPos1, trackNeg, trackPos2);
      bool topolDsToPiKK = selectionPiKK(candidate, trackPos1, trackNeg, trackPos2);
      if (!topolDsToKKPi && !topolDsToPiKK) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        continue;
      }
      if (topolDsToKKPi) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoTopol);
      }
      if (topolDsToPiKK) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoTopol);
      }

      // track-level PID selection
      int pidTrackPos1Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos1);
      int pidTrackPos1Kaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackPos1);
      int pidTrackPos2Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos2);
      int pidTrackPos2Kaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackPos2);
      int pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackNeg);

      bool pidDsToKKPi = !(pidTrackPos1Kaon == TrackSelectorPID::Status::PIDRejected ||
                           pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                           pidTrackPos2Pion == TrackSelectorPID::Status::PIDRejected);

      bool pidDsToPiKK = !(pidTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
                           pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                           pidTrackPos2Kaon == TrackSelectorPID::Status::PIDRejected);

      if (!pidDsToKKPi && !pidDsToPiKK) {
        hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
        continue;
      }
      if (pidDsToKKPi) {
        SETBIT(statusDsToKKPi, aod::SelectionStep::RecoPID);
      }
      if (pidDsToPiKK) {
        SETBIT(statusDsToPiKK, aod::SelectionStep::RecoPID);
      }

      hfSelDsToKKPiCandidate(statusDsToKKPi, statusDsToPiKK);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorDsToKKPi>(cfgc)};
}
