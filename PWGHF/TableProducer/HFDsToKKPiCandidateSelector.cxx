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

/// \file HFDsToKKPiCandidateSelector.cxx
/// \brief Ds± → K± K∓ π± selection task
///
/// \author Stefano Politano <stefano.politano@cern.ch>, Politecnico and INFN Torino

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::analysis::hf_cuts_ds_tokkpi;

/// Struct for applying Ds to KKpi selection cuts
struct HfDsToKKPiCandidateSelector {
  Produces<aod::HFSelDsToKKPiCandidate> hfSelDsToKKPiCandidate;

  Configurable<double> pTCandMin{"pTCandMin", 1., "Lower bound of candidate pT"};
  Configurable<double> pTCandMax{"pTCandMax", 36., "Upper bound of candidate pT"};
  // TPC
  Configurable<double> pidTPCMinpT{"pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> pidTPCMaxpT{"pidTPCMaxpT", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTPC{"nSigmaTPC", 3., "Nsigma cut on TPC"};
  // Configurable<double> TPCNClsFindablePIDCut{"TPCNClsFindablePIDCut", 50., "Lower bound of TPC findable clusters for good PID"};
  //  TOF
  Configurable<double> pidTOFMinpT{"pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> pidTOFMaxpT{"pidTOFMaxpT", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTOF{"nSigmaTOF", 3., "Nsigma cut on TOF"};
  // topological cuts
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_ds_tokkpi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Ds_to_K_K_Pi_cuts", {hf_cuts_ds_tokkpi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Ds candidate selection per pT bin"};

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
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }
    // check that the candidate pT is within the analysis range
    if (candpT < pTCandMin || candpT > pTCandMax) {
      return false;
    }
    // cut on daughter pT
    if (trackKaon1.pt() < cuts->get(pTBin, "pT K") || trackKaon2.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    // invariant-mass cut
    if (std::abs(InvMassDsKKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kDs)) > cuts->get(pTBin, "m") && (std::abs(InvMassDsPiKK(candidate) - RecoDecay::getMassPDG(pdg::Code::kDs)) > cuts->get(pTBin, "m"))) {
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

  void process(aod::HfCandProng3 const& candidates, aod::BigTracksPID const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(pidTPCMinpT, pidTPCMaxpT);
    selectorPion.setRangeNSigmaTPC(-nSigmaTPC, nSigmaTPC);
    selectorPion.setRangePtTOF(pidTOFMinpT, pidTOFMaxpT);
    selectorPion.setRangeNSigmaTOF(-nSigmaTOF, nSigmaTOF);

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

      auto trackPos1 = candidate.index0_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.index1_as<aod::BigTracksPID>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.index2_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)

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
  return WorkflowSpec{adaptAnalysisTask<HfDsToKKPiCandidateSelector>(cfgc)};
}