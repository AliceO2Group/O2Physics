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

/// \file candidateSelectorXicToPKPi.cxx
/// \brief Ξc± → p± K∓ π± selection task
/// \note Inspired from candidateSelectorLc.cxx
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN PADOVA
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Common/Core/TrackSelectorPID.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_xic_to_p_k_pi;

/// Struct for applying Xic selection cuts
struct HfCandidateSelectorXicToPKPi {
  Produces<aod::HfSelXicToPKPi> hfSelXicToPKPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> usePid{"usePid", true, "Bool to use or not the PID at filtering level"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 1., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.5, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 4., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<double> decayLengthXYNormalisedMin{"decayLengthXYNormalisedMin", 3., "Min. normalised decay length XY"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xic_to_p_k_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "Xic candidate selection per pT bin"};

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

    // candidate chi2PC
    if (candidate.chi2PCA() > cuts->get(pTBin, "chi2PCA")) {
      return false;
    }

    // candidate decay length
    if (candidate.decayLength() <= cuts->get(pTBin, "decay length")) {
      return false;
    }

    // candidate decay length XY
    if (candidate.decayLengthXY() <= cuts->get(pTBin, "decLengthXY")) {
      return false;
    }

    // candidate normalized decay length XY
    if (candidate.decayLengthXYNormalised() < cuts->get(pTBin, "normDecLXY")) {
      return false;
    }

    // candidate normalised decay length (Inspired from Lc selector)
    if (candidate.decayLengthXYNormalised() < decayLengthXYNormalisedMin) {
      return false;
    }

    // candidate ct
    if (ctXic(candidate) > cuts->get(pTBin, "ct")) {
      return false;
    }

    // candidate impact parameter XY
    if (std::abs(candidate.impactParameterXY()) > cuts->get(pTBin, "impParXY")) {
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
      if (std::abs(invMassXicToPKPi(candidate) - RecoDecay::getMassPDG(pdg::Code::kXiCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(invMassXicToPiKP(candidate) - RecoDecay::getMassPDG(pdg::Code::kXiCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    return true;
  }

  void process(aod::HfCand3Prong const& candidates, aod::BigTracksPID const&)
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

    TrackSelectorPID selectorProton(selectorPion);
    selectorProton.setPDG(kProton);

    // looping over 3-prong candidates
    for (auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      auto statusXicToPKPi = 0;
      auto statusXicToPiKP = 0;

      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        continue;
      }

      auto trackPos1 = candidate.prong0_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.prong1_as<aod::BigTracksPID>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.prong2_as<aod::BigTracksPID>(); // positive daughter (negative for the antiparticles)

      /*
      // daughter track validity selection
      if (!daughterSelection(trackPos1) || !daughterSelection(trackNeg) || !daughterSelection(trackPos2)) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        continue;
      }
      */

      // implement filter bit 4 cut - should be done before this task at the track selection level

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        continue;
      }

      // conjugate-dependent topplogical selection for Xic

      bool topolXicToPKPi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolXicToPiKP = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolXicToPKPi && !topolXicToPiKP) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        continue;
      }

      auto pidXicToPKPi = -1;
      auto pidXicToPiKP = -1;

      if (!usePid) {
        // PID non applied
        pidXicToPKPi = 1;
        pidXicToPiKP = 1;
      } else {
        // track-level PID selection
        auto pidTrackPos1Proton = selectorProton.getStatusTrackPIDTpcOrTof(trackPos1);
        auto pidTrackPos2Proton = selectorProton.getStatusTrackPIDTpcOrTof(trackPos2);
        auto pidTrackPos1Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos1);
        auto pidTrackPos2Pion = selectorPion.getStatusTrackPIDTpcOrTof(trackPos2);
        auto pidTrackNegKaon = selectorKaon.getStatusTrackPIDTpcOrTof(trackNeg);

        if (pidTrackPos1Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidXicToPKPi = 1; // accept XicpKpi
        } else if (pidTrackPos1Proton == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos2Pion == TrackSelectorPID::Status::PIDRejected) {
          pidXicToPKPi = 0; // exclude XicpKpi
        }
        if (pidTrackPos2Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidXicToPiKP = 1; // accept XicpiKp
        } else if (pidTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos2Proton == TrackSelectorPID::Status::PIDRejected) {
          pidXicToPiKP = 0; // exclude XicpiKp
        }
      }

      if (pidXicToPKPi == 0 && pidXicToPiKP == 0) {
        hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
        continue;
      }

      if ((pidXicToPKPi == -1 || pidXicToPKPi == 1) && topolXicToPKPi) {
        statusXicToPKPi = 1; // identified as Xic->pKpi
      }
      if ((pidXicToPiKP == -1 || pidXicToPiKP == 1) && topolXicToPiKP) {
        statusXicToPiKP = 1; // identified as Xic->piKp
      }

      hfSelXicToPKPiCandidate(statusXicToPKPi, statusXicToPiKP);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorXicToPKPi>(cfgc)};
}
