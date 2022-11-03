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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::analysis::hf_cuts_lc_topkpi;

/// Struct for applying Lc selection cuts
struct HfCandidateSelectorLc {
  Produces<aod::HFSelLcCandidate> hfSelLcCandidate;

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> d_FilterPID{"d_FilterPID", true, "Bool to use or not the PID based on nSigma cut at filtering level"};
  // TPC
  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 0.1, "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 1., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 3., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTPCCombined{"d_nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};
  //Configurable<double> d_TPCNClsFindablePIDCut{"d_TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};
  // TOF
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.5, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 2.5, "Upper bound of track pT for TOF PID"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaTOFCombined{"d_nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};
  // Bayesian
  Configurable<bool> d_BayesPID{"d_BayesPID", true, "Bool to use or not the PID based on Bayesian probability cut at filtering level"};
  Configurable<double> d_pidBayesMinpT{"d_pidBayesMinpT", 0., "Lower bound of track pT for Bayesian PID"};
  Configurable<double> d_pidBayesMaxpT{"d_pidBayesMaxpT", 100, "Upper bound of track pT for Bayesian PID"};
  // topological cuts
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_lc_topkpi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Lc_to_p_K_pi_cuts", {hf_cuts_lc_topkpi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Lc candidate selection per pT bin"};

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

    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    //candidate chi2PCA
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
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // cut on daughter pT
    if (trackProton.pt() < cuts->get(pTBin, "pT p") || trackKaon.pt() < cuts->get(pTBin, "pT K") || trackPion.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    if (trackProton.globalIndex() == candidate.index0Id()) {
      if (std::abs(InvMassLcpKpi(candidate) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    } else {
      if (std::abs(InvMassLcpiKp(candidate) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "m")) {
        return false;
      }
    }

    return true;
  }

  using TrksPID = soa::Join<aod::BigTracksPID, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidBayes>;

  void process(aod::HfCandProng3 const& candidates, TrksPID const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPion.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPion.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPion.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPion.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPion.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);
    selectorPion.setRangePtBayes(d_pidBayesMinpT, d_pidBayesMaxpT);

    TrackSelectorPID selectorKaon(selectorPion);
    selectorKaon.setPDG(kKPlus);

    TrackSelectorPID selectorProton(selectorPion);
    selectorProton.setPDG(kProton);

    // looping over 3-prong candidates
    for (auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      auto statusLcpKpi = 0;
      auto statusLcpiKp = 0;

      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        hfSelLcCandidate(statusLcpKpi, statusLcpiKp);
        continue;
      }

      auto trackPos1 = candidate.index0_as<TrksPID>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.index1_as<TrksPID>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.index2_as<TrksPID>(); // positive daughter (negative for the antiparticles)

      /*
      // daughter track validity selection
      if (!daughterSelection(trackPos1) || !daughterSelection(trackNeg) || !daughterSelection(trackPos2)) {
        hfSelLcCandidate(statusLcpKpi, statusLcpiKp);
        continue;
      }
      */

      // implement filter bit 4 cut - should be done before this task at the track selection level

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelLcCandidate(statusLcpKpi, statusLcpiKp);
        continue;
      }

      // conjugate-dependent topological selection for Lc

      bool topolLcpKpi = selectionTopolConjugate(candidate, trackPos1, trackNeg, trackPos2);
      bool topolLcpiKp = selectionTopolConjugate(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolLcpKpi && !topolLcpiKp) {
        hfSelLcCandidate(statusLcpKpi, statusLcpiKp);
        continue;
      }

      auto pidLcpKpi = -1;
      auto pidLcpiKp = -1;
      auto pidBayesLcpKpi = -1;
      auto pidBayesLcpiKp = -1;

      if (!d_FilterPID) {
        // PID non applied
        pidLcpKpi = 1;
        pidLcpiKp = 1;
      } else {
        // track-level PID selection
        int pidTrackPos1Proton = selectorProton.getStatusTrackPIDAll(trackPos1);
        int pidTrackPos2Proton = selectorProton.getStatusTrackPIDAll(trackPos2);
        int pidTrackPos1Pion = selectorPion.getStatusTrackPIDAll(trackPos1);
        int pidTrackPos2Pion = selectorPion.getStatusTrackPIDAll(trackPos2);
        int pidTrackNegKaon = selectorKaon.getStatusTrackPIDAll(trackNeg);

        if (pidTrackPos1Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidLcpKpi = 1; // accept LcpKpi
        } else if (pidTrackPos1Proton == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos2Pion == TrackSelectorPID::Status::PIDRejected) {
          pidLcpKpi = 0; // exclude LcpKpi
        }
        if (pidTrackPos2Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidLcpiKp = 1; // accept LcpiKp
        } else if (pidTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidTrackPos2Proton == TrackSelectorPID::Status::PIDRejected) {
          pidLcpiKp = 0; // exclude LcpiKp
        }
      }

      if (!d_BayesPID) {
        // PID non applied
        pidBayesLcpKpi = 1;
        pidBayesLcpiKp = 1;
      } else {
        int pidBayesTrackPos1Proton = selectorProton.getStatusTrackBayesPID(trackPos1);
        int pidBayesTrackPos2Proton = selectorProton.getStatusTrackBayesPID(trackPos2);
        int pidBayesTrackPos1Pion = selectorPion.getStatusTrackBayesPID(trackPos1);
        int pidBayesTrackPos2Pion = selectorPion.getStatusTrackBayesPID(trackPos2);
        int pidBayesTrackNegKaon = selectorKaon.getStatusTrackBayesPID(trackNeg);

        if (pidBayesTrackPos1Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackPos2Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidBayesLcpKpi = 1; // accept LcpKpi
        } else if (pidBayesTrackPos1Proton == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackPos2Pion == TrackSelectorPID::Status::PIDRejected) {
          pidBayesLcpKpi = 0; // exclude LcpKpi
        }
        if (pidBayesTrackPos2Proton == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDAccepted &&
            pidBayesTrackPos1Pion == TrackSelectorPID::Status::PIDAccepted) {
          pidBayesLcpiKp = 1; // accept LcpiKp
        } else if (pidBayesTrackPos1Pion == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                   pidBayesTrackPos2Proton == TrackSelectorPID::Status::PIDRejected) {
          pidBayesLcpiKp = 0; // exclude LcpiKp
        }
      }

      if (pidLcpKpi == 0 && pidLcpiKp == 0) {
        hfSelLcCandidate(statusLcpKpi, statusLcpiKp);
        continue;
      }

      if (pidBayesLcpKpi == 0 && pidBayesLcpiKp == 0) {
        hfSelLcCandidate(statusLcpKpi, statusLcpiKp);
        continue;
      }

      if ((pidLcpKpi == -1 || pidLcpKpi == 1) && (pidBayesLcpKpi == -1 || pidBayesLcpKpi == 1) && topolLcpKpi) {
        statusLcpKpi = 1; // identified as LcpKpi
      }
      if ((pidLcpiKp == -1 || pidLcpiKp == 1) && (pidBayesLcpiKp == -1 || pidBayesLcpiKp == 1) && topolLcpiKp) {
        statusLcpiKp = 1; // identified as LcpiKp
      }

      hfSelLcCandidate(statusLcpKpi, statusLcpiKp);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorLc>(cfgc)};
}
