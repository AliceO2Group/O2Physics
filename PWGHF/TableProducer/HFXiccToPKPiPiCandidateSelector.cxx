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

/// \file HFXiccToPKPiPiCandidateSelector.cxx
/// \brief Ξcc± → p± K∓ π± π±selection task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::framework;
//using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_xicc;
using namespace o2::analysis::hf_cuts_xicc_topkpipi;

/// Struct for applying Xicc selection cuts
struct HfXiccToPKPiPiCandidateSelector {
  Produces<aod::HFSelXiccToPKPiPiCandidate> hfSelXiccToPKPiPiCandidate;

  Configurable<double> d_pTCandMin{"d_pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> d_pTCandMax{"d_pTCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> d_FilterPID{"d_FilterPID", true, "Bool to use or not the PID at filtering level"};
  // TPC
  Configurable<double> d_pidTPCMinpT{"d_pidTPCMinpT", 9999., "Lower bound of track pT for TPC PID"};
  Configurable<double> d_pidTPCMaxpT{"d_pidTPCMaxpT", 99999., "Upper bound of track pT for TPC PID"};
  Configurable<double> d_nSigmaTPC{"d_nSigmaTPC", 9999., "Nsigma cut on TPC only"};
  Configurable<double> d_nSigmaTPCCombined{"d_nSigmaTPCCombined", 9999., "Nsigma cut on TPC combined with TOF"};
  // TOF
  Configurable<double> d_pidTOFMinpT{"d_pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> d_pidTOFMaxpT{"d_pidTOFMaxpT", 4., "Upper bound of track pT for TOF PID"};
  Configurable<double> d_nSigmaTOF{"d_nSigmaTOF", 3., "Nsigma cut on TOF only"};
  Configurable<double> d_nSigmaTOFCombined{"d_nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_xicc_topkpipi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Xicc_to_p_K_pi_pi_cuts", {hf_cuts_xicc_topkpipi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Xicc candidate selection per pT bin"};

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandXicc, const T2& hfCandXic, const T3& trackPi)
  {
    auto candpT = hfCandXicc.pt();
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < d_pTCandMin || candpT >= d_pTCandMax) {
      return false;
    }

    // check candidate mass is within a defined mass window
    if (std::abs(InvMassXiccToXicPi(hfCandXicc) - RecoDecay::getMassPDG(pdg::Code::kXiCCPlusPlus)) > cuts->get(pTBin, "m")) {
      return false;
    }

    // impact parameter product
    if (hfCandXicc.impactParameterProduct() > cuts->get(pTBin, "d0d0")) {
      return false;
    }

    // cosine of pointing angle
    if (hfCandXicc.cpa() <= cuts->get(pTBin, "cos pointing angle")) {
      return false;
    }

    // cosine of pointing angle XY
    if (hfCandXicc.cpaXY() <= cuts->get(pTBin, "cos pointing angle XY")) {
      return false;
    }

    // candidate maximum decay length
    if (hfCandXicc.decayLength() > cuts->get(pTBin, "max decay length")) {
      return false;
    }

    // candidate maximum decay length XY
    if (hfCandXicc.decayLengthXY() > cuts->get(pTBin, "max decay length XY")) {
      return false;
    }

    // candidate minimum decay length
    if (hfCandXicc.decayLength() <= cuts->get(pTBin, "min decay length")) {
      return false;
    }

    // candidate chi2PC
    if (hfCandXicc.chi2PCA() > cuts->get(pTBin, "chi2PCA")) {
      return false;
    }

    // minimum DCA of daughters
    if ((std::abs(hfCandXicc.impactParameter0()) <= cuts->get(pTBin, "min d0 Xic")) ||
        (std::abs(hfCandXicc.impactParameter1()) <= cuts->get(pTBin, "min d0 Pi"))) {
      return false;
    }

    // maximum DCA of daughters
    if ((std::abs(hfCandXicc.impactParameter0()) > cuts->get(pTBin, "max d0 Xic")) ||
        (std::abs(hfCandXicc.impactParameter1()) > cuts->get(pTBin, "max d0 Pi"))) {
      return false;
    }

    // cut on daughter pT
    if (hfCandXic.pt() < cuts->get(pTBin, "pT Xic") || trackPi.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    return true;
  }

  void process(aod::HfCandXicc const& hfCandXiccs, aod::HfCandProng3 const&, aod::BigTracksPID const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(d_pidTPCMinpT, d_pidTPCMaxpT);
    selectorPion.setRangeNSigmaTPC(-d_nSigmaTPC, d_nSigmaTPC);
    selectorPion.setRangeNSigmaTPCCondTOF(-d_nSigmaTPCCombined, d_nSigmaTPCCombined);
    selectorPion.setRangePtTOF(d_pidTOFMinpT, d_pidTOFMaxpT);
    selectorPion.setRangeNSigmaTOF(-d_nSigmaTOF, d_nSigmaTOF);
    selectorPion.setRangeNSigmaTOFCondTPC(-d_nSigmaTOFCombined, d_nSigmaTOFCombined);

    // looping over 3-prong candidates
    for (auto& hfCandXicc : hfCandXiccs) {
      auto hfCandXic = hfCandXicc.index0();
      auto trackPi = hfCandXicc.index1_as<aod::BigTracksPID>();
      // final selection flag: 0 - rejected, 1 - accepted
      auto statusXiccToPKPiPi = 0;

      if (!(hfCandXicc.hfflag() & 1 << DecayType::XiccToXicPi)) {
        hfSelXiccToPKPiPiCandidate(statusXiccToPKPiPi);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(hfCandXicc, hfCandXic, trackPi)) {
        hfSelXiccToPKPiPiCandidate(statusXiccToPKPiPi);
        continue;
      }

      auto pidPi = 0;

      if (!d_FilterPID) {
        // PID non applied
        pidPi = 1;
      } else {
        // track-level PID selection
        int pidPion = selectorPion.getStatusTrackPIDAll(trackPi);
        if (pidPion == TrackSelectorPID::Status::PIDAccepted) {
          pidPi = 1;
        }
      }

      if (pidPi == 0) {
        hfSelXiccToPKPiPiCandidate(statusXiccToPKPiPi);
        continue;
      }

      hfSelXiccToPKPiPiCandidate(1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfXiccToPKPiPiCandidateSelector>(cfgc)};
}
