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

/// \file candidateSelectorXiccToPKPiPi.cxx
/// \brief Ξcc± → p± K∓ π± π±selection task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN

#include <vector>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Struct for applying Xicc selection cuts
struct HfCandidateSelectorXiccToPKPiPi {
  Produces<aod::HfSelXiccToPKPiPi> hfSelXiccToPKPiPiCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> usePid{"usePid", true, "Bool to use or not the PID at filtering level"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 9999., "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 99999., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 9999., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 9999., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 4., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xicc_to_p_k_pi_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_xicc_to_p_k_pi_pi::cuts[0], hf_cuts_xicc_to_p_k_pi_pi::nBinsPt, hf_cuts_xicc_to_p_k_pi_pi::nCutVars, hf_cuts_xicc_to_p_k_pi_pi::labelsPt, hf_cuts_xicc_to_p_k_pi_pi::labelsCutVar}, "Xicc candidate selection per pT bin"};

  HfHelper hfHelper;
  TrackSelectorPi selectorPion;

  using TracksSel = soa::Join<aod::Tracks, aod::TracksPidPi>;

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
  }

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandXicc, const T2& hfCandXic, const T3& trackPi)
  {
    auto candpT = hfCandXicc.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < ptCandMin || candpT >= ptCandMax) {
      return false;
    }

    // check candidate mass is within a defined mass window
    if (std::abs(hfHelper.invMassXiccToXicPi(hfCandXicc) - o2::constants::physics::MassXiCCPlusPlus) > cuts->get(pTBin, "m")) {
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

  void process(aod::HfCandXicc const& hfCandXiccs,
               aod::HfCand3Prong const&,
               TracksSel const&)
  {
    // looping over 3-prong candidates
    for (const auto& hfCandXicc : hfCandXiccs) {
      auto hfCandXic = hfCandXicc.prong0();
      auto trackPi = hfCandXicc.prong1_as<TracksSel>();
      // final selection flag: 0 - rejected, 1 - accepted
      auto statusXiccToPKPiPi = 0;

      if (!(hfCandXicc.hfflag() & 1 << aod::hf_cand_xicc::DecayType::XiccToXicPi)) {
        hfSelXiccToPKPiPiCandidate(statusXiccToPKPiPi);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol(hfCandXicc, hfCandXic, trackPi)) {
        hfSelXiccToPKPiPiCandidate(statusXiccToPKPiPi);
        continue;
      }

      auto pidPi = 0;

      if (!usePid) {
        // PID non applied
        pidPi = 1;
      } else {
        // track-level PID selection
        int pidPion = selectorPion.statusTpcOrTof(trackPi);
        if (pidPion == TrackSelectorPID::Accepted) {
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
    adaptAnalysisTask<HfCandidateSelectorXiccToPKPiPi>(cfgc)};
}
