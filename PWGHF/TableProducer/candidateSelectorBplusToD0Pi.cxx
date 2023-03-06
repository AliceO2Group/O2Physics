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

/// \file candidateSelectorBplusToD0Pi.cxx
/// \brief B± → D0bar(D0) π± candidate selector
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari
/// \author Deepa Thomas <deepa.thomas@cern.ch>, UT Austin

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_bplus;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_bplus_to_d0_pi;

struct HfCandidateSelectorBplusToD0Pi {
  Produces<aod::HfSelBplusToD0Pi> hfSelBplusToD0PiCandidate;

  // Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  // Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // Enable PID
  Configurable<bool> usePid{"usePid", true, "Bool to use or not the PID at filtering level"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 999, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 9999, "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 50., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 999., "Nsigma cut on TOF combined with TPC"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bplus_to_d0_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bplus_to_d0_pi::cuts[0], nBinsPt, nCutVars, labelsPt, labelsCutVar}, "B+ candidate selection per pT bin"};

  // Apply topological cuts as defined in SelectorCuts.h; return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandBplus, const T2& hfCandD0, const T3& trackPi)
  {
    auto candpT = hfCandBplus.pt();
    int pTBin = findBin(binsPt, candpT);
    if (pTBin == -1) {
      // Printf("B+ topol selection failed at getpTBin");
      return false;
    }

    // pi pt
    if (trackPi.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    // d0(D0)xd0(pi)
    if (hfCandBplus.impactParameterProduct() > cuts->get(pTBin, "Imp. Par. Product")) {
      return false;
    }

    // D0 mass
    if (trackPi.sign() > 0) {
      if (std::abs(invMassD0barToKPi(hfCandD0) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "DeltaMD0")) {
        return false;
      }
    }
    if (trackPi.sign() < 0) {
      if (std::abs(invMassD0ToPiK(hfCandD0) - RecoDecay::getMassPDG(pdg::Code::kD0)) > cuts->get(pTBin, "DeltaMD0")) {
        return false;
      }
    }

    // B Decay length
    if (hfCandBplus.decayLength() < cuts->get(pTBin, "B decLen")) {
      return false;
    }

    // B Decay length XY
    if (hfCandBplus.decayLengthXY() < cuts->get(pTBin, "B decLenXY")) {
      return false;
    }

    // B+ CPA cut
    if (hfCandBplus.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    // if (candpT < ptCandMin || candpT >= ptCandMax) {
    //  Printf("B+ topol selection failed at cand pT check");
    //  return false;
    //  }

    // B+ mass cut
    // if (std::abs(invMassBplusToD0Pi(hfCandBplus) - RecoDecay::getMassPDG(521)) > cuts->get(pTBin, "m")) {
    //  Printf("B+ topol selection failed at mass diff check");
    //   return false;
    //  }

    // d0 of D0 and pi
    if (std::abs(hfCandBplus.impactParameter0()) < cuts->get(pTBin, "d0 D0")) {
      return false;
    }
    if (std::abs(hfCandBplus.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }
    // D0 CPA
    //  if (std::abs(hfCandD0.cpa()) < cuts->get(pTBin, "CPA D0")){
    //   return false;
    // }
    return true;
  }

  void process(aod::HfCandBplus const& hfCandBs, soa::Join<aod::HfCand2Prong, aod::HfSelD0> const&, aod::BigTracksPID const& tracks)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTPCCondTOF(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTOF(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTOF(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTOFCondTPC(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);

    for (auto& hfCandB : hfCandBs) { // looping over Bplus candidates

      int statusBplus = 0;

      // check if flagged as B+ --> D0bar Pi
      if (!(hfCandB.hfflag() & 1 << hf_cand_bplus::DecayType::BplusToD0Pi)) {
        hfSelBplusToD0PiCandidate(statusBplus);
        // Printf("B+ candidate selection failed at hfflag check");
        continue;
      }

      // D0 is always index0 and pi is index1 by default
      auto candD0 = hfCandB.prong0_as<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>();
      auto trackPi = hfCandB.prong1_as<aod::BigTracksPID>();

      // topological cuts
      if (!selectionTopol(hfCandB, candD0, trackPi)) {
        hfSelBplusToD0PiCandidate(statusBplus);
        // Printf("B+ candidate selection failed at selection topology");
        continue;
      }

      if (usePid) {
        // PID applied
        if (selectorPion.getStatusTrackPIDTpcOrTof(trackPi) != TrackSelectorPID::Status::PIDAccepted) {
          // Printf("PID not successful");
          hfSelBplusToD0PiCandidate(statusBplus);
          continue;
        }
      }

      statusBplus = 1; // Successful topological and PID

      hfSelBplusToD0PiCandidate(statusBplus);
      // Printf("B+ candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfCandidateSelectorBplusToD0Pi>(cfgc));
  return workflow;
}
