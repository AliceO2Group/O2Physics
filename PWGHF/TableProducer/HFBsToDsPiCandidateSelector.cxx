// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFLbToLcPiCandidateSelector.cxx
/// \brief Bs → Ds+ π- candidate selector
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/Core/HFSelectorCuts.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "PWGHF/Core/HFSelectorCuts.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_bs;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::analysis::hf_cuts_bs_todspi;

struct HfBsToDsPiCandidateSelector {
  Produces<aod::HFSelBsToDsPiCandidate> hfSelBsToDsPiCandidate;

  Configurable<double> pTCandMin{"pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> pTCandMax{"pTCandMax", 50., "Upper bound of candidate pT"};

  //Track quality
  Configurable<double> TPCNClsFindablePIDCut{"TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};

  //TPC PID
  Configurable<double> pidTPCMinpT{"pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> pidTPCMaxpT{"pidTPCMaxpT", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTPC{"nSigmaTPC", 5., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTPCCombined{"nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};

  //TOF PID
  Configurable<double> pidTOFMinpT{"pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> pidTOFMaxpT{"pidTOFMaxpT", 10., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTOF{"nSigmaTOF", 5., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTOFCombined{"nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};

  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_bs_todspi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Bs_to_dspi_cuts", {hf_cuts_bs_todspi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Bs0 candidate selection per pT bin"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 1, "Selection Flag for Ds+"};

  // Apply topological cuts as defined in HFSelectorCuts.h; return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandBs, const T2& hfCandDs, const T3& trackPi)
  {
    auto candpT = hfCandBs.pt();
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      // Printf("Bs topol selection failed at getpTBin");
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (candpT < pTCandMin || candpT >= pTCandMax) {
      return false;
    }

    //Λb0 mass cut
    if (std::abs(InvMassBsToDsPi(hfCandBs) - RecoDecay::getMassPDG(pdg::Code::kLambdaB0)) > cuts->get(pTBin, "m")) {
      //Printf("Bs topol selection failed at mass diff check");
      return false;
    }

    //pion pt
    if (trackPi.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }

    //Ds+ pt
    if (hfCandDs.pt() < cuts->get(pTBin, "pT Ds+")) {
      return false;
    }

    //Ds mass
    //if (trackPi.sign() < 0) {
    //if (std::abs(InvMassDspKpi(hfCandDs) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "DeltaMDs")) {
    //return false;
    //}
    //}

    //Bs Decay length
    if (hfCandBs.decayLength() < cuts->get(pTBin, "Bs decLen")) {
      return false;
    }

    //Bs Decay length XY
    if (hfCandBs.decayLengthXY() < cuts->get(pTBin, "Bs decLenXY")) {
      return false;
    }

    //Bs chi2PCA cut
    if (hfCandBs.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      //Printf("Bs selection failed at chi2PCA");
      return false;
    }

    //Bs CPA cut
    if (hfCandBs.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }

    //d0 of pi
    if (std::abs(hfCandBs.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    //d0 of Ds+
    if (std::abs(hfCandBs.impactParameter0()) < cuts->get(pTBin, "d0 Ds+")) {
      return false;
    }

    return true;
  }

  void process(aod::HfCandBs const& hfCandBss, soa::Join<aod::HfCandProng3, aod::HFSelDsCandidate>, aod::BigTracksPID const&)
  {
    for (auto& hfCandBs : hfCandBss) { //looping over Bs candidates

      int statusBs = 0;

      // check if flagged as Λb --> Λc+ π-
      if (!(hfCandBs.hfflag() & 1 << hf_cand_bs::DecayType::BsToDsPi)) {
        hfSelBsToDsPiCandidate(statusBs);
        //Printf("Bs candidate selection failed at hfflag check");
        continue;
      }

      // Ds is always index0 and pi is index1 by default
      //auto candDs = hfCandBs.index0();
      auto candDs = hfCandBs.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelDsCandidate>>();
      auto trackPi = hfCandBs.index1_as<aod::BigTracksPID>();

      //topological cuts
      if (!selectionTopol(hfCandBs, candDs, trackPi)) {
        hfSelBsToDsPiCandidate(statusBs);
        // Printf("Bs candidate selection failed at selection topology");
        continue;
      }

      hfSelBsToDsPiCandidate(1);
      //Printf("Bs candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfBsToDsPiCandidateSelector>(cfgc));
  return workflow;
}
