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

/// \file CandidateSelectionTables.h
/// \brief Definitions of tables produced by candidate selectors

#ifndef PWGHF_DATAMODEL_CANDIDATESELECTIONTABLES_H_
#define PWGHF_DATAMODEL_CANDIDATESELECTIONTABLES_H_

#include <vector>

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

namespace o2::aod
{
// selection steps
enum SelectionStep {
  RecoSkims = 0,
  RecoTopol,
  RecoPID,
  RecoMl,
  NSelectionSteps
};

namespace hf_sel_candidate_d0
{
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);           //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int);     //!
DECLARE_SOA_COLUMN(IsRecoHfFlag, isRecoHfFlag, int); //!
DECLARE_SOA_COLUMN(IsRecoTopol, isRecoTopol, int);   //!
DECLARE_SOA_COLUMN(IsRecoCand, isRecoCand, int);     //!
DECLARE_SOA_COLUMN(IsRecoPid, isRecoPid, int);
DECLARE_SOA_COLUMN(MlProbD0, mlProbD0, std::vector<float>); //!
} // namespace hf_sel_candidate_d0
DECLARE_SOA_TABLE(HfSelD0, "AOD", "HFSELD0", //!
                  hf_sel_candidate_d0::IsSelD0,
                  hf_sel_candidate_d0::IsSelD0bar,
                  hf_sel_candidate_d0::IsRecoHfFlag,
                  hf_sel_candidate_d0::IsRecoTopol,
                  hf_sel_candidate_d0::IsRecoCand,
                  hf_sel_candidate_d0::IsRecoPid);
DECLARE_SOA_TABLE(HfMlD0, "AOD", "HFMLD0", //!
                  hf_sel_candidate_d0::MlProbD0);

namespace hf_sel_candidate_d0_parametrized_pid
{
DECLARE_SOA_COLUMN(IsSelD0NoPid, isSelD0NoPid, int);                 //!
DECLARE_SOA_COLUMN(IsSelD0PerfectPid, isSelD0PerfectPid, int);       //!
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);                           //!
DECLARE_SOA_COLUMN(IsSelD0barNoPid, isSelD0barNoPid, int);           //!
DECLARE_SOA_COLUMN(IsSelD0barPerfectPid, isSelD0barPerfectPid, int); //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int);                     //!
} // namespace hf_sel_candidate_d0_parametrized_pid
DECLARE_SOA_TABLE(HfSelD0ParametrizedPid, "AOD", "HFSELD0P", //!
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0NoPid,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0PerfectPid,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0barNoPid,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0barPerfectPid,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0bar);

namespace hf_sel_candidate_d0_alice3_barrel
{
DECLARE_SOA_COLUMN(IsSelHfFlag, isSelHfFlag, int);                           //!
DECLARE_SOA_COLUMN(IsSelD0NoPid, isSelD0NoPid, int);                         //!
DECLARE_SOA_COLUMN(IsSelD0PerfectPid, isSelD0PerfectPid, int);               //!
DECLARE_SOA_COLUMN(IsSelD0TofPid, isSelD0TofPid, int);                       //!
DECLARE_SOA_COLUMN(IsSelD0RichPid, isSelD0RichPid, int);                     //!
DECLARE_SOA_COLUMN(IsSelD0TofPlusRichPid, isSelD0TofPlusRichPid, int);       //!
DECLARE_SOA_COLUMN(IsSelD0barTofPlusRichPid, isSelD0barTofPlusRichPid, int); //!
} // namespace hf_sel_candidate_d0_alice3_barrel
DECLARE_SOA_TABLE(HfSelD0Alice3Barrel, "AOD", "HFSELD0A3B", //!
                  hf_sel_candidate_d0_alice3_barrel::IsSelHfFlag,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0NoPid,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0PerfectPid,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0TofPid,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0RichPid,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0TofPlusRichPid,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0barTofPlusRichPid);

namespace hf_sel_candidate_d0_alice3_forward
{
DECLARE_SOA_COLUMN(IsSelHfFlag, isSelHfFlag, int);         //!
DECLARE_SOA_COLUMN(IsSelD0FNoPid, isSelD0FNoPid, int);     //!
DECLARE_SOA_COLUMN(IsSelD0FRichPid, isSelD0FRichPid, int); //!
} // namespace hf_sel_candidate_d0_alice3_forward
DECLARE_SOA_TABLE(HfSelD0Alice3Forward, "AOD", "HFSELD0A3F", //!
                  hf_sel_candidate_d0_alice3_forward::IsSelHfFlag,
                  hf_sel_candidate_d0_alice3_forward::IsSelD0FNoPid,
                  hf_sel_candidate_d0_alice3_forward::IsSelD0FRichPid);

namespace hf_sel_candidate_dplus
{
DECLARE_SOA_COLUMN(IsSelDplusToPiKPi, isSelDplusToPiKPi, int); //!
DECLARE_SOA_COLUMN(MlProbDplusToPiKPi, mlProbDplusToPiKPi, std::vector<float>); //!
} // namespace hf_sel_candidate_dplus
DECLARE_SOA_TABLE(HfSelDplusToPiKPi, "AOD", "HFSELDPLUS", //!
                  hf_sel_candidate_dplus::IsSelDplusToPiKPi);
DECLARE_SOA_TABLE(HfMlDplusToPiKPi, "AOD", "HFMLDPLUS", //!
                  hf_sel_candidate_dplus::MlProbDplusToPiKPi);

namespace hf_sel_candidate_ds
{
DECLARE_SOA_COLUMN(IsSelDsToKKPi, isSelDsToKKPi, int); //!
DECLARE_SOA_COLUMN(IsSelDsToPiKK, isSelDsToPiKK, int); //!
DECLARE_SOA_COLUMN(MlProbDsToKKPi, mlProbDsToKKPi, std::vector<float>); //!
} // namespace hf_sel_candidate_ds
DECLARE_SOA_TABLE(HfSelDsToKKPi, "AOD", "HFSELDS", //!
                  hf_sel_candidate_ds::IsSelDsToKKPi, hf_sel_candidate_ds::IsSelDsToPiKK);
DECLARE_SOA_TABLE(HfMlDsToKKPi, "AOD", "HFMLDS", //!
                  hf_sel_candidate_ds::MlProbDsToKKPi);

namespace hf_sel_candidate_lc
{
DECLARE_SOA_COLUMN(IsSelLcToPKPi, isSelLcToPKPi, int); //!
DECLARE_SOA_COLUMN(IsSelLcToPiKP, isSelLcToPiKP, int); //!
} // namespace hf_sel_candidate_lc
DECLARE_SOA_TABLE(HfSelLc, "AOD", "HFSELLC", //!
                  hf_sel_candidate_lc::IsSelLcToPKPi, hf_sel_candidate_lc::IsSelLcToPiKP);

namespace hf_sel_candidate_lc_alice3
{
DECLARE_SOA_COLUMN(IsSelLcToPKPiNoPid, isSelLcToPKPiNoPid, int);                   //!
DECLARE_SOA_COLUMN(IsSelLcToPKPiPerfectPid, isSelLcToPKPiPerfectPid, int);         //!
DECLARE_SOA_COLUMN(IsSelLcToPKPiTofPid, isSelLcToPKPiTofPid, int);                 //!
DECLARE_SOA_COLUMN(IsSelLcToPKPiTofPlusRichPid, isSelLcToPKPiTofPlusRichPid, int); //!
DECLARE_SOA_COLUMN(IsSelLcToPiKPNoPid, isSelLcToPiKPNoPid, int);                   //!
DECLARE_SOA_COLUMN(IsSelLcToPiKPPerfectPid, isSelLcToPiKPPerfectPid, int);         //!
DECLARE_SOA_COLUMN(IsSelLcToPiKPTofPid, isSelLcToPiKPTofPid, int);                 //!
DECLARE_SOA_COLUMN(IsSelLcToPiKPTofPlusRichPid, isSelLcToPiKPTofPlusRichPid, int); //!
} // namespace hf_sel_candidate_lc_alice3
DECLARE_SOA_TABLE(HfSelLcAlice3, "AOD", "HFSELLCA3B", //!
                  hf_sel_candidate_lc_alice3::IsSelLcToPKPiNoPid,
                  hf_sel_candidate_lc_alice3::IsSelLcToPKPiPerfectPid,
                  hf_sel_candidate_lc_alice3::IsSelLcToPKPiTofPid,
                  hf_sel_candidate_lc_alice3::IsSelLcToPKPiTofPlusRichPid,
                  hf_sel_candidate_lc_alice3::IsSelLcToPiKPNoPid,
                  hf_sel_candidate_lc_alice3::IsSelLcToPiKPPerfectPid,
                  hf_sel_candidate_lc_alice3::IsSelLcToPiKPTofPid,
                  hf_sel_candidate_lc_alice3::IsSelLcToPiKPTofPlusRichPid);

namespace hf_sel_candidate_lc_parametrized_pid
{
DECLARE_SOA_COLUMN(IsSelLcToPKPiNoPid, isSelLcToPKPiNoPid, int);           //!
DECLARE_SOA_COLUMN(IsSelLcToPKPiPerfectPid, isSelLcToPKPiPerfectPid, int); //!
DECLARE_SOA_COLUMN(IsSelLcToPKPi, isSelLcToPKPi, int);                     //!
DECLARE_SOA_COLUMN(IsSelLcToPiKPNoPid, isSelLcToPiKPNoPid, int);           //!
DECLARE_SOA_COLUMN(IsSelLcToPiKPPerfectPid, isSelLcToPiKPPerfectPid, int); //!
DECLARE_SOA_COLUMN(IsSelLcToPiKP, isSelLcToPiKP, int);                     //!
} // namespace hf_sel_candidate_lc_parametrized_pid
DECLARE_SOA_TABLE(HfSelLcParametrizedPid, "AOD", "HFSELLCP", //!
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcToPKPiNoPid,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcToPKPiPerfectPid,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcToPKPi,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcToPiKPNoPid,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcToPiKPPerfectPid,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcToPiKP);

namespace hf_sel_candidate_jpsi
{
DECLARE_SOA_COLUMN(IsSelJpsiToEE, isSelJpsiToEE, int);               //!
DECLARE_SOA_COLUMN(IsSelJpsiToMuMu, isSelJpsiToMuMu, int);           //!
DECLARE_SOA_COLUMN(IsSelJpsiToEETopol, isSelJpsiToEETopol, int);     //!
DECLARE_SOA_COLUMN(IsSelJpsiToEETpc, isSelJpsiToEETpc, int);         //!
DECLARE_SOA_COLUMN(IsSelJpsiToEETof, isSelJpsiToEETof, int);         //!
DECLARE_SOA_COLUMN(IsSelJpsiToEERich, isSelJpsiToEERich, int);       //!
DECLARE_SOA_COLUMN(IsSelJpsiToEETofRich, isSelJpsiToEETofRich, int); //!
DECLARE_SOA_COLUMN(IsSelJpsiToMuMuTopol, isSelJpsiToMuMuTopol, int); //!
DECLARE_SOA_COLUMN(IsSelJpsiToMuMuMid, isSelJpsiToMuMuMid, int);     //!
} // namespace hf_sel_candidate_jpsi
DECLARE_SOA_TABLE(HfSelJpsi, "AOD", "HFSELJPSI", //!
                  hf_sel_candidate_jpsi::IsSelJpsiToEE,
                  hf_sel_candidate_jpsi::IsSelJpsiToMuMu,
                  hf_sel_candidate_jpsi::IsSelJpsiToEETopol,
                  hf_sel_candidate_jpsi::IsSelJpsiToEETpc,
                  hf_sel_candidate_jpsi::IsSelJpsiToEETof,
                  hf_sel_candidate_jpsi::IsSelJpsiToEERich,
                  hf_sel_candidate_jpsi::IsSelJpsiToEETofRich,
                  hf_sel_candidate_jpsi::IsSelJpsiToMuMuTopol,
                  hf_sel_candidate_jpsi::IsSelJpsiToMuMuMid);

namespace hf_sel_candidate_lc_to_k0s_p
{
DECLARE_SOA_COLUMN(IsSelLcToK0sP, isSelLcToK0sP, int);
} // namespace hf_sel_candidate_lc_to_k0s_p
DECLARE_SOA_TABLE(HfSelLcToK0sP, "AOD", "HFSELLCK0SP", //!
                  hf_sel_candidate_lc_to_k0s_p::IsSelLcToK0sP);

namespace hf_sel_candidate_b0
{
DECLARE_SOA_COLUMN(IsSelB0ToDPi, isSelB0ToDPi, int); //!

/// Apply topological cuts as defined in SelectorCuts.h
/// \param candB0 B0 candidate
/// \param cuts B0 candidate selection per pT bin"
/// \param binsPt pT bin limits
/// \return true if candidate passes all selections
template <typename T1, typename T2, typename T3>
bool selectionTopol(const T1& candB0, const T2& cuts, const T3& binsPt)
{
  auto ptCandB0 = candB0.pt();
  auto ptD = RecoDecay::pt(candB0.pxProng0(), candB0.pyProng0());
  auto ptPi = RecoDecay::pt(candB0.pxProng1(), candB0.pyProng1());

  int pTBin = findBin(binsPt, ptCandB0);
  if (pTBin == -1) {
    // LOGF(info, "B0 topol selection failed at getpTBin");
    return false;
  }

  // // check that the candidate pT is within the analysis range
  // if (ptCandB0 < ptCandMin || ptCandB0 >= ptCandMax) {
  //   return false;
  // }

  // B0 mass cut
  if (std::abs(o2::aod::hf_cand_b0::invMassB0ToDPi(candB0) - RecoDecay::getMassPDG(o2::analysis::pdg::Code::kB0)) > cuts->get(pTBin, "m")) {
    // Printf("B0 topol selection failed at mass diff check");
    return false;
  }

  // pion pt
  if (ptPi < cuts->get(pTBin, "pT Pi")) {
    return false;
  }

  // D- pt
  if (ptD < cuts->get(pTBin, "pT D")) {
    return false;
  }

  /*
  // D mass cut | already applied in candidateSelectorDplusToPiKPi.cxx
  if (std::abs(o2::aod::hf_cand_3prong::invMassDplusToPiKPi(hfCandD) - RecoDecay::getMassPDG(o2::analysis::pdg::Code::kDMinus)) > cuts->get(pTBin, "DeltaMD")) {
    return false;
  }
  */

  // B0 Decay length
  if (candB0.decayLength() < cuts->get(pTBin, "B0 decLen")) {
    return false;
  }

  // B0 Decay length XY
  if (candB0.decayLengthXY() < cuts->get(pTBin, "B0 decLenXY")) {
    return false;
  }

  // B0 chi2PCA cut
  if (candB0.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
    return false;
  }

  // B0 CPA cut
  if (candB0.cpa() < cuts->get(pTBin, "CPA")) {
    return false;
  }

  // d0 of pi
  if (std::abs(candB0.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
    return false;
  }

  // d0 of D
  if (std::abs(candB0.impactParameter0()) < cuts->get(pTBin, "d0 D")) {
    return false;
  }

  return true;
}

/// Apply PID selection
/// \param pidTrackPi PID status of trackPi (prong1 of B0 candidate)
/// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
/// \return true if prong1 of B0 candidate passes all selections
template <typename T1 = int, typename T2 = bool>
bool selectionPID(const T1& pidTrackPi, const T2& acceptPIDNotApplicable)
{
  if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
    return false;
  }
  if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
    return false;
  }

  return true;
}

} // namespace hf_sel_candidate_b0
DECLARE_SOA_TABLE(HfSelB0ToDPi, "AOD", "HFSELB0", //!
                  hf_sel_candidate_b0::IsSelB0ToDPi);

namespace hf_sel_candidate_bs
{
DECLARE_SOA_COLUMN(IsSelBsToDsPi, isSelBsToDsPi, int); //!

/// Apply topological cuts as defined in SelectorCuts.h
/// \param candBs Bs candidate
/// \param cuts Bs candidate selections
/// \param binsPt pT bin limits
/// \return true if candidate passes all selections
template <typename T1, typename T2, typename T3>
bool selectionTopol(const T1& candBs, const T2& cuts, const T3& binsPt)
{
  auto ptCandBs = candBs.pt();
  auto ptDs = RecoDecay::pt(candBs.pxProng0(), candBs.pyProng0());
  auto ptPi = RecoDecay::pt(candBs.pxProng1(), candBs.pyProng1());

  int pTBin = findBin(binsPt, ptCandBs);
  if (pTBin == -1) {
    return false;
  }

  // Bs mass cut
  if (std::abs(o2::aod::hf_cand_bs::invMassBsToDsPi(candBs) - RecoDecay::getMassPDG(o2::analysis::pdg::Code::kBS)) > cuts->get(pTBin, "m")) {
    return false;
  }

  // pion pt
  if (ptPi < cuts->get(pTBin, "pT Pi")) {
    return false;
  }

  // Ds pt
  if (ptDs < cuts->get(pTBin, "pT Ds")) {
    return false;
  }

  // Bs Decay length
  if (candBs.decayLength() < cuts->get(pTBin, "Bs decLen")) {
    return false;
  }

  // Bs Decay length XY
  if (candBs.decayLengthXY() < cuts->get(pTBin, "Bs decLenXY")) {
    return false;
  }

  // Bs chi2PCA cut
  if (candBs.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
    return false;
  }

  // Bs CPA cut
  if (candBs.cpa() < cuts->get(pTBin, "CPA")) {
    return false;
  }

  // d0 of pi
  if (std::abs(candBs.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
    return false;
  }

  // d0 of Ds
  if (std::abs(candBs.impactParameter0()) < cuts->get(pTBin, "d0 Ds")) {
    return false;
  }

  // d0(Ds)xd0(pi)
  if (candBs.impactParameterProduct() > cuts->get(pTBin, "Imp. Par. Product")) {
    return false;
  }

  return true;
}

/// Apply PID selection
/// \param pidTrackPi PID status of trackPi (prong1 of Bs candidate)
/// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
/// \return true if prong1 of Bs candidate passes all selections
template <typename T1 = int, typename T2 = bool>
bool selectionPID(const T1& pidTrackPi, const T2& acceptPIDNotApplicable)
{
  if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
    return false;
  }
  if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
    return false;
  }

  return true;
}

} // namespace hf_sel_candidate_bs
DECLARE_SOA_TABLE(HfSelBsToDsPi, "AOD", "HFSELBS", //!
                  hf_sel_candidate_bs::IsSelBsToDsPi);

namespace hf_sel_candidate_bplus
{
DECLARE_SOA_COLUMN(IsSelBplusToD0Pi, isSelBplusToD0Pi, int); //!
// Apply topological cuts as defined in SelectorCuts.h
/// \param candBp B+ candidate
/// \param cuts B+ candidate selection per pT bin"
/// \param binsPt pT bin limits
/// \return true if candidate passes all selections
template <typename T1, typename T2, typename T3>
bool selectionTopol(const T1& candBp, const T2& cuts, const T3& binsPt)
{
  auto ptcandBp = candBp.pt();
  auto ptPi = RecoDecay::pt(candBp.pxProng1(), candBp.pyProng1());

  int pTBin = findBin(binsPt, ptcandBp);
  if (pTBin == -1) {
    return false;
  }

  // B+ mass cut
  if (std::abs(invMassBplusToD0Pi(candBp) - RecoDecay::getMassPDG(521)) > cuts->get(pTBin, "m")) {
    return false;
  }

  // pion pt
  if (ptPi < cuts->get(pTBin, "pT Pi")) {
    return false;
  }

  // d0(D0)xd0(pi)
  if (candBp.impactParameterProduct() > cuts->get(pTBin, "Imp. Par. Product")) {
    return false;
  }

  // B Decay length
  if (candBp.decayLength() < cuts->get(pTBin, "B decLen")) {
    return false;
  }

  // B Decay length XY
  if (candBp.decayLengthXY() < cuts->get(pTBin, "B decLenXY")) {
    return false;
  }

  // B+ CPA cut
  if (candBp.cpa() < cuts->get(pTBin, "CPA")) {
    return false;
  }

  return true;
}

/// Apply PID selection
/// \param pidTrackPi PID status of trackPi (prong1 of B+ candidate)
/// \param acceptPIDNotApplicable switch to accept Status::NotApplicable
/// \return true if prong1 of B+ candidate passes all selections
template <typename T1 = int, typename T2 = bool>
bool selectionPID(const T1& pidTrackPi, const T2& acceptPIDNotApplicable)
{
  if (!acceptPIDNotApplicable && pidTrackPi != TrackSelectorPID::Accepted) {
    return false;
  }
  if (acceptPIDNotApplicable && pidTrackPi == TrackSelectorPID::Rejected) {
    return false;
  }

  return true;
}

} // namespace hf_sel_candidate_bplus
DECLARE_SOA_TABLE(HfSelBplusToD0Pi, "AOD", "HFSELBPLUS", //!
                  hf_sel_candidate_bplus::IsSelBplusToD0Pi);

namespace hf_sel_candidate_lb
{
DECLARE_SOA_COLUMN(IsSelLbToLcPi, isSelLbToLcPi, int); //!
} // namespace hf_sel_candidate_lb
DECLARE_SOA_TABLE(HfSelLbToLcPi, "AOD", "HFSELLB", //!
                  hf_sel_candidate_lb::IsSelLbToLcPi);

namespace hf_sel_candidate_x
{
DECLARE_SOA_COLUMN(IsSelXToJpsiToEEPiPi, isSelXToJpsiToEEPiPi, int);     //!
DECLARE_SOA_COLUMN(IsSelXToJpsiToMuMuPiPi, isSelXToJpsiToMuMuPiPi, int); //!
} // namespace hf_sel_candidate_x
DECLARE_SOA_TABLE(HfSelXToJpsiPiPi, "AOD", "HFSELX", //!
                  hf_sel_candidate_x::IsSelXToJpsiToEEPiPi, hf_sel_candidate_x::IsSelXToJpsiToMuMuPiPi);

namespace hf_sel_candidate_chic
{
DECLARE_SOA_COLUMN(IsSelChicToJpsiToEEGamma, isSelChicToJpsiToEEGamma, int);     //!
DECLARE_SOA_COLUMN(IsSelChicToJpsiToMuMuGamma, isSelChicToJpsiToMuMuGamma, int); //!
} // namespace hf_sel_candidate_chic
DECLARE_SOA_TABLE(HfSelChicToJpsiGamma, "AOD", "HFSELCHIC", //!
                  hf_sel_candidate_chic::IsSelChicToJpsiToEEGamma, hf_sel_candidate_chic::IsSelChicToJpsiToMuMuGamma);

namespace hf_sel_candidate_xic
{
DECLARE_SOA_COLUMN(IsSelXicToPKPi, isSelXicToPKPi, int); //!
DECLARE_SOA_COLUMN(IsSelXicToPiKP, isSelXicToPiKP, int); //!
} // namespace hf_sel_candidate_xic
DECLARE_SOA_TABLE(HfSelXicToPKPi, "AOD", "HFSELXIC", //!
                  hf_sel_candidate_xic::IsSelXicToPKPi, hf_sel_candidate_xic::IsSelXicToPiKP);

namespace hf_sel_candidate_xicc
{
DECLARE_SOA_COLUMN(IsSelXiccToPKPiPi, isSelXiccToPKPiPi, int); //!
} // namespace hf_sel_candidate_xicc
DECLARE_SOA_TABLE(HfSelXiccToPKPiPi, "AOD", "HFSELXICC", //!
                  hf_sel_candidate_xicc::IsSelXiccToPKPiPi);

namespace hf_sel_toxipi
{
DECLARE_SOA_COLUMN(StatusPidLambda, statusPidLambda, bool);
DECLARE_SOA_COLUMN(StatusPidCascade, statusPidCascade, bool);
DECLARE_SOA_COLUMN(StatusPidOmegac, statusPidOmegac, bool);
DECLARE_SOA_COLUMN(StatusInvMassLambda, statusInvMassLambda, bool);
DECLARE_SOA_COLUMN(StatusInvMassCascade, statusInvMassCascade, bool);
DECLARE_SOA_COLUMN(StatusInvMassOmegac, statusInvMassOmegac, bool);
DECLARE_SOA_COLUMN(ResultSelections, resultSelections, bool);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromOmega, tpcNSigmaPiFromOmega, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCasc, tpcNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromOmega, tofNSigmaPiFromOmega, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCasc, tofNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);

} // namespace hf_sel_toxipi
DECLARE_SOA_TABLE(HfSelToXiPi, "AOD", "HFSELTOXIPI",
                  hf_sel_toxipi::StatusPidLambda, hf_sel_toxipi::StatusPidCascade, hf_sel_toxipi::StatusPidOmegac,
                  hf_sel_toxipi::StatusInvMassLambda, hf_sel_toxipi::StatusInvMassCascade, hf_sel_toxipi::StatusInvMassOmegac,
                  hf_sel_toxipi::ResultSelections, hf_sel_toxipi::PidTpcInfoStored, hf_sel_toxipi::PidTofInfoStored,
                  hf_sel_toxipi::TpcNSigmaPiFromOmega, hf_sel_toxipi::TpcNSigmaPiFromCasc, hf_sel_toxipi::TpcNSigmaPiFromLambda, hf_sel_toxipi::TpcNSigmaPrFromLambda,
                  hf_sel_toxipi::TofNSigmaPiFromOmega, hf_sel_toxipi::TofNSigmaPiFromCasc, hf_sel_toxipi::TofNSigmaPiFromLambda, hf_sel_toxipi::TofNSigmaPrFromLambda);

} // namespace o2::aod
#endif // PWGHF_DATAMODEL_CANDIDATESELECTIONTABLES_H_
