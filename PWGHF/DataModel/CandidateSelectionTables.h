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

#ifndef O2_ANALYSIS_CANDIDATESELECTIONTABLES_H_
#define O2_ANALYSIS_CANDIDATESELECTIONTABLES_H_

namespace o2::aod
{
// selection steps
enum SelectionStep {
  RecoSkims = 0,
  RecoTopol,
  RecoPID,
  NSelectionSteps
};

namespace hf_sel_candidate_d0
{
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);           //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int);     //!
DECLARE_SOA_COLUMN(IsRecoHFFlag, isRecoHFFlag, int); //!
DECLARE_SOA_COLUMN(IsRecoTopol, isRecoTopol, int);   //!
DECLARE_SOA_COLUMN(IsRecoCand, isRecoCand, int);     //!
DECLARE_SOA_COLUMN(IsRecoPID, isRecoPID, int);
} // namespace hf_sel_candidate_d0
DECLARE_SOA_TABLE(HfSelD0, "AOD", "HFSELD0", //!
                  hf_sel_candidate_d0::IsSelD0,
                  hf_sel_candidate_d0::IsSelD0bar,
                  hf_sel_candidate_d0::IsRecoHFFlag,
                  hf_sel_candidate_d0::IsRecoTopol,
                  hf_sel_candidate_d0::IsRecoCand,
                  hf_sel_candidate_d0::IsRecoPID);

namespace hf_sel_candidate_d0_parametrized_pid
{
DECLARE_SOA_COLUMN(IsSelD0NoPID, isSelD0NoPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelD0PerfectPID, isSelD0PerfectPID, int);       //!
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);                           //!
DECLARE_SOA_COLUMN(IsSelD0barNoPID, isSelD0barNoPID, int);           //!
DECLARE_SOA_COLUMN(IsSelD0barPerfectPID, isSelD0barPerfectPID, int); //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int);                     //!
} // namespace hf_sel_candidate_d0_parametrized_pid
DECLARE_SOA_TABLE(HfSelD0ParametrizedPid, "AOD", "HFSELD0P", //!
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0NoPID,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0PerfectPID,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0barNoPID,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0barPerfectPID,
                  hf_sel_candidate_d0_parametrized_pid::IsSelD0bar);

namespace hf_sel_candidate_d0_alice3_barrel
{
DECLARE_SOA_COLUMN(IsSelHFFlag, isSelHFFlag, int);                     //!
DECLARE_SOA_COLUMN(IsSelD0NoPID, isSelD0NoPID, int);                   //!
DECLARE_SOA_COLUMN(IsSelD0PerfectPID, isSelD0PerfectPID, int);         //!
DECLARE_SOA_COLUMN(IsSelD0TOFPID, isSelD0TOFPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelD0RICHPID, isSelD0RICHPID, int);               //!
DECLARE_SOA_COLUMN(IsSelD0TOFplusRICHPID, isSelD0TOFplusRICHPID, int); //!
DECLARE_SOA_COLUMN(IsSelD0barTOFplusRICHPID, isSelD0barTOFplusRICHPID, int); //!
} // namespace hf_sel_candidate_d0_alice3_barrel
DECLARE_SOA_TABLE(HfSelD0Alice3Barrel, "AOD", "HFSELD0A3B", //!
                  hf_sel_candidate_d0_alice3_barrel::IsSelHFFlag,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0NoPID,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0PerfectPID,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0TOFPID,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0RICHPID,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0TOFplusRICHPID,
                  hf_sel_candidate_d0_alice3_barrel::IsSelD0barTOFplusRICHPID);

namespace hf_sel_candidate_d0_alice3_forward
{
DECLARE_SOA_COLUMN(IsSelHFFFlag, isSelHFFFlag, int);       //!
DECLARE_SOA_COLUMN(IsSelD0FNoPID, isSelD0FNoPID, int);     //!
DECLARE_SOA_COLUMN(IsSelD0FRICHPID, isSelD0FRICHPID, int); //!
} // namespace hf_sel_candidate_d0_alice3_forward
DECLARE_SOA_TABLE(HfSelD0Alice3Forward, "AOD", "HFSELD0A3F", //!
                  hf_sel_candidate_d0_alice3_forward::IsSelHFFFlag,
                  hf_sel_candidate_d0_alice3_forward::IsSelD0FNoPID,
                  hf_sel_candidate_d0_alice3_forward::IsSelD0FRICHPID);

namespace hf_sel_candidate_dplus
{
DECLARE_SOA_COLUMN(IsSelDplusToPiKPi, isSelDplusToPiKPi, int); //!
} // namespace hf_sel_candidate_dplus
DECLARE_SOA_TABLE(HfSelDplusToPiKPi, "AOD", "HFSELDPLUS", //!
                  hf_sel_candidate_dplus::IsSelDplusToPiKPi);

namespace hf_sel_candidate_ds
{
DECLARE_SOA_COLUMN(IsSelDsToKKPi, isSelDsToKKPi, int); //!
DECLARE_SOA_COLUMN(IsSelDsToPiKK, isSelDsToPiKK, int); //!
} // namespace hf_sel_candidate_ds
DECLARE_SOA_TABLE(HfSelDsToKKPi, "AOD", "HFSELDS", //!
                  hf_sel_candidate_ds::IsSelDsToKKPi, hf_sel_candidate_ds::IsSelDsToPiKK);

namespace hf_sel_candidate_lc
{
DECLARE_SOA_COLUMN(IsSelLcpKpi, isSelLcpKpi, int); //!
DECLARE_SOA_COLUMN(IsSelLcpiKp, isSelLcpiKp, int); //!
} // namespace hf_sel_candidate_lc
DECLARE_SOA_TABLE(HfSelLc, "AOD", "HFSELLC", //!
                  hf_sel_candidate_lc::IsSelLcpKpi, hf_sel_candidate_lc::IsSelLcpiKp);

namespace hf_sel_candidate_lc_alice3
{
DECLARE_SOA_COLUMN(IsSelLcPKPiNoPID, isSelLcPKPiNoPID, int);                   //!
DECLARE_SOA_COLUMN(IsSelLcPKPiPerfectPID, isSelLcPKPiPerfectPID, int);         //!
DECLARE_SOA_COLUMN(IsSelLcPKPiTOFPID, isSelLcPKPiTOFPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelLcPKPiTOFplusRICHPID, isSelLcPKPiTOFplusRICHPID, int); //!
DECLARE_SOA_COLUMN(IsSelLcPiKPNoPID, isSelLcPiKPNoPID, int);                   //!
DECLARE_SOA_COLUMN(IsSelLcPiKPPerfectPID, isSelLcPiKPPerfectPID, int);         //!
DECLARE_SOA_COLUMN(IsSelLcPiKPTOFPID, isSelLcPiKPTOFPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelLcPiKPTOFplusRICHPID, isSelLcPiKPTOFplusRICHPID, int); //!
} // namespace hf_sel_candidate_lc_alice3
DECLARE_SOA_TABLE(HfSelLcAlice3, "AOD", "HFSELLCA3B", //!
                  hf_sel_candidate_lc_alice3::IsSelLcPKPiNoPID,
                  hf_sel_candidate_lc_alice3::IsSelLcPKPiPerfectPID,
                  hf_sel_candidate_lc_alice3::IsSelLcPKPiTOFPID,
                  hf_sel_candidate_lc_alice3::IsSelLcPKPiTOFplusRICHPID,
                  hf_sel_candidate_lc_alice3::IsSelLcPiKPNoPID,
                  hf_sel_candidate_lc_alice3::IsSelLcPiKPPerfectPID,
                  hf_sel_candidate_lc_alice3::IsSelLcPiKPTOFPID,
                  hf_sel_candidate_lc_alice3::IsSelLcPiKPTOFplusRICHPID);

namespace hf_sel_candidate_lc_parametrized_pid
{
DECLARE_SOA_COLUMN(IsSelLcPKPiNoPID, isSelLcPKPiNoPID, int);           //!
DECLARE_SOA_COLUMN(IsSelLcPKPiPerfectPID, isSelLcPKPiPerfectPID, int); //!
DECLARE_SOA_COLUMN(IsSelLcPKPi, isSelLcPKPi, int);                     //!
DECLARE_SOA_COLUMN(IsSelLcPiKPNoPID, isSelLcPiKPNoPID, int);           //!
DECLARE_SOA_COLUMN(IsSelLcPiKPPerfectPID, isSelLcPiKPPerfectPID, int); //!
DECLARE_SOA_COLUMN(IsSelLcPiKP, isSelLcPiKP, int);                     //!
} // namespace hf_sel_candidate_lc_parametrized_pid
DECLARE_SOA_TABLE(HfSelLcParametrizedPid, "AOD", "HFSELLCP", //!
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcPKPiNoPID,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcPKPiPerfectPID,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcPKPi,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcPiKPNoPID,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcPiKPPerfectPID,
                  hf_sel_candidate_lc_parametrized_pid::IsSelLcPiKP);

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
DECLARE_SOA_COLUMN(IsSelLcK0sP, isSelLcK0sP, int);
} // namespace hf_sel_candidate_lc_to_k0s_p
DECLARE_SOA_TABLE(HfSelLcK0sP, "AOD", "HFSELLCK0SP", //!
                  hf_sel_candidate_lc_to_k0s_p::IsSelLcK0sP);

namespace hf_sel_candidate_b0
{
DECLARE_SOA_COLUMN(IsSelB0ToDPi, isSelB0ToDPi, int); //!
} // namespace hf_sel_candidate_b0
DECLARE_SOA_TABLE(HfSelB0ToDPi, "AOD", "HFSELB0", //!
                  hf_sel_candidate_b0::IsSelB0ToDPi);

namespace hf_sel_candidate_bplus
{
DECLARE_SOA_COLUMN(IsSelBPlusToD0Pi, isSelBPlusToD0Pi, int); //!
} // namespace hf_sel_candidate_bplus
DECLARE_SOA_TABLE(HfSelBplusToD0Pi, "AOD", "HFSELBPLUS", //!
                  hf_sel_candidate_bplus::IsSelBPlusToD0Pi);

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
} // namespace o2::aod
#endif // O2_ANALYSIS_CANDIDATESELECTIONTABLES_H_
