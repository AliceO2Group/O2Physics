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

#ifndef O2_ANALYSIS_HFCANDIDATESELECTIONTABLES_H_
#define O2_ANALYSIS_HFCANDIDATESELECTIONTABLES_H_

namespace o2::aod
{
// selection steps
enum SelectionStep {
  RecoSkims = 0,
  RecoTopol,
  RecoPID,
  NSelectionSteps
};

namespace hf_selcandidate_d0
{
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);           //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int);     //!
DECLARE_SOA_COLUMN(IsRecoHFFlag, isRecoHFFlag, int); //!
DECLARE_SOA_COLUMN(IsRecoTopol, isRecoTopol, int);   //!
DECLARE_SOA_COLUMN(IsRecoCand, isRecoCand, int);     //!
DECLARE_SOA_COLUMN(IsRecoPID, isRecoPID, int);
} // namespace hf_selcandidate_d0
DECLARE_SOA_TABLE(HFSelD0Candidate, "AOD", "HFSELD0CAND", //!
                  hf_selcandidate_d0::IsSelD0,
                  hf_selcandidate_d0::IsSelD0bar,
                  hf_selcandidate_d0::IsRecoHFFlag,
                  hf_selcandidate_d0::IsRecoTopol,
                  hf_selcandidate_d0::IsRecoCand,
                  hf_selcandidate_d0::IsRecoPID);

namespace hf_selcandidate_d0_parametrizedPID
{
DECLARE_SOA_COLUMN(IsSelD0NoPID, isSelD0NoPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelD0PerfectPID, isSelD0PerfectPID, int);       //!
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);                           //!
DECLARE_SOA_COLUMN(IsSelD0barNoPID, isSelD0barNoPID, int);           //!
DECLARE_SOA_COLUMN(IsSelD0barPerfectPID, isSelD0barPerfectPID, int); //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int);                     //!
} // namespace hf_selcandidate_d0_parametrizedPID
DECLARE_SOA_TABLE(HFSelD0CandidateparametrizedPID, "AOD", "HFSELD0CANDP", //!
                  hf_selcandidate_d0_parametrizedPID::IsSelD0NoPID,
                  hf_selcandidate_d0_parametrizedPID::IsSelD0PerfectPID,
                  hf_selcandidate_d0_parametrizedPID::IsSelD0,
                  hf_selcandidate_d0_parametrizedPID::IsSelD0barNoPID,
                  hf_selcandidate_d0_parametrizedPID::IsSelD0barPerfectPID,
                  hf_selcandidate_d0_parametrizedPID::IsSelD0bar);

namespace hf_selcandidate_d0_ALICE3_Barrel
{
DECLARE_SOA_COLUMN(IsSelHFFlag, isSelHFFlag, int);                     //!
DECLARE_SOA_COLUMN(IsSelD0NoPID, isSelD0NoPID, int);                   //!
DECLARE_SOA_COLUMN(IsSelD0PerfectPID, isSelD0PerfectPID, int);         //!
DECLARE_SOA_COLUMN(IsSelD0TOFPID, isSelD0TOFPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelD0RICHPID, isSelD0RICHPID, int);               //!
DECLARE_SOA_COLUMN(IsSelD0TOFplusRICHPID, isSelD0TOFplusRICHPID, int); //!
DECLARE_SOA_COLUMN(IsSelD0barTOFplusRICHPID, isSelD0barTOFplusRICHPID, int); //!
} // namespace hf_selcandidate_d0_ALICE3_Barrel
DECLARE_SOA_TABLE(HFSelD0CandidateALICE3Barrel, "AOD", "HFSELD0CANDA3B", //!
                  hf_selcandidate_d0_ALICE3_Barrel::IsSelHFFlag,
                  hf_selcandidate_d0_ALICE3_Barrel::IsSelD0NoPID,
                  hf_selcandidate_d0_ALICE3_Barrel::IsSelD0PerfectPID,
                  hf_selcandidate_d0_ALICE3_Barrel::IsSelD0TOFPID,
                  hf_selcandidate_d0_ALICE3_Barrel::IsSelD0RICHPID,
                  hf_selcandidate_d0_ALICE3_Barrel::IsSelD0TOFplusRICHPID,
                  hf_selcandidate_d0_ALICE3_Barrel::IsSelD0barTOFplusRICHPID);

namespace hf_selcandidate_d0_ALICE3_Forward
{
DECLARE_SOA_COLUMN(IsSelHFFFlag, isSelHFFFlag, int);       //!
DECLARE_SOA_COLUMN(IsSelD0FNoPID, isSelD0FNoPID, int);     //!
DECLARE_SOA_COLUMN(IsSelD0FRICHPID, isSelD0FRICHPID, int); //!
} // namespace hf_selcandidate_d0_ALICE3_Forward
DECLARE_SOA_TABLE(HFSelD0CandidateALICE3Forward, "AOD", "HFSELD0CANDA3F", //!
                  hf_selcandidate_d0_ALICE3_Forward::IsSelHFFFlag,
                  hf_selcandidate_d0_ALICE3_Forward::IsSelD0FNoPID,
                  hf_selcandidate_d0_ALICE3_Forward::IsSelD0FRICHPID);

namespace hf_selcandidate_dplus
{
DECLARE_SOA_COLUMN(IsSelDplusToPiKPi, isSelDplusToPiKPi, int); //!
} // namespace hf_selcandidate_dplus
DECLARE_SOA_TABLE(HFSelDplusToPiKPiCandidate, "AOD", "HFSELDPLUSCAND", //!
                  hf_selcandidate_dplus::IsSelDplusToPiKPi);

namespace hf_selcandidate_ds
{
DECLARE_SOA_COLUMN(IsSelDsToKKPi, isSelDsToKKPi, int); //!
DECLARE_SOA_COLUMN(IsSelDsToPiKK, isSelDsToPiKK, int); //!
} // namespace hf_selcandidate_ds
DECLARE_SOA_TABLE(HFSelDsToKKPiCandidate, "AOD", "HFSELDSCAND", //!
                  hf_selcandidate_ds::IsSelDsToKKPi, hf_selcandidate_ds::IsSelDsToPiKK);

namespace hf_selcandidate_lc
{
DECLARE_SOA_COLUMN(IsSelLcpKpi, isSelLcpKpi, int); //!
DECLARE_SOA_COLUMN(IsSelLcpiKp, isSelLcpiKp, int); //!
} // namespace hf_selcandidate_lc
DECLARE_SOA_TABLE(HFSelLcCandidate, "AOD", "HFSELLCCAND", //!
                  hf_selcandidate_lc::IsSelLcpKpi, hf_selcandidate_lc::IsSelLcpiKp);

namespace hf_selcandidate_lc_ALICE3
{
DECLARE_SOA_COLUMN(IsSelLcPKPiNoPID, isSelLcPKPiNoPID, int);                   //!
DECLARE_SOA_COLUMN(IsSelLcPKPiPerfectPID, isSelLcPKPiPerfectPID, int);         //!
DECLARE_SOA_COLUMN(IsSelLcPKPiTOFPID, isSelLcPKPiTOFPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelLcPKPiTOFplusRICHPID, isSelLcPKPiTOFplusRICHPID, int); //!
DECLARE_SOA_COLUMN(IsSelLcPiKPNoPID, isSelLcPiKPNoPID, int);                   //!
DECLARE_SOA_COLUMN(IsSelLcPiKPPerfectPID, isSelLcPiKPPerfectPID, int);         //!
DECLARE_SOA_COLUMN(IsSelLcPiKPTOFPID, isSelLcPiKPTOFPID, int);                 //!
DECLARE_SOA_COLUMN(IsSelLcPiKPTOFplusRICHPID, isSelLcPiKPTOFplusRICHPID, int); //!
} // namespace hf_selcandidate_lc_ALICE3
DECLARE_SOA_TABLE(HFSelLcCandidateALICE3, "AOD", "HFSELLCCANDA3B", //!
                  hf_selcandidate_lc_ALICE3::IsSelLcPKPiNoPID,
                  hf_selcandidate_lc_ALICE3::IsSelLcPKPiPerfectPID,
                  hf_selcandidate_lc_ALICE3::IsSelLcPKPiTOFPID,
                  hf_selcandidate_lc_ALICE3::IsSelLcPKPiTOFplusRICHPID,
                  hf_selcandidate_lc_ALICE3::IsSelLcPiKPNoPID,
                  hf_selcandidate_lc_ALICE3::IsSelLcPiKPPerfectPID,
                  hf_selcandidate_lc_ALICE3::IsSelLcPiKPTOFPID,
                  hf_selcandidate_lc_ALICE3::IsSelLcPiKPTOFplusRICHPID);

namespace hf_selcandidate_lc_parametrizedPID
{
DECLARE_SOA_COLUMN(IsSelLcPKPiNoPID, isSelLcPKPiNoPID, int);           //!
DECLARE_SOA_COLUMN(IsSelLcPKPiPerfectPID, isSelLcPKPiPerfectPID, int); //!
DECLARE_SOA_COLUMN(IsSelLcPKPi, isSelLcPKPi, int);                     //!
DECLARE_SOA_COLUMN(IsSelLcPiKPNoPID, isSelLcPiKPNoPID, int);           //!
DECLARE_SOA_COLUMN(IsSelLcPiKPPerfectPID, isSelLcPiKPPerfectPID, int); //!
DECLARE_SOA_COLUMN(IsSelLcPiKP, isSelLcPiKP, int);                     //!
} // namespace hf_selcandidate_lc_parametrizedPID
DECLARE_SOA_TABLE(HFSelLcCandidateparametrizedPID, "AOD", "HFSELLCCANDP", //!
                  hf_selcandidate_lc_parametrizedPID::IsSelLcPKPiNoPID,
                  hf_selcandidate_lc_parametrizedPID::IsSelLcPKPiPerfectPID,
                  hf_selcandidate_lc_parametrizedPID::IsSelLcPKPi,
                  hf_selcandidate_lc_parametrizedPID::IsSelLcPiKPNoPID,
                  hf_selcandidate_lc_parametrizedPID::IsSelLcPiKPPerfectPID,
                  hf_selcandidate_lc_parametrizedPID::IsSelLcPiKP);

namespace hf_selcandidate_jpsi
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
} // namespace hf_selcandidate_jpsi
DECLARE_SOA_TABLE(HFSelJpsiCandidate, "AOD", "HFSELJPSICAND", //!
                  hf_selcandidate_jpsi::IsSelJpsiToEE,
                  hf_selcandidate_jpsi::IsSelJpsiToMuMu,
                  hf_selcandidate_jpsi::IsSelJpsiToEETopol,
                  hf_selcandidate_jpsi::IsSelJpsiToEETpc,
                  hf_selcandidate_jpsi::IsSelJpsiToEETof,
                  hf_selcandidate_jpsi::IsSelJpsiToEERich,
                  hf_selcandidate_jpsi::IsSelJpsiToEETofRich,
                  hf_selcandidate_jpsi::IsSelJpsiToMuMuTopol,
                  hf_selcandidate_jpsi::IsSelJpsiToMuMuMid);

namespace hf_selcandidate_lc_k0sp
{
DECLARE_SOA_COLUMN(IsSelLcK0sP, isSelLcK0sP, int);
} // namespace hf_selcandidate_lc_k0sp
DECLARE_SOA_TABLE(HFSelLcK0sPCandidate, "AOD", "HFSELLCK0SPCAND", //!
                  hf_selcandidate_lc_k0sp::IsSelLcK0sP);

namespace hf_selcandidate_bplus
{
DECLARE_SOA_COLUMN(IsSelBPlusToD0Pi, isSelBPlusToD0Pi, int); //!
} // namespace hf_selcandidate_bplus
DECLARE_SOA_TABLE(HFSelBPlusToD0PiCandidate, "AOD", "HFSELBPLUSCAND", //!
                  hf_selcandidate_bplus::IsSelBPlusToD0Pi);

namespace hf_selcandidate_lb
{
DECLARE_SOA_COLUMN(IsSelLbToLcPi, isSelLbToLcPi, int); //!
} // namespace hf_selcandidate_lb
DECLARE_SOA_TABLE(HFSelLbToLcPiCandidate, "AOD", "HFSELLBCAND", //!
                  hf_selcandidate_lb::IsSelLbToLcPi);

namespace hf_selcandidate_x
{
DECLARE_SOA_COLUMN(IsSelXToJpsiToEEPiPi, isSelXToJpsiToEEPiPi, int);     //!
DECLARE_SOA_COLUMN(IsSelXToJpsiToMuMuPiPi, isSelXToJpsiToMuMuPiPi, int); //!
} // namespace hf_selcandidate_x
DECLARE_SOA_TABLE(HFSelXToJpsiPiPiCandidate, "AOD", "HFSELXCAND", //!
                  hf_selcandidate_x::IsSelXToJpsiToEEPiPi, hf_selcandidate_x::IsSelXToJpsiToMuMuPiPi);

namespace hf_selcandidate_chic
{
DECLARE_SOA_COLUMN(IsSelChicToJpsiToEEGamma, isSelChicToJpsiToEEGamma, int);     //!
DECLARE_SOA_COLUMN(IsSelChicToJpsiToMuMuGamma, isSelChicToJpsiToMuMuGamma, int); //!
} // namespace hf_selcandidate_chic
DECLARE_SOA_TABLE(HFSelChicToJpsiGammaCandidate, "AOD", "HFSELCHICCAND", //!
                  hf_selcandidate_chic::IsSelChicToJpsiToEEGamma, hf_selcandidate_chic::IsSelChicToJpsiToMuMuGamma);

namespace hf_selcandidate_xic
{
DECLARE_SOA_COLUMN(IsSelXicToPKPi, isSelXicToPKPi, int); //!
DECLARE_SOA_COLUMN(IsSelXicToPiKP, isSelXicToPiKP, int); //!
} // namespace hf_selcandidate_xic
DECLARE_SOA_TABLE(HFSelXicToPKPiCandidate, "AOD", "HFSELXICCAND", //!
                  hf_selcandidate_xic::IsSelXicToPKPi, hf_selcandidate_xic::IsSelXicToPiKP);

namespace hf_selcandidate_xicc
{
DECLARE_SOA_COLUMN(IsSelXiccToPKPiPi, isSelXiccToPKPiPi, int); //!
} // namespace hf_selcandidate_xicc
DECLARE_SOA_TABLE(HFSelXiccToPKPiPiCandidate, "AOD", "HFSELXICCCAND", //!
                  hf_selcandidate_xicc::IsSelXiccToPKPiPi);
} // namespace o2::aod
#endif // O2_ANALYSIS_HFCANDIDATESELECTIONTABLES_H_
