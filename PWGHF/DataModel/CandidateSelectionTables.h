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
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN

#ifndef PWGHF_DATAMODEL_CANDIDATESELECTIONTABLES_H_
#define PWGHF_DATAMODEL_CANDIDATESELECTIONTABLES_H_

#include <Framework/ASoA.h>

#include <vector>

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
DECLARE_SOA_COLUMN(MlProbD0, mlProbD0, std::vector<float>);       //!
DECLARE_SOA_COLUMN(MlProbD0bar, mlProbD0bar, std::vector<float>); //!
} // namespace hf_sel_candidate_d0

DECLARE_SOA_TABLE(HfSelD0, "AOD", "HFSELD0", //!
                  hf_sel_candidate_d0::IsSelD0,
                  hf_sel_candidate_d0::IsSelD0bar,
                  hf_sel_candidate_d0::IsRecoHfFlag,
                  hf_sel_candidate_d0::IsRecoTopol,
                  hf_sel_candidate_d0::IsRecoCand,
                  hf_sel_candidate_d0::IsRecoPid);

DECLARE_SOA_TABLE(HfMlD0, "AOD", "HFMLD0", //!
                  hf_sel_candidate_d0::MlProbD0,
                  hf_sel_candidate_d0::MlProbD0bar);

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
DECLARE_SOA_COLUMN(IsSelDplusToPiKPi, isSelDplusToPiKPi, int);                  //!
DECLARE_SOA_COLUMN(MlProbDplusToPiKPi, mlProbDplusToPiKPi, std::vector<float>); //!
} // namespace hf_sel_candidate_dplus

DECLARE_SOA_TABLE(HfSelDplusToPiKPi, "AOD", "HFSELDPLUS", //!
                  hf_sel_candidate_dplus::IsSelDplusToPiKPi);

DECLARE_SOA_TABLE(HfMlDplusToPiKPi, "AOD", "HFMLDPLUS", //!
                  hf_sel_candidate_dplus::MlProbDplusToPiKPi);

namespace hf_sel_candidate_ds
{
DECLARE_SOA_COLUMN(IsSelDsToKKPi, isSelDsToKKPi, int);                  //!
DECLARE_SOA_COLUMN(IsSelDsToPiKK, isSelDsToPiKK, int);                  //!
DECLARE_SOA_COLUMN(MlProbDsToKKPi, mlProbDsToKKPi, std::vector<float>); //!
DECLARE_SOA_COLUMN(MlProbDsToPiKK, mlProbDsToPiKK, std::vector<float>); //!
} // namespace hf_sel_candidate_ds

DECLARE_SOA_TABLE(HfSelDsToKKPi, "AOD", "HFSELDS", //!
                  hf_sel_candidate_ds::IsSelDsToKKPi, hf_sel_candidate_ds::IsSelDsToPiKK);

DECLARE_SOA_TABLE(HfMlDsToKKPi, "AOD", "HFMLDS", //!
                  hf_sel_candidate_ds::MlProbDsToKKPi, hf_sel_candidate_ds::MlProbDsToPiKK);

namespace hf_sel_candidate_dstar
{
DECLARE_SOA_COLUMN(IsSelDstarToD0Pi, isSelDstarToD0Pi, bool);                 //! checking if all four of following check pass
DECLARE_SOA_COLUMN(IsRecoD0Flag, isRecoD0Flag, bool);                         //! checking DecayType::D0ToPiK of D0prong
DECLARE_SOA_COLUMN(IsRecoTopol, isRecoTopol, bool);                           //! checking conjugate independent Topological selection on Dstar
DECLARE_SOA_COLUMN(IsRecoCand, isRecoCand, bool);                             //! checking conjugate dependent Topological selecton on Dstar
DECLARE_SOA_COLUMN(IsRecoPid, isRecoPid, bool);                               //! checking PID selection on daughters of D0Prong
DECLARE_SOA_COLUMN(MlProbDstarToD0Pi, mlProbDstarToD0Pi, std::vector<float>); //! ML probability for Dstar to D0Pi
} // namespace hf_sel_candidate_dstar

DECLARE_SOA_TABLE(HfSelDstarToD0Pi, "AOD", "HFSELDSTAR", //! Table stores information about selection flag on Dstar candidate
                  hf_sel_candidate_dstar::IsSelDstarToD0Pi,
                  hf_sel_candidate_dstar::IsRecoD0Flag,
                  hf_sel_candidate_dstar::IsRecoTopol,
                  hf_sel_candidate_dstar::IsRecoCand,
                  hf_sel_candidate_dstar::IsRecoPid);

DECLARE_SOA_TABLE(HfMlDstarToD0Pi, "AOD", "HFMLDSTAR", //! Table stores ML probability for Dstar to D0Pi
                  hf_sel_candidate_dstar::MlProbDstarToD0Pi);

namespace hf_sel_candidate_lc
{
DECLARE_SOA_COLUMN(IsSelLcToPKPi, isSelLcToPKPi, int);                  //!
DECLARE_SOA_COLUMN(IsSelLcToPiKP, isSelLcToPiKP, int);                  //!
DECLARE_SOA_COLUMN(MlProbLcToPKPi, mlProbLcToPKPi, std::vector<float>); //!
DECLARE_SOA_COLUMN(MlProbLcToPiKP, mlProbLcToPiKP, std::vector<float>); //!
} // namespace hf_sel_candidate_lc
DECLARE_SOA_TABLE(HfSelLc, "AOD", "HFSELLC", //!
                  hf_sel_candidate_lc::IsSelLcToPKPi, hf_sel_candidate_lc::IsSelLcToPiKP);
DECLARE_SOA_TABLE(HfMlLcToPKPi, "AOD", "HFMLLc", //!
                  hf_sel_candidate_lc::MlProbLcToPKPi, hf_sel_candidate_lc::MlProbLcToPiKP);

namespace hf_sel_candidate_cd
{
DECLARE_SOA_COLUMN(IsSelCdToDeKPi, isSelCdToDeKPi, int); //!
DECLARE_SOA_COLUMN(IsSelCdToPiKDe, isSelCdToPiKDe, int); //!
} // namespace hf_sel_candidate_cd
DECLARE_SOA_TABLE(HfSelCd, "AOD", "HFSELCD", //!
                  hf_sel_candidate_cd::IsSelCdToDeKPi, hf_sel_candidate_cd::IsSelCdToPiKDe);

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
DECLARE_SOA_COLUMN(MlProbLcToK0sP, mlProbLcToK0sP, std::vector<float>); //!
} // namespace hf_sel_candidate_lc_to_k0s_p
DECLARE_SOA_TABLE(HfSelLcToK0sP, "AOD", "HFSELLCK0SP", //!
                  hf_sel_candidate_lc_to_k0s_p::IsSelLcToK0sP);
DECLARE_SOA_TABLE(HfMlLcToK0sP, "AOD", "HFMLLcK0sP", //!
                  hf_sel_candidate_lc_to_k0s_p::MlProbLcToK0sP);

namespace hf_sel_candidate_b0
{
DECLARE_SOA_COLUMN(IsSelB0ToDPi, isSelB0ToDPi, int);     //! selection flag on B0 candidate
DECLARE_SOA_COLUMN(MlProbB0ToDPi, mlProbB0ToDPi, float); //! ML score of B0 candidate for signal class
} // namespace hf_sel_candidate_b0

DECLARE_SOA_TABLE(HfSelB0ToDPi, "AOD", "HFSELB0", //!
                  hf_sel_candidate_b0::IsSelB0ToDPi);

DECLARE_SOA_TABLE(HfMlB0ToDPi, "AOD", "HFMLB0", //!
                  hf_sel_candidate_b0::MlProbB0ToDPi);

namespace hf_sel_candidate_bs
{
DECLARE_SOA_COLUMN(IsSelBsToDsPi, isSelBsToDsPi, int);                        //!
DECLARE_SOA_COLUMN(MlProbBsToDsPi, mlProbBsToDsPi, std::vector<float>);       //!
DECLARE_SOA_COLUMN(MlProbBsToJpsiPhi, mlProbBsToJpsiPhi, std::vector<float>); //!
} // namespace hf_sel_candidate_bs

DECLARE_SOA_TABLE(HfSelBsToDsPi, "AOD", "HFSELBS", //!
                  hf_sel_candidate_bs::IsSelBsToDsPi);
DECLARE_SOA_TABLE(HfMlBsToDsPi, "AOD", "HFMLBS", //!
                  hf_sel_candidate_bs::MlProbBsToDsPi);
DECLARE_SOA_TABLE(HfMlBsToJpsiPhi, "AOD", "HFMLBSTOJPSIPHI", //!
                  hf_sel_candidate_bs::MlProbBsToDsPi);

namespace hf_sel_candidate_bplus
{
DECLARE_SOA_COLUMN(IsSelBplusToD0Pi, isSelBplusToD0Pi, int);       //! selection flag on B+ candidate
DECLARE_SOA_COLUMN(MlProbBplusToD0Pi, mlProbBplusToD0Pi, float);   //! ML score of B+ candidate for signal class
DECLARE_SOA_COLUMN(MlProbBplusToJpsiK, mlProbBplusToJpsiK, float); //! ML score of B+ candidate for signal class
} // namespace hf_sel_candidate_bplus

DECLARE_SOA_TABLE(HfSelBplusToD0Pi, "AOD", "HFSELBPLUS", //!
                  hf_sel_candidate_bplus::IsSelBplusToD0Pi);

DECLARE_SOA_TABLE(HfMlBplusToD0Pi, "AOD", "HFMLBPLUS", //!
                  hf_sel_candidate_bplus::MlProbBplusToD0Pi);

DECLARE_SOA_TABLE(HfMlBplusToJpsiK, "AOD", "HFMLBPLUSTOJPSIK", //!
                  hf_sel_candidate_bplus::MlProbBplusToJpsiK);
namespace hf_sel_candidate_lb
{
DECLARE_SOA_COLUMN(IsSelLbToLcPi, isSelLbToLcPi, int);     //! selection flag on Lb candidate
DECLARE_SOA_COLUMN(MlProbLbToLcPi, mlProbLbToLcPi, float); //! ML score of Lb candidate for signal class

} // namespace hf_sel_candidate_lb

DECLARE_SOA_TABLE(HfSelLbToLcPi, "AOD", "HFSELLB", //!
                  hf_sel_candidate_lb::IsSelLbToLcPi);
DECLARE_SOA_TABLE(HfMlLbToLcPi, "AOD", "HFMLLB", //!
                  hf_sel_candidate_lb::MlProbLbToLcPi);

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
// XicPlus to P K Pi
DECLARE_SOA_COLUMN(IsSelXicToPKPi, isSelXicToPKPi, int);                  //!
DECLARE_SOA_COLUMN(IsSelXicToPiKP, isSelXicToPiKP, int);                  //!
DECLARE_SOA_COLUMN(MlProbXicToPKPi, mlProbXicToPKPi, std::vector<float>); //!
DECLARE_SOA_COLUMN(MlProbXicToPiKP, mlProbXicToPiKP, std::vector<float>); //!
// XicPlus to Xi Pi Pi
DECLARE_SOA_COLUMN(IsSelXicToXiPiPi, isSelXicToXiPiPi, int);                  //!
DECLARE_SOA_COLUMN(MlProbXicToXiPiPi, mlProbXicToXiPiPi, std::vector<float>); //!
enum XicToXiPiPiSelectionStep {
  RecoTotal = 0,
  RecoKinTopol,
  RecoTrackQuality,
  RecoPID,
  RecoMl,
  NSelectionSteps
};
} // namespace hf_sel_candidate_xic

DECLARE_SOA_TABLE(HfSelXicToPKPi, "AOD", "HFSELXIC", //!
                  hf_sel_candidate_xic::IsSelXicToPKPi, hf_sel_candidate_xic::IsSelXicToPiKP);
DECLARE_SOA_TABLE(HfMlXicToPKPi, "AOD", "HFMLXIC", //!
                  hf_sel_candidate_xic::MlProbXicToPKPi, hf_sel_candidate_xic::MlProbXicToPiKP);
// XicPlus to Xi Pi Pi
DECLARE_SOA_TABLE(HfSelXicToXiPiPi, "AOD", "HFSELXICTOXI2PI", //!
                  hf_sel_candidate_xic::IsSelXicToXiPiPi);
DECLARE_SOA_TABLE(HfMlXicToXiPiPi, "AOD", "HFMLXICTOXI2PI", //!
                  hf_sel_candidate_xic::MlProbXicToXiPiPi);

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
DECLARE_SOA_COLUMN(StatusPidCharmBaryon, statusPidCharmBaryon, bool);
DECLARE_SOA_COLUMN(StatusInvMassLambda, statusInvMassLambda, bool);
DECLARE_SOA_COLUMN(StatusInvMassCascade, statusInvMassCascade, bool);
DECLARE_SOA_COLUMN(StatusInvMassCharmBaryon, statusInvMassCharmBaryon, bool);
DECLARE_SOA_COLUMN(ResultSelections, resultSelections, bool);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCharmBaryon, tpcNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCasc, tpcNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCharmBaryon, tofNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCasc, tofNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);
DECLARE_SOA_COLUMN(MlProbToXiPi, mlProbToXiPi, std::vector<float>);
} // namespace hf_sel_toxipi

DECLARE_SOA_TABLE(HfSelToXiPi, "AOD", "HFSELTOXIPI",
                  hf_sel_toxipi::StatusPidLambda, hf_sel_toxipi::StatusPidCascade, hf_sel_toxipi::StatusPidCharmBaryon,
                  hf_sel_toxipi::StatusInvMassLambda, hf_sel_toxipi::StatusInvMassCascade, hf_sel_toxipi::StatusInvMassCharmBaryon,
                  hf_sel_toxipi::ResultSelections, hf_sel_toxipi::PidTpcInfoStored, hf_sel_toxipi::PidTofInfoStored,
                  hf_sel_toxipi::TpcNSigmaPiFromCharmBaryon, hf_sel_toxipi::TpcNSigmaPiFromCasc, hf_sel_toxipi::TpcNSigmaPiFromLambda, hf_sel_toxipi::TpcNSigmaPrFromLambda,
                  hf_sel_toxipi::TofNSigmaPiFromCharmBaryon, hf_sel_toxipi::TofNSigmaPiFromCasc, hf_sel_toxipi::TofNSigmaPiFromLambda, hf_sel_toxipi::TofNSigmaPrFromLambda);

DECLARE_SOA_TABLE(HfSelToXiPiKf, "AOD", "HFSELTOXIPIKF",
                  hf_sel_toxipi::ResultSelections,
                  hf_sel_toxipi::TpcNSigmaPiFromCharmBaryon, hf_sel_toxipi::TpcNSigmaPiFromCasc, hf_sel_toxipi::TpcNSigmaPiFromLambda, hf_sel_toxipi::TpcNSigmaPrFromLambda,
                  hf_sel_toxipi::TofNSigmaPiFromCharmBaryon, hf_sel_toxipi::TofNSigmaPiFromCasc, hf_sel_toxipi::TofNSigmaPiFromLambda, hf_sel_toxipi::TofNSigmaPrFromLambda);

DECLARE_SOA_TABLE(HfMlToXiPi, "AOD", "HFMLSELTOXIPI",
                  hf_sel_toxipi::MlProbToXiPi);

namespace hf_sel_toomegapi
{
DECLARE_SOA_COLUMN(StatusPidLambda, statusPidLambda, bool);
DECLARE_SOA_COLUMN(StatusPidCascade, statusPidCascade, bool);
DECLARE_SOA_COLUMN(StatusPidCharmBaryon, statusPidCharmBaryon, bool);
DECLARE_SOA_COLUMN(StatusInvMassLambda, statusInvMassLambda, bool);
DECLARE_SOA_COLUMN(StatusInvMassCascade, statusInvMassCascade, bool);
DECLARE_SOA_COLUMN(StatusInvMassCharmBaryon, statusInvMassCharmBaryon, bool);
DECLARE_SOA_COLUMN(ResultSelections, resultSelections, bool);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCharmBaryon, tpcNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCharmBaryon, tpcNSigmaKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCasc, tpcNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCharmBaryon, tofNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCharmBaryon, tofNSigmaKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCasc, tofNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);
// Machine learning column for omegac0 to omega pi
DECLARE_SOA_COLUMN(MlProbOmegac, mlProbOmegac, std::vector<float>);
} // namespace hf_sel_toomegapi

DECLARE_SOA_TABLE(HfSelToOmegaPi, "AOD", "HFSELTOOMEPI",
                  hf_sel_toomegapi::StatusPidLambda, hf_sel_toomegapi::StatusPidCascade, hf_sel_toomegapi::StatusPidCharmBaryon,
                  hf_sel_toomegapi::StatusInvMassLambda, hf_sel_toomegapi::StatusInvMassCascade, hf_sel_toomegapi::StatusInvMassCharmBaryon,
                  hf_sel_toomegapi::ResultSelections, hf_sel_toomegapi::PidTpcInfoStored, hf_sel_toomegapi::PidTofInfoStored,
                  hf_sel_toomegapi::TpcNSigmaPiFromCharmBaryon, hf_sel_toomegapi::TpcNSigmaKaFromCasc, hf_sel_toomegapi::TpcNSigmaPiFromLambda, hf_sel_toomegapi::TpcNSigmaPrFromLambda,
                  hf_sel_toomegapi::TofNSigmaPiFromCharmBaryon, hf_sel_toomegapi::TofNSigmaKaFromCasc, hf_sel_toomegapi::TofNSigmaPiFromLambda, hf_sel_toomegapi::TofNSigmaPrFromLambda);

DECLARE_SOA_TABLE(HfSelToOmegaKaKf, "AOD", "HFSELTOOMEGAKAKF",
                  hf_sel_toomegapi::StatusPidLambda, hf_sel_toomegapi::StatusPidCascade, hf_sel_toomegapi::StatusPidCharmBaryon,
                  hf_sel_toomegapi::StatusInvMassLambda, hf_sel_toomegapi::StatusInvMassCascade, hf_sel_toomegapi::StatusInvMassCharmBaryon,
                  hf_sel_toomegapi::ResultSelections, hf_sel_toomegapi::PidTpcInfoStored, hf_sel_toomegapi::PidTofInfoStored,
                  hf_sel_toomegapi::TpcNSigmaKaFromCharmBaryon, hf_sel_toomegapi::TpcNSigmaKaFromCasc, hf_sel_toomegapi::TpcNSigmaPiFromLambda, hf_sel_toomegapi::TpcNSigmaPrFromLambda,
                  hf_sel_toomegapi::TofNSigmaKaFromCharmBaryon, hf_sel_toomegapi::TofNSigmaKaFromCasc, hf_sel_toomegapi::TofNSigmaPiFromLambda, hf_sel_toomegapi::TofNSigmaPrFromLambda);

DECLARE_SOA_TABLE(HfMlSelOmegacToOmegaPi, "AOD", "HFMLOMEGAC", //!
                  hf_sel_toomegapi::MlProbOmegac);
namespace hf_sel_toomegaka
{
DECLARE_SOA_COLUMN(StatusPidLambda, statusPidLambda, bool);
DECLARE_SOA_COLUMN(StatusPidCascade, statusPidCascade, bool);
DECLARE_SOA_COLUMN(StatusPidCharmBaryon, statusPidCharmBaryon, bool);
DECLARE_SOA_COLUMN(StatusInvMassLambda, statusInvMassLambda, bool);
DECLARE_SOA_COLUMN(StatusInvMassCascade, statusInvMassCascade, bool);
DECLARE_SOA_COLUMN(StatusInvMassCharmBaryon, statusInvMassCharmBaryon, bool);
DECLARE_SOA_COLUMN(ResultSelections, resultSelections, bool);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCharmBaryon, tpcNSigmaKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCasc, tpcNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCharmBaryon, tofNSigmaKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCasc, tofNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);
} // namespace hf_sel_toomegaka

DECLARE_SOA_TABLE(HfSelToOmegaKa, "AOD", "HFSELTOOMEKA",
                  hf_sel_toomegaka::StatusPidLambda, hf_sel_toomegaka::StatusPidCascade, hf_sel_toomegaka::StatusPidCharmBaryon,
                  hf_sel_toomegaka::StatusInvMassLambda, hf_sel_toomegaka::StatusInvMassCascade, hf_sel_toomegaka::StatusInvMassCharmBaryon,
                  hf_sel_toomegaka::ResultSelections, hf_sel_toomegaka::PidTpcInfoStored, hf_sel_toomegaka::PidTofInfoStored,
                  hf_sel_toomegaka::TpcNSigmaKaFromCharmBaryon, hf_sel_toomegaka::TpcNSigmaKaFromCasc, hf_sel_toomegaka::TpcNSigmaPiFromLambda, hf_sel_toomegaka::TpcNSigmaPrFromLambda,
                  hf_sel_toomegaka::TofNSigmaKaFromCharmBaryon, hf_sel_toomegaka::TofNSigmaKaFromCasc, hf_sel_toomegaka::TofNSigmaPiFromLambda, hf_sel_toomegaka::TofNSigmaPrFromLambda)

} // namespace o2::aod

#endif // PWGHF_DATAMODEL_CANDIDATESELECTIONTABLES_H_
