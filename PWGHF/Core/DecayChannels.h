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

/// \file DecayChannels.h
/// \brief Definitions of constants for MC flagging of HF decay channels.
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University
/// \note DecayChannelMain enums define unique combinations of the mother and the daughters for main channels.
/// \note DecayChannelResonant enums define unique combinations of the mother and the daughters for resonant channels.
/// \note Value 0 is reserved to indicate no match.
/// \note Daughter ordering convention: (charm|strange|π±|K±|π0), (baryon|meson), (+|−)

#ifndef PWGHF_CORE_DECAYCHANNELS_H_
#define PWGHF_CORE_DECAYCHANNELS_H_

#include <cstdint>

namespace o2::hf_decay
{

// TODO
// - HF cascades (Λc+ → p K0short)
// - HF cascades to LF cascades (Ωc0/Ξc0 → Ξ+ π−, Ξc+ → Ξ+ π− π+)
// - Σc

namespace hf_cand_2prong
{
/// @brief 2-prong candidates: main channels
enum DecayChannelMain : int8_t {
  // D0
  D0ToPiK = 1, // π+ K−
  D0ToPiKPi0,  // π+ K− π0
  D0ToPiPi,    // π+ π−
  D0ToPiPiPi0, // π+ π− π0
  D0ToKK,      // K+ K−
  //
  LastChannelMain
};
/// @brief 2-prong candidates: resonant channels
enum DecayChannelResonant : int8_t {
  // D0
  D0ToRhoplusPi = 1, // ρ+ π−
  D0ToRhoplusK,      // ρ+ K−
  D0ToKstar0Pi0,     // anti-K*0 π0
  D0ToKstarPi,       // K*− π+
  //
  LastChannelResonant
};
} // namespace hf_cand_2prong

namespace hf_cand_3prong
{
/// @brief 3-prong candidates: main channels
enum DecayChannelMain : int8_t {
  // D+
  DplusToPiKPi = 1, // π+ K− π+
  DplusToPiKPiPi0,  // π+ K− π+ π0
  DplusToPiPiPi,    // π+ π− π+
  DplusToPiKK,      // π+ K− K+
  // Ds+
  DsToPiKK,      // π+ K− K+
  DsToPiKKPi0,   // π+ K− K+ π0
  DsToPiPiK,     // π+ π− K+
  DsToPiPiPi,    // π+ π− π+
  DsToPiPiPiPi0, // π+ π− π+ π0
  // D*+
  DstarToPiKPi, // π+ K− π+ (from [(D0 → π+ K−) π+])
  // Λc+
  LcToPKPi,    // p K− π+
  LcToPKPiPi0, // p K− π+ π0
  LcToPPiPi,   // p π− π+
  LcToPKK,     // p K− K+
  // Ξc+
  XicToPKPi,  // p K− π+
  XicToPKK,   // p K− K+
  XicToSPiPi, // Σ+ π− π+
  //
  LastChannelMain
};
/// @brief 3-prong candidates: resonant channels
enum DecayChannelResonant : int8_t {
  // D+
  DplusToPhiPi = 1,    // φ π+
  DplusToKstar0K,      // anti-K*0 K+
  DplusToKstar1430_0K, // anti-K*0(1430) K+
  DplusToRho0Pi,       // ρ0 π+
  DplusToF2_1270Pi,    // f2(1270) π+
  // Ds+
  DsToPhiPi,      // φ π+
  DsToPhiRhoplus, // φ ρ+
  DsToKstar0K,    // anti-K*0 K+
  DsToKstar0Pi,   // anti-K*0 π+
  DsToRho0Pi,     // ρ0 π+
  DsToRho0K,      // ρ0 K+
  DsToF2_1270Pi,  // f2(1270) π+
  DsToF0_1370K,   // f0(1370) K+
  DsToEtaPi,      // η π+
  // Λc+
  LcToPKstar0,        // p K*0(892)
  LcToDeltaplusplusK, // Δ++ K−
  LcToL1520Pi,        // Λ(1520) π+
  // Ξc+
  XicToPKstar0, // p anti-K*0(892)
  XicToPPhi,    // p φ
  //
  LastChannelResonant
};
} // namespace hf_cand_3prong

namespace hf_cand_dstar
{
/// @brief D*+ candidates: main channels
enum DecayChannelMain : int8_t {
  // D*+
  DstarToPiKPi = 1, // π+ K− π+ (from [(D0 → π+ K−) π+])
  DstarToPiKPiPi0,  // π+ K− π+ π0 (from [(D0 → π+ K− π0) π+] or [(D+ → π+ K− π+) π0])
  //
  LastChannelMain
};
} // namespace hf_cand_dstar

namespace hf_cand_beauty
{
/// @brief beauty candidates: main channels
enum DecayChannelMain : int8_t {
  // B0
  B0ToDminusPi = 1,  // D− π+
  B0ToDminusPiPi0,   // D− π+ π0
  B0ToDminusPiGamma, // D− π+ γ0
  B0ToDminusK,       // D− K+
  B0ToD0PiPi,        // anti-D0 π+ π−
  // Bs0
  BsToDsPi,      // Ds− π+
  BsToDsPiPi0,   // Ds− π+ π0
  BsToDsPiGamma, // Ds− π+ γ0
  BsToDsK,       // Ds− K+
  // Λb0
  LbToLcPi,      // Λc+ π−
  LbToLcPiPi0,   // Λc+ π− π0
  LbToLcPiGamma, // Λc+ π− γ0
  LbToLcK,       // Λc+ K−
  LbToLcKPi0,    // Λc+ K− π0
  // B+
  BplusToD0Pi,      // anti-D0 π+
  BplusToD0PiPi0,   // anti-D0 π+ π0
  BplusToD0PiGamma, // anti-D0 π+ γ0
  BplusToD0K,       // anti-D0 K+
  //
  LastChannelMain
};
/// @brief beauty candidates: resonant channels
enum DecayChannelResonant : int8_t {
  // B0
  B0ToDminusRhoplus = 1, // D− ρ+
  B0ToDstarminusPi,      // D*− π+
  // Bs0
  BsToDsRhoplus, // Ds− ρ+
  BsToDsstarPi,  // Ds*− π+
  // Λb0
  LbToLcRhoplus, // Λc+ ρ−
  LbToScPi,      // Σc+ π−
  LbToScK,       // Σc+ K−
  LbToSc0Pi0,    // Σc0 π0
  // B+
  BplusToD0Rhoplus, // anti-D0 ρ+
  BplusToDstar0Pi,  // anti-D*0 π+
  //
  LastChannelResonant
};
/// @brief beauty candidates: beauty to J/ψ decay channels
enum DecayChannelToJpsiMain : int8_t {
  // B0
  B0ToJpsiPiK = 1, // J/ψ π- K+
  // Bs0
  BsToJpsiKK, // J/ψ K+ K-
  // Λb0
  LbToJpsiPK, // J/ψ p K-
  // B+
  BplusToJpsiK, // J/ψ K+
  // Bc+
  BcToJpsiPi, // J/ψ π+
  //
  LastChannelToJpsiMain
};
/// @brief beauty candidates: beauty to J/ψ resonant decay channels
enum DecayChannelToJpsiResonant : int8_t {
  // B0
  B0ToJpsiKstar0 = 1, // J/ψ K*0(892)
  // Bs0
  BsToJpsiPhi, // J/ψ φ
  //
  LastChannelToJpsiResonant
};
} // namespace hf_cand_beauty
} // namespace o2::hf_decay

#endif // PWGHF_CORE_DECAYCHANNELS_H_
