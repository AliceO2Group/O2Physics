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
  D0ToPiK = 1,     // π+ K−
  D0ToPiKPi0 = 2,  // π+ K− π0
  D0ToPiPi = 3,    // π+ π−
  D0ToPiPiPi0 = 4, // π+ π− π0
  D0ToKK = 5,      // K+ K−
  //
  NChannelsMain = D0ToKK // last channel
};
/// @brief 2-prong candidates: resonant channels
enum DecayChannelResonant : int8_t {
  // D0
  D0ToRhoplusPi = 1, // ρ+ π−
  D0ToRhoplusK = 2,  // ρ+ K−
  D0ToKstar0Pi0 = 3, // anti-K*0 π0
  D0ToKstarPi = 4,   // K*− π+
  //
  NChannelsResonant = D0ToKstarPi // last channel
};
} // namespace hf_cand_2prong

namespace hf_cand_3prong
{
/// @brief 3-prong candidates: main channels
enum DecayChannelMain : int8_t {
  // D+
  DplusToPiKPi = 1,    // π+ K− π+
  DplusToPiKPiPi0 = 2, // π+ K− π+ π0
  DplusToPiPiPi = 3,   // π+ π− π+
  DplusToPiKK = 4,     // π+ K− K+
  // Ds+
  DsToPiKK = 5,      // π+ K− K+
  DsToPiKKPi0 = 6,   // π+ K− K+ π0
  DsToPiPiK = 7,     // π+ π− K+
  DsToPiPiPi = 8,    // π+ π− π+
  DsToPiPiPiPi0 = 9, // π+ π− π+ π0
  // D*+
  DstarToPiKPi = 10,       // π+ K− π+ (from [(D0 → π+ K−) π+])
  DstarToPiKPiPi0 = 11,    // π+ K− π+ π0
  DstarToPiKPiPi0Pi0 = 12, // π+ K− π+ π0 π0
  DstarToPiKK = 13,        // π+ K− K+
  DstarToPiKKPi0 = 14,     // π+ K− K+ π0
  DstarToPiPiPi = 15,      // π+ π− π+
  DstarToPiPiPiPi0 = 16,   // π+ π− π+ π0
  // Λc+
  LcToPKPi = 17,    // p K− π+
  LcToPKPiPi0 = 18, // p K− π+ π0
  LcToPPiPi = 19,   // p π− π+
  LcToPKK = 20,     // p K− K+
  // Ξc+
  XicToPKPi = 21,  // p K− π+
  XicToPKK = 22,   // p K− K+
  XicToSPiPi = 23, // Σ+ π− π+
  //
  NChannelsMain = XicToSPiPi // last channel
};
/// @brief 3-prong candidates: resonant channels
enum DecayChannelResonant : int8_t {
  // D+
  DplusToPhiPi = 1,        // φ π+
  DplusToKstar0K = 2,      // anti-K*0 K+
  DplusToKstar1430_0K = 3, // anti-K*0(1430) K+
  DplusToRho0Pi = 4,       // ρ0 π+
  DplusToF2_1270Pi = 5,    // f2(1270) π+
  // Ds+
  DsToPhiPi = 6,      // φ π+
  DsToPhiRhoplus = 7, // φ ρ+
  DsToKstar0K = 8,    // anti-K*0 K+
  DsToKstar0Pi = 9,   // anti-K*0 π+
  DsToRho0Pi = 10,    // ρ0 π+
  DsToRho0K = 11,     // ρ0 K+
  DsToF2_1270Pi = 12, // f2(1270) π+
  DsToF0_1370K = 13,  // f0(1370) K+
  DsToEtaPi = 14,     // η π+
  // D*+
  DstarToD0ToRhoplusPi = 15,       // ρ+ π−
  DstarToD0ToRhoplusK = 16,        // ρ+ K−
  DstarToD0ToKstar0Pi0 = 17,       // anti-K*0 π0
  DstarToD0ToKstarPi = 18,         // K*− π+
  DstarToDplusToPhiPi = 19,        // φ π+
  DstarToDplusToKstar0K = 20,      // anti-K*0 K+
  DstarToDplusToKstar1430_0K = 21, // anti-K*0(1430) K+
  DstarToDplusToRho0Pi = 22,       // ρ0 π+
  DstarToDplusToF2_1270Pi = 23,    // f2(1270) π+
  // Λc+
  LcToPKstar0 = 24,        // p K*0(892)
  LcToDeltaplusplusK = 25, // Δ++ K−
  LcToL1520Pi = 26,        // Λ(1520) π+
  // Ξc+
  XicToPKstar0 = 27, // p anti-K*0(892)
  XicToPPhi = 28,    // p φ
  //
  NChannelsResonant = XicToPPhi // last channel
};
} // namespace hf_cand_3prong

namespace hf_cand_dstar
{
/// @brief D*+ candidates: main channels
enum DecayChannelMain : int8_t {
  // D*+
  DstarToPiKPi = 1,    // π+ K− π+ (from [(D0 → π+ K−) π+])
  DstarToPiKPiPi0 = 2, // π+ K− π+ π0 (from [(D0 → π+ K− π0) π+] or [(D+ → π+ K− π+) π0])
  //
  NChannelsMain = DstarToPiKPiPi0 // last channel
};
} // namespace hf_cand_dstar

namespace hf_cand_beauty
{
/// @brief beauty candidates: main channels
enum DecayChannelMain : int8_t {
  // B0
  B0ToDminusPi = 1,      // D− π+
  B0ToDminusPiPi0 = 2,   // D− π+ π0
  B0ToDminusPiGamma = 3, // D− π+ γ0
  B0ToDminusK = 4,       // D− K+
  B0ToD0PiPi = 5,        // anti-D0 π+
  // Bs0
  BsToDsPi = 6,      // Ds− π+
  BsToDsPiPi0 = 7,   // Ds− π+ π0
  BsToDsPiGamma = 8, // Ds− π+ γ0
  BsToDsK = 9,       // Ds− K+
  // Λb0
  LbToLcPi = 10,      // Λc+ π−
  LbToLcPiPi0 = 11,   // Λc+ π− π0
  LbToLcPiGamma = 12, // Λc+ π− γ0
  LbToLcK = 13,       // Λc+ K−
  LbToLcKPi0 = 14,    // Λc+ K− π0
  // B+
  BplusToD0Pi = 15,      // anti-D0 π+
  BplusToD0PiPi0 = 16,   // anti-D0 π+ π0
  BplusToD0PiGamma = 17, // anti-D0 π+ γ0
  BplusToD0K = 18,       // anti-D0 K+
  //
  NChannelsMain = BplusToD0K // last channel
};
/// @brief beauty candidates: resonant channels
enum DecayChannelResonant : int8_t {
  // B0
  B0ToDminusRhoplus = 1, // D− ρ+
  B0ToDstarminusPi = 2,  // D*− π+
  // Bs0
  BsToDsRhoplus = 3, // Ds− ρ+
  BsToDsstarPi = 4,  // Ds*− π+
  // Λb0
  LbToLcRhoplus = 5, // Λc+ ρ−
  LbToScPi = 6,      // Σc+ π−
  LbToScK = 7,       // Σc+ K−
  LbToSc0Pi0 = 8,    // Σc0 π0
  // B+
  BplusToD0Rhoplus = 9, // anti-D0 ρ+
  BplusToDstar0Pi = 10, // anti-D*0 π+
  //
  NChannelsResonant = BplusToDstar0Pi // last channel
};
/// @brief beauty candidates: beauty to J/ψ decay channels
enum DecayChannelToJpsiMain : int8_t {
  // B0
  B0ToJpsiPiK = 1, // J/ψ π− K+
  // Bs0
  BsToJpsiKK = 2, // J/ψ K+ K−
  // Λb0
  LbToJpsiPK = 3, // J/ψ p K−
  // B+
  BplusToJpsiK = 4, // J/ψ K+
  // Bc+
  BcToJpsiPi = 5, // J/ψ π+
  //
  NChannelsToJpsiMain = BcToJpsiPi // last channel
};
/// @brief beauty candidates: beauty to J/ψ resonant decay channels
enum DecayChannelToJpsiResonant : int8_t {
  // B0
  B0ToJpsiKstar0 = 1, // J/ψ K*0(892)
  // Bs0
  BsToJpsiPhi = 2, // J/ψ φ
  //
  NChannelsToJpsiResonant = BsToJpsiPhi // last channel
};
} // namespace hf_cand_beauty
} // namespace o2::hf_decay

#endif // PWGHF_CORE_DECAYCHANNELS_H_
