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
/// \brief Definitions of constants for MC flagging of decay channels.
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University
/// \note DecayChannel enums define unique combinations of the mother particle and final daughters.
/// \note DecayChannelResonant enums define unique combinations of the mother particle and intermediate daughters.
/// \note Daughter ordering convention: (charm|strange|π±|K±|π0), (baryon|meson), (+|−)

#ifndef PWGHF_CORE_DECAYCHANNELS_H_
#define PWGHF_CORE_DECAYCHANNELS_H_

#include <cstdint>

namespace o2::hf_decay
{

// TODO
// - HF cascades (Λc+ → p K0short)
// - HF cascades to LF cascades (Ωc0/Ξc0 → Ξ+ π−, Ξc+ → Ξ+ π− π+)

namespace hf_cand_2prong
{
/// @brief 2-prong candidates
enum DecayChannel : int8_t {
  // D0
  D0ToPiK = 1, // π+ K−
  D0ToPiKPi0,  // π+ K− π0
  D0ToPiPi,    // π+ π−
  D0ToPiPiPi0, // π+ π− π0
  D0ToKK,      // K+ K−
  //
  LastChannel
};
} // namespace hf_cand_2prong

namespace hf_cand_3prong
{
/// @brief 3-prong candidates
enum DecayChannel : int8_t {
  // D+
  DplusToPiKPi = 1, // π+ K− π+
  DplusToPiKPiPi0,  // π+ K− π+ π0
  DplusToPiPiPi,    // π+ π− π+
  DplusToPiKK,      // π+ K− K+
  // Ds+
  DsToPiKK,      // π+ K− K+
  DsToPiKKPi0,   // π+ K− K+ π0
  DsToPiPiK,     // π+ π− K+
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
  LastChannel
};
/// @brief 3-prong candidates resonant
enum DecayChannelResonant : int8_t {
  // D+
  DplusToPhiPi,   // φ π+
  DplusToK0starK, // K0* K+
  // Ds+
  DsToPhiPi,   // φ π+
  DsToK0starK, // K0* K+
  //
  LastChannelResonant
};
} // namespace hf_cand_3prong

namespace hf_cand_dstar
{
/// @brief D*+ candidates
enum DecayChannel : int8_t {
  // D*+
  DstarToPiKPi = 1, // π+ K− π+ (from [(D0 → π+ K−) π+])
  DstarToPiKPiPi0,  // π+ K− π+ π0 (from [(D0 → π+ K− π0) π+] or [(D+ → π+ K− π+) π0])
  //
  LastChannel
};
} // namespace hf_cand_dstar

namespace hf_cand_beauty_neutral
{
/// @brief neutral beauty candidates
enum DecayChannel : int8_t {
  // B0
  B0ToDminusPi = 1, // D− π+
  B0ToDminusK,      // D− K+
  B0ToD0PiPi,       // anti-D0 π+ π−
  // Bs0
  BsToDsPi, // Ds− π+
  BsToDsK,  // Ds− K+
  // Λb0
  LbToLcPi, // Λc+ π−
  LbToLcK,  // Λc+ K−
  LbToScPi, // Σc+ π−
  LbToScK,  // Σc+ K−
  //
  LastChannel
};
enum DecayChannelResonant : int8_t {
  // B0
  B0ToDminusRhoplus = 1, // D− ρ+
  B0ToDstarminusPi,      // D*− π+
  // Bs0
  BsToDsRhoplus, // Ds− ρ+
  BsToDsstarPi,  // Ds*− π+
  // Λb0
  LbToLcRhoplus, // Λc+ ρ−
  //
  LastChannelResonant
};
} // namespace hf_cand_beauty_neutral

namespace hf_cand_bplus
{
/// @brief B+ candidates
enum DecayChannel : int8_t {
  // B+
  BplusToD0Pi = 1, // anti-D0 π+
  BplusToD0K,      // anti-D0 K+
  //
  LastChannel
};
/// @brief B+ candidates resonant
enum DecayChannelResonant : int8_t {
  // B+
  BplusToD0Rhoplus = 1, // anti-D0 ρ+
  BplusToD0starPi,      // anti-D0* π+
  //
  LastChannelResonant
};
} // namespace hf_cand_bplus
} // namespace o2::hf_decay

#endif // PWGHF_CORE_DECAYCHANNELS_H_
