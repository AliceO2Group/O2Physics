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
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN
/// \note DecayChannelMain enums define unique combinations of the mother and the daughters for main channels.
/// \note DecayChannelResonant enums define unique combinations of the mother and the daughters for resonant channels.
/// \note Value 0 is reserved to indicate no match.
/// \note Daughter ordering convention: (charm|strange|π±|K±|π0), (baryon|meson), (+|−)

#ifndef PWGHF_CORE_DECAYCHANNELS_H_
#define PWGHF_CORE_DECAYCHANNELS_H_

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

using namespace o2::constants::physics;

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

std::unordered_map<DecayChannelMain, std::vector<int>> finalStates2Prongs =
  {
    {DecayChannelMain::D0ToPiK, std::vector<int>{+kKMinus, +kPiPlus}},
    {DecayChannelMain::D0ToPiKPi0, std::vector<int>{+kKMinus, +kPiPlus, +kPi0}},
    {DecayChannelMain::D0ToPiPi, std::vector<int>{+kPiMinus, +kPiPlus}},
    {DecayChannelMain::D0ToPiPiPi0, std::vector<int>{+kPiMinus, +kPiPlus, +kPi0}},
    {DecayChannelMain::D0ToKK, std::vector<int>{+kKMinus, +kKPlus}},
};

std::unordered_map<DecayChannelResonant, std::array<int, 2>> resoStatesD0 =
  {
    {DecayChannelResonant::D0ToRhoplusPi, std::array<int, 2>{213, +kPiMinus}},
    {DecayChannelResonant::D0ToRhoplusK, std::array<int, 2>{213, +kKMinus}},
    {DecayChannelResonant::D0ToKstar0Pi0, std::array<int, 2>{-kK0Star892, +kPi0}},
    {DecayChannelResonant::D0ToKstarPi, std::array<int, 2>{-kKPlusStar892, +kPiPlus}},
};

} // namespace hf_cand_2prong

namespace hf_cand_3prong
{
/// @brief 3-prong candidates: main channels
enum DecayChannelMain : int8_t {
  // D+
  DplusToPiKPi = 1, // 1, π+ K− π+
  DplusToPiKPiPi0,  // 2, π+ K− π+ π0
  DplusToPiPiPi,    // 3, π+ π− π+
  DplusToPiKK,      // 4, π+ K− K+
  // Ds+
  DsToPiKK,      // 5, π+ K− K+
  DsToPiKKPi0,   // 6, π+ K− K+ π0
  DsToPiPiK,     // 7, π+ π− K+
  DsToPiPiPi,    // 8, π+ π− π+
  DsToPiPiPiPi0, // 9, π+ π− π+ π0
  // D*+
  DstarToPiKPi,       // 10, π+ K− π+ (from [(D0 → π+ K−) π+])
  DstarToPiKPiPi0,    // 11, π+ K− π+ π0
  DstarToPiKPiPi0Pi0, // 12, π+ K− π+ π0 π0
  DstarToPiKK,        // 13, π+ K− K+
  DstarToPiKKPi0,     // 14, π+ K− K+ π0
  DstarToPiPiPi,      // 15, π+ π− π+
  DstarToPiPiPiPi0,   // 16, π+ π− π+ π0
  // Λc+
  LcToPKPi,    // 17, p K− π+
  LcToPKPiPi0, // 18, p K− π+ π0
  LcToPPiPi,   // 19, p π− π+
  LcToPKK,     // 20, p K− K+
  // Ξc+
  XicToPKPi,  // 21, p K− π+
  XicToPKK,   // 22, p K− K+
  XicToSPiPi, // 23, Σ+ π− π+
  //
  LastChannelMain
};
/// @brief 3-prong candidates: resonant channels
enum DecayChannelResonant : int8_t {
  // D+
  DplusToPhiPi = 1,    // 1, φ π+
  DplusToKstar0K,      // 2, anti-K*0 K+
  DplusToKstar1430_0K, // 3, anti-K*0(1430) K+
  DplusToRho0Pi,       // 4, ρ0 π+
  DplusToF2_1270Pi,    // 5, f2(1270) π+
  // Ds+
  DsToPhiPi,      // 6, φ π+
  DsToPhiRhoplus, // 7, φ ρ+
  DsToKstar0K,    // 8, anti-K*0 K+
  DsToKstar0Pi,   // 9, anti-K*0 π+
  DsToRho0Pi,     // 10, ρ0 π+
  DsToRho0K,      // 11, ρ0 K+
  DsToF2_1270Pi,  // 12, f2(1270) π+
  DsToF0_1370K,   // 13, f0(1370) K+
  DsToEtaPi,      // 14, η π+
  // D*+
  DstarToD0ToRhoplusPi,       // 15, ρ+ π−
  DstarToD0ToRhoplusK,        // 16, ρ+ K−
  DstarToD0ToKstar0Pi0,       // 17, anti-K*0 π0
  DstarToD0ToKstarPi,         // 18, K*− π+
  DstarToDplusToPhiPi,        // 19, φ π+
  DstarToDplusToKstar0K,      // 20, anti-K*0 K+
  DstarToDplusToKstar1430_0K, // 21, anti-K*0(1430) K+
  DstarToDplusToRho0Pi,       // 22, ρ0 π+
  DstarToDplusToF2_1270Pi,    // 23, f2(1270) π+
  // Λc+
  LcToPKstar0,        // 24, p K*0(892)
  LcToDeltaplusplusK, // 25, Δ++ K−
  LcToL1520Pi,        // 26, Λ(1520) π+
  // Ξc+
  XicToPKstar0, // 27, p anti-K*0(892)
  XicToPPhi,    // 28, p φ
  //
  LastChannelResonant
};

std::unordered_map<DecayChannelResonant, std::array<int, 2>> resoStatesDPlus =
  {
    {DecayChannelResonant::DplusToPhiPi, std::array<int, 2>{+kPhi, +kPiPlus}},
    {DecayChannelResonant::DplusToKstar0K, std::array<int, 2>{-kK0Star892, +kKPlus}},
    {DecayChannelResonant::DplusToKstar1430_0K, std::array<int, 2>{+10311, +kKPlus}},
    {DecayChannelResonant::DplusToRho0Pi, std::array<int, 2>{+113, +kPiPlus}},
    {DecayChannelResonant::DplusToF2_1270Pi, std::array<int, 2>{+225, +kPiPlus}},
};

std::unordered_map<DecayChannelMain, std::vector<int>> finalStatesDPlus =
  {
    {DecayChannelMain::DplusToPiKPi, std::vector<int>{+kKMinus, +kKPlus, +kPiPlus}},
    {DecayChannelMain::DplusToPiKK, std::vector<int>{+kKMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DplusToPiKPiPi0, std::vector<int>{+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DplusToPiPiPi, std::vector<int>{+kPiMinus, +kPiPlus, +kPiPlus}},
};

// Ds± → K± K∓ π±
std::unordered_map<DecayChannelResonant, std::array<int, 2>> resoStatesDs =
  {
    {DecayChannelResonant::DsToPhiPi, std::array<int, 2>{+kPhi, +kPiPlus}},
    {DecayChannelResonant::DsToPhiRhoplus, std::array<int, 2>{+kPhi, +213}},
    {DecayChannelResonant::DsToKstar0K, std::array<int, 2>{-kK0Star892, +kKPlus}},
    {DecayChannelResonant::DsToKstar0Pi, std::array<int, 2>{+kK0Star892, +kPiPlus}},
    {DecayChannelResonant::DsToRho0Pi, std::array<int, 2>{113, +kPiPlus}},
    {DecayChannelResonant::DsToRho0K, std::array<int, 2>{113, +kKPlus}},
    {DecayChannelResonant::DsToF2_1270Pi, std::array<int, 2>{225, +kPiPlus}},
    {DecayChannelResonant::DsToF0_1370K, std::array<int, 2>{10221, +kKPlus}},
    {DecayChannelResonant::DsToEtaPi, std::array<int, 2>{221, +kPiPlus}},
};

std::unordered_map<DecayChannelMain, std::vector<int>> finalStatesDs =
  {
    {DecayChannelMain::DsToPiKK, std::vector<int>{+kKMinus, +kKPlus, +kPiPlus}},
    {DecayChannelMain::DsToPiKKPi0, std::vector<int>{+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DsToPiPiK, std::vector<int>{+kKPlus, +kPiPlus, +kPiMinus}},
    {DecayChannelMain::DsToPiPiPi, std::vector<int>{+kPiMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DsToPiPiPiPi0, std::vector<int>{+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

// Dstar → K± K∓ π±
std::unordered_map<DecayChannelResonant, std::array<int, 2>> resoStatesDstar =
  {
    {DecayChannelResonant::DstarToD0ToRhoplusPi, std::array<int, 2>{213, +kPiMinus}},
    {DecayChannelResonant::DstarToD0ToRhoplusK, std::array<int, 2>{213, +kKMinus}},
    {DecayChannelResonant::DstarToD0ToKstar0Pi0, std::array<int, 2>{-kK0Star892, +kPi0}},
    {DecayChannelResonant::DstarToD0ToKstarPi, std::array<int, 2>{-kKPlusStar892, +kPiPlus}},
    {DecayChannelResonant::DstarToDplusToPhiPi, std::array<int, 2>{+kPhi, +kPiPlus}},
    {DecayChannelResonant::DstarToDplusToKstar0K, std::array<int, 2>{-kK0Star892, +kKPlus}},
    {DecayChannelResonant::DstarToDplusToKstar1430_0K, std::array<int, 2>{+10311, +kKPlus}},
    {DecayChannelResonant::DstarToDplusToRho0Pi, std::array<int, 2>{+113, +kPiPlus}},
    {DecayChannelResonant::DstarToDplusToF2_1270Pi, std::array<int, 2>{+225, +kPiPlus}},
};

std::unordered_map<DecayChannelMain, std::vector<int>> finalStatesDstar =
  {
    {DecayChannelMain::DstarToPiKPi, std::vector<int>{+kKMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DstarToPiKPiPi0, std::vector<int>{+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DstarToPiKPiPi0Pi0, std::vector<int>{+kKMinus, +kPiPlus, +kPiPlus, +kPi0, +kPi0}},
    {DecayChannelMain::DstarToPiKK, std::vector<int>{+kKMinus, +kKPlus, +kPiPlus}},
    {DecayChannelMain::DstarToPiKKPi0, std::vector<int>{+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DstarToPiPiPi, std::vector<int>{+kPiMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DstarToPiPiPiPi0, std::vector<int>{+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

// Lc → p K∓ π±
std::unordered_map<DecayChannelResonant, std::array<int, 2>> resoStatesLambdaC =
  {
    {DecayChannelResonant::LcToPKstar0, std::array<int, 2>{+kK0Star892, +kProton}},
    {DecayChannelResonant::LcToDeltaplusplusK, std::array<int, 2>{+2224, +kKMinus}},
    {DecayChannelResonant::LcToL1520Pi, std::array<int, 2>{+102134, +kPiPlus}},
};

std::unordered_map<DecayChannelMain, std::vector<int>> finalStatesLc =
  {
    {DecayChannelMain::LcToPKPi, std::vector<int>{+kProton, +kKMinus, +kPiPlus}},
    {DecayChannelMain::LcToPKPiPi0, std::vector<int>{+kProton, +kKMinus, +kPiPlus, +kPi0}},
    {DecayChannelMain::LcToPPiPi, std::vector<int>{+kProton, +kPiMinus, +kPiPlus}},
    {DecayChannelMain::LcToPKK, std::vector<int>{+kProton, +kKMinus, +kKPlus}}};

// Xic → p K∓ π±
std::unordered_map<DecayChannelResonant, std::array<int, 2>> resoStatesXiC =
  {
    {DecayChannelResonant::XicToPKstar0, std::array<int, 2>{-kK0Star892, +kProton}},
    {DecayChannelResonant::XicToPPhi, std::array<int, 2>{+kProton, +kPhi}},
};

std::unordered_map<DecayChannelMain, std::vector<int>> finalStatesXic =
  {
    {DecayChannelMain::XicToPKPi, std::vector<int>{+kProton, +kKMinus, +kPiPlus}},
    {DecayChannelMain::XicToPKK, std::vector<int>{+kProton, +kKMinus, +kKPlus}},
    {DecayChannelMain::XicToSPiPi, std::vector<int>{+kSigmaPlus, +kPiMinus, +kPiPlus}},
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

namespace o2::hf_corrbkg
{
using namespace o2::hf_decay;

/// Returns a map of the possible final states for a specific 3-prong particle specie
/// \param pdgMother PDG code of the mother particle
/// \return a map of final states with their corresponding PDG codes
std::unordered_map<hf_cand_3prong::DecayChannelMain, std::vector<int>> getDecayChannel3Prong(int pdgMother)
{
  switch (pdgMother) {
    case Pdg::kDPlus:
      return hf_cand_3prong::finalStatesDPlus;
    case Pdg::kDS:
      return hf_cand_3prong::finalStatesDs;
    case Pdg::kDStar:
      return hf_cand_3prong::finalStatesDstar;
    case Pdg::kLambdaCPlus:
      return hf_cand_3prong::finalStatesLc;
    case Pdg::kXiCPlus:
      return hf_cand_3prong::finalStatesXic;
    default:
      LOG(error) << "Unknown PDG code for 3-prong final states: " << pdgMother;
      return {};
  }
}

/// Returns a map of the resonant decay channels for a specific 3-prong particle specie
/// \param pdgMother PDG code of the mother particle
/// \return a map of resonant decay channels with their corresponding PDG codes
std::unordered_map<hf_cand_3prong::DecayChannelResonant, std::array<int, 2>> getResoChannels3Prong(int pdgMother)
{
  switch (pdgMother) {
    case Pdg::kDPlus:
      return hf_cand_3prong::resoStatesDPlus;
    case Pdg::kDS:
      return hf_cand_3prong::resoStatesDs;
    case Pdg::kDStar:
      return hf_cand_3prong::resoStatesDstar;
    case Pdg::kLambdaCPlus:
      return hf_cand_3prong::resoStatesLambdaC;
    case Pdg::kXiCPlus:
      return hf_cand_3prong::resoStatesXiC;
    default:
      LOG(error) << "Unknown PDG code for 3-prong final states: " << pdgMother;
      return {};
  }
}

/// Perform the matching for a single resonant channel
/// \tparam N size of the array of daughter PDG codes
/// \param arrPdgResoChn array of daughter indices
/// \param arrPdgDaugs array of PDG codes for the resonant decay
/// \return true if the resonant channel is matched, false otherwise
template <std::size_t N>
bool checkResonantDecay(std::array<int, N> arrPdgResoChn, std::array<int, N> arrPdgDaugs)
{
  for (size_t i = 0; i < N; i++) {
    bool findDaug = false;
    for (size_t j = 0; j < N; j++) {
      if (std::abs(arrPdgResoChn[i]) == std::abs(arrPdgDaugs[j])) {
        arrPdgDaugs[j] = -1; // Mark as found
        findDaug = true;
        break;
      }
    }
    if (!findDaug) {
      return false;
    }
  }
  return true;
}

/// Flag the resonant decays
/// Function documentation:
/// \tparam is3Prong bool to specify if the mother decays with a 3-prong decay
/// \tparam N size of the array of daughter PDG codes
/// \param motherPdg PDG code of the mother particle
/// \param channel decay channel flag to be set
/// \param arrDaughPdgs array of daughter PDG codes
template <bool is3Prong = false, std::size_t N>
void flagResonantDecay(int motherPdg, int8_t* channel, std::array<int, N> arrDaughPdgs)
{
  if constexpr (is3Prong) {
    std::unordered_map<hf_cand_3prong::DecayChannelResonant, std::array<int, 2>> resoStates = getResoChannels3Prong(motherPdg);
    for (const auto& [flag, pdgCodes] : resoStates) {
      if (checkResonantDecay(arrDaughPdgs, pdgCodes)) {
        *channel = flag;
        break;
      }
    }
  } else {
    if (motherPdg != Pdg::kD0) {
      return;
    }
    for (const auto& [flag, pdgCodes] : hf_cand_2prong::resoStatesD0) {
      if (checkResonantDecay(arrDaughPdgs, pdgCodes)) {
        *channel = flag;
        break;
      }
    }
  }
}
} // namespace o2::hf_corrbkg

#endif // PWGHF_CORE_DECAYCHANNELS_H_
