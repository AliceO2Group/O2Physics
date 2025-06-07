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

#include <cstdint>
#include <unordered_map>
#include <array>
#include <variant>

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

std::unordered_map<DecayChannelMain, std::vector<int> > finalStates2Prongs = 
{
    {DecayChannelMain::D0ToPiK,          std::vector<int>{+kKMinus,   +kPiPlus}},
    {DecayChannelMain::D0ToKK,           std::vector<int>{+kKMinus,   +kKPlus}},
    {DecayChannelMain::D0ToPiKPi0,       std::vector<int>{+kKMinus,   +kPiPlus, +kPi0}},
    {DecayChannelMain::D0ToPiPi,         std::vector<int>{+kPiMinus,  +kPiPlus}},
    {DecayChannelMain::D0ToPiPiPi0,      std::vector<int>{+kPiMinus,  +kPiPlus, +kPi0}}
};

std::unordered_map<DecayChannelResonant, std::array<int, 2> > resoStatesD0 = 
{
  {DecayChannelResonant::D0ToRhoplusPi,  std::array<int, 2>{213,            +kPiMinus}},
  {DecayChannelResonant::D0ToRhoplusK,   std::array<int, 2>{213,            +kKMinus}},
  {DecayChannelResonant::D0ToKstar0Pi0,  std::array<int, 2>{-kK0Star892,    +kPi0}},
  {DecayChannelResonant::D0ToKstarPi,    std::array<int, 2>{-kKPlusStar892, +kPiPlus}},
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


std::unordered_map<DecayChannelResonant, std::array<int, 2> > resoStatesDPlus = 
{
   {DecayChannelResonant::DplusToPhiPi,         std::array<int, 2>{+kPhi,       +kPiPlus}},
   {DecayChannelResonant::DplusToKstar0K,       std::array<int, 2>{-kK0Star892, +kKPlus}},
   {DecayChannelResonant::DplusToKstar1430_0K,  std::array<int, 2>{+10311,      +kKPlus}},
   {DecayChannelResonant::DplusToRho0Pi,        std::array<int, 2>{+113,        +kPiPlus}},
   {DecayChannelResonant::DplusToF2_1270Pi,     std::array<int, 2>{+225,        +kPiPlus}},
};

std::unordered_map<DecayChannelMain, std::vector<int> > finalStatesDPlus = 
{
  {DecayChannelMain::DplusToPiKPi,      std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus}},
  {DecayChannelMain::DplusToPiKK,       std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus}},
  {DecayChannelMain::DplusToPiKPiPi0,   std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus, +kPi0}},
  {DecayChannelMain::DplusToPiPiPi,     std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus}},
  // {DecayChannelMain::PiPiPiPi0,     std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus, +kPi0}},
};
    
// Ds± → K± K∓ π±
std::unordered_map<DecayChannelResonant, std::array<int, 2> > resoStatesDs = 
{
  {DecayChannelResonant::DsToPhiPi,       std::array<int, 2>{+kPhi,        +kPiPlus}},
  {DecayChannelResonant::DsToPhiRhoplus,  std::array<int, 2>{+kPhi,        +213}},
  {DecayChannelResonant::DsToKstar0K,     std::array<int, 2>{-kK0Star892,  +kKPlus}},
  {DecayChannelResonant::DsToKstar0Pi,    std::array<int, 2>{+kK0Star892,  +kPiPlus}},
  {DecayChannelResonant::DsToRho0Pi,      std::array<int, 2>{113,          +kPiPlus}},
  {DecayChannelResonant::DsToRho0K,       std::array<int, 2>{113,          +kKPlus}},
  {DecayChannelResonant::DsToF2_1270Pi,   std::array<int, 2>{225,          +kPiPlus}},
  {DecayChannelResonant::DsToF0_1370K,    std::array<int, 2>{10221,        +kKPlus}},
  {DecayChannelResonant::DsToEtaPi,       std::array<int, 2>{221,          +kPiPlus}},
};

std::unordered_map<DecayChannelMain, std::vector<int> > finalStatesDs = 
{
  {DecayChannelMain::DsToPiKK,       std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus}},
  {DecayChannelMain::DsToPiKKPi0,    std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus, +kPi0}},
  {DecayChannelMain::DsToPiPiK,      std::vector<int>{+kKPlus,    +kPiPlus,  +kPiMinus}},
  {DecayChannelMain::DsToPiPiPi,     std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus}},
  {DecayChannelMain::DsToPiPiPiPi0,  std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus, +kPi0}},
};

// // Dstar → K± K∓ π±
// std::unordered_map<DecayChannelResonant, std::array<int, 2> > resoStatesDStarD0 = 
// {
//   {DecayChannelResonant::RhoPi,         std::array<int, 2>{213,            +kPiMinus}},
//   {DecayChannelResonant::RhoK,          std::array<int, 2>{213,            +kKMinus}},
//   {DecayChannelResonant::K0starPi0,     std::array<int, 2>{-kK0Star892,    +kPi0}},
//   {DecayChannelResonant::K0starPiPlus,  std::array<int, 2>{-kKPlusStar892, +kPiPlus}},
// };

// std::unordered_map<DecayChannelMain, std::vector<int> > finalStatesDStar = 
// {
//   {DecayChannelMain::KKPi,           std::vector<int>{+kKMinus,   +kKPlus,   +kPiPlus}},
//   {DecayChannelMain::KMinusPiPi,     std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus}},
//   {DecayChannelMain::KMinusPiPiPi0,  std::vector<int>{+kKMinus,   +kPiPlus,  +kPiPlus, +kPi0}},
//   {DecayChannelMain::PiPiPi,         std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus}},
//   {DecayChannelMain::PiPiPiPi0,      std::vector<int>{+kPiMinus,  +kPiPlus,  +kPiPlus, +kPi0}},
// };

    // Lc → p K∓ π±
std::unordered_map<DecayChannelResonant, std::array<int, 2> > resoStatesLambdaC = 
{
  {DecayChannelResonant::LcToPKstar0,        std::array<int, 2>{+kK0Star892, +kProton}},
  {DecayChannelResonant::LcToDeltaplusplusK, std::array<int, 2>{+2224,       +kKMinus}},
  {DecayChannelResonant::LcToL1520Pi,        std::array<int, 2>{+102134,     +kPiPlus}},
};

std::unordered_map<DecayChannelMain, std::vector<int> > finalStatesLc = 
{
  {DecayChannelMain::LcToPKPi,     std::vector<int>{+kProton,   +kKMinus,  +kPiPlus}},
  {DecayChannelMain::LcToPKPiPi0,  std::vector<int>{+kProton,   +kKMinus,  +kPiPlus, +kPi0}},
  {DecayChannelMain::LcToPPiPi,    std::vector<int>{+kProton,   +kPiMinus, +kPiPlus}},
  {DecayChannelMain::LcToPKK,      std::vector<int>{+kProton,   +kKMinus,  +kKPlus}}
};
    
// Xic → p K∓ π±
std::unordered_map<DecayChannelResonant, std::array<int, 2> > resoStatesXiC = 
{
  {DecayChannelResonant::XicToPKstar0,  std::array<int, 2>{-kK0Star892, +kProton}},
  {DecayChannelResonant::XicToPPhi,     std::array<int, 2>{+kProton,    +kPhi}},
};

std::unordered_map<DecayChannelMain, std::vector<int> > finalStatesXic = 
{
  {DecayChannelMain::XicToPKPi,   std::vector<int>{+kProton,   +kKMinus,  +kPiPlus}},
  {DecayChannelMain::XicToPKK,    std::vector<int>{+kProton,   +kKMinus,  +kKPlus}},
  {DecayChannelMain::XicToSPiPi,  std::vector<int>{+kSigmaPlus,   +kPiMinus,  +kPiPlus}},
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
} // namespace hf_cand_beauty
} // namespace o2::hf_decay

// using namespace o2::hf_decay;

namespace o2::hf_corrbkg
{
  std::unordered_map<o2::hf_decay::hf_cand_3prong::DecayChannelMain, std::vector<int> > getDecayChannel3Prong(int pdgMother) {
    switch (pdgMother) {
      case Pdg::kDPlus:
        return o2::hf_decay::hf_cand_3prong::finalStatesDPlus;
      case Pdg::kDS:
        return o2::hf_decay::hf_cand_3prong::finalStatesDs;
      case Pdg::kDStar:
        return o2::hf_decay::hf_cand_3prong::finalStatesDPlus;
      case Pdg::kLambdaCPlus:
        return o2::hf_decay::hf_cand_3prong::finalStatesLc;
      case Pdg::kXiCPlus:
        return o2::hf_decay::hf_cand_3prong::finalStatesXic;
      default:
        LOG(error) << "Unknown PDG code for 3-prong final states: " << pdgMother;
        return {};
    }
  }

  bool checkResonantDecay(std::vector<int> arrDaughIndex, std::array<int, 2> arrPDGResonant) {
      // LOG(info) << "Entered checkResonantDecay with daughters: " << arrDaughIndex[0] << ", " << arrDaughIndex[1] << " and resonant PDG codes: " << arrPDGResonant[0] << ", " << arrPDGResonant[1];
      // LOG(info) << "arrDaughIndex.size(): " << arrDaughIndex.size() << ", arrPDGResonant.size(): " << arrPDGResonant.size();
      std::array<int, 2> arrPDGResonantAbs = {std::abs(arrPDGResonant[0]), std::abs(arrPDGResonant[1])};
      LOG(info) << "Testing: " << arrDaughIndex[0] << ", " << arrDaughIndex[1] << " matching PDG codes: " << arrPDGResonant[0] << ", " << arrPDGResonant[1];
      for (int i = 0; i < 2; i++) {
        LOG(info) << "Checking daughter index: " << arrDaughIndex[i];
        bool findDaug = false;
        for (int j = 0; j < 2; j++) {
          LOG(info) << "Checking daughter PDG: " << arrDaughIndex[i] << " against resonant PDG: " << arrPDGResonantAbs[j];
          if (std::abs(arrDaughIndex[i]) == arrPDGResonantAbs[j]) {
            arrPDGResonantAbs[j] = -1; // Mark as found
            LOG(info) << "Matched!";
            findDaug = true;
            break; // If a daughter matches, break the inner loop
          }
        }
        if (!findDaug) {
          LOG(info) << "Returning false";
          return false; // If any daughter does not match, return false
        }
      }
      LOG(info) << "Resonant decay found with daughters: " << arrDaughIndex[0] << ", " << arrDaughIndex[1] << " matching PDG codes: " << arrPDGResonant[0] << ", " << arrPDGResonant[1];
      return true;
    }

    /// Check if the decay is resonant
    /// \tparam arrDaughIndex index of the particle daughters at resonance level
    /// \tparam arrPDGResonant PDG code of the resonant decay
    /// \return true if the decay is resonant
    void flagResonantDecay(int motherPdg, int8_t* channel, std::vector<int> arrDaughIndex) {
      switch (motherPdg) {
        case Pdg::kD0:
          for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_2prong::resoStatesD0) {
            std::cout << "Checking D0 resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << " vs " << arrDaughIndex[0] << " " << arrDaughIndex[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "D0 resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kDPlus:
          for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::resoStatesDPlus) {
            // std::cout << "Checking DPlus resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "D+ resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kDS:
          for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::resoStatesDs) {
            // std::cout << "Checking DS resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "Ds resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kDStar:
          // for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::resoStatesDStarD0) {
          //   std::cout << "Checking DStar resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
          //   if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
          //     *channel = flag;
          //     LOG(info) << "Dstar resonant decay found with channel: " << static_cast<int>(*channel);
          //     break;
          //   }
          // }
          LOG(info) << "Dstar resonant decay not found, checking D0 resonances";
          break;
        case Pdg::kLambdaCPlus:
          for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::resoStatesLambdaC) {
            // std::cout << "Checking LambdaC resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "Lc resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
        case Pdg::kXiCPlus:
          for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::resoStatesXiC) {
            // std::cout << "Checking XiC resonant decay with flag: " << flag << ", pdgCodes: " << pdgCodes[0] << ", " << pdgCodes[1] << std::endl;
            if (checkResonantDecay(arrDaughIndex, pdgCodes)) {
              *channel = flag;
              LOG(info) << "Xic resonant decay found with channel: " << static_cast<int>(*channel);
              break;
            }
          }
          break;
      }
      LOG(info) << "Leaving function with channel: " << static_cast<int>(*channel);
    }

} // namespace o2::hf_corrbkg

#endif // PWGHF_CORE_DECAYCHANNELS_H_
