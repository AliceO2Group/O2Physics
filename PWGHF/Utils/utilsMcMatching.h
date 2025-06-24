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

/// \file utilsMcMatching.h
/// \brief Mapping of MC flags contained in DecayChannels.h
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN

#ifndef PWGHF_UTILS_UTILSMCMATCHING_H_
#define PWGHF_UTILS_UTILSMCMATCHING_H_

#include "PWGHF/Core/DecayChannels.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>

#include <TPDGCode.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace o2::hf_decay
{

namespace hf_cand_2prong
{

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersD0Main{
  {DecayChannelMain::D0ToPiK, {+kKMinus, +kPiPlus}},
  {DecayChannelMain::D0ToPiKPi0, {+kKMinus, +kPiPlus, +kPi0}},
  {DecayChannelMain::D0ToPiPi, {+kPiMinus, +kPiPlus}},
  {DecayChannelMain::D0ToPiPiPi0, {+kPiMinus, +kPiPlus, +kPi0}},
  {DecayChannelMain::D0ToKK, {+kKMinus, +kKPlus}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersD0Resonant{
  {DecayChannelResonant::D0ToRhoplusPi, {+kRho770Plus, +kPiMinus}},
  {DecayChannelResonant::D0ToRhoplusK, {+kRho770Plus, +kKMinus}},
  {DecayChannelResonant::D0ToKstar0Pi0, {-o2::constants::physics::kK0Star892, +kPi0}},
  {DecayChannelResonant::D0ToKstarPi, {-o2::constants::physics::kKPlusStar892, +kPiPlus}},
};

} // namespace hf_cand_2prong

namespace hf_cand_3prong
{

// D±
static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersDplusMain{
  {DecayChannelMain::DplusToPiKPi, {+kKMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DplusToPiKK, {+kKMinus, +kKPlus, +kPiPlus}},
  {DecayChannelMain::DplusToPiKPiPi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DplusToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersDplusResonant{
  {DecayChannelResonant::DplusToPhiPi, {+o2::constants::physics::kPhi, +kPiPlus}},
  {DecayChannelResonant::DplusToKstar0K, {-o2::constants::physics::kK0Star892, +kKPlus}},
  {DecayChannelResonant::DplusToKstar1430_0K, {-10311, +kKPlus}},
  {DecayChannelResonant::DplusToRho0Pi, {+kRho770_0, +kPiPlus}},
  {DecayChannelResonant::DplusToF2_1270Pi, {+225, +kPiPlus}},
};

// Ds±
static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersDsMain{
  {DecayChannelMain::DsToPiKK, {+kKMinus, +kKPlus, +kPiPlus}},
  {DecayChannelMain::DsToPiKKPi0, {+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DsToPiPiK, {+kKPlus, +kPiPlus, +kPiMinus}},
  {DecayChannelMain::DsToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DsToPiPiPiPi0, {+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersDsResonant{
  {DecayChannelResonant::DsToPhiPi, {+o2::constants::physics::kPhi, +kPiPlus}},
  {DecayChannelResonant::DsToPhiRhoplus, {+o2::constants::physics::kPhi, +kRho770Plus}},
  {DecayChannelResonant::DsToKstar0K, {-o2::constants::physics::kK0Star892, +kKPlus}},
  {DecayChannelResonant::DsToKstar0Pi, {+o2::constants::physics::kK0Star892, +kPiPlus}},
  {DecayChannelResonant::DsToRho0Pi, {+kRho770_0, +kPiPlus}},
  {DecayChannelResonant::DsToRho0K, {+kRho770_0, +kKPlus}},
  {DecayChannelResonant::DsToF2_1270Pi, {225, +kPiPlus}},
  {DecayChannelResonant::DsToF0_1370K, {10221, +kKPlus}},
  {DecayChannelResonant::DsToEtaPi, {221, +kPiPlus}},
};

// Dstar
static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersDstarMain{
  {DecayChannelMain::DstarToPiKPi, {+kKMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DstarToPiKPiPi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DstarToPiKPiPi0Pi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0, +kPi0}},
  {DecayChannelMain::DstarToPiKK, {+kKMinus, +kKPlus, +kPiPlus}},
  {DecayChannelMain::DstarToPiKKPi0, {+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DstarToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DstarToPiPiPiPi0, {+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersDstarResonant{
  {DecayChannelResonant::DstarToD0ToRhoplusPi, {+kRho770Plus, +kPiMinus}},
  {DecayChannelResonant::DstarToD0ToRhoplusK, {+kRho770Plus, +kKMinus}},
  {DecayChannelResonant::DstarToD0ToKstar0Pi0, {-o2::constants::physics::kK0Star892, +kPi0}},
  {DecayChannelResonant::DstarToD0ToKstarPi, {-o2::constants::physics::kKPlusStar892, +kPiPlus}},
  {DecayChannelResonant::DstarToDplusToPhiPi, {+o2::constants::physics::kPhi, +kPiPlus}},
  {DecayChannelResonant::DstarToDplusToKstar0K, {-o2::constants::physics::kK0Star892, +kKPlus}},
  {DecayChannelResonant::DstarToDplusToKstar1430_0K, {-10311, +kKPlus}},
  {DecayChannelResonant::DstarToDplusToRho0Pi, {+kRho770_0, +kPiPlus}},
  {DecayChannelResonant::DstarToDplusToF2_1270Pi, {+225, +kPiPlus}},
};

// Lc
static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersLcMain{
  {DecayChannelMain::LcToPKPi, {+kProton, +kKMinus, +kPiPlus}},
  {DecayChannelMain::LcToPKPiPi0, {+kProton, +kKMinus, +kPiPlus, +kPi0}},
  {DecayChannelMain::LcToPPiPi, {+kProton, +kPiMinus, +kPiPlus}},
  {DecayChannelMain::LcToPKK, {+kProton, +kKMinus, +kKPlus}}};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersLcResonant{
  {DecayChannelResonant::LcToPKstar0, {-o2::constants::physics::kK0Star892, +kProton}},
  {DecayChannelResonant::LcToDeltaplusplusK, {+2224, +kKMinus}},
  {DecayChannelResonant::LcToL1520Pi, {+102134, +kPiPlus}},
};

// Xic
static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersXicMain{
  {DecayChannelMain::XicToPKPi, {+kProton, +kKMinus, +kPiPlus}},
  {DecayChannelMain::XicToPKK, {+kProton, +kKMinus, +kKPlus}},
  {DecayChannelMain::XicToSPiPi, {+kSigmaPlus, +kPiMinus, +kPiPlus}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersXicResonant{
  {DecayChannelResonant::XicToPKstar0, {-o2::constants::physics::kK0Star892, +kProton}},
  {DecayChannelResonant::XicToPPhi, {+kProton, +o2::constants::physics::kPhi}},
};

/// Returns a map of the possible final states for a specific 3-prong particle specie
/// \param pdgMother PDG code of the mother particle
/// \return a map of final states with their corresponding PDG codes
inline std::unordered_map<DecayChannelMain, const std::vector<int>> getDecayChannelMain(int pdgMother)
{
  switch (pdgMother) {
    case o2::constants::physics::Pdg::kDPlus:
      return daughtersDplusMain;
    case o2::constants::physics::Pdg::kDS:
      return daughtersDsMain;
    case o2::constants::physics::Pdg::kDStar:
      return daughtersDstarMain;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      return daughtersLcMain;
    case o2::constants::physics::Pdg::kXiCPlus:
      return daughtersXicMain;
    default:
      LOG(fatal) << "Unknown PDG code for 3-prong final states: " << pdgMother;
      return {};
  }
}
} // namespace hf_cand_3prong

/// Compare an array of PDG codes with an expected array
/// \tparam N size of the arrays to be compared
/// \param arrPdgTested array of PDG codes to be tested
/// \param arrPdgExpected array of the expected PDG codes
/// \return true if the arrays are equal, false otherwise
template <std::size_t N>
inline bool checkDecayChannel(std::array<int, N> const& arrPdgTested, std::array<int, N> arrPdgExpected)
{
  for (std::size_t i = 0; i < N; i++) {
    bool foundPdg = false;
    for (std::size_t j = 0; j < N; j++) {
      if (std::abs(arrPdgTested[i]) == std::abs(arrPdgExpected[j])) {
        arrPdgExpected[j] = -1; // Mark as found
        foundPdg = true;
        break;
      }
    }
    if (!foundPdg) {
      return false;
    }
  }
  return true;
}

/// Flag the resonant decays
/// \tparam N size of the array of daughter PDG codes
/// \param motherPdg PDG code of the mother particle
/// \param arrDaughPdgs array of daughter PDG codes
/// \return the flag for the matched resonant decay channel
template <std::size_t N>
inline int8_t flagResonantDecay(const int motherPdg, std::array<int, N> const& arrDaughPdgs)
{
  switch (motherPdg) {
    case o2::constants::physics::Pdg::kD0:
      for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_2prong::daughtersD0Resonant) {
        if (o2::hf_decay::checkDecayChannel(arrDaughPdgs, pdgCodes)) {
          return flag;
        }
      }
      break;
    case o2::constants::physics::Pdg::kDPlus:
      for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::daughtersDplusResonant) {
        if (o2::hf_decay::checkDecayChannel(arrDaughPdgs, pdgCodes)) {
          return flag;
        }
      }
      break;
    case o2::constants::physics::Pdg::kDS:
      for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::daughtersDsResonant) {
        if (o2::hf_decay::checkDecayChannel(arrDaughPdgs, pdgCodes)) {
          return flag;
        }
      }
      break;
    case o2::constants::physics::Pdg::kDStar:
      for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::daughtersDstarResonant) {
        if (o2::hf_decay::checkDecayChannel(arrDaughPdgs, pdgCodes)) {
          return flag;
        }
      }
      break;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::daughtersLcResonant) {
        if (o2::hf_decay::checkDecayChannel(arrDaughPdgs, pdgCodes)) {
          return flag;
        }
      }
      break;
    case o2::constants::physics::Pdg::kXiCPlus:
      for (const auto& [flag, pdgCodes] : o2::hf_decay::hf_cand_3prong::daughtersXicResonant) {
        if (o2::hf_decay::checkDecayChannel(arrDaughPdgs, pdgCodes)) {
          return flag;
        }
      }
      break;
    default:
      LOG(fatal) << "Unknown PDG code for 3-prong final states: " << motherPdg;
      return -1;
  }
  return 0;
}

/// Flip the sign of a specific PDG code in an array
/// of PDG codes associated to an antiparticle.
/// \tparam N size of the array of PDG codes
/// \param motherPdgCode PDG code of the mother particle
/// \param partPdgCode PDG code to be flipped
/// \param arrFinalStatePdgs array of PDG codes to be modified
template <std::size_t N>
inline void changeFinalStatePdgSign(const int motherPdgCode, const int partPdgCode, std::array<int, N>& arrFinalStatePdgs)
{
  if (motherPdgCode >= 0) {
    return;
  }
  for (auto& part : arrFinalStatePdgs) { // o2-linter: disable=const-ref-in-for-loop (arrFinalStatePdgs entries are modified)
    if (part == partPdgCode) {
      part = -part;
    }
  }
}
} // namespace o2::hf_decay

#endif // PWGHF_UTILS_UTILSMCMATCHING_H_
