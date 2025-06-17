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

static const std::unordered_map<DecayChannelMain, const std::vector<int>> DaughtersD0Main{
  {DecayChannelMain::D0ToPiK, {+kKMinus, +kPiPlus}},
  {DecayChannelMain::D0ToPiKPi0, {+kKMinus, +kPiPlus, +kPi0}},
  {DecayChannelMain::D0ToPiPi, {+kPiMinus, +kPiPlus}},
  {DecayChannelMain::D0ToPiPiPi0, {+kPiMinus, +kPiPlus, +kPi0}},
  {DecayChannelMain::D0ToKK, {+kKMinus, +kKPlus}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> DaughtersD0Resonant{
  {DecayChannelResonant::D0ToRhoplusPi, {+kRho770Plus, +kPiMinus}},
  {DecayChannelResonant::D0ToRhoplusK, {+kRho770Plus, +kKMinus}},
  {DecayChannelResonant::D0ToKstar0Pi0, {-o2::constants::physics::kK0Star892, +kPi0}},
  {DecayChannelResonant::D0ToKstarPi, {-o2::constants::physics::kKPlusStar892, +kPiPlus}},
};

} // namespace hf_cand_2prong

namespace hf_cand_3prong
{

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> DaughtersDPlusResonant{
  {DecayChannelResonant::DplusToPhiPi, {+o2::constants::physics::kPhi, +kPiPlus}},
  {DecayChannelResonant::DplusToKstar0K, {-o2::constants::physics::kK0Star892, +kKPlus}},
  {DecayChannelResonant::DplusToKstar1430_0K, {+10311, +kKPlus}},
  {DecayChannelResonant::DplusToRho0Pi, {+kRho770_0, +kPiPlus}},
  {DecayChannelResonant::DplusToF2_1270Pi, {+225, +kPiPlus}},
};

static const std::unordered_map<DecayChannelMain, const std::vector<int>> DaughtersDPlusMain{
  {DecayChannelMain::DplusToPiKPi, {+kKMinus, +kKPlus, +kPiPlus}},
  {DecayChannelMain::DplusToPiKK, {+kKMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DplusToPiKPiPi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DplusToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
};

// Ds± → K± K∓ π±
static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> DaughtersDsResonant{
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

static const std::unordered_map<DecayChannelMain, const std::vector<int>> DaughtersDsMain{
  {DecayChannelMain::DsToPiKK, {+kKMinus, +kKPlus, +kPiPlus}},
  {DecayChannelMain::DsToPiKKPi0, {+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DsToPiPiK, {+kKPlus, +kPiPlus, +kPiMinus}},
  {DecayChannelMain::DsToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DsToPiPiPiPi0, {+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

// Dstar → K± K∓ π±
static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> DaughtersDstarResonant{
  {DecayChannelResonant::DstarToD0ToRhoplusPi, {+kRho770Plus, +kPiMinus}},
  {DecayChannelResonant::DstarToD0ToRhoplusK, {+kRho770Plus, +kKMinus}},
  {DecayChannelResonant::DstarToD0ToKstar0Pi0, {-o2::constants::physics::kK0Star892, +kPi0}},
  {DecayChannelResonant::DstarToD0ToKstarPi, {-o2::constants::physics::kKPlusStar892, +kPiPlus}},
  {DecayChannelResonant::DstarToDplusToPhiPi, {+o2::constants::physics::kPhi, +kPiPlus}},
  {DecayChannelResonant::DstarToDplusToKstar0K, {-o2::constants::physics::kK0Star892, +kKPlus}},
  {DecayChannelResonant::DstarToDplusToKstar1430_0K, {+10311, +kKPlus}},
  {DecayChannelResonant::DstarToDplusToRho0Pi, {+kRho770_0, +kPiPlus}},
  {DecayChannelResonant::DstarToDplusToF2_1270Pi, {+225, +kPiPlus}},
};

static const std::unordered_map<DecayChannelMain, const std::vector<int>> DaughtersDstarMain{
  {DecayChannelMain::DstarToPiKPi, {+kKMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DstarToPiKPiPi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DstarToPiKPiPi0Pi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0, +kPi0}},
  {DecayChannelMain::DstarToPiKK, {+kKMinus, +kKPlus, +kPiPlus}},
  {DecayChannelMain::DstarToPiKKPi0, {+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
  {DecayChannelMain::DstarToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
  {DecayChannelMain::DstarToPiPiPiPi0, {+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

// Lc → p K∓ π±
static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> DaughtersLcResonant{
  {DecayChannelResonant::LcToPKstar0, {+o2::constants::physics::kK0Star892, +kProton}},
  {DecayChannelResonant::LcToDeltaplusplusK, {+2224, +kKMinus}},
  {DecayChannelResonant::LcToL1520Pi, {+102134, +kPiPlus}},
};

static const std::unordered_map<DecayChannelMain, const std::vector<int>> DaughtersLcMain{
  {DecayChannelMain::LcToPKPi, {+kProton, +kKMinus, +kPiPlus}},
  {DecayChannelMain::LcToPKPiPi0, {+kProton, +kKMinus, +kPiPlus, +kPi0}},
  {DecayChannelMain::LcToPPiPi, {+kProton, +kPiMinus, +kPiPlus}},
  {DecayChannelMain::LcToPKK, {+kProton, +kKMinus, +kKPlus}}};

// Xic → p K∓ π±
static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> DaughtersXiCResonant{
  {DecayChannelResonant::XicToPKstar0, {-o2::constants::physics::kK0Star892, +kProton}},
  {DecayChannelResonant::XicToPPhi, {+kProton, +o2::constants::physics::kPhi}},
};

static const std::unordered_map<DecayChannelMain, const std::vector<int>> DaughtersXiCMain{
  {DecayChannelMain::XicToPKPi, {+kProton, +kKMinus, +kPiPlus}},
  {DecayChannelMain::XicToPKK, {+kProton, +kKMinus, +kKPlus}},
  {DecayChannelMain::XicToSPiPi, {+kSigmaPlus, +kPiMinus, +kPiPlus}},
};

/// Returns a map of the possible final states for a specific 3-prong particle specie
/// \param pdgMother PDG code of the mother particle
/// \return a map of final states with their corresponding PDG codes
inline std::unordered_map<DecayChannelMain, const std::vector<int>> getDecayChannel3Prong(int pdgMother)
{
  switch (pdgMother) {
    case o2::constants::physics::Pdg::kDPlus:
      return DaughtersDPlusMain;
    case o2::constants::physics::Pdg::kDS:
      return DaughtersDsMain;
    case o2::constants::physics::Pdg::kDStar:
      return DaughtersDstarMain;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      return DaughtersLcMain;
    case o2::constants::physics::Pdg::kXiCPlus:
      return DaughtersXiCMain;
    default:
      LOG(error) << "Unknown PDG code for 3-prong final states: " << pdgMother;
      return {};
  }
}

/// Returns a map of the resonant decay channels for a specific 3-prong particle specie
/// \param pdgMother PDG code of the mother particle
/// \return a map of resonant decay channels with their corresponding PDG codes
inline std::unordered_map<DecayChannelResonant, const std::array<int, 2>> getResoChannels3Prong(int pdgMother)
{
  switch (pdgMother) {
    case o2::constants::physics::Pdg::kDPlus:
      return DaughtersDPlusResonant;
    case o2::constants::physics::Pdg::kDS:
      return DaughtersDsResonant;
    case o2::constants::physics::Pdg::kDStar:
      return DaughtersDstarResonant;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      return DaughtersLcResonant;
    case o2::constants::physics::Pdg::kXiCPlus:
      return DaughtersXiCResonant;
    default:
      LOG(error) << "Unknown PDG code for 3-prong final states: " << pdgMother;
      return {};
  }
}

} // namespace hf_cand_3prong

/// Perform the matching for a single resonant channel
/// \tparam N size of the array of daughter PDG codes
/// \param arrPdgResoChn array of daughter indices
/// \param arrPdgDaugs array of PDG codes for the resonant decay
/// \return true if the resonant channel is matched, false otherwise
template <std::size_t N>
inline bool checkResonantDecay(std::array<int, N> const& arrPdgResoChn, std::array<int, N> arrPdgDaugs)
{
  for (std::size_t i = 0; i < N; i++) {
    bool findDaug = false;
    for (std::size_t j = 0; j < N; j++) {
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
inline void flagResonantDecay(int motherPdg, int8_t* channel, std::array<int, N> const& arrDaughPdgs)
{
  if constexpr (is3Prong) {
    std::unordered_map<o2::hf_decay::hf_cand_3prong::DecayChannelResonant, const std::array<int, 2>> resoStates = o2::hf_decay::hf_cand_3prong::getResoChannels3Prong(motherPdg);
    for (const auto& [flag, pdgCodes] : resoStates) {
      if (o2::hf_decay::checkResonantDecay(arrDaughPdgs, pdgCodes)) {
        *channel = flag;
        break;
      }
    }
  } else {
    if (motherPdg != o2::constants::physics::Pdg::kD0) {
      return;
    }
    for (const auto& [flag, pdgCodes] : hf_cand_2prong::DaughtersD0Resonant) {
      if (o2::hf_decay::checkResonantDecay(arrDaughPdgs, pdgCodes)) {
        *channel = flag;
        break;
      }
    }
  }
}
} // namespace o2::hf_decay

#endif // PWGHF_UTILS_UTILSMCMATCHING_H_
