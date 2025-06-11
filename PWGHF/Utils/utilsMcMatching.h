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

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "PWGHF/Core/DecayChannels.h"

namespace o2::hf_corrbkg
{
using namespace o2::hf_decay
using namespace o2::constants::physics;

namespace 2prong
{

std::unordered_map<hf_cand_2prong::DecayChannelMain, const std::vector<int>> finalStates2Prongs =
  {
    {DecayChannelMain::D0ToPiK, {+kKMinus, +kPiPlus}},
    {DecayChannelMain::D0ToPiKPi0, {+kKMinus, +kPiPlus, +kPi0}},
    {DecayChannelMain::D0ToPiPi, {+kPiMinus, +kPiPlus}},
    {DecayChannelMain::D0ToPiPiPi0, {+kPiMinus, +kPiPlus, +kPi0}},
    {DecayChannelMain::D0ToKK, {+kKMinus, +kKPlus}},
};

std::unordered_map<hf_cand_2prong::DecayChannelResonant, std::array<int, 2>> resoStatesD0 =
  {
    {DecayChannelResonant::D0ToRhoplusPi, {213, +kPiMinus}},
    {DecayChannelResonant::D0ToRhoplusK, {213, +kKMinus}},
    {DecayChannelResonant::D0ToKstar0Pi0, {-kK0Star892, +kPi0}},
    {DecayChannelResonant::D0ToKstarPi, {-kKPlusStar892, +kPiPlus}},
};

} // namespace 2prong

namespace 3prong
{

std::unordered_map<hf_cand_3prong::DecayChannelResonant, const std::array<int, 2>> resoStatesDPlus =
  {
    {DecayChannelResonant::DplusToPhiPi, {+kPhi, +kPiPlus}},
    {DecayChannelResonant::DplusToKstar0K, {-kK0Star892, +kKPlus}},
    {DecayChannelResonant::DplusToKstar1430_0K, {+10311, +kKPlus}},
    {DecayChannelResonant::DplusToRho0Pi, {+113, +kPiPlus}},
    {DecayChannelResonant::DplusToF2_1270Pi, {+225, +kPiPlus}},
};

std::unordered_map<hf_cand_3prong::DecayChannelMain, const std::vector<int>> finalStatesDPlus =
  {
    {DecayChannelMain::DplusToPiKPi, {+kKMinus, +kKPlus, +kPiPlus}},
    {DecayChannelMain::DplusToPiKK, {+kKMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DplusToPiKPiPi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DplusToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
};

// Ds± → K± K∓ π±
std::unordered_map<hf_cand_3prong::DecayChannelResonant, const std::array<int, 2>> resoStatesDs =
  {
    {DecayChannelResonant::DsToPhiPi, {+kPhi, +kPiPlus}},
    {DecayChannelResonant::DsToPhiRhoplus, {+kPhi, +213}},
    {DecayChannelResonant::DsToKstar0K, {-kK0Star892, +kKPlus}},
    {DecayChannelResonant::DsToKstar0Pi, {+kK0Star892, +kPiPlus}},
    {DecayChannelResonant::DsToRho0Pi, {113, +kPiPlus}},
    {DecayChannelResonant::DsToRho0K, {113, +kKPlus}},
    {DecayChannelResonant::DsToF2_1270Pi, {225, +kPiPlus}},
    {DecayChannelResonant::DsToF0_1370K, {10221, +kKPlus}},
    {DecayChannelResonant::DsToEtaPi, {221, +kPiPlus}},
};

std::unordered_map<hf_cand_3prong::DecayChannelMain, const std::vector<int>> finalStatesDs =
  {
    {DecayChannelMain::DsToPiKK, {+kKMinus, +kKPlus, +kPiPlus}},
    {DecayChannelMain::DsToPiKKPi0, {+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DsToPiPiK, {+kKPlus, +kPiPlus, +kPiMinus}},
    {DecayChannelMain::DsToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DsToPiPiPiPi0, {+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

// Dstar → K± K∓ π±
std::unordered_map<hf_cand_3prong::DecayChannelResonant, const std::array<int, 2>> resoStatesDstar =
  {
    {DecayChannelResonant::DstarToD0ToRhoplusPi, {213, +kPiMinus}},
    {DecayChannelResonant::DstarToD0ToRhoplusK, {213, +kKMinus}},
    {DecayChannelResonant::DstarToD0ToKstar0Pi0, {-kK0Star892, +kPi0}},
    {DecayChannelResonant::DstarToD0ToKstarPi, {-kKPlusStar892, +kPiPlus}},
    {DecayChannelResonant::DstarToDplusToPhiPi, {+kPhi, +kPiPlus}},
    {DecayChannelResonant::DstarToDplusToKstar0K, {-kK0Star892, +kKPlus}},
    {DecayChannelResonant::DstarToDplusToKstar1430_0K, {+10311, +kKPlus}},
    {DecayChannelResonant::DstarToDplusToRho0Pi, {+113, +kPiPlus}},
    {DecayChannelResonant::DstarToDplusToF2_1270Pi, {+225, +kPiPlus}},
};

std::unordered_map<hf_cand_3prong::DecayChannelMain, const std::vector<int>> finalStatesDstar =
  {
    {DecayChannelMain::DstarToPiKPi, {+kKMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DstarToPiKPiPi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DstarToPiKPiPi0Pi0, {+kKMinus, +kPiPlus, +kPiPlus, +kPi0, +kPi0}},
    {DecayChannelMain::DstarToPiKK, {+kKMinus, +kKPlus, +kPiPlus}},
    {DecayChannelMain::DstarToPiKKPi0, {+kKMinus, +kKPlus, +kPiPlus, +kPi0}},
    {DecayChannelMain::DstarToPiPiPi, {+kPiMinus, +kPiPlus, +kPiPlus}},
    {DecayChannelMain::DstarToPiPiPiPi0, {+kPiMinus, +kPiPlus, +kPiPlus, +kPi0}},
};

// Lc → p K∓ π±
std::unordered_map<hf_cand_3prong::DecayChannelResonant, const std::array<int, 2>> resoStatesLambdaC =
  {
    {DecayChannelResonant::LcToPKstar0, {+kK0Star892, +kProton}},
    {DecayChannelResonant::LcToDeltaplusplusK, {+2224, +kKMinus}},
    {DecayChannelResonant::LcToL1520Pi, {+102134, +kPiPlus}},
};

std::unordered_map<hf_cand_3prong::DecayChannelMain, const std::vector<int>> finalStatesLc =
  {
    {DecayChannelMain::LcToPKPi, {+kProton, +kKMinus, +kPiPlus}},
    {DecayChannelMain::LcToPKPiPi0, {+kProton, +kKMinus, +kPiPlus, +kPi0}},
    {DecayChannelMain::LcToPPiPi, {+kProton, +kPiMinus, +kPiPlus}},
    {DecayChannelMain::LcToPKK, {+kProton, +kKMinus, +kKPlus}}};

// Xic → p K∓ π±
std::unordered_map<hf_cand_3prong::DecayChannelResonant, const std::array<int, 2>> resoStatesXiC =
  {
    {DecayChannelResonant::XicToPKstar0, {-kK0Star892, +kProton}},
    {DecayChannelResonant::XicToPPhi, {+kProton, +kPhi}},
};

std::unordered_map<hf_cand_3prong::DecayChannelMain, const std::vector<int>> finalStatesXic =
  {
    {DecayChannelMain::XicToPKPi, {+kProton, +kKMinus, +kPiPlus}},
    {DecayChannelMain::XicToPKK, {+kProton, +kKMinus, +kKPlus}},
    {DecayChannelMain::XicToSPiPi, {+kSigmaPlus, +kPiMinus, +kPiPlus}},
};
}

/// Returns a map of the possible final states for a specific 3-prong particle specie
/// \param pdgMother PDG code of the mother particle
/// \return a map of final states with their corresponding PDG codes
std::unordered_map<hf_cand_3prong::DecayChannelMain, std::vector<int>> getDecayChannel3Prong(int pdgMother)
{
  switch (pdgMother) {
    case Pdg::kDPlus:
      return 3prong::finalStatesDPlus;
    case Pdg::kDS:
      return 3prong::finalStatesDs;
    case Pdg::kDStar:
      return 3prong::finalStatesDstar;
    case Pdg::kLambdaCPlus:
      return 3prong::finalStatesLc;
    case Pdg::kXiCPlus:
      return 3prong::finalStatesXic;
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
      return 3prong::resoStatesDPlus;
    case Pdg::kDS:
      return 3prong::resoStatesDs;
    case Pdg::kDStar:
      return 3prong::resoStatesDstar;
    case Pdg::kLambdaCPlus:
      return 3prong::resoStatesLambdaC;
    case Pdg::kXiCPlus:
      return 3prong::resoStatesXiC;
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
bool checkResonantDecay(std::array<int, N> const& arrPdgResoChn, std::array<int, N> const& arrPdgDaugs)
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
    for (const auto& [flag, pdgCodes] : 3prong::resoStatesD0) {
      if (checkResonantDecay(arrDaughPdgs, pdgCodes)) {
        *channel = flag;
        break;
      }
    }
  }
}
} // namespace o2::hf_corrbkg

#endif // PWGHF_UTILS_UTILSMCMATCHING_H_
