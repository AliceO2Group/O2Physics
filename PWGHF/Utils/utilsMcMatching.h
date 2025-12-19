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

// D0

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersD0Main{
  {DecayChannelMain::D0ToPiK, {+PDG_t::kKMinus, +PDG_t::kPiPlus}},
  {DecayChannelMain::D0ToPiKPi0, {+PDG_t::kKMinus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
  {DecayChannelMain::D0ToPiPi, {+PDG_t::kPiMinus, +PDG_t::kPiPlus}},
  {DecayChannelMain::D0ToPiPiPi0, {+PDG_t::kPiMinus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
  {DecayChannelMain::D0ToKK, {+PDG_t::kKMinus, +PDG_t::kKPlus}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersD0Resonant{
  {DecayChannelResonant::D0ToRhoplusPi, {+PDG_t::kRho770Plus, +PDG_t::kPiMinus}},
  {DecayChannelResonant::D0ToRhoplusK, {+PDG_t::kRho770Plus, +PDG_t::kKMinus}},
  {DecayChannelResonant::D0ToKstar0Pi0, {-o2::constants::physics::Pdg::kK0Star892, +PDG_t::kPi0}},
  {DecayChannelResonant::D0ToKstarPi, {-o2::constants::physics::Pdg::kKPlusStar892, +PDG_t::kPiPlus}},
};

// J/ψ

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersJpsiMain{
  {DecayChannelMain::JpsiToEE, {+PDG_t::kElectron, +PDG_t::kPositron}},
  {DecayChannelMain::JpsiToMuMu, {+PDG_t::kMuonMinus, +PDG_t::kMuonPlus}},
};
} // namespace hf_cand_2prong

namespace hf_cand_3prong
{

// D±

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersDplusMain{
  {DecayChannelMain::DplusToPiKPi, {+PDG_t::kKMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus}},
  {DecayChannelMain::DplusToPiKK, {+PDG_t::kKMinus, +PDG_t::kKPlus, +PDG_t::kPiPlus}},
  {DecayChannelMain::DplusToPiKPiPi0, {+PDG_t::kKMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
  {DecayChannelMain::DplusToPiPiPi, {+PDG_t::kPiMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersDplusResonant{
  {DecayChannelResonant::DplusToPhiPi, {+o2::constants::physics::Pdg::kPhi, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DplusToKstar0K, {-o2::constants::physics::Pdg::kK0Star892, +PDG_t::kKPlus}},
  {DecayChannelResonant::DplusToKstar1430_0K, {-10311, +PDG_t::kKPlus}},
  {DecayChannelResonant::DplusToRho0Pi, {+PDG_t::kRho770_0, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DplusToF2_1270Pi, {+225, +PDG_t::kPiPlus}},
};

// Ds±

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersDsMain{
  {DecayChannelMain::DsToPiKK, {+PDG_t::kKMinus, +PDG_t::kKPlus, +PDG_t::kPiPlus}},
  {DecayChannelMain::DsToPiKKPi0, {+PDG_t::kKMinus, +PDG_t::kKPlus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
  {DecayChannelMain::DsToPiPiK, {+PDG_t::kKPlus, +PDG_t::kPiPlus, +PDG_t::kPiMinus}},
  {DecayChannelMain::DsToPiPiPi, {+PDG_t::kPiMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus}},
  {DecayChannelMain::DsToPiPiPiPi0, {+PDG_t::kPiMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersDsResonant{
  {DecayChannelResonant::DsToPhiPi, {+o2::constants::physics::Pdg::kPhi, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DsToPhiRhoplus, {+o2::constants::physics::Pdg::kPhi, +PDG_t::kRho770Plus}},
  {DecayChannelResonant::DsToKstar0K, {-o2::constants::physics::Pdg::kK0Star892, +PDG_t::kKPlus}},
  {DecayChannelResonant::DsToKstar0Pi, {+o2::constants::physics::Pdg::kK0Star892, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DsToRho0Pi, {+PDG_t::kRho770_0, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DsToRho0K, {+PDG_t::kRho770_0, +PDG_t::kKPlus}},
  {DecayChannelResonant::DsToF2_1270Pi, {225, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DsToF0_1370K, {10221, +PDG_t::kKPlus}},
  {DecayChannelResonant::DsToEtaPi, {221, +PDG_t::kPiPlus}},
};

// D*+

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersDstarMain{
  {DecayChannelMain::DstarToPiKPi, {+PDG_t::kKMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus}},
  {DecayChannelMain::DstarToPiKPiPi0, {+PDG_t::kKMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
  {DecayChannelMain::DstarToPiKPiPi0Pi0, {+PDG_t::kKMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus, +PDG_t::kPi0, +PDG_t::kPi0}},
  {DecayChannelMain::DstarToPiKK, {+PDG_t::kKMinus, +PDG_t::kKPlus, +PDG_t::kPiPlus}},
  {DecayChannelMain::DstarToPiKKPi0, {+PDG_t::kKMinus, +PDG_t::kKPlus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
  {DecayChannelMain::DstarToPiPiPi, {+PDG_t::kPiMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus}},
  {DecayChannelMain::DstarToPiPiPiPi0, {+PDG_t::kPiMinus, +PDG_t::kPiPlus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersDstarResonant{
  {DecayChannelResonant::DstarToD0ToRhoplusPi, {+PDG_t::kRho770Plus, +PDG_t::kPiMinus}},
  {DecayChannelResonant::DstarToD0ToRhoplusK, {+PDG_t::kRho770Plus, +PDG_t::kKMinus}},
  {DecayChannelResonant::DstarToD0ToKstar0Pi0, {-o2::constants::physics::Pdg::kK0Star892, +PDG_t::kPi0}},
  {DecayChannelResonant::DstarToD0ToKstarPi, {-o2::constants::physics::Pdg::kKPlusStar892, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DstarToDplusToPhiPi, {+o2::constants::physics::Pdg::kPhi, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DstarToDplusToKstar0K, {-o2::constants::physics::Pdg::kK0Star892, +PDG_t::kKPlus}},
  {DecayChannelResonant::DstarToDplusToKstar1430_0K, {-10311, +PDG_t::kKPlus}},
  {DecayChannelResonant::DstarToDplusToRho0Pi, {+PDG_t::kRho770_0, +PDG_t::kPiPlus}},
  {DecayChannelResonant::DstarToDplusToF2_1270Pi, {+225, +PDG_t::kPiPlus}},
};

// Λc+

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersLcMain{
  {DecayChannelMain::LcToPKPi, {+PDG_t::kProton, +PDG_t::kKMinus, +PDG_t::kPiPlus}},
  {DecayChannelMain::LcToPKPiPi0, {+PDG_t::kProton, +PDG_t::kKMinus, +PDG_t::kPiPlus, +PDG_t::kPi0}},
  {DecayChannelMain::LcToPPiPi, {+PDG_t::kProton, +PDG_t::kPiMinus, +PDG_t::kPiPlus}},
  {DecayChannelMain::LcToPKK, {+PDG_t::kProton, +PDG_t::kKMinus, +PDG_t::kKPlus}}};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersLcResonant{
  {DecayChannelResonant::LcToPKstar0, {-o2::constants::physics::Pdg::kK0Star892, +PDG_t::kProton}},
  {DecayChannelResonant::LcToDeltaplusplusK, {+2224, +PDG_t::kKMinus}},
  {DecayChannelResonant::LcToL1520Pi, {+102134, +PDG_t::kPiPlus}},
  {DecayChannelResonant::LcToPPhi, {+PDG_t::kProton, +o2::constants::physics::Pdg::kPhi}},
};

// Ξc+

static const std::unordered_map<DecayChannelMain, const std::vector<int>> daughtersXicMain{
  {DecayChannelMain::XicToPKPi, {+PDG_t::kProton, +PDG_t::kKMinus, +PDG_t::kPiPlus}},
  {DecayChannelMain::XicToPKK, {+PDG_t::kProton, +PDG_t::kKMinus, +PDG_t::kKPlus}},
  {DecayChannelMain::XicToSPiPi, {+PDG_t::kSigmaPlus, +PDG_t::kPiMinus, +PDG_t::kPiPlus}},
};

static const std::unordered_map<DecayChannelResonant, const std::array<int, 2>> daughtersXicResonant{
  {DecayChannelResonant::XicToPKstar0, {-o2::constants::physics::Pdg::kK0Star892, +PDG_t::kProton}},
  {DecayChannelResonant::XicToPPhi, {+PDG_t::kProton, +o2::constants::physics::Pdg::kPhi}},
};

/// Returns a map of the possible final states for a specific 3-prong particle specie
/// \param pdgMother PDG code of the mother particle
/// \return a map of final states with their corresponding PDG codes
inline std::unordered_map<DecayChannelMain, const std::vector<int>> getDecayChannelsMain(int pdgMother)
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

namespace hf_cand_reso
{
const std::unordered_map<int, int> particlesToDstarK0s = {
  {DecayChannelMain::Ds1ToDstarK0s, constants::physics::Pdg::kDS1},
  {DecayChannelMain::Ds2starToDstarK0s, constants::physics::Pdg::kDS2Star},
  {DecayChannelMain::Ds1star2700ToDstarK0s, constants::physics::Pdg::kDS1Star2700},
  {DecayChannelMain::Ds1star2860ToDstarK0s, constants::physics::Pdg::kDS1Star2860},
  {DecayChannelMain::Ds3star2860ToDstarK0s, constants::physics::Pdg::kDS3Star2860}};
const std::unordered_map<int, int> particlesToDplusK0s = {
  {DecayChannelMain::Ds2starToDplusK0s, constants::physics::Pdg::kDS2Star}};
const std::unordered_map<int, int> particlesToDplusLambda = {
  {DecayChannelMain::Xic3055plusToDplusLambda, constants::physics::Pdg::kXiC3055Plus},
  {DecayChannelMain::Xic3080plusToDplusLambda, constants::physics::Pdg::kXiC3080Plus}};
const std::unordered_map<int, int> particlesToD0Lambda = {
  {DecayChannelMain::Xic3055zeroToD0Lambda, constants::physics::Pdg::kXiC3055_0},
  {DecayChannelMain::Xic3080zeroToD0Lambda, constants::physics::Pdg::kXiC3080_0}};
const std::unordered_map<int, int> particlesToDstarPi = {
  {DecayChannelMain::D1zeroToDstarPi, constants::physics::Pdg::kD10},
  {DecayChannelMain::D2starzeroToDstarPi, constants::physics::Pdg::kD2Star0}};
const std::unordered_map<int, int> particlesToDplusPi = {
  {DecayChannelMain::D2starzeroToDplusPi, constants::physics::Pdg::kD2Star0}};
const std::unordered_map<int, int> particlesToD0Pi = {
  {DecayChannelMain::D2starplusToD0Pi, constants::physics::Pdg::kD2StarPlus},
  {DecayChannelMain::DstarToD0Pi, constants::physics::Pdg::kDStar}};
const std::unordered_map<int, int> particlesToD0Kplus = {
  {DecayChannelMain::Ds2starToD0Kplus, constants::physics::Pdg::kDS2Star}};

enum PartialMatchMc : uint8_t {
  D0Matched = 0,
  DstarMatched,
  DplusMatched,
  K0Matched,
  LambdaMatched,
  PionMatched,
  KaonMatched,
  ProtonMatched,
  ResoPartlyMatched
};
} // namespace hf_cand_reso

/// Compare an array of PDG codes with an expected array
/// \tparam N size of the arrays to be compared
/// \param arrPdgTested array of PDG codes to be tested
/// \param arrPdgExpected array of the expected PDG codes
/// \return true if the arrays are equal, false otherwise
template <std::size_t N>
inline bool areSamePdgArrays(std::array<int, N> const& arrPdgTested, std::array<int, N> arrPdgExpected)
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

/// Flag the resonant decay channel
/// \tparam N size of the array of daughter PDG codes
/// \param pdgMother PDG code of the mother particle
/// \param arrPdgDaughters array of daughter PDG codes
/// \return the channel for the matched resonant decay channel
template <std::size_t N>
inline int8_t getDecayChannelResonant(const int pdgMother, std::array<int, N> const& arrPdgDaughters)
{
  switch (pdgMother) {
    case o2::constants::physics::Pdg::kD0:
      for (const auto& [channelResonant, arrPdgDaughtersResonant] : hf_cand_2prong::daughtersD0Resonant) {
        if (areSamePdgArrays(arrPdgDaughters, arrPdgDaughtersResonant)) {
          return channelResonant;
        }
      }
      break;
    case o2::constants::physics::Pdg::kDPlus:
      for (const auto& [channelResonant, arrPdgDaughtersResonant] : hf_cand_3prong::daughtersDplusResonant) {
        if (areSamePdgArrays(arrPdgDaughters, arrPdgDaughtersResonant)) {
          return channelResonant;
        }
      }
      break;
    case o2::constants::physics::Pdg::kDS:
      for (const auto& [channelResonant, arrPdgDaughtersResonant] : hf_cand_3prong::daughtersDsResonant) {
        if (areSamePdgArrays(arrPdgDaughters, arrPdgDaughtersResonant)) {
          return channelResonant;
        }
      }
      break;
    case o2::constants::physics::Pdg::kDStar:
      for (const auto& [channelResonant, arrPdgDaughtersResonant] : hf_cand_3prong::daughtersDstarResonant) {
        if (areSamePdgArrays(arrPdgDaughters, arrPdgDaughtersResonant)) {
          return channelResonant;
        }
      }
      break;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      for (const auto& [channelResonant, arrPdgDaughtersResonant] : hf_cand_3prong::daughtersLcResonant) {
        if (areSamePdgArrays(arrPdgDaughters, arrPdgDaughtersResonant)) {
          return channelResonant;
        }
      }
      break;
    case o2::constants::physics::Pdg::kXiCPlus:
      for (const auto& [channelResonant, arrPdgDaughtersResonant] : hf_cand_3prong::daughtersXicResonant) {
        if (areSamePdgArrays(arrPdgDaughters, arrPdgDaughtersResonant)) {
          return channelResonant;
        }
      }
      break;
    default:
      LOG(fatal) << "Unknown PDG code for 3-prong final states: " << pdgMother;
      return -1;
  }
  return 0;
}

/// Flip the sign of a specific PDG code in an array
/// of PDG codes associated to an antiparticle.
/// \tparam N size of the array of PDG codes
/// \param pdgMother PDG code of the mother particle
/// \param pdgToFlip PDG code to be flipped
/// \param arrPdg array of PDG codes to be modified
template <std::size_t N>
inline void flipPdgSign(const int pdgMother, const int pdgToFlip, std::array<int, N>& arrPdg)
{
  if (pdgMother >= 0) {
    return;
  }
  for (auto& pdg : arrPdg) { // o2-linter: disable=const-ref-in-for-loop (arrPdg entries are modified)
    if (pdg == pdgToFlip) {
      pdg = -pdg;
    }
  }
}
} // namespace o2::hf_decay

#endif // PWGHF_UTILS_UTILSMCMATCHING_H_
