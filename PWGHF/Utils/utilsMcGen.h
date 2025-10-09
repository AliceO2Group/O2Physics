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

/// \file utilsMcGen.h
/// \brief utility functions for HF MC gen. workflows
///
/// \author Nima Zardoshti, nima.zardoshti@cern.ch, CERN

#ifndef PWGHF_UTILS_UTILSMCGEN_H_
#define PWGHF_UTILS_UTILSMCGEN_H_

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Utils/utilsMcMatching.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>

#include <TPDGCode.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hf_mc_gen
{

template <typename TMcParticles, typename TMcParticlesPerColl, typename TCursor>
void fillMcMatchGen2Prong(TMcParticles const& mcParticles,
                          TMcParticlesPerColl const& mcParticlesPerMcColl,
                          TCursor& rowMcMatchGen,
                          const bool rejectBackground,
                          const bool matchCorrelatedBackground)
{
  using namespace o2::constants::physics;
  using namespace o2::hf_decay::hf_cand_2prong;

  constexpr std::size_t NDaughtersResonant{2u};

  // Match generated particles.
  for (const auto& particle : mcParticlesPerMcColl) {
    int8_t flagChannelMain = 0;
    int8_t flagChannelResonant = 0;
    int8_t origin = 0;
    int8_t sign = 0;
    std::vector<int> idxBhadMothers{};
    // Reject particles from background events
    if (particle.fromBackgroundEvent() && rejectBackground) {
      rowMcMatchGen(flagChannelMain, origin, flagChannelResonant, -1);
      continue;
    }
    if (matchCorrelatedBackground) {
      constexpr int DepthMainMax = 2; // Depth for final state matching
      constexpr int DepthResoMax = 1; // Depth for resonant decay matching
      bool matched = false;

      // TODO: J/ψ
      for (const auto& [channelMain, finalState] : daughtersD0Main) {
        if (finalState.size() == 3) { // o2-linter: disable=magic-number (partially reconstructed 3-prong decays)
          std::array<int, 3> arrPdgDaughtersMain3Prongs = std::array{finalState[0], finalState[1], finalState[2]};
          o2::hf_decay::flipPdgSign(particle.pdgCode(), +kPi0, arrPdgDaughtersMain3Prongs);
          matched = RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kD0, arrPdgDaughtersMain3Prongs, true, &sign, DepthMainMax);
        } else if (finalState.size() == 2) { // o2-linter: disable=magic-number (fully reconstructed 2-prong decays)
          std::array<int, 2> arrPdgDaughtersMain2Prongs = std::array{finalState[0], finalState[1]};
          matched = RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kD0, arrPdgDaughtersMain2Prongs, true, &sign, DepthMainMax);
        } else {
          LOG(fatal) << "Final state size not supported: " << finalState.size();
          return;
        }
        if (matched) {
          flagChannelMain = sign * channelMain;

          // Flag the resonant decay channel
          std::vector<int> arrResoDaughIndex = {};
          RecoDecay::getDaughters(particle, &arrResoDaughIndex, std::array{0}, DepthResoMax);
          std::array<int, NDaughtersResonant> arrPdgDaughters = {};
          if (arrResoDaughIndex.size() == NDaughtersResonant) {
            for (auto iProng = 0u; iProng < arrResoDaughIndex.size(); ++iProng) {
              auto daughI = mcParticles.rawIteratorAt(arrResoDaughIndex[iProng]);
              arrPdgDaughters[iProng] = daughI.pdgCode();
            }
            flagChannelResonant = o2::hf_decay::getDecayChannelResonant(Pdg::kD0, arrPdgDaughters);
          }
          break;
        }
      }
    } else {
      // D0(bar) → π± K∓
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign)) {
        flagChannelMain = sign * DecayChannelMain::D0ToPiK;
      }

      // J/ψ → e+ e−
      if (flagChannelMain == 0) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kJPsi, std::array{+kElectron, -kElectron}, true)) {
          flagChannelMain = DecayChannelMain::JpsiToEE;
        }
      }

      // J/ψ → μ+ μ−
      if (flagChannelMain == 0) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kJPsi, std::array{+kMuonPlus, -kMuonPlus}, true)) {
          flagChannelMain = DecayChannelMain::JpsiToMuMu;
        }
      }
    }

    // Check whether the particle is non-prompt (from a b quark).
    if (flagChannelMain != 0) {
      origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
    }
    if (origin == RecoDecay::OriginType::NonPrompt) {
      rowMcMatchGen(flagChannelMain, origin, flagChannelResonant, idxBhadMothers[0]);
    } else {
      rowMcMatchGen(flagChannelMain, origin, flagChannelResonant, -1);
    }
  }
}

template <typename TMcParticles, typename TMcParticlesPerColl, typename TCursor>
void fillMcMatchGen3Prong(TMcParticles const& mcParticles,
                          TMcParticlesPerColl const& mcParticlesPerMcColl,
                          TCursor& rowMcMatchGen,
                          const bool rejectBackground,
                          std::vector<int> const& pdgMothersCorrelBkg = {})
{
  using namespace o2::constants::physics;
  using namespace o2::hf_decay::hf_cand_3prong;

  constexpr std::size_t NDaughtersResonant{2u};

  // Match generated particles.
  for (const auto& particle : mcParticlesPerMcColl) {
    int8_t flagChannelMain = 0;
    int8_t flagChannelResonant = 0;
    int8_t origin = 0;
    int8_t sign = 0;
    std::vector<int> arrDaughIndex;
    std::vector<int> idxBhadMothers{};
    std::array<int, NDaughtersResonant> arrPdgDaugResonant{};
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantLcToPKstar0{daughtersLcResonant.at(DecayChannelResonant::LcToPKstar0)};               // Λc± → p± K*
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantLcToDeltaplusplusK{daughtersLcResonant.at(DecayChannelResonant::LcToDeltaplusplusK)}; // Λc± → Δ(1232)±± K∓
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantLcToL1520Pi{daughtersLcResonant.at(DecayChannelResonant::LcToL1520Pi)};               // Λc± → Λ(1520) π±
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantDToPhiPi{daughtersDsResonant.at(DecayChannelResonant::DsToPhiPi)};                    // Ds± → φ π± and D± → φ π±
    const std::array<int, NDaughtersResonant> arrPdgDaugResonantDToKstar0K{daughtersDsResonant.at(DecayChannelResonant::DsToKstar0K)};                // Ds± → anti-K*(892)0 K± and D± → anti-K*(892)0 K±

    // Reject particles from background events
    if (particle.fromBackgroundEvent() && rejectBackground) {
      rowMcMatchGen(flagChannelMain, origin, flagChannelResonant, -1);
      continue;
    }

    if (!pdgMothersCorrelBkg.empty()) {
      for (const auto& pdgMother : pdgMothersCorrelBkg) {
        if (std::abs(particle.pdgCode()) != pdgMother) {
          continue; // Skip if the particle PDG code does not match the mother PDG code
        }
        const auto finalStates = getDecayChannelsMain(pdgMother);
        constexpr int DepthMainMax = 2; // Depth for final state matching
        constexpr int DepthResoMax = 1; // Depth for resonant decay matching

        int depthMainMax = DepthMainMax;
        bool matched = false;
        if (pdgMother == Pdg::kDStar) {
          depthMainMax = DepthMainMax + 1; // D0 resonant decays are switched on
        }

        std::vector<int> arrAllDaughtersIndex;
        for (const auto& [channelMain, finalState] : finalStates) {
          if (finalState.size() == 5) { // o2-linter: disable=magic-number (partially reconstructed 3-prong decays from 5-prong decays)
            std::array<int, 5> arrPdgDaughtersMain5Prongs = std::array{finalState[0], finalState[1], finalState[2], finalState[3], finalState[4]};
            o2::hf_decay::flipPdgSign(particle.pdgCode(), +kPi0, arrPdgDaughtersMain5Prongs);
            RecoDecay::getDaughters<false>(particle, &arrAllDaughtersIndex, arrPdgDaughtersMain5Prongs, depthMainMax);
            matched = RecoDecay::isMatchedMCGen(mcParticles, particle, pdgMother, arrPdgDaughtersMain5Prongs, true, &sign, -1);
          } else if (finalState.size() == 4) { // o2-linter: disable=magic-number (partially reconstructed 3-prong decays from 4-prong decays)
            std::array<int, 4> arrPdgDaughtersMain4Prongs = std::array{finalState[0], finalState[1], finalState[2], finalState[3]};
            o2::hf_decay::flipPdgSign(particle.pdgCode(), +kPi0, arrPdgDaughtersMain4Prongs);
            RecoDecay::getDaughters<false>(particle, &arrAllDaughtersIndex, arrPdgDaughtersMain4Prongs, depthMainMax);
            matched = RecoDecay::isMatchedMCGen(mcParticles, particle, pdgMother, arrPdgDaughtersMain4Prongs, true, &sign, -1);
          } else if (finalState.size() == 3) { // o2-linter: disable=magic-number (fully reconstructed 3-prong decays)
            std::array<int, 3> arrPdgDaughtersMain3Prongs = std::array{finalState[0], finalState[1], finalState[2]};
            RecoDecay::getDaughters<false>(particle, &arrAllDaughtersIndex, arrPdgDaughtersMain3Prongs, depthMainMax);
            matched = RecoDecay::isMatchedMCGen(mcParticles, particle, pdgMother, arrPdgDaughtersMain3Prongs, true, &sign, depthMainMax);
          } else {
            LOG(fatal) << "Final state size not supported: " << finalState.size();
            return;
          }
          if (matched) {
            flagChannelMain = sign * channelMain;
            // Flag the resonant decay channel
            std::vector<int> arrResoDaughIndex = {};
            if (std::abs(pdgMother) == Pdg::kDStar) {
              std::vector<int> arrResoDaughIndexDStar = {};
              RecoDecay::getDaughters(particle, &arrResoDaughIndexDStar, std::array{0}, DepthResoMax);
              for (const int iDaug : arrResoDaughIndexDStar) {
                auto daughDstar = mcParticles.rawIteratorAt(iDaug);
                if (std::abs(daughDstar.pdgCode()) == Pdg::kD0 || std::abs(daughDstar.pdgCode()) == Pdg::kDPlus) {
                  RecoDecay::getDaughters(daughDstar, &arrResoDaughIndex, std::array{0}, DepthResoMax);
                  break;
                }
              }
            } else {
              RecoDecay::getDaughters(particle, &arrResoDaughIndex, std::array{0}, DepthResoMax);
            }
            std::array<int, NDaughtersResonant> arrPdgDaughters = {};
            if (arrResoDaughIndex.size() == NDaughtersResonant) {
              for (auto iProng = 0u; iProng < NDaughtersResonant; ++iProng) {
                auto daughI = mcParticles.rawIteratorAt(arrResoDaughIndex[iProng]);
                arrPdgDaughters[iProng] = daughI.pdgCode();
              }
              flagChannelResonant = o2::hf_decay::getDecayChannelResonant(pdgMother, arrPdgDaughters);
            }
            break; // Exit loop if a match is found
          }
        }
        if (matched) {
          break; // Exit loop if a match is found
        }
      }
    } else {

      // D± → π± K∓ π±
      if (flagChannelMain == 0) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flagChannelMain = sign * DecayChannelMain::DplusToPiKPi;
        }
      }

      // Ds± → K± K∓ π± and D± → K± K∓ π±
      if (flagChannelMain == 0) {
        bool isDplus = false;
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
          // TODO: move to different and explicit flags
          flagChannelMain = sign * DecayChannelMain::DsToPiKK;
        } else if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
          // TODO: move to different and explicit flags
          flagChannelMain = sign * DecayChannelMain::DplusToPiKK;
          isDplus = true;
        }
        if (flagChannelMain != 0) {
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
          if (arrDaughIndex.size() == NDaughtersResonant) {
            for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
              auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
              arrPdgDaugResonant[iProng] = std::abs(daughI.pdgCode());
            }
            if ((arrPdgDaugResonant[0] == arrPdgDaugResonantDToPhiPi[0] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToPhiPi[1]) || (arrPdgDaugResonant[0] == arrPdgDaugResonantDToPhiPi[1] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToPhiPi[0])) {
              flagChannelResonant = isDplus ? DecayChannelResonant::DplusToPhiPi : DecayChannelResonant::DsToPhiPi;
            } else if ((arrPdgDaugResonant[0] == arrPdgDaugResonantDToKstar0K[0] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToKstar0K[1]) || (arrPdgDaugResonant[0] == arrPdgDaugResonantDToKstar0K[1] && arrPdgDaugResonant[1] == arrPdgDaugResonantDToKstar0K[0])) {
              flagChannelResonant = isDplus ? DecayChannelResonant::DplusToKstar0K : DecayChannelResonant::DsToKstar0K;
            }
          }
        }
      }

      // D*± → D0(bar) π±
      if (flagChannelMain == 0) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &sign, 2)) {
          flagChannelMain = sign * DecayChannelMain::DstarToPiKPi;
        }
      }

      // Λc± → p± K∓ π±
      if (flagChannelMain == 0) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flagChannelMain = sign * DecayChannelMain::LcToPKPi;

          // Flagging the different Λc± → p± K∓ π± decay channels
          RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
          if (arrDaughIndex.size() == NDaughtersResonant) {
            for (auto iProng = 0u; iProng < arrDaughIndex.size(); ++iProng) {
              auto daughI = mcParticles.rawIteratorAt(arrDaughIndex[iProng]);
              arrPdgDaugResonant[iProng] = std::abs(daughI.pdgCode());
            }
            if ((arrPdgDaugResonant[0] == arrPdgDaugResonantLcToPKstar0[0] && arrPdgDaugResonant[1] == arrPdgDaugResonantLcToPKstar0[1]) || (arrPdgDaugResonant[0] == arrPdgDaugResonantLcToPKstar0[1] && arrPdgDaugResonant[1] == arrPdgDaugResonantLcToPKstar0[0])) {
              flagChannelResonant = DecayChannelResonant::LcToPKstar0;
            } else if ((arrPdgDaugResonant[0] == arrPdgDaugResonantLcToDeltaplusplusK[0] && arrPdgDaugResonant[1] == arrPdgDaugResonantLcToDeltaplusplusK[1]) || (arrPdgDaugResonant[0] == arrPdgDaugResonantLcToDeltaplusplusK[1] && arrPdgDaugResonant[1] == arrPdgDaugResonantLcToDeltaplusplusK[0])) {
              flagChannelResonant = DecayChannelResonant::LcToDeltaplusplusK;
            } else if ((arrPdgDaugResonant[0] == arrPdgDaugResonantLcToL1520Pi[0] && arrPdgDaugResonant[1] == arrPdgDaugResonantLcToL1520Pi[1]) || (arrPdgDaugResonant[0] == arrPdgDaugResonantLcToL1520Pi[1] && arrPdgDaugResonant[1] == arrPdgDaugResonantLcToL1520Pi[0])) {
              flagChannelResonant = DecayChannelResonant::LcToL1520Pi;
            }
          }
        }
      }

      // Ξc± → p± K∓ π±
      if (flagChannelMain == 0) {
        if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          flagChannelMain = sign * DecayChannelMain::XicToPKPi;
        }
      }
    }

    // Check whether the particle is non-prompt (from a b quark).
    if (flagChannelMain != 0) {
      origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
    }
    if (origin == RecoDecay::OriginType::NonPrompt) {
      rowMcMatchGen(flagChannelMain, origin, flagChannelResonant, idxBhadMothers[0]);
    } else {
      rowMcMatchGen(flagChannelMain, origin, flagChannelResonant, -1);
    }
  }
}

template <typename TMcParticles, typename TCursor>
void fillMcMatchGenBplus(TMcParticles const& mcParticles, TCursor& rowMcMatchGen)
{
  using namespace o2::constants::physics;
  using namespace o2::hf_decay::hf_cand_beauty;

  // Match generated particles.
  for (const auto& particle : mcParticles) {
    int8_t flagChannelMain = 0;
    int8_t flagChannelReso = 0;
    int8_t origin = 0;
    int8_t signB = 0;
    int8_t signD0 = 0;
    int indexGenD0 = -1;

    // B± → D0bar(D0) π± → (K± π∓) π±
    std::vector<int> arrayDaughterB;
    if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kBPlus, std::array{-Pdg::kD0, +kPiPlus}, true, &signB, 1, &arrayDaughterB)) {
      // D0(bar) → π± K∓
      for (const auto iD : arrayDaughterB) { // o2-linter: disable=const-ref-in-for-loop (int values)
        auto candDaughterMC = mcParticles.rawIteratorAt(iD);
        if (std::abs(candDaughterMC.pdgCode()) == Pdg::kD0) {
          indexGenD0 = RecoDecay::isMatchedMCGen(mcParticles, candDaughterMC, Pdg::kD0, std::array{-kKPlus, +kPiPlus}, true, &signD0, 1);
        }
      }
      if (indexGenD0 > -1) {
        flagChannelMain = signB * DecayChannelMain::BplusToD0Pi;
      }
    }
    rowMcMatchGen(flagChannelMain, flagChannelReso, origin);
  } // B candidate
}

template <typename TMcParticles, typename TCursor>
void fillMcMatchGenB0(TMcParticles const& mcParticles, TCursor& rowMcMatchGen)
{
  using namespace o2::constants::physics;
  using namespace o2::hf_decay::hf_cand_beauty;

  // Match generated particles.
  for (const auto& particle : mcParticles) {
    int8_t flagChannelMain = 0;
    int8_t flagChannelReso = 0;
    int8_t origin = 0;
    int8_t sign = 0;
    // B0 → D- π+
    if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kB0, std::array{-static_cast<int>(Pdg::kDPlus), +kPiPlus}, true)) {
      // D- → π- K+ π-
      auto candDMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
      if (RecoDecay::isMatchedMCGen(mcParticles, candDMC, -static_cast<int>(Pdg::kDPlus), std::array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign)) {
        flagChannelMain = sign * DecayChannelMain::B0ToDminusPi;
      }
    }
    rowMcMatchGen(flagChannelMain, flagChannelReso, origin);
  } // gen
}

} // namespace hf_mc_gen

#endif // PWGHF_UTILS_UTILSMCGEN_H_
