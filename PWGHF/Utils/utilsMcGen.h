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
/// \brief utility functions for HF Mc gen. workflows
///
/// \author Nima Zardoshti, nima.zardoshti@cern.ch, CERN

#ifndef PWGHF_UTILS_UTILSMCGEN_H_
#define PWGHF_UTILS_UTILSMCGEN_H_

#include <vector>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"

namespace hf_mc_gen
{

template <typename T, typename U, typename V>
void fillMcMatchGen2Prong(T const& mcParticles, U const& mcParticlesPerMcColl, V& rowMcMatchGen, bool rejectBackground)
{
  using namespace o2::constants::physics;

  // Match generated particles.
  for (const auto& particle : mcParticlesPerMcColl) {
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t sign = 0;
    std::vector<int> idxBhadMothers{};
    // Reject particles from background events
    if (particle.fromBackgroundEvent() && rejectBackground) {
      rowMcMatchGen(flag, origin, -1);
      continue;
    }

    // D0(bar) → π± K∓
    if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign)) {
      flag = sign * (1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK);
    }

    // J/ψ → e+ e−
    if (flag == 0) {
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kJPsi, std::array{+kElectron, -kElectron}, true)) {
        flag = 1 << o2::aod::hf_cand_2prong::DecayType::JpsiToEE;
      }
    }

    // J/ψ → μ+ μ−
    if (flag == 0) {
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kJPsi, std::array{+kMuonPlus, -kMuonPlus}, true)) {
        flag = 1 << o2::aod::hf_cand_2prong::DecayType::JpsiToMuMu;
      }
    }

    // Check whether the particle is non-prompt (from a b quark).
    if (flag != 0) {
      origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
    }
    if (origin == RecoDecay::OriginType::NonPrompt) {
      rowMcMatchGen(flag, origin, idxBhadMothers[0]);
    } else {
      rowMcMatchGen(flag, origin, -1);
    }
  }
}

template <typename T, typename U, typename V>
void fillMcMatchGen3Prong(T const& mcParticles, U const& mcParticlesPerMcColl, V& rowMcMatchGen, bool rejectBackground)
{
  using namespace o2::constants::physics;

  // Match generated particles.
  for (const auto& particle : mcParticlesPerMcColl) {
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t channel = 0;
    int8_t sign = 0;
    std::vector<int> arrDaughIndex;
    std::vector<int> idxBhadMothers{};
    std::array<int, 2> arrPDGDaugh;
    std::array<int, 2> arrPDGResonant1 = {kProton, 313};      // Λc± → p± K*
    std::array<int, 2> arrPDGResonant2 = {2224, kKPlus};      // Λc± → Δ(1232)±± K∓
    std::array<int, 2> arrPDGResonant3 = {102134, kPiPlus};   // Λc± → Λ(1520) π±
    std::array<int, 2> arrPDGResonantDPhiPi = {333, kPiPlus}; // Ds± → Phi π± and D± → Phi π±
    std::array<int, 2> arrPDGResonantDKstarK = {313, kKPlus}; // Ds± → K*(892)0bar K± and D± → K*(892)0bar K±
    // Reject particles from background events
    if (particle.fromBackgroundEvent() && rejectBackground) {
      rowMcMatchGen(flag, origin, channel, -1);
      continue;
    }

    // D± → π± K∓ π±
    if (flag == 0) {
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
        flag = sign * (1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi);
      }
    }

    // Ds± → K± K∓ π± and D± → K± K∓ π±
    if (flag == 0) {
      bool isDplus = false;
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
        // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
        // TODO: move to different and explicit flags
        flag = sign * (1 << o2::aod::hf_cand_3prong::DecayType::DsToKKPi);
      } else if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDPlus, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2)) {
        // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
        // TODO: move to different and explicit flags
        flag = sign * (1 << o2::aod::hf_cand_3prong::DecayType::DsToKKPi);
        isDplus = true;
      }
      if (flag != 0) {
        RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
        if (arrDaughIndex.size() == 2) {
          for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) {
            auto daughJ = mcParticles.rawIteratorAt(arrDaughIndex[jProng]);
            arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
          }
          if ((arrPDGDaugh[0] == arrPDGResonantDPhiPi[0] && arrPDGDaugh[1] == arrPDGResonantDPhiPi[1]) || (arrPDGDaugh[0] == arrPDGResonantDPhiPi[1] && arrPDGDaugh[1] == arrPDGResonantDPhiPi[0])) {
            channel = isDplus ? o2::aod::hf_cand_3prong::DecayChannelDToKKPi::DplusToPhiPi : o2::aod::hf_cand_3prong::DecayChannelDToKKPi::DsToPhiPi;
          } else if ((arrPDGDaugh[0] == arrPDGResonantDKstarK[0] && arrPDGDaugh[1] == arrPDGResonantDKstarK[1]) || (arrPDGDaugh[0] == arrPDGResonantDKstarK[1] && arrPDGDaugh[1] == arrPDGResonantDKstarK[0])) {
            channel = isDplus ? o2::aod::hf_cand_3prong::DecayChannelDToKKPi::DplusToK0starK : o2::aod::hf_cand_3prong::DecayChannelDToKKPi::DsToK0starK;
          }
        }
      }
    }

    // D*± → D0(bar) π±
    if (flag == 0) {
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kDStar, std::array{+kPiPlus, +kPiPlus, -kKPlus}, true, &sign, 2)) {
        flag = sign * (1 << o2::aod::hf_cand_3prong::DstarToPiKPiBkg);
      }
    }

    // Λc± → p± K∓ π±
    if (flag == 0) {
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
        flag = sign * (1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi);

        // Flagging the different Λc± → p± K∓ π± decay channels
        RecoDecay::getDaughters(particle, &arrDaughIndex, std::array{0}, 1);
        if (arrDaughIndex.size() == 2) {
          for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) {
            auto daughJ = mcParticles.rawIteratorAt(arrDaughIndex[jProng]);
            arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
          }
          if ((arrPDGDaugh[0] == arrPDGResonant1[0] && arrPDGDaugh[1] == arrPDGResonant1[1]) || (arrPDGDaugh[0] == arrPDGResonant1[1] && arrPDGDaugh[1] == arrPDGResonant1[0])) {
            channel = 1;
          } else if ((arrPDGDaugh[0] == arrPDGResonant2[0] && arrPDGDaugh[1] == arrPDGResonant2[1]) || (arrPDGDaugh[0] == arrPDGResonant2[1] && arrPDGDaugh[1] == arrPDGResonant2[0])) {
            channel = 2;
          } else if ((arrPDGDaugh[0] == arrPDGResonant3[0] && arrPDGDaugh[1] == arrPDGResonant3[1]) || (arrPDGDaugh[0] == arrPDGResonant3[1] && arrPDGDaugh[1] == arrPDGResonant3[0])) {
            channel = 3;
          }
        }
      }
    }

    // Ξc± → p± K∓ π±
    if (flag == 0) {
      if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
        flag = sign * (1 << o2::aod::hf_cand_3prong::DecayType::XicToPKPi);
      }
    }

    // Check whether the particle is non-prompt (from a b quark).
    if (flag != 0) {
      origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);
    }
    if (origin == RecoDecay::OriginType::NonPrompt) {
      rowMcMatchGen(flag, origin, channel, idxBhadMothers[0]);
    } else {
      rowMcMatchGen(flag, origin, channel, -1);
    }
  }
}

template <typename T, typename U>
void fillMcMatchGenBplus(T const& mcParticles, U& rowMcMatchGen)
{
  using namespace o2::constants::physics;

  // Match generated particles.
  for (const auto& particle : mcParticles) {
    int8_t flag = 0;
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
        flag = signB * (1 << o2::aod::hf_cand_bplus::DecayType::BplusToD0Pi);
      }
    }
    rowMcMatchGen(flag, origin);
  } // B candidate
}

template <typename T, typename U>
void fillMcMatchGenB0(T const& mcParticles, U& rowMcMatchGen)
{
  using namespace o2::constants::physics;

  // Match generated particles.
  for (const auto& particle : mcParticles) {
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t sign = 0;
    // B0 → D- π+
    if (RecoDecay::isMatchedMCGen(mcParticles, particle, Pdg::kB0, std::array{-static_cast<int>(Pdg::kDPlus), +kPiPlus}, true)) {
      // D- → π- K+ π-
      auto candDMC = mcParticles.rawIteratorAt(particle.daughtersIds().front());
      if (RecoDecay::isMatchedMCGen(mcParticles, candDMC, -static_cast<int>(Pdg::kDPlus), std::array{-kPiPlus, +kKPlus, -kPiPlus}, true, &sign)) {
        flag = sign * BIT(o2::aod::hf_cand_b0::DecayType::B0ToDPi);
      }
    }
    rowMcMatchGen(flag, origin);
  } // gen
}

} // namespace hf_mc_gen

#endif // PWGHF_UTILS_UTILSMCGEN_H_
