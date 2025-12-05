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

/// \file utilsCorrelations.h
/// \brief Utilities for HFC analyses
/// \author Xu Wang <waxu@cern.ch>

#ifndef PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_
#define PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"

#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>

#include <TPDGCode.h>

#include <Rtypes.h>

#include <cmath>
#include <cstddef>

namespace o2::analysis::hf_correlations
{
enum Region {
  Default = 0,
  Toward,
  Away,
  Transverse
};

enum PairSign {
  SignNotDefined = 0,
  LcPosTrkPos,
  LcPosTrkNeg,
  LcNegTrkPos,
  LcNegTrkNeg
};

constexpr float PhiTowardMax{o2::constants::math::PIThird};
constexpr float PhiAwayMin{2.f * o2::constants::math::PIThird};
constexpr float PhiAwayMax{4.f * o2::constants::math::PIThird};

template <typename T>
Region getRegion(T const deltaPhi)
{
  if (std::abs(deltaPhi) < PhiTowardMax) {
    return Toward;
  }
  if (deltaPhi > PhiAwayMin && deltaPhi < PhiAwayMax) {
    return Away;
  }
  return Transverse;
}

// Pair Sign Calculation
template <typename TrgPt, typename TrkPt>
int signCalulation(TrgPt const& trigPt, TrkPt const& assocPt)
{
  int sign = 0;
  if (trigPt > 0. && assocPt > 0.) {
    sign = LcPosTrkPos;
  } else if (trigPt > 0. && assocPt < 0.) {
    sign = LcPosTrkNeg;
  } else if (trigPt < 0. && assocPt > 0.) {
    sign = LcNegTrkPos;
  } else if (trigPt < 0. && assocPt < 0.) {
    sign = LcNegTrkNeg;
  } else {
    sign = SignNotDefined;
  }
  return sign;
}

template <typename Atrack, typename SpeciesContainer, typename T1, typename T2>
bool passPIDSelection(Atrack const& track, SpeciesContainer const mPIDspecies,
                      T1 const maxTPC, T2 const maxTOF, double ptThreshold = 0.75, bool tofForced = false)
{
  // Ensure size consistency
  if (mPIDspecies.value.size() != maxTPC.value.size() || mPIDspecies.value.size() != maxTOF.value.size()) {
    LOGF(error, "Size of particle species and corresponding nSigma selection arrays should be the same");
    return false; // Early exit on error
  }

  for (size_t speciesIndex = 0; speciesIndex < mPIDspecies.value.size(); ++speciesIndex) {
    auto const& pid = mPIDspecies->at(speciesIndex);
    auto nSigmaTPC = o2::aod::pidutils::tpcNSigma(pid, track);

    if (tofForced && !track.hasTOF()) {
      return false;
    }

    if (speciesIndex == 0) { // First species logic
      if (std::abs(nSigmaTPC) > maxTPC->at(speciesIndex)) {
        return false; // TPC check failed
      }
      if (tofForced || (track.pt() > ptThreshold && track.hasTOF())) {
        auto nSigmaTOF = o2::aod::pidutils::tofNSigma(pid, track);
        if (std::abs(nSigmaTOF) > maxTOF->at(speciesIndex)) {
          return false; // TOF check failed
        }
      }
    } else {                                                // Other species logic
      if (std::abs(nSigmaTPC) < maxTPC->at(speciesIndex)) { // Check TPC nSigma  first
        if (track.hasTOF()) {
          auto nSigmaTOF = o2::aod::pidutils::tofNSigma(pid, track);
          if (std::abs(nSigmaTOF) < maxTOF->at(speciesIndex)) {
            return false; // Reject if both TPC and TOF are within thresholds
          }
        } else {
          return false; // Reject if only TPC is within threshold and TOF is unavailable
        }
      }
    }
  }
  return true; // Passed all checks
}

/// @brief Selects a candidate based on its PDG code, decay channel, and assigns the corresponding mass.
///
/// @tparam isScCandidate  Boolean template parameter:
///                   - `true` to check for Sigma_c candidates
///                   - `false` to check for Lambda_c candidates
/// @tparam McParticleType Type representing the MC particle, must provide `pdgCode()` and `flagMcMatchGen()`
///
/// @param[in] particle  MC particle whose PDG code and decay flag are evaluated
/// @param[out] massCand Mass of the matched candidate is set here, if a valid match is found
///
/// @return `true` if candidate matches expected PDG and decay flag, and mass is set; `false` otherwise
template <bool IsScCandidate, typename McParticleType>
bool matchCandAndMass(McParticleType const& particle, double& massCand)
{
  const auto pdgCand = std::abs(particle.pdgCode());
  const auto matchGenFlag = std::abs(particle.flagMcMatchGen());

  // Validate PDG code based on candidate type
  if (IsScCandidate) {
    if (!(pdgCand == o2::constants::physics::Pdg::kSigmaC0 ||
          pdgCand == o2::constants::physics::Pdg::kSigmaCPlusPlus ||
          pdgCand == o2::constants::physics::Pdg::kSigmaCStar0 ||
          pdgCand == o2::constants::physics::Pdg::kSigmaCStarPlusPlus)) {
      return false;
    }
  } else {
    if (pdgCand != o2::constants::physics::Pdg::kLambdaCPlus) {
      return false;
    }
  }

  // Map decay type to mass
  switch (matchGenFlag) {
    case o2::hf_decay::hf_cand_sigmac::DecayChannelMain::Sc0ToPKPiPi: {
      massCand = o2::constants::physics::MassSigmaC0;
      return true;
    }

    case o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStar0ToPKPiPi: {
      massCand = o2::constants::physics::MassSigmaCStar0;
      return true;
    }

    case o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScplusplusToPKPiPi: {
      massCand = o2::constants::physics::MassSigmaCPlusPlus;
      return true;
    }

    case o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStarPlusPlusToPKPiPi: {
      massCand = o2::constants::physics::MassSigmaCStarPlusPlus;
      return true;
    }

    case hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi: {
      massCand = o2::constants::physics::MassLambdaCPlus;
      return true;
    }

    default: {
      return false;
    }
  }
}

// ========= Find Leading Particle ==============
template <typename TTracks, typename T1> //// FIXME: 14 days
int findLeadingParticle(TTracks const& tracks, T1 const etaTrackMax)
{
  auto leadingParticle = tracks.begin();
  for (auto const& track : tracks) {
    if (std::abs(track.eta()) > etaTrackMax) {
      continue;
    }
    if (track.pt() > leadingParticle.pt()) {
      leadingParticle = track;
    }
  }
  return leadingParticle.globalIndex();
}

// ======= Find Leading Particle for McGen ============
template <typename TMcParticles, typename T1, typename T2>
int findLeadingParticleMcGen(TMcParticles const& mcParticles, T1 const etaTrackMax, T2 const ptTrackMin)
{
  auto leadingParticle = mcParticles.begin();
  for (auto const& mcParticle : mcParticles) {
    if (std::abs(mcParticle.eta()) > etaTrackMax) {
      continue;
    }
    if (mcParticle.pt() < ptTrackMin) {
      continue;
    }
    if ((std::abs(mcParticle.pdgCode()) != kElectron) && (std::abs(mcParticle.pdgCode()) != kMuonMinus) && (std::abs(mcParticle.pdgCode()) != kPiPlus) && (std::abs(mcParticle.pdgCode()) != kKPlus) && (std::abs(mcParticle.pdgCode()) != kProton)) {
      continue;
    }
    if (mcParticle.pt() > leadingParticle.pt()) {
      leadingParticle = mcParticle;
    }
  }
  return leadingParticle.globalIndex();
}
} // namespace o2::analysis::hf_correlations
#endif // PWGHF_HFC_UTILS_UTILSCORRELATIONS_H_
