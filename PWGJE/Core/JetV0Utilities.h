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

/// \file JetV0Utilities.h
/// \brief Jet V0 related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETV0UTILITIES_H_
#define PWGJE_CORE_JETV0UTILITIES_H_

#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>

#include <TPDGCode.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <string>
#include <type_traits>

namespace jetv0utilities
{

/**
 * returns true if the candidate is from a V0 table
 */
template <typename T>
constexpr bool isV0Candidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesV0Data::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesV0Data::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesV0MCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesV0MCD::filtered_iterator>;
}

/**
 * returns true if the candidate is from a MC V0 table
 */
template <typename T>
constexpr bool isV0McCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesV0MCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesV0MCP::filtered_iterator>;
}

/**
 * returns true if the table is a V0 table
 */
template <typename T>
constexpr bool isV0Table()
{

  return isV0Candidate<typename T::iterator>() || isV0Candidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a V0 MC table
 */
template <typename T>
constexpr bool isV0McTable()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesV0MCP> || std::is_same_v<std::decay_t<T>, o2::soa::Filtered<o2::aod::CandidatesV0MCP>>; // note not optimal way but needed for jetfindingutilities::analyseParticles()
}

/**
 * returns true if the track is a daughter of the V0 candidate
 *
 * @param track track that is being checked
 * @param candidate V0 candidate that is being checked
 */
template <typename T, typename U>
bool isV0DaughterTrack(T& track, U& candidate)
{
  if constexpr (isV0Candidate<U>()) {
    if (candidate.posTrackId() == track.globalIndex() || candidate.negTrackId() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

/**
 * returns the index of the JMcParticle matched to the V0 candidate
 *
 * @param candidate V0 candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedV0ParticleId(const T& candidate, const U& /*tracks*/, const V& /*particles*/)
{
  const auto candidateDaughterParticle = candidate.template posTrack_as<U>().template mcParticle_as<V>();
  return candidateDaughterParticle.template mothers_first_as<V>().globalIndex(); // can we get the Id directly?
}

/**
 * returns the JMcParticle matched to the V0 candidate
 *
 * @param candidate v0 candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedV0Particle(const T& candidate, const U& /*tracks*/, const V& /*particles*/)
{
  const auto candidateDaughterParticle = candidate.template posTrack_as<U>().template mcParticle_as<V>();
  return candidateDaughterParticle.template mothers_first_as<V>();
}

/**
 * returns a slice of the table depending on the index of the V0 candidate
 *
 * @param candidate v0 candidate that is being checked
 * @param table the table to be sliced
 */
template <typename T, typename U, typename V>
auto slicedPerV0Candidate(T const& table, U const& candidate, V const& perV0Candidate)
{
  if constexpr (isV0Candidate<U>()) {
    return table.sliceBy(perV0Candidate, candidate.globalIndex());
  } else {
    return table;
  }
}

template <typename T, typename U>
bool isV0Particle(T const& particles, U const& particle, bool v0ChargedDecaysOnly)
{
  if (v0ChargedDecaysOnly) {
    return RecoDecay::isMatchedMCGen(particles, particle, +kK0Short, std::array{+kPiPlus, -kPiPlus}, true) || RecoDecay::isMatchedMCGen(particles, particle, +kLambda0, std::array{+kProton, -kPiPlus}, true) || RecoDecay::isMatchedMCGen(particles, particle, -kLambda0, std::array{-kProton, +kPiPlus}, true);
  } else {
    return RecoDecay::isMatchedMCGen(particles, particle, +kK0Short, std::array{+kPiPlus, -kPiPlus}, true) || RecoDecay::isMatchedMCGen(particles, particle, +kK0Short, std::array{+kPi0, +kPi0}, true) || RecoDecay::isMatchedMCGen(particles, particle, +kLambda0, std::array{+kProton, -kPiPlus}, true) || RecoDecay::isMatchedMCGen(particles, particle, +kLambda0, std::array{+kNeutron, +kPi0}, true) || RecoDecay::isMatchedMCGen(particles, particle, -kLambda0, std::array{-kProton, +kPiPlus}, true) || RecoDecay::isMatchedMCGen(particles, particle, -kLambda0, std::array{+kNeutron, +kPi0}, true);
  }
}

enum JV0ParticleDecays {
  K0sToPiPi = 0,
  LambdaToPPi = 1
};

template <typename T>
bool selectV0ParticleDecay(T const& v0Particle, int v0ParticleDecaySelection = -1)
{
  if (v0ParticleDecaySelection == -1) {
    return true;
  }
  return (v0Particle.decayFlag() & (1 << v0ParticleDecaySelection));
}

int initialiseV0ParticleDecaySelection(const std::string& v0ParticleDecaySelection)
{
  if (v0ParticleDecaySelection == "K0sToPiPi") {
    return JV0ParticleDecays::K0sToPiPi;
  } else if (v0ParticleDecaySelection == "LambdaToPPi") {
    return JV0ParticleDecays::LambdaToPPi;
  }
  return -1;
}

template <typename T, typename U>
uint8_t setV0ParticleDecayBit(T const& particles, U const& particle)
{
  uint8_t bit = 0;
  if (RecoDecay::isMatchedMCGen(particles, particle, +kK0Short, std::array{+kPiPlus, -kPiPlus}, true)) {
    SETBIT(bit, JV0ParticleDecays::K0sToPiPi);
  }
  if (RecoDecay::isMatchedMCGen(particles, particle, +kLambda0, std::array{+kProton, -kPiPlus}, true) || RecoDecay::isMatchedMCGen(particles, particle, -kLambda0, std::array{-kProton, +kPiPlus}, true)) {
    SETBIT(bit, JV0ParticleDecays::LambdaToPPi);
  }
  return bit;
}

}; // namespace jetv0utilities

#endif // PWGJE_CORE_JETV0UTILITIES_H_
