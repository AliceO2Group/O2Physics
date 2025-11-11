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

/// \file JetHFUtilities.h
/// \brief Jet HF related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETHFUTILITIES_H_
#define PWGJE_CORE_JETHFUTILITIES_H_

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/DecayChannelsLegacy.h"
#include "PWGJE/DataModel/Jet.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <vector>

namespace jethfutilities
{

/**
 * returns true if the candidate is from a D0 table
 */
template <typename T>
constexpr bool isD0Candidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesD0Data::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesD0Data::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesD0MCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesD0MCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a D0 MC table
 */
template <typename T>
constexpr bool isD0McCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesD0MCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesD0MCP::filtered_iterator>;
}

/**
 * returns true if the table is a D0 table
 */
template <typename T>
constexpr bool isD0Table()
{
  return isD0Candidate<typename T::iterator>() || isD0Candidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a D0 MC table
 */
template <typename T>
constexpr bool isD0McTable()
{
  return isD0McCandidate<typename T::iterator>() || isD0McCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a D+ table
 */
template <typename T>
constexpr bool isDplusCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDplusData::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDplusData::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDplusMCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDplusMCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a D+ MC table
 */
template <typename T>
constexpr bool isDplusMcCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDplusMCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDplusMCP::filtered_iterator>;
}

/**
 * returns true if the table is a D+ table
 */
template <typename T>
constexpr bool isDplusTable()
{
  return isDplusCandidate<typename T::iterator>() || isDplusCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a D+ MC table
 */
template <typename T>
constexpr bool isDplusMcTable()
{
  return isDplusMcCandidate<typename T::iterator>() || isDplusMcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a Ds table
 */
template <typename T>
constexpr bool isDsCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDsData::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDsData::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDsMCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDsMCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a Ds MC table
 */
template <typename T>
constexpr bool isDsMcCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDsMCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDsMCP::filtered_iterator>;
}

/**
 * returns true if the table is a Ds table
 */
template <typename T>
constexpr bool isDsTable()
{
  return isDsCandidate<typename T::iterator>() || isDsCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a Ds MC table
 */
template <typename T>
constexpr bool isDsMcTable()
{
  return isDsMcCandidate<typename T::iterator>() || isDsMcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a D* table
 */
template <typename T>
constexpr bool isDstarCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDstarData::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDstarData::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDstarMCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDstarMCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a D* MC table
 */
template <typename T>
constexpr bool isDstarMcCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDstarMCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesDstarMCP::filtered_iterator>;
}

/**
 * returns true if the table is a D* table
 */
template <typename T>
constexpr bool isDstarTable()
{
  return isDstarCandidate<typename T::iterator>() || isDstarCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a D* MC table
 */
template <typename T>
constexpr bool isDstarMcTable()
{
  return isDstarMcCandidate<typename T::iterator>() || isDstarMcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a Lc table
 */
template <typename T>
constexpr bool isLcCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesLcData::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesLcData::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesLcMCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesLcMCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a Lc MC table
 */
template <typename T>
constexpr bool isLcMcCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesLcMCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesLcMCP::filtered_iterator>;
}

/**
 * returns true if the table is a Lc table
 */
template <typename T>
constexpr bool isLcTable()
{
  return isLcCandidate<typename T::iterator>() || isLcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a Lc MC table
 */
template <typename T>
constexpr bool isLcMcTable()
{
  return isLcMcCandidate<typename T::iterator>() || isLcMcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a B0 table
 */
template <typename T>
constexpr bool isB0Candidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesB0Data::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesB0Data::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesB0MCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesB0MCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a B0 MC table
 */
template <typename T>
constexpr bool isB0McCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesB0MCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesB0MCP::filtered_iterator>;
}

/**
 * returns true if the table is a B0 table
 */
template <typename T>
constexpr bool isB0Table()
{
  return isB0Candidate<typename T::iterator>() || isB0Candidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a B0 MC table
 */
template <typename T>
constexpr bool isB0McTable()
{
  return isB0McCandidate<typename T::iterator>() || isB0McCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a Bplus table
 */
template <typename T>
constexpr bool isBplusCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesBplusData::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesBplusData::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesBplusMCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesBplusMCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a Bplus MC table
 */
template <typename T>
constexpr bool isBplusMcCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesBplusMCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesBplusMCP::filtered_iterator>;
}

/**
 * returns true if the table is a Bplus table
 */
template <typename T>
constexpr bool isBplusTable()
{
  return isBplusCandidate<typename T::iterator>() || isBplusCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a Bplus MC table
 */
template <typename T>
constexpr bool isBplusMcTable()
{
  return isBplusMcCandidate<typename T::iterator>() || isBplusMcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a XicToXiPiPi table
 */
template <typename T>
constexpr bool isXicToXiPiPiCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesXicToXiPiPiData::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesXicToXiPiPiData::filtered_iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesXicToXiPiPiMCD::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesXicToXiPiPiMCD::filtered_iterator>;
}

/**
 * returns true if the particle is from a XicToXiPiPi MC table
 */
template <typename T>
constexpr bool isXicToXiPiPiMcCandidate()
{
  return std::is_same_v<std::decay_t<T>, o2::aod::CandidatesXicToXiPiPiMCP::iterator> || std::is_same_v<std::decay_t<T>, o2::aod::CandidatesXicToXiPiPiMCP::filtered_iterator>;
}

/**
 * returns true if the table is a XicToXiPiPi table
 */
template <typename T>
constexpr bool isXicToXiPiPiTable()
{
  return isXicToXiPiPiCandidate<typename T::iterator>() || isXicToXiPiPiCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a XicToXiPiPi MC table
 */
template <typename T>
constexpr bool isXicToXiPiPiMcTable()
{
  return isXicToXiPiPiMcCandidate<typename T::iterator>() || isXicToXiPiPiMcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is from a HF table
 * * @param candidate candidate that is being checked
 */
template <typename T>
constexpr bool isHFCandidate()
{
  if constexpr (isD0Candidate<T>()) {
    return true;
  } else if constexpr (isDplusCandidate<T>()) {
    return true;
  } else if constexpr (isDsCandidate<T>()) {
    return true;
  } else if constexpr (isDstarCandidate<T>()) {
    return true;
  } else if constexpr (isLcCandidate<T>()) {
    return true;
  } else if constexpr (isB0Candidate<T>()) {
    return true;
  } else if constexpr (isBplusCandidate<T>()) {
    return true;
  } else if constexpr (isXicToXiPiPiCandidate<T>()) {
    return true;
  } else {
    return false;
  }
}

/**
 * returns true if the candidate is from a HF MC candidate table
 * * @param candidate candidate that is being checked
 */
template <typename T>
constexpr bool isHFMcCandidate()
{
  if constexpr (isD0McCandidate<T>()) {
    return true;
  } else if constexpr (isDplusMcCandidate<T>()) {
    return true;
  } else if constexpr (isDsMcCandidate<T>()) {
    return true;
  } else if constexpr (isDstarMcCandidate<T>()) {
    return true;
  } else if constexpr (isLcMcCandidate<T>()) {
    return true;
  } else if constexpr (isB0McCandidate<T>()) {
    return true;
  } else if constexpr (isXicToXiPiPiMcCandidate<T>()) {
    return true;
  } else if constexpr (isDstarMcCandidate<T>()) {
    return true;
  } else {
    return false;
  }
}

/**
 * returns true if the table type is a HF table
 */
template <typename T>
constexpr bool isHFTable()
{
  if constexpr (isD0Candidate<typename T::iterator>() || isD0Candidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isDplusCandidate<typename T::iterator>() || isDplusCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isDsCandidate<typename T::iterator>() || isDsCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isDstarCandidate<typename T::iterator>() || isDstarCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isLcCandidate<typename T::iterator>() || isLcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isB0Candidate<typename T::iterator>() || isB0Candidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isBplusCandidate<typename T::iterator>() || isBplusCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isXicToXiPiPiCandidate<typename T::iterator>() || isXicToXiPiPiCandidate<typename T::filtered_iterator>()) {
    return true;
  } else {
    return false;
  }
}

/**
 * returns true if the table type is a HF MC table
 */
template <typename T>
constexpr bool isHFMcTable()
{
  if constexpr (isD0McCandidate<typename T::iterator>() || isD0McCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isDplusMcCandidate<typename T::iterator>() || isDplusMcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isDsMcCandidate<typename T::iterator>() || isDsMcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isDstarMcCandidate<typename T::iterator>() || isDstarMcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isLcMcCandidate<typename T::iterator>() || isLcMcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isB0McCandidate<typename T::iterator>() || isB0McCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isBplusMcCandidate<typename T::iterator>() || isBplusMcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isXicToXiPiPiMcCandidate<typename T::iterator>() || isXicToXiPiPiMcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else {
    return false;
  }
}

/**
 * returns true if the candidate is matched to a reconstructed level candidate with the correct decay
 * * @param candidate candidate that is being checked
 */
template <typename T>
constexpr bool isMatchedHFCandidate(T const& candidate)
{
  if constexpr (isD0Candidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDplusCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDsCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDstarCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isLcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isB0Candidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_beauty::DecayChannelMain::B0ToDminusPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isBplusCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::hf_decay::hf_cand_beauty::DecayChannelMain::BplusToD0Pi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isXicToXiPiPiCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isD0McCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDplusMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDsMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDstarMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isLcMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isB0McCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_beauty::DecayChannelMain::B0ToDminusPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isBplusMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::hf_decay::hf_cand_beauty::DecayChannelMain::BplusToD0Pi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isXicToXiPiPiMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == o2::aod::hf_cand_xic_to_xi_pi_pi::DecayType::XicToXiPiPi) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

/**
 * returns true if the track is a daughter of the HF candidate
 *
 * @param track track that is being checked
 * @param candidate HF candidate that is being checked
 */
template <typename T, typename U>
bool isHFDaughterTrack(T& track, U& candidate)
{
  if constexpr (isD0Candidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDplusCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDsCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDstarCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isLcCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isB0Candidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex() || candidate.prong3Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isBplusCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isXicToXiPiPiCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex() || candidate.prong3Id() == track.globalIndex() || candidate.prong4Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

/**
 * returns the index of the JMcParticle matched to the HF candidate
 *
 * @param candidate hf candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedHFParticleId(const T& candidate, const U& /*tracks*/, const V& /*particles*/)
{
  const auto candidateDaughterParticle = candidate.template prong1_as<U>().template mcParticle_as<V>();
  return candidateDaughterParticle.template mothers_first_as<V>().globalIndex(); // can we get the Id directly?
}

/**
 * returns the JMcParticle matched to the HF candidate
 *
 * @param candidate hf candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedHFParticle(const T& candidate, const U& /*tracks*/, const V& /*particles*/)
{
  const auto candidateDaughterParticle = candidate.template prong1_as<U>().template mcParticle_as<V>();
  return candidateDaughterParticle.template mothers_first_as<V>();
}

/**
 * returns a slice of the table depending on the index of the HF candidate
 *
 * @param candidate HF candidate that is being checked
 * @param table the table to be sliced
 */
template <typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename R, typename S>
auto slicedPerHFCandidate(T const& table, U const& candidate, V const& perD0Candidate, M const& perDplusCandidate, N const& perDsCandidate, O const& perDstarCandidate, P const& perLcCandidate, Q const& perB0Candidate, R const& perBplusCandidate, S const& perXicToXiPiPiCandidate)
{
  if constexpr (isD0Candidate<U>()) {
    return table.sliceBy(perD0Candidate, candidate.globalIndex());
  } else if constexpr (isDplusCandidate<U>()) {
    return table.sliceBy(perDplusCandidate, candidate.globalIndex());
  } else if constexpr (isDstarCandidate<U>()) {
    return table.sliceBy(perDstarCandidate, candidate.globalIndex());
  } else if constexpr (isDsCandidate<U>()) {
    return table.sliceBy(perDsCandidate, candidate.globalIndex());
  } else if constexpr (isLcCandidate<U>()) {
    return table.sliceBy(perLcCandidate, candidate.globalIndex());
  } else if constexpr (isB0Candidate<U>()) {
    return table.sliceBy(perB0Candidate, candidate.globalIndex());
  } else if constexpr (isBplusCandidate<U>()) {
    return table.sliceBy(perBplusCandidate, candidate.globalIndex());
  } else if constexpr (isXicToXiPiPiCandidate<U>()) {
    return table.sliceBy(perXicToXiPiPiCandidate, candidate.globalIndex());
  } else {
    return table;
  }
}

/**
 * returns a slice of the table depending on the index of the HF candidate
 *
 * @param HFTable HF table type
 * @param jet jet that is being sliced based on
 * @param table the table to be sliced
 */
template <typename HFTable, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename R, typename S>
auto slicedPerHFJet(T const& table, U const& jet, V const& perD0Jet, M const& perDplusJet, N const& perDsJet, O const& perDstarJet, P const& perLcJet, Q const& perB0Jet, R const& perBplusJet, S const& perXicToXiPiPiJet)
{
  if constexpr (isD0Table<HFTable>() || isD0McTable<HFTable>()) {
    return table.sliceBy(perD0Jet, jet.globalIndex());
  } else if constexpr (isDplusTable<HFTable>() || isDplusMcTable<HFTable>()) {
    return table.sliceBy(perDplusJet, jet.globalIndex());
  } else if constexpr (isDstarTable<HFTable>() || isDstarMcTable<HFTable>()) {
    return table.sliceBy(perDstarJet, jet.globalIndex());
  } else if constexpr (isDsTable<HFTable>() || isDsMcTable<HFTable>()) {
    return table.sliceBy(perDsJet, jet.globalIndex());
  } else if constexpr (isLcTable<HFTable>() || isLcMcTable<HFTable>()) {
    return table.sliceBy(perLcJet, jet.globalIndex());
  } else if constexpr (isB0Table<HFTable>() || isB0McTable<HFTable>()) {
    return table.sliceBy(perB0Jet, jet.globalIndex());
  } else if constexpr (isBplusTable<HFTable>() || isBplusMcTable<HFTable>()) {
    return table.sliceBy(perBplusJet, jet.globalIndex());
  } else if constexpr (isXicToXiPiPiTable<HFTable>() || isXicToXiPiPiMcTable<HFTable>()) {
    return table.sliceBy(perXicToXiPiPiJet, jet.globalIndex());
  } else {
    return table;
  }
}

/**
 * returns the HF collision Id of candidate based on type of HF candidate
 *
 * @param candidate HF candidate that is being checked
 */
template <typename T>
int getHFCandidateCollisionId(T const& candidate)
{
  return candidate.hfCollBaseId();
}

/**
 * returns the HF Mc collision Id of candidate based on type of HF candidate
 *
 * @param candidate HF candidate that is being checked
 */
template <typename T>
int getHFMcCandidateCollisionId(T const& candidate)
{
  return candidate.hfMcCollBaseId();
}

/**
 * returns the PDG of the candidate based on HF Table
 *
 * @param candidate HF candidate that is being checked
 */
template <typename T>
int getHFCandidatePDG(T const& /*candidate*/)
{
  if constexpr (isD0Candidate<T>() || isD0McCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kD0);
  } else if constexpr (isDplusCandidate<T>() || isDplusMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kDPlus);
  } else if constexpr (isDsCandidate<T>() || isDsMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kDS);
  } else if constexpr (isDstarCandidate<T>() || isDstarMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kDStar);
  } else if constexpr (isLcCandidate<T>() || isLcMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kLambdaCPlus);
  } else if constexpr (isB0Candidate<T>() || isB0McCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kB0);
  } else if constexpr (isBplusCandidate<T>() || isBplusMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kBPlus);
  } else if constexpr (isXicToXiPiPiCandidate<T>() || isXicToXiPiPiMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kXiCPlus);
  } else {
    return 0;
  }
}

/**
 * returns the PDG of the candidates in the table type
 */
template <typename T>
int getHFTablePDG()
{
  if constexpr (isD0Table<T>() || isD0McTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kD0);
  } else if constexpr (isDplusTable<T>() || isDplusMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kDPlus);
  } else if constexpr (isDsTable<T>() || isDsMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kDS);
  } else if constexpr (isDstarTable<T>() || isDstarMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kDStar);
  } else if constexpr (isLcTable<T>() || isLcMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kLambdaCPlus);
  } else if constexpr (isB0Table<T>() || isB0McTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kB0);
  } else if constexpr (isBplusTable<T>() || isBplusMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kBPlus);
  } else if constexpr (isXicToXiPiPiTable<T>() || isXicToXiPiPiMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kXiCPlus);
  } else {
    return 0;
  }
}

/**
 * returns the mass of the candidate based on HF Table
 *
 * @param candidate HF candidate that is being checked
 */
template <typename T>
float getHFCandidatePDGMass(T const& /*candidate*/)
{
  if constexpr (isD0Candidate<T>() || isD0McCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassD0);
  } else if constexpr (isDplusCandidate<T>() || isDplusMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassDPlus);
  } else if constexpr (isDsCandidate<T>() || isDsMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassDS);
  } else if constexpr (isDstarCandidate<T>() || isDstarMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassDStar);
  } else if constexpr (isLcCandidate<T>() || isLcMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassLambdaCPlus);
  } else if constexpr (isB0Candidate<T>() || isB0McCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassB0);
  } else if constexpr (isBplusCandidate<T>() || isBplusMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassBPlus);
  } else if constexpr (isXicToXiPiPiCandidate<T>() || isXicToXiPiPiMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassXiCPlus);
  } else {
    return -1.0;
  }
}

/**
 * returns the mass of the candidates in the table type
 *
 */
template <typename T>
float getHFTablePDGMass()
{
  if constexpr (isD0Table<T>() || isD0McTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassD0);
  } else if constexpr (isDplusTable<T>() || isDplusMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassDPlus);
  } else if constexpr (isDsTable<T>() || isDsMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassDS);
  } else if constexpr (isDstarTable<T>() || isDstarMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassDStar);
  } else if constexpr (isLcTable<T>() || isLcMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassLambdaCPlus);
  } else if constexpr (isB0Table<T>() || isB0McTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassB0);
  } else if constexpr (isBplusTable<T>() || isBplusMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassBPlus);
  } else if constexpr (isXicToXiPiPiTable<T>() || isXicToXiPiPiMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassXiCPlus);
  } else {
    return -1.0;
  }
}

/**
 * returns the mass of the candidate based on HF Table
 *
 * @param candidate HF candidate that is being checked
 */
template <typename T>
float getHFCandidateInvariantMass(T const& candidate)
{
  return candidate.m();
}

template <typename T, typename U>
void fillHFCollisionTable(T const& collision, U& HFCollisionTable)
{
  HFCollisionTable(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.centFT0A(), collision.centFT0C(), collision.centFT0M(), collision.centFV0A(), collision.multZeqNTracksPV());
}

template <typename T, typename U>
void fillHFMcCollisionTable(T const& mcCollision, U& HFMcCollisionTable)
{
  HFMcCollisionTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.centFT0M());
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N>
void fillD0CandidateTable(T const& candidate, U& D0ParTable, V& D0ParETable, M& D0MlTable, N& D0MCDTable)
{
  D0ParTable(
    candidate.chi2PCA(),
    candidate.cpa(),
    candidate.cpaXY(),
    candidate.decayLength(),
    candidate.decayLengthXY(),
    candidate.decayLengthNormalised(),
    candidate.decayLengthXYNormalised(),
    candidate.ptProng0(),
    candidate.ptProng1(),
    candidate.impactParameter0(),
    candidate.impactParameter1(),
    candidate.impactParameterNormalised0(),
    candidate.impactParameterNormalised1(),
    candidate.nSigTpcPiExpPi(),
    candidate.nSigTofPiExpPi(),
    candidate.nSigTpcTofPiExpPi(),
    candidate.nSigTpcKaExpPi(),
    candidate.nSigTofKaExpPi(),
    candidate.nSigTpcTofKaExpPi(),
    candidate.nSigTpcPiExpKa(),
    candidate.nSigTofPiExpKa(),
    candidate.nSigTpcTofPiExpKa(),
    candidate.nSigTpcKaExpKa(),
    candidate.nSigTofKaExpKa(),
    candidate.nSigTpcTofKaExpKa(),
    candidate.maxNormalisedDeltaIP(),
    candidate.impactParameterProduct());

  D0ParETable(
    candidate.xSecondaryVertex(),
    candidate.ySecondaryVertex(),
    candidate.zSecondaryVertex(),
    candidate.errorDecayLength(),
    candidate.errorDecayLengthXY(),
    candidate.kfTopolChi2OverNdf(),
    candidate.rSecondaryVertex(),
    candidate.pProng0(),
    candidate.pProng1(),
    candidate.pxProng0(),
    candidate.pyProng0(),
    candidate.pzProng0(),
    candidate.pxProng1(),
    candidate.pyProng1(),
    candidate.pzProng1(),
    candidate.errorImpactParameter0(),
    candidate.errorImpactParameter1(),
    candidate.cosThetaStar(),
    candidate.ct());

  std::vector<float> mlScoresVector;
  auto mlScoresSpan = candidate.mlScores();
  std::copy(mlScoresSpan.begin(), mlScoresSpan.end(), std::back_inserter(mlScoresVector));
  D0MlTable(mlScoresVector);

  if constexpr (isMc) {
    D0MCDTable(candidate.flagMcMatchRec(), candidate.originMcRec());
  }
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N>
void fillDplusCandidateTable(T const& candidate, U& DplusParTable, V& DplusParETable, M& DplusMlTable, N& DplusMCDTable)
{

  DplusParTable(
    candidate.chi2PCA(),
    candidate.nProngsContributorsPV(),
    candidate.cpa(),
    candidate.cpaXY(),
    candidate.decayLength(),
    candidate.decayLengthXY(),
    candidate.decayLengthNormalised(),
    candidate.decayLengthXYNormalised(),
    candidate.ptProng0(),
    candidate.ptProng1(),
    candidate.ptProng2(),
    candidate.impactParameter0(),
    candidate.impactParameter1(),
    candidate.impactParameter2(),
    candidate.impactParameterNormalised0(),
    candidate.impactParameterNormalised1(),
    candidate.impactParameterNormalised2(),
    candidate.nSigTpcPi0(),
    candidate.nSigTofPi0(),
    candidate.nSigTpcTofPi0(),
    candidate.nSigTpcKa1(),
    candidate.nSigTofKa1(),
    candidate.nSigTpcTofKa1(),
    candidate.nSigTpcPi2(),
    candidate.nSigTofPi2(),
    candidate.nSigTpcTofPi2());

  DplusParETable(
    candidate.xSecondaryVertex(),
    candidate.ySecondaryVertex(),
    candidate.zSecondaryVertex(),
    candidate.errorDecayLength(),
    candidate.errorDecayLengthXY(),
    candidate.rSecondaryVertex(),
    candidate.pProng0(),
    candidate.pProng1(),
    candidate.pProng2(),
    candidate.pxProng0(),
    candidate.pyProng0(),
    candidate.pzProng0(),
    candidate.pxProng1(),
    candidate.pyProng1(),
    candidate.pzProng1(),
    candidate.pxProng2(),
    candidate.pyProng2(),
    candidate.pzProng2(),
    candidate.errorImpactParameter0(),
    candidate.errorImpactParameter1(),
    candidate.errorImpactParameter2(),
    candidate.ct());

  std::vector<float> mlScoresVector;
  auto mlScoresSpan = candidate.mlScores();
  std::copy(mlScoresSpan.begin(), mlScoresSpan.end(), std::back_inserter(mlScoresVector));
  DplusMlTable(mlScoresVector);

  if constexpr (isMc) {
    DplusMCDTable(candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.isCandidateSwapped(), candidate.flagMcDecayChanRec());
  }
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N>
void fillDsCandidateTable(T const& candidate, U& DsParTable, V& DsParETable, M& DsMlTable, N& DsMCDTable)
{

  DsParTable(
    candidate.chi2PCA(),
    candidate.nProngsContributorsPV(),
    candidate.cpa(),
    candidate.cpaXY(),
    candidate.decayLength(),
    candidate.decayLengthXY(),
    candidate.decayLengthNormalised(),
    candidate.decayLengthXYNormalised(),
    candidate.ptProng0(),
    candidate.ptProng1(),
    candidate.ptProng2(),
    candidate.impactParameter0(),
    candidate.impactParameter1(),
    candidate.impactParameter2(),
    candidate.impactParameterNormalised0(),
    candidate.impactParameterNormalised1(),
    candidate.impactParameterNormalised2(),
    candidate.nSigTpcPi0(),
    candidate.nSigTpcKa0(),
    candidate.nSigTofPi0(),
    candidate.nSigTofKa0(),
    candidate.nSigTpcTofPi0(),
    candidate.nSigTpcTofKa0(),
    candidate.nSigTpcKa1(),
    candidate.nSigTofKa1(),
    candidate.nSigTpcTofKa1(),
    candidate.nSigTpcPi2(),
    candidate.nSigTpcKa2(),
    candidate.nSigTofPi2(),
    candidate.nSigTofKa2(),
    candidate.nSigTpcTofPi2(),
    candidate.nSigTpcTofKa2());

  DsParETable(
    candidate.xSecondaryVertex(),
    candidate.ySecondaryVertex(),
    candidate.zSecondaryVertex(),
    candidate.errorDecayLength(),
    candidate.errorDecayLengthXY(),
    candidate.rSecondaryVertex(),
    candidate.pProng0(),
    candidate.pProng1(),
    candidate.pProng2(),
    candidate.pxProng0(),
    candidate.pyProng0(),
    candidate.pzProng0(),
    candidate.pxProng1(),
    candidate.pyProng1(),
    candidate.pzProng1(),
    candidate.pxProng2(),
    candidate.pyProng2(),
    candidate.pzProng2(),
    candidate.errorImpactParameter0(),
    candidate.errorImpactParameter1(),
    candidate.errorImpactParameter2(),
    candidate.ct());

  std::vector<float> mlScoresVector;
  auto mlScoresSpan = candidate.mlScores();
  std::copy(mlScoresSpan.begin(), mlScoresSpan.end(), std::back_inserter(mlScoresVector));
  DsMlTable(mlScoresVector);

  if constexpr (isMc) {
    DsMCDTable(candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.isCandidateSwapped(), candidate.flagMcDecayChanRec());
  }
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N>
void fillDstarCandidateTable(T const& candidate, U& DstarParTable, V& DstarParDaughterTable, M& DstarMlTable, N& DstarMCDTable)
{

  DstarParTable(
    candidate.pxProng0(),
    candidate.pyProng0(),
    candidate.pzProng0(),
    candidate.pxProng1(),
    candidate.pyProng1(),
    candidate.pzProng1(),
    candidate.signProng1(),
    candidate.impactParameter1(),
    candidate.impactParameterNormalised1(),
    candidate.nSigTpcPi1(),
    candidate.nSigTofPi1(),
    candidate.nSigTpcTofPi1());

  DstarParDaughterTable(
    candidate.chi2PCACharm(),
    candidate.cpaCharm(),
    candidate.cpaXYCharm(),
    candidate.decayLengthCharm(),
    candidate.decayLengthXYCharm(),
    candidate.decayLengthNormalisedCharm(),
    candidate.decayLengthXYNormalisedCharm(),
    candidate.pxProng0Charm(),
    candidate.pyProng0Charm(),
    candidate.pzProng0Charm(),
    candidate.pxProng1Charm(),
    candidate.pyProng1Charm(),
    candidate.pzProng1Charm(),
    candidate.invMassCharm(),
    candidate.impactParameter0Charm(),
    candidate.impactParameter1Charm(),
    candidate.impactParameterNormalised0Charm(),
    candidate.impactParameterNormalised1Charm(),
    candidate.nSigTpcPi0Charm(),
    candidate.nSigTofPi0Charm(),
    candidate.nSigTpcTofPi0Charm(),
    candidate.nSigTpcKa0Charm(),
    candidate.nSigTofKa0Charm(),
    candidate.nSigTpcTofKa0Charm(),
    candidate.nSigTpcPi1Charm(),
    candidate.nSigTofPi1Charm(),
    candidate.nSigTpcTofPi1Charm(),
    candidate.nSigTpcKa1Charm(),
    candidate.nSigTofKa1Charm(),
    candidate.nSigTpcTofKa1Charm());

  std::vector<float> mlScoresVector;
  auto mlScoresSpan = candidate.mlScores();
  std::copy(mlScoresSpan.begin(), mlScoresSpan.end(), std::back_inserter(mlScoresVector));
  DstarMlTable(mlScoresVector);

  if constexpr (isMc) {
    DstarMCDTable(
      candidate.flagMcMatchRec(),
      candidate.flagMcMatchRecCharm(),
      candidate.originMcRec(),
      candidate.ptBhadMotherPart(),
      candidate.pdgBhadMotherPart(),
      candidate.nTracksDecayed());
  }
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N>
void fillLcCandidateTable(T const& candidate, U& LcParTable, V& LcParETable, M& LcMlTable, N& LcMCDTable)
{

  LcParTable(
    candidate.chi2PCA(),
    candidate.nProngsContributorsPV(),
    candidate.cpa(),
    candidate.cpaXY(),
    candidate.decayLength(),
    candidate.decayLengthXY(),
    candidate.decayLengthNormalised(),
    candidate.decayLengthXYNormalised(),
    candidate.ptProng0(),
    candidate.ptProng1(),
    candidate.ptProng2(),
    candidate.impactParameter0(),
    candidate.impactParameter1(),
    candidate.impactParameter2(),
    candidate.impactParameterNormalised0(),
    candidate.impactParameterNormalised1(),
    candidate.impactParameterNormalised2(),
    candidate.nSigTpcPi0(),
    candidate.nSigTpcPr0(),
    candidate.nSigTofPi0(),
    candidate.nSigTofPr0(),
    candidate.nSigTpcTofPi0(),
    candidate.nSigTpcTofPr0(),
    candidate.nSigTpcKa1(),
    candidate.nSigTofKa1(),
    candidate.nSigTpcTofKa1(),
    candidate.nSigTpcPi2(),
    candidate.nSigTpcPr2(),
    candidate.nSigTofPi2(),
    candidate.nSigTofPr2(),
    candidate.nSigTpcTofPi2(),
    candidate.nSigTpcTofPr2());

  LcParETable(
    candidate.xSecondaryVertex(),
    candidate.ySecondaryVertex(),
    candidate.zSecondaryVertex(),
    candidate.errorDecayLength(),
    candidate.errorDecayLengthXY(),
    candidate.rSecondaryVertex(),
    candidate.pProng0(),
    candidate.pProng1(),
    candidate.pProng2(),
    candidate.pxProng0(),
    candidate.pyProng0(),
    candidate.pzProng0(),
    candidate.pxProng1(),
    candidate.pyProng1(),
    candidate.pzProng1(),
    candidate.pxProng2(),
    candidate.pyProng2(),
    candidate.pzProng2(),
    candidate.errorImpactParameter0(),
    candidate.errorImpactParameter1(),
    candidate.errorImpactParameter2(),
    candidate.ct());

  std::vector<float> mlScoresVector;
  auto mlScoresSpan = candidate.mlScores();
  std::copy(mlScoresSpan.begin(), mlScoresSpan.end(), std::back_inserter(mlScoresVector));
  LcMlTable(mlScoresVector);

  if constexpr (isMc) {
    LcMCDTable(candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.isCandidateSwapped());
  }
}

// need to update this
template <bool isMc, typename T, typename U, typename V, typename M, typename N, typename O, typename P>
void fillB0CandidateTable(T const& candidate, U& B0ParTable, V& B0ParETable, M& B0ParD0Table, N& B0MlTable, O& B0MlD0Table, P& B0MCDTable)
{

  B0ParTable(
    candidate.chi2PCA(),
    candidate.cpa(),
    candidate.cpaXY(),
    candidate.decayLength(),
    candidate.decayLengthXY(),
    candidate.decayLengthNormalised(),
    candidate.decayLengthXYNormalised(),
    candidate.ptProng0(),
    candidate.ptProng1(),
    candidate.impactParameter0(),
    candidate.impactParameter1(),
    candidate.impactParameterNormalised0(),
    candidate.impactParameterNormalised1(),
    candidate.nSigTpcPiExpPi(),
    candidate.nSigTofPiExpPi(),
    candidate.nSigTpcTofPiExpPi(),
    candidate.nSigTpcKaExpPi(),
    candidate.nSigTofKaExpPi(),
    candidate.nSigTpcTofKaExpPi(),
    candidate.maxNormalisedDeltaIP(),
    candidate.impactParameterProduct());

  B0ParETable(
    candidate.xSecondaryVertex(),
    candidate.ySecondaryVertex(),
    candidate.zSecondaryVertex(),
    candidate.errorDecayLength(),
    candidate.errorDecayLengthXY(),
    candidate.rSecondaryVertex(),
    candidate.pProng1(),
    candidate.pxProng1(),
    candidate.pyProng1(),
    candidate.pzProng1(),
    candidate.errorImpactParameter1(),
    candidate.cosThetaStar(),
    candidate.ct());

  B0ParD0Table(
    candidate.chi2PCACharm(),
    candidate.nProngsContributorsPVCharm(),
    candidate.cpaCharm(),
    candidate.cpaXYCharm(),
    candidate.decayLengthCharm(),
    candidate.decayLengthXYCharm(),
    candidate.decayLengthNormalisedCharm(),
    candidate.decayLengthXYNormalisedCharm(),
    candidate.ptProng0Charm(),
    candidate.ptProng1Charm(),
    candidate.ptProng2Charm(),
    candidate.impactParameter0Charm(),
    candidate.impactParameter1Charm(),
    candidate.impactParameter2Charm(),
    candidate.impactParameterNormalised0Charm(),
    candidate.impactParameterNormalised1Charm(),
    candidate.impactParameterNormalised2Charm(),
    candidate.nSigTpcPi0Charm(),
    candidate.nSigTofPi0Charm(),
    candidate.nSigTpcTofPi0Charm(),
    candidate.nSigTpcKa1Charm(),
    candidate.nSigTofKa1Charm(),
    candidate.nSigTpcTofKa1Charm(),
    candidate.nSigTpcPi2Charm(),
    candidate.nSigTofPi2Charm(),
    candidate.nSigTpcTofPi2Charm());

  // B0SelectionFlagTable(candidate.candidateSelFlag());

  B0MlTable(candidate.mlScoreSig());

  std::vector<float> mlScoresCharmVector;
  auto mlScoresCharmSpan = candidate.mlScoresCharm();
  std::copy(mlScoresCharmSpan.begin(), mlScoresCharmSpan.end(), std::back_inserter(mlScoresCharmVector));
  B0MlD0Table(mlScoresCharmVector);

  if constexpr (isMc) {
    B0MCDTable(candidate.flagMcMatchRec(), candidate.originMcRec());
  }
}

// need to update this
template <bool isMc, typename T, typename U, typename V, typename M, typename N, typename O, typename P>
void fillBplusCandidateTable(T const& candidate, U& BplusParTable, V& BplusParETable, M& BplusParD0Table, N& BplusMlTable, O& BplusMlD0Table, P& BplusMCDTable)
{

  BplusParTable(
    candidate.chi2PCA(),
    candidate.cpa(),
    candidate.cpaXY(),
    candidate.decayLength(),
    candidate.decayLengthXY(),
    candidate.decayLengthNormalised(),
    candidate.decayLengthXYNormalised(),
    candidate.ptProng0(),
    candidate.ptProng1(),
    candidate.impactParameter0(),
    candidate.impactParameter1(),
    candidate.impactParameterNormalised0(),
    candidate.impactParameterNormalised1(),
    candidate.nSigTpcPiExpPi(),
    candidate.nSigTofPiExpPi(),
    candidate.nSigTpcTofPiExpPi(),
    candidate.nSigTpcKaExpPi(),
    candidate.nSigTofKaExpPi(),
    candidate.nSigTpcTofKaExpPi(),
    candidate.maxNormalisedDeltaIP(),
    candidate.impactParameterProduct());

  BplusParETable(
    candidate.xSecondaryVertex(),
    candidate.ySecondaryVertex(),
    candidate.zSecondaryVertex(),
    candidate.errorDecayLength(),
    candidate.errorDecayLengthXY(),
    candidate.rSecondaryVertex(),
    candidate.pProng1(),
    candidate.pxProng1(),
    candidate.pyProng1(),
    candidate.pzProng1(),
    candidate.errorImpactParameter1(),
    candidate.cosThetaStar(),
    candidate.ct());

  BplusParD0Table(
    candidate.cpaCharm(),
    candidate.decayLengthCharm(),
    candidate.impactParameter0Charm(),
    candidate.impactParameter1Charm(),
    candidate.impactParameterProductCharm(),
    candidate.nSigTpcPiExpPiCharm(),
    candidate.nSigTofPiExpPiCharm(),
    candidate.nSigTpcTofPiExpPiCharm(),
    candidate.nSigTpcKaExpPiCharm(),
    candidate.nSigTofKaExpPiCharm(),
    candidate.nSigTpcTofKaExpPiCharm(),
    candidate.nSigTpcPiExpKaCharm(),
    candidate.nSigTofPiExpKaCharm(),
    candidate.nSigTpcTofPiExpKaCharm(),
    candidate.nSigTpcKaExpKaCharm(),
    candidate.nSigTofKaExpKaCharm(),
    candidate.nSigTpcTofKaExpKaCharm());

  // BplusSelectionFlagTable(candidate.candidateSelFlag());

  BplusMlTable(candidate.mlScoreSig());

  std::vector<float> mlScoresCharmVector;
  auto mlScoresCharmSpan = candidate.mlScoresCharm();
  std::copy(mlScoresCharmSpan.begin(), mlScoresCharmSpan.end(), std::back_inserter(mlScoresCharmVector));
  BplusMlD0Table(mlScoresCharmVector);

  if constexpr (isMc) {
    BplusMCDTable(candidate.flagMcMatchRec(), candidate.originMcRec());
  }
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N>
void fillXicToXiPiPiCandidateTable(T const& candidate, U& XicToXiPiPiParTable, V& XicToXiPiPiParETable, M& XicToXiPiPiMlTable, N& XicToXiPiPiMCDTable)
{

  XicToXiPiPiParTable(
    candidate.sign(),
    candidate.ptProngXi(),
    candidate.ptProngPi0(),
    candidate.ptProngPi1(),
    candidate.invMassXi(),
    candidate.invMassLambda(),
    candidate.invMassXiPi0(),
    candidate.invMassXiPi1(),
    candidate.chi2PCA(),
    candidate.ct(),
    candidate.decayLength(),
    candidate.decayLengthXY(),
    candidate.decayLengthNormalised(),
    candidate.decayLengthXYNormalised(),
    candidate.cpa(),
    candidate.cpaXY(),
    candidate.cpaXi(),
    candidate.cpaXYXi(),
    candidate.cpaLambda(),
    candidate.cpaXYLambda(),
    candidate.impactParameterXi(),
    candidate.impactParameterNormalisedXi(),
    candidate.impactParameterPi0(),
    candidate.impactParameterNormalisedPi0(),
    candidate.impactParameterPi1(),
    candidate.impactParameterNormalisedPi1(),
    candidate.maxNormalisedDeltaIP());

  XicToXiPiPiParETable(
    candidate.cpaLambdaToXi(),
    candidate.cpaXYLambdaToXi(),
    candidate.pProngPi0(),
    candidate.pProngPi1(),
    candidate.pBachelorPi(),
    candidate.pPiFromLambda(),
    candidate.pPrFromLambda(),
    candidate.dcaXiDaughters(),
    candidate.dcaV0Daughters(),
    candidate.dcaPosToPV(),
    candidate.dcaNegToPV(),
    candidate.dcaBachelorToPV(),
    candidate.dcaXYCascToPV(),
    candidate.dcaZCascToPV(),
    candidate.nSigTpcPiFromXicPlus0(),
    candidate.nSigTpcPiFromXicPlus1(),
    candidate.nSigTpcBachelorPi(),
    candidate.nSigTpcPiFromLambda(),
    candidate.nSigTpcPrFromLambda(),
    candidate.nSigTofPiFromXicPlus0(),
    candidate.nSigTofPiFromXicPlus1(),
    candidate.nSigTofBachelorPi(),
    candidate.nSigTofPiFromLambda(),
    candidate.nSigTofPrFromLambda());

  std::vector<float> mlScoresVector;
  auto mlScoresSpan = candidate.mlScores();
  std::copy(mlScoresSpan.begin(), mlScoresSpan.end(), std::back_inserter(mlScoresVector));
  XicToXiPiPiMlTable(mlScoresVector);

  if constexpr (isMc) {
    XicToXiPiPiMCDTable(candidate.flagMcMatchRec(), candidate.originMcRec());
  }
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename S>
void fillHFCandidateTable(T const& candidate, int32_t collisionIndex, U& HFBaseTable, V& HFParTable, M& HFParETable, N& HFParDaughterTable, O& HFSelectionFlagTable, P& HFMlTable, Q& HFMlDaughterTable, S& HFMCDTable)
{
  HFBaseTable(collisionIndex, candidate.pt(), candidate.eta(), candidate.phi(), candidate.m(), candidate.y());
  HFSelectionFlagTable(candidate.candidateSelFlag());

  if constexpr (isD0Candidate<T>()) {
    fillD0CandidateTable<isMc>(candidate, HFParTable, HFParETable, HFMlTable, HFMCDTable);
  }
  if constexpr (isDplusCandidate<T>()) {
    fillDplusCandidateTable<isMc>(candidate, HFParTable, HFParETable, HFMlTable, HFMCDTable);
  }
  if constexpr (isDsCandidate<T>()) {
    fillDsCandidateTable<isMc>(candidate, HFParTable, HFParETable, HFMlTable, HFMCDTable);
  }
  if constexpr (isDstarCandidate<T>()) {
    fillDstarCandidateTable<isMc>(candidate, HFParTable, HFParDaughterTable, HFMlTable, HFMCDTable);
  }
  if constexpr (isLcCandidate<T>()) {
    fillLcCandidateTable<isMc>(candidate, HFParTable, HFParETable, HFMlTable, HFMCDTable);
  }
  if constexpr (isB0Candidate<T>()) {
    fillB0CandidateTable<isMc>(candidate, HFParTable, HFParETable, HFParDaughterTable, HFMlTable, HFMlDaughterTable, HFMCDTable);
  }
  if constexpr (isBplusCandidate<T>()) {
    fillBplusCandidateTable<isMc>(candidate, HFParTable, HFParETable, HFParDaughterTable, HFMlTable, HFMlDaughterTable, HFMCDTable);
  }
  if constexpr (isXicToXiPiPiCandidate<T>()) {
    fillXicToXiPiPiCandidateTable<isMc>(candidate, HFParTable, HFParETable, HFMlTable, HFMCDTable);
  }
}

template <typename T, typename U>
void fillHFCandidateMcTable(T const& candidate, int32_t mcCollisionIndex, U& BaseMcTable)
{
  BaseMcTable(mcCollisionIndex, candidate.pt(), candidate.eta(), candidate.phi(), candidate.y(), candidate.flagMcMatchGen(), candidate.originMcGen());
}

}; // namespace jethfutilities

#endif // PWGJE_CORE_JETHFUTILITIES_H_
