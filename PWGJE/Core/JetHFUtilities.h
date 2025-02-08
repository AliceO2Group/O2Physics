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

#include <array>
#include <vector>
#include <string>
#include <optional>
#include <algorithm>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Framework/Logger.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

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
  } else if constexpr (isLcCandidate<T>()) {
    return true;
  } else if constexpr (isBplusCandidate<T>()) {
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
  } else if constexpr (isLcMcCandidate<T>()) {
    return true;
  } else if constexpr (isBplusMcCandidate<T>()) {
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
  } else if constexpr (isLcCandidate<typename T::iterator>() || isLcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isBplusCandidate<typename T::iterator>() || isBplusCandidate<typename T::filtered_iterator>()) {
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
  } else if constexpr (isLcMcCandidate<typename T::iterator>() || isLcMcCandidate<typename T::filtered_iterator>()) {
    return true;
  } else if constexpr (isBplusMcCandidate<typename T::iterator>() || isBplusMcCandidate<typename T::filtered_iterator>()) {
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
    if (std::abs(candidate.flagMcMatchRec()) == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDplusCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == 1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isLcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == 1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isBplusCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchRec()) == 1 << o2::aod::hf_cand_bplus::DecayType::BplusToD0Pi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isD0McCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isDplusMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == 1 << o2::aod::hf_cand_3prong::DecayType::DplusToPiKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isLcMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == 1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi) {
      return true;
    } else {
      return false;
    }
  } else if constexpr (isBplusMcCandidate<T>()) {
    if (std::abs(candidate.flagMcMatchGen()) == 1 << o2::aod::hf_cand_bplus::DecayType::BplusToD0Pi) {
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
 * @param tracks the track table
 */
template <typename T, typename U, typename V>
bool isHFDaughterTrack(T& track, U& candidate, V const& /*tracks*/)
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
  } else if constexpr (isLcCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex() || candidate.prong2Id() == track.globalIndex()) {
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
template <typename T, typename U, typename V, typename M, typename N, typename O>
auto slicedPerHFCandidate(T const& table, U const& candidate, V const& perD0Candidate, M const& perDplusCandidate, N const& perLcCandidate, O const& perBplusCandidate)
{
  if constexpr (isD0Candidate<U>()) {
    return table.sliceBy(perD0Candidate, candidate.globalIndex());
  } else if constexpr (isDplusCandidate<U>()) {
    return table.sliceBy(perDplusCandidate, candidate.globalIndex());
  } else if constexpr (isLcCandidate<U>()) {
    return table.sliceBy(perLcCandidate, candidate.globalIndex());
  } else if constexpr (isBplusCandidate<U>()) {
    return table.sliceBy(perBplusCandidate, candidate.globalIndex());
  } else {
    return table;
  }
}

/**
 * returns a slice of the table depending on the type of the HF candidate and index of the collision
 *
 * @param candidate HF candidate that is being checked
 * @param table the table to be sliced
 */
template <typename T, typename U, typename V, typename M, typename N, typename O, typename P>
auto slicedPerHFCollision(T const& table, U const& /*candidates*/, V const& collision, M const& D0CollisionPerCollision, N const& DplusCollisionPerCollision, O const& LcCollisionPerCollision, P const& BplusCollisionPerCollision)
{
  if constexpr (isD0Table<U>() || isD0McTable<U>()) {
    return table.sliceBy(D0CollisionPerCollision, collision.globalIndex());
  } else if constexpr (isDplusTable<U>() || isDplusMcTable<U>()) {
    return table.sliceBy(DplusCollisionPerCollision, collision.globalIndex());
  } else if constexpr (isLcTable<U>() || isLcMcTable<U>()) {
    return table.sliceBy(LcCollisionPerCollision, collision.globalIndex());
  } else if constexpr (isBplusTable<U>() || isBplusMcTable<U>()) {
    return table.sliceBy(BplusCollisionPerCollision, collision.globalIndex());
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
  } else if constexpr (isLcCandidate<T>() || isLcMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kLambdaCPlus);
  } else if constexpr (isBplusCandidate<T>() || isBplusMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kBPlus);
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
  } else if constexpr (isLcTable<T>() || isLcMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kLambdaCPlus);
  } else if constexpr (isBplusTable<T>() || isBplusMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kBPlus);
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
  } else if constexpr (isLcCandidate<T>() || isLcMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassLambdaCPlus);
  } else if constexpr (isBplusCandidate<T>() || isBplusMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassBPlus);
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
  } else if constexpr (isLcTable<T>() || isLcMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassLambdaCPlus);
  } else if constexpr (isBplusTable<T>() || isBplusMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassBPlus);
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
    DplusMCDTable(candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.isCandidateSwapped());
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
  if constexpr (isLcCandidate<T>()) {
    fillLcCandidateTable<isMc>(candidate, HFParTable, HFParETable, HFMlTable, HFMCDTable);
  }
  if constexpr (isBplusCandidate<T>()) {
    fillBplusCandidateTable<isMc>(candidate, HFParTable, HFParETable, HFParDaughterTable, HFMlTable, HFMlDaughterTable, HFMCDTable);
  }
}

template <typename T, typename U>
void fillHFCandidateMcTable(T const& candidate, int32_t mcCollisionIndex, U& BaseMcTable)
{
  BaseMcTable(mcCollisionIndex, candidate.pt(), candidate.eta(), candidate.phi(), candidate.y(), candidate.flagMcMatchGen(), candidate.originMcGen());
}

}; // namespace jethfutilities

#endif // PWGJE_CORE_JETHFUTILITIES_H_
