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

/// \file JetDQUtilities.h
/// \brief Jet DQ related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETDQUTILITIES_H_
#define PWGJE_CORE_JETDQUTILITIES_H_

#include <array>
#include <vector>
#include <string>
#include <optional>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Framework/Logger.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

namespace jetdqutilities
{

/**
 * returns true if the candidate is from a Dielectron table
 */
template <typename T>
constexpr bool isDielectronCandidate()
{
  return std::is_same_v<std::decay_t<T>, CandidatesDielectronData::iterator> || std::is_same_v<std::decay_t<T>, CandidatesDielectronData::filtered_iterator> || std::is_same_v<std::decay_t<T>, CandidatesDielectronMCD::iterator> || std::is_same_v<std::decay_t<T>, CandidatesDielectronMCD::filtered_iterator>;
}

/**
 * returns true if the candidate is from a MC Dielectron table
 */
template <typename T>
constexpr bool isDielectronMcCandidate()
{
  return std::is_same_v<std::decay_t<T>, CandidatesDielectronMCP::iterator> || std::is_same_v<std::decay_t<T>, CandidatesDielectronMCP::filtered_iterator>;
}

/**
 * returns true if the table is a Dielectron table
 */
template <typename T>
constexpr bool isDielectronTable()
{
  return isDielectronCandidate<typename T::iterator>() || isDielectronCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the table is a Dielectron MC table
 */
template <typename T>
constexpr bool isDielectronMcTable()
{
  return isDielectronMcCandidate<typename T::iterator>() || isDielectronMcCandidate<typename T::filtered_iterator>();
}

/**
 * returns true if the candidate is matched to a reconstructed level candidate with the correct decay
 * * @param candidate candidate that is being checked
 */
template <typename T>
constexpr bool isMatchedDielectronCandidate(T const& /*candidate*/)
{
  if constexpr (isDielectronCandidate<T>()) {
    // For now the decision to select signals is done in the DQ framework
    return true;
  } else if constexpr (isDielectronMcCandidate<T>()) {
    return true; // this is true because we always only select Jpsi in our decay channel. If more channels are added this needs to be expanded
  } else {
    return false;
  }
}

/**
 * returns true if the track is a daughter of the dielectron candidate
 *
 * @param track track that is being checked
 * @param candidate Dielectron candidate that is being checked
 * @param tracks the track table
 */
template <typename T, typename U, typename V>
bool isDielectronDaughterTrack(T& track, U& candidate, V const& /*tracks*/)
{
  if constexpr (isDielectronCandidate<U>()) {
    if (candidate.prong0Id() == track.globalIndex() || candidate.prong1Id() == track.globalIndex()) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

/**
 * returns the index of the JMcParticle matched to the Dielectron candidate
 *
 * @param candidate dielectron candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedDielectronParticleId(const T& candidate, const U& /*tracks*/, const V& /*particles*/)
{
  const auto candidateDaughterParticle = candidate.template prong1_as<U>().template mcParticle_as<V>();
  return candidateDaughterParticle.template mothers_first_as<V>().globalIndex(); // can we get the Id directly?
}

/**
 * returns the JMcParticle matched to the Dielectron candidate
 *
 * @param candidate dielectron candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedDielectronParticle(const T& candidate, const U& /*tracks*/, const V& /*particles*/)
{
  const auto candidateDaughterParticle = candidate.template prong1_as<U>().template mcParticle_as<V>();
  return candidateDaughterParticle.template mothers_first_as<V>();
}

/**
 * returns a slice of the table depending on the index of the Dielectron candidate
 *
 * @param candidate dielectron candidate that is being checked
 * @param table the table to be sliced
 */
template <typename T, typename U, typename V>
auto slicedPerDielectronCandidate(T const& table, U const& candidate, V const& perDielectronCandidate)
{
  if constexpr (isDielectronCandidate<U>()) {
    return table.sliceBy(perDielectronCandidate, candidate.globalIndex());
  } else {
    return table;
  }
}

/**
 * returns a slice of the table depending on the type of the Dielectron candidate and index of the collision
 *
 * @param candidate dielectron candidate that is being checked
 * @param table the table to be sliced
 */
template <typename T, typename U, typename V, typename M>
auto slicedPerDielectronCollision(T const& table, U const& /*candidates*/, V const& collision, M const& DielectronCollisionPerCollision)
{
  if constexpr (isDielectronTable<U>() || isDielectronMcTable<U>()) {
    return table.sliceBy(DielectronCollisionPerCollision, collision.globalIndex());
  } else {
    return table;
  }
}

/**
 * returns the Dielectron collision Id of candidate based on type of Dielectron candidate
 *
 * @param candidate dielectron candidate that is being checked
 */
template <typename T>
int getDielectronCandidateCollisionId(T const& candidate)
{
  if constexpr (isDielectronCandidate<T>()) {
    return candidate.reducedeventId();
  } else {
    return -1;
  }
}

/**
 * returns the Dielectron Mc collision Id of candidate based on type of Dielectron candidate
 *
 * @param candidate dielectron candidate that is being checked
 */
template <typename T>
int getDielectronMcCandidateCollisionId(T const& candidate)
{
  if constexpr (isDielectronMcCandidate<T>()) {
    return candidate.dielectronmccollisionId();
  } else {
    return -1;
  }
}

/**
 * returns the PDG of the candidate based on Dielectron Table
 *
 * @param candidate dielectron candidate that is being checked
 */
template <typename T>
int getDielectronCandidatePDG(T const& /*candidate*/)
{
  if constexpr (isDielectronCandidate<T>() || isDielectronMcCandidate<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kJPsi);
  } else {
    return 0;
  }
}

/**
 * returns the PDG of the candidates in the table type
 */
template <typename T>
int getDielectronTablePDG()
{
  if constexpr (isDielectronTable<T>() || isDielectronMcTable<T>()) {
    return static_cast<int>(o2::constants::physics::Pdg::kJPsi);
  } else {
    return 0;
  }
}

/**
 * returns the mass of the candidate based on Dielectron Table
 *
 * @param candidate dielectron candidate that is being checked
 */
template <typename T>
float getDielectronCandidatePDGMass(T const& /*candidate*/)
{
  if constexpr (isDielectronCandidate<T>() || isDielectronMcCandidate<T>()) {
    return static_cast<float>(o2::constants::physics::MassJPsi);
  } else {
    return -1.0;
  }
}

/**
 * returns the mass of the candidates in the table type
 *
 */
template <typename T>
float getDielectronTablePDGMass()
{
  if constexpr (isDielectronTable<T>() || isDielectronMcTable<T>()) {
    return static_cast<float>(o2::constants::physics::MassJPsi);
  } else {
    return -1.0;
  }
}

/**
 * returns the mass of the candidate based on Dielectron Table
 *
 * @param candidate dielectron candidate that is being checked
 */
template <typename T>
float getDielectronCandidateInvariantMass(T const& candidate)
{
  if constexpr (isDielectronCandidate<T>()) {
    return candidate.mass();
  } else {
    return -1.0;
  }
}

template <typename T, typename U>
bool isDielectronParticle(T const& particles, U const& particle)
{
  return RecoDecay::isMatchedMCGen(particles, particle, o2::constants::physics::Pdg::kJPsi, std::array{+kElectron, -kElectron}, true);
}

enum JDielectronParticleDecays {
  JPsiToEE = 0,
};

template <typename T>
bool selectDielectronParticleDecay(T const& dielectronParticle, int dielectronParticleDecaySelection = -1)
{
  if (dielectronParticleDecaySelection == -1) {
    return true;
  }
  return (dielectronParticle.decayFlag() & (1 << dielectronParticleDecaySelection));
}

int initialiseDielectronParticleDecaySelection(std::string dielectronParticleDecaySelection)
{
  if (dielectronParticleDecaySelection == "JPsiToEE") {
    return JDielectronParticleDecays::JPsiToEE;
  }
  return -1;
}

template <typename T, typename U>
uint8_t setDielectronParticleDecayBit(T const& particles, U const& particle)
{
  uint8_t bit = 0;
  if (RecoDecay::isMatchedMCGen(particles, particle, o2::constants::physics::Pdg::kJPsi, std::array{+kElectron, -kElectron}, true)) {
    SETBIT(bit, JDielectronParticleDecays::JPsiToEE);
  }
  return bit;
}

template <typename T, typename U>
void fillDielectronCollisionTable(T const& collision, U& DielectronCollisionTable, int32_t& DielectronCollisionTableIndex)
{
  DielectronCollisionTable(collision.tag_raw(), collision.runNumber(), collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes());
  DielectronCollisionTableIndex = DielectronCollisionTable.lastIndex();
}

template <typename T, typename U>
void fillDielectronMcCollisionTable(T const& mcCollision, U& DielectronMcCollisionTable, int32_t& DielectronMcCollisionTableIndex)
{
  DielectronMcCollisionTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
  DielectronMcCollisionTableIndex = DielectronMcCollisionTable.lastIndex();
}

template <typename T, typename U>
void fillDielectronCandidateTable(T const& candidate, int32_t collisionIndex, U& DielectronTable, int32_t& DielectronCandidateTableIndex)
{
  DielectronTable(collisionIndex, candidate.mass(), candidate.pt(), candidate.eta(), candidate.phi(), candidate.sign(), candidate.filterMap_raw(), candidate.mcDecision());
  DielectronCandidateTableIndex = DielectronTable.lastIndex();
}

template <typename T, typename U>
void fillDielectronCandidateMcTable(T const& candidate, int32_t mcCollisionIndex, U& DielectronMcTable, int32_t& DielectronCandidateTableIndex)
{
  DielectronMcTable(mcCollisionIndex, candidate.pt(), candidate.eta(), candidate.phi(), candidate.y(), candidate.e(), candidate.m(), candidate.pdgCode(), candidate.getGenStatusCode(), candidate.getHepMCStatusCode(), candidate.isPhysicalPrimary(), candidate.decayFlag(), candidate.origin());
  DielectronCandidateTableIndex = DielectronMcTable.lastIndex();
}

}; // namespace jetdqutilities

#endif // PWGJE_CORE_JETDQUTILITIES_H_
