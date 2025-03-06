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

/// \file JetCandidateUtilities.h
/// \brief Jet candidate related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETCANDIDATEUTILITIES_H_
#define PWGJE_CORE_JETCANDIDATEUTILITIES_H_

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
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetDQUtilities.h"
#include "PWGJE/Core/JetV0Utilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

namespace jetcandidateutilities
{

/**
 * returns true if the candidate is from a candidate table
 * * @param candidate candidate that is being checked
 */
template <typename T>
constexpr bool isCandidate()
{
  if constexpr (jethfutilities::isHFCandidate<T>()) {
    return true;
  } else if constexpr (jetv0utilities::isV0Candidate<T>()) {
    return true;
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>()) {
    return true;
  } else {
    return false;
  }
}

/**
 * returns true if the candidate is from a MC candidate table
 * * @param candidate candidate that is being checked
 */
template <typename T>
constexpr bool isMcCandidate()
{
  if constexpr (jethfutilities::isHFMcCandidate<T>()) {
    return true;
  } else if constexpr (jetv0utilities::isV0McCandidate<T>()) {
    return true;
  } else if constexpr (jetdqutilities::isDielectronMcCandidate<T>()) {
    return true;
  } else {
    return false;
  }
}

/**
 * returns true if the table type is a candidate table
 */
template <typename T>
constexpr bool isCandidateTable()
{
  if constexpr (jethfutilities::isHFTable<T>()) {
    return true;
  } else if constexpr (jetv0utilities::isV0Table<T>()) {
    return true;
  } else if constexpr (jetdqutilities::isDielectronTable<T>()) {
    return true;
  } else {
    return false;
  }
}

/**
 * returns true if the table type is a candidate table
 */
template <typename T>
constexpr bool isCandidateMcTable()
{
  if constexpr (jethfutilities::isHFMcTable<T>()) {
    return true;
  } else if constexpr (jetv0utilities::isV0McTable<T>()) {
    return true;
  } else if constexpr (jetdqutilities::isDielectronMcTable<T>()) {
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
constexpr bool isMatchedCandidate(T const& candidate)
{
  if constexpr (jethfutilities::isHFCandidate<T>() || jethfutilities::isHFMcCandidate<T>()) {
    return jethfutilities::isMatchedHFCandidate(candidate);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>() || jetdqutilities::isDielectronMcCandidate<T>()) {
    return jetdqutilities::isMatchedDielectronCandidate(candidate);
  } else {
    return false;
  }
}

/**
 * returns true if the track is a daughter of the candidate
 *
 * @param track track that is being checked
 * @param candidate candidate that is being checked
 * @param tracks the track table
 */
template <typename T, typename U, typename V>
bool isDaughterTrack(T& track, U& candidate, V const& tracks)
{
  if constexpr (jethfutilities::isHFCandidate<U>()) {
    return jethfutilities::isHFDaughterTrack(track, candidate, tracks);
  } else if constexpr (jetv0utilities::isV0Candidate<U>()) {
    return jetv0utilities::isV0DaughterTrack(track, candidate, tracks);
  } else if constexpr (jetdqutilities::isDielectronCandidate<U>()) {
    return jetdqutilities::isDielectronDaughterTrack(track, candidate, tracks);
  } else {
    return false;
  }
}

/**
 * returns true if the particle has any daughters with the given global index
 *
 * @param candidate mother hf particle that is being checked
 * @param globalIndex global index of potnetial daughter particle
 */
template <typename T>
bool isDaughterParticle(const T& particle, int globalIndex)
{
  for (auto daughter : particle.template daughters_as<typename std::decay_t<T>::parent_t>()) {
    if (daughter.globalIndex() == globalIndex) {
      return true;
    }
    if (isDaughterParticle(daughter, globalIndex)) {
      return true;
    }
  }
  return false;
}

/**
 * returns the index of the JMcParticle matched to the candidate
 *
 * @param candidate hf candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedParticleId(const T& candidate, const U& tracks, const V& particles)
{
  if constexpr (jethfutilities::isHFCandidate<T>()) {
    return jethfutilities::matchedHFParticleId(candidate, tracks, particles);
  } else if constexpr (jetv0utilities::isV0Candidate<T>()) {
    return jetv0utilities::matchedV0ParticleId(candidate, tracks, particles);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>()) {
    return jetdqutilities::matchedDielectronParticleId(candidate, tracks, particles);
  } else {
    return -1;
  }
}

/**
 * returns the JMcParticle matched to the candidate
 *
 * @param candidate hf candidate that is being checked
 * @param tracks track table
 * @param particles particle table
 */
template <typename T, typename U, typename V>
auto matchedParticle(const T& candidate, const U& tracks, const V& particles)
{
  if constexpr (jethfutilities::isHFCandidate<T>()) {
    return jethfutilities::matchedHFParticle(candidate, tracks, particles);
  } else if constexpr (jetv0utilities::isV0Candidate<T>()) {
    return jetv0utilities::matchedV0Particle(candidate, tracks, particles);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>()) {
    return jetdqutilities::matchedDielectronParticle(candidate, tracks, particles);
  } else {
    return jethfutilities::matchedHFParticle(candidate, tracks, particles); // this is a dummy output which should never be triggered
  }
}

/**
 * returns a slice of the table depending on the index of the candidate
 *
 * @param candidate candidate that is being checked
 * @param table the table to be sliced
 */
template <typename T, typename U, typename V, typename M, typename N, typename O, typename P>
auto slicedPerCandidate(T const& table, U const& candidate, V const& perD0Candidate, M const& perDplusCandidate, N const& perLcCandidate, O const& perBplusCandidate, P const& perDielectronCandidate)
{
  if constexpr (jethfutilities::isHFCandidate<U>()) {
    return jethfutilities::slicedPerHFCandidate(table, candidate, perD0Candidate, perDplusCandidate, perLcCandidate, perBplusCandidate);
  } else if constexpr (jetdqutilities::isDielectronCandidate<U>()) {
    return jetdqutilities::slicedPerDielectronCandidate(table, candidate, perDielectronCandidate);
  } else {
    return table;
  }
}

/**
 * returns a slice of the table depending on the index of the candidate
 * @param CandidateTable candidtae table type
 * @param jet jet that the slice is based on
 * @param table the table to be sliced
 */
template <typename CandidateTable, typename T, typename U, typename V, typename M, typename N, typename O, typename P>
auto slicedPerJet(T const& table, U const& jet, V const& perD0Jet, M const& perDplusJet, N const& perLcJet, O const& perBplusJet, P const& perDielectronJet)
{
  if constexpr (jethfutilities::isHFTable<CandidateTable>() || jethfutilities::isHFMcTable<CandidateTable>()) {
    return jethfutilities::slicedPerHFJet<CandidateTable>(table, jet, perD0Jet, perDplusJet, perLcJet, perBplusJet);
  } else if constexpr (jetdqutilities::isDielectronTable<CandidateTable>() || jetdqutilities::isDielectronMcTable<CandidateTable>()) {
    return jetdqutilities::slicedPerDielectronJet<CandidateTable>(table, jet, perDielectronJet);
  } else {
    return table;
  }
}

/**
 * returns the candidate collision Id of candidate based on type of candidate
 *
 * @param candidate candidate that is being checked
 */
template <typename T>
int getCandidateCollisionId(T const& candidate)
{
  if constexpr (jethfutilities::isHFCandidate<T>()) {
    return jethfutilities::getHFCandidateCollisionId(candidate);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>()) {
    return jetdqutilities::getDielectronCandidateCollisionId(candidate);
  } else {
    return -1;
  }
}

/**
 * returns the candidate Mc collision Id of candidate based on type of candidate
 *
 * @param candidate candidate that is being checked
 */
template <typename T>
int getMcCandidateCollisionId(T const& candidate)
{
  if constexpr (jethfutilities::isHFMcCandidate<T>()) {
    return jethfutilities::getHFMcCandidateCollisionId(candidate);
  } else if constexpr (jetdqutilities::isDielectronMcCandidate<T>()) {
    return jetdqutilities::getDielectronMcCandidateCollisionId(candidate);
  } else {
    return -1;
  }
}

/**
 * returns the PDG of the candidate based on Table
 *
 * @param candidate candidate that is being checked
 */

template <typename T>
int getCandidatePDG(T const& candidate)
{
  if constexpr (jethfutilities::isHFCandidate<T>() || jethfutilities::isHFMcCandidate<T>()) {
    return jethfutilities::getHFCandidatePDG(candidate);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>() || jetdqutilities::isDielectronMcCandidate<T>()) {
    return jetdqutilities::getDielectronCandidatePDG(candidate);
  } else {
    return 0;
  }
}

/**
 * returns the PDG of the candidates in the table type
 */
template <typename T>
int getTablePDG()
{
  if constexpr (jethfutilities::isHFTable<T>() || jethfutilities::isHFMcTable<T>()) {
    return jethfutilities::getHFTablePDG<T>();
  } else if constexpr (jetdqutilities::isDielectronTable<T>() || jetdqutilities::isDielectronMcTable<T>()) {
    return jetdqutilities::getDielectronTablePDG<T>();
  } else {
    return 0;
  }
}

/**
 * returns the pdg mass of the candidate based on Table
 *
 * @param candidate candidate that is being checked
 */
template <typename T>
float getCandidatePDGMass(T const& candidate)
{
  if constexpr (jethfutilities::isHFCandidate<T>() || jethfutilities::isHFMcCandidate<T>()) {
    return jethfutilities::getHFCandidatePDGMass(candidate);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>() || jetdqutilities::isDielectronMcCandidate<T>()) {
    return jetdqutilities::getDielectronCandidatePDGMass(candidate);
  } else {
    return -1.0;
  }
}

/**
 * returns the pdg mass of the candidates in the table type
 *
 */
template <typename T>
float getTablePDGMass()
{
  if constexpr (jethfutilities::isHFTable<T>() || jethfutilities::isHFMcTable<T>()) {
    return jethfutilities::getHFTablePDGMass<T>();
  } else if constexpr (jetdqutilities::isDielectronTable<T>() || jetdqutilities::isDielectronMcTable<T>()) {
    return jetdqutilities::getDielectronTablePDGMass<T>();
  } else {
    return -1.0;
  }
}

/**
 * returns the invariant mass of the candidate based on Table
 *
 * @param candidate  candidate that is being checked
 */
template <typename T>
float getCandidateInvariantMass(T const& candidate)
{
  if constexpr (jethfutilities::isHFCandidate<T>()) {
    return jethfutilities::getHFCandidateInvariantMass(candidate);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>()) {
    return jetdqutilities::getDielectronCandidateInvariantMass(candidate);
  } else {
    return -1.0;
  }
}

template <typename T, typename U, typename V>
void fillCandidateCollisionTable(T const& collision, U const& /*candidates*/, V& CandiateCollisionTable)
{
  if constexpr (jethfutilities::isHFTable<U>()) {
    jethfutilities::fillHFCollisionTable(collision, CandiateCollisionTable);
  } else if constexpr (jetdqutilities::isDielectronTable<U>()) {
    jetdqutilities::fillDielectronCollisionTable(collision, CandiateCollisionTable); // if more dilepton tables are added we would need a fillDQCollisionTable
  }
}

template <typename T, typename U, typename V>
void fillCandidateMcCollisionTable(T const& mcCollision, U const& /*candidates*/, V& CandiateMcCollisionTable)
{
  if constexpr (jethfutilities::isHFMcTable<U>()) {
    jethfutilities::fillHFMcCollisionTable(mcCollision, CandiateMcCollisionTable);
  } else if constexpr (jetdqutilities::isDielectronMcTable<U>()) {
    jetdqutilities::fillDielectronMcCollisionTable(mcCollision, CandiateMcCollisionTable);
  }
}

template <bool isMc, typename T, typename U, typename V, typename M, typename N, typename O, typename P, typename Q, typename S>
void fillCandidateTable(T const& candidate, int32_t collisionIndex, U& BaseTable, V& HFParTable, M& HFParETable, N& HFParDaughterTable, O& HFSelectionFlagTable, P& HFMlTable, Q& HFMlDaughterTable, S& HFMCDTable)
{
  if constexpr (jethfutilities::isHFCandidate<T>()) {
    jethfutilities::fillHFCandidateTable<isMc>(candidate, collisionIndex, BaseTable, HFParTable, HFParETable, HFParDaughterTable, HFSelectionFlagTable, HFMlTable, HFMlDaughterTable, HFMCDTable);
  } else if constexpr (jetdqutilities::isDielectronCandidate<T>()) {
    jetdqutilities::fillDielectronCandidateTable(candidate, collisionIndex, BaseTable);
  }
}

template <typename T, typename U>
void fillCandidateMcTable(T const& candidate, int32_t mcCollisionIndex, U& BaseMcTable)
{
  if constexpr (jethfutilities::isHFMcCandidate<T>()) {
    jethfutilities::fillHFCandidateMcTable(candidate, mcCollisionIndex, BaseMcTable);
  } else if constexpr (jetdqutilities::isDielectronMcCandidate<T>()) {
    jetdqutilities::fillDielectronCandidateMcTable(candidate, mcCollisionIndex, BaseMcTable);
  }
}

}; // namespace jetcandidateutilities

#endif // PWGJE_CORE_JETCANDIDATEUTILITIES_H_
