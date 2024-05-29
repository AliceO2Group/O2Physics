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

#include <array>
#include <vector>
#include <string>
#include <optional>

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

namespace jetv0utilities
{

/**
 * returns true if the candidate is from a V0 table
 */
template <typename T>
constexpr bool isV0Candidate()
{
  return std::is_same_v<std::decay_t<T>, CandidatesV0Data::iterator> || std::is_same_v<std::decay_t<T>, CandidatesV0Data::filtered_iterator> || std::is_same_v<std::decay_t<T>, CandidatesV0MCD::iterator> || std::is_same_v<std::decay_t<T>, CandidatesV0MCD::filtered_iterator>;
}

/**
 * returns true if the candidate is from a MC V0 table
 */
template <typename T>
constexpr bool isV0McCandidate()
{
  return std::is_same_v<std::decay_t<T>, CandidatesV0MCP::iterator> || std::is_same_v<std::decay_t<T>, CandidatesV0MCP::filtered_iterator>;
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
  return std::is_same_v<std::decay_t<T>, CandidatesV0MCP> || std::is_same_v<std::decay_t<T>, o2::soa::Filtered<CandidatesV0MCP>>; // note not optimal way but needed for jetfindingutilities::analyseParticles()
}

template <typename T>
bool isV0Particle(T const& particle)
{

  if (TMath::Abs(particle.pdgCode()) == 310) { // k0s
    return true;
  }
  if (TMath::Abs(particle.pdgCode()) == 3122) { // Lambda
    return true;
  }
  return false;
}

enum JV0ParticleDecays {
  K0sToPiPi = 0,
  LambdaPPi = 1
};

template <typename T>
bool selectV0ParticleDecay(T const& v0Particle, int v0ParticleDecaySelection = -1)
{
  if (v0ParticleDecaySelection == -1) {
    return true;
  }
  return (v0Particle.decayFlag() & (1 << v0ParticleDecaySelection));
}

int initialiseV0ParticleDecaySelection(std::string v0ParticleDecaySelection)
{
  if (v0ParticleDecaySelection == "K0sToPiPi") {
    return JV0ParticleDecays::K0sToPiPi;
  } else if (v0ParticleDecaySelection == "LambdaPPi") {
    return JV0ParticleDecays::LambdaPPi;
  }
  return -1;
}

template <typename U, typename T>
uint8_t setV0ParticleDecayBit(T const& particle)
{

  uint8_t bit = 0;

  if (particle.has_daughters()) {
    int daughter1PdgCode = 0;
    int daughter2PdgCode = 0;
    int i = 0;
    for (auto daughter : particle.template daughters_as<U>()) {
      if (i == 0) {
        daughter1PdgCode = daughter.pdgCode();
      }
      if (i == 1) {
        daughter2PdgCode = daughter.pdgCode();
      }
      i++;
    }
    if (TMath::Abs(daughter1PdgCode) == 211 && TMath::Abs(daughter2PdgCode) == 211) {
      SETBIT(bit, JV0ParticleDecays::K0sToPiPi);
    }
    if ((TMath::Abs(daughter1PdgCode) == 211 && TMath::Abs(daughter2PdgCode) == 2212) || (TMath::Abs(daughter2PdgCode) == 211 && TMath::Abs(daughter1PdgCode) == 2212)) {
      SETBIT(bit, JV0ParticleDecays::LambdaPPi);
    }
  }
  return bit;
}

}; // namespace jetv0utilities

#endif // PWGJE_CORE_JETV0UTILITIES_H_
