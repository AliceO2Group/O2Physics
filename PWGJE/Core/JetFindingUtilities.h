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

/// \file JetFindingUtilities.h
/// \brief Jet Finding related utilities
///
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#ifndef PWGJE_CORE_JETFINDINGUTILITIES_H_
#define PWGJE_CORE_JETFINDINGUTILITIES_H_

#include <array>
#include <vector>
#include <string>
#include <optional>
#include <cmath>

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

// #include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/DataModel/Jet.h"

namespace jetfindingutilities
{

/**
 * Adds tracks to a fastjet inputParticles list
 *
 * @param inputParticles fastjet container
 * @param tracks track table to be added
 * @param trackSelection track selection to be applied to tracks
 * @param candidate optional HF candidiate
 */

template <typename T, typename U>
void analyseTracks(std::vector<fastjet::PseudoJet>& inputParticles, T const& tracks, int trackSelection, std::optional<U> const& candidate = std::nullopt)
{
  for (auto& track : tracks) {
    if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
      continue;
    }
    if (candidate != std::nullopt) {
      auto cand = candidate.value();
      if (jethfutilities::isDaughterTrack(track, cand, tracks)) {
        continue;
      }
    }
    fastjetutilities::fillTracks(track, inputParticles, track.globalIndex());
  }
}

/**
 * Adds tracks to a fastjet inputParticles list for the case where there are multiple candidates per event
 *
 * @param inputParticles fastjet container
 * @param tracks track table to be added
 * @param trackSelection track selection to be applied to tracks
 * @param candidates candidiates
 */

template <typename T, typename U>
void analyseTracksMultipleCandidates(std::vector<fastjet::PseudoJet>& inputParticles, T const& tracks, int trackSelection, U const& candidates)
{
  for (auto& track : tracks) {
    if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
      continue;
    }
    for (auto& candidate : candidates) {
      if (jethfutilities::isDaughterTrack(track, candidate, tracks)) {
        continue;
      }
    }
    fastjetutilities::fillTracks(track, inputParticles, track.globalIndex());
  }
}

/**
 * Adds clusters to a fastjet inputParticles list
 *
 * @param inputParticles fastjet container
 * @param clusters track table to be added
 */
template <typename T>
void analyseClusters(std::vector<fastjet::PseudoJet>& inputParticles, T const& clusters)
{
  for (auto& cluster : *clusters) {
    // add cluster selections
    fastjetutilities::fillClusters(cluster, inputParticles, cluster.globalIndex());
  }
}

/**
 * Adds hf candidates to a fastjet inputParticles list (for data)
 *
 * @param inputParticles fastjet container
 * @param candPtMin minimum pT of hf candidate
 * @param candPtMax maximum pT of hf candidate
 * @param candYMin minimum Y of hf candidate
 * @param candYMax maximum Y of hf candidate
 * @param candidate hf candidate
 */
template <typename T>
bool analyseCandidate(std::vector<fastjet::PseudoJet>& inputParticles, T const& candidate, float candPtMin, float candPtMax, float candYMin, float candYMax)
{
  auto candMass = jethfutilities::getCandidatePDGMass(candidate);
  if (isnan(candidate.y())) {
    return false;
  }
  if (candidate.y() < candYMin || candidate.y() > candYMax) {
    return false;
  }
  if (candidate.pt() < candPtMin || candidate.pt() >= candPtMax) {
    return false;
  }
  fastjetutilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), candMass);
  return true;
}

/**
 * Adds hf candidates to a fastjet inputParticles list (for MC det)
 *
 * @param inputParticles fastjet container
 * @param candPtMin minimum pT of hf candidate
 * @param candPtMax maximum pT of hf candidate
 * @param candYMin minimum Y of hf candidate
 * @param candYMax maximum Y of hf candidate
 * @param candidate hf candidate
 * @param rejectBackgroundMCCandidates choose whether to accept background hf candidates as defined by the selection flag
 */
template <typename T>
bool analyseCandidateMC(std::vector<fastjet::PseudoJet>& inputParticles, T const& candidate, float candPtMin, float candPtMax, float candYMin, float candYMax, bool rejectBackgroundMCCandidates)
{
  if (rejectBackgroundMCCandidates && !jethfutilities::isMatchedHFCandidate(candidate)) {
    return false;
  }
  return analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax);
}

/**
 * Adds hf candidates to a fastjet inputParticles list (for data)
 *
 * @param inputParticles fastjet container
 * @param v0Mass pdg mass of v0 candidate
 * @param v0PtMin minimum pT of v0 candidate
 * @param v0PtMax maximum pT of v candidate
 * @param v0EtaMin minimum eta of v0 candidate
 * @param v0EtaMax maximum eta of v0 candidate
 * @param v0s V0 candidates
 */
template <typename T>
bool analyseV0s(std::vector<fastjet::PseudoJet>& inputParticles, T const& v0s, float v0PtMin, float v0PtMax, float v0YMin, float v0YMax, int v0Index)
{
  float v0Mass = 0;
  float v0Y = -10.0;

  int nSelectedV0s = 0;
  for (auto const& v0 : v0s) {
    // if (isnan(candidate.y())) {
    //   continue;
    // }
    if constexpr (jetv0utilities::isV0McTable<T>()) {
      v0Mass = v0.m();
      v0Y = v0.y();
    } else {
      if (v0Index == 0) {
        v0Mass = o2::constants::physics::MassKaonNeutral;
      }
      if (v0Index == 1) {
        v0Mass = o2::constants::physics::MassLambda0;
      }
      v0Y = v0.rapidity(v0Index);
    }
    if (v0Y < v0YMin || v0Y > v0YMax) {
      continue;
    }
    if (v0.pt() < v0PtMin || v0.pt() >= v0PtMax) {
      continue;
    }
    fastjetutilities::fillTracks(v0, inputParticles, v0.globalIndex(), static_cast<int>(JetConstituentStatus::v0), v0Mass);
    nSelectedV0s++;
  }
  if (nSelectedV0s > 0) {
    return true;
  } else {
    return false;
  }
}

/**
 * Performs jet finding and fills jet tables
 *
 * @param jetFinder JetFinder object which carries jet finding parameters
 * @param inputParticles fastjet container
 * @param jetRadius jet finding radii
 * @param collision the collision within which jets are being found
 * @param jetsTable output table of jets
 * @param constituentsTable output table of jet constituents
 * @param doHFJetFinding set whether only jets containing a HF candidate are saved
 */
template <typename T, typename U, typename V>
void findJets(JetFinder& jetFinder, std::vector<fastjet::PseudoJet>& inputParticles, float jetPtMin, float jetPtMax, std::vector<double> jetRadius, float jetAreaFractionMin, T const& collision, U& jetsTable, V& constituentsTable, bool doHFJetFinding = false)
{
  auto jetRValues = static_cast<std::vector<double>>(jetRadius);
  jetFinder.jetPtMin = jetPtMin;
  jetFinder.jetPtMax = jetPtMax;
  for (auto R : jetRValues) {
    jetFinder.jetR = R;
    std::vector<fastjet::PseudoJet> jets;
    fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
    for (const auto& jet : jets) {
      bool isHFJet = false;
      if (doHFJetFinding) {
        for (const auto& constituent : jet.constituents()) {
          if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
            isHFJet = true;
            break;
          }
        }
        if (!isHFJet) {
          continue;
        }
      }
      if (jet.has_area() && jet.area() < jetAreaFractionMin * M_PI * R * R) {
        continue;
      }
      std::vector<int> trackconst;
      std::vector<int> candconst;
      std::vector<int> clusterconst;
      std::vector<int> v0const;
      jetsTable(collision.globalIndex(), jet.pt(), jet.eta(), jet.phi(),
                jet.E(), jet.rapidity(), jet.m(), jet.has_area() ? jet.area() : 0., std::round(R * 100));
      for (const auto& constituent : sorted_by_pt(jet.constituents())) {
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
          trackconst.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::cluster)) {
          clusterconst.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
          candconst.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
        if (constituent.template user_info<fastjetutilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::v0)) {
          v0const.push_back(constituent.template user_info<fastjetutilities::fastjet_user_info>().getIndex());
        }
      }
      constituentsTable(jetsTable.lastIndex(), trackconst, clusterconst, candconst, v0const);
    }
  }
}

/**
 * Adds particles to a fastjet inputParticles list
 *
 * @param inputParticles fastjet container
 * @param particleSelection particle selection to be applied to particles
 * @param jetTypeParticleLevel set whether charged particles, neutral particles or both are accepted
 * @param particles particle table to be added
 * @param pdgDatabase database of pdg codes
 * @param candidate optional hf candidiate
 */
template <typename T, typename U>
void analyseParticles(std::vector<fastjet::PseudoJet>& inputParticles, std::string particleSelection, int jetTypeParticleLevel, T const& particles, o2::framework::Service<o2::framework::O2DatabasePDG> pdgDatabase, std::optional<U> const& candidate = std::nullopt)
{
  for (auto& particle : particles) {
    if (particleSelection == "PhysicalPrimary" && !particle.isPhysicalPrimary()) { // CHECK : Does this exclude the HF hadron?
      continue;
    } else if (particleSelection == "HepMCStatus" && particle.getHepMCStatusCode() != 1) { // do we need isPhysicalPrimary as well? Note: Might give unforseen results if the generator isnt PYTHIA
      continue;
    } else if (particleSelection == "GenStatus" && particle.getGenStatusCode() != 1) {
      continue;
    } else if (particleSelection == "PhysicalPrimaryAndHepMCStatus" && (!particle.isPhysicalPrimary() || particle.getHepMCStatusCode() != 1)) {
      continue;
    }
    if (isinf(particle.eta())) {
      continue;
    }
    auto pdgParticle = pdgDatabase->GetParticle(particle.pdgCode());
    auto pdgCharge = pdgParticle ? std::abs(pdgParticle->Charge()) : -1.0;
    if (jetTypeParticleLevel == static_cast<int>(JetType::charged) && pdgCharge < 3.0) {
      continue;
    }
    if (jetTypeParticleLevel == static_cast<int>(JetType::neutral) && pdgCharge != 0.0) {
      continue;
    }
    if constexpr (jethfutilities::isHFMcCandidate<U>()) {
      if (candidate != std::nullopt) {
        auto cand = candidate.value();
        if (cand.mcParticleId() == particle.globalIndex()) {
          continue;
        }
        auto hfParticle = cand.template mcParticle_as<T>();
        if (jethfutilities::isDaughterParticle(hfParticle, particle.globalIndex())) {
          continue;
        }
      }
    }
    if constexpr (jetv0utilities::isV0McTable<U>()) { // note that for V0s the candidate table is given to this function, not a single candidate
      if (candidate != std::nullopt) {
        auto cands = candidate.value();
        for (auto const& cand : cands) {
          if (cand.mcParticleId() == particle.globalIndex()) {
            continue;
          }
          auto V0Particle = cand.template mcParticle_as<T>();
          if (jethfutilities::isDaughterParticle(V0Particle, particle.globalIndex())) {
            continue;
          }
        }
      }
    }
    fastjetutilities::fillTracks(particle, inputParticles, particle.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdgParticle->Mass());
  }
}

template <typename T>
bool isInEtaAcceptance(T const& jet, float jetEtaMin, float jetEtaMax, float etaMin = -0.9, float etaMax = 0.9)
{
  if (jetEtaMin < -98.0) {
    return (jet.eta() >= etaMin + (jet.r() / 100.0) && jet.eta() <= etaMax - (jet.r() / 100.0));
  } else {
    return (jet.eta() >= jetEtaMin && jet.eta() <= jetEtaMax);
  }
}

}; // namespace jetfindingutilities

#endif // PWGJE_CORE_JETFINDINGUTILITIES_H_
