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
 * @param candMass pdg mass of hf candidate
 * @param candPtMin minimum pT of hf candidate
 * @param candPtMax maximum pT of hf candidate
 * @param candYMin minimum Y of hf candidate
 * @param candYMax maximum Y of hf candidate
 * @param candidate hf candidate
 */
template <typename T>
bool analyseCandidate(std::vector<fastjet::PseudoJet>& inputParticles, T const& candidate, float candPtMin, float candPtMax, float candYMin, float candYMax, float candMass)
{
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
 * @param candMass pdg mass of hf candidate
 * @param candPtMin minimum pT of hf candidate
 * @param candPtMax maximum pT of hf candidate
 * @param candYMin minimum Y of hf candidate
 * @param candYMax maximum Y of hf candidate
 * @param candidate hf candidate
 * @param rejectBackgroundMCCandidates choose whether to accept background hf candidates as defined by the selection flag
 */
template <typename T>
bool analyseCandidateMC(std::vector<fastjet::PseudoJet>& inputParticles, T const& candidate, float candPtMin, float candPtMax, float candYMin, float candYMax, float candMass, bool rejectBackgroundMCCandidates)
{
  if (rejectBackgroundMCCandidates && !jethfutilities::isMatchedHFCandidate(candidate)) {
    return false;
  }
  return analyseCandidate(inputParticles, candidate, candPtMin, candPtMax, candYMin, candYMax, candMass);
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
void findJets(JetFinder& jetFinder, std::vector<fastjet::PseudoJet>& inputParticles, std::vector<double> jetRadius, float jetAreaFractionMin, T const& collision, U& jetsTable, V& constituentsTable, bool doHFJetFinding = false)
{
  auto jetRValues = static_cast<std::vector<double>>(jetRadius);
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
      jetsTable(collision.globalIndex(), jet.pt(), jet.eta(), jet.phi(),
                jet.E(), jet.m(), jet.has_area() ? jet.area() : 0., std::round(R * 100));
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
      }
      constituentsTable(jetsTable.lastIndex(), trackconst, clusterconst, candconst);
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
        auto hfParticle = cand.template mcParticle_as<T>();
        if (jethfutilities::isDaughterParticle(hfParticle, particle.globalIndex()) || (hfParticle.globalIndex() == particle.globalIndex())) {
          continue;
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
