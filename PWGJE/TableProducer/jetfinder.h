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

// jet finder task header file
//
// Authors: Nima Zardoshti, Jochen Klein

#ifndef PWGJE_TABLEPRODUCER_JETFINDER_H_
#define PWGJE_TABLEPRODUCER_JETFINDER_H_

#include <array>
#include <vector>
#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"

using JetTracks = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection>>;
using JetClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

using ParticlesD0 = o2::soa::Filtered<o2::soa::Join<o2::aod::McParticles, o2::aod::HfCand2ProngMcGen>>;
using ParticlesLc = o2::soa::Filtered<o2::soa::Join<o2::aod::McParticles, o2::aod::HfCand3ProngMcGen>>;
using ParticlesBplus = o2::soa::Filtered<o2::soa::Join<o2::aod::McParticles, o2::aod::HfCandBplusMcGen>>;

using CandidatesD0Data = o2::soa::Filtered<o2::soa::Join<o2::aod::HfCand2Prong, o2::aod::HfSelD0>>;
using CandidatesD0MCD = o2::soa::Filtered<o2::soa::Join<o2::aod::HfCand2Prong, o2::aod::HfSelD0, o2::aod::HfCand2ProngMcRec>>;

using CandidatesBplusData = o2::soa::Filtered<o2::soa::Join<o2::aod::HfCandBplus, o2::aod::HfSelBplusToD0Pi>>;
using CandidatesBplusMCD = o2::soa::Filtered<o2::soa::Join<o2::aod::HfCandBplus, o2::aod::HfSelBplusToD0Pi, o2::aod::HfCandBplusMcRec>>;

using CandidatesLcData = o2::soa::Filtered<o2::soa::Join<o2::aod::HfCand3Prong, o2::aod::HfSelLc>>;
using CandidatesLcMCD = o2::soa::Filtered<o2::soa::Join<o2::aod::HfCand3Prong, o2::aod::HfSelLc, o2::aod::HfCand3ProngMcRec>>;

// functions for track, cluster and candidate selection

// function that performs track selections on each track
template <typename T>
bool selectTrack(T const& track, std::string trackSelection)
{
  if (trackSelection == "globalTracks") {
    return track.isGlobalTrackWoPtEta();
  } else if (trackSelection == "QualityTracks") {
    return track.isQualityTrack();
  } else if (trackSelection == "hybridTracksJE") {
    return track.trackCutFlagFb5();
  } else {
    return true;
  }
}

// function that adds tracks to the fastjet list, removing daughters of 2Prong candidates
template <typename T, typename U>
void analyseTracks(std::vector<fastjet::PseudoJet>& inputParticles, T const& tracks, std::string trackSelection, std::optional<U> const& candidate = std::nullopt)
{
  for (auto& track : tracks) {
    if (!selectTrack(track, trackSelection)) {
      continue;
    }
    if (candidate != std::nullopt) {
      auto cand = candidate.value();
      if constexpr (std::is_same_v<std::decay_t<U>, CandidatesD0Data::iterator> || std::is_same_v<std::decay_t<U>, CandidatesD0Data::filtered_iterator> || std::is_same_v<std::decay_t<U>, CandidatesD0MCD::iterator> || std::is_same_v<std::decay_t<U>, CandidatesD0MCD::filtered_iterator>) {
        if (cand.template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
      }

      if constexpr (std::is_same_v<std::decay_t<U>, CandidatesLcData::iterator> || std::is_same_v<std::decay_t<U>, CandidatesLcData::filtered_iterator> || std::is_same_v<std::decay_t<U>, CandidatesLcMCD::iterator> || std::is_same_v<std::decay_t<U>, CandidatesLcMCD::filtered_iterator>) {
        if (cand.template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong1_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong2_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
      }

      if constexpr (std::is_same_v<std::decay_t<U>, CandidatesBplusData::iterator> || std::is_same_v<std::decay_t<U>, CandidatesBplusData::filtered_iterator> || std::is_same_v<std::decay_t<U>, CandidatesBplusMCD::iterator> || std::is_same_v<std::decay_t<U>, CandidatesBplusMCD::filtered_iterator>) {
        if (cand.template prong0_as<o2::aod::HfCand2Prong>().template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong0_as<o2::aod::HfCand2Prong>().template prong1_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
      }
    }
    FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex());
  }
}

// function that adds clusters to the fastjet list
template <typename T>
void analyseClusters(std::vector<fastjet::PseudoJet>& inputParticles, T const& clusters)
{
  for (auto& cluster : *clusters) {
    // add cluster selections
    FastJetUtilities::fillClusters(cluster, inputParticles, cluster.globalIndex());
  }
}

// function that takes any generic candidate, performs selections and adds the candidate to the fastjet list
template <typename T>
bool analyseCandidate(std::vector<fastjet::PseudoJet>& inputParticles, int candPDG, float candPtMin, float candPtMax, float candYMin, float candYMax, T const& candidate)
{
  if (candidate.y(RecoDecay::getMassPDG(candPDG)) < candYMin || candidate.y(RecoDecay::getMassPDG(candPDG)) > candYMax) {
    return false;
  }
  if (candidate.pt() < candPtMin || candidate.pt() >= candPtMax) {
    return false;
  }
  FastJetUtilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
  return true;
}

// function that checks the MC status of a candidate and then calls the function to analyseCandidates
template <typename T>
bool analyseCandidateMC(std::vector<fastjet::PseudoJet>& inputParticles, int candPDG, int candDecay, float candPtMin, float candPtMax, float candYMin, float candYMax, T const& candidate, bool rejectBackgroundMCCandidates)
{
  if (rejectBackgroundMCCandidates && !(std::abs(candidate.flagMcMatchRec()) == 1 << candDecay)) {
    return false;
  }
  return analyseCandidate(inputParticles, candPDG, candPtMin, candPtMax, candYMin, candYMax, candidate);
}

// function that calls the jet finding and fills the relevant tables
template <typename T, typename U, typename V, typename W>
void findJets(JetFinder& jetFinder, std::vector<fastjet::PseudoJet>& inputParticles, std::vector<double> jetRadius, T const& collision, U& jetsTable, V& constituentsTable, W& constituentsSubTable, bool DoConstSub, bool doHFJetFinding = false)
{
  // auto candidatepT = 0.0;
  auto jetRValues = static_cast<std::vector<double>>(jetRadius);
  for (auto R : jetRValues) {
    jetFinder.jetR = R;
    std::vector<fastjet::PseudoJet> jets;
    fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
    for (const auto& jet : jets) {
      bool isHFJet = false;
      if (doHFJetFinding) {
        for (const auto& constituent : jet.constituents()) {
          if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
            isHFJet = true;
            // candidatepT = constituent.pt();
            break;
          }
        }
        if (!isHFJet) {
          continue;
        }
      }
      std::vector<int> trackconst;
      std::vector<int> candconst;
      std::vector<int> clusterconst;
      jetsTable(collision.globalIndex(), jet.pt(), jet.eta(), jet.phi(),
                jet.E(), jet.m(), jet.area(), std::round(R * 100));
      for (const auto& constituent : sorted_by_pt(jet.constituents())) {
        // need to add seperate thing for constituent subtraction
        if (DoConstSub) { // FIXME: needs to be addressed in Haadi's PR
          constituentsSubTable(jetsTable.lastIndex(), constituent.pt(), constituent.eta(), constituent.phi(),
                               constituent.E(), constituent.m(), constituent.user_index());
        }

        if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
          trackconst.push_back(constituent.template user_info<FastJetUtilities::fastjet_user_info>().getIndex());
        }
        if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::cluster)) {
          clusterconst.push_back(constituent.template user_info<FastJetUtilities::fastjet_user_info>().getIndex());
        }
        if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
          candconst.push_back(constituent.template user_info<FastJetUtilities::fastjet_user_info>().getIndex());
        }
      }
      constituentsTable(jetsTable.lastIndex(), trackconst, clusterconst, candconst);
    }
  }
}

// function that checks if a candidate has any daughters that need to be removed from the event at gen level
template <typename T>
bool checkDaughters(T const& particle, int globalIndex)
{
  for (auto daughter : particle.template daughters_as<o2::aod::McParticles>()) {
    if (daughter.globalIndex() == globalIndex) {
      return true;
    }
    if (checkDaughters(daughter, globalIndex)) {
      return true;
    }
  }
  return false;
}

template <typename T, typename U>
void analyseParticles(std::vector<fastjet::PseudoJet>& inputParticles, std::string particleSelection, int jetTypeParticleLevel, T const& particles, TDatabasePDG* pdg, std::optional<U> const& candidate = std::nullopt)
{
  inputParticles.clear();
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
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    auto pdgCharge = pdgParticle ? std::abs(pdgParticle->Charge()) : -1.0;
    if (jetTypeParticleLevel == static_cast<int>(JetType::charged) && pdgCharge < 3.0) {
      continue;
    }
    if (jetTypeParticleLevel == static_cast<int>(JetType::neutral) && pdgCharge != 0.0) {
      continue;
    }
    if (candidate != std::nullopt) {
      auto cand = candidate.value();
      if (checkDaughters(cand, particle.globalIndex())) {
        continue;
      }
    }
    FastJetUtilities::fillTracks(particle, inputParticles, particle.globalIndex(), static_cast<int>(JetConstituentStatus::track), RecoDecay::getMassPDG(particle.pdgCode()));
  }
}

template <typename T>
bool selectCollision(T const& collision, std::string eventSelection)
{
  if (eventSelection == "sel8") {
    return collision.sel8();
  } else if (eventSelection == "sel7") {
    return collision.sel7();
  } else {
    return true;
  }
}

template <typename T>
bool hasEMCAL(T const& collision)
{
  // Check for EMCAL in readout requires any of the EMCAL trigger classes (including EMC in MB trigger) to fire
  std::array<triggerAliases, 11> selectAliases = {{triggerAliases::kTVXinEMC, triggerAliases::kEMC7, triggerAliases::kDMC7, triggerAliases::kEG1, triggerAliases::kEG2, triggerAliases::kDG1, triggerAliases::kDG2, triggerAliases::kEJ1, triggerAliases::kEJ2, triggerAliases::kDJ1, triggerAliases::kDJ2}};
  bool found = false;
  for (auto alias : selectAliases) {
    if (collision.alias_bit(alias)) {
      found = true;
      break;
    }
  }
  return found;
}

#endif // PWGJE_TABLEPRODUCER_JETFINDER_H_
