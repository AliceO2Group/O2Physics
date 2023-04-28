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

using JetTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
using JetClusters = o2::soa::Filtered<o2::aod::EMCALClusters>;

using JetParticles2Prong = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;
using JetParticles3Prong = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
using JetParticlesBplus = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandBplusMcGen>>;

using CandidateD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
using CandidateD0MC = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>>;

using CandidateBplusData = soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>>;
using CandidateBplusMC = soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>>;

using CandidateLcData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
using CandidateLcMC = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;

// functions for track, cluster and candidate selection

// function that performs track selections on each track
template <typename T>
bool selectTrack(T const& track, std::string trackSelection)
{
  if (trackSelection == "globalTracks" && !track.isGlobalTrackWoPtEta()) {
    return false;
  } else if (trackSelection == "QualityTracks" && !track.isQualityTrack()) {
    return false;
  } else if (trackSelection == "hybridTracksJE" && !track.trackCutFlagFb5()) {//isQualityTrack
    return false;
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
      if constexpr (std::is_same_v<std::decay_t<U>, CandidateD0Data::iterator> || std::is_same_v<std::decay_t<U>, CandidateD0Data::filtered_iterator> || std::is_same_v<std::decay_t<U>, CandidateD0MC::iterator> || std::is_same_v<std::decay_t<U>, CandidateD0MC::filtered_iterator>) {
        if (cand.template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
      }

      if constexpr (std::is_same_v<std::decay_t<U>, CandidateLcData::iterator> || std::is_same_v<std::decay_t<U>, CandidateLcData::filtered_iterator> || std::is_same_v<std::decay_t<U>, CandidateLcMC::iterator> || std::is_same_v<std::decay_t<U>, CandidateLcMC::filtered_iterator>) {
        if (cand.template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong1_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong2_as<JetTracks>().globalIndex() == track.globalIndex()) {
          continue;
        }
      }

      if constexpr (std::is_same_v<std::decay_t<U>, CandidateBplusData::iterator> || std::is_same_v<std::decay_t<U>, CandidateBplusData::filtered_iterator> || std::is_same_v<std::decay_t<U>, CandidateBplusMC::iterator> || std::is_same_v<std::decay_t<U>, CandidateBplusMC::filtered_iterator>) {
        if (cand.template prong0_as<aod::HfCand2Prong>().template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong0_as<aod::HfCand2Prong>().template prong1_as<JetTracks>().globalIndex() == track.globalIndex() || cand.template prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
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
      jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
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
      break;
    }
  }
}

// function that checks if a candidate has any daughters that need to be removed from the event at gen level
template <typename T>
bool checkDaughters(T const& particle, int globalIndex)
{
  for (auto daughter : particle.template daughters_as<aod::McParticles>()) {
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
void analyseParticles(std::vector<fastjet::PseudoJet>& inputParticles, float particleEtaMin, float particleEtaMax, int jetTypeParticleLevel, T const& particles, TDatabasePDG* pdg, std::optional<U> const& candidate = std::nullopt)
{
  inputParticles.clear();
  for (auto& particle : particles) {
    // TODO: can we do this through the filter?
    if (particle.eta() < particleEtaMin || particle.eta() > particleEtaMax) {
      continue;
    }
    if (particle.getGenStatusCode() != 1) { // CHECK : Does this exclude the HF hadron?
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
bool selectCollision(T const& collision)
{
  if (!collision.sel8()) {
    return false;
  }
  return true;
}

#endif // PWGJE_TABLEPRODUCER_JETFINDER_H_
